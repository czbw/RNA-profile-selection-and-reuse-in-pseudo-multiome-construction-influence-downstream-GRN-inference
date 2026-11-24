#!/usr/bin/env python
"""The shared chromatin chain, step 1: peak matrix -> cisTopic object -> LDA topic model."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

# A region name pycisTopic can parse. Deliberately permissive about the contig token (GRCh38
# scaffolds such as KI270713.1 are legal) and strict about the ':' / '-' punctuation, which is
# the thing that actually breaks region_names_to_coordinates.
REGION_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")
# The deposited form, per config `dataset.peak_name_deposited`.
DEPOSITED_RE = re.compile(r"^([^:\s]+)-(\d+)-(\d+)$")
NNZ_SAMPLE = 1_000_000        # how many stored values to check for "raw counts, not normalised"

def log(m):
    print(f"[b_01] {m}", flush=True)

def die(msg):
    """Abort loudly. Nothing in this script continues past something impossible."""
    print(f"[b_01] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

# config

def load_config():
    """config.yaml, with `paths:` resolved against AMLGRN_ROOT."""
    name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    p = HERE / name
    if not p.is_file() or p.stat().st_size == 0:
        die(f"config {p} is missing or empty (AMLGRN_CONFIG={name}); run_bmmc.sh cds into "
            f"scripts/ and this script resolves the config relative to its own directory")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or "paths" not in cfg:
        die(f"{p} did not parse into a mapping with a `paths:` block")
    cfg["P"] = {k: (ROOT / v) for k, v in cfg["paths"].items()}
    cfg["CFG_PATH"] = p
    return cfg

# atomic write helpers — `.part.<ext>` then rename, exactly as run_bmmc.sh's partof()/commit()

def partof(final: Path) -> Path:
    """out/x.pkl -> out/x.part.pkl. The extension is preserved because anndata/mudata/polars
    all dispatch on it, and because run_bmmc.sh's own partof() does the same."""
    return final.with_suffix(".part" + final.suffix)

def commit_bytes(final: Path, payload: bytes):
    tmp = partof(final)
    with open(tmp, "wb") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)

def commit_pickle(final: Path, obj):
    tmp = partof(final)
    try:
        with open(tmp, "wb") as fh:
            pickle.dump(obj, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, final)
    except BaseException:
        # a killed write must leave nothing a later resubmit's `-s` guard would accept
        if tmp.exists():
            tmp.unlink()
        raise

# h5ad readers

def _index_of(g) -> np.ndarray:
    """String index of an AnnData dataframe group, honouring its `_index` attribute."""
    key = g.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[str(key)][:]],
                    dtype=object)

def _categorical(f, path) -> np.ndarray:
    """Read a categorical column in either layout: the older __categories layout (used by the
    deposit and read by census.py/export.py) or the modern obs/col/{codes,categories}."""
    grp, col = path.rsplit("/", 1)
    node = f[path]
    if isinstance(node, h5py.Group):                       # modern: obs/col/{codes,categories}
        codes = node["codes"][:]
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]]
    else:
        legacy = f"{grp}/__categories/{col}"
        if legacy not in f:
            return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in node[:]],
                            dtype=object)
        codes = node[:]
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f[legacy][:]]
    return np.array([cats[c] for c in codes], dtype=object)

def _sparse_from_group(g) -> sp.csr_matrix:
    """Materialise an AnnData X/layer group as CSR, whatever orientation it was stored in."""
    enc = g.attrs.get("encoding-type", "")
    if isinstance(enc, bytes):
        enc = enc.decode()
    shape = tuple(int(x) for x in g.attrs["shape"])
    data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
    if enc == "csc_matrix":
        return sp.csc_matrix((data, indices, indptr), shape=shape).tocsr()
    if enc not in ("csr_matrix", ""):
        die(f"unsupported sparse encoding-type {enc!r}")
    return sp.csr_matrix((data, indices, indptr), shape=shape)

def assert_raw_counts(M: sp.spmatrix, what: str):
    """Abort unless M looks like raw counts: X in this deposit is library-size-normalised + log1p
    and must not be summed, so check the values rather than the name of the dataset they came from."""
    if M.nnz == 0:
        die(f"{what}: matrix has zero stored values")
    v = M.data[:NNZ_SAMPLE]
    if float(v.min()) < 0:
        die(f"{what}: negative values present -- this is not a count matrix")
    if not np.all(np.equal(np.mod(v.astype("float64"), 1.0), 0.0)):
        die(f"{what}: non-integral values in the first {len(v):,} stored entries -- this looks "
            f"like the NORMALISED X, not raw counts. config.yaml's closing note: 'X in this "
            f"deposit is library-size-normalised + log1p and must NEVER be summed.'")
    log(f"{what}: raw-count check OK on {len(v):,} stored values "
        f"(min {float(v.min()):.0f}, max {float(v.max()):.0f})")

def normalise_peak_names(names: np.ndarray) -> np.ndarray:
    """Return peak names in the only form pycisTopic parses, or abort.

    Accepts the already-rewritten `chr1:9776-10668` (what export.py emits) and the
    deposited `chr1-9776-10668`, rewriting the first delimiter only. A mixture, or anything
    matching neither, is impossible for this deposit and therefore fatal rather than repaired.
    """
    n_ok = sum(1 for x in names if REGION_RE.match(str(x)))
    if n_ok == len(names):
        log(f"peak names already in pycisTopic form ({names[0]} ... {names[-1]})")
        return np.asarray(names, dtype=object)
    n_dep = sum(1 for x in names if DEPOSITED_RE.match(str(x)))
    if n_dep == len(names):
        out = np.array([str(x).replace("-", ":", 1) for x in names], dtype=object)
        bad = [x for x in out if not REGION_RE.match(str(x))]
        if bad:
            die(f"peak rewrite produced {len(bad)} unparseable names, e.g. {bad[:3]}")
        log(f"rewrote all {len(out):,} peak names  {names[0]} -> {out[0]}")
        return out
    bad = [str(x) for x in names if not REGION_RE.match(str(x)) and not DEPOSITED_RE.match(str(x))]
    die(f"peak names are neither `chr:start-end` ({n_ok}/{len(names)}) nor the deposited "
        f"`chr-start-end` ({n_dep}/{len(names)}); {len(bad)} match neither, e.g. {bad[:5]}. "
        f"config.yaml records the two legal forms as dataset.peak_name_deposited / "
        f"dataset.peak_name_required.")

def load_atac_presplit(path: Path):
    """out/bmmc_atac_counts.h5ad -> (cells x peaks CSR, cell ids, peak names)."""
    log(f"reading pre-split ATAC counts from {path}")
    with h5py.File(path, "r") as f:
        obs_names = _index_of(f["obs"])
        var_names = _index_of(f["var"])
        X = f["X"]
        if isinstance(X, h5py.Group):
            M = _sparse_from_group(X)
        else:
            log(f"WARNING: X is DENSE {X.shape} -- materialising "
                f"{np.prod(X.shape) * 4 / 1e9:.0f} GB before converting to CSR")
            M = sp.csr_matrix(X[:])
    if M.shape != (len(obs_names), len(var_names)):
        die(f"{path.name}: X {M.shape} does not match obs/var ({len(obs_names)}, {len(var_names)})")
    return M, obs_names, normalise_peak_names(var_names)

def load_atac_from_multiome(path: Path, block: int = 5000):
    """Derive the ATAC half from the raw GSE194122 multiome h5ad.

    A faithful port of export.py: the GEX/ATAC blocks are asserted contiguous (if they
    are ever interleaved the column slice below would emit a scrambled matrix, so it aborts),
    `layers/counts` is read rather than X, and rows are streamed in blocks so peak RSS stays a
    small multiple of one block rather than of the whole 325M-nonzero matrix.
    """
    log(f"deriving the ATAC half from the raw multiome {path}")
    with h5py.File(path, "r") as f:
        obs_names = _index_of(f["obs"])
        var_names = _index_of(f["var"])
        n_cells, n_feat = len(obs_names), len(var_names)
        if "var/feature_types" not in f:
            die(f"{path.name} has no var/feature_types -- cannot tell GEX from ATAC columns")
        ft = _categorical(f, "var/feature_types")
        gex_i = np.where(ft == "GEX")[0]
        atac_i = np.where(ft == "ATAC")[0]
        if len(gex_i) == 0 or len(atac_i) == 0:
            die(f"var/feature_types has {len(gex_i)} GEX and {len(atac_i)} ATAC columns")
        if not (gex_i.max() < atac_i.min()
                and np.array_equal(gex_i, np.arange(gex_i.min(), gex_i.max() + 1))
                and np.array_equal(atac_i, np.arange(atac_i.min(), atac_i.max() + 1))):
            die("GEX/ATAC blocks are NOT contiguous -- the column slice would scramble features")
        split = int(atac_i.min())
        log(f"{n_cells:,} cells x {n_feat:,} features | GEX 0..{split - 1} | "
            f"ATAC {split}..{n_feat - 1}")

        if "layers" in f and "counts" in f["layers"]:
            src = f["layers/counts"]
            log("reading layers/counts (NOT X: X is normalised + log1p)")
        else:
            die(f"{path.name} has no layers/counts. X in this deposit is library-size-normalised "
                f"+ log1p and must never be summed, so there is no safe fallback.")
        indptr, data, indices = src["indptr"], src["data"], src["indices"]
        blocks = []
        for s in range(0, n_cells, block):
            e = min(s + block, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[lo:hi], indices[lo:hi], indptr[s:e + 1] - lo),
                                shape=(e - s, n_feat))
            blocks.append(blk[:, split:])
            if (s // block) % 5 == 0:
                log(f"  rows {s:,}-{e:,}")
    M = sp.vstack(blocks).tocsr()
    del blocks
    return M, obs_names, normalise_peak_names(var_names[atac_i])

def digest_strings(arr) -> str:
    h = hashlib.sha256()
    for x in arr:
        h.update(str(x).encode())
    return h.hexdigest()[:16]

def verify_existing(out_pkl: Path, n_cells_exp: int, n_regions_exp: int, group_by: str):
    """Resumability, done honestly: an artifact counts as complete only if its content is
    complete. `-s` proves a file is non-empty, not that a pickle finished writing, so the
    existing object is unpickled and its shape, its grouping column and its fitted LDA model
    are all checked. A truncated pickle raises here instead of being frozen into the chain."""
    log(f"{out_pkl.name} exists ({out_pkl.stat().st_size / 1e9:.2f} GB) -- verifying content "
        f"rather than trusting its size")
    try:
        with open(out_pkl, "rb") as fh:
            obj = pickle.load(fh)
    except Exception as e:
        die(f"{out_pkl} exists but does NOT unpickle ({type(e).__name__}: {e}). That is a "
            f"TRUNCATED artifact -- delete it and resubmit; do not let chain_guard freeze it.")
    n_c, n_r = len(obj.cell_names), len(obj.region_names)
    if (n_c, n_r) != (n_cells_exp, n_regions_exp):
        die(f"{out_pkl.name} is {n_c:,} x {n_r:,}, config declares {n_cells_exp:,} x "
            f"{n_regions_exp:,}. Refusing to reuse a chain of the wrong shape.")
    if getattr(obj, "selected_model", None) is None:
        die(f"{out_pkl.name} carries NO fitted LDA model (selected_model is None) -- it was "
            f"written before add_LDA_model. Delete it and resubmit.")
    if group_by not in obj.cell_data.columns:
        die(f"{out_pkl.name} cell_data has no `{group_by}` column -- region_sets.py calls "
            f"find_diff_features on it and would fail hours later. Delete it and resubmit.")
    log(f"verified: {n_c:,} cells x {n_r:,} regions, LDA model present, `{group_by}` present")
    return obj

def write_cell_names(path: Path, names) -> str:
    payload = ("\n".join(str(x) for x in names) + "\n").encode("utf-8")
    commit_bytes(path, payload)
    return digest_strings(names)

def main():
    ap = argparse.ArgumentParser(description="BMMC shared chromatin chain: cisTopic + LDA")
    ap.add_argument("--h5ad", required=True,
                    help="raw GSE194122 multiome h5ad (data/multiome_BMMC.h5ad). Used only if "
                         "out/bmmc_atac_counts.h5ad is absent; always recorded as provenance.")
    ap.add_argument("--metadata", required=True,
                    help="uploads/b00_cell_metadata.csv -- the verified per-cell manifest")
    ap.add_argument("--stage", default=os.environ.get("AMLGRN_STAGE", "bmmc"))
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if a verified cistopic object already exists")
    a = ap.parse_args()

    cfg = load_config()
    STAGE = a.stage
    OUTP = cfg["P"]["out"]
    WORK = cfg["P"]["work"]
    OUTP.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    SEP = str(cfg["cell_id_sep"])
    SEED = int(os.environ.get("AMLGRN_SEED", cfg["seed"]))
    GROUP_BY = cfg["pseudobulk"]["group_by"]
    DS = cfg["dataset"]
    N_CELLS, N_PEAKS = int(DS["n_cells"]), int(DS["n_atac_peaks"])
    TM = cfg["topic_model"]
    TOPIC_NCPU = int(cfg.get("topic_n_cpu", 3))

    out_pkl = OUTP / f"cistopic_obj_{STAGE}.pkl"
    out_ids = OUTP / f"cistopic_cell_names_{STAGE}.txt"
    out_rep = OUTP / f"b_01_chain_report_{STAGE}.json"

    log(f"stage {STAGE} | seed {SEED} | config {cfg['CFG_PATH'].name} | root {ROOT}")
    log(f"topic model: backend={TM['backend']} n_topics={TM['n_topics']} n_iter={TM['n_iter']} "
        f"alpha={TM['alpha']} eta={TM['eta']} | {TOPIC_NCPU} models concurrently")

    if TM["backend"] != "lda":
        die(f"config topic_model.backend is {TM['backend']!r}; only the VALIDATED `lda` "
            f"(collapsed Gibbs) backend is wired here. Configuration required: port the Mallet path "
            f"from the reference bundles before setting anything else.")

    # ---- already done? ----------------------------------------------------------------------
    if out_pkl.exists() and out_pkl.stat().st_size > 0 and not a.force:
        obj = verify_existing(out_pkl, N_CELLS, N_PEAKS, GROUP_BY)
        d = write_cell_names(out_ids, obj.cell_names)
        log(f"SKIP: the chromatin chain's cisTopic object is already built and complete "
            f"({out_pkl.name}); rewrote {out_ids.name} (digest {d}) and did nothing else.")
        return
    for stale in (partof(out_pkl), partof(out_ids)):
        if stale.exists():
            log(f"removing stale sidecar {stale.name} (a killed write leaves no valid output)")
            stale.unlink()

    # ---- the peak x cell matrix -------------------------------------------------------------
    presplit = OUTP / "bmmc_atac_counts.h5ad"
    src_h5ad = Path(a.h5ad)
    if presplit.is_file() and presplit.stat().st_size > 0:
        M, cell_ids, region_names = load_atac_presplit(presplit)
        source = str(presplit)
    else:
        if not src_h5ad.is_file() or src_h5ad.stat().st_size == 0:
            die(f"neither {presplit} nor --h5ad {src_h5ad} is present and non-empty. "
                f"The raw multiome is expected at data/multiome_BMMC.h5ad; "
                f"the pre-split ATAC matrix is kept local.")
        M, cell_ids, region_names = load_atac_from_multiome(src_h5ad)
        source = str(src_h5ad)

    log(f"ATAC {M.shape[0]:,} cells x {M.shape[1]:,} peaks | nnz {M.nnz:,} | dtype {M.dtype}")
    assert_raw_counts(M, "ATAC counts")
    # `lda`, the collapsed-Gibbs backend, requires an integer dtype: lda.utils.matrix_to_lists
    # raises "expected sparse matrix with integer values, found float values" on float input even
    # when every stored value is integral. The fragment-based path (AML, PBMC) yields int
    # natively; this count-matrix path inherits float32 from the deposit's layers/counts. The cast
    # is lossless and placed after assert_raw_counts on purpose: that check has just shown the
    # values are non-negative and integral, so the cast cannot hide a normalised matrix. Observed max is 68 — nowhere near int32's range.
    if M.dtype.kind != "i":
        log(f"casting ATAC counts {M.dtype} -> int32 for the `lda` backend "
            f"(values already verified integral above)")
        M = M.astype(np.int32)
    if M.shape != (N_CELLS, N_PEAKS):
        die(f"ATAC matrix is {M.shape[0]:,} x {M.shape[1]:,}; config.yaml declares "
            f"{N_CELLS:,} cells x {N_PEAKS:,} peaks. Those are MEASURED counts, not targets -- "
            f"a different number means a different deposit, so this aborts rather than "
            f"silently building a chain nothing else in the bundle agrees with.")
    if len(set(map(str, cell_ids))) != len(cell_ids):
        die("cell ids in the ATAC matrix are not unique -- the h5ad barcode is supposed to be "
            "globally unique (the bare fragment barcode is NOT, and is never keyed on)")
    if len(set(map(str, region_names))) != len(region_names):
        die("peak names are not unique after normalisation")

    # ---- the per-cell manifest --------------------------------------------------------------
    meta_p = Path(a.metadata)
    if not meta_p.is_file() or meta_p.stat().st_size == 0:
        die(f"--metadata {meta_p} is missing or empty")
    meta = pd.read_csv(meta_p, dtype=str)
    need = ["h5ad_barcode", "fragment_barcode", "batch", "donor", "site", "cell_type"]
    missing_cols = [c for c in need if c not in meta.columns]
    if missing_cols:
        die(f"{meta_p.name} is missing column(s) {missing_cols}; it has {list(meta.columns)}")
    if meta["h5ad_barcode"].duplicated().any():
        die(f"{meta_p.name} has duplicate h5ad_barcode values -- it must be one row per cell")
    meta = meta.set_index("h5ad_barcode")
    only_matrix = set(map(str, cell_ids)) - set(meta.index)
    only_meta = set(meta.index) - set(map(str, cell_ids))
    if only_matrix or only_meta:
        die(f"CELL-ID MISMATCH between the ATAC matrix and {meta_p.name}: "
            f"{len(only_matrix)} in the matrix only (e.g. {sorted(only_matrix)[:3]}), "
            f"{len(only_meta)} in the manifest only (e.g. {sorted(only_meta)[:3]}). Both must "
            f"be the same {N_CELLS:,} cells, keyed on the h5ad barcode.")
    log(f"manifest {len(meta):,} rows covers the ATAC matrix exactly "
        f"({meta['batch'].nunique()} batches, {meta['donor'].nunique()} donors, "
        f"{meta['site'].nunique()} sites, {meta['cell_type'].nunique()} cell types)")
    for key, want in (("batch", "n_batches"), ("donor", "n_donors"),
                      ("site", "n_sites"), ("cell_type", "n_cell_types")):
        if int(meta[key].nunique()) != int(DS[want]):
            die(f"{meta_p.name} has {meta[key].nunique()} distinct {key} values; config "
                f"dataset.{want} = {DS[want]}")

    # ---- cisTopic object --------------------------------------------------------------------
    from pycisTopic.cistopic_class import create_cistopic_object
    from pycisTopic.lda_models import run_cgs_models, evaluate_models

    clock = {}
    t0 = time.time()
    # pycisTopic wants regions x cells. The transpose is the one large allocation here.
    log("transposing to regions x cells for create_cistopic_object")
    MT = M.T.tocsr()
    del M
    obj = create_cistopic_object(
        fragment_matrix=MT,
        cell_names=[str(x) for x in cell_ids],
        region_names=[str(x) for x in region_names],
        project=STAGE,
        split_pattern=SEP,
    )
    del MT
    clock["cistopic_object"] = round(time.time() - t0, 1)
    log(f"cisTopic object: {len(obj.cell_names):,} cells x {len(obj.region_names):,} regions "
        f"({clock['cistopic_object']}s)")

    # create_cistopic_object filters cells (min_frag) and regions (min_cell). On this deposit
    # nothing should be dropped: it is already filtered and annotated by the NeurIPS organisers
    # and the consensus peak set is fixed. If anything is dropped the analysis set and the peak
    # set both move, so it is fatal, not a warning.
    if len(obj.cell_names) != N_CELLS:
        die(f"create_cistopic_object kept {len(obj.cell_names):,} of {N_CELLS:,} cells. The "
            f"deposit is pre-filtered and no cell may be dropped: every arm's assignment table "
            f"covers all {N_CELLS:,}. Investigate before proceeding (Configuration required: decide "
            f"whether the analysis set shrinks or the input is wrong).")
    if len(obj.region_names) != N_PEAKS:
        die(f"create_cistopic_object kept {len(obj.region_names):,} of {N_PEAKS:,} peaks. The "
            f"NeurIPS consensus peak set is fixed and must be used as shipped.")

    # ---- cell ids: read back what the object actually named them ----------------------------
    cis_names = [str(x) for x in obj.cell_names]
    base = [n.split(SEP)[0] for n in cis_names]
    if len(set(base)) != len(base):
        die(f"splitting cisTopic cell names on {SEP!r} does not give unique barcodes -- the "
            f"cell_id_sep is wrong (a '-' separator silently matches zero cells)")
    if set(base) != set(map(str, cell_ids)):
        die(f"cisTopic cell names do not map back onto the input barcodes when split on {SEP!r}; "
            f"example {cis_names[0]!r} -> {base[0]!r}")
    tagged = sum(1 for n in cis_names if SEP in n)
    log(f"cisTopic cell ids: {cis_names[0]} ({tagged:,}/{len(cis_names):,} carry the "
        f"{SEP!r}{STAGE} tag) -- this is the ANALYSIS SET b_03 must reproduce exactly")

    # ---- cell_data, incl. the grouping key run_grn_bmmc.sh needs ----------------------------
    cd = meta.loc[base].copy()
    cd.index = pd.Index(cis_names, name=None)
    # The census emits the label as `cell_type`; b_02's find_diff_features and run_grn_bmmc.sh's
    # --key_to_group_by both want Classified_Celltype. Copying rather than renaming keeps the
    # raw column too.
    cd[GROUP_BY] = cd["cell_type"].astype(str)
    obj.add_cell_data(cd, split_pattern=SEP)
    if GROUP_BY not in obj.cell_data.columns:
        die(f"add_cell_data did not attach `{GROUP_BY}`; columns are "
            f"{list(obj.cell_data.columns)}")
    n_na = int(obj.cell_data[GROUP_BY].isna().sum())
    if n_na:
        die(f"{n_na:,} cells have no `{GROUP_BY}` after add_cell_data -- the join key is wrong")
    log(f"cell_data attached: {obj.cell_data.shape[0]:,} rows, `{GROUP_BY}` complete "
        f"({obj.cell_data[GROUP_BY].nunique()} levels)")

    try:                      # a headless node has no display; evaluate_models draws
        import matplotlib
        matplotlib.use("Agg")
    except Exception as e:    # matplotlib absent or already configured -- not fatal
        log(f"(matplotlib backend not forced: {e})")

    models_dir = OUTP / f"models_{STAGE}"
    models_dir.mkdir(parents=True, exist_ok=True)

    # Whatever matrix run_cgs_models hands to `lda` must be integer-typed. The input was cast to
    # int32 above, but create_cistopic_object builds its own matrix and its dtype is not controlled
    # here, so check it rather than assume the cast survived (a float matrix would otherwise fail
    # deep inside the model-fitting workers).
    for _attr in ("binary_matrix", "fragment_matrix"):
        _mtx = getattr(obj, _attr, None)
        _dt = getattr(_mtx, "dtype", None)
        if _dt is not None:
            log(f"  {_attr}.dtype = {_dt}")
            if _dt.kind not in "iub":
                die(f"cisTopic object's `{_attr}` is {_dt}, but the `lda` backend requires an "
                    f"integer dtype and will raise 'expected sparse matrix with integer values, "
                    f"found float values'. The ATAC counts WERE cast to int32 before "
                    f"create_cistopic_object, so pycisTopic re-cast them internally. Fix that "
                    f"conversion -- do NOT binarise here, because the topic model is fitted on "
                    f"counts and silently swapping in a 0/1 matrix would change the result "
                    f"without changing the log.")

    t0 = time.time()
    log(f"run_cgs_models: {len(TM['n_topics'])} models -> {models_dir}")
    models = run_cgs_models(
        obj,
        n_topics=TM["n_topics"],
        n_cpu=TOPIC_NCPU,
        n_iter=int(TM["n_iter"]),
        random_state=SEED,
        alpha=TM["alpha"],
        eta=TM["eta"],
        save_path=str(models_dir),
    )
    clock["topic_models"] = round(time.time() - t0, 1)
    if not models:
        die("run_cgs_models returned no models")
    log(f"{len(models)} topic models fitted ({clock['topic_models']}s)")

    best = evaluate_models(models, return_model=True)
    if best is None:
        die("evaluate_models returned no model")
    obj.add_LDA_model(best)
    if getattr(obj, "selected_model", None) is None:
        die("add_LDA_model did not set selected_model")
    n_topics_sel = getattr(best, "n_topic", None) or getattr(best, "n_topics", None)
    log(f"selected model: {n_topics_sel} topics")

    # ---- commit: ids first, object last (run_bmmc.sh guards on the object) -------------------
    id_digest = write_cell_names(out_ids, cis_names)
    commit_pickle(out_pkl, obj)
    log(f"committed {out_ids.name} (digest {id_digest}) and {out_pkl.name} "
        f"({out_pkl.stat().st_size / 1e9:.2f} GB)")

    rep = {
        "stage": STAGE,
        "seed": SEED,
        "source_matrix": source,
        "provenance_h5ad": str(src_h5ad),
        "metadata": str(meta_p),
        "n_cells": len(cis_names),
        "n_regions": len(obj.region_names),
        "cell_id_sep": SEP,
        "cell_id_example": cis_names[0],
        "cell_id_sha256_16": id_digest,
        "region_sha256_16": digest_strings(obj.region_names),
        "group_by": GROUP_BY,
        "topic_model": {"backend": TM["backend"], "n_topics": TM["n_topics"],
                        "n_iter": TM["n_iter"], "alpha": TM["alpha"], "eta": TM["eta"],
                        "n_cpu_models": TOPIC_NCPU, "selected_n_topics": n_topics_sel},
        "timings_s": clock,
    }
    commit_bytes(out_rep, (json.dumps(rep, indent=1) + "\n").encode("utf-8"))
    log(f"=== timings (s): {clock} ===")
    log(f"b_01 COMPLETE -> {out_pkl}")

if __name__ == "__main__":
    main()
