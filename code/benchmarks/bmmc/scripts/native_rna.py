#!/usr/bin/env python
"""The native (unpaired) RNA object, the input the metacell arms consume."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

NNZ_SAMPLE = 1_000_000

def log(m):
    print(f"[b_10] {m}", flush=True)

def die(msg):
    print(f"[b_10] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def load_config():
    """See atac_chain.py: this bundle ships no common.py, so load_config()'s contract
    (read AMLGRN_CONFIG from scripts/, resolve `paths:` against AMLGRN_ROOT) is inlined."""
    name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    p = HERE / name
    if not p.is_file() or p.stat().st_size == 0:
        die(f"config {p} is missing or empty (AMLGRN_CONFIG={name})")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or "paths" not in cfg:
        die(f"{p} did not parse into a mapping with a `paths:` block")
    cfg["P"] = {k: (ROOT / v) for k, v in cfg["paths"].items()}
    cfg["CFG_PATH"] = p
    return cfg

def partof(final: Path) -> Path:
    """out/native_rna_bmmc.h5ad -> out/native_rna_bmmc.part.h5ad. The extension is preserved
    because anndata dispatches on it (identical to run_bmmc.sh's partof())."""
    return final.with_suffix(".part" + final.suffix)

def commit_bytes(final: Path, payload: bytes):
    tmp = partof(final)
    with open(tmp, "wb") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()

def digest_strings(arr) -> str:
    h = hashlib.sha256()
    for x in arr:
        h.update(str(x).encode())
    return h.hexdigest()[:16]

# h5ad readers — the same helpers b_01/b_03 carry; this bundle ships no shared module

def _index_of(g) -> np.ndarray:
    key = g.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[str(key)][:]],
                    dtype=object)

def _categorical(f, path) -> np.ndarray:
    grp, col = path.rsplit("/", 1)
    node = f[path]
    if isinstance(node, h5py.Group):
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

def assert_raw_counts(M: sp.spmatrix, what: str) -> dict:
    """X in this deposit is library-size-normalised + log1p and must never be summed. Check the
    values, not the name of the dataset they were read from."""
    if M.nnz == 0:
        die(f"{what}: matrix has zero stored values")
    v = M.data[:NNZ_SAMPLE]
    if float(v.min()) < 0:
        die(f"{what}: negative values present -- this is not a count matrix")
    if not np.all(np.equal(np.mod(v.astype("float64"), 1.0), 0.0)):
        die(f"{what}: non-integral values in the first {len(v):,} stored entries -- this looks "
            f"like the NORMALISED X, not raw counts. Averaging a log1p matrix into metacells "
            f"would silently compute a geometric mean.")
    log(f"{what}: raw-count check OK on {len(v):,} stored values "
        f"(min {float(v.min()):.0f}, max {float(v.max()):.0f})")
    return {"n_checked": int(len(v)), "min": float(v.min()), "max": float(v.max()),
            "all_integral": True, "any_negative": False}

def load_rna_presplit(path: Path):
    log(f"reading pre-split RNA counts from {path}")
    with h5py.File(path, "r") as f:
        obs_names = _index_of(f["obs"])
        var_names = _index_of(f["var"])
        X = f["X"]
        if isinstance(X, h5py.Group):
            M = _sparse_from_group(X)
        else:
            log(f"WARNING: X is DENSE {X.shape} -- materialising it")
            M = sp.csr_matrix(X[:])
    if M.shape != (len(obs_names), len(var_names)):
        die(f"{path.name}: X {M.shape} does not match obs/var "
            f"({len(obs_names)}, {len(var_names)})")
    return M, obs_names, var_names, "X (pre-split; export.py wrote layers/counts here)"

def load_rna_from_multiome(path: Path, block: int = 5000):
    """Derive the GEX half from the raw GSE194122 multiome h5ad -- the same port of
    export.py that build_arm_rna.py carries, including the contiguity assertion
    and the layers/counts requirement."""
    log(f"deriving the GEX half from the raw multiome {path}")
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
        log(f"{n_cells:,} cells x {n_feat:,} features | GEX 0..{split - 1}")
        if "layers" in f and "counts" in f["layers"]:
            src = f["layers/counts"]
            log("reading layers/counts (NOT X: X is normalised + log1p)")
        else:
            die(f"{path.name} has no layers/counts, and X must never be summed -- no fallback")
        indptr, data, indices = src["indptr"], src["data"], src["indices"]
        blocks = []
        for s in range(0, n_cells, block):
            e = min(s + block, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[lo:hi], indices[lo:hi], indptr[s:e + 1] - lo),
                                shape=(e - s, n_feat))
            blocks.append(blk[:, :split])
            if (s // block) % 5 == 0:
                log(f"  rows {s:,}-{e:,}")
    M = sp.vstack(blocks).tocsr()
    del blocks
    return M, obs_names, var_names[gex_i], "layers/counts of the raw multiome deposit"

def resolve_cell_ids(outp: Path, stage: str, sep: str, barcodes: np.ndarray):
    """The ids to write, in the order of `barcodes` (the RNA matrix's own order).

    Returns (ids, base, source, verified). `verified` is True only when the ids came from the
    frozen chain; when it is False the caller must record that, and the assertion against the
    real analysis set is deferred to metacell_floor.py --analysis-set on CUBIC.
    """
    ids_p = outp / f"cistopic_cell_names_{stage}.txt"
    if ids_p.is_file() and ids_p.stat().st_size > 0:
        names = [ln.strip() for ln in ids_p.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if len(set(names)) != len(names):
            die(f"{ids_p.name} contains duplicate cell ids")
        base = [n.split(sep)[0] for n in names]
        if len(set(base)) != len(base):
            die(f"splitting {ids_p.name} on cell_id_sep {sep!r} does not give unique barcodes")
        pos = {b: i for i, b in enumerate(base)}
        missing = [b for b in barcodes if b not in pos]
        extra = sorted(set(base) - set(map(str, barcodes)))
        if missing or extra:
            die(f"the analysis set in {ids_p.name} does not cover the RNA matrix: "
                f"{len(missing)} RNA-only (e.g. {list(missing)[:3]}), "
                f"{len(extra)} chain-only (e.g. {extra[:3]})")
        ids = [names[pos[b]] for b in barcodes]
        log(f"cell ids: {len(ids):,} taken FROM THE FROZEN CHAIN ({ids_p.name}); "
            f"example {ids[0]}")
        return ids, [str(b) for b in barcodes], str(ids_p), True

    ids = [f"{b}{sep}{stage}" for b in barcodes]
    log(f"cell ids: {ids_p.name} is ABSENT -- composing <h5ad_barcode>{sep!r}{stage} instead, "
        f"which is what create_cistopic_object(project={stage!r}, split_pattern={sep!r}) "
        f"produces and what b_01 asserts it produced. Example {ids[0]}")
    log("cell ids: analysis_set_verified=FALSE. Assert this on CUBIC, where the chain exists:")
    log(f"cell ids:   python metacell_floor.py --gex out/native_rna_{stage}.metacell.h5ad "
        f"--stage {stage} --k <k> --analysis-set out/{ids_p.name} "
        f"--metadata uploads/b00_cell_metadata.csv")
    return ids, [str(b) for b in barcodes], f"convention <h5ad_barcode>{sep}{stage} " \
                                            f"({ids_p.name} absent)", False

def main():
    ap = argparse.ArgumentParser(description="Build the BMMC native (unpaired) RNA object")
    ap.add_argument("--stage", default=os.environ.get("AMLGRN_STAGE", "bmmc"))
    ap.add_argument("--h5ad", required=True,
                    help="raw GSE194122 multiome h5ad; used only if out/bmmc_rna_counts.h5ad is "
                         "absent, always recorded as provenance")
    ap.add_argument("--metadata", required=True, help="uploads/b00_cell_metadata.csv "
                                                      "(or out/b00_cell_metadata.csv locally)")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_config()
    STAGE = a.stage
    OUTP = cfg["P"]["out"]
    OUTP.mkdir(parents=True, exist_ok=True)
    SEP = str(cfg["cell_id_sep"])
    GROUP_BY = cfg["pseudobulk"]["group_by"]
    DS = cfg["dataset"]

    out_h5ad = OUTP / f"native_rna_{STAGE}.h5ad"
    out_rep = OUTP / f"b_10_native_rna_report_{STAGE}.json"
    if out_h5ad.is_file() and out_h5ad.stat().st_size > 0 and not a.force:
        log(f"SKIP: {out_h5ad.name} already exists ({out_h5ad.stat().st_size / 1e6:.0f} MB)")
        return
    stale = partof(out_h5ad)
    if stale.exists():
        log(f"removing stale sidecar {stale.name} (a killed write leaves no valid output)")
        stale.unlink()

    log(f"stage {STAGE} | config {cfg['CFG_PATH'].name} | root {ROOT}")

    # ---- 1) the per-cell manifest -----------------------------------------------------------
    meta_p = Path(a.metadata)
    if not meta_p.is_file() or meta_p.stat().st_size == 0:
        die(f"--metadata {meta_p} is missing or empty")
    meta = pd.read_csv(meta_p, dtype=str)
    need = ["h5ad_barcode", "fragment_barcode", "batch", "donor", "site", "cell_type"]
    missing_cols = [c for c in need if c not in meta.columns]
    if missing_cols:
        die(f"{meta_p.name} is missing column(s) {missing_cols}; it has {list(meta.columns)}")
    if meta["h5ad_barcode"].duplicated().any():
        die(f"{meta_p.name} has duplicate h5ad_barcode values -- one row per cell is required")
    meta = meta.set_index("h5ad_barcode")
    for key, want in (("batch", "n_batches"), ("donor", "n_donors"),
                      ("site", "n_sites"), ("cell_type", "n_cell_types")):
        if int(meta[key].nunique()) != int(DS[want]):
            die(f"{meta_p.name} has {meta[key].nunique()} distinct {key} values; config "
                f"dataset.{want} = {DS[want]}")
    log(f"manifest {len(meta):,} rows | {meta['batch'].nunique()} batches, "
        f"{meta['donor'].nunique()} donors, {meta['site'].nunique()} sites, "
        f"{meta['cell_type'].nunique()} cell types")

    # ---- 2) the RNA matrix ------------------------------------------------------------------
    presplit = OUTP / "bmmc_rna_counts.h5ad"
    src_h5ad = Path(a.h5ad)
    if presplit.is_file() and presplit.stat().st_size > 0:
        RNA, rna_ids, gene_names, layer = load_rna_presplit(presplit)
        source = presplit
    else:
        if not src_h5ad.is_file() or src_h5ad.stat().st_size == 0:
            die(f"neither {presplit} nor --h5ad {src_h5ad} is present and non-empty")
        RNA, rna_ids, gene_names, layer = load_rna_from_multiome(src_h5ad)
        source = src_h5ad
    rna_ids = np.array([str(x) for x in rna_ids], dtype=object)
    log(f"RNA {RNA.shape[0]:,} cells x {RNA.shape[1]:,} genes | nnz {RNA.nnz:,} | layer: {layer}")
    counts_check = assert_raw_counts(RNA, "RNA counts")
    if RNA.shape != (int(DS["n_cells"]), int(DS["n_gex_features"])):
        die(f"RNA matrix is {RNA.shape[0]:,} x {RNA.shape[1]:,}; config declares "
            f"{int(DS['n_cells']):,} x {int(DS['n_gex_features']):,}")
    if len(set(rna_ids)) != len(rna_ids):
        die("RNA cell ids are not unique")
    if set(rna_ids) != set(meta.index):
        die(f"CELL-ID MISMATCH between the RNA matrix and {meta_p.name}: "
            f"{len(set(rna_ids) - set(meta.index))} RNA-only, "
            f"{len(set(meta.index) - set(rna_ids))} manifest-only. Both modalities are keyed on "
            f"the h5ad barcode.")

    # ---- 3) cell ids: from the chain when it exists, else the asserted convention ------------
    cell_ids, base, id_source, id_verified = resolve_cell_ids(OUTP, STAGE, SEP, rna_ids)

    # ---- 4) obs, var, write -----------------------------------------------------------------
    import anndata as ad

    rows = meta.loc[base]
    obs = pd.DataFrame(index=pd.Index(cell_ids, name=None))
    for c in ("fragment_barcode", "batch", "donor", "site", "cell_type"):
        obs[c] = rows[c].values
    obs["atac_barcode"] = np.array(base, dtype=object)
    obs["sample_id"] = STAGE
    # The grouping key run_grn_bmmc.sh passes as --key_to_group_by and config pseudobulk.group_by
    # declares. Copied from `cell_type` rather than renamed, as b_01 does when it builds the
    # cisTopic cell_data, so the two modalities carry the same 22 labels, spelled the same way.
    obs[GROUP_BY] = rows["cell_type"].astype(str).values
    if obs[GROUP_BY].isna().any():
        die(f"`{GROUP_BY}` has {int(obs[GROUP_BY].isna().sum())} missing values")

    var = pd.DataFrame(index=pd.Index([str(g) for g in gene_names], name=None))
    adata = ad.AnnData(X=RNA.tocsr(), obs=obs, var=var)
    adata.uns["bmmc_native"] = {
        "stage": STAGE, "rna_source": str(source), "rna_layer": layer,
        "cell_id_source": id_source, "analysis_set_verified": bool(id_verified),
        "cell_id_sep": SEP, "group_by": GROUP_BY,
        "note": "unpaired native RNA -- the AVERAGING family's input. No assignment exists.",
    }

    tmp = partof(out_h5ad)
    try:
        adata.write_h5ad(tmp, compression="gzip")
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise
    if not tmp.is_file() or tmp.stat().st_size == 0:
        die("write_h5ad returned but produced no output")
    os.replace(tmp, out_h5ad)
    log(f"committed {out_h5ad.name} ({out_h5ad.stat().st_size / 1e6:.0f} MB)")

    vc = obs[GROUP_BY].value_counts().sort_values()
    log(f"per-{GROUP_BY} counts: {len(vc)} levels, smallest "
        f"'{vc.index[0]}' n={int(vc.iloc[0]):,}, largest '{vc.index[-1]}' n={int(vc.iloc[-1]):,}")

    # ---- 5) report --------------------------------------------------------------------------
    inputs = {str(meta_p): sha256_file(meta_p), str(source): sha256_file(source)}
    if source != src_h5ad and src_h5ad.is_file():
        inputs[str(src_h5ad)] = "not read (pre-split GEX present); recorded as provenance only"
    rep = {
        "script": "native_rna.py",
        "stage": STAGE,
        "config": str(cfg["CFG_PATH"]),
        "output": str(out_h5ad),
        "output_sha256": sha256_file(out_h5ad),
        "output_bytes": int(out_h5ad.stat().st_size),
        "inputs_sha256": inputs,
        "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
        "rna_layer": layer,
        "raw_count_check": counts_check,
        "cell_id_source": id_source,
        "cell_id_example": cell_ids[0],
        "cell_id_sep": SEP,
        "cell_id_sha256_16": digest_strings(cell_ids),
        "analysis_set_verified": bool(id_verified),
        "analysis_set_verification_note": (
            "ids came from the frozen chain" if id_verified else
            f"cistopic_cell_names_{STAGE}.txt was ABSENT at build time; the ids follow the "
            f"convention b_01 asserts create_cistopic_object produces. The assertion against "
            f"the real analysis set is made on CUBIC by "
            f"metacell_floor.py --analysis-set."),
        "group_by": GROUP_BY,
        "celltype_counts": {str(k): int(v) for k, v in
                            obs[GROUP_BY].value_counts().sort_index().items()},
    }
    commit_bytes(out_rep, (json.dumps(rep, indent=1) + "\n").encode("utf-8"))
    log(f"wrote {out_rep.name}")
    log("b_10 COMPLETE")

if __name__ == "__main__":
    main()
