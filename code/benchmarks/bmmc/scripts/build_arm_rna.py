#!/usr/bin/env python
"""Build one arm's paired RNA object: the RNA matrix, permuted onto the ATAC cells."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

# The assignment table's column names, in preference order: atac_barcode / assigned_rna_barcode
# (written by arms.py) or cistopic_id / Linked_scRNA_barcode (PBMC bundle). Exactly one candidate from each list must be present --
# anything else is reported with the actual header rather than guessed at.
ATAC_COLS = ["atac_barcode", "atac_cell_id", "cistopic_id", "cell_id"]
RNA_COLS = ["assigned_rna_barcode", "rna_barcode", "assigned_rna_cell_id",
            "Linked_scRNA_barcode"]
NNZ_SAMPLE = 1_000_000

def log(m):
    print(f"[b_03] {m}", flush=True)

def die(msg):
    print(f"[b_03] FATAL: {msg}", file=sys.stderr, flush=True)
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
    "out/paired_rna_bmmc.native reference.h5ad -> out/paired_rna_bmmc.native reference.part.h5ad. Identical to\n    run_bmmc.sh's partof(): the extension is preserved because anndata dispatches on it."
    return final.with_suffix(".part" + final.suffix)

def commit_bytes(final: Path, payload: bytes):
    tmp = partof(final)
    with open(tmp, "wb") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)

def digest_strings(arr) -> str:
    h = hashlib.sha256()
    for x in arr:
        h.update(str(x).encode())
    return h.hexdigest()[:16]

# h5ad readers (same helpers as b_01; this bundle ships no shared module)

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

def assert_raw_counts(M: sp.spmatrix, what: str):
    """X in this deposit is library-size-normalised + log1p and must not be summed. Check the
    values, not the name of the dataset they were read from."""
    if M.nnz == 0:
        die(f"{what}: matrix has zero stored values")
    v = M.data[:NNZ_SAMPLE]
    if float(v.min()) < 0:
        die(f"{what}: negative values present -- this is not a count matrix")
    if not np.all(np.equal(np.mod(v.astype("float64"), 1.0), 0.0)):
        die(f"{what}: non-integral values in the first {len(v):,} stored entries -- this looks "
            f"like the NORMALISED X, not raw counts.")
    log(f"{what}: raw-count check OK on {len(v):,} stored values "
        f"(min {float(v.min()):.0f}, max {float(v.max()):.0f})")

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
    return M, obs_names, var_names

def load_rna_from_multiome(path: Path, block: int = 5000):
    """Derive the GEX half from the raw GSE194122 multiome h5ad -- a faithful port of
    export.py, including its contiguity assertion and its layers/counts requirement."""
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
    return M, obs_names, var_names[gex_i]

def analysis_set(outp: Path, stage: str, sep: str):
    """The cisTopic cell ids, in cisTopic order. From b_01's sidecar when it exists, otherwise
    from the pickle itself -- never from a naming convention."""
    ids_p = outp / f"cistopic_cell_names_{stage}.txt"
    cis_p = outp / f"cistopic_obj_{stage}.pkl"
    if ids_p.is_file() and ids_p.stat().st_size > 0:
        names = [ln.strip() for ln in ids_p.read_text(encoding="utf-8").splitlines()
                 if ln.strip()]
        log(f"analysis set: {len(names):,} cell ids from {ids_p.name}")
    elif cis_p.is_file() and cis_p.stat().st_size > 0:
        log(f"{ids_p.name} absent -- reading cell names from {cis_p.name} instead "
            f"({cis_p.stat().st_size / 1e9:.2f} GB)")
        try:
            with open(cis_p, "rb") as fh:
                obj = pickle.load(fh)
        except Exception as e:
            die(f"{cis_p} does not unpickle ({type(e).__name__}: {e}) -- a TRUNCATED chain")
        names = [str(x) for x in obj.cell_names]
        del obj
        log(f"analysis set: {len(names):,} cell ids from {cis_p.name}")
    else:
        die(f"neither {ids_p} nor {cis_p} exists. The analysis set is defined by the chromatin "
            f"chain, which run_bmmc.sh phase [2] builds before phase [3] runs this script. "
            f"Refusing to invent cell ids: prepare_GEX_ACC joins the modalities BY NAME, and a "
            f"name that differs by one separator matches zero cells.")
    if len(set(names)) != len(names):
        die("the analysis set contains duplicate cell ids")
    base = [n.split(sep)[0] for n in names]
    if len(set(base)) != len(base):
        die(f"splitting the analysis set on cell_id_sep {sep!r} does not give unique barcodes")
    return names, base

def resolve_column(df: pd.DataFrame, candidates, what: str, path: Path) -> str:
    found = [c for c in candidates if c in df.columns]
    if len(found) == 1:
        return found[0]
    if not found:
        die(f"{path.name} has no {what} column. Looked for {candidates}; the file's columns are "
            f"{list(df.columns)}.")
    die(f"{path.name} has MORE THAN ONE candidate {what} column ({found}) -- ambiguous, so this "
        f"aborts rather than picking one.")

def main():
    ap = argparse.ArgumentParser(description="Build one BMMC arm's paired RNA object")
    ap.add_argument("--stage", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--h5ad", required=True,
                    help="raw GSE194122 multiome h5ad; used only if out/bmmc_rna_counts.h5ad "
                         "is absent, always recorded as provenance")
    ap.add_argument("--metadata", required=True, help="uploads/b00_cell_metadata.csv")
    ap.add_argument("--assignment", default=None,
                    help="per-arm assignment table. Omitted for `truth` (barcode identity).")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_config()
    STAGE, ARM = a.stage, a.arm
    OUTP = cfg["P"]["out"]
    OUTP.mkdir(parents=True, exist_ok=True)
    SEP = str(cfg["cell_id_sep"])
    GROUP_BY = cfg["pseudobulk"]["group_by"]
    DS = cfg["dataset"]

    if ARM == "pseudomultiome" or ARM.startswith("pseudomultiome_k"):
        die(f"arm {ARM!r} is a METACELL arm: SCENIC+ builds the metacells itself inside "
            f"prepare_GEX_ACC, so it has no per-cell assignment and no paired RNA object. "
            f"run_bmmc.sh skips it in phase [3]; this script must never be called for it.")
    if int(a.seed) != int(cfg["seed"]):
        log(f"WARNING: --seed {a.seed} differs from config seed {cfg['seed']}. The seed only "
            f"selects which assignment table run_bmmc.sh passes; it is recorded, not applied.")

    out_h5ad = OUTP / f"paired_rna_{STAGE}.{ARM}.h5ad"
    out_rep = OUTP / f"b_03_arm_report_{STAGE}.{ARM}.json"
    if out_h5ad.is_file() and out_h5ad.stat().st_size > 0 and not a.force:
        log(f"SKIP: {out_h5ad.name} already exists ({out_h5ad.stat().st_size / 1e6:.0f} MB)")
        return
    stale = partof(out_h5ad)
    if stale.exists():
        log(f"removing stale sidecar {stale.name} (a killed write leaves no valid output)")
        stale.unlink()

    log(f"stage {STAGE} | arm {ARM} | seed {a.seed} | config {cfg['CFG_PATH'].name}")

    # ---- 1) the analysis set, from the chain ------------------------------------------------
    cis_names, cis_base = analysis_set(OUTP, STAGE, SEP)
    if len(cis_names) != int(DS["n_cells"]):
        die(f"the analysis set holds {len(cis_names):,} cells; config dataset.n_cells is "
            f"{int(DS['n_cells']):,}. Every arm must cover exactly the cells the chain was "
            f"built on.")
    base_pos = {b: i for i, b in enumerate(cis_base)}

    # ---- 2) the per-cell manifest -----------------------------------------------------------
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
    if set(cis_base) != set(meta.index):
        only_chain = set(cis_base) - set(meta.index)
        only_meta = set(meta.index) - set(cis_base)
        die(f"CELL-ID MISMATCH between the chain and {meta_p.name}: {len(only_chain)} chain-only "
            f"(e.g. {sorted(only_chain)[:3]}), {len(only_meta)} manifest-only "
            f"(e.g. {sorted(only_meta)[:3]})")
    log(f"manifest covers the analysis set exactly ({len(meta):,} cells, "
        f"{meta['cell_type'].nunique()} cell types)")

    # ---- 3) the RNA matrix ------------------------------------------------------------------
    presplit = OUTP / "bmmc_rna_counts.h5ad"
    src_h5ad = Path(a.h5ad)
    if presplit.is_file() and presplit.stat().st_size > 0:
        RNA, rna_ids, gene_names = load_rna_presplit(presplit)
        source = str(presplit)
    else:
        if not src_h5ad.is_file() or src_h5ad.stat().st_size == 0:
            die(f"neither {presplit} nor --h5ad {src_h5ad} is present and non-empty")
        RNA, rna_ids, gene_names = load_rna_from_multiome(src_h5ad)
        source = str(src_h5ad)
    rna_ids = np.array([str(x) for x in rna_ids], dtype=object)
    log(f"RNA {RNA.shape[0]:,} cells x {RNA.shape[1]:,} genes | nnz {RNA.nnz:,}")
    assert_raw_counts(RNA, "RNA counts")
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
    rna_pos = {b: i for i, b in enumerate(rna_ids)}

    # ---- 4) the assignment: ATAC cell -> RNA cell -------------------------------------------
    # A function over the analysis set: total (every ATAC cell assigned), single-valued (no ATAC
    # cell twice) and closed (no ATAC cell outside the set). Anything else is a bug in the table.
    if a.assignment:
        asg_p = Path(a.assignment)
        if not asg_p.is_file() or asg_p.stat().st_size == 0:
            die(f"--assignment {asg_p} is missing or empty")
        asg = pd.read_csv(asg_p, dtype=str)
        acol = resolve_column(asg, ATAC_COLS, "ATAC cell id", asg_p)
        rcol = resolve_column(asg, RNA_COLS, "assigned RNA cell id", asg_p)
        log(f"assignment {len(asg):,} rows from {asg_p.name} | ATAC key `{acol}` -> RNA `{rcol}`")
        atac_key = asg[acol].astype(str).values
        rna_key = asg[rcol].astype(str).values
        if pd.isna(asg[rcol]).any():
            die(f"{asg_p.name} has {int(pd.isna(asg[rcol]).sum())} empty `{rcol}` value(s) -- "
                f"BMMC is a true multiome, so no arm can legitimately leave a cell unpaired")

        dup = pd.Index(atac_key)[pd.Index(atac_key).duplicated()].unique().tolist()
        if dup:
            die(f"{asg_p.name} assigns {len(dup)} ATAC cell(s) more than once (e.g. {dup[:5]}) "
                f"-- an assignment must be a FUNCTION over the analysis set")
        aset = set(atac_key)
        missing = [b for b in cis_base if b not in aset]
        extra = sorted(aset - set(cis_base))
        if missing or extra:
            hint = ""
            if extra and set(extra) <= set(meta["fragment_barcode"].astype(str)):
                hint = (" Those values are FRAGMENT barcodes, which are NOT unique across the 13 "
                        "batches. Key on the h5ad barcode (`AAACAGCCAATTAAGG-1-s1d1`) instead.")
            die(f"ASSIGNMENT IS NOT A FUNCTION OVER THE ANALYSIS SET -- refusing to build.\n"
                f"  in the analysis set but not assigned : {len(missing)} "
                f"(e.g. {missing[:3]})\n"
                f"  assigned but not in the analysis set : {len(extra)} (e.g. {extra[:3]})\n"
                f"  BMMC is a true multiome: n_ATAC == n_RNA identically and a bijection is "
                f"always feasible, so a drop here is a bug in the table, not a property of the "
                f"data.{hint}")
        unknown = sorted(set(rna_key) - set(rna_ids))
        if unknown:
            die(f"{len(unknown)} assigned RNA barcode(s) are absent from the RNA matrix, e.g. "
                f"{unknown[:5]}")
        # align the table to cisTopic order: row i of the output is ATAC cell i's RNA profile
        order = pd.Series(np.arange(len(atac_key)), index=pd.Index(atac_key)).loc[cis_base].values
        assigned = rna_key[order]
        assignment_src = str(asg_p)
    else:
        # `truth` is barcode identity -- the one pairing only a real multiome gives, and the
        # yardstick for every other arm. run_bmmc.sh exempts it from --assignment by design.
        if ARM != "truth":
            log(f"WARNING: no --assignment for arm {ARM!r}; falling back to barcode identity, "
                f"which is the `truth` pairing. If that is not what {ARM!r} means, this arm was "
                f"dispatched without its table.")
        assigned = np.array(cis_base, dtype=object)
        assignment_src = "barcode identity (no assignment table by design)"
        log(f"arm {ARM}: barcode-identity pairing over {len(assigned):,} cells")

    rows = np.array([rna_pos[b] for b in assigned], dtype=np.int64)
    if len(rows) != len(cis_names):
        die(f"built {len(rows):,} assignments for {len(cis_names):,} ATAC cells")

    # ---- 5) permute, label, write -----------------------------------------------------------
    import anndata as ad

    sub = RNA[rows]                       # duplicate rows intentional -- reuse is the phenomenon
    del RNA
    if sub.shape[0] != len(cis_names):
        die(f"permuted matrix has {sub.shape[0]:,} rows, expected {len(cis_names):,}")

    atac_meta = meta.loc[cis_base]
    obs = pd.DataFrame(index=pd.Index(cis_names, name=None))
    for c in ("fragment_barcode", "batch", "donor", "site", "cell_type"):
        obs[c] = atac_meta[c].values
    obs["atac_barcode"] = np.array(cis_base, dtype=object)
    # the grouping key run_grn_bmmc.sh passes as --key_to_group_by and pseudobulk.group_by
    # declares. It describes the ATAC cell, so it agrees with the cisTopic object's cell_data.
    obs[GROUP_BY] = atac_meta["cell_type"].astype(str).values
    obs["assigned_rna_barcode"] = assigned
    obs["assigned_rna_cell_type"] = meta.loc[assigned, "cell_type"].values
    obs["assigned_rna_donor"] = meta.loc[assigned, "donor"].values
    obs["arm"] = ARM
    if obs[GROUP_BY].isna().any():
        die(f"`{GROUP_BY}` has {int(obs[GROUP_BY].isna().sum())} missing values")

    var = pd.DataFrame(index=pd.Index([str(g) for g in gene_names], name=None))
    adata = ad.AnnData(X=sub.tocsr(), obs=obs, var=var)
    adata.uns["bmmc_arm"] = {
        "stage": STAGE, "arm": ARM, "seed": int(a.seed),
        "assignment": assignment_src, "rna_source": source,
        "cell_id_sep": SEP, "group_by": GROUP_BY,
    }

    tmp = partof(out_h5ad)
    try:
        adata.write_h5ad(tmp)
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise
    if not tmp.is_file() or tmp.stat().st_size == 0:
        die("write_h5ad returned but produced no output")
    os.replace(tmp, out_h5ad)

    n_distinct = len(set(assigned.tolist()))
    n_self = int(sum(1 for a_, r_ in zip(cis_base, assigned) if a_ == r_))
    same_ct = float(np.mean(obs[GROUP_BY].values == obs["assigned_rna_cell_type"].values))
    same_donor = float(np.mean(obs["donor"].values == obs["assigned_rna_donor"].values))
    log(f"{ARM}: {adata.n_obs:,} x {adata.n_vars:,} | distinct RNA cells {n_distinct:,} "
        f"| reuse {adata.n_obs / n_distinct:.2f}x | self-pairs {n_self:,} "
        f"({n_self / adata.n_obs:.1%}) | same cell type {same_ct:.1%} "
        f"| same donor {same_donor:.1%}")
    log(f"committed {out_h5ad.name} ({out_h5ad.stat().st_size / 1e6:.0f} MB)")

    rep = {
        "stage": STAGE, "arm": ARM, "seed": int(a.seed),
        "assignment": assignment_src, "rna_source": source,
        "provenance_h5ad": str(src_h5ad), "metadata": str(meta_p),
        "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
        "n_distinct_rna_cells": n_distinct,
        "reuse_factor": round(adata.n_obs / n_distinct, 4),
        "n_self_pairs": n_self,
        "frac_same_cell_type": round(same_ct, 6),
        "frac_same_donor": round(same_donor, 6),
        "cell_id_example": cis_names[0],
        "cell_id_sha256_16": digest_strings(cis_names),
        "group_by": GROUP_BY,
    }
    commit_bytes(out_rep, (json.dumps(rep, indent=1) + "\n").encode("utf-8"))
    log("b_03 COMPLETE")

if __name__ == "__main__":
    main()
