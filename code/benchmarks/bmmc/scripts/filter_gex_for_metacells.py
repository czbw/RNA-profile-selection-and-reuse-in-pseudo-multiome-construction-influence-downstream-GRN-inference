#!/usr/bin/env python
"""Filter the native GEX object to the cell types that clear the metacell floor in both modalities.

The floor is max(min_cells_per_celltype, k): a cell type below it makes SCENIC+ metacell sampling
raise "Sample larger than population". For every k <= min_cells_per_celltype the filtered object
is identical, so one file serves all metacell arms (checked with a content digest).
Adapted from the AML bundle's filter script.
"""
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
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

def log(m):
    print(f"[b_11] {m}", flush=True)

def die(msg):
    print(f"[b_11] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def load_config():
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

def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)

def _index_of(g) -> np.ndarray:
    key = g.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return _decode(g[str(key)][:])

def _obs_col(f, col) -> np.ndarray:
    """obs column as strings, whether it is stored categorical (group), legacy-categorical, or
    a plain string dataset."""
    path = f"obs/{col}"
    node = f[path]
    if isinstance(node, h5py.Group):
        cats = _decode(node["categories"][:])
        return np.array([cats[c] for c in node["codes"][:]], dtype=object)
    legacy = f"obs/__categories/{col}"
    if legacy in f:
        cats = _decode(f[legacy][:])
        return np.array([cats[c] for c in node[:]], dtype=object)
    return _decode(node[:])

def content_digest(obs_ids, var_ids, kept_labels, X) -> str:
    """A digest of what the object is, not of how h5py happened to serialise it."""
    h = hashlib.sha256()
    for tag, seq in (("obs", obs_ids), ("var", var_ids), ("kept", sorted(kept_labels))):
        h.update(tag.encode())
        for x in seq:
            h.update(str(x).encode())
            h.update(b"\x00")
    h.update(b"X")
    h.update(str(X.shape).encode())
    for name in ("data", "indices", "indptr"):
        h.update(np.ascontiguousarray(getattr(X, name)).tobytes())
    return h.hexdigest()

def atac_counts(a, cfg, group_by, gex_barcodes):
    """(counts Series indexed by label, source description). Uses the first available ATAC-side
    source, in order: the ATAC counts h5ad, a cell-data table, cisTopic cell names + metadata."""
    OUTP = cfg["P"]["out"]
    if a.atac_h5ad:
        p = Path(a.atac_h5ad)
    else:
        p = OUTP / "bmmc_atac_counts.h5ad"
    if p.is_file() and p.stat().st_size > 0:
        with h5py.File(p, "r") as f:
            ids = _index_of(f["obs"])
            cols = [c for c in (group_by, "cell_type") if f"obs/{c}" in f]
            if not cols:
                die(f"{p.name} obs has neither `{group_by}` nor `cell_type`; it has "
                    f"{list(f['obs'].keys())}")
            lab = _obs_col(f, cols[0])
        if len(set(ids)) != len(ids):
            die(f"{p.name} has duplicate cell ids")
        if set(ids) != set(gex_barcodes):
            die(f"ATAC/GEX CELL-SET MISMATCH: {p.name} covers {len(set(ids) - set(gex_barcodes))} "
                f"cells the GEX object does not, and misses "
                f"{len(set(gex_barcodes) - set(ids))} that it has. b_01 asserts both modalities "
                f"are the same cells keyed on the h5ad barcode; a disjoint ATAC side would make "
                f"every count below meaningless.")
        log(f"ATAC side: {p.name} obs/{cols[0]} ({len(ids):,} cells)")
        return pd.Series(lab).astype(str).value_counts(), f"{p} (obs/{cols[0]})"

    if a.cell_data:
        p = Path(a.cell_data)
        if not p.is_file() or p.stat().st_size == 0:
            die(f"--cell-data {p} is missing or empty")
        cd = pd.read_csv(p, index_col=0, dtype=str)
        col = next((c for c in (group_by, "cell_type") if c in cd.columns), None)
        if col is None:
            die(f"{p.name} has neither `{group_by}` nor `cell_type`; it has {list(cd.columns)}")
        log(f"ATAC side: {p.name} column {col} ({len(cd):,} rows)")
        return cd[col].astype(str).value_counts(), f"{p} ({col})"

    names_p = Path(a.cell_names) if a.cell_names else OUTP / f"cistopic_cell_names_{a.stage}.txt"
    meta_p = Path(a.metadata) if a.metadata else None
    if names_p.is_file() and names_p.stat().st_size > 0 and meta_p and meta_p.is_file():
        sep = str(cfg["cell_id_sep"])
        names = [ln.strip() for ln in names_p.read_text(encoding="utf-8").splitlines()
                 if ln.strip()]
        base = [n.split(sep)[0] for n in names]
        meta = pd.read_csv(meta_p, dtype=str)
        if "h5ad_barcode" not in meta.columns or "cell_type" not in meta.columns:
            die(f"{meta_p.name} needs h5ad_barcode and cell_type; it has {list(meta.columns)}")
        meta = meta.set_index("h5ad_barcode")
        if set(base) != set(gex_barcodes):
            die(f"ATAC/GEX CELL-SET MISMATCH between {names_p.name} and the GEX object's "
                f"atac_barcode column")
        missing = [b for b in base if b not in meta.index]
        if missing:
            die(f"{len(missing)} chain cell id(s) have no row in {meta_p.name} "
                f"(e.g. {missing[:3]})")
        log(f"ATAC side: {names_p.name} joined to {meta_p.name} ({len(base):,} cells)")
        return (meta.loc[base, "cell_type"].astype(str).value_counts(),
                f"{names_p} joined to {meta_p} (cell_type)")

    die("no ATAC-side cell-type source. Give --atac-h5ad, or --cell-data, or "
        "--cell-names together with --metadata. SCENIC+ samples k cells from BOTH modalities, "
        "so filtering on the RNA side alone would not measure the thing this script claims.")

def main() -> int:
    ap = argparse.ArgumentParser(description="Apply the metacell cell-type floor to the native "
                                             "BMMC RNA object")
    ap.add_argument("--stage", default=os.environ.get("AMLGRN_STAGE", "bmmc"))
    ap.add_argument("--k", type=int, default=None,
                    help="nr_cells_per_metacells; default = config "
                         "pairing_params.pseudomultiome.nr_cells_per_metacell")
    ap.add_argument("--gex", default=None, help="default out/native_rna_<STAGE>.h5ad")
    ap.add_argument("--atac-h5ad", default=None, help="default out/bmmc_atac_counts.h5ad")
    ap.add_argument("--cell-data", default=None)
    ap.add_argument("--cell-names", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing metacell object even if its content digest "
                         "differs -- see the k-independence block in the module docstring")
    a = ap.parse_args()

    cfg = load_config()
    STAGE = a.stage
    OUTP = cfg["P"]["out"]
    pm = cfg["pairing_params"]["pseudomultiome"]
    GROUP_BY = str(pm.get("group_by", cfg["pseudobulk"]["group_by"]))
    K = int(a.k) if a.k is not None else int(pm["nr_cells_per_metacell"])
    if K < 2:
        die(f"metacell size must be >= 2 (got {K}) -- run_grn_bmmc.sh asserts the same")
    FLOOR = max(int(pm.get("min_cells_per_celltype", 25)), K)

    gex_fp = Path(a.gex) if a.gex else OUTP / f"native_rna_{STAGE}.h5ad"
    out_fp = OUTP / f"native_rna_{STAGE}.metacell.h5ad"
    rep_fp = OUTP / f"b_11_metacell_filter_report_{STAGE}_k{K}.json"
    if not gex_fp.is_file() or gex_fp.stat().st_size == 0:
        die(f"missing {gex_fp} -- native_rna.py builds it")

    log(f"stage {STAGE} | k {K} | floor max(min_cells_per_celltype="
        f"{int(pm.get('min_cells_per_celltype', 25))}, k={K}) = {FLOOR} | group_by {GROUP_BY}")

    import anndata as ad
    A = ad.read_h5ad(gex_fp)
    if GROUP_BY not in A.obs.columns:
        die(f"{gex_fp.name} has no `{GROUP_BY}` column; it has {list(A.obs.columns)}")
    if "atac_barcode" not in A.obs.columns:
        die(f"{gex_fp.name} has no `atac_barcode` column -- b_10 writes it, and without it the "
            f"ATAC/GEX cell-set assertion below cannot be made")
    log(f"GEX {A.n_obs:,} cells x {A.n_vars:,} genes from {gex_fp.name}")
    gex_n = A.obs[GROUP_BY].astype(str).value_counts()
    atac_n, atac_src = atac_counts(a, cfg, GROUP_BY, A.obs["atac_barcode"].astype(str).values)

    only_one = sorted(set(gex_n.index) ^ set(atac_n.index))
    for ct in only_one:
        log(f"'{ct}' is in ONE modality only (RNA={int(gex_n.get(ct, 0))}, "
            f"ATAC={int(atac_n.get(ct, 0))}) -- SCENIC+ drops it regardless")

    shared = sorted(set(gex_n.index) & set(atac_n.index))
    keep, drop, table = [], [], []
    for ct in shared:
        r, t = int(gex_n[ct]), int(atac_n[ct])
        worst = min(r, t)
        (keep if worst >= FLOOR else drop).append(ct)
        table.append({"cell_type": ct, "n_rna": r, "n_atac": t, "n_min": worst,
                      "kept": bool(worst >= FLOOR)})
    for row in sorted(table, key=lambda d: d["n_min"]):
        log(f"  {row['cell_type']:24s} RNA={row['n_rna']:>6,}  ATAC={row['n_atac']:>6,}  "
            f"min={row['n_min']:>6,}  {'keep' if row['kept'] else 'DROP (below floor)'}")
    if not keep:
        die(f"no cell type clears the {FLOOR}-cell floor in both modalities")

    if drop:
        log(f"DROPPING {len(drop)} cell type(s) below the {FLOOR}-cell floor: {drop}")
    else:
        log(f"FILTER IS A NO-OP at k={K}: all {len(shared)} shared cell types clear the "
            f"{FLOOR}-cell floor (thinnest is "
            f"'{min(table, key=lambda d: d['n_min'])['cell_type']}' at "
            f"{min(d['n_min'] for d in table):,}). Recorded as a measured result.")

    B = A[A.obs[GROUP_BY].astype(str).isin(keep)].copy()
    B.obs[GROUP_BY] = B.obs[GROUP_BY].astype(str)      # shed unused categorical levels
    n_dropped_cells = int(A.n_obs - B.n_obs)
    log(f"metacell GEX {B.n_obs:,} cells x {B.n_vars:,} genes | kept {len(keep)} cell types, "
        f"dropped {len(drop)} ({n_dropped_cells:,} cells)")

    dig = content_digest(list(B.obs_names), list(B.var_names), keep, B.X.tocsr())
    log(f"content digest {dig[:16]}")

    # ---- k-independence: the same file must serve every metacell arm ------------------------
    prev = None
    prev_rep = sorted(OUTP.glob(f"b_11_metacell_filter_report_{STAGE}_k*.json"))
    for p in prev_rep:
        if p == rep_fp:
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if d.get("content_digest"):
            prev = (p.name, int(d.get("k", -1)), d["content_digest"])
            break
    if prev and prev[2] != dig and not a.force:
        die(f"K-INDEPENDENCE VIOLATED: this run at k={K} produces content digest {dig[:16]}, but "
            f"{prev[0]} recorded {prev[2][:16]} at k={prev[1]}. The floor is "
            f"max(min_cells_per_celltype, k), so every k <= "
            f"{int(pm.get('min_cells_per_celltype', 25))} must give the SAME object and both "
            f"metacell arms must share one file. Two arms differing in cell set as well as in k "
            f"is the confound this ladder exists to exclude. Investigate before overwriting "
            f"(--force overrides).")
    if prev and prev[2] == dig:
        log(f"k-independence CONFIRMED against {prev[0]} (k={prev[1]}): identical content digest")

    if out_fp.is_file() and out_fp.stat().st_size > 0 and prev and prev[2] == dig and not a.force:
        log(f"{out_fp.name} already exists with the same content -- not rewritten")
    else:
        tmp = partof(out_fp)
        try:
            B.write_h5ad(tmp, compression="gzip")
        except BaseException:
            if tmp.exists():
                tmp.unlink()
            raise
        if not tmp.is_file() or tmp.stat().st_size == 0:
            die("write_h5ad returned but produced no output")
        os.replace(tmp, out_fp)
        log(f"committed {out_fp.name} ({out_fp.stat().st_size / 1e6:.0f} MB)")

    rep = {
        "script": "filter_gex_for_metacells.py",
        "ported_from": "paml/corrected/scripts/filter_gex_for_metacells.py",
        "stage": STAGE, "k": K, "floor": FLOOR,
        "min_cells_per_celltype": int(pm.get("min_cells_per_celltype", 25)),
        "group_by": GROUP_BY,
        "gex_input": str(gex_fp),
        "gex_input_sha256": sha256_file(gex_fp),
        "atac_source": atac_src,
        "output": str(out_fp),
        "output_sha256": sha256_file(out_fp),
        "output_bytes": int(out_fp.stat().st_size),
        "content_digest": dig,
        "n_obs_in": int(A.n_obs), "n_obs_out": int(B.n_obs),
        "n_cells_dropped": n_dropped_cells,
        "n_celltypes_shared": len(shared),
        "n_celltypes_kept": len(keep), "n_celltypes_dropped": len(drop),
        "celltypes_dropped": drop,
        "celltypes_in_one_modality_only": only_one,
        "filter_is_noop": bool(not drop and not only_one),
        "per_celltype": sorted(table, key=lambda d: d["n_min"]),
    }
    commit_bytes(rep_fp, (json.dumps(rep, indent=1) + "\n").encode("utf-8"))
    log(f"wrote {rep_fp.name}")
    log("b_11 COMPLETE")
    return 0

if __name__ == "__main__":
    sys.exit(main())
