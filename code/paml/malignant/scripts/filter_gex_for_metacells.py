"""Preflight for Option B (pseudomultiome): drop cell types with too few cells in either modality
for metacell construction, and write the filtered RNA h5ad."""
import os
import sys

os.environ.setdefault("AMLGRN_CONFIG", "config.yaml")

import anndata as ad
import pandas as pd

from common import load_config, log

cfg = load_config()
OUTP = cfg["P"]["out"]
WORK = cfg["P"]["work"]

STAGE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("AMLGRN_STAGE", "")
if not STAGE:
    raise SystemExit("usage: python filter_gex_for_metacells.py <Stage>")

pm = cfg["pairing_params"]["pseudomultiome"]
GROUP_BY = pm["group_by"]
NR_CELLS = int(pm["nr_cells_per_metacell"])          # must match run_grn.sh --nr_cells_per_metacells
FLOOR = max(int(pm.get("min_cells_per_celltype", 25)), NR_CELLS)

gex_fp = OUTP / f"native_rna_{STAGE}.h5ad"           # built by build_native_rna.py
atac_fp = WORK / f"cell_data_{STAGE}.csv"            # built by make_cell_data.py (= cisTopic cells)
out_fp = OUTP / f"native_rna_{STAGE}.metacell.h5ad"
for fp in (gex_fp, atac_fp):
    if not fp.exists():
        raise SystemExit(f"{STAGE}: missing {fp}")

A = ad.read_h5ad(gex_fp)
gex_n = A.obs[GROUP_BY].astype(str).value_counts()
atac_n = pd.read_csv(atac_fp, index_col=0, dtype=str)[GROUP_BY].astype(str).value_counts()

for ct in sorted(set(gex_n.index) ^ set(atac_n.index)):
    log(f"{STAGE}: '{ct}' is in one modality only (RNA={int(gex_n.get(ct, 0))}, "
        f"ATAC={int(atac_n.get(ct, 0))}) — SCENIC+ drops it regardless")

keep, drop = [], []
for ct in sorted(set(gex_n.index) & set(atac_n.index)):
    (keep if min(int(gex_n[ct]), int(atac_n[ct])) >= FLOOR else drop).append(ct)
for ct in drop:
    log(f"{STAGE}: DROP '{ct}' — RNA={int(gex_n[ct])}, ATAC={int(atac_n[ct])} "
        f"(floor {FLOOR}: too thin for {NR_CELLS}-cell metacells)")

if not keep:
    raise SystemExit(f"{STAGE}: no cell type clears the {FLOOR}-cell floor in both modalities.")

A = A[A.obs[GROUP_BY].astype(str).isin(keep)].copy()
A.obs[GROUP_BY] = A.obs[GROUP_BY].astype(str)        # shed unused categorical levels
A.write_h5ad(out_fp)
log(f"{STAGE}: metacell GEX {A.n_obs:,} cells x {A.n_vars:,} genes | kept {len(keep)} cell types, "
    f"dropped {len(drop)} -> {out_fp.name}")
