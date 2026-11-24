#!/usr/bin/env python
"""Preflight assertion for the metacell (non-multiome) arms."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

GROUP_BY = "Classified_Celltype"   # hardcoded in run_grn.sh's --key_to_group_by; keep in step

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gex", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--cell-data", default=None,
                    help="ATAC cell_data csv; defaults to $AMLGRN_ROOT/work/cell_data_<STAGE>.csv")
    a = ap.parse_args()

    gex_fp = Path(a.gex)
    if not gex_fp.is_file():
        print(f"[p14] MISSING {gex_fp}", file=sys.stderr)
        return 2

    A = ad.read_h5ad(gex_fp)
    if GROUP_BY not in A.obs.columns:
        print(f"[p14] {gex_fp.name} has no '{GROUP_BY}' column; columns = {list(A.obs.columns)}",
              file=sys.stderr)
        return 2
    rna_n = A.obs[GROUP_BY].astype(str).value_counts()

    root = Path(os.environ.get("AMLGRN_ROOT", "."))
    cd_fp = Path(a.cell_data) if a.cell_data else root / "work" / f"cell_data_{a.stage}.csv"
    atac_n = None
    if cd_fp.is_file():
        cd = pd.read_csv(cd_fp, index_col=0, dtype=str)
        if GROUP_BY in cd.columns:
            atac_n = cd[GROUP_BY].astype(str).value_counts()
        else:
            print(f"[p14] WARNING: {cd_fp.name} has no '{GROUP_BY}'; ATAC side NOT checked",
                  file=sys.stderr)
    else:
        print(f"[p14] WARNING: {cd_fp} not found; ATAC side NOT checked", file=sys.stderr)

    labels = sorted(set(rna_n.index) if atac_n is None
                    else set(rna_n.index) & set(atac_n.index))
    print(f"[p14] {a.stage}: k={a.k}  GEX {A.n_obs:,} cells x {A.n_vars:,} genes  "
          f"| {len(labels)} shared cell types"
          f"{' (RNA only -- ATAC unchecked)' if atac_n is None else ''}")

    bad = []
    for ct in labels:
        r = int(rna_n.get(ct, 0))
        t = int(atac_n.get(ct, 0)) if atac_n is not None else None
        worst = r if t is None else min(r, t)
        flag = "  <-- BELOW FLOOR" if worst < a.k else ""
        print(f"[p14]   {ct:34s} RNA={r:>6,}"
              + (f"  ATAC={t:>6,}" if t is not None else "") + f"  min={worst:>6,}{flag}")
        if worst < a.k:
            bad.append((ct, r, t))

    if bad:
        print(f"\n[p14] REFUSING to run at k={a.k}: {len(bad)} cell type(s) hold fewer than {a.k} "
              f"cells and would raise 'Sample larger than population' inside SCENIC+:",
              file=sys.stderr)
        for ct, r, t in bad:
            print(f"[p14]   {ct}: RNA={r}" + (f", ATAC={t}" if t is not None else ""),
                  file=sys.stderr)
        print("[p14] Rebuild the metacell GEX with a floor of at least this k, or lower k.",
              file=sys.stderr)
        return 3

    print(f"[p14] OK: every shared cell type clears the {a.k}-cell floor.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
