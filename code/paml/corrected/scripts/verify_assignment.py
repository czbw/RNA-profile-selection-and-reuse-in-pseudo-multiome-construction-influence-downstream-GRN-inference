#!/usr/bin/env python
"""Verify the repaired assignment independently of the script that produced it."""
from __future__ import annotations

import argparse
import glob
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STAGES = ["Diagnosis", "Relapse"]
EXCLUDE_UNITS = ["AML7_DX"]
BUILDER_RNA_POOL = 220_149   # must equal p_a01.BUILDER_RNA_POOL and the builder's log

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    a = ap.parse_args()
    root = Path(a.root)
    outdir = root / "paml/corrected" / "out"
    fails = []

    # RNA deposit, re-read from source
    files = sorted(glob.glob(str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")))
    R = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in files], ignore_index=True)
    R["_pid"] = R["Patient_ID"].str.strip(); R["_tp"] = R["Patient_Sample"].str.strip()
    R["_ct"] = R["Classified_Celltype"].str.strip(); R["_unit"] = R["Library_ID"].str.strip()
    R = R[R["_tp"].isin(STAGES)].copy()
    # The pool definition must match p_a01's exactly; asserting the builder's pool size keeps
    # assignment, verification and the downstream builder on one definition.
    R = R[R["Malignant"].str.strip().str.lower() == "malignant"].copy()
    if len(R) != BUILDER_RNA_POOL:
        raise SystemExit(f"[a02] RNA pool {len(R):,} != builder's {BUILDER_RNA_POOL:,}")
    R = R[~R["_unit"].isin(EXCLUDE_UNITS)].copy()
    R["_key"] = R["_unit"] + "||" + R["Cell_Barcode"].str.strip()
    rna_keys = set(R["_key"])
    meta = R.set_index("_key")[["_pid", "_tp", "_ct"]].to_dict("index")
    pool_n = R.groupby(["_pid", "_tp", "_ct"]).size()
    print(f"[a02] RNA deposit: {len(R):,} cells, {len(rna_keys):,} unique (library, barcode) keys")

    for stage in STAGES:
        rep = outdir / f"cell_data_{stage}.repaired.csv"
        orig = root / "paml/malignant" / "out" / f"cell_data_{stage}.csv"
        if not rep.is_file():
            fails.append(f"{stage}: MISSING {rep.name}")
            continue
        D = pd.read_csv(rep, index_col=0, dtype=str)
        O = pd.read_csv(orig, index_col=0, dtype=str)
        print(f"\n[a02] --- {stage}: {len(D):,} repaired rows (original {len(O):,}) ---")

        D["_pid"] = D["Patient_ID"].str.strip(); D["_tp"] = D["timepoint"].str.strip()
        D["_ct"] = D["Classified_Celltype"].str.strip()
        D["_key"] = (D["Linked_scRNA_Sample"].str.strip() + "||"
                     + D["Linked_scRNA_barcode"].str.strip())

        # (1) complete, single-valued
        if D["Linked_scRNA_barcode"].isna().any():
            fails.append(f"{stage}: unassigned cells present")
        if D.index.duplicated().any():
            fails.append(f"{stage}: duplicate ATAC cell ids")

        # (3) the join key must resolve -- this is what build_paired_rna.py does
        unresolved = ~D["_key"].isin(rna_keys)
        if unresolved.any():
            fails.append(f"{stage}: {int(unresolved.sum()):,} donor keys absent from the RNA deposit "
                         f"-- build_paired_rna.py would silently drop them")
        else:
            print(f"[a02]   join keys resolve: {len(D):,}/{len(D):,} OK")

        # (2) no crossing, re-derived from the RNA deposit rather than from p_a01's maps
        ok = ~unresolved
        got = pd.DataFrame([meta[k] for k in D.loc[ok, "_key"]], index=D.index[ok])
        for field, label in (("_pid", "patient"), ("_tp", "timepoint"), ("_ct", "cell type")):
            bad = int((got[field].to_numpy() != D.loc[ok, field].to_numpy()).sum())
            if bad:
                fails.append(f"{stage}: {bad:,} assignments CROSS {label}")
            else:
                print(f"[a02]   no {label} crossing: {int(ok.sum()):,} checked")

        # (5) per-stratum reuse never exceeds the independently computed floor
        worst = []
        for key, grp in D.groupby(["_pid", "_tp", "_ct"]):
            n_rna = int(pool_n.get(key, 0))
            if n_rna == 0:
                fails.append(f"{stage}: stratum {key} kept but has no donors")
                continue
            floor = math.ceil(len(grp) / n_rna)
            used = grp["Linked_scRNA_barcode"].value_counts().max()
            if used > floor:
                worst.append((key, int(used), floor))
        if worst:
            fails.append(f"{stage}: {len(worst)} strata exceed their reuse floor, e.g. {worst[:3]}")
        else:
            print(f"[a02]   every stratum at or below its reuse floor")

        # (4) dropped cells are exactly the no-donor ones
        dropped = O.index.difference(D.index)
        Od = O.loc[dropped].copy()
        Od["_pid"] = Od["Patient_ID"].str.strip(); Od["_tp"] = Od["timepoint"].str.strip()
        Od["_ct"] = Od["Classified_Celltype"].str.strip()
        wrongly = 0
        for key, grp in Od.groupby(["_pid", "_tp", "_ct"]):
            if int(pool_n.get(key, 0)) > 0 and key[0] not in ("",):
                wrongly += len(grp)
        if wrongly:
            fails.append(f"{stage}: {wrongly:,} cells dropped despite having donors available")
        else:
            print(f"[a02]   {len(dropped):,} dropped ({len(dropped)/len(O):.2%}), all no-donor")

        reuse = len(D) / D["Linked_scRNA_barcode"].nunique()
        print(f"[a02]   achieved reuse {reuse:.4f}x on {D['Linked_scRNA_barcode'].nunique():,} donors")

    print()
    if fails:
        print("[a02] ===== VERIFICATION FAILED =====", file=sys.stderr)
        for f in fails:
            print(f"[a02]   {f}", file=sys.stderr)
        return 1
    print("[a02] ===== ALL CHECKS PASSED — safe to upload =====")
    return 0

if __name__ == "__main__":
    sys.exit(main())
