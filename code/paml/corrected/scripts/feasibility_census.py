#!/usr/bin/env python
"""Feasibility census: how much RNA-donor reuse is forced by per-stratum cell counts."""
from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXCLUDE_UNITS_DEFAULT = ["AML7_DX"]
ATAC_REQ = ["barcode", "sample_id", "timepoint", "Classified_Celltype", "Patient_ID",
            "Linked_scRNA_barcode", "Linked_scRNA_Sample"]
RNA_REQ = ["Cell_Barcode", "Library_ID", "Patient_Sample", "Classified_Celltype",
           "Malignant", "Patient_ID"]

def load_atac(root: Path) -> pd.DataFrame:
    """The ATAC set actually used for SCENIC+ (not the raw GEO metadata)."""
    frames = []
    for tp in ("Diagnosis", "Relapse"):
        fp = root / "paml/malignant" / "out" / f"cell_data_{tp}.csv"
        if not fp.is_file():
            raise SystemExit(f"MISSING {fp}")
        d = pd.read_csv(fp, index_col=0, dtype=str)
        missing = [c for c in ATAC_REQ if c not in d.columns]
        if missing:
            raise SystemExit(f"{fp.name}: missing columns {missing}; found {list(d.columns)}")
        d["_src"] = fp.name
        frames.append(d)
        print(f"[a00] ATAC {fp.name}: {len(d):,} cells")
    return pd.concat(frames, ignore_index=True)

def load_rna(root: Path, malignant_only: bool) -> pd.DataFrame:
    """RNA side from the GEO deposit.

    The malignant-only filter is OFF by default. The ATAC "malignant compartment" contains cells
    labelled CD4.Naive, CD8.Naive, B.Cell, NK, Late.Erythrocyte and pDC, which are called
    non-malignant on the RNA side; requiring malignant RNA would leave those strata with no donors
    because of a label disagreement between modalities. The effect of the filter is printed.
    """
    pat = str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")
    files = sorted(glob.glob(pat))
    if not files:
        raise SystemExit(f"no RNA metadata matched {pat}")
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t", dtype=str, compression="gzip")
        missing = [c for c in RNA_REQ if c not in d.columns]
        if missing:
            raise SystemExit(f"{Path(f).name}: missing columns {missing}")
        frames.append(d)
    R = pd.concat(frames, ignore_index=True)
    mal = R["Malignant"].astype(str).str.strip().str.lower() == "malignant"
    print(f"[a00] RNA: {len(R):,} cells across {len(files)} libraries "
          f"({int(mal.sum()):,} malignant, {mal.mean():.1%})")
    if malignant_only:
        print("[a00]   RNA pool = MALIGNANT ONLY (--malignant-rna-only)")
        return R[mal].copy()
    print("[a00]   RNA pool = ALL cells (default; see docstring for why)")
    return R.copy()

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument("--include-aml7-dx", action="store_true",
                    help="override the declared AML7_DX exclusion (documented 2026-08-16)")
    ap.add_argument("--malignant-rna-only", action="store_true",
                    help="restrict the RNA pool to Malignant cells. OFF by default: it creates "
                         "artificial scarcity in normal-lineage strata (see load_rna docstring).")
    a = ap.parse_args()
    root = Path(a.root)
    excl = [] if a.include_aml7_dx else list(EXCLUDE_UNITS_DEFAULT)

    A = load_atac(root)
    R = load_rna(root, a.malignant_rna_only)

    # normalise the join keys. timepoint on the ATAC side; Patient_Sample on the RNA side.
    A["_tp"] = A["timepoint"].astype(str).str.strip()
    R["_tp"] = R["Patient_Sample"].astype(str).str.strip()
    for df in (A, R):
        df["_pid"] = df["Patient_ID"].astype(str).str.strip()
        df["_ct"] = df["Classified_Celltype"].astype(str).str.strip()
    R["_unit"] = R["Library_ID"].astype(str).str.strip()
    A["_unit"] = A["sample_id"].astype(str).str.strip()

    A = A[A["_tp"].isin(["Diagnosis", "Relapse"])].copy()
    R = R[R["_tp"].isin(["Diagnosis", "Relapse"])].copy()
    print(f"[a00] after Dx/Rel restriction: ATAC {len(A):,} | RNA {len(R):,}")

    if excl:
        nA0, nR0 = len(A), len(R)
        A = A[~A["_unit"].isin(excl)].copy()
        R = R[~R["_unit"].isin(excl)].copy()
        print(f"[a00] DECLARED EXCLUSION {excl}: dropped {nA0-len(A):,} ATAC, {nR0-len(R):,} RNA")

    # ---- per-stratum census ------------------------------------------------
    na = A.groupby(["_pid", "_tp", "_ct"]).size().rename("n_atac")
    nr = R.groupby(["_pid", "_tp", "_ct"]).size().rename("n_rna")
    C = pd.concat([na, nr], axis=1).fillna(0).astype(int)
    C = C[C["n_atac"] > 0].copy()
    C["min_reuse"] = [
        (math.inf if r == 0 else max(1.0, math.ceil(x / r)))
        for x, r in zip(C["n_atac"], C["n_rna"])
    ]
    C["short"] = C["n_rna"] < C["n_atac"]
    C["unassignable"] = C["n_rna"] == 0

    n_strata = len(C)
    n_short = int(C["short"].sum())
    atac_total = int(C["n_atac"].sum())
    atac_short = int(C.loc[C["short"], "n_atac"].sum())
    atac_unassignable = int(C.loc[C["unassignable"], "n_atac"].sum())

    assignable = C[~C["unassignable"]]
    # cell-weighted global floor: total ATAC / total distinct donors actually needed
    donors_needed = int(np.minimum(assignable["n_atac"], assignable["n_rna"]).sum())
    global_floor = float(assignable["n_atac"].sum() / donors_needed) if donors_needed else float("inf")

    print("\n[a00] ================= FEASIBILITY CENSUS =================")
    print(f"[a00] strata                       : {n_strata:,}")
    print(f"[a00] strata with n_RNA < n_ATAC   : {n_short:,}  ({n_short/n_strata:.1%})")
    print(f"[a00] ATAC cells in short strata   : {atac_short:,} / {atac_total:,} "
          f"({atac_short/atac_total:.2%})")
    print(f"[a00] ATAC cells UNASSIGNABLE (n_RNA=0): {atac_unassignable:,} "
          f"({atac_unassignable/atac_total:.2%})")
    print(f"[a00] MINIMUM FORCED GLOBAL REUSE  : {global_floor:.3f}x")

    # ---- what the original study's link did, for contrast ------------------
    link_donors = A["Linked_scRNA_barcode"].astype(str).nunique()
    link_reuse = len(A) / link_donors if link_donors else float("nan")
    rna_pool = len(R)
    print(f"\n[a00] Lambo link, same cells       : {link_reuse:.2f}x reuse")
    print(f"[a00]   distinct RNA donors used   : {link_donors:,} of {rna_pool:,} malignant RNA "
          f"({link_donors/rna_pool:.2%})")
    print(f"[a00]   attributable to SCARCITY   : {global_floor:.3f}x")
    print(f"[a00]   attributable to the RULE   : {link_reuse/global_floor:.2f}x on top "
          f"({1 - (global_floor-1)/(link_reuse-1):.1%} of the excess)")

    # ---- the within-timepoint violation, counted ---------------------------
    same_tp = None
    if "Linked_scRNA_Sample" in A.columns:
        tp_of_unit = dict(zip(R["_unit"], R["_tp"]))
        linked_tp = A["Linked_scRNA_Sample"].astype(str).str.strip().map(tp_of_unit)
        ok = linked_tp.notna()
        same_tp = float((linked_tp[ok] == A.loc[ok, "_tp"]).mean())
        print(f"\n[a00] links staying WITHIN timepoint: {same_tp:.2%} "
              f"(resolvable for {int(ok.sum()):,} of {len(A):,})")
        print("[a00]   -> a hard within-timepoint constraint is the single largest correction")

    worst = C.sort_values("min_reuse", ascending=False).head(12)
    print("\n[a00] worst strata by forced reuse:")
    for (pid, tp, ct), row in worst.iterrows():
        mr = "INF (no RNA)" if not np.isfinite(row["min_reuse"]) else f"{row['min_reuse']:.2f}x"
        print(f"[a00]   {pid:>10} {tp:<10} {ct:<24} ATAC {row['n_atac']:>6,}  "
              f"RNA {row['n_rna']:>6,}  floor {mr}")

    out = Path(a.root) / "paml/corrected" / "out"
    out.mkdir(parents=True, exist_ok=True)
    C.reset_index().rename(columns={"_pid": "patient", "_tp": "timepoint", "_ct": "cell_type"}) \
     .to_csv(out / "a00_stratum_census.csv", index=False)
    rep = {
        "scope": "malignant Diagnosis + Relapse",
        "excluded_units": excl,
        "n_strata": n_strata, "n_strata_short": n_short,
        "atac_total": atac_total, "atac_in_short_strata": atac_short,
        "pct_atac_in_short_strata": atac_short / atac_total,
        "atac_unassignable": atac_unassignable,
        "minimum_forced_global_reuse": global_floor,
        "lambo_link_reuse": link_reuse,
        "lambo_link_distinct_donors": int(link_donors),
        "malignant_rna_pool": int(rna_pool),
        "pct_rna_pool_used_by_link": link_donors / rna_pool,
        "pct_links_within_timepoint": same_tp,
    }
    (out / "a00_feasibility_report.json").write_text(json.dumps(rep, indent=2))
    print(f"\n[a00] wrote {out/'a00_stratum_census.csv'} and a00_feasibility_report.json")
    return 0

if __name__ == "__main__":
    sys.exit(main())
