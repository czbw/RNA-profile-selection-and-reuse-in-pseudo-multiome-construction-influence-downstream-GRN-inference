#!/usr/bin/env python
"""Assign eligible RNA profiles within patient, stage and cell type with balanced reuse.

The eligible malignant RNA pool must match the downstream RNA builder. Empty strata
are excluded. Link confidence is a schema placeholder, not matching accuracy."""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXCLUDE_UNITS = ["AML7_DX"]
STAGES = ["Diagnosis", "Relapse"]

BUILDER_RNA_POOL = 220_149   # RNA pool size reported by build_paired_rna.py

def load_rna(root: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")))
    if not files:
        raise SystemExit("no RNA metadata found")
    R = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in files], ignore_index=True)
    R["_pid"] = R["Patient_ID"].str.strip()
    R["_tp"] = R["Patient_Sample"].str.strip()
    R["_ct"] = R["Classified_Celltype"].str.strip()
    R["_unit"] = R["Library_ID"].str.strip()
    R = R[R["_tp"].isin(STAGES)].copy()
    mal = R["Malignant"].str.strip().str.lower() == "malignant"
    R = R[mal].copy()
    n_before_excl = len(R)
    if n_before_excl != BUILDER_RNA_POOL:
        raise SystemExit(
            f"[a01] RNA pool is {n_before_excl:,} but build_paired_rna.py supplies "
            f"{BUILDER_RNA_POOL:,}. These MUST match or donors will silently fail to join. "
            f"Investigate before proceeding -- do not adjust this constant to make it pass.")
    R = R[~R["_unit"].isin(EXCLUDE_UNITS)].copy()
    print(f"[a01] RNA pool: {len(R):,} malignant cells (Dx/Rel; matches the builder's "
          f"{BUILDER_RNA_POOL:,} before the AML7_DX exclusion)")
    return R

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument("--seed", type=int, default=555)
    a = ap.parse_args()
    root = Path(a.root)
    outdir = root / "paml/corrected" / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    R = load_rna(root)
    # donor pools keyed by the full stratum (patient, timepoint, cell type), so no assignment can cross strata
    pools = {k: v for k, v in R.groupby(["_pid", "_tp", "_ct"], sort=True)}

    summary, per_stratum = {}, []
    for stage in STAGES:
        fp = root / "paml/malignant" / "out" / f"cell_data_{stage}.csv"
        A = pd.read_csv(fp, index_col=0, dtype=str)
        A["_pid"] = A["Patient_ID"].str.strip()
        A["_tp"] = A["timepoint"].str.strip()
        A["_ct"] = A["Classified_Celltype"].str.strip()
        n0 = len(A)
        A = A[~A["sample_id"].str.strip().isin(EXCLUDE_UNITS)].copy()

        new_bc = pd.Series(index=A.index, dtype=object)
        new_lib = pd.Series(index=A.index, dtype=object)
        dropped_idx = []

        for key, grp in A.groupby(["_pid", "_tp", "_ct"], sort=True):
            pool = pools.get(key)
            n_atac = len(grp)
            if pool is None or len(pool) == 0:
                dropped_idx.extend(grp.index.tolist())
                per_stratum.append({"patient": key[0], "timepoint": key[1], "cell_type": key[2],
                                    "n_atac": n_atac, "n_rna": 0, "achieved_max_reuse": None,
                                    "status": "DROPPED_no_donor"})
                continue
            n_rna = len(pool)
            # shuffle both sides, then assign cyclically -> max use = ceil(n_atac/n_rna), the floor
            atac_order = rng.permutation(n_atac)
            rna_order = rng.permutation(n_rna)
            pick = rna_order[np.arange(n_atac) % n_rna]
            idx = grp.index.to_numpy()[atac_order]
            new_bc.loc[idx] = pool["Cell_Barcode"].to_numpy()[pick]
            new_lib.loc[idx] = pool["_unit"].to_numpy()[pick]
            per_stratum.append({
                "patient": key[0], "timepoint": key[1], "cell_type": key[2],
                "n_atac": n_atac, "n_rna": n_rna,
                "achieved_max_reuse": int(np.ceil(n_atac / n_rna)), "status": "OK"})

        keep = ~pd.Index(A.index).isin(dropped_idx)
        out = A.loc[keep].copy()
        out["Linked_scRNA_barcode"] = new_bc.loc[out.index]
        out["Linked_scRNA_Sample"] = new_lib.loc[out.index]
        out["Linked_scRNA_confidence_score"] = 1.0   # placeholder; original was a measured no-op

        assert out["Linked_scRNA_barcode"].notna().all(), "unassigned cell survived"
        # the constraint, asserted rather than trusted
        tp_of_unit = dict(zip(R["_unit"], R["_tp"]))
        pid_of_unit = dict(zip(R["_unit"], R["_pid"]))
        lt = out["Linked_scRNA_Sample"].map(tp_of_unit)
        lp = out["Linked_scRNA_Sample"].map(pid_of_unit)
        assert (lt == out["_tp"]).all(), "TIMEPOINT CROSSED — constraint violated"
        assert (lp == out["_pid"]).all(), "PATIENT CROSSED — constraint violated"

        reuse = len(out) / out["Linked_scRNA_barcode"].nunique()
        donors = out["Linked_scRNA_barcode"].nunique()
        old = A["Linked_scRNA_barcode"]
        self_hit = float((out["Linked_scRNA_barcode"].to_numpy()
                          == A.loc[out.index, "Linked_scRNA_barcode"].to_numpy()).mean())

        dest = outdir / f"cell_data_{stage}.repaired.csv"
        out.drop(columns=["_pid", "_tp", "_ct"]).to_csv(dest)
        summary[stage] = {
            "n_in": n0, "n_out": len(out), "n_dropped_no_donor": len(dropped_idx),
            "pct_dropped": len(dropped_idx) / n0,
            "achieved_reuse": reuse, "distinct_donors": int(donors),
            "old_reuse": len(A) / old.nunique(), "old_distinct_donors": int(old.nunique()),
            "agreement_with_old_link": self_hit, "file": dest.name}
        print(f"[a01] {stage}: {len(out):,}/{n0:,} cells kept "
              f"({len(dropped_idx):,} dropped, {len(dropped_idx)/n0:.2%}) | "
              f"reuse {reuse:.3f}x (was {len(A)/old.nunique():.2f}x) | "
              f"donors {donors:,} (was {old.nunique():,}) | "
              f"agreement with old link {self_hit:.2%}")

    pd.DataFrame(per_stratum).to_csv(outdir / "a01_per_stratum.csv", index=False)
    (outdir / "a01_assignment_report.json").write_text(json.dumps(
        {"seed": a.seed, "scope": "malignant Diagnosis + Relapse",
         "excluded_units": EXCLUDE_UNITS, "malignant_rna_filter": False,
         "unassignable_policy": "DROP (no cell-type relaxation)",
         "constraints": ["within-patient", "within-timepoint", "within-cell-type"],
         "stages": summary}, indent=2))
    print(f"[a01] wrote {outdir/'a01_assignment_report.json'} and a01_per_stratum.csv")
    return 0

if __name__ == "__main__":
    sys.exit(main())
