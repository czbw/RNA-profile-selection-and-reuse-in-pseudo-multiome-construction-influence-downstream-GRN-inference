#!/usr/bin/env python
"""Build the corrected pairing for healthy cells: minimum reuse, within patient, timepoint and cell type.

The confidence score is written as a 1.0 placeholder; it is not a measured matching accuracy.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ATAC_MANIFEST = ("preprocessing", "out")
ALL_STAGES = ["Diagnosis", "Remission", "Relapse"]
DID_STAGES = ["Diagnosis", "Relapse"]
HEALTHY_LABEL = "Normal"
CONF_GATE = 0.5          # build_paired_rna.py's linkage.min_confidence

# healthy RNA cells per stage in the downstream h5ads; asserted against the metadata in load_rna
BUILDER_RNA_POOL = {"Diagnosis": 18_910, "Remission": 78_432, "Relapse": 35_170}

def log(m: str) -> None:
    print(f"[h01] {m}", flush=True)

def atac_unit(sample_id: pd.Series) -> pd.Series:
    """'GSM7498708_AML16_DX' -> 'AML16_DX'. The RNA side carries no GSM prefix."""
    return sample_id.astype(str).str.strip().str.split("_", n=1).str[1]

def load_rna(root: Path, exclude: list[str]) -> pd.DataFrame:
    files = sorted(glob.glob(str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")))
    if not files:
        raise SystemExit("[h01] no RNA metadata found")
    R = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in files], ignore_index=True)
    R["_pid"] = R["Patient_ID"].str.strip()
    R["_tp"] = R["Patient_Sample"].str.strip()          # authoritative; NOT the library token
    R["_ct"] = R["Classified_Celltype"].str.strip()
    R["_unit"] = R["Library_ID"].str.strip()
    R = R[R["Malignant"].astype(str).str.strip() == HEALTHY_LABEL].copy()

    counts = R["_tp"].value_counts().to_dict()
    for tp in ALL_STAGES:
        got, want = int(counts.get(tp, 0)), BUILDER_RNA_POOL[tp]
        if got != want:
            raise SystemExit(
                f"[h01] healthy RNA at {tp} is {got:,} but the consumer's "
                f"data/paml/healthy_{tp.lower()}.h5ad supplies {want:,}. These MUST match or donors "
                f"will silently fail to join -- do not adjust the constant to make it pass.")
    log(f"healthy RNA pool: {len(R):,} cells "
        + ", ".join(f"{t} {counts[t]:,}" for t in ALL_STAGES)
        + "  (locked to the consumer's h5ads)")
    if exclude:
        n0 = len(R)
        R = R[~R["_unit"].isin(exclude)].copy()
        log(f"excluded units {exclude}: dropped {n0-len(R):,} RNA cells")
    return R

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--include-remission", action="store_true")
    ap.add_argument("--include-aml7-dx", action="store_true")
    a = ap.parse_args()
    root = Path(a.root)
    stages = ALL_STAGES if a.include_remission else DID_STAGES
    exclude = [] if a.include_aml7_dx else ["AML7_DX"]
    outdir = root / "paml/healthy_corrected" / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    log(f"scope: healthy {' + '.join(stages)} | seed {a.seed} | excluded {exclude or 'none'}")

    R = load_rna(root, exclude)
    # Donor pools keyed by the full stratum (patient, timepoint, cell type), so no assignment can cross strata.
    pools = {k: v for k, v in R.groupby(["_pid", "_tp", "_ct"], sort=True)}

    summary, per_stratum = {}, []
    for stage in stages:
        fp = root.joinpath(*ATAC_MANIFEST) / f"cell_data_{stage}.csv"
        if not fp.is_file():
            raise SystemExit(f"[h01] MISSING {fp}")
        A = pd.read_csv(fp, index_col=0, dtype=str)
        seen_tp = set(A["timepoint"].astype(str).str.strip())
        if seen_tp != {stage}:
            raise SystemExit(f"[h01] {fp.name}: expected timepoint {{{stage}}}, got {sorted(seen_tp)}")
        A["_pid"] = A["Patient_ID"].str.strip()
        A["_tp"] = A["timepoint"].str.strip()
        A["_ct"] = A["Classified_Celltype"].str.strip()
        A["_unit"] = atac_unit(A["sample_id"])
        n0 = len(A)
        A = A[~A["_unit"].isin(exclude)].copy()
        n_excluded = n0 - len(A)

        # the confidence gate the original build applied, counted before it is discarded
        conf = pd.to_numeric(A["Linked_scRNA_confidence_score"], errors="coerce")
        if conf.isna().any():
            raise SystemExit(f"[h01] {stage}: non-numeric Linked_scRNA_confidence_score")
        n_below_gate = int((conf < CONF_GATE).sum())

        new_bc = pd.Series(index=A.index, dtype=object)
        new_lib = pd.Series(index=A.index, dtype=object)
        dropped_idx = []

        for key, grp in A.groupby(["_pid", "_tp", "_ct"], sort=True):
            if key[1] != stage:
                raise SystemExit(f"[h01] stratum {key} carries a timepoint that is not {stage}")
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
        out["Linked_scRNA_confidence_score"] = 1.0   # placeholder; see the module docstring

        assert out["Linked_scRNA_barcode"].notna().all(), "unassigned cell survived"
        # the constraints, asserted rather than trusted
        tp_of_unit = dict(zip(R["_unit"], R["_tp"]))
        pid_of_unit = dict(zip(R["_unit"], R["_pid"]))
        lt = out["Linked_scRNA_Sample"].map(tp_of_unit)
        lp = out["Linked_scRNA_Sample"].map(pid_of_unit)
        assert lt.notna().all(), "a donor library is absent from the healthy RNA pool"
        assert (lt == out["_tp"]).all(), "TIMEPOINT CROSSED - constraint violated"
        assert (lp == out["_pid"]).all(), "PATIENT CROSSED - constraint violated"

        donors = out["Linked_scRNA_barcode"].nunique()
        reuse = len(out) / donors
        old = A["Linked_scRNA_barcode"]
        self_hit = float((out["Linked_scRNA_barcode"].to_numpy()
                          == A.loc[out.index, "Linked_scRNA_barcode"].to_numpy()).mean())
        old_within = float((A["Linked_scRNA_Sample"].map(tp_of_unit) == A["_tp"]).mean())

        dest = outdir / f"cell_data_{stage}.healthy_repaired.csv"
        out.drop(columns=["_pid", "_tp", "_ct", "_unit"]).to_csv(dest)
        summary[stage] = {
            "n_in": n0, "n_excluded_units": n_excluded, "n_out": len(out),
            "n_dropped_no_donor": len(dropped_idx), "pct_dropped": len(dropped_idx) / n0,
            "achieved_reuse": reuse, "distinct_donors": int(donors),
            "old_reuse": len(A) / old.nunique(), "old_distinct_donors": int(old.nunique()),
            "old_pct_within_timepoint": old_within, "new_pct_within_timepoint": 1.0,
            "agreement_with_old_link": self_hit,
            "n_below_conf_gate_now_admitted": n_below_gate,
            "file": dest.name}
        log(f"{stage}: {len(out):,}/{n0:,} cells kept "
            f"({len(dropped_idx):,} dropped no-donor {len(dropped_idx)/n0:.2%}"
            + (f", {n_excluded:,} excluded-unit" if n_excluded else "") + ") | "
            f"reuse {reuse:.3f}x (was {len(A)/old.nunique():.2f}x) | "
            f"donors {donors:,} (was {old.nunique():,}) | "
            f"within-timepoint 100.00% (was {old_within:.2%}) | "
            f"agreement with old link {self_hit:.2%}")
        log(f"{stage}:   {n_below_gate:,} cells were below the original {CONF_GATE} confidence "
            f"gate and are now ADMITTED (the placeholder is not behaviour-neutral here)")

    # Outputs are scope-tagged so a three-stage exploratory run does not overwrite the two-stage
    # outputs. The per-stage cell_data files are keyed by stage and need no suffix.
    sfx = ".withRemission" if a.include_remission else ""
    pd.DataFrame(per_stratum).to_csv(outdir / f"h01_per_stratum{sfx}.csv", index=False)
    (outdir / f"h01_healthy_assignment_report{sfx}.json").write_text(json.dumps(
        {"seed": a.seed, "scope": f"healthy {' + '.join(stages)}", "stages_built": stages,
         "excluded_units": exclude,
         "rna_pool": f"{HEALTHY_LABEL} only -- forced by the consumer (healthy h5ads), not chosen",
         "unassignable_policy": "DROP (no cell-type relaxation)",
         "confidence_placeholder": 1.0,
         "confidence_placeholder_is_behaviour_neutral": False,
         "timepoint_label_authority": "Patient_Sample (RNA) vs timepoint (ATAC); library token "
                                      "is NEVER used -- see the PAULSL swap in p_h00",
         "constraints": ["within-patient", "within-timepoint", "within-cell-type"],
         "stages": summary}, indent=2))
    log(f"wrote h01_healthy_assignment_report{sfx}.json and h01_per_stratum{sfx}.csv -> {outdir}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
