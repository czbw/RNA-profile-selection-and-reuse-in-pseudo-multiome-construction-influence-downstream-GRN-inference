#!/usr/bin/env python
"L6: native reference's frozen regulon vocabulary scored on every BMMC arm's paired data."
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import activity_common as C

TAG = "b_06"
log = C.make_log(TAG)

def spearman_cols(A, B):
    """Per-column Spearman between two (n x p) matrices. NaN where a column is constant.

    Byte-for-byte the p_08 estimator, including the `np.ptp == 0` guard: a constant column has no
    defined rank correlation and writing 0.0 there would invent a measured zero.
    """
    out = np.full(A.shape[1], np.nan)
    for j in range(A.shape[1]):
        a, b = A[:, j], B[:, j]
        if np.ptp(a) == 0 or np.ptp(b) == 0:
            continue
        out[j] = stats.spearmanr(a, b).statistic
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default=C.STAGE_DEFAULT)
    ap.add_argument("--seed", type=int, default=C.SEED_DEFAULT)
    ap.add_argument("--uploads", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args()

    t0 = time.time()
    out = Path(a.out_dir) if a.out_dir else C.BUNDLE / "out"
    uploads = Path(a.uploads) if a.uploads else C.BUNDLE / "uploads"
    meta_path = Path(a.metadata) if a.metadata else uploads / "b00_cell_metadata.csv"
    h5mu = out / f"eRegulon_AUC_{a.stage}_truth.h5mu"

    arms = C.LADDER_ARMS + C.NONLADDER_ARMS
    tables = {arm: C.resolve_table(TAG, uploads, a.stage, a.seed, arm) for arm in arms}
    for arm, p in tables.items():
        log(f"arm {arm:28s} -> {p.name}")

    X, cells, regs, mods = C.load_gene_based(TAG, h5mu)
    meta = C.load_metadata(TAG, meta_path)
    premise, row_of = C.verify_premise(TAG, h5mu, X, cells, regs, mods, meta, tables)
    ct = C.celltype_labels(TAG, meta)
    log(f"truth activity matrix: {X.shape[0]:,} cells x {X.shape[1]} regulons "
        f"({X.nbytes / 2**20:.1f} MiB)")

    results, per_reg, per_ct = {}, [], []
    for arm in arms:
        df = C.load_assignment(TAG, tables[arm], arm)
        own = np.array([row_of[b] for b in df["true_rna_barcode"]])
        got = np.array([row_of[b] for b in df["assigned_rna_barcode"]])
        T, A = X[own], X[got]                    # native reference vs arm: same cells, same vocabulary
        labels = np.array([ct[b] for b in df["true_rna_barcode"]])
        order = C.celltype_order(labels)

        rho_all = spearman_cols(A, T)
        within = {c: spearman_cols(A[labels == c], T[labels == c]) for c in order}
        W = np.vstack([within[c] for c in order])          # (n_celltypes x n_regulons)

        # ---- the identity control -------------------------------------------------
        if arm == C.TRUTH:
            nonconst = ~np.isnan(rho_all)
            worst = float(np.nanmin(rho_all)) if nonconst.any() else float("nan")
            if not np.allclose(rho_all[nonconst], 1.0, atol=1e-12):
                C.die(TAG, f"identity control FAILED: truth-vs-truth rho is not 1.0 for every "
                           f"non-constant regulon (min {worst!r}). The barcode join is wrong.")
            log(f"identity control PASS: truth vs truth rho == 1.0 for all "
                f"{int(nonconst.sum())} non-constant regulons "
                f"({int((~nonconst).sum())} constant, rho undefined)")

        res = {
            "arm": arm,
            "ladder_arm": arm in C.LADDER_ARMS,
            "seed": a.seed,
            "n_cells": int(len(df)),
            "n_cells_dropped_unscoreable": premise["per_arm"][arm]["n_cells_dropped_unscoreable"],
            "n_regulons": int(X.shape[1]),
            "n_regulons_rho_all_undefined": int(np.isnan(rho_all).sum()),
            "n_celltypes": len(order),
            "rho_all_median": float(np.nanmedian(rho_all)),
            "rho_within_median": float(np.nanmedian(W)),
            "rho_within_q25": float(np.nanquantile(W, .25)),
            "rho_within_q75": float(np.nanquantile(W, .75)),
            "frac_within_above_0.5": float(np.nanmean(W > .5)),
            "frac_within_below_0.2": float(np.nanmean(W < .2)),
            "attenuation": float(np.nanmedian(rho_all) - np.nanmedian(W)),
            "self_hit_rate": float((df["assigned_rna_barcode"] == df["true_rna_barcode"]).mean()),
        }
        results[arm] = res
        per_reg.append(pd.DataFrame({"arm": arm, "seed": a.seed, "regulon": regs,
                                     "rho_all": rho_all,
                                     "rho_within_median": np.nanmedian(W, axis=0)}))
        for i, c in enumerate(order):
            per_ct.append({"arm": arm, "ladder_arm": arm in C.LADDER_ARMS, "seed": a.seed,
                           "cell_type": c, "n_cells": int((labels == c).sum()),
                           "rho_within_median": float(np.nanmedian(W[i])),
                           "rho_within_q25": float(np.nanquantile(W[i], .25)),
                           "rho_within_q75": float(np.nanquantile(W[i], .75)),
                           "frac_above_0.5": float(np.nanmean(W[i] > .5))})
        log(f"{arm:28s} rho_all {res['rho_all_median']:+.4f} | within {res['rho_within_median']:+.4f} "
            f"[{res['rho_within_q25']:+.4f},{res['rho_within_q75']:+.4f}] | "
            f"atten {res['attenuation']:+.4f} | selfhit {res['self_hit_rate']:.6f}"
            + ("" if arm in C.LADDER_ARMS else "   [non-ladder]"))

    # ---- the donor contrast: one designed variable, read at both resolutions -----------------
    per_ct_df = pd.DataFrame(per_ct)
    blind, aware = results[C.DONOR_BLIND], results[C.DONOR_AWARE]
    bct = per_ct_df[per_ct_df.arm == C.DONOR_BLIND].set_index("cell_type")["rho_within_median"]
    act = per_ct_df[per_ct_df.arm == C.DONOR_AWARE].set_index("cell_type")["rho_within_median"]
    d_ct = (act - bct).sort_values(ascending=False)
    donor_axis = {
        "donor_blind_arm": C.DONOR_BLIND,
        "donor_aware_arm": C.DONOR_AWARE,
        "contrast_note": ("the two arms are both exact bijections within cell type over the same "
                          "69,249 cells and differ in exactly one designed variable, "
                          "donor-awareness. b_05 measures 80.26% of the donor-blind arm's cells "
                          "receiving a different donor's transcriptome. NOT a pure-donor "
                          "contrast: site is nested under donor for 8 of 10 BMMC donors, so "
                          "constraining within donor drags site concordance 0.3303 -> 0.7918 and "
                          "batch 0.1433 -> 0.7918 with it. Read every delta "
                          "below as donor+site-composite."),
        "delta_rho_all_median": aware["rho_all_median"] - blind["rho_all_median"],
        "delta_rho_within_median": aware["rho_within_median"] - blind["rho_within_median"],
        "delta_attenuation": aware["attenuation"] - blind["attenuation"],
        "delta_frac_within_above_0.5": (aware["frac_within_above_0.5"]
                                        - blind["frac_within_above_0.5"]),
        "per_celltype_delta_rho_within_median": {k: float(v) for k, v in d_ct.items()},
        "n_celltypes_aware_higher": int((d_ct > 0).sum()),
        "n_celltypes_blind_higher": int((d_ct < 0).sum()),
    }
    log("")
    log("=" * 92)
    log("DONOR CONTRAST (one designed variable; donor+site-composite -- see the note)")
    log("=" * 92)
    log(f"  {C.DONOR_BLIND:30s} within-celltype rho {blind['rho_within_median']:+.4f}")
    log(f"  {C.DONOR_AWARE:30s} within-celltype rho {aware['rho_within_median']:+.4f}")
    log(f"  delta (aware - blind)                         {donor_axis['delta_rho_within_median']:+.4f}"
        f"   [aware higher in {donor_axis['n_celltypes_aware_higher']} of {len(d_ct)} cell types]")

    report = {
        "stage": a.stage, "seed": a.seed, "track": "A", "rung": "L6",
        "arms_scored": arms,
        "ladder_arms": C.LADDER_ARMS,
        "nonladder_arms": C.NONLADDER_ARMS,
        "declarations": C.declarations(),
        "premise": premise,
        "arms": results,
        "donor_axis": donor_axis,
        "pbmc_counterpart": {
            "file": str(C.PBMC_OUT / "trackA_L6_10k.json"),
            "read_not_rerun": True,
            "values_seed555": {
                k.split("|")[0]: {f: v[f] for f in
                                  ("n_cells", "n_cells_dropped_unscoreable", "n_regulons",
                                   "n_celltypes", "rho_all_median", "rho_within_median",
                                   "attenuation", "self_hit_rate")}
                for k, v in C.read_pbmc_committed("trackA_L6_10k.json").items()
                if isinstance(v, dict) and k.endswith("|555")},
            "note": ("PBMC 10k is a DIFFERENT DATASET with DIFFERENT ARMS -- one donor, 16 cell "
                     "types, 245 regulons, 10,881 analysis cells, 3 seeds. No BMMC arm is the "
                     "same construction as any PBMC rand_* arm, so any arm-to-arm pairing across "
                     "the two is an INTERPRETIVE CHOICE and is not made here. The one directly "
                     "comparable pair is link vs link (same construction family), and even that "
                     "differs in that PBMC's link dropped 650 unscoreable cells and BMMC's drops "
                     "none."),
        },
        "provenance": C.provenance(
            "activity_agreement.py", "benchmarks/pbmc/scripts/activity_agreement.py",
            inputs={
                "truth_auc_h5mu": {"path": str(h5mu), "sha256": premise["h5mu_sha256"]},
                "cell_metadata": {"path": str(meta_path), "sha256": C.sha256(meta_path)},
                "arm_tables": {arm: premise["per_arm"][arm]["table"] for arm in arms},
                "arm_table_sha256": {arm: premise["per_arm"][arm]["table_sha256"] for arm in arms},
            },
            extra={"runtime_s": None},
        ),
    }
    report["provenance"]["runtime_s"] = round(time.time() - t0, 1)
    report["provenance"]["peak_rss_bytes"] = C.peak_rss_bytes()

    if a.no_write:
        log("--no-write: nothing written")
        return
    pd.concat(per_reg).to_csv(out / f"b06_trackA_L6_per_regulon_{a.stage}.csv", index=False)
    per_ct_df.to_csv(out / f"b06_trackA_L6_per_celltype_{a.stage}.csv", index=False)
    C.write_json(TAG, out / f"b06_trackA_L6_{a.stage}.json", report)
    log(f"wrote {out}/b06_trackA_L6_per_regulon_{a.stage}.csv and "
        f"b06_trackA_L6_per_celltype_{a.stage}.csv")
    log(f"b_06 COMPLETE in {report['provenance']['runtime_s']}s, peak RSS "
        f"{(report['provenance']['peak_rss_bytes'] or 0) / 2**20:.0f} MiB")

if __name__ == "__main__":
    main()
