#!/usr/bin/env python
'Cell-type mean-profile correlation, each BMMC arm vs native reference.'
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import activity_common as C

TAG = "b_07"
log = C.make_log(TAG)

SEED_BOOT, B = 555, 2000          # identical to p_09

def profiles(M, labels, order):
    """(cells x regulons) -> (regulons x n_celltypes) of within-cell-type mean activity."""
    return np.vstack([M[labels == c].mean(axis=0) for c in order]).T

def corr_matrix(Pa, Pt):
    'C[x, y] = Pearson(arm profile of regulon x, native reference profile of regulon y). Row-centred.'
    A = Pa - Pa.mean(axis=1, keepdims=True)
    T = Pt - Pt.mean(axis=1, keepdims=True)
    A = A / np.linalg.norm(A, axis=1, keepdims=True)
    T = T / np.linalg.norm(T, axis=1, keepdims=True)
    return A @ T.T

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
    X, cells, regs, mods = C.load_gene_based(TAG, h5mu)
    meta = C.load_metadata(TAG, meta_path)
    premise, row_of = C.verify_premise(TAG, h5mu, X, cells, regs, mods, meta, tables)
    ct = C.celltype_labels(TAG, meta)
    log(f"truth matrix {X.shape[0]:,} x {X.shape[1]}")

    rng = np.random.default_rng(SEED_BOOT)
    results, rows = {}, []
    for arm in arms:
        df = C.load_assignment(TAG, tables[arm], arm)
        own = np.array([row_of[b] for b in df["true_rna_barcode"]])
        got = np.array([row_of[b] for b in df["assigned_rna_barcode"]])
        lab = np.array([ct[b] for b in df["true_rna_barcode"]])
        order = C.celltype_order(lab)

        Pt = profiles(X[own], lab, order)          # native reference's cell-type profiles
        Pa = profiles(X[got], lab, order)          # the arm's
        keep = (Pt.std(axis=1) > 0) & (Pa.std(axis=1) > 0)
        Cm = corr_matrix(Pa[keep], Pt[keep])
        n = Cm.shape[0]
        matched = np.diag(Cm)
        off = Cm[~np.eye(n, dtype=bool)]

        pct = np.array([stats.percentileofscore(Cm[i], Cm[i, i], kind="rank") for i in range(n)])
        bs = np.array([np.median(matched[rng.integers(0, n, n)]) for _ in range(B)])

        if arm == C.TRUTH and not np.allclose(matched, 1.0, atol=1e-10):
            C.die(TAG, f"identity control FAILED: truth-vs-truth matched profile correlation is "
                       f"not 1.0 (min {matched.min()!r}); the barcode join is wrong")

        # ---- the cell-type-mean degeneracy -- detected from the data, not assumed -------------
        sct, degenerate = C.within_celltype_bijection(
            df, ct, premise["per_arm"][arm]["is_exact_bijection"])
        max_profile_gap = float(np.abs(Pa - Pt).max())
        if degenerate and not np.allclose(matched, 1.0, atol=1e-9):
            C.die(TAG, f"arm {arm!r} is an exact within-cell-type bijection, so its cell-type mean "
                       f"profile MUST equal truth's, yet matched correlation min is "
                       f"{matched.min()!r}. The derivation or the join is wrong -- one of them "
                       f"must be fixed before any number here is quoted.")
        if degenerate:
            log(f"    [{arm}] DEGENERATE at cell-type-mean resolution: same-cell-type rate "
                f"{sct:.6f}, exact bijection -> profiles identical to truth "
                f"(max |Pa-Pt| = {max_profile_gap:.3e}). Its matched correlation is an IDENTITY.")

        r = {"arm": arm, "ladder_arm": arm in C.LADDER_ARMS, "seed": a.seed,
             "same_cell_type_rate": sct,
             "is_exact_bijection": premise["per_arm"][arm]["is_exact_bijection"],
             "degenerate_by_construction": degenerate,
             "max_abs_profile_difference_vs_truth": max_profile_gap,
             "n_regulons_used": int(n), "n_regulons_dropped_constant": int((~keep).sum()),
             "n_celltypes": len(order), "n_cells": int(len(df)),
             "matched_median": float(np.median(matched)),
             "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
             "mismatched_null_median": float(np.median(off)),
             "null_p95": float(np.percentile(off, 95)),
             "identifiability_median_pct": float(np.median(pct)),
             "own_is_rank1_frac": float(np.mean(pct >= 100 * (n - 1) / n)),
             "frac_matched_above_0.5": float(np.mean(matched > 0.5))}
        results[arm] = r
        rows.append(pd.DataFrame({"arm": arm, "seed": a.seed,
                                  "regulon": np.array(regs)[keep], "r_profile": matched}))
        log(f"{arm:28s} matched {r['matched_median']:+.4f}"
            + (" (IDENTITY)" if degenerate else " ") + " "
            f"[{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}] | null {r['mismatched_null_median']:+.4f} "
            f"| ident {r['identifiability_median_pct']:.1f}% | rank1 {r['own_is_rank1_frac']:.1%} "
            f"| >0.5 {r['frac_matched_above_0.5']:.1%}"
            + ("" if arm in C.LADDER_ARMS else "   [non-ladder]"))

    blind, aware = results[C.DONOR_BLIND], results[C.DONOR_AWARE]
    both_degenerate = blind["degenerate_by_construction"] and aware["degenerate_by_construction"]
    donor_axis = {
        "donor_blind_arm": C.DONOR_BLIND, "donor_aware_arm": C.DONOR_AWARE,
        "contrast_note": ("one designed variable (donor-awareness), but donor+site-composite: "
                          "site is nested under donor for 8 of 10 BMMC donors, so the donor "
                          "constraint drags site concordance 0.3303 -> 0.7918 with it."),
        "both_arms_degenerate_by_construction": both_degenerate,
        "statistic_can_separate_the_two_arms": not both_degenerate,
        "degeneracy_note": (C.CELLTYPE_MEAN_DEGENERACY + " Both BMMC ladder random-within-cell-type "
                            "arms are exact within-cell-type bijections, so at THIS resolution the "
                            "donor delta below is identically zero and is NOT evidence that donor "
                            "awareness does not matter. It is measured at cell resolution by b_06 "
                            "(L6) and b_08 (L7).") if both_degenerate else None,
        "delta_matched_median": aware["matched_median"] - blind["matched_median"],
        "delta_identifiability_median_pct": (aware["identifiability_median_pct"]
                                             - blind["identifiability_median_pct"]),
        "delta_own_is_rank1_frac": aware["own_is_rank1_frac"] - blind["own_is_rank1_frac"],
        "delta_frac_matched_above_0.5": (aware["frac_matched_above_0.5"]
                                         - blind["frac_matched_above_0.5"]),
        "ci95_overlap": not (blind["ci95"][1] < aware["ci95"][0]
                             or aware["ci95"][1] < blind["ci95"][0]),
        "ci95_overlap_note": ("the two CIs are regulon-clustered bootstraps of two medians over "
                              "the SAME regulons on the SAME cells, so non-overlap is suggestive "
                              "of a difference but is not a paired test of one. No p-value is "
                              "quoted because none was pre-specified at this rung."),
    }

    log("")
    log("=" * 96)
    log("STATISTIC 3 -- cell-type mean-profile correlation vs TRUTH (the resolution used in the AML analysis)")
    log("=" * 96)
    log(f"{'arm':30s}{'matched':>10}{'95% CI':>21}{'null':>9}{'ident%':>9}{'rank1':>8}{'>0.5':>8}")
    for arm in arms:
        r = results[arm]
        ci = "[{:+.4f},{:+.4f}]".format(r["ci95"][0], r["ci95"][1])
        log(f"{arm:30s}{r['matched_median']:>+10.4f}{ci:>21}"
            f"{r['mismatched_null_median']:>+9.4f}{r['identifiability_median_pct']:>8.1f}%"
            f"{r['own_is_rank1_frac']:>8.1%}{r['frac_matched_above_0.5']:>8.1%}"
            + ("   <- IDENTITY, degenerate by construction"
               if r["degenerate_by_construction"] else ""))
    log("")
    log(f"DONOR CONTRAST  delta matched median (aware - blind) = "
        f"{donor_axis['delta_matched_median']:+.4f}")

    report = {
        "stage": a.stage, "seed": a.seed, "track": "A", "rung": "statistic_3_celltype_profile",
        "arms_scored": arms, "ladder_arms": C.LADDER_ARMS, "nonladder_arms": C.NONLADDER_ARMS,
        "estimator": {"correlation": "Pearson on cell-type mean-activity profiles",
                      "summary": "median over regulons",
                      "bootstrap": "regulon-clustered", "B": B, "bootstrap_seed": SEED_BOOT,
                      "n_points_per_profile": results[C.TRUTH]["n_celltypes"]},
        "declarations": C.declarations(),
        "declared_degeneracy": {
            "statement": C.CELLTYPE_MEAN_DEGENERACY,
            "detected_from_data": True,
            "degenerate_arms": [k for k, v in results.items()
                                if v["degenerate_by_construction"]],
            "informative_arms": [k for k, v in results.items()
                                 if not v["degenerate_by_construction"]],
            "consequence": ("the matched-profile correlation of a degenerate arm is an identity "
                            "and must not be quoted as a measurement of that arm's fidelity"),
        },
        "premise": premise,
        "arms": results,
        "donor_axis": donor_axis,
        "pbmc_counterpart": {
            "file": str(C.PBMC_OUT / "stat3_celltype_profile_10k.json"),
            "read_not_rerun": True,
            "values_seed555": {
                k.split("|")[0]: v
                for k, v in C.read_pbmc_committed("stat3_celltype_profile_10k.json").items()
                if isinstance(v, dict) and k.endswith("|555")},
            "note": ("different dataset, different arms, 16 cell types not 22, 245 regulons not "
                     "221, 3 seeds not 1. Arm-to-arm pairing across the datasets is an "
                     "interpretive choice and is not made here."),
        },
        "aml_comparator_note": ("p_09 logs the AML arm's activity-level figures (Option A vs Option B, "
                                "NO ground truth): healthy 0.844 / 0.900 / 0.862, malignant "
                                "0.820 / 0.710. Those compare two constructions to each other; "
                                "everything here compares a construction to KNOWN truth. Not the "
                                "same quantity -- quoted only as the number this rung exists to "
                                "put in context."),
        "provenance": C.provenance(
            "group_mean_activity.py",
            "benchmarks/pbmc/scripts/group_mean_activity.py",
            inputs={
                "truth_auc_h5mu": {"path": str(h5mu), "sha256": premise["h5mu_sha256"]},
                "cell_metadata": {"path": str(meta_path), "sha256": C.sha256(meta_path)},
                "arm_tables": {arm: premise["per_arm"][arm]["table"] for arm in arms},
                "arm_table_sha256": {arm: premise["per_arm"][arm]["table_sha256"] for arm in arms},
            }),
    }
    report["provenance"]["runtime_s"] = round(time.time() - t0, 1)
    report["provenance"]["peak_rss_bytes"] = C.peak_rss_bytes()

    if a.no_write:
        log("--no-write: nothing written")
        return
    pd.concat(rows).to_csv(out / f"b07_stat3_per_regulon_{a.stage}.csv", index=False)
    C.write_json(TAG, out / f"b07_stat3_celltype_profile_{a.stage}.json", report)
    log(f"b_07 COMPLETE in {report['provenance']['runtime_s']}s")

if __name__ == "__main__":
    main()
