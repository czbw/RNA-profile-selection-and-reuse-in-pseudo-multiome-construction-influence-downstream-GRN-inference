"""FIG1 panel d, the accompanying report - is a regulon's activity agreement related to its"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
STAGES = ["Diagnosis", "Relapse"]

def log(m):
    print(f"[f1-actvtgt] {m}", flush=True)

act = pd.read_csv(OUT / "f1_celltype_agreement_allshared.csv")
act = act[(act.level == "regulon")][["stage", "unit", "pearson_r", "spearman_r"]]
act = act.rename(columns={"unit": "regulon", "pearson_r": "activity_pearson_r",
                          "spearman_r": "activity_spearman_r"})

tgt = pd.read_csv(OUT / "f1_shared_regulon_jaccard.csv")
tgt = tgt[(tgt.arm_A == "link") & (tgt.arm_B == "pseudomultiome")][
    ["stage", "regulon", "TF", "sign_mode", "jaccard_genes", "jaccard_regions",
     "n_genes_A", "n_genes_B"]]

def grid(v, label):
    """The percentile grid, computed by one function for both statistics so the two rows
    of the supplementary table are summarised identically."""
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    return {
        "statistic": label, "n": int(v.size),
        "min": float(v.min()), "p10": float(np.percentile(v, 10)),
        "q25": float(np.percentile(v, 25)), "median": float(np.median(v)),
        "q75": float(np.percentile(v, 75)), "p90": float(np.percentile(v, 90)),
        "max": float(v.max()), "mean": float(v.mean()),
        "frac_gt_0.5": float((v > 0.5).mean()), "frac_lt_0": float((v < 0).mean()),
        "n_gt_0.5": int((v > 0.5).sum()), "n_lt_0": int((v < 0).sum()),
    }

res = {"status": "DIAGNOSTIC report accompanying Figure 1 panel d",
       "inputs": {"activity": "analysis/out/f1_celltype_agreement_allshared.csv "
                              "(level == 'regulon')",
                  "target": "analysis/out/f1_shared_regulon_jaccard.csv "
                            "(arm_A == 'link', arm_B == 'pseudomultiome')"},
       "shared_substrate_caveat": (
           "The two currencies are NOT measured on independent substrates. AUCell scores a "
           "regulon FROM its target genes, so both arms' activity profiles for a regulon are "
           "computed from the very gene sets whose overlap is the target statistic. A regulon "
           "whose two arms share more targets will therefore tend to have more similar activity "
           "profiles for a partly MECHANICAL reason, and the baseline expectation for the "
           "correlation below is positive rather than zero. It is bounded, though: the largest "
           "target Jaccard reached by any shared regulon is reported in marginals.target_max, "
           "so even the most-overlapping regulon shares only that fraction of its target union."),
       "set_identity": {}, "grids": {}, "by_stage": {}}

joined = []
for stage in STAGES:
    a = act[act.stage == stage]
    t = tgt[tgt.stage == stage]
    sa, st_ = set(a.regulon), set(t.regulon)
    ident = {
        "n_activity": len(sa), "n_target": len(st_),
        "n_intersection": len(sa & st_),
        "sets_identical": bool(sa == st_),
        "in_activity_not_target": sorted(sa - st_),
        "in_target_not_activity": sorted(st_ - sa),
    }
    res["set_identity"][stage] = ident
    if not ident["sets_identical"]:
        raise SystemExit(
            f"REGULON-SET IDENTITY ASSERTION FAILED at {stage}: "
            f"{len(sa)} activity vs {len(st_)} target, "
            f"{len(sa - st_)} / {len(st_ - sa)} asymmetric. Panel d must not be built.")
    log(f"{stage}: set identity OK - {len(sa)} regulons, sets equal")

    m = a.merge(t, on=["stage", "regulon"], how="inner", validate="one_to_one")
    assert len(m) == len(sa), f"{stage}: join produced {len(m)} rows, expected {len(sa)}"
    joined.append(m)

    x = m.activity_pearson_r.values.astype(float)
    y = m.jaccard_genes.values.astype(float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]

    res["grids"][stage] = {
        "activity_pearson": grid(x, "cell-type activity agreement (Pearson r)"),
        "activity_spearman": grid(m.activity_spearman_r.values.astype(float)[ok],
                                  "cell-type activity agreement (Spearman rho)"),
        "target_jaccard_genes": grid(y, "target-set agreement (gene Jaccard)"),
        "target_jaccard_regions": grid(m.jaccard_regions.values.astype(float)[ok],
                                       "target-set agreement (region Jaccard)"),
    }
    pr = stats.pearsonr(x, y)
    sr = stats.spearmanr(x, y)
    # the same question on the Spearman activity statistic, as a sensitivity
    x2 = m.activity_spearman_r.values.astype(float)[ok]
    pr2 = stats.pearsonr(x2, y)
    sr2 = stats.spearmanr(x2, y)

    lo = m[m.activity_pearson_r < np.percentile(x, 25)]
    hi = m[m.activity_pearson_r > np.percentile(x, 75)]
    res["by_stage"][stage] = {
        "n": int(ok.sum()),
        "activity_pearson_r_vs_target_jaccard": {
            "pearson_r": float(pr.statistic), "pearson_p": float(pr.pvalue),
            "spearman_rho": float(sr.statistic), "spearman_p": float(sr.pvalue),
            "r_squared": float(pr.statistic ** 2)},
        "sensitivity_using_activity_spearman": {
            "pearson_r": float(pr2.statistic), "pearson_p": float(pr2.pvalue),
            "spearman_rho": float(sr2.statistic), "spearman_p": float(sr2.pvalue)},
        "marginals": {
            "activity_median": float(np.median(x)), "activity_min": float(x.min()),
            "activity_max": float(x.max()),
            "target_median": float(np.median(y)), "target_min": float(y.min()),
            "target_max": float(y.max())},
        "quartile_contrast": {
            "note": "target-set agreement in the bottom vs top quartile of ACTIVITY agreement; "
                    "if the two currencies were one quantity these would separate",
            "bottom_quartile_activity_n": int(len(lo)),
            "bottom_quartile_target_median": float(lo.jaccard_genes.median()),
            "top_quartile_activity_n": int(len(hi)),
            "top_quartile_target_median": float(hi.jaccard_genes.median()),
            "difference": float(hi.jaccard_genes.median() - lo.jaccard_genes.median())},
        "the_dissociation_in_one_line": {
            "note": "regulons in the TOP decile of activity agreement that are STILL below the "
                    "0.50 target bar - the panel's claim, counted",
            "n_top_decile_activity": int((m.activity_pearson_r
                                          >= np.percentile(x, 90)).sum()),
            "n_of_those_below_0.50_target": int(((m.activity_pearson_r
                                                  >= np.percentile(x, 90))
                                                 & (m.jaccard_genes < 0.5)).sum()),
            "max_target_jaccard_among_them": float(
                m.loc[m.activity_pearson_r >= np.percentile(x, 90), "jaccard_genes"].max())},
    }
    b = res["by_stage"][stage]
    log(f"{stage}: n={b['n']} | Pearson(activity, target) = "
        f"{b['activity_pearson_r_vs_target_jaccard']['pearson_r']:+.4f} "
        f"(p={b['activity_pearson_r_vs_target_jaccard']['pearson_p']:.3g}, "
        f"R2={b['activity_pearson_r_vs_target_jaccard']['r_squared']:.4f}) | Spearman = "
        f"{b['activity_pearson_r_vs_target_jaccard']['spearman_rho']:+.4f} "
        f"(p={b['activity_pearson_r_vs_target_jaccard']['spearman_p']:.3g})")
    log(f"    target median, bottom vs top activity quartile: "
        f"{b['quartile_contrast']['bottom_quartile_target_median']:.4f} vs "
        f"{b['quartile_contrast']['top_quartile_target_median']:.4f} "
        f"(diff {b['quartile_contrast']['difference']:+.4f})")
    log(f"    top-decile activity regulons still under the 0.50 target bar: "
        f"{b['the_dissociation_in_one_line']['n_of_those_below_0.50_target']} of "
        f"{b['the_dissociation_in_one_line']['n_top_decile_activity']}")

pd.concat(joined, ignore_index=True).to_csv(OUT / "f1_activity_vs_target.csv", index=False)
res["set_identity"]["ALL_IDENTICAL"] = bool(
    all(v["sets_identical"] for v in res["set_identity"].values() if isinstance(v, dict)))
with open(OUT / "f1_activity_vs_target.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote f1_activity_vs_target.{csv,json}")
