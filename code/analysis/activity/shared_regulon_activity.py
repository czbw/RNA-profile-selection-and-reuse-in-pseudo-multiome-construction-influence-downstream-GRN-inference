"""FIG1 panel c diagnostic - cell-type activity agreement on all shared regulators."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
L3 = ROOT / "longitudinal/published" / "out"
R3 = ROOT / "longitudinal/corrected" / "out"

STAGES = ["Diagnosis", "Relapse"]
A_ARM, B_ARM = "link", "pseudomultiome"
MIN_CELLTYPES_FLAG = 5

def log(m):
    print(f"[f1-ctagree] {m}", flush=True)

def corr_pair(x, y):
    """Pearson and Spearman with degenerate cases returned as NaN, never as 0."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan, np.nan, n
    return float(np.corrcoef(x, y)[0, 1]), float(stats.spearmanr(x, y).statistic), n

def dist(v, label):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return {"statistic": label, "n": 0}
    return {
        "statistic": label, "n": int(v.size),
        "min": float(v.min()), "p10": float(np.percentile(v, 10)),
        "q25": float(np.percentile(v, 25)), "median": float(np.median(v)),
        "q75": float(np.percentile(v, 75)), "p90": float(np.percentile(v, 90)),
        "max": float(v.max()), "mean": float(v.mean()),
        "frac_gt_0.5": float((v > 0.5).mean()), "frac_lt_0": float((v < 0).mean()),
    }

def fisher_mean(rs):
    rs = np.asarray([r for r in rs if np.isfinite(r)], float)
    if rs.size == 0:
        return np.nan
    rs = np.clip(rs, -0.999999, 0.999999)
    return float(np.tanh(np.arctanh(rs).mean()))

prof = pd.read_csv(OUT / "f1_celltype_profile.csv")
ctj = json.load(open(OUT / "f1_celltype.json"))

# module definitions, read from the funnels rather than restated
fl = pd.read_csv(L3 / "nomination_funnel.csv", index_col=0)
fr = pd.read_csv(R3 / "nomination_funnel.csv", index_col=0)
MODULES = {
    "link_T4_26": list(fl.index[fl.survives_4]),
    "repair_core34": list(fr.index[fr.ax4_tier1 & fr.ax1_directional & fr.ax2_selective
                                   & fr.ax3_celltype_robust]),
}

res = {
    "status": "DIAGNOSTIC - report only; no figure was built and nothing in figures/ was touched",
    "question": ("does the all-shared-regulator cell-type agreement hold at the level the "
                 "module-level number implies?"),
    "inputs": {
        "profile": "analysis/out/f1_celltype_profile.csv",
        "floor_and_celltypes": "analysis/out/f1_celltype.json",
        "modules": ["longitudinal/published/out/nomination_funnel.csv",
                    "longitudinal/corrected/out/nomination_funnel.csv"],
        "signmode_mismatch": "analysis/out/f1_signmode_mismatch.csv",
    },
    "design": {}, "distributions": {}, "modules": {}, "parity_gate": {},
    "checks": {},
}
rows = []

for stage in STAGES:
    P = prof[prof.stage == stage]
    wa = P[P.arm == A_ARM].pivot(index="celltype", columns="regulon",
                                 values="z_within_arm_across_celltypes")
    wb = P[P.arm == B_ARM].pivot(index="celltype", columns="regulon",
                                 values="z_within_arm_across_celltypes")
    ra = P[P.arm == A_ARM].pivot(index="celltype", columns="regulon", values="mean_AUC_raw")
    rb = P[P.arm == B_ARM].pivot(index="celltype", columns="regulon", values="mean_AUC_raw")

    ct_shared = sorted(set(wa.index) & set(wb.index))
    reg_shared = sorted(set(wa.columns) & set(wb.columns))
    res["design"][stage] = {
        "celltypes_retained_link": sorted(wa.index),
        "celltypes_retained_pseudomultiome": sorted(wb.index),
        "celltypes_shared": ct_shared,
        "n_celltypes_link": int(len(wa.index)),
        "n_celltypes_pseudomultiome": int(len(wb.index)),
        "n_celltypes_shared": int(len(ct_shared)),
        "celltypes_link_only": sorted(set(wa.index) - set(wb.index)),
        "celltypes_pseudomultiome_only": sorted(set(wb.index) - set(wa.index)),
        "excluded_below_20_unit_floor_link": ctj["excluded_below_floor"][f"{stage}_{A_ARM}"],
        "excluded_below_20_unit_floor_pseudomultiome":
            ctj["excluded_below_floor"][f"{stage}_{B_ARM}"],
        "n_regulons_link": int(wa.shape[1]),
        "n_regulons_pseudomultiome": int(wb.shape[1]),
        "n_regulons_shared": int(len(reg_shared)),
        "n_TFs_shared": int(len({r.split("_")[0] for r in reg_shared})),
    }
    log(f"{stage}: {len(ct_shared)} shared cell types, {len(reg_shared)} shared regulons, "
        f"{len({r.split('_')[0] for r in reg_shared})} shared TFs")

    # ---------------------------------------------------------------- (1) per regulon
    za, zb = wa.loc[ct_shared], wb.loc[ct_shared]
    Ra, Rb = ra.loc[ct_shared], rb.loc[ct_shared]
    max_abs_z_vs_raw = 0.0
    for r in reg_shared:
        pr, sr, n = corr_pair(za[r].values, zb[r].values)
        # check: the z transform cannot move a per-regulon correlation
        pr_raw, sr_raw, _ = corr_pair(Ra[r].values, Rb[r].values)
        if np.isfinite(pr) and np.isfinite(pr_raw):
            max_abs_z_vs_raw = max(max_abs_z_vs_raw, abs(pr - pr_raw))
        rows.append({"stage": stage, "level": "regulon", "unit": r,
                     "TF": r.split("_")[0], "sign_mode": r.split("_")[-1],
                     "aggregation": "none",
                     "pearson_r": pr, "spearman_r": sr, "n_celltypes": n,
                     "n_celltypes_below_5": bool(n < MIN_CELLTYPES_FLAG),
                     "pearson_r_on_raw_AUC": pr_raw})
    res["checks"][f"{stage}_max_abs_diff_pearson_z_vs_raw_per_regulon"] = max_abs_z_vs_raw

    reg_rows = [x for x in rows if x["stage"] == stage and x["level"] == "regulon"]
    res["distributions"].setdefault(stage, {})
    res["distributions"][stage]["regulon_pearson"] = dist(
        [x["pearson_r"] for x in reg_rows], "per-regulon Pearson r")
    res["distributions"][stage]["regulon_spearman"] = dist(
        [x["spearman_r"] for x in reg_rows], "per-regulon Spearman rho")
    res["distributions"][stage]["regulon_undefined"] = int(
        sum(1 for x in reg_rows if not np.isfinite(x["pearson_r"])))

    tf_shared_any = sorted(set(c.split("_")[0] for c in wa.columns)
                           & set(c.split("_")[0] for c in wb.columns))
    tf_with_shared_regulon = sorted({r.split("_")[0] for r in reg_shared})
    res["design"][stage]["n_TFs_shared_at_TF_level"] = int(len(tf_shared_any))
    res["design"][stage]["n_TFs_with_at_least_one_shared_regulon"] = int(len(tf_with_shared_regulon))
    res["design"][stage]["TFs_shared_but_sharing_no_regulon"] = sorted(
        set(tf_shared_any) - set(tf_with_shared_regulon))

    by_tf_shared = {t: [r for r in reg_shared if r.split("_")[0] == t]
                    for t in tf_shared_any}
    by_tf_a = {t: [c for c in wa.columns if c.split("_")[0] == t] for t in tf_shared_any}
    by_tf_b = {t: [c for c in wb.columns if c.split("_")[0] == t] for t in tf_shared_any}
    rmap = {x["unit"]: x["pearson_r"] for x in reg_rows}
    smap = {x["unit"]: x["spearman_r"] for x in reg_rows}

    for t in tf_shared_any:
        sh = by_tf_shared[t]
        if not sh:
            # no shared regulon: metrics 1 and 3 are undefined for this TF; recorded as such,
            # not as zero and not skipped.
            for tag in ("T1_mean_of_per_regulon_r", "T3_profile_mean_over_shared_regulons"):
                rows.append({"stage": stage, "level": "TF", "unit": t, "TF": t,
                             "sign_mode": "", "aggregation": tag,
                             "pearson_r": np.nan, "spearman_r": np.nan,
                             "n_celltypes": len(ct_shared), "n_celltypes_below_5": False,
                             "n_regulons_used_link": len(by_tf_a[t]),
                             "n_regulons_used_pm": len(by_tf_b[t]),
                             "note": "no regulon shared between the arms for this TF; "
                                     "this aggregation does not exist for it"})
            pa = za[by_tf_a[t]].mean(axis=1).values
            pb = zb[by_tf_b[t]].mean(axis=1).values
            pr2, sr2, n2 = corr_pair(pa, pb)
            rows.append({"stage": stage, "level": "TF", "unit": t, "TF": t, "sign_mode": "",
                         "aggregation": "T2_profile_mean_over_each_arms_own_regulons",
                         "pearson_r": pr2, "spearman_r": sr2, "n_celltypes": n2,
                         "n_celltypes_below_5": bool(n2 < MIN_CELLTYPES_FLAG),
                         "n_regulons_used_link": len(by_tf_a[t]),
                         "n_regulons_used_pm": len(by_tf_b[t]),
                         "signmode_count_mismatch": int(abs(len(by_tf_a[t]) - len(by_tf_b[t]))),
                         "DA26_exposed": True,
                         "note": "no shared regulon: the two profiles are averages over "
                                 "entirely disjoint sign modes"})
            continue
        # Metric 1: summary of the TF's per-regulon agreements (no profile averaged).
        rs = [rmap[r] for r in sh]
        ss = [smap[r] for r in sh]
        rs_f = [x for x in rs if np.isfinite(x)]
        rows.append({"stage": stage, "level": "TF", "unit": t, "TF": t, "sign_mode": "",
                     "aggregation": "T1_mean_of_per_regulon_r",
                     "pearson_r": float(np.mean(rs_f)) if rs_f else np.nan,
                     "spearman_r": float(np.mean([x for x in ss if np.isfinite(x)])) if rs_f else np.nan,
                     "n_celltypes": len(ct_shared),
                     "n_celltypes_below_5": bool(len(ct_shared) < MIN_CELLTYPES_FLAG),
                     "n_regulons_used_link": len(sh), "n_regulons_used_pm": len(sh),
                     "pearson_r_fisher_mean": fisher_mean(rs),
                     "pearson_r_median_of_regulons": float(np.median(rs_f)) if rs_f else np.nan})
        
        pa = za[by_tf_a[t]].mean(axis=1).values if by_tf_a[t] else None
        pb = zb[by_tf_b[t]].mean(axis=1).values if by_tf_b[t] else None
        pr2, sr2, n2 = corr_pair(pa, pb) if pa is not None and pb is not None else (np.nan, np.nan, 0)
        rows.append({"stage": stage, "level": "TF", "unit": t, "TF": t, "sign_mode": "",
                     "aggregation": "T2_profile_mean_over_each_arms_own_regulons",
                     "pearson_r": pr2, "spearman_r": sr2, "n_celltypes": n2,
                     "n_celltypes_below_5": bool(n2 < MIN_CELLTYPES_FLAG),
                     "n_regulons_used_link": len(by_tf_a[t]),
                     "n_regulons_used_pm": len(by_tf_b[t]),
                     "signmode_count_mismatch": int(abs(len(by_tf_a[t]) - len(by_tf_b[t]))),
                     "DA26_exposed": True})
        # Metric 3: profiles averaged over the shared regulons only (mismatch removed).
        pr3, sr3, n3 = corr_pair(za[sh].mean(axis=1).values, zb[sh].mean(axis=1).values)
        rows.append({"stage": stage, "level": "TF", "unit": t, "TF": t, "sign_mode": "",
                     "aggregation": "T3_profile_mean_over_shared_regulons",
                     "pearson_r": pr3, "spearman_r": sr3, "n_celltypes": n3,
                     "n_celltypes_below_5": bool(n3 < MIN_CELLTYPES_FLAG),
                     "n_regulons_used_link": len(sh), "n_regulons_used_pm": len(sh),
                     "DA26_exposed": True})

    for tag in ("T1_mean_of_per_regulon_r", "T2_profile_mean_over_each_arms_own_regulons",
                "T3_profile_mean_over_shared_regulons"):
        sub = [x for x in rows if x["stage"] == stage and x["level"] == "TF"
               and x["aggregation"] == tag]
        res["distributions"][stage][f"TF_pearson_{tag}"] = dist(
            [x["pearson_r"] for x in sub], f"per-TF Pearson r ({tag})")
        res["distributions"][stage][f"TF_spearman_{tag}"] = dist(
            [x["spearman_r"] for x in sub], f"per-TF Spearman rho ({tag})")

    # ---------------------------------------------------------------- (3) modules
    for mname, mregs in MODULES.items():
        pa_c = [r for r in mregs if r in wa.columns]
        pb_c = [r for r in mregs if r in wb.columns]
        if not pa_c or not pb_c:
            continue
        prm, srm, nm = corr_pair(za[pa_c].mean(axis=1).values, zb[pb_c].mean(axis=1).values)
        rows.append({"stage": stage, "level": "module", "unit": mname, "TF": "", "sign_mode": "",
                     "aggregation": "module_profile_mean_over_module_regulons",
                     "pearson_r": prm, "spearman_r": srm, "n_celltypes": nm,
                     "n_celltypes_below_5": bool(nm < MIN_CELLTYPES_FLAG),
                     "n_regulons_used_link": len(pa_c), "n_regulons_used_pm": len(pb_c),
                     "DA26_exposed": True})
        # the module's own members' per-regulon agreements, so the module-averaged value can be
        # read against its parts. Averaging regulons with heterogeneous cell-type profiles gives
        # a flatter mean profile than the parts, so this difference measures that loss.
        mem = [rmap[r] for r in mregs if r in rmap and np.isfinite(rmap[r])]
        res["modules"][f"{stage}|{mname}"] = {
            "pearson_r": prm, "spearman_r": srm, "n_celltypes": nm,
            "n_regulons_in_link": len(pa_c), "n_regulons_in_pseudomultiome": len(pb_c),
            "member_per_regulon_median": float(np.median(mem)) if mem else None,
            "n_members_with_a_per_regulon_value": int(len(mem)),
            "module_minus_member_median": (float(prm - np.median(mem)) if mem else None)}
        committed = ctj["module_profile_shape_agreement"].get(
            f"{stage}|{mname}|link_vs_pseudomultiome", {})
        if "pearson_r" in committed:
            res["parity_gate"][f"{stage}|{mname}"] = {
                "committed": committed["pearson_r"], "recomputed": prm,
                "abs_diff": abs(committed["pearson_r"] - prm),
                "PASS": bool(abs(committed["pearson_r"] - prm) < 1e-12)}

D = pd.DataFrame(rows)
D.to_csv(OUT / "f1_celltype_agreement_allshared.csv", index=False)

# sign-mode mismatch context, read not restated
sm = pd.read_csv(OUT / "f1_signmode_mismatch.csv")
res["signmode_mismatch_context"] = {
    s: {"n_shared_TFs": int((sm.stage == s).sum()),
        "frac_mode_count_differs": float((sm[sm.stage == s].mode_count_mismatch > 0).mean())}
    for s in STAGES}

res["parity_gate"]["ALL_PASS"] = bool(
    all(v["PASS"] for k, v in res["parity_gate"].items() if isinstance(v, dict)))

with open(OUT / "f1_celltype_agreement_allshared.json", "w") as f:
    json.dump(res, f, indent=1)

for stage in STAGES:
    dd = res["distributions"][stage]
    log(f"{stage}: per-regulon Pearson median {dd['regulon_pearson']['median']:+.4f} "
        f"(n={dd['regulon_pearson']['n']}, IQR {dd['regulon_pearson']['q25']:+.3f}.."
        f"{dd['regulon_pearson']['q75']:+.3f}, >0.5 {dd['regulon_pearson']['frac_gt_0.5']:.3f}, "
        f"<0 {dd['regulon_pearson']['frac_lt_0']:.3f})")
log(f"parity check: {res['parity_gate']['ALL_PASS']}")
log(f"max |Pearson(z) - Pearson(raw)| per regulon: "
    f"{max(res['checks'][k] for k in res['checks']):.3e}")
log("wrote f1_celltype_agreement_allshared.{csv,json}")
