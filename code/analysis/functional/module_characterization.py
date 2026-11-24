#!/usr/bin/env python
"""The corrected arm's module, recomputed on the construction-independent 55 regulons."""
import json
import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
R3 = os.path.join(ROOT, "longitudinal/corrected", "out")
R3B = os.path.join(ROOT, "longitudinal/healthy_corrected", "out")
OUT = os.path.join(ROOT, "analysis", "out")
REP = os.path.join(ROOT, "adult/corrected", "out")

def log(m):
    print(f"[set55] {m}", flush=True)

fun = pd.read_csv(os.path.join(R3, "nomination_funnel.csv"))
sel = fun.ax1_directional & fun.ax2_selective & fun.ax3_celltype_robust
S55 = fun[sel].copy()
assert len(S55) == 55, len(S55)
log(f"the 55: {len(S55)} regulons over {S55.TF.nunique()} TFs")

rec = json.load(open(os.path.join(R3, "s3_recalibrated.json"), encoding="utf-8"))
LMM26 = rec["results"]["E1_random_intercept"]["confirmed_regulons"]

fam = pd.read_csv(os.path.join(OUT, "f1_family_assignment.csv")).set_index("TF")
S55 = S55.join(fam, on="TF")

# ---- effect size and the ax4 counterfactual -------------------------------------------
S55["abs_beta_max"] = S55[["beta_mal_Dx", "beta_mal_Rel"]].abs().max(axis=1)
S55["ax4_would_remove"] = ~S55.ax4_tier1
S55["n_axes_passed_of_3"] = 3

for voc in ("Diagnosis", "Relapse"):
    p = os.path.join(REP, f"per_regulon_external_T3_repair_{voc}.csv")
    e = pd.read_csv(p)
    key = "regulon" if "regulon" in e.columns else e.columns[0]
    val = next((c for c in e.columns if c.lower() in
                ("oriented_mean_delta", "mean_delta", "delta", "mean", "effect")), None)
    if val is None:
        log(f"external {voc}: no recognised effect column in {os.path.basename(p)}; columns are "
            f"{list(e.columns)} -- left blank")
        S55[f"external_delta_{voc}"] = pd.NA
    else:
        keep = [key, val] + [c for c in ("frac_positive", "p", "agrees_with_internal_sign")
                             if c in e.columns]
        ren = {key: "regulon", val: f"external_delta_{voc}"}
        ren.update({c: f"external_{c}_{voc}" for c in keep[2:]})
        S55 = S55.merge(e[keep].rename(columns=ren), on="regulon", how="left")
        log(f"external {voc}: joined {S55[f'external_delta_{voc}'].notna().sum()} of {len(S55)}")

# ---- cell-type distribution of module activity ----------------------------------------
ct = pd.read_csv(os.path.join(OUT, "f1_celltype_profile.csv"))
ctr = ct[(ct.arm == "repair") & (ct.regulon.isin(set(S55.regulon)))]
ct_summary = (ctr.groupby(["stage", "celltype"])
                 .agg(n_regulons=("regulon", "nunique"), mean_AUC=("mean_AUC_raw", "mean"))
                 .reset_index())
log(f"cell-type profile rows for the 55: {len(ctr)} over {ctr.celltype.nunique()} cell types")

# ---- sign-mode survival: TFs present in all four modes, within the 55 ------------------
def modes(regs):
    t = {}
    for r in regs:
        t.setdefault(r.split("_")[0], set()).add(r.split("_")[-1])
    return t

m55 = modes(list(S55.regulon))
four55 = sorted(tf for tf, ms in m55.items() if len(ms) == 4)
m_univ = modes(list(fun.regulon))
four_univ = sorted(tf for tf, ms in m_univ.items() if len(ms) == 4)
S55["n_sign_modes_in_set"] = S55.TF.map({tf: len(ms) for tf, ms in m55.items()})
log(f"all-four-sign-mode TFs WITHIN the 55: {four55}")

# ---- set relations against the LMM-confirmed 26 ----------------------------------------
r55, l26 = set(S55.regulon), set(LMM26)
t55 = set(S55.TF)
t26 = {r.split("_")[0] for r in LMM26}
inter_r, inter_t = sorted(r55 & l26), sorted(t55 & t26)
jac_r = len(r55 & l26) / len(r55 | l26)
jac_t = len(t55 & t26) / len(t55 | t26)
ieg_def = set(json.load(open(os.path.join(ROOT, "figures", "data", "manifest.json"),
                             encoding="utf-8"))["s3b.ap1.ieg_definition"]["value"])
n_ieg55 = int(S55.TF.isin(ieg_def).sum())
n_ieg26 = sum(1 for r in LMM26 if r.split("_")[0] in ieg_def)
log(f"55 vs LMM26: {len(inter_r)} regulons, {len(inter_t)} TFs, Jaccard {jac_r:.4f}")
log(f"AP-1/IEG: 55 -> {n_ieg55}/{len(S55)} ({n_ieg55/len(S55):.1%}) | "
    f"LMM26 -> {n_ieg26}/{len(LMM26)} ({n_ieg26/len(LMM26):.1%})")

res = {
    "_provenance": {
        "definition": "ax1_directional & ax2_selective & ax3_celltype_robust, corrected arm",
        "note": "ax4_tier1 dropped; see analysis/T4_DEFINITION_STRESS_TEST.md",
        "funnel": "longitudinal/corrected/out/nomination_funnel.csv",
        "family_source": "analysis/out/f1_family_assignment.csv "
                         "(Lambert DBD level A, project family level B, project IEG list level C)",
        "external": "adult/corrected/out/per_regulon_external_T3_repair_{vocab}.csv",
        "celltype": "analysis/out/f1_celltype_profile.csv (arm == repair)"},
    "set55": {"n_regulons": len(S55), "n_TFs": int(S55.TF.nunique()),
              "regulons": sorted(r55), "TFs": sorted(t55)},
    "vs_LMM26": {"n_intersection_regulons": len(inter_r), "intersection_regulons": inter_r,
                 "n_intersection_TFs": len(inter_t), "intersection_TFs": inter_t,
                 "jaccard_regulon": jac_r, "jaccard_TF": jac_t},
    "ieg": {"n_in_set55": n_ieg55, "frac_in_set55": n_ieg55 / len(S55),
            "n_in_LMM26": n_ieg26, "frac_in_LMM26": n_ieg26 / len(LMM26),
            "ieg_list_size": len(ieg_def)},
    "ax4_counterfactual": {
        "n_would_be_removed": int(S55.ax4_would_remove.sum()),
        "removed_regulons": sorted(S55.regulon[S55.ax4_would_remove]),
        "median_abs_beta_kept": float(S55.abs_beta_max[~S55.ax4_would_remove].median()),
        "median_abs_beta_removed": float(S55.abs_beta_max[S55.ax4_would_remove].median())},
    "sign_modes": {"all_four_within_set55": four55,
                   "all_four_within_repair_universe201": four_univ},
    "family_composition_level_B": S55.level_B_family.value_counts().to_dict(),
    "family_composition_TF_level": (S55.drop_duplicates("TF")
                                    .level_B_family.value_counts().to_dict()),
    "celltype_summary": ct_summary.to_dict(orient="records"),
}

cols = ["regulon", "TF", "sign_mode_lead", "direction", "n_axes_passed_of_3",
        "ax1_directional", "ax2_selective", "ax3_celltype_robust", "ax4_tier1",
        "ax4_would_remove", "n_sign_modes_in_set", "beta_mal_Dx", "beta_mal_Rel",
        "abs_beta_max", "beta_healthy_Diagnosis", "beta_healthy_Relapse",
        "did_Diagnosis", "did_Relapse", "q_did_Diagnosis", "q_did_Relapse",
        "n_celltypes_same_sign_Diagnosis", "n_celltypes_same_sign_Relapse",
        "level_A_Lambert_DBD", "level_B_family", "level_C_project_IEG", "tier"]
cols += [c for c in S55.columns if c.startswith("external_")]
S55[[c for c in cols if c in S55.columns]].sort_values("abs_beta_max", ascending=False) \
    .to_csv(os.path.join(OUT, "f6_set55_per_regulon.csv"), index=False)
json.dump(res, open(os.path.join(OUT, "f6_set55.json"), "w", encoding="utf-8"), indent=1)
log("wrote analysis/out/f6_set55_per_regulon.csv and f6_set55.json")
