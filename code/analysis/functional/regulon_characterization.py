"""FIG6 (b)(d)(e) - the per-regulon table for the corrected arm, and the selectivity / external axes."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
R3 = ROOT / "longitudinal/corrected" / "out"
R3B = ROOT / "longitudinal/healthy_corrected" / "out"
RR = ROOT / "adult/corrected" / "out"
OUT = ROOT / "analysis" / "out"

def log(m):
    print(f"[f6-reg] {m}", flush=True)

fun = pd.read_csv(R3 / "nomination_funnel.csv")
rec = json.load(open(R3 / "s3_recalibrated.json"))
LMM26 = rec["results"]["E1_random_intercept"]["confirmed_regulons"]
e1 = pd.read_csv(R3 / "e1_random_intercept_reconciliation.csv").set_index("regulon")
did = pd.read_csv(R3B / "did_reconciliation.csv").set_index("regulon")
fam = pd.read_csv(OUT / "f1_family_assignment.csv").set_index("TF")
ct = pd.read_csv(OUT / "f1_celltype_profile.csv")

core = fun[fun.ax4_tier1 & fun.ax1_directional & fun.ax2_selective & fun.ax3_celltype_robust]
CORE34 = list(core.regulon)
assert len(CORE34) == 34, len(CORE34)

# ---- four-sign-mode survivors, at two scopes, computed not assumed
def modes(regs):
    t = {}
    for r in regs:
        t.setdefault(r.split("_")[0], set()).add(r.split("_")[-1])
    return t
m_core = modes(CORE34)
m_univ = modes(list(fun.regulon))
ALL4_CORE = sorted(k for k, v in m_core.items() if len(v) == 4)
ALL4_UNIV = sorted(k for k, v in m_univ.items() if len(v) == 4)
log(f"all four sign modes WITHIN the corrected core: {ALL4_CORE}")
log(f"all four sign modes within the corrected arm's 201-regulon universe: {ALL4_UNIV}")

# ---- external replication, per regulon
ext = {}
for d in ("Diagnosis", "Relapse"):
    e = pd.read_csv(RR / f"per_regulon_external_CORE34_{d}.csv").set_index("regulon")
    ext[d] = e
    log(f"external per-regulon {d}: {len(e)} regulons")

# ---- cell-type distribution summary, corrected arm
ctr = ct[ct.arm == "repair"]
top_ct = {}
for d in ("Diagnosis", "Relapse"):
    s = ctr[ctr.stage == d]
    for r, g in s.groupby("regulon"):
        g = g.sort_values("z_within_arm_across_celltypes", ascending=False)
        top_ct[(d, r)] = ("; ".join(f"{a}({b:+.2f})" for a, b in
                                    zip(g.celltype.head(3), g.z_within_arm_across_celltypes.head(3))),
                          "; ".join(f"{a}({b:+.2f})" for a, b in
                                    zip(g.celltype.tail(2), g.z_within_arm_across_celltypes.tail(2))))

rows = []
for r in sorted(set(CORE34) | set(LMM26)):
    tf = r.split("_")[0]
    f = fun[fun.regulon == r]
    f = f.iloc[0] if len(f) else None
    row = {
        "regulon": r, "TF": tf, "sign_mode": r.split("_")[-1],
        "in_corrected_core34": r in CORE34,
        "in_LMM_confirmed26": r in LMM26,
        "in_both": (r in CORE34) and (r in LMM26),
        "TF_all_four_sign_modes_in_core34": tf in ALL4_CORE,
        "TF_all_four_sign_modes_in_repair_universe201": tf in ALL4_UNIV,
        "TF_n_sign_modes_in_core34": len(m_core.get(tf, set())),
        "TF_n_sign_modes_in_universe201": len(m_univ.get(tf, set())),
        "family_A_Lambert_DBD": fam.level_A_Lambert_DBD.get(tf),
        "family_B_refined": fam.level_B_family.get(tf),
        "family_C_project_IEG": bool(fam.level_C_project_IEG.get(tf, False)),
    }
    if f is not None:
        row.update({
            "in_repair_universe201": True,
            "ax1_directional": bool(f.ax1_directional), "ax2_selective": bool(f.ax2_selective),
            "ax3_celltype_robust": bool(f.ax3_celltype_robust), "ax4_tier1": bool(f.ax4_tier1),
            "ax5_chip": "UNAVAILABLE BY CONSTRUCTION (ChIP gate computed against the PUBLISHED "
                        "network's targets; axis dropped)",
            "n_axes_passed_of_4": int(bool(f.ax1_directional) + bool(f.ax2_selective)
                                      + bool(f.ax3_celltype_robust) + bool(f.ax4_tier1)),
            "tier_imported_from_published_catalogue": f.tier,
            "direction": f.direction, "sign_mode_lead": f.sign_mode_lead, "tf_output": f.tf_output,
            "beta_malignant_Diagnosis": float(f.beta_malignant_Diagnosis),
            "beta_malignant_Relapse": float(f.beta_malignant_Relapse),
            "beta_healthy_Diagnosis": float(f.beta_healthy_Diagnosis),
            "beta_healthy_Relapse": float(f.beta_healthy_Relapse),
            "did_Diagnosis": float(f.did_Diagnosis), "did_Relapse": float(f.did_Relapse),
            "p_did_Diagnosis": float(f.p_did_Diagnosis), "p_did_Relapse": float(f.p_did_Relapse),
            "q_did_Diagnosis": float(f.q_did_Diagnosis), "q_did_Relapse": float(f.q_did_Relapse),
            "selectivity_verdict": f.verdict,
            "n_celltypes_same_sign_Diagnosis_of6": f.n_celltypes_same_sign_Diagnosis,
            "n_celltypes_same_sign_Relapse_of6": f.n_celltypes_same_sign_Relapse,
        })
    else:
        row["in_repair_universe201"] = False
    if r in e1.index:
        row.update({"E1_LMM_beta_Diagnosis": float(e1.loc[r, "beta_Diagnosis"]),
                    "E1_LMM_beta_Relapse": float(e1.loc[r, "beta_Relapse"]),
                    "E1_LMM_q_Diagnosis": float(e1.loc[r, "q_Diagnosis"]),
                    "E1_LMM_q_Relapse": float(e1.loc[r, "q_Relapse"]),
                    "E1_LMM_sign_agree": bool(e1.loc[r, "sign_agree"]),
                    "E1_LMM_confirmed": bool(e1.loc[r, "confirmed"])})
    else:
        row["E1_LMM_q_Diagnosis"] = None
        row["E1_LMM_note"] = "no E1 LMM fit exists for this regulon (not in the 136-regulon shared set)"
    for d in ("Diagnosis", "Relapse"):
        if r in ext[d].index:
            row[f"external59_oriented_delta_{d}"] = float(ext[d].loc[r, "oriented_mean_delta"])
            row[f"external59_frac_positive_{d}"] = float(ext[d].loc[r, "frac_positive"])
            row[f"external59_p_{d}"] = float(ext[d].loc[r, "p"])
            row[f"external59_agrees_internal_sign_{d}"] = bool(ext[d].loc[r, "agrees_with_internal_sign"])
        else:
            row[f"external59_oriented_delta_{d}"] = None
            row[f"external59_note_{d}"] = ("no external per-regulon value: the 59-patient external "
                                           "scoring was run for the corrected core only")
        tc = top_ct.get((d, r))
        row[f"top3_celltypes_by_z_{d}"] = tc[0] if tc else None
        row[f"bottom2_celltypes_by_z_{d}"] = tc[1] if tc else None
    rows.append(row)

T = pd.DataFrame(rows)
T.to_csv(OUT / "f6_per_regulon_table.csv", index=False)

# ---- does a module-level test / selectivity exist for the 34? verified, not asserted
hits = []
for p in ROOT.rglob("*.json"):
    if any(x in str(p) for x in (".git", "node_modules")):
        continue
    try:
        txt = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        continue
    if ("null_C" in txt or "primary_LOPO" in txt or "selectivity_DiD" in txt):
        hits.append(str(p.relative_to(ROOT)).replace("\\", "/"))
modules_named = {}
for h in hits:
    try:
        j = json.load(open(ROOT / h))
    except Exception:
        continue
    names = list(j.get("modules", {}).keys()) if isinstance(j.get("modules"), dict) else []
    sizes = {n: j["modules"][n].get("n_regulons") for n in names} if names else {}
    modules_named[h] = sizes

res = {
    "status": "extraction with provenance; nothing here is a new statistic except the two "
              "four-sign-mode lists and the family join",
    "prespec": "analysis/FIG1_FIG6_PRESPEC.md 4(b), 4(d), 4(e)",
    "sets": {
        "corrected_core34": {"n_regulons": len(CORE34), "n_TFs": len({r.split('_')[0] for r in CORE34}),
                             "TFs": sorted({r.split('_')[0] for r in CORE34}), "regulons": sorted(CORE34)},
        "LMM_confirmed26": {"n_regulons": len(LMM26), "n_TFs": len({r.split('_')[0] for r in LMM26}),
                            "TFs": sorted({r.split('_')[0] for r in LMM26}), "regulons": sorted(LMM26)},
        "intersection": {"n_regulons": len(set(CORE34) & set(LMM26)),
                         "regulons": sorted(set(CORE34) & set(LMM26)),
                         "n_TFs": len({r.split('_')[0] for r in CORE34} & {r.split('_')[0] for r in LMM26}),
                         "TFs": sorted({r.split('_')[0] for r in CORE34} & {r.split('_')[0] for r in LMM26}),
                         "jaccard_regulon": len(set(CORE34) & set(LMM26)) / len(set(CORE34) | set(LMM26))},
    },
    "four_sign_mode_survivors": {"within_corrected_core34": ALL4_CORE,
                                 "within_repair_universe201": ALL4_UNIV,
                                 "blueprint_DA57_7_says": "SPI1 only",
                                 "verdict": ("the analysis protocol statement is INCOMPLETE: "
                                             + ", ".join(ALL4_CORE) + " both survive in all four "
                                             "sign modes within the corrected core")},
    "ax5_status": "UNAVAILABLE BY CONSTRUCTION - not a failure",
    "module_level_artifacts_search": {
        "files_containing_null_C_or_primary_LOPO_or_selectivity_DiD": modules_named,
    },
    "IEG_check_on_LMM26": {
        "n_LMM26_regulons_whose_TF_is_in_the_project_IEG_list":
            int(sum(1 for r in LMM26 if bool(fam.level_C_project_IEG.get(r.split("_")[0], False)))),
        "n_core34_regulons_whose_TF_is_in_the_project_IEG_list":
            int(sum(1 for r in CORE34 if bool(fam.level_C_project_IEG.get(r.split("_")[0], False)))),
        "note": "an earlier note warns the 26 are dominated by the family axis 2 exists to remove; this "
                "counts it rather than repeating it on trust",
    },
}
with open(OUT / "f6_per_regulon.json", "w") as f:
    json.dump(res, f, indent=1)
log(f"core34 n IEG regulons = {res['IEG_check_on_LMM26']['n_core34_regulons_whose_TF_is_in_the_project_IEG_list']}, "
    f"LMM26 n IEG regulons = {res['IEG_check_on_LMM26']['n_LMM26_regulons_whose_TF_is_in_the_project_IEG_list']}")
log(f"core34 n LMM26 intersection = {res['sets']['intersection']['n_regulons']}")
log("wrote f6_per_regulon_table.csv, f6_per_regulon.json")
