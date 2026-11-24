"""Module test (see s3_23) re-run on frozen module definitions read from modules/out/ (arm-specific T3 and the cross-arm intersection T4)."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--arm", required=True, choices=["link", "repair"])
ARM = ap.parse_args().arm

REPO = Path(__file__).resolve().parents[3]
ARMS = {
    "link":   dict(root=REPO / "longitudinal/published",  s3b=REPO / "longitudinal/healthy_published" / "out",
                   t3=REPO / "modules" / "out" / "published_module.json",   t3name="T3_link"),
    "repair": dict(root=REPO / "longitudinal/corrected", s3b=REPO / "longitudinal/healthy_corrected" / "out",
                   t3=REPO / "modules" / "out" / "corrected_module.json", t3name="T3_repair"),
}
CFG = ARMS[ARM]
ROOT = CFG["root"]
OUT = ROOT / "out"
S3B = CFG["s3b"]
T4_JSON = REPO / "modules" / "out" / "shared_module.json"

SEED = 555
MIN_CELLS = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
N_PERM = 2000
N_RAND = 2000
rng = np.random.default_rng(SEED)

def log(m):
    print(f"[s3-module-{ARM}] {m}", flush=True)

def die(m):
    raise SystemExit(f"[s3-module-{ARM}] ABORT — {m}")

# Module definitions are read from file, never re-derived

fun = pd.read_csv(OUT / "nomination_funnel.csv", index_col=0)
UNIVERSE = list(fun.index)

t3_def = json.loads(CFG["t3"].read_text(encoding="utf-8"))
t4_def = json.loads(T4_JSON.read_text(encoding="utf-8"))
MODULES = {t3_def["name"]: list(t3_def["regulons"]),
           t4_def["name"]: list(t4_def["regulons"])}
MODULE_FILES = {t3_def["name"]: CFG["t3"], t4_def["name"]: T4_JSON}
MODULE_FROZEN = {t3_def["name"]: t3_def["frozen"], t4_def["name"]: t4_def["frozen"]}

# Entry assertions: the module files must agree with the funnel
derived_t3 = set(fun.index[fun.ax1_directional & fun.ax2_selective & fun.ax3_celltype_robust])
if set(MODULES[t3_def["name"]]) != derived_t3:
    die(f"{t3_def['name']} on disk ({len(MODULES[t3_def['name']])}) != ax1&ax2&ax3 on this arm's "
        f"funnel ({len(derived_t3)}); the module file and the funnel disagree")
for k, v in MODULES.items():
    missing = [r for r in v if r not in fun.index]
    if missing:
        die(f"{k}: {len(missing)} regulon(s) absent from the {ARM} funnel universe "
            f"({len(UNIVERSE)}), e.g. {missing[:3]}")
if not set(MODULES[t4_def["name"]]) <= set(MODULES[t3_def["name"]]):
    die("T4_intersection is not a subset of this arm's T3 — the intersection is inconsistent")

log(f"arm={ARM}  funnel universe {len(UNIVERSE)}  (module background for NULL B)")
for k, v in MODULES.items():
    log(f"{k}: {len(v)} regulons / {fun.loc[v, 'TF'].nunique()} TFs "
        f"[{MODULE_FILES[k].name}, frozen {MODULE_FROZEN[k][:10]}]")
# T4 is arm-keyed: record which arm's vocabulary is used; the two arms are never merged
log(f"T4 orientation/targets are arm-keyed; this run uses key '{ARM}' only")

# Data — unchanged from s3_23

def load(path, paired_path, admitted):
    t = pd.read_csv(path)
    pu = set(map(tuple, pd.read_csv(paired_path).values))
    m = (t.n_cells >= MIN_CELLS).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(admitted).values
    f = t[m].copy().reset_index(drop=True)
    regs = [c for c in f.columns if c not in META]
    A = f[regs].values
    return f, A - A.mean(axis=1, keepdims=True), regs        # offset-centred

HEAL_ADM = sorted(pd.read_csv(S3B / "paired_units_healthy.csv").celltype
                  .value_counts().pipe(lambda s: s[s >= 5]).index)
arms = {}
for d in DEFS:
    mf, mc, mr = load(OUT / f"pseudobulk_def{d}.csv", OUT / "paired_units.csv", ADMITTED)
    hf, hc, hr = load(S3B / f"pseudobulk_healthy_def{d}.csv",
                      S3B / "paired_units_healthy.csv", HEAL_ADM)
    assert mr == hr
    arms[d] = dict(mf=mf, mc=mc, hf=hf, hc=hc, regs=mr,
                   ix={r: j for j, r in enumerate(mr)})
    log(f"def{d}: malignant {len(mf)} units / {mf.patient.nunique()} patients | "
        f"healthy {len(hf)} units / {hf.patient.nunique()} patients | "
        f"{len(mr)} regulon columns")

# Core machinery — unchanged from s3_23

def unit_deltas(f, M, stage_num):
    d = pd.DataFrame(M, columns=[f"r{j}" for j in range(M.shape[1])])
    d["patient"], d["celltype"], d["s"] = f.patient.values, f.celltype.values, stage_num
    dx = d[d.s == 0.0].set_index(["patient", "celltype"]).drop(columns="s")
    rl = d[d.s == 1.0].set_index(["patient", "celltype"]).drop(columns="s")
    ix = dx.index.intersection(rl.index)
    return (rl.loc[ix] - dx.loc[ix])

def patient_matrix(delta_units):
    return delta_units.groupby(level="patient").mean()

def lopo_scores(pm):
    out = {}
    pats = list(pm.index)
    V = pm.values
    for i, p in enumerate(pats):
        others = np.delete(V, i, axis=0)
        s = np.sign(others.mean(axis=0))
        s[s == 0] = 1.0
        out[p] = float(np.mean(V[i] * s))
    return pd.Series(out)

def test_stat(scores):
    v = scores.values
    if len(v) < 5 or np.allclose(v, 0):
        return None
    return dict(n_patients=int(len(v)), mean=float(v.mean()), median=float(np.median(v)),
                frac_up=float((v > 0).mean()),
                p_wilcoxon=float(stats.wilcoxon(v, zero_method="wilcox").pvalue))

res = {"arm": ARM,
       "seed": SEED, "n_perm": N_PERM, "n_random_modules": N_RAND,
       "orientation_primary": "leave-one-patient-out",
       "status": "fork of s3_23_module_test.py for the 2026-08-20 module scheme; "
                 "the modules are pre-defined and frozen in modules/out/, this step tests them",
       "provenance": {
           "forked_from": "longitudinal/published/scripts/s3_23_module_test.py",
           "committed_outputs_not_touched": ["out/module_test.csv", "out/s3_module_test.json"],
           "arm": ARM,
           "funnel": str((OUT / "nomination_funnel.csv").relative_to(REPO)),
           "funnel_universe_size": len(UNIVERSE),
           "malignant_pseudobulk": [str((OUT / f"pseudobulk_def{d}.csv").relative_to(REPO))
                                    for d in DEFS],
           "healthy_pseudobulk": [str((S3B / f"pseudobulk_healthy_def{d}.csv").relative_to(REPO))
                                  for d in DEFS],
           "paired_units": str((OUT / "paired_units.csv").relative_to(REPO)),
           "paired_units_healthy": str((S3B / "paired_units_healthy.csv").relative_to(REPO)),
           "module_files": {k: str(v.relative_to(REPO)) for k, v in MODULE_FILES.items()},
           "module_frozen": MODULE_FROZEN,
           "T4_arm_key_used": ARM,
           "T4_arm_key_note": "T4's orientation/targets are keyed [arm][vocabulary]; only the "
                              "'%s' key is read in this run. The two arms are never merged or "
                              "averaged." % ARM,
           "admitted_malignant_celltypes": ADMITTED,
           "admitted_malignant_celltypes_note":
               "held at the committed six on both arms because the repair funnel's ax3 (>=4 of 6) "
               "was computed over these six; the repair arm's seventh >=5-patient cell type (CLP) "
               "is therefore excluded here, as it was there",
           "admitted_healthy_celltypes": HEAL_ADM,
           "min_cells_per_unit": MIN_CELLS,
           "nulls_redrawn": "A, B and C are all recomputed here at the size of the module being "
                            "tested; no null on disk is reused (the committed .npy/.json nulls are "
                            "matched to the OLD 26-regulon module)",
       },
       "null_C_pool_diagnostics": {
           "what_this_is": "candidate pool sizes at every funnel depth on THIS arm, and the "
                           "fraction of the pool each module would occupy if that depth were used "
                           "as the null C background. Reported so a decision about redefining or "
                           "retiring null C has numbers. NO judgement is made here and no pool "
                           "other than the committed s3_23 one (ax1&ax2&ax3) is used for any "
                           "p-value in this file.",
           "pools": {name: int(mask.sum()) for name, mask in {
               "universe": pd.Series(True, index=fun.index),
               "ax1": fun.ax1_directional,
               "ax1_and_ax2": fun.ax1_directional & fun.ax2_selective,
               "ax1_and_ax2_and_ax3": fun.ax1_directional & fun.ax2_selective
                                      & fun.ax3_celltype_robust,
           }.items()},
           "draw_fraction_module_over_pool": {
               k: {name: round(len(v) / int(mask.sum()), 4) for name, mask in {
                   "universe": pd.Series(True, index=fun.index),
                   "ax1": fun.ax1_directional,
                   "ax1_and_ax2": fun.ax1_directional & fun.ax2_selective,
                   "ax1_and_ax2_and_ax3": fun.ax1_directional & fun.ax2_selective
                                          & fun.ax3_celltype_robust,
               }.items()} for k, v in MODULES.items()},
       },
       "modules": {k: {"n_regulons": len(v), "n_TFs": int(fun.loc[v, "TF"].nunique()),
                       "TFs": sorted(fun.loc[v, "TF"].unique().tolist()),
                       "regulons": sorted(v)}
                   for k, v in MODULES.items()},
       "results": {}}

for mod_name, mod_regs in MODULES.items():
    res["results"][mod_name] = {}
    for d in DEFS:
        a = arms[d]
        cols = [a["ix"][r] for r in mod_regs if r in a["ix"]]
        if len(cols) != len(mod_regs):
            die(f"{mod_name} / def{d}: only {len(cols)} of {len(mod_regs)} module regulons are "
                f"columns of pseudobulk_def{d}.csv")
        sn_m = (a["mf"].stage == "Relapse").astype(float).values
        sn_h = (a["hf"].stage == "Relapse").astype(float).values

        # ---- PRIMARY: malignant, LOPO orientation ----
        du = unit_deltas(a["mf"], a["mc"][:, cols], sn_m)
        pm = patient_matrix(du)
        obs = lopo_scores(pm)
        st = test_stat(obs)

        # ---- NULL A: stage-label permutation, LOPO inside the loop ----
        key = (a["mf"].patient + "|" + a["mf"].celltype).values
        units = pd.unique(key)
        r2 = np.random.default_rng(SEED)
        null_a = np.empty(N_PERM)
        for b in range(N_PERM):
            sn = sn_m.copy()
            for u, flip in zip(units, r2.integers(0, 2, len(units))):
                if flip:
                    idx = np.where(key == u)[0]
                    sn[idx] = 1.0 - sn[idx]
            pmb = patient_matrix(unit_deltas(a["mf"], a["mc"][:, cols], sn))
            null_a[b] = lopo_scores(pmb).mean()
        p_a = float((np.sum(np.abs(null_a) >= abs(st["mean"])) + 1) / (N_PERM + 1))

        # ---- NULL B: size-matched random modules from the funnel universe ----
        all_shared = [a["ix"][r] for r in fun.index if r in a["ix"]]
        du_all = unit_deltas(a["mf"], a["mc"][:, all_shared], sn_m)
        pm_all = patient_matrix(du_all)
        r3 = np.random.default_rng(SEED + 1)
        null_b = np.empty(N_RAND)
        for b in range(N_RAND):
            pick = r3.choice(pm_all.shape[1], size=len(cols), replace=False)
            null_b[b] = lopo_scores(pm_all.iloc[:, pick]).mean()
        p_b = float((np.sum(np.abs(null_b) >= abs(st["mean"])) + 1) / (N_RAND + 1))

        # ---- NULL C: random modules drawn from the funnel survivors at step 3 (strictest) ----
        pool3 = [a["ix"][r] for r in fun.index[fun.ax1_directional & fun.ax2_selective
                                               & fun.ax3_celltype_robust] if r in a["ix"]]
        n_draw = min(len(cols), len(pool3))
        # Degeneracy is decided before drawing: if the module is the whole pool the null is not
        # computed, so no tautological p-value is produced.
        c_degenerate = bool(n_draw >= len(pool3))
        if c_degenerate:
            null_c, p_c = None, None
            null_c_block = {
                "background": f"{ARM} ax1&ax2&ax3 survivors (s3_23 definition, unchanged)",
                "pool_size": len(pool3), "draw_size": n_draw, "n_draws": 0,
                "status": "undefined_by_construction_module_is_the_pool",
                "why": ("the module under test IS the null C pool (draw size == pool size), so "
                        "every draw is the whole pool, the null is a constant equal to the "
                        "observation, and any p-value would be 1.0 as an arithmetic identity, not "
                        "as evidence. NOT COMPUTED. No replacement pool is substituted here."),
                "p_empirical": None, "mean": None, "sd": None, "ci95": None,
                "percentile_of_observed": None}
        else:
            du3 = unit_deltas(a["mf"], a["mc"][:, pool3], sn_m)
            pm3 = patient_matrix(du3)
            r4 = np.random.default_rng(SEED + 2)
            null_c = np.empty(N_RAND)
            for b in range(N_RAND):
                pick = r4.choice(pm3.shape[1], size=n_draw, replace=False)
                null_c[b] = lopo_scores(pm3.iloc[:, pick]).mean()
            p_c = float((np.sum(np.abs(null_c) >= abs(st["mean"])) + 1) / (N_RAND + 1))
            null_c_block = {
                "background": f"{ARM} ax1&ax2&ax3 survivors (s3_23 definition, unchanged)",
                "pool_size": len(pool3), "draw_size": n_draw, "n_draws": N_RAND,
                "status": "computed",
                "mean": float(null_c.mean()), "sd": float(null_c.std()),
                "ci95": [float(np.percentile(null_c, 2.5)), float(np.percentile(null_c, 97.5))],
                "p_empirical": p_c,
                "percentile_of_observed": float(stats.percentileofscore(np.abs(null_c),
                                                                        abs(st["mean"])))}

        # ---- SELECTIVITY: same module on the patient's own healthy cells ----
        s_full = np.sign(pm.values.mean(axis=0))
        s_full[s_full == 0] = 1.0
        duh = unit_deltas(a["hf"], a["hc"][:, cols], sn_h)
        pmh = patient_matrix(duh)
        heal = pd.Series((pmh.values * s_full).mean(axis=1), index=pmh.index)
        mal_full = pd.Series((pm.values * s_full).mean(axis=1), index=pm.index)
        sth = test_stat(heal)
        shared_pat = [p for p in mal_full.index if p in heal.index]
        did = mal_full.loc[shared_pat].values - heal.loc[shared_pat].values
        std = dict(n_patients=len(shared_pat), mean=float(did.mean()),
                   frac_malignant_larger=float((did > 0).mean()),
                   p_wilcoxon=float(stats.wilcoxon(did, zero_method="wilcox").pvalue))

        # ---- SECONDARY: cross-vocabulary orientation ----
        other = "Relapse" if d == "Diagnosis" else "Diagnosis"
        ao = arms[other]
        cols_o = [ao["ix"][r] for r in mod_regs if r in ao["ix"]]
        pm_o = patient_matrix(unit_deltas(ao["mf"], ao["mc"][:, cols_o],
                                          (ao["mf"].stage == "Relapse").astype(float).values))
        s_o = np.sign(pm_o.values.mean(axis=0))
        s_o[s_o == 0] = 1.0
        xv = pd.Series((pm.values * s_o).mean(axis=1), index=pm.index)
        stx = test_stat(xv)

        res["results"][mod_name][d] = {
            "primary_LOPO": st,
            "null_A_stage_permutation": {
                "n_draws": N_PERM, "module_size": len(cols),
                "mean": float(null_a.mean()), "sd": float(null_a.std()),
                "ci95": [float(np.percentile(null_a, 2.5)), float(np.percentile(null_a, 97.5))],
                "p_empirical": p_a},
            "null_B_random_modules": {
                "background": f"{ARM} funnel universe", "background_size": pm_all.shape[1],
                "draw_size": len(cols), "n_draws": N_RAND,
                "mean": float(null_b.mean()), "sd": float(null_b.std()),
                "ci95": [float(np.percentile(null_b, 2.5)), float(np.percentile(null_b, 97.5))],
                "p_empirical": p_b,
                "percentile_of_observed": float(stats.percentileofscore(np.abs(null_b),
                                                                        abs(st["mean"])))},
            "null_C_random_from_filtered_pool": null_c_block,
            "healthy_same_module": sth,
            "selectivity_DiD": std,
            "secondary_cross_vocabulary_orientation": stx,
        }
        log(f"{mod_name} / def{d}: LOPO mean {st['mean']:+.5f} "
            f"({st['frac_up']:.0%} of {st['n_patients']} patients same direction, "
            f"Wilcoxon p={st['p_wilcoxon']:.4g})")
        log(f"    NULL A stage-permutation p={p_a:.4g} [{N_PERM} draws, size {len(cols)}] | "
            f"NULL B random-module p={p_b:.4g} "
            f"({res['results'][mod_name][d]['null_B_random_modules']['percentile_of_observed']:.1f}th pct)"
            f" [{N_RAND} draws of {len(cols)} from {pm_all.shape[1]}]"
            + (f" | NULL C NOT COMPUTED — undefined by construction, the module IS the pool "
               f"({n_draw} of {len(pool3)})"
               if c_degenerate else
               f" | NULL C random-from-filtered p={p_c:.4g} "
               f"({res['results'][mod_name][d]['null_C_random_from_filtered_pool']['percentile_of_observed']:.1f}th pct)"
               f" [{N_RAND} draws of {n_draw} from {len(pool3)}]"))
        log(f"    healthy same module: mean {sth['mean']:+.5f} (p={sth['p_wilcoxon']:.3g}) | "
            f"DiD {std['mean']:+.5f} (p={std['p_wilcoxon']:.4g}, "
            f"{std['frac_malignant_larger']:.0%} of patients malignant>healthy)")
        log(f"    cross-vocabulary orientation: mean {stx['mean']:+.5f} (p={stx['p_wilcoxon']:.4g})")

    both = res["results"][mod_name]
    c_deg = any(both[d]["null_C_random_from_filtered_pool"]["status"] != "computed" for d in DEFS)
    res["results"][mod_name]["agreement"] = {
        "same_sign_both_vocabularies":
            bool(np.sign(both["Diagnosis"]["primary_LOPO"]["mean"]) ==
                 np.sign(both["Relapse"]["primary_LOPO"]["mean"])),
        "both_clear_null_A": bool(both["Diagnosis"]["null_A_stage_permutation"]["p_empirical"] < 0.05
                                  and both["Relapse"]["null_A_stage_permutation"]["p_empirical"] < 0.05),
        "both_clear_null_B": bool(both["Diagnosis"]["null_B_random_modules"]["p_empirical"] < 0.05
                                  and both["Relapse"]["null_B_random_modules"]["p_empirical"] < 0.05),
        "both_clear_null_C": (None if c_deg else bool(
            both["Diagnosis"]["null_C_random_from_filtered_pool"]["p_empirical"] < 0.05
            and both["Relapse"]["null_C_random_from_filtered_pool"]["p_empirical"] < 0.05)),
        "null_C_status": ("undefined_by_construction_module_is_the_pool" if c_deg else "computed"),
        "both_selective": bool(both["Diagnosis"]["selectivity_DiD"]["p_wilcoxon"] < 0.05
                               and both["Relapse"]["selectivity_DiD"]["p_wilcoxon"] < 0.05),
    }
    ag = res["results"][mod_name]["agreement"]
    log(f"{mod_name} AGREEMENT: same sign {ag['same_sign_both_vocabularies']}, "
        f"null A both {ag['both_clear_null_A']}, null B both {ag['both_clear_null_B']}, "
        f"null C both {ag['both_clear_null_C']} (status={ag['null_C_status']}), "
        f"selective both {ag['both_selective']}")

rows = []
for m, byd in res["results"].items():
    for d in DEFS:
        x = byd[d]
        rows.append(dict(arm=ARM, module=m, vocabulary=d,
                         n_regulons=len(MODULES[m]),
                         n_TFs=int(fun.loc[MODULES[m], "TF"].nunique()),
                         n_patients=x["primary_LOPO"]["n_patients"],
                         lopo_mean=x["primary_LOPO"]["mean"],
                         frac_patients_same_dir=x["primary_LOPO"]["frac_up"],
                         p_wilcoxon=x["primary_LOPO"]["p_wilcoxon"],
                         p_null_A=x["null_A_stage_permutation"]["p_empirical"],
                         null_A_n_draws=x["null_A_stage_permutation"]["n_draws"],
                         p_null_B=x["null_B_random_modules"]["p_empirical"],
                         null_B_background_size=x["null_B_random_modules"]["background_size"],
                         null_B_draw_size=x["null_B_random_modules"]["draw_size"],
                         null_B_n_draws=x["null_B_random_modules"]["n_draws"],
                         p_null_C=x["null_C_random_from_filtered_pool"]["p_empirical"],
                         null_C_status=x["null_C_random_from_filtered_pool"]["status"],
                         null_C_pool_size=x["null_C_random_from_filtered_pool"]["pool_size"],
                         null_C_draw_size=x["null_C_random_from_filtered_pool"]["draw_size"],
                         null_C_n_draws=x["null_C_random_from_filtered_pool"]["n_draws"],
                         healthy_mean=x["healthy_same_module"]["mean"],
                         healthy_p=x["healthy_same_module"]["p_wilcoxon"],
                         did_mean=x["selectivity_DiD"]["mean"],
                         did_p=x["selectivity_DiD"]["p_wilcoxon"],
                         xvocab_mean=x["secondary_cross_vocabulary_orientation"]["mean"],
                         xvocab_p=x["secondary_cross_vocabulary_orientation"]["p_wilcoxon"]))
tag = "_".join(MODULES.keys())
csv_path = OUT / f"module_test_{tag}.csv"
json_path = OUT / f"s3_module_test_{tag}.json"
pd.DataFrame(rows).to_csv(csv_path, index=False)
with open(json_path, "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log(f"wrote {json_path.name} + {csv_path.name}")
