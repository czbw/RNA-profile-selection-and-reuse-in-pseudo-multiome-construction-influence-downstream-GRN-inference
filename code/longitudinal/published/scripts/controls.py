"""Controls for the step-11 mixed-model results: sign-agreement null, permutation calibration, minimum detectable effect."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/published')
OUT = ROOT / "out"

SEED = 555
N_DRAWS_C1 = 2000
N_PERM_C2 = 20
N_REG_C2 = 100
N_REG_C3 = 60
DELTAS = [0.002, 0.005, 0.01, 0.02, 0.04, 0.08]
MIN_CELLS_PER_UNIT = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")

rng = np.random.default_rng(SEED)

def log(m):
    print(f"[s3-ctrl] {m}", flush=True)

def load(defs):
    t = pd.read_csv(OUT / f"pseudobulk_def{defs}.csv")
    pu = set(map(tuple, pd.read_csv(OUT / "paired_units.csv").values))
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(ADMITTED).values
    f = t[m].copy()
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    return f.sort_values(["patient", "celltype", "stage"]).reset_index(drop=True)

def fit_stage(d):
    """Return p for stage_num under the step-11 primary; None on failure."""
    for re_f, kind in (("~stage_num", "slope"), (None, "int")):
        try:
            r = smf.mixedlm("y ~ stage_num + C(celltype)", d, groups=d.patient,
                            re_formula=re_f).fit(reml=True, method="lbfgs", maxiter=200)
            b, se = r.params.get("stage_num"), r.bse.get("stage_num")
            if b is None or not np.isfinite(b) or not np.isfinite(se) or se == 0:
                continue
            if kind == "slope" and not r.converged:
                continue
            return float(r.pvalues["stage_num"]), float(b)
        except Exception:
            continue
    return None

frames = {d: load(d) for d in DEFS}
regs = {d: [c for c in frames[d].columns if c not in META + ["stage_num"]] for d in DEFS}
fits = pd.read_csv(OUT / "lmm_fits_by_definition.csv")
res = {"seed": SEED}

# C1 — mismatched-name sign-agreement null

SOURCES = {"raw": OUT / "lmm_fits_by_definition.csv",
           "offset_centred": OUT / "c4_offset_centred_fits.csv"}
res["C1_sign_agreement"] = {}
for label, path in SOURCES.items():
    if not path.exists():
        log(f"C1 {label}: {path.name} absent — skipped")
        continue
    src = pd.read_csv(path)
    w = src.pivot(index="regulon", columns="defs", values="beta")
    shared = [r for r in set(regs["Diagnosis"]) & set(regs["Relapse"]) if r in w.index]
    w = w.loc[shared].dropna()
    bD, bR = w["Diagnosis"].values, w["Relapse"].values
    obs = float(np.mean(np.sign(bD) == np.sign(bR)))

    fracD = float(np.mean(bD < 0))
    fracR = float(np.mean(bR < 0))
    analytic = fracD * fracR + (1 - fracD) * (1 - fracR)

    draws = np.empty(N_DRAWS_C1)
    for i in range(N_DRAWS_C1):
        draws[i] = np.mean(np.sign(bD) == np.sign(rng.permutation(bR)))
    p_emp = float((np.sum(draws >= obs) + 1) / (N_DRAWS_C1 + 1))

    res["C1_sign_agreement"][label] = {
        "n_shared_regulons": int(len(w)),
        "observed": obs,
        "marginal_frac_negative": {"Diagnosis_defs": fracD, "Relapse_defs": fracR},
        "analytic_chance_from_marginals": analytic,
        "permutation_null_mean": float(draws.mean()),
        "permutation_null_ci95": [float(np.percentile(draws, 2.5)),
                                  float(np.percentile(draws, 97.5))],
        "p_empirical": p_emp,
        "excess_over_null": obs - float(draws.mean()),
    }
    log(f"C1 [{label}] n={len(w)} observed {obs:.3f} | mismatched-name null "
        f"{draws.mean():.3f} [{np.percentile(draws, 2.5):.3f}, "
        f"{np.percentile(draws, 97.5):.3f}] | p={p_emp:.4g}")
    log(f"C1 [{label}] marginal negative: Dx defs {fracD:.3f}, Rel defs {fracR:.3f} "
        f"-> chance from marginals alone = {analytic:.3f}")
    np.save(OUT / f"c1_sign_agreement_null_{label}.npy", draws)

# C2 — within-patient stage-label permutation calibration

c2 = []
for d in DEFS:
    f = frames[d]
    pick = list(rng.choice(regs[d], size=min(N_REG_C2, len(regs[d])), replace=False))
    # unit = (patient, celltype); exactly two rows, one per stage
    key = f.patient + "|" + f.celltype
    units = key.unique()
    for perm in range(N_PERM_C2):
        flip = {u: bool(v) for u, v in zip(units, rng.integers(0, 2, len(units)))}
        sn = f.stage_num.values.copy()
        for u in units:
            if flip[u]:
                idx = np.where(key.values == u)[0]
                sn[idx] = 1.0 - sn[idx]
        base = f[["patient", "celltype"]].copy()
        base["stage_num"] = sn
        for reg in pick:
            base["y"] = f[reg].values
            r = fit_stage(base)
            if r:
                c2.append(dict(defs=d, perm=perm, regulon=reg, p=r[0], beta=r[1]))
    log(f"C2 def{d}: {N_PERM_C2} permutations x {len(pick)} regulons done")

c2 = pd.DataFrame(c2)
c2.to_csv(OUT / "c2_permutation_calibration.csv", index=False)
res["C2_permutation"] = {}
for d in DEFS:
    g = c2[c2.defs == d]
    res["C2_permutation"][d] = {
        "n_fits": int(len(g)),
        "frac_p_below_05": float((g.p < 0.05).mean()),
        "frac_p_below_01": float((g.p < 0.01).mean()),
        "median_p": float(g.p.median()),
        "ks_uniform_p": float(stats.kstest(g.p.values, "uniform").pvalue),
    }
    log(f"C2 def{d}: p<0.05 in {(g.p < 0.05).mean():.4f} of permuted fits "
        f"(nominal 0.05), median p {g.p.median():.3f}")

# C3 — minimum detectable effect

c3 = []
for d in DEFS:
    f = frames[d]
    n_tests = int((fits.defs == d).sum())
    p_first = 0.05 / n_tests          # p a regulon must clear to be the first BH discovery
    pick = list(rng.choice(regs[d], size=min(N_REG_C3, len(regs[d])), replace=False))
    for delta in DELTAS:
        hits, tot = 0, 0
        for reg in pick:
            base = f[["patient", "celltype", "stage_num"]].copy()
            base["y"] = f[reg].values + delta * f.stage_num.values
            r = fit_stage(base)
            if r:
                tot += 1
                hits += r[0] < p_first
        c3.append(dict(defs=d, delta=delta, n=tot, frac_detected=hits / tot if tot else np.nan,
                       p_threshold=p_first))
        log(f"C3 def{d} delta={delta:<6} detected {hits}/{tot} at p<{p_first:.2g}")
c3 = pd.DataFrame(c3)
c3.to_csv(OUT / "c3_minimum_detectable_effect.csv", index=False)

res["C3_mde"] = {}
for d in DEFS:
    g = c3[c3.defs == d].sort_values("delta")
    ok = g[g.frac_detected >= 0.80]
    sd = float(np.std([frames[d][r].mean() for r in regs[d]]))
    res["C3_mde"][d] = {
        "p_threshold_first_discovery": float(g.p_threshold.iloc[0]),
        "mde_80pct_power_auc_units": float(ok.delta.min()) if len(ok) else None,
        "sweep": g[["delta", "frac_detected"]].to_dict("records"),
        "between_regulon_sd_of_mean_activity": sd,
        "mde_in_sd_units": float(ok.delta.min() / sd) if len(ok) and sd else None,
    }

with open(OUT / "s3_controls.json", "w") as fh:
    json.dump(res, fh, indent=2)
log("wrote s3_controls.json")
