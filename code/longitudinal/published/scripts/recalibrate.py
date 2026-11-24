"""Post hoc re-estimation, run after the permutation calibration (control C2) failed for the primary model."""
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
MIN_CELLS_PER_UNIT = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
N_PERM = 20
N_REG_PERM = 100

rng = np.random.default_rng(SEED)

def log(m):
    print(f"[s3-recal] {m}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def load(defs):
    t = pd.read_csv(OUT / f"pseudobulk_def{defs}.csv")
    pu = set(map(tuple, pd.read_csv(OUT / "paired_units.csv").values))
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(ADMITTED).values
    f = t[m].copy()
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    return f.reset_index(drop=True)

def e1(d):
    """(1 | patient) LMM -> (beta, p)."""
    try:
        r = smf.mixedlm("y ~ stage_num + C(celltype)", d, groups=d.patient).fit(
            reml=True, method="lbfgs", maxiter=200)
        b, se = r.params.get("stage_num"), r.bse.get("stage_num")
        if b is None or not np.isfinite(b) or not np.isfinite(se) or se == 0:
            return None
        return float(b), float(r.pvalues["stage_num"])
    except Exception:
        return None

def e2(d):
    """Per-patient paired Wilcoxon on within-patient deltas -> (mean delta, p)."""
    p = d.pivot_table(index=["patient", "celltype"], columns="stage_num", values="y")
    if 0.0 not in p or 1.0 not in p:
        return None
    p = p.dropna()
    per = (p[1.0] - p[0.0]).groupby(level="patient").mean()
    v = per.values
    if len(v) < 5 or np.allclose(v, 0):
        return None
    return float(np.mean(v)), float(stats.wilcoxon(v, zero_method="wilcox").pvalue)

frames, centred, regs = {}, {}, {}
for d in DEFS:
    f = load(d)
    r = [c for c in f.columns if c not in META + ["stage_num"]]
    A = f[r].values
    frames[d], regs[d] = f, r
    centred[d] = A - A.mean(axis=1, keepdims=True)
    log(f"def{d}: {len(f)} rows, {len(r)} regulons (offset-centred)")

res = {"seed": SEED, "status": "POST-HOC — triggered by the C2 calibration failure",
       "predeclared_primary_result_stands": "0 regulons CONFIRMED (s3_results.json)"}

# Calibrate first: an estimator is read only after its permutation null is measured

res["calibration"] = {}
for name, fn in (("E1_random_intercept", e1), ("E2_paired_wilcoxon", e2)):
    res["calibration"][name] = {}
    for d in DEFS:
        f = frames[d]
        key = (f.patient + "|" + f.celltype).values
        units = pd.unique(key)
        pick = rng.choice(len(regs[d]), size=min(N_REG_PERM, len(regs[d])), replace=False)
        ps = []
        for perm in range(N_PERM):
            sn = f.stage_num.values.copy()
            for u, flip in zip(units, rng.integers(0, 2, len(units))):
                if flip:
                    idx = np.where(key == u)[0]
                    sn[idx] = 1.0 - sn[idx]
            base = f[["patient", "celltype"]].copy()
            base["stage_num"] = sn
            for j in pick:
                base["y"] = centred[d][:, j]
                r = fn(base)
                if r:
                    ps.append(r[1])
        ps = np.asarray(ps)
        res["calibration"][name][d] = {
            "n_fits": int(len(ps)),
            "frac_p_below_05": float((ps < 0.05).mean()),
            "frac_p_below_01": float((ps < 0.01).mean()),
            "median_p": float(np.median(ps)),
            "ks_uniform_p": float(stats.kstest(ps, "uniform").pvalue),
        }
        c = res["calibration"][name][d]
        log(f"calibration {name} def{d}: p<0.05 in {c['frac_p_below_05']:.4f} "
            f"(nominal 0.05), median p {c['median_p']:.3f}, KS p {c['ks_uniform_p']:.2g}")

# Readability bar: the realised type-I rate at 0.05 must lie in [0.02, 0.10] and the median
# permutation p in [0.40, 0.60], under both vocabularies.
def readable(name):
    return all(0.02 <= res["calibration"][name][d]["frac_p_below_05"] <= 0.10
               and 0.40 <= res["calibration"][name][d]["median_p"] <= 0.60 for d in DEFS)

res["readable"] = {n: readable(n) for n in res["calibration"]}
for n, ok in res["readable"].items():
    log(f"{n}: calibration {'PASSES' if ok else 'FAILS'} the readability bar")

# Fit whichever estimators are readable

res["results"] = {}
for name, fn in (("E1_random_intercept", e1), ("E2_paired_wilcoxon", e2)):
    if not res["readable"][name]:
        log(f"{name}: NOT READ (failed calibration) — no hits are counted from it")
        res["results"][name] = {"read": False}
        continue
    rows = []
    for d in DEFS:
        f = frames[d]
        base = f[["patient", "celltype", "stage_num"]].copy()
        for j, reg in enumerate(regs[d]):
            base["y"] = centred[d][:, j]
            r = fn(base)
            if r:
                rows.append(dict(defs=d, regulon=reg, beta=r[0], p=r[1]))
        log(f"{name} def{d}: {sum(x['defs'] == d for x in rows)} fitted")
    fits = pd.DataFrame(rows)
    fits["q"] = np.nan
    for d in DEFS:
        m = fits.defs == d
        fits.loc[m, "q"] = bh(fits.loc[m, "p"].values)
    fits.to_csv(OUT / f"{name.lower()}_fits.csv", index=False)

    w = fits.pivot(index="regulon", columns="defs", values=["beta", "q"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.dropna()
    w["sign_agree"] = np.sign(w.beta_Diagnosis) == np.sign(w.beta_Relapse)
    w["confirmed"] = w.sign_agree & (w.q_Diagnosis < Q_CUT) & (w.q_Relapse < Q_CUT)
    w.sort_values("q_Diagnosis").to_csv(OUT / f"{name.lower()}_reconciliation.csv")
    res["results"][name] = {
        "read": True,
        "per_vocabulary": {d: {"fitted": int((fits.defs == d).sum()),
                               "q05": int(((fits.defs == d) & (fits.q < Q_CUT)).sum()),
                               "min_q": float(fits[fits.defs == d].q.min())} for d in DEFS},
        "n_shared": int(len(w)),
        "sign_agreement": float(w.sign_agree.mean()),
        "confirmed": int(w.confirmed.sum()),
        "confirmed_regulons": sorted(w[w.confirmed].index.tolist()),
    }
    for d in DEFS:
        x = res["results"][name]["per_vocabulary"][d]
        log(f"{name} def{d}: {x['q05']} at q<{Q_CUT} (min q {x['min_q']:.4g})")
    log(f"{name}: {int(w.confirmed.sum())} CONFIRMED under the cross-definition agreement rule")

with open(OUT / "s3_recalibrated.json", "w") as fh:
    json.dump(res, fh, indent=2)
log("wrote s3_recalibrated.json")
