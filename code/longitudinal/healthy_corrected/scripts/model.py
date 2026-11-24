"""Healthy-compartment Dx->Rel contrast, fitted with the same models as the malignant-cell analysis."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/healthy_corrected')
OUT = ROOT / "out"

SEED = 555
MIN_CELLS_PER_UNIT = 20
MIN_PATIENTS_PER_CELLTYPE = 5
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
N_PERM = 20                        # as the malignant calibration
N_REG_PERM = 100
N_DRAWS_C1 = 2000
rng = np.random.default_rng(SEED)

def log(m):
    print(f"[s3b-model] {m}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

# Design

paired = pd.read_csv(OUT / "paired_units_healthy.csv")
pu = set(map(tuple, paired.values))
ct_n = paired.celltype.value_counts()
ADMITTED = sorted(ct_n[ct_n >= MIN_PATIENTS_PER_CELLTYPE].index)
log(f"H3 admitted healthy cell types ({len(ADMITTED)}): " +
    ", ".join(f"{c} ({ct_n[c]})" for c in ADMITTED))
log(f"H3 excluded ({int((ct_n < MIN_PATIENTS_PER_CELLTYPE).sum())}): " +
    ", ".join(f"{c} ({v})" for c, v in ct_n[ct_n < MIN_PATIENTS_PER_CELLTYPE].items()))

def load(defs):
    t = pd.read_csv(OUT / f"pseudobulk_healthy_def{defs}.csv")
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(ADMITTED).values
    f = t[m].copy()
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    return f.reset_index(drop=True)

def e1(d):
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

frames, centred, raw, regs = {}, {}, {}, {}
for d in DEFS:
    f = load(d)
    r = [c for c in f.columns if c not in META + ["stage_num"]]
    A = f[r].values
    frames[d], regs[d] = f, r
    raw[d] = A
    centred[d] = A - A.mean(axis=1, keepdims=True)          # H4
    log(f"def{d}: {len(f)} rows = {f.patient.nunique()} patients x {f.celltype.nunique()} "
        f"cell types x 2 stages, {len(r)} regulons")

res = {"seed": SEED,
       "arm": "HEALTHY bystander compartment, scored under the MALIGNANT vocabularies",
       "design": {d: {"rows": int(len(frames[d])), "patients": int(frames[d].patient.nunique()),
                      "celltypes": int(frames[d].celltype.nunique()),
                      "regulons": len(regs[d])} for d in DEFS},
       "admitted_celltypes": ADMITTED,
       "admitted_overlap_with_malignant": sorted(
           set(ADMITTED) & {"Early.Basophil", "Early.Erythrocyte", "GMP", "HSC",
                            "Monocytes", "Progenitor"})}

# Calibration (run before the fits; the healthy data must pass its own calibration)

res["calibration"] = {}
for name, fn in (("E1_random_intercept", e1), ("E2_paired_wilcoxon", e2)):
    res["calibration"][name] = {}
    for d in DEFS:
        f = frames[d]
        key = (f.patient + "|" + f.celltype).values
        units = pd.unique(key)
        pick = rng.choice(len(regs[d]), size=min(N_REG_PERM, len(regs[d])), replace=False)
        ps = []
        for _ in range(N_PERM):
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
            "n_fits": int(len(ps)), "frac_p_below_05": float((ps < 0.05).mean()),
            "frac_p_below_01": float((ps < 0.01).mean()), "median_p": float(np.median(ps)),
            "ks_uniform_p": float(stats.kstest(ps, "uniform").pvalue)}
        c = res["calibration"][name][d]
        log(f"calibration {name} def{d}: p<0.05 in {c['frac_p_below_05']:.4f} "
            f"(nominal 0.05), median p {c['median_p']:.3f}, KS p {c['ks_uniform_p']:.2g}")

def readable(name):                                          # H2
    return all(0.02 <= res["calibration"][name][d]["frac_p_below_05"] <= 0.10
               and 0.40 <= res["calibration"][name][d]["median_p"] <= 0.60 for d in DEFS)

res["readable"] = {n: readable(n) for n in res["calibration"]}
for n, ok in res["readable"].items():
    log(f"{n}: calibration {'PASSES' if ok else 'FAILS'} the readability bar")

# Fits — offset-centred (primary read) and raw (drift)

res["results"] = {}
for scale, MAT in (("centred", centred), ("raw", raw)):
    rows = []
    for d in DEFS:
        f = frames[d]
        base = f[["patient", "celltype", "stage_num"]].copy()
        for j, reg in enumerate(regs[d]):
            base["y"] = MAT[d][:, j]
            r = e2(base)
            if r:
                rows.append(dict(defs=d, regulon=reg, beta=r[0], p=r[1]))
        log(f"E2 {scale} def{d}: {sum(x['defs'] == d for x in rows)} fitted")
    fits = pd.DataFrame(rows)
    fits["q"] = np.nan
    for d in DEFS:
        m = fits.defs == d
        fits.loc[m, "q"] = bh(fits.loc[m, "p"].values)
    fits.to_csv(OUT / f"e2_healthy_fits_{scale}.csv", index=False)

    w = fits.pivot(index="regulon", columns="defs", values=["beta", "q"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.dropna()
    w["sign_agree"] = np.sign(w.beta_Diagnosis) == np.sign(w.beta_Relapse)
    w["confirmed"] = w.sign_agree & (w.q_Diagnosis < Q_CUT) & (w.q_Relapse < Q_CUT)
    w.to_csv(OUT / f"e2_healthy_reconciliation_{scale}.csv")

    # mismatched-name null for the agreement statistic
    bD, bR = w.beta_Diagnosis.values, w.beta_Relapse.values
    draws = np.empty(N_DRAWS_C1)
    r2 = np.random.default_rng(SEED)
    for i in range(N_DRAWS_C1):
        draws[i] = np.mean(np.sign(bD) == np.sign(r2.permutation(bR)))
    if scale == "centred":
        np.save(OUT / "c1_healthy_sign_null.npy", draws)
    obs = float(w.sign_agree.mean())

    res["results"][scale] = {
        "per_vocabulary": {d: {"fitted": int((fits.defs == d).sum()),
                               "q05": int(((fits.defs == d) & (fits.q < Q_CUT)).sum()),
                               "min_q": float(fits.loc[fits.defs == d, "q"].min()),
                               "frac_down": float((fits.loc[fits.defs == d, "beta"] < 0).mean())}
                           for d in DEFS},
        "n_shared": int(len(w)),
        "sign_agreement": obs,
        "sign_null_mean": float(draws.mean()),
        "sign_null_ci95": [float(np.percentile(draws, 2.5)),
                           float(np.percentile(draws, 97.5))],
        "sign_p_empirical": float((np.sum(draws >= obs) + 1) / (N_DRAWS_C1 + 1)),
        "confirmed": int(w.confirmed.sum()),
        "confirmed_regulons": list(w.index[w.confirmed]),
    }
    x = res["results"][scale]
    log(f"E2 {scale}: agreement {x['sign_agreement']:.3f} vs null "
        f"{x['sign_null_mean']:.3f} {x['sign_null_ci95']}, empirical p {x['sign_p_empirical']:.4g}; "
        f"confirmed {x['confirmed']}")
    log(f"E2 {scale} drift: frac moving DOWN = " +
        ", ".join(f"def{d} {x['per_vocabulary'][d]['frac_down']:.3f}" for d in DEFS))

with open(OUT / "s3b_recalibrated.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote s3b_recalibrated.json")
