"""The within-patient Dx->Rel mixed model, fitted twice and reconciled."""
import json
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/published')
OUT = ROOT / "out"

MIN_CELLS_PER_UNIT = 20
MIN_PATIENTS_PER_CELLTYPE = 5  # minimum paired patients for a cell type to be admitted
Q_CUT = 0.05                   # FDR cutoff
LOGIT_TRIGGER_FRAC = 0.01      # use the logit scale if more than this fraction of values lie at the edges
LOGIT_EDGE = 1e-3              # distance from 0 or 1 that counts as an edge value
DEFS = ("Diagnosis", "Relapse")

tf_of = lambda s: re.split(r"[_(]", str(s))[0]
META = ["stage", "patient", "celltype", "driver", "n_cells"]

# A singular random-effects covariance means the random-slope variance collapsed to ~0, i.e.
# the stage effect is essentially identical across patients. That is a result, not an error,
# and it is captured per regulon in `slope_var` / `icc_slope`. Only the console warning is silenced.
warnings.filterwarnings("ignore", message="Random effects covariance is singular")

def log(m):
    print(f"[s3-model] {m}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def driver_class(d):
    if d in ("UNLABELLED", "AMBIGUOUS") or not isinstance(d, str):
        return "unlabelled"
    if d.startswith("KMT2A"):
        return "KMT2A-r"
    if d in ("RUNX1-RUNX1T1", "CBFB-MYH11"):
        return "CBF-r"
    if d == "FLT3-ITD":
        return "FLT3-ITD"
    return "other"

# Load + restrict to the paired, admitted design

tabs, regs = {}, {}
for d in DEFS:
    t = pd.read_csv(OUT / f"pseudobulk_def{d}.csv")
    tabs[d] = t
    regs[d] = [c for c in t.columns if c not in META]

paired = pd.read_csv(OUT / "paired_units.csv")
pu = set(map(tuple, paired.values))

ct_n = paired.celltype.value_counts()
admitted = sorted(ct_n[ct_n >= MIN_PATIENTS_PER_CELLTYPE].index)
dropped_ct = sorted(ct_n[ct_n < MIN_PATIENTS_PER_CELLTYPE].index)
log(f"M4 admitted cell types ({len(admitted)}): " +
    ", ".join(f"{c} ({ct_n[c]})" for c in admitted))
log(f"M4 excluded cell types ({len(dropped_ct)}): " +
    ", ".join(f"{c} ({ct_n[c]})" for c in dropped_ct))

frames = {}
for d in DEFS:
    t = tabs[d]
    m = (t.n_cells >= MIN_CELLS_PER_UNIT) & \
        [tuple(x) in pu for x in t[["patient", "celltype"]].values] & \
        t.celltype.isin(admitted)
    f = t[m].copy()
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    f["dclass"] = f.driver.map(driver_class)
    frames[d] = f
    log(f"def{d}: {len(f)} rows = {f.patient.nunique()} patients x "
        f"{f.celltype.nunique()} cell types x 2 stages, {len(regs[d])} regulons")

# Outcome-scale trigger (logit vs raw), evaluated once on the admitted rows

edge_frac = {}
for d in DEFS:
    V = frames[d][regs[d]].values
    edge_frac[d] = float(np.mean((V <= LOGIT_EDGE) | (V >= 1 - LOGIT_EDGE)))
USE_LOGIT = max(edge_frac.values()) > LOGIT_TRIGGER_FRAC
log(f"M5 edge fraction: " + ", ".join(f"def{d}={edge_frac[d]:.4f}" for d in DEFS) +
    f"  -> trigger {LOGIT_TRIGGER_FRAC} => scale = {'logit' if USE_LOGIT else 'raw'}")

def transform(y, n):
    if not USE_LOGIT:
        return y
    off = 1.0 / (2.0 * np.maximum(n, 1))
    p = np.clip(y, off, 1 - off)
    return np.log(p / (1 - p))

# The primary fit

def fit_one(df, reg):
    """MixedLM `y ~ stage + C(celltype)` with (stage|patient); fall back to (1|patient)."""
    d = df[["patient", "celltype", "stage_num", "n_cells"]].copy()
    d["y"] = transform(df[reg].values, df.n_cells.values)
    if not np.isfinite(d.y).all() or d.y.std() == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        for re_f, kind in (("~stage_num", "slope"), (None, "intercept")):
            try:
                md = smf.mixedlm("y ~ stage_num + C(celltype)", d, groups=d.patient,
                                 re_formula=re_f)
                r = md.fit(reml=True, method="lbfgs", maxiter=200)
                b, se = r.params.get("stage_num"), r.bse.get("stage_num")
                if b is None or not np.isfinite(b) or not np.isfinite(se) or se == 0:
                    continue
                # A non-converged random-slope fit is not accepted; fall back to (1|patient)
                # rather than report a coefficient the optimiser did not settle on.
                if kind == "slope" and not r.converged:
                    continue
                slope_var = float(r.cov_re.iloc[1, 1]) if kind == "slope" else np.nan
                return dict(beta=float(b), se=float(se), p=float(r.pvalues["stage_num"]),
                            re_kind=kind, converged=bool(r.converged),
                            slope_var=slope_var, resid_var=float(r.scale),
                            icc_slope=float(slope_var / (slope_var + r.scale))
                            if kind == "slope" and np.isfinite(slope_var) else np.nan)
            except Exception:
                continue
    return None

FIT_CSV = OUT / "lmm_fits_by_definition.csv"
REUSE = os.environ.get("S3_REUSE_FITS") == "1" and FIT_CSV.exists()
if REUSE:
    fits = pd.read_csv(FIT_CSV)
    log(f"reusing {FIT_CSV.name} ({len(fits)} fits) — S3_REUSE_FITS=1")
else:
    rows = []
    for d in DEFS:
        f = frames[d]
        for i, reg in enumerate(regs[d]):
            r = fit_one(f, reg)
            if r is None:
                continue
            r.update(defs=d, regulon=reg, TF=tf_of(reg))
            rows.append(r)
            if (i + 1) % 200 == 0:
                log(f"  def{d}: {i + 1}/{len(regs[d])} fitted")
        log(f"def{d}: {sum(x['defs'] == d for x in rows)} regulons fitted")

    fits = pd.DataFrame(rows)
    fits["q"] = np.nan
    for d in DEFS:
        m = fits.defs == d
        fits.loc[m, "q"] = bh(fits.loc[m, "p"].values)
    fits.to_csv(FIT_CSV, index=False)

# Sensitivity analysis: precision-weighted within-patient paired delta

def paired_delta(df, reg):
    d = df[["patient", "celltype", "stage", "n_cells"]].copy()
    d["y"] = transform(df[reg].values, df.n_cells.values)
    w = d.pivot_table(index=["patient", "celltype"], columns="stage",
                      values=["y", "n_cells"])
    if ("y", "Diagnosis") not in w or ("y", "Relapse") not in w:
        return None
    w = w.dropna()
    if not len(w):
        return None
    delta = w[("y", "Relapse")] - w[("y", "Diagnosis")]
    nh = 2 / (1 / w[("n_cells", "Diagnosis")] + 1 / w[("n_cells", "Relapse")])
    pw = pd.DataFrame({"d": delta.values, "w": nh.values},
                      index=delta.index).groupby(level="patient").apply(
        lambda g: np.average(g.d, weights=g.w))
    v = pw.values
    if len(v) < 5:
        return None
    return dict(n_patients=int(len(v)), wmean_delta=float(np.mean(v)),
                wmedian_delta=float(np.median(v)),
                p_wilcoxon=float(stats.wilcoxon(v, zero_method="wilcox").pvalue),
                frac_patients_same_sign=float(max((v > 0).mean(), (v < 0).mean())))

srows = []
for d in DEFS:
    f = frames[d]
    for reg in regs[d]:
        r = paired_delta(f, reg)
        if r:
            r.update(defs=d, regulon=reg)
            srows.append(r)
sens = pd.DataFrame(srows)
sens["q_wilcoxon"] = np.nan
for d in DEFS:
    m = sens.defs == d
    sens.loc[m, "q_wilcoxon"] = bh(sens.loc[m, "p_wilcoxon"].values)
sens.to_csv(OUT / "sensitivity_weighted_paired_delta.csv", index=False)
log(f"M6 sensitivity: {len(sens)} regulon x definition tests")

# Reconcile the two definition sets

shared = sorted(set(regs["Diagnosis"]) & set(regs["Relapse"]))
priv = {d: sorted(set(regs[d]) - set(regs["Relapse" if d == "Diagnosis" else "Diagnosis"]))
        for d in DEFS}
log(f"vocabulary: {len(shared)} shared regulon names | "
    f"{len(priv['Diagnosis'])} Dx-private | {len(priv['Relapse'])} Rel-private")

wide = fits.pivot(index="regulon", columns="defs",
                  values=["beta", "se", "p", "q", "re_kind", "converged", "icc_slope"])
wide.columns = [f"{a}_{b}" for a, b in wide.columns]
for c in wide.columns:
    if c.split("_")[0] in ("beta", "se", "p", "q", "icc") or c.startswith("icc_slope"):
        wide[c] = pd.to_numeric(wide[c], errors="coerce")
w = wide.loc[[r for r in shared if r in wide.index]].copy()
w = w.dropna(subset=["beta_Diagnosis", "beta_Relapse", "q_Diagnosis", "q_Relapse"])
w["TF"] = [tf_of(r) for r in w.index]
w["sign_agree"] = np.sign(w.beta_Diagnosis) == np.sign(w.beta_Relapse)
w["both_sig"] = (w.q_Diagnosis < Q_CUT) & (w.q_Relapse < Q_CUT)
w["verdict"] = np.where(w.sign_agree & w.both_sig, "CONFIRMED",
                 np.where(w.sign_agree, "DIRECTIONAL", "NOT CONFIRMED"))
w["beta_min_abs"] = np.minimum(w.beta_Diagnosis.abs(), w.beta_Relapse.abs()) * \
                    np.sign(w.beta_Diagnosis)
sens_w = sens.pivot(index="regulon", columns="defs",
                    values=["wmean_delta", "q_wilcoxon", "frac_patients_same_sign"])
sens_w.columns = [f"{a}_{b}" for a, b in sens_w.columns]
w = w.join(sens_w)
w["sensitivity_concurs"] = (np.sign(w.wmean_delta_Diagnosis) == np.sign(w.beta_Diagnosis)) & \
                           (np.sign(w.wmean_delta_Relapse) == np.sign(w.beta_Relapse))
w = w.sort_values(["verdict", "q_Diagnosis"])
w.to_csv(OUT / "cross_definition_reconciliation.csv")

vc = w.verdict.value_counts().to_dict()
log("M7 verdicts over shared-name regulons: " +
    ", ".join(f"{k} {v}" for k, v in vc.items()))
conf = w[w.verdict == "CONFIRMED"]
log(f"M7 CONFIRMED: {len(conf)} regulons over {conf.TF.nunique()} TFs "
    f"({int((conf.beta_min_abs > 0).sum())} up at relapse, "
    f"{int((conf.beta_min_abs < 0).sum())} down)")

# stage-private regulons: single definition, reported separately
prows = []
for d in DEFS:
    sub = fits[(fits.defs == d) & (fits.regulon.isin(priv[d]))].copy()
    sub["private_to"] = d
    prows.append(sub)
private = pd.concat(prows, ignore_index=True).sort_values("q")
private.to_csv(OUT / "stage_private_regulons.csv", index=False)
log(f"M7 stage-private: {len(private)} fitted, "
    f"{int((private.q < Q_CUT).sum())} at q<{Q_CUT} (single definition — NOT cross-validated)")

# Composition vs within-cell-type activity change: stage x cell-type interaction on the CONFIRMED set

decomp = []
for reg in conf.index:
    for d in DEFS:
        f = frames[d]
        y = transform(f[reg].values, f.n_cells.values)
        dd = f[["patient", "celltype", "stage_num", "n_cells"]].copy()
        dd["y"] = y
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                r = smf.mixedlm("y ~ stage_num * C(celltype)", dd, groups=dd.patient,
                                re_formula="~stage_num").fit(reml=False, method="lbfgs",
                                                             maxiter=200)
                r0 = smf.mixedlm("y ~ stage_num + C(celltype)", dd, groups=dd.patient,
                                 re_formula="~stage_num").fit(reml=False, method="lbfgs",
                                                              maxiter=200)
                k = len(r.params) - len(r0.params)
                lr = 2 * (r.llf - r0.llf)
                pval = float(stats.chi2.sf(lr, k)) if k > 0 and np.isfinite(lr) and lr > 0 else 1.0
            except Exception:
                continue
        # per-cell-type within-stage effect, from the interaction model
        per_ct = {}
        base = r.params.get("stage_num", np.nan)
        for c in admitted:
            key = f"stage_num:C(celltype)[T.{c}]"
            per_ct[c] = float(base + r.params.get(key, 0.0))
        decomp.append(dict(regulon=reg, defs=d, TF=tf_of(reg),
                           lrt_p_stage_x_celltype=pval,
                           n_celltypes_same_sign=int(sum(
                               np.sign(v) == np.sign(conf.loc[reg, "beta_min_abs"])
                               for v in per_ct.values())),
                           n_celltypes=len(per_ct), **{f"beta_{c}": v for c, v in per_ct.items()}))
dec = pd.DataFrame(decomp)
if len(dec):
    dec["q_lrt"] = bh(dec.lrt_p_stage_x_celltype.values)
    dec.to_csv(OUT / "decomposition_stage_by_celltype.csv", index=False)
    log(f"rewiring-vs-composition decomposition: {len(dec)} regulon x definition fits, "
        f"{int((dec.q_lrt < Q_CUT).sum())} with a significant stage:cell_type interaction")

# Driver class: descriptive only (the stage x driver term is not estimable)

ddrows = []
for reg in conf.index:
    for d in DEFS:
        f = frames[d]
        dd = f[["patient", "celltype", "stage", "dclass", "n_cells"]].copy()
        dd["y"] = transform(f[reg].values, f.n_cells.values)
        pv = dd.pivot_table(index=["patient", "celltype", "dclass"], columns="stage", values="y")
        if pv.shape[1] < 2:
            continue
        pv = pv.dropna()
        delta = (pv.Relapse - pv.Diagnosis).groupby(level=["patient", "dclass"]).mean()
        for dc, g in delta.groupby(level="dclass"):
            ddrows.append(dict(regulon=reg, defs=d, driver_class=dc, n_patients=len(g),
                               mean_delta=float(g.mean())))
if ddrows:
    pd.DataFrame(ddrows).to_csv(OUT / "driver_descriptive.csv", index=False)

results = {
    "decisions": {"M1": "stage:driver dropped (not estimable)",
                  "M2": "driver main effect dropped (collinear with patient RE)",
                  "M3": "activity ~ stage + C(celltype), (stage|patient), fallback (1|patient)",
                  "M4": {"min_patients_per_celltype": MIN_PATIENTS_PER_CELLTYPE,
                         "admitted": admitted, "excluded": dropped_ct},
                  "M5": {"trigger_frac": LOGIT_TRIGGER_FRAC, "edge": LOGIT_EDGE,
                         "realised_edge_frac": edge_frac, "scale": "logit" if USE_LOGIT else "raw"},
                  "M6": "precision weights as declared sensitivity (MixedLM takes no weights)",
                  "M7": {"q_cut": Q_CUT},
                  "M8": "BH within definition set"},
    "design": {d: {"rows": int(len(frames[d])), "patients": int(frames[d].patient.nunique()),
                   "celltypes": int(frames[d].celltype.nunique()),
                   "regulons": len(regs[d])} for d in DEFS},
    "vocabulary": {"shared": len(shared), "dx_private": len(priv["Diagnosis"]),
                   "rel_private": len(priv["Relapse"])},
    "fits": {d: {"fitted": int((fits.defs == d).sum()),
                 "converged": int(((fits.defs == d) & fits.converged).sum()),
                 "random_slope": int(((fits.defs == d) & (fits.re_kind == "slope")).sum()),
                 "q05": int(((fits.defs == d) & (fits.q < Q_CUT)).sum())} for d in DEFS},
    "verdicts": vc,
    "confirmed": {"n": int(len(conf)), "n_TF": int(conf.TF.nunique()),
                  "up": int((conf.beta_min_abs > 0).sum()),
                  "down": int((conf.beta_min_abs < 0).sum()),
                  "sensitivity_concurs": int(conf.sensitivity_concurs.fillna(False).sum())},
    "stage_private": {"fitted": int(len(private)), "q05": int((private.q < Q_CUT).sum())},
}
with open(OUT / "s3_results.json", "w") as fh:
    json.dump(results, fh, indent=2, default=str)
log("wrote s3_results.json")
