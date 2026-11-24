"""Control C4: is the Dx->Rel signal a global per-unit offset or regulon-specific change?"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/published')
OUT = ROOT / "out"

MIN_CELLS_PER_UNIT = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
SEED = 555

def log(m):
    print(f"[s3-offset] {m}", flush=True)

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

def fit_stage(d):
    for re_f, kind in (("~stage_num", "slope"), (None, "int")):
        try:
            r = smf.mixedlm("y ~ stage_num + C(celltype)", d, groups=d.patient,
                            re_formula=re_f).fit(reml=True, method="lbfgs", maxiter=200)
            b, se = r.params.get("stage_num"), r.bse.get("stage_num")
            if b is None or not np.isfinite(b) or not np.isfinite(se) or se == 0:
                continue
            if kind == "slope" and not r.converged:
                continue
            return float(b), float(se), float(r.pvalues["stage_num"]), kind
        except Exception:
            continue
    return None

res = {"seed": SEED, "note": "activity centred within pseudobulk unit across the vocabulary"}
rows, offsets = [], {}

for d in DEFS:
    f = load(d)
    regs = [c for c in f.columns if c not in META + ["stage_num"]]
    A = f[regs].values                       # units x regulons
    unit_mean = A.mean(axis=1)               # the per-unit global offset being removed
    C = A - unit_mean[:, None]

    off = pd.DataFrame({"stage": f.stage.values, "offset": unit_mean})
    offsets[d] = {
        "mean_offset_Diagnosis": float(off[off.stage == "Diagnosis"].offset.mean()),
        "mean_offset_Relapse": float(off[off.stage == "Relapse"].offset.mean()),
        "delta": float(off[off.stage == "Relapse"].offset.mean()
                       - off[off.stage == "Diagnosis"].offset.mean()),
    }
    # is the offset itself a significant stage effect?
    dd = f[["patient", "celltype", "stage_num"]].copy()
    dd["y"] = unit_mean
    r = fit_stage(dd)
    offsets[d]["offset_stage_beta"] = r[0] if r else None
    offsets[d]["offset_stage_p"] = r[2] if r else None
    log(f"def{d}: per-unit global offset Dx {offsets[d]['mean_offset_Diagnosis']:.4f} -> "
        f"Rel {offsets[d]['mean_offset_Relapse']:.4f} "
        f"(delta {offsets[d]['delta']:+.4f}, p={offsets[d]['offset_stage_p']})")

    for i, reg in enumerate(regs):
        dd = f[["patient", "celltype", "stage_num"]].copy()
        dd["y"] = C[:, i]
        r = fit_stage(dd)
        if r:
            rows.append(dict(defs=d, regulon=reg, beta=r[0], se=r[1], p=r[2], re_kind=r[3]))
        if (i + 1) % 200 == 0:
            log(f"  def{d}: {i + 1}/{len(regs)} centred fits")

fits = pd.DataFrame(rows)
fits["q"] = np.nan
for d in DEFS:
    m = fits.defs == d
    fits.loc[m, "q"] = bh(fits.loc[m, "p"].values)
fits.to_csv(OUT / "c4_offset_centred_fits.csv", index=False)

w = fits.pivot(index="regulon", columns="defs", values=["beta", "q"])
w.columns = [f"{a}_{b}" for a, b in w.columns]
w = w.dropna()
w["sign_agree"] = np.sign(w.beta_Diagnosis) == np.sign(w.beta_Relapse)
w["confirmed"] = w.sign_agree & (w.q_Diagnosis < Q_CUT) & (w.q_Relapse < Q_CUT)
w.sort_values("q_Diagnosis").to_csv(OUT / "c4_offset_centred_reconciliation.csv")

fD = float((w.beta_Diagnosis < 0).mean())
fR = float((w.beta_Relapse < 0).mean())
res["offsets"] = offsets
res["centred"] = {
    d: {"fitted": int((fits.defs == d).sum()),
        "q05": int(((fits.defs == d) & (fits.q < Q_CUT)).sum()),
        "min_q": float(fits[fits.defs == d].q.min())} for d in DEFS}
res["centred_reconciliation"] = {
    "n_shared": int(len(w)),
    "sign_agreement": float(w.sign_agree.mean()),
    "marginal_frac_negative": {"Diagnosis_defs": fD, "Relapse_defs": fR},
    "analytic_chance_from_marginals": fD * fR + (1 - fD) * (1 - fR),
    "confirmed": int(w.confirmed.sum()),
    "confirmed_regulons": sorted(w[w.confirmed].index.tolist()),
}
for d in DEFS:
    log(f"C4 def{d}: {res['centred'][d]['q05']} regulons at q<{Q_CUT} "
        f"(min q {res['centred'][d]['min_q']:.4g})")
log(f"C4 reconciliation: sign agreement {w.sign_agree.mean():.3f} "
    f"(marginal chance {res['centred_reconciliation']['analytic_chance_from_marginals']:.3f}), "
    f"{int(w.confirmed.sum())} CONFIRMED after centring")

with open(OUT / "s3_offset.json", "w") as fh:
    json.dump(res, fh, indent=2)
log("wrote s3_offset.json")
