"""The directional candidate list, with its error rate estimated."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/published')
OUT = ROOT / "out"
S2A = ROOT.parents[1] / "analysis/network_comparison" / "out"

CHANCE = 0.50               # sign agreement expected by chance: mismatched-name null after offset centring (0.512 measured)
A_SWEEP = [1.00, 0.95, 0.90, 0.85]
MIN_CELLS_PER_UNIT = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]

def log(m):
    print(f"[s3-cand] {m}", flush=True)

def precision_band(agreement):
    """Probability an agreeing regulon's shared sign is the true one, swept over `a`."""
    out = {}
    for a in A_SWEEP:
        p = (agreement - CHANCE) / (a - CHANCE)
        p = min(max(p, 0.0), 1.0)
        out[a] = float(1 - (1 - p) * CHANCE / agreement) if agreement > 0 else None
    return out

cat = pd.read_csv(S2A / "regulon_catalogue.csv")
tier = {r: t.split()[0] for r, t in zip(cat.regulon, cat.tier)}

w = pd.read_csv(OUT / "c4_offset_centred_reconciliation.csv", index_col=0)
w["TF"] = [re.split(r"[_(]", str(i))[0] for i in w.index]
w["sign_mode"] = [re.search(r"_([+-]/[+-])$", str(i)).group(1)
                  if re.search(r"_([+-]/[+-])$", str(i)) else "?" for i in w.index]
w["tier"] = [tier.get(i, "none") for i in w.index]
w["minabs"] = np.minimum(w.beta_Diagnosis.abs(), w.beta_Relapse.abs())
w["direction"] = np.where(w.beta_Diagnosis > 0, "UP at relapse", "DOWN at relapse")

# per-patient consistency, computed on the same offset-centred activity the betas came from
cons = {}
for defs in ("Diagnosis", "Relapse"):
    t = pd.read_csv(OUT / f"pseudobulk_def{defs}.csv")
    pu = set(map(tuple, pd.read_csv(OUT / "paired_units.csv").values))
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(ADMITTED).values
    f = t[m].copy()
    regs = [c for c in f.columns if c not in META]
    A = f[regs].values
    Cn = A - A.mean(axis=1, keepdims=True)
    key = pd.MultiIndex.from_arrays([f.patient.values, f.celltype.values])
    isrel = (f.stage == "Relapse").values
    for j, reg in enumerate(regs):
        s = pd.Series(Cn[:, j], index=key)
        dl = s[isrel].groupby(level=[0, 1]).mean() - s[~isrel].groupby(level=[0, 1]).mean()
        per_pat = dl.groupby(level=0).mean().dropna()
        if len(per_pat):
            cons.setdefault(reg, {})[defs] = per_pat

for defs in ("Diagnosis", "Relapse"):
    w[f"frac_pat_{defs}"] = [
        float(max((cons[r][defs] > 0).mean(), (cons[r][defs] < 0).mean()))
        if r in cons and defs in cons[r] else np.nan for r in w.index]
w["n_patients"] = [len(cons[r]["Diagnosis"]) if r in cons and "Diagnosis" in cons[r] else np.nan
                   for r in w.index]

res = {"chance_agreement": CHANCE, "a_sweep": A_SWEEP,
       "status": "EXPLORATORY / hypothesis-generating — post-hoc in ordering",
       "strata": {}}

w["q_bin"] = pd.qcut(w.minabs, 4, labels=["q1", "q2", "q3", "q4"])
STRATA = {
    "all shared": w,
    "tier 1": w[w.tier == "T1"],
    "top |beta| quartile": w[w.q_bin == "q4"],
    "tier 1 AND top quartile": w[(w.tier == "T1") & (w.q_bin == "q4")],
}
for name, sub in STRATA.items():
    if not len(sub):
        continue
    ag = float(sub.sign_agree.mean())
    band = precision_band(ag)
    res["strata"][name] = {
        "n_regulons": int(len(sub)), "agreement": ag,
        "n_agreeing": int(sub.sign_agree.sum()),
        "precision_band": band,
        "implied_fdr_band": {str(k): (None if v is None else float(1 - v))
                             for k, v in band.items()},
    }
    lo, hi = min(band.values()), max(band.values())
    log(f"{name:26s} n={len(sub):3d} agreement={ag:.3f} agreeing={int(sub.sign_agree.sum()):3d} "
        f"-> precision {lo:.2f}-{hi:.2f} (FDR {1 - hi:.2f}-{1 - lo:.2f})")

# ---- the list itself: agreeing regulons, ranked by the conservative effect size ----
cand = w[w.sign_agree].copy()
cand = cand.reindex(cand.minabs.sort_values(ascending=False).index)
cols = ["TF", "sign_mode", "tier", "direction", "beta_Diagnosis", "beta_Relapse", "minabs",
        "q_Diagnosis", "q_Relapse", "frac_pat_Diagnosis", "frac_pat_Relapse", "n_patients",
        "q_bin"]
cand[cols].to_csv(OUT / "directional_candidates.csv")
log(f"wrote directional_candidates.csv — {len(cand)} sign-agreeing regulons")

top = cand[cand.q_bin == "q4"]
res["top_quartile_list"] = {
    "n": int(len(top)),
    "up": sorted(top[top.direction.str.startswith("UP")].index.tolist()),
    "down": sorted(top[top.direction.str.startswith("DOWN")].index.tolist()),
    "TFs_up": sorted(set(top[top.direction.str.startswith("UP")].TF)),
    "TFs_down": sorted(set(top[top.direction.str.startswith("DOWN")].TF)),
}
res["tf_appears_both_directions"] = sorted(
    set(res["top_quartile_list"]["TFs_up"]) & set(res["top_quartile_list"]["TFs_down"]))
log(f"top quartile & agreeing: {len(top)} regulons | "
    f"{len(res['top_quartile_list']['up'])} up, {len(res['top_quartile_list']['down'])} down | "
    f"{len(res['top_quartile_list']['TFs_up'])} TFs up, "
    f"{len(res['top_quartile_list']['TFs_down'])} TFs down")
if res["tf_appears_both_directions"]:
    log(f"  TFs appearing in BOTH directions (different sign modes): "
        f"{', '.join(res['tf_appears_both_directions'])}")

with open(OUT / "s3_candidates.json", "w") as fh:
    json.dump(res, fh, indent=2)
log("wrote s3_candidates.json")
