"""The bystander comparison: which Dx->Rel movement is malignant-specific?"""
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/healthy_published')
REPO = ROOT.parents[1]
OUT = ROOT / "out"
MAL = REPO / "longitudinal/published" / "out"

SEED = 555
MIN_CELLS_PER_UNIT = 20
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
N_DRAWS = 2000
MAL_ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
rng = np.random.default_rng(SEED)

tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[s3b-sel] {m}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def load_arm(pb_path, paired_path, admitted, defs):
    """-> (frame restricted to the paired admitted design, offset-centred matrix, regulon list)"""
    t = pd.read_csv(pb_path)
    pu = set(map(tuple, pd.read_csv(paired_path).values))
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(admitted).values
    f = t[m].copy().reset_index(drop=True)
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    regs = [c for c in f.columns if c not in META + ["stage_num"]]
    A = f[regs].values
    return f, A - A.mean(axis=1, keepdims=True), regs

def per_patient_delta(f, col):
    """within-patient Rel-Dx, averaged over that patient's cell types -> Series indexed by patient"""
    d = f[["patient", "celltype", "stage_num"]].copy()
    d["y"] = col
    p = d.pivot_table(index=["patient", "celltype"], columns="stage_num", values="y")
    if 0.0 not in p or 1.0 not in p:
        return None
    p = p.dropna()
    return (p[1.0] - p[0.0]).groupby(level="patient").mean()

# Load both arms

arms = {}
for d in DEFS:
    hf, hc, hr = load_arm(OUT / f"pseudobulk_healthy_def{d}.csv",
                          OUT / "paired_units_healthy.csv",
                          sorted(pd.read_csv(OUT / "paired_units_healthy.csv")
                                 .celltype.value_counts().pipe(lambda s: s[s >= 5]).index), d)
    mf, mc, mr = load_arm(MAL / f"pseudobulk_def{d}.csv", MAL / "paired_units.csv",
                          MAL_ADMITTED, d)
    assert hr == mr, f"vocabulary mismatch between arms for def{d}"
    arms[d] = dict(hf=hf, hc=hc, mf=mf, mc=mc, regs=hr)
    log(f"def{d}: healthy {len(hf)} rows / {hf.patient.nunique()} patients | "
        f"malignant {len(mf)} rows / {mf.patient.nunique()} patients | {len(hr)} regulons")

shared_patients = sorted(set(arms["Diagnosis"]["hf"].patient) &
                         set(arms["Diagnosis"]["mf"].patient))
log(f"patients contributing to BOTH compartments: {len(shared_patients)} "
    f"(healthy-only {len(set(arms['Diagnosis']['hf'].patient) - set(shared_patients))}, "
    f"malignant-only {len(set(arms['Diagnosis']['mf'].patient) - set(shared_patients))})")
if len(shared_patients) < 5:
    raise SystemExit("[s3b-sel] ABORT — fewer than 5 patients in both compartments; "
                     "the paired DiD is not estimable")

# The DiD, per regulon per vocabulary

rows = []
for d in DEFS:
    a = arms[d]
    for j, reg in enumerate(a["regs"]):
        dh = per_patient_delta(a["hf"], a["hc"][:, j])
        dm = per_patient_delta(a["mf"], a["mc"][:, j])
        if dh is None or dm is None:
            continue
        idx = [p for p in shared_patients if p in dh.index and p in dm.index]
        if len(idx) < 5:
            continue
        v = dm.loc[idx].values - dh.loc[idx].values
        if np.allclose(v, 0):
            continue
        rows.append(dict(
            defs=d, regulon=reg, TF=tf_of(reg), n_patients=len(idx),
            beta_malignant=float(dm.loc[idx].mean()), beta_healthy=float(dh.loc[idx].mean()),
            did=float(np.mean(v)), did_median=float(np.median(v)),
            p_did=float(stats.wilcoxon(v, zero_method="wilcox").pvalue),
            frac_patients_same_sign=float(max((v > 0).mean(), (v < 0).mean())),
            sd_did=float(np.std(v, ddof=1))))
    log(f"def{d}: {sum(r['defs'] == d for r in rows)} regulons with an estimable DiD")

did = pd.DataFrame(rows)
did["q_did"] = np.nan
for d in DEFS:
    m = did.defs == d
    did.loc[m, "q_did"] = bh(did.loc[m, "p_did"].values)
did.to_csv(OUT / "did_fits.csv", index=False)

res = {"seed": SEED,
       "design": {"shared_patients": len(shared_patients),
                  "healthy_admitted": sorted(pd.read_csv(OUT / "paired_units_healthy.csv")
                                             .celltype.value_counts()
                                             .pipe(lambda s: s[s >= 5]).index),
                  "malignant_admitted": MAL_ADMITTED,
                  "admitted_overlap": []},
       "caveat": "compartment-level bystander control; the two arms share NO cell types"}

# Aggregate: do the two compartments move together at all?

res["cross_compartment_agreement"] = {}
for d in DEFS:
    g = did[did.defs == d]
    bm, bh_ = g.beta_malignant.values, g.beta_healthy.values
    pear = stats.pearsonr(bm, bh_)
    spear = stats.spearmanr(bm, bh_)
    obs_sign = float(np.mean(np.sign(bm) == np.sign(bh_)))
    draws = np.empty(N_DRAWS)
    r2 = np.random.default_rng(SEED)
    for i in range(N_DRAWS):
        draws[i] = np.mean(np.sign(bm) == np.sign(r2.permutation(bh_)))
    res["cross_compartment_agreement"][d] = {
        "n_regulons": int(len(g)),
        "pearson_r": float(pear[0]), "pearson_p": float(pear[1]),
        "spearman_r": float(spear[0]), "spearman_p": float(spear[1]),
        "sign_concordance": obs_sign,
        "sign_null_mean": float(draws.mean()),
        "sign_null_ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        "sign_p_empirical": float((np.sum(draws >= obs_sign) + 1) / (N_DRAWS + 1)),
        "median_abs_beta_malignant": float(np.median(np.abs(bm))),
        "median_abs_beta_healthy": float(np.median(np.abs(bh_))),
        "magnitude_ratio": float(np.median(np.abs(bm)) / np.median(np.abs(bh_))),
        # how much of the DiD is simply the malignant arm? (circularity audit)
        "corr_did_vs_beta_malignant": float(np.corrcoef(g.did.values, bm)[0, 1]),
        "corr_did_vs_beta_healthy": float(np.corrcoef(g.did.values, bh_)[0, 1]),
    }
    x = res["cross_compartment_agreement"][d]
    log(f"S4 def{d}: malignant-vs-healthy beta  Pearson r={x['pearson_r']:+.3f} "
        f"(p={x['pearson_p']:.3g}), Spearman {x['spearman_r']:+.3f}; sign concordance "
        f"{x['sign_concordance']:.3f} vs null {x['sign_null_mean']:.3f} "
        f"{x['sign_null_ci95']}, p={x['sign_p_empirical']:.4g}")
    log(f"      median |beta|: malignant {x['median_abs_beta_malignant']:.5f} | "
        f"healthy {x['median_abs_beta_healthy']:.5f}")

# Reconcile the two vocabularies and assign verdicts

w = did.pivot(index="regulon", columns="defs",
              values=["beta_malignant", "beta_healthy", "did", "p_did", "q_did",
                      "frac_patients_same_sign"])
w.columns = [f"{a}_{b}" for a, b in w.columns]
w = w.dropna(subset=["did_Diagnosis", "did_Relapse"])
w["TF"] = [tf_of(r) for r in w.index]
w["did_sign_agree"] = np.sign(w.did_Diagnosis) == np.sign(w.did_Relapse)
w["did_both_sig"] = (w.q_did_Diagnosis < Q_CUT) & (w.q_did_Relapse < Q_CUT)
w["compartments_agree"] = (
    (np.sign(w.beta_malignant_Diagnosis) == np.sign(w.beta_healthy_Diagnosis)) &
    (np.sign(w.beta_malignant_Relapse) == np.sign(w.beta_healthy_Relapse)))
w["verdict"] = np.where(w.compartments_agree, "SHARED",
                 np.where(w.did_sign_agree & w.did_both_sig, "MALIGNANT-SPECIFIC",
                   np.where(w.did_sign_agree, "DIRECTIONAL-SELECTIVE", "NOT CONFIRMED")))
w["did_min_abs"] = np.minimum(w.did_Diagnosis.abs(), w.did_Relapse.abs()) * \
                   np.sign(w.did_Diagnosis)
w = w.sort_values(["verdict", "q_did_Diagnosis"])
w.to_csv(OUT / "did_reconciliation.csv")

vc = w.verdict.value_counts().to_dict()
log("S3 verdicts over shared-name regulons: " + ", ".join(f"{k} {v}" for k, v in vc.items()))
res["verdicts"] = vc
res["n_shared_name"] = int(len(w))

# DiD sign agreement across vocabularies, against the mismatched-name null
dD, dR = w.did_Diagnosis.values, w.did_Relapse.values
obs = float(np.mean(np.sign(dD) == np.sign(dR)))
draws = np.empty(N_DRAWS)
r2 = np.random.default_rng(SEED)
for i in range(N_DRAWS):
    draws[i] = np.mean(np.sign(dD) == np.sign(r2.permutation(dR)))
res["did_cross_vocabulary"] = {
    "sign_agreement": obs, "null_mean": float(draws.mean()),
    "null_ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
    "p_empirical": float((np.sum(draws >= obs) + 1) / (N_DRAWS + 1))}
log(f"S2 DiD cross-vocabulary sign agreement {obs:.3f} vs null "
    f"{res['did_cross_vocabulary']['null_mean']:.3f}, p="
    f"{res['did_cross_vocabulary']['p_empirical']:.4g}")

# Minimum detectable DiD, so "healthy did not move" carries its own floor

mde = {}
for d in DEFS:
    g = did[did.defs == d]
    sd = float(np.median(g.sd_did))
    n = int(np.median(g.n_patients))
    # two-sided 0.05, 80% power, paired t approximation (Wilcoxon ARE ~0.955 vs t)
    mde[d] = {"median_sd_did": sd, "n_patients": n,
              "mde_80pct_power": float((1.96 + 0.84) * sd / np.sqrt(n) / 0.955),
              "median_abs_observed_did": float(np.median(g.did.abs()))}
    log(f"S5 def{d}: median SD(DiD)={sd:.5f}, n={n} -> minimum detectable DiD "
        f"{mde[d]['mde_80pct_power']:.5f}; median |observed DiD| "
        f"{mde[d]['median_abs_observed_did']:.5f}")
res["minimum_detectable_did"] = mde

# The malignant candidate lists, re-read through the bystander control

cand = json.load(open(MAL / "s3_candidates.json"))
tq = cand["top_quartile_list"]
cand_regs = [r for r in tq.get("up", []) + tq.get("down", []) if r in w.index]
sub = w.loc[cand_regs]
res["malignant_top_quartile_candidates"] = {
    "n_in_list": len(tq.get("up", [])) + len(tq.get("down", [])),
    "n_with_did": int(len(sub)),
    "verdicts": sub.verdict.value_counts().to_dict(),
    "n_compartments_agree_SHARED": int(sub.compartments_agree.sum()),
    "n_did_sign_agree": int(sub.did_sign_agree.sum()),
    "healthy_median_abs_beta": float(sub.beta_healthy_Diagnosis.abs().median()),
    "malignant_median_abs_beta": float(sub.beta_malignant_Diagnosis.abs().median()),
}

# per-vocabulary healthy-vs-malignant sign concordance, the non-circular selectivity read
res["healthy_vs_malignant_concordance"] = {}
for label, g in (("all_shared", w), ("malignant_top_quartile", sub)):
    res["healthy_vs_malignant_concordance"][label] = {
        "n": int(len(g)),
        **{f"concordance_{d}": float(np.mean(np.sign(g[f"beta_malignant_{d}"]) ==
                                             np.sign(g[f"beta_healthy_{d}"]))) for d in DEFS},
        "both_vocabularies": float(g.compartments_agree.mean()),
        "both_vocabularies_independent_chance": 0.25,
        "chance_per_vocabulary": 0.50}
    x = res["healthy_vs_malignant_concordance"][label]
    log(f"concordance [{label}] n={x['n']}: healthy moves WITH malignant in "
        f"{x['concordance_Diagnosis']:.3f} (Dx vocab) / {x['concordance_Relapse']:.3f} "
        f"(Rel vocab) of regulons vs chance 0.50; both vocabularies "
        f"{x['both_vocabularies']:.3f} vs 0.25")
sub.to_csv(OUT / "did_malignant_candidates.csv")
log("candidates: " + ", ".join(f"{k} {v}" for k, v in
                               res["malignant_top_quartile_candidates"]["verdicts"].items()))

with open(OUT / "s3b_selectivity.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote s3b_selectivity.json")
