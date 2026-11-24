"""Within-cell-type activity and cell-type composition effects in the directional regulon set."""
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/published')
REPO = ROOT.parents[1]
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

tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[s3-decomp] {m}", flush=True)

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
    return t[m].copy().reset_index(drop=True)

# The directional set

rec = pd.read_csv(OUT / "e2_paired_wilcoxon_reconciliation.csv", index_col=0)
DIRECTIONAL = list(rec.index[rec.sign_agree])
log(f"D1 directional set: {len(DIRECTIONAL)} of {len(rec)} shared-name regulons "
    f"(sign agreement {rec.sign_agree.mean():.3f} vs a 0.516 mismatched-name null)")

frames, centred, raw, regs = {}, {}, {}, {}
for d in DEFS:
    f = load(d)
    r = [c for c in f.columns if c not in META]
    A = f[r].values
    frames[d], regs[d] = f, r
    raw[d] = A
    centred[d] = A - A.mean(axis=1, keepdims=True)
    log(f"def{d}: {len(f)} rows = {f.patient.nunique()} patients x {f.celltype.nunique()} "
        f"cell types x 2 stages")

col_ix = {d: {reg: j for j, reg in enumerate(regs[d])} for d in DEFS}

# The two statistics

def unit_table(f, y, stage_num):
    """-> tidy frame of one value per (patient, celltype, stage)"""
    return pd.DataFrame({"patient": f.patient.values, "celltype": f.celltype.values,
                         "stage_num": stage_num, "n_cells": f.n_cells.values, "y": y})

def decompose(tab):
    '(A) per-patient naive / within-cell-type activity change / composition, then Wilcoxon over patients.'
    out = []
    for p, g in tab.groupby("patient"):
        dx = g[g.stage_num == 0.0].set_index("celltype")
        rl = g[g.stage_num == 1.0].set_index("celltype")
        cts = dx.index.intersection(rl.index)
        if len(cts) < 2:
            continue
        wdx = dx.loc[cts, "n_cells"].values.astype(float)
        wrl = rl.loc[cts, "n_cells"].values.astype(float)
        wdx, wrl = wdx / wdx.sum(), wrl / wrl.sum()
        mdx, mrl = dx.loc[cts, "y"].values, rl.loc[cts, "y"].values
        naive = float((wrl * mrl).sum() - (wdx * mdx).sum())
        rewire = float((wdx * (mrl - mdx)).sum())
        out.append((naive, rewire, naive - rewire))
    if len(out) < 5:
        return None
    a = np.asarray(out)
    naive, rewire, comp = a[:, 0], a[:, 1], a[:, 2]
    res = dict(n_patients=len(a),
               naive_mean=float(naive.mean()), rewiring_mean=float(rewire.mean()),
               composition_mean=float(comp.mean()),
               p_rewiring=float(stats.wilcoxon(rewire, zero_method="wilcox").pvalue)
               if not np.allclose(rewire, 0) else 1.0,
               p_composition=float(stats.wilcoxon(comp, zero_method="wilcox").pvalue)
               if not np.allclose(comp, 0) else 1.0)
    denom = abs(res["rewiring_mean"]) + abs(res["composition_mean"])
    res["composition_share"] = float(abs(res["composition_mean"]) / denom) if denom > 0 else np.nan
    res["same_sign_naive_vs_rewiring"] = bool(np.sign(res["naive_mean"]) ==
                                              np.sign(res["rewiring_mean"]))
    return res

def heterogeneity(tab):
    """(B) are the per-(patient x cell-type) deltas different across cell types?"""
    p = tab.pivot_table(index=["patient", "celltype"], columns="stage_num", values="y")
    if 0.0 not in p or 1.0 not in p:
        return None
    p = p.dropna()
    delta = (p[1.0] - p[0.0]).rename("d").reset_index()
    groups = [g.d.values for _, g in delta.groupby("celltype") if len(g) >= 3]
    if len(groups) < 3:
        return None
    try:
        kw = stats.kruskal(*groups)
    except ValueError:
        return None
    per_ct = delta.groupby("celltype").d.mean()
    overall = float(delta.groupby("patient").d.mean().mean())
    same = int((np.sign(per_ct.values) == np.sign(overall)).sum()) if overall != 0 else 0
    return dict(p_heterogeneity=float(kw.pvalue), n_celltypes=len(per_ct),
                n_celltypes_same_sign=same,
                **{f"beta_{c}": float(v) for c, v in per_ct.items()})

# Calibration before any p-value is read

res = {"seed": SEED, "status": "EXPLORATORY / post-hoc in ordering",
       "predeclared_primary_stands": "0 regulons CONFIRMED (s3_results.json); this step does not "
                                     "alter that result, it reports the rewiring-vs-composition decomposition on the directional set",
       "directional_set_size": len(DIRECTIONAL),
       "scale": "offset-centred (raw run as declared sensitivity)",
       "admitted_celltypes": ADMITTED,
       "composition_component_is_a_lower_bound": True}

res["calibration"] = {}
for name in ("rewiring", "composition", "heterogeneity"):
    res["calibration"][name] = {}
for d in DEFS:
    f = frames[d]
    key = (f.patient + "|" + f.celltype).values
    units = pd.unique(key)
    pick = rng.choice(len(regs[d]), size=min(N_REG_PERM, len(regs[d])), replace=False)
    ps = {"rewiring": [], "composition": [], "heterogeneity": []}
    sn0 = (f.stage == "Relapse").astype(float).values
    for _ in range(N_PERM):
        sn = sn0.copy()
        for u, flip in zip(units, rng.integers(0, 2, len(units))):
            if flip:
                idx = np.where(key == u)[0]
                sn[idx] = 1.0 - sn[idx]
        for j in pick:
            tab = unit_table(f, centred[d][:, j], sn)
            a = decompose(tab)
            if a:
                ps["rewiring"].append(a["p_rewiring"])
                ps["composition"].append(a["p_composition"])
            b = heterogeneity(tab)
            if b:
                ps["heterogeneity"].append(b["p_heterogeneity"])
    for name, v in ps.items():
        v = np.asarray(v)
        res["calibration"][name][d] = {
            "n_fits": int(len(v)), "frac_p_below_05": float((v < 0.05).mean()),
            "median_p": float(np.median(v)),
            "ks_uniform_p": float(stats.kstest(v, "uniform").pvalue)}
        c = res["calibration"][name][d]
        log(f"calibration {name} def{d}: p<0.05 in {c['frac_p_below_05']:.4f} "
            f"(nominal 0.05), median p {c['median_p']:.3f}")

def readable(name):
    return all(0.02 <= res["calibration"][name][d]["frac_p_below_05"] <= 0.10
               and 0.40 <= res["calibration"][name][d]["median_p"] <= 0.60 for d in DEFS)

res["readable"] = {n: readable(n) for n in res["calibration"]}
for n, ok in res["readable"].items():
    log(f"{n}: calibration {'PASSES' if ok else 'FAILS'} the readability bar")

# The real fits, over the directional set

rows = []
for d in DEFS:
    f = frames[d]
    sn = (f.stage == "Relapse").astype(float).values
    for reg in DIRECTIONAL:
        j = col_ix[d].get(reg)
        if j is None:
            continue
        tab = unit_table(f, centred[d][:, j], sn)
        a, b = decompose(tab), heterogeneity(tab)
        if a is None:
            continue
        r = dict(defs=d, regulon=reg, TF=tf_of(reg), **a)
        if b:
            r.update(b)
        rows.append(r)
    log(f"def{d}: {sum(x['defs'] == d for x in rows)} directional regulons decomposed")

dec = pd.DataFrame(rows)
for d in DEFS:
    m = dec.defs == d
    for src, dst in (("p_rewiring", "q_rewiring"), ("p_composition", "q_composition"),
                     ("p_heterogeneity", "q_heterogeneity")):
        if src in dec.columns:
            dec.loc[m, dst] = bh(dec.loc[m, src].fillna(1.0).values)
dec.to_csv(OUT / "decomposition_by_regulon.csv", index=False)

# Reconcile across vocabularies

w = dec.pivot(index="regulon", columns="defs",
              values=["naive_mean", "rewiring_mean", "composition_mean", "composition_share",
                      "q_rewiring", "q_composition", "q_heterogeneity",
                      "n_celltypes_same_sign", "n_celltypes"])
w.columns = [f"{a}_{b}" for a, b in w.columns]
w = w.dropna(subset=["rewiring_mean_Diagnosis", "rewiring_mean_Relapse"])
w["TF"] = [tf_of(r) for r in w.index]
w["rewiring_sign_agree"] = np.sign(w.rewiring_mean_Diagnosis) == np.sign(w.rewiring_mean_Relapse)
w["composition_sign_agree"] = (np.sign(w.composition_mean_Diagnosis) ==
                               np.sign(w.composition_mean_Relapse))
w["heterogeneous_both"] = (w.q_heterogeneity_Diagnosis < Q_CUT) & \
                          (w.q_heterogeneity_Relapse < Q_CUT)
w["uniform_both"] = (w.q_heterogeneity_Diagnosis >= Q_CUT) & (w.q_heterogeneity_Relapse >= Q_CUT)
w["verdict"] = np.where(w.heterogeneous_both, "CELL-TYPE-RESTRICTED",
                 np.where(w.uniform_both, "UNIFORM ACROSS HIERARCHY", "NOT CONFIRMED"))
w.to_csv(OUT / "decomposition_reconciliation.csv")

res["n_decomposed"] = int(len(w))
res["verdicts"] = w.verdict.value_counts().to_dict()
log("verdicts: " + ", ".join(f"{k} {v}" for k, v in res["verdicts"].items()))

# headline decomposition numbers
res["decomposition"] = {}
for d in DEFS:
    g = dec[dec.defs == d]
    res["decomposition"][d] = {
        "n": int(len(g)),
        "median_abs_naive": float(g.naive_mean.abs().median()),
        "median_abs_rewiring": float(g.rewiring_mean.abs().median()),
        "median_abs_composition": float(g.composition_mean.abs().median()),
        "median_composition_share": float(g.composition_share.median()),
        "frac_composition_share_above_50pct": float((g.composition_share > 0.5).mean()),
        "frac_naive_and_rewiring_same_sign": float(g.same_sign_naive_vs_rewiring.mean()),
        "q05_rewiring": int((g.q_rewiring < Q_CUT).sum()),
        "q05_composition": int((g.q_composition < Q_CUT).sum()),
        "q05_heterogeneity": int((g.q_heterogeneity < Q_CUT).sum())
        if "q_heterogeneity" in g else None,
        "mean_n_celltypes_same_sign": float(g.n_celltypes_same_sign.mean())
        if "n_celltypes_same_sign" in g else None,
    }
    x = res["decomposition"][d]
    log(f"def{d}: median |naive| {x['median_abs_naive']:.5f} = |rewiring| "
        f"{x['median_abs_rewiring']:.5f} + |composition| {x['median_abs_composition']:.5f}; "
        f"median composition share {x['median_composition_share']:.3f}; "
        f"composition dominates in {x['frac_composition_share_above_50pct']:.3f} of regulons")
    log(f"      naive and rewiring agree in sign in "
        f"{x['frac_naive_and_rewiring_same_sign']:.3f}; cell types moving with the overall "
        f"direction: {x['mean_n_celltypes_same_sign']:.2f} of 6")

# which cell types carry the within-cell-type activity change
ctcols = [c for c in dec.columns if c.startswith("beta_")]
res["per_celltype_mean_rewiring"] = {
    d: {c.replace("beta_", ""): float(dec.loc[dec.defs == d, c].mean()) for c in ctcols}
    for d in DEFS}
for d in DEFS:
    s = sorted(res["per_celltype_mean_rewiring"][d].items(), key=lambda kv: -abs(kv[1]))
    log(f"def{d} per-cell-type mean delta: " + ", ".join(f"{k} {v:+.5f}" for k, v in s))

# cross-reference with the Stage-3B selectivity filter
sl = REPO / "longitudinal/healthy_published" / "out" / "s3b_shortlist.csv"
if sl.exists():
    keep = set(pd.read_csv(sl, index_col=0).index)
    sub = w.loc[[i for i in w.index if i in keep]]
    res["selectivity_gated_subset"] = {
        "n": int(len(sub)), "verdicts": sub.verdict.value_counts().to_dict()}
    log(f"Stage-3B gated subset ({len(sub)} regulons): " +
        ", ".join(f"{k} {v}" for k, v in sub.verdict.value_counts().items()))

with open(OUT / "s3_decomposition.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote s3_decomposition.json")
