'Decompose module activity change into within-type and composition components.\n\nThe calculation is exploratory and restricted to the six declared cell types.\nThe output field named within-cell-type activity change denotes within-cell-type activity change; it\ndoes not estimate changes in regulatory edges.'
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
S3 = REPO / "longitudinal/corrected" / "out"
OUT = HERE.parent / "out"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 555                                                                    
MIN_CELLS_PER_UNIT = 20
ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
N_PERM_MODULE = 2000          
rng = np.random.default_rng(SEED)

tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[f6-decomp55] {m}", flush=True)

def bh(p):                                                    
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def load(defs):                                               
    t = pd.read_csv(S3 / f"pseudobulk_def{defs}.csv")
    pu = set(map(tuple, pd.read_csv(S3 / "paired_units.csv").values))
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in pu for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(ADMITTED).values
    return t[m].copy().reset_index(drop=True)

def unit_table(f, y, stage_num):                              
    return pd.DataFrame({"patient": f.patient.values, "celltype": f.celltype.values,
                         "stage_num": stage_num, "n_cells": f.n_cells.values, "y": y})

def decompose(tab, keep_patient=False):                       
    '(A) per-patient naive / within-cell-type activity change / composition, then Wilcoxon over patients.'
    out, pids = [], []
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
        pids.append((p, len(cts)))
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
    if keep_patient:
        res["_per_patient"] = pd.DataFrame(
            {"patient": [x[0] for x in pids], "n_celltypes": [x[1] for x in pids],
             "naive": naive, "rewiring": rewire, "composition": comp})
    return res

def heterogeneity(tab):                                       
    """(B/C) are the per-(patient x cell-type) deltas different across cell types?"""
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

mod = json.loads((REPO / "modules" / "out" / "corrected_module.json").read_text())
MODULE = list(mod["regulons"])
MODULE_TFS = sorted(set(mod["TFs"]))
assert len(MODULE) == 55 == mod["n_regulons"], f"module is not 55: {len(MODULE)}"
assert len(MODULE_TFS) == 31 == mod["n_TFs"], f"module is not 31 TFs: {len(MODULE_TFS)}"
assert len(set(MODULE)) == 55, "duplicate regulons in module definition"
log(f"module of record: {len(MODULE)} regulons / {len(MODULE_TFS)} TFs, "
    f"rule = ax1_directional & ax2_selective & ax3_celltype_robust, frozen {mod['frozen'][:10]}")

committed = pd.read_csv(S3 / "decomposition_by_regulon.csv")
DIR160 = sorted(set(committed.regulon))
missing = sorted(set(MODULE) - set(DIR160))
assert not missing, f"module NOT contained in the committed 160: {missing}"
log(f"containment asserted: all 55 module regulons are a strict subset of the committed "
    f"{len(DIR160)} directional regulons (0 missing)")

ORIENT = mod["orientation"]
assert all(all(r in ORIENT[d] for r in MODULE) for d in DEFS), "orientation missing a module regulon"

frames, M = {}, {}
substrate = {}
for d in DEFS:
    f = load(d)
    regs = [c for c in f.columns if c not in META]
    miss = [r for r in MODULE if r not in regs]
    assert not miss, f"module regulons missing from pseudobulk def{d}: {miss}"
    A = f[regs].values
    C = A - A.mean(axis=1, keepdims=True)                     # offset-centred, full-panel offset
    ix = [regs.index(r) for r in MODULE]
    sgn = np.array([ORIENT[d][r] for r in MODULE], float)
    frames[d] = f
    M[d] = {"offsetcentred": C[:, ix].mean(axis=1),           # sensitivity (unoriented)
            "raw": A[:, ix].mean(axis=1),                     # sensitivity
            "oriented_offsetcentred": (C[:, ix] * sgn).mean(axis=1)}  # primary (orientation-signed, offset-centred)
    substrate[d] = {"file": f"longitudinal/corrected/out/pseudobulk_def{d}.csv",
                    "n_rows_admitted": int(len(f)), "n_patients": int(f.patient.nunique()),
                    "n_units": int(f.groupby(["patient", "celltype"]).ngroups),
                    "n_regulon_columns_in_panel": int(len(regs)),
                    "n_module_regulons_present": 55,
                    "n_orientation_positive": int((sgn > 0).sum())}
    log(f"def{d}: {len(f)} admitted rows = {f.patient.nunique()} patients x "
        f"{f.groupby(['patient','celltype']).ngroups} paired units x 2 stages; "
        f"panel {len(regs)} regulons; module 55/55 present")

def build_fast(d):
    f = frames[d]
    key = (f.patient + "|" + f.celltype).values
    order = np.lexsort((f.stage.values, key))                 # within unit: Diagnosis then Relapse
    k, pat, ct = key[order], f.patient.values[order], f.celltype.values[order]
    n = f.n_cells.values[order].astype(float)
    assert (k[0::2] == k[1::2]).all(), "units are not exactly paired"
    assert (f.stage.values[order][0::2] == "Diagnosis").all()
    return dict(order=order, upat=pat[0::2], uct=ct[0::2],
                n0=n[0::2], n1=n[1::2])

def fast_stats(fs, y, flip):
    """Exactly decompose()+heterogeneity() for a permuted stage assignment, vectorised."""
    y0, y1 = y[fs["order"]][0::2], y[fs["order"]][1::2]
    mdx = np.where(flip, y1, y0); mrl = np.where(flip, y0, y1)
    wdx = np.where(flip, fs["n1"], fs["n0"]); wrl = np.where(flip, fs["n0"], fs["n1"])
    naive, rewire = [], []
    for p in np.unique(fs["upat"]):
        m = fs["upat"] == p
        if m.sum() < 2:
            continue
        a, b = wdx[m] / wdx[m].sum(), wrl[m] / wrl[m].sum()
        naive.append((b * mrl[m]).sum() - (a * mdx[m]).sum())
        rewire.append((a * (mrl[m] - mdx[m])).sum())
    if len(naive) < 5:
        return None
    naive = np.asarray(naive); rewire = np.asarray(rewire); comp = naive - rewire
    pr = float(stats.wilcoxon(rewire, zero_method="wilcox").pvalue) if not np.allclose(rewire, 0) else 1.0
    pc = float(stats.wilcoxon(comp, zero_method="wilcox").pvalue) if not np.allclose(comp, 0) else 1.0
    delta = mrl - mdx
    groups = [delta[fs["uct"] == c] for c in np.unique(fs["uct"])]
    groups = [g for g in groups if len(g) >= 3]
    ph = float(stats.kruskal(*groups).pvalue) if len(groups) >= 3 else np.nan
    return pr, pc, ph

FAST = {d: build_fast(d) for d in DEFS}

for d in DEFS:                                                # equivalence assertion, no permutation
    y = M[d]["offsetcentred"]
    sn = (frames[d].stage == "Relapse").astype(float).values
    ref_a = decompose(unit_table(frames[d], y, sn))
    ref_b = heterogeneity(unit_table(frames[d], y, sn))
    got = fast_stats(FAST[d], y, np.zeros(len(FAST[d]["upat"]), bool))
    assert np.isclose(ref_a["p_rewiring"], got[0]) and np.isclose(ref_a["p_composition"], got[1]) \
        and np.isclose(ref_b["p_heterogeneity"], got[2]), f"fast path diverges on def{d}"
log("fast permutation path ASSERTED identical to the verbatim reference estimators")

res = {
    "analysis": "F6 decomposition on the module of record",
    "status": "POST-HOC AND DESCRIPTIVE",
    "status_note": "The 55-regulon module was itself defined post-hoc (2026-08-30); "
                   "this analysis is conditioned on it and is not confirmatory. The "
                   "pre-declared Stage-3 primary (0 regulons CONFIRMED, s3_results.json) is "
                   "unchanged by anything here.",
    "prespec": "analysis/F6_DECOMPOSITION_PRESPEC.md (written before any number below was computed)",
    "reference_implementation": "longitudinal/corrected/scripts/decomposition.py",
    "seed": SEED,
    "module": {"file": "modules/out/corrected_module.json",
               "rule": "ax1_directional & ax2_selective & ax3_celltype_robust",
               "arm": "corrected (repair)", "n_regulons": 55, "n_TFs": 31,
               "frozen": mod["frozen"][:10],
               "containment_in_committed_160_asserted": True, "n_missing_from_160": 0},
    "scale": "offset-centred, orientation-signed (unoriented offset-centred and raw run as sensitivities)",
    "admitted_celltypes": ADMITTED,
    "admitted_celltypes_source": "longitudinal/corrected/out/s3_decomposition.json -> "
                                 "admitted_celltypes",
    "composition_component_is_a_lower_bound": True,
    "composition_lower_bound_note": "measured over 6 admitted malignant cell types only; the "
                                    "label-level shift spans all 20 types, so the "
                                    "composition component and its share are LOWER BOUNDS",
    "substrate": substrate,
    "n_permutations_calibration": N_PERM_MODULE,
}

res["calibration"] = {}
for scale in ("offsetcentred", "raw", "oriented_offsetcentred"):
    res["calibration"][scale] = {n: {} for n in ("rewiring", "composition", "heterogeneity")}
    for d in DEFS:
        fs, y = FAST[d], M[d][scale]
        nu = len(fs["upat"])
        ps = {"rewiring": [], "composition": [], "heterogeneity": []}
        for _ in range(N_PERM_MODULE):
            got = fast_stats(fs, y, rng.integers(0, 2, nu).astype(bool))
            if got is None:
                continue
            ps["rewiring"].append(got[0]); ps["composition"].append(got[1])
            if np.isfinite(got[2]):
                ps["heterogeneity"].append(got[2])
        for name, v in ps.items():
            v = np.asarray(v, float)
            res["calibration"][scale][name][d] = {
                "n_fits": int(len(v)), "frac_p_below_05": float((v < 0.05).mean()),
                "median_p": float(np.median(v)),
                "ks_uniform_p": float(stats.kstest(v, "uniform").pvalue)}
            c = res["calibration"][scale][name][d]
            log(f"calibration [{scale}] {name} def{d}: p<0.05 in {c['frac_p_below_05']:.4f} "
                f"(nominal 0.05), median p {c['median_p']:.3f}")

def readable(scale, name):                                                            
    return all(0.02 <= res["calibration"][scale][name][d]["frac_p_below_05"] <= 0.10
               and 0.40 <= res["calibration"][scale][name][d]["median_p"] <= 0.60 for d in DEFS)

res["readable"] = {s: {n: readable(s, n) for n in res["calibration"][s]}
                   for s in res["calibration"]}
for s, dd in res["readable"].items():
    for n, ok in dd.items():
        log(f"[{s}] {n}: calibration {'PASSES' if ok else 'FAILS'} the readability bar")

# (A) + (C) — the real module-level fits

res["module_decomposition"] = {}
pp_rows = []
for scale in ("offsetcentred", "raw", "oriented_offsetcentred"):
    res["module_decomposition"][scale] = {}
    for d in DEFS:
        f = frames[d]
        sn = (f.stage == "Relapse").astype(float).values
        tab = unit_table(f, M[d][scale], sn)
        a = decompose(tab, keep_patient=True)
        b = heterogeneity(tab)
        pp = a.pop("_per_patient")
        pp.insert(0, "scale", scale); pp.insert(0, "vocabulary", d)
        pp_rows.append(pp)
        v = pp[["naive", "rewiring", "composition"]]
        summ = dict(a)
        for col in ("naive", "rewiring", "composition"):
            x = pp[col].values
            summ[f"{col}_median"] = float(np.median(x))
            summ[f"{col}_sd"] = float(x.std(ddof=1))
            summ[f"{col}_iqr_lo"] = float(np.percentile(x, 25))
            summ[f"{col}_iqr_hi"] = float(np.percentile(x, 75))
            summ[f"frac_patients_{col}_positive"] = float((x > 0).mean())
        summ["frac_patients_abs_composition_gt_abs_rewiring"] = \
            float((pp.composition.abs() > pp.rewiring.abs()).mean())
        summ["frac_patients_naive_rewiring_same_sign"] = \
            float((np.sign(pp.naive) == np.sign(pp.rewiring)).mean())
        summ["pearson_naive_rewiring"] = float(np.corrcoef(pp.naive, pp.rewiring)[0, 1])
        summ["verdict_primary_cut_0.50"] = ("MOSTLY COMPOSITION" if summ["composition_share"] > 0.5
                                            else "MOSTLY REWIRING" if summ["composition_share"] < 0.5
                                            else "UNDETERMINED")
        summ["verdict_decisive_band"] = ("DECISIVELY REWIRING" if summ["composition_share"] <= 1/3
                                         else "DECISIVELY COMPOSITION" if summ["composition_share"] >= 2/3
                                         else "MIXED")
        if b:
            summ.update(b)
        res["module_decomposition"][scale][d] = summ
        log(f"[{scale}] def{d}: naive {summ['naive_mean']:+.6f} = rewiring "
            f"{summ['rewiring_mean']:+.6f} + composition {summ['composition_mean']:+.6f}; "
            f"share {summ['composition_share']:.4f} -> {summ['verdict_primary_cut_0.50']}")

pp_all = pd.concat(pp_rows, ignore_index=True)
pp_all["abs_composition_gt_abs_rewiring"] = pp_all.composition.abs() > pp_all.rewiring.abs()
pp_all["naive_rewiring_same_sign"] = np.sign(pp_all.naive) == np.sign(pp_all.rewiring)
pp_all.to_csv(OUT / "f6_decomposition_module_by_patient.csv", index=False)
log(f"wrote f6_decomposition_module_by_patient.csv ({len(pp_all)} rows)")

res["module_verdict"] = {}
for scale in res["module_decomposition"]:
    v = [res["module_decomposition"][scale][d]["verdict_primary_cut_0.50"] for d in DEFS]
    b = [res["module_decomposition"][scale][d]["verdict_decisive_band"] for d in DEFS]
    res["module_verdict"][scale] = {
        "per_vocabulary": dict(zip(DEFS, v)),
        "both_vocabularies_agree": bool(v[0] == v[1]),
        "verdict": v[0] if v[0] == v[1] else "NOT CONFIRMED (vocabularies disagree)",
        "decisive_band_per_vocabulary": dict(zip(DEFS, b)),
        "decisive_band_verdict": b[0] if b[0] == b[1] else "NOT CONFIRMED (vocabularies disagree)"}

# (B) — per-regulon on the 55, BH recomputed within the family

sub = committed[committed.regulon.isin(MODULE)].copy()
assert sub.regulon.nunique() == 55 and len(sub) == 110, f"expected 110 rows, got {len(sub)}"
sub = sub.rename(columns={"q_rewiring": "q_rewiring_committed_160",
                          "q_composition": "q_composition_committed_160",
                          "q_heterogeneity": "q_heterogeneity_committed_160"})
status_changes = {}
for d in DEFS:
    m = (sub.defs == d).values
    for test in ("rewiring", "composition", "heterogeneity"):
        p = sub.loc[m, f"p_{test}"].fillna(1.0).values
        q55 = bh(p)
        sub.loc[m, f"q_{test}_set55"] = q55
        q160 = sub.loc[m, f"q_{test}_committed_160"].fillna(1.0).values
        gain = int(((q55 < Q_CUT) & (q160 >= Q_CUT)).sum())
        lose = int(((q55 >= Q_CUT) & (q160 < Q_CUT)).sum())
        status_changes[f"{test}|{d}"] = {
            "n_family": int(m.sum()),
            "n_sig_committed_160_wide_q": int((q160 < Q_CUT).sum()),
            "n_sig_recomputed_set55_q": int((q55 < Q_CUT).sum()),
            "n_became_significant": gain, "n_lost_significance": lose,
            "n_status_changes": gain + lose}
        log(f"BH {test} def{d}: sig {int((q160<Q_CUT).sum())} (160-wide q) -> "
            f"{int((q55<Q_CUT).sum())} (55-wide q); {gain+lose} status change(s)")
res["per_regulon_bh_status_changes"] = status_changes
res["per_regulon_bh_total_status_changes"] = int(sum(v["n_status_changes"]
                                                     for v in status_changes.values()))
sub["in_module_55"] = True
sub.to_csv(OUT / "f6_decomposition_set55_per_regulon.csv", index=False)
log(f"wrote f6_decomposition_set55_per_regulon.csv ({len(sub)} rows)")

res["per_regulon_summary_set55"] = {}
ctcols = [c for c in sub.columns if c.startswith("beta_")]
for d in DEFS:
    g = sub[sub.defs == d]
    res["per_regulon_summary_set55"][d] = {
        "n": int(len(g)),
        "median_abs_naive": float(g.naive_mean.abs().median()),
        "median_abs_rewiring": float(g.rewiring_mean.abs().median()),
        "median_abs_composition": float(g.composition_mean.abs().median()),
        "median_composition_share": float(g.composition_share.median()),
        "mean_composition_share": float(g.composition_share.mean()),
        "iqr_composition_share": [float(np.percentile(g.composition_share, 25)),
                                  float(np.percentile(g.composition_share, 75))],
        "frac_composition_share_above_50pct": float((g.composition_share > 0.5).mean()),
        "n_composition_share_above_50pct": int((g.composition_share > 0.5).sum()),
        "frac_naive_and_rewiring_same_sign": float(g.same_sign_naive_vs_rewiring.mean()),
        "n_naive_and_rewiring_same_sign": int(g.same_sign_naive_vs_rewiring.sum()),
        "q05_rewiring_set55": int((g.q_rewiring_set55 < Q_CUT).sum()),
        "q05_composition_set55": int((g.q_composition_set55 < Q_CUT).sum()),
        "q05_heterogeneity_set55": int((g.q_heterogeneity_set55 < Q_CUT).sum()),
        "q05_rewiring_committed_160": int((g.q_rewiring_committed_160 < Q_CUT).sum()),
        "q05_composition_committed_160": int((g.q_composition_committed_160 < Q_CUT).sum()),
        "q05_heterogeneity_committed_160": int((g.q_heterogeneity_committed_160 < Q_CUT).sum()),
        "mean_n_celltypes_same_sign": float(g.n_celltypes_same_sign.mean()),
        "n_patients_range": [int(g.n_patients.min()), int(g.n_patients.max())],
        "per_celltype_mean_delta": {c.replace("beta_", ""): float(g[c].mean()) for c in ctcols},
        "regulons_q05_rewiring_set55": sorted(g.loc[g.q_rewiring_set55 < Q_CUT, "regulon"]),
    }

# cross-vocabulary reconciliation on the 55
w = sub.pivot(index="regulon", columns="defs",
              values=["naive_mean", "rewiring_mean", "composition_mean", "composition_share",
                      "q_rewiring_set55", "q_composition_set55", "q_heterogeneity_set55"])
w.columns = [f"{a}_{b}" for a, b in w.columns]
res["per_regulon_cross_vocabulary_set55"] = {
    "n": int(len(w)),
    "rewiring_sign_agree": int((np.sign(w.rewiring_mean_Diagnosis) ==
                                np.sign(w.rewiring_mean_Relapse)).sum()),
    "composition_sign_agree": int((np.sign(w.composition_mean_Diagnosis) ==
                                   np.sign(w.composition_mean_Relapse)).sum()),
    "heterogeneous_both_vocabularies": int(((w.q_heterogeneity_set55_Diagnosis < Q_CUT) &
                                            (w.q_heterogeneity_set55_Relapse < Q_CUT)).sum()),
    "uniform_both_vocabularies": int(((w.q_heterogeneity_set55_Diagnosis >= Q_CUT) &
                                      (w.q_heterogeneity_set55_Relapse >= Q_CUT)).sum()),
    "rewiring_q05_both_vocabularies": int(((w.q_rewiring_set55_Diagnosis < Q_CUT) &
                                           (w.q_rewiring_set55_Relapse < Q_CUT)).sum()),
    "composition_share_above_50pct_both": int(((w.composition_share_Diagnosis > 0.5) &
                                               (w.composition_share_Relapse > 0.5)).sum()),
}

comm = json.loads((S3 / "s3_decomposition.json").read_text())
res["per_regulon_calibration_committed"] = {
    "source": "longitudinal/corrected/out/s3_decomposition.json -> calibration / readable",
    "note": "governs the per-regulon p-values, which are re-read not re-inferred; NOT re-run here",
    "calibration": comm["calibration"], "readable": comm["readable"],
    "directional_set_size": comm["directional_set_size"]}

with open(OUT / "f6_decomposition.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote f6_decomposition.json")
