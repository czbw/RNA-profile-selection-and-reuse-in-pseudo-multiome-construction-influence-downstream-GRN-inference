'SENSITIVITY (Supplementary Data 6): the 55-regulon nomination and the pediatric decomposition re-run with seven cell types (CLP added) and a cell-type bar of 5 of 7.'
import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True, help="working repository root")
ap.add_argument("--out", required=True, help="output directory (created)")
ap.add_argument("--bar", type=int, default=5, help="required cell types out of the admitted set")
args = ap.parse_args()
REPO = Path(args.repo)
OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)
S3 = REPO / "longitudinal/corrected" / "out"

MIN_CELLS_PER_UNIT = 20
MIN_PATIENTS_PER_CELLTYPE = 5
SIX = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
META = ["stage", "patient", "celltype", "driver", "n_cells"]
DEFS = ("Diagnosis", "Relapse")
Q_CUT = 0.05
tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[sens-clp7] {m}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

paired = pd.read_csv(S3 / "paired_units.csv")
PU = set(map(tuple, paired.values))
ct_n = paired.celltype.value_counts()
SEVEN = sorted(ct_n[ct_n >= MIN_PATIENTS_PER_CELLTYPE].index)
assert set(SIX) < set(SEVEN) and set(SEVEN) - set(SIX) == {"CLP"}, \
    f"admission rule does not give six + CLP: {SEVEN}"
log(f"admitted by rule (>= {MIN_PATIENTS_PER_CELLTYPE} paired patients): {SEVEN}")

# --------------------------------------------------------------------------- #
# substrate
# --------------------------------------------------------------------------- #
def load(d, adm):
    t = pd.read_csv(S3 / f"pseudobulk_def{d}.csv")
    m = (t.n_cells >= MIN_CELLS_PER_UNIT).values & \
        np.array([tuple(x) in PU for x in t[["patient", "celltype"]].values]) & \
        t.celltype.isin(adm).values
    f = t[m].copy().reset_index(drop=True)
    f["stage_num"] = (f.stage == "Relapse").astype(float)
    regs = [c for c in f.columns if c not in META + ["stage_num"]]
    A = f[regs].values
    return f, regs, A, A - A.mean(axis=1, keepdims=True)

def unit_deltas(f, C):
    """paired (patient, cell type) units -> unit deltas (relapse minus diagnosis), per regulon."""
    idx = pd.DataFrame({"patient": f.patient, "celltype": f.celltype,
                        "st": f.stage_num == 1.0, "i": np.arange(len(f))})
    p0 = idx[~idx.st].set_index(["patient", "celltype"]).i
    p1 = idx[idx.st].set_index(["patient", "celltype"]).i
    both = p0.index.intersection(p1.index)
    D = C[p1.loc[both].values] - C[p0.loc[both].values]
    return D, np.array([b[0] for b in both]), np.array([b[1] for b in both])

def e2_and_signs(d, adm):
    """per regulon: E2 mean change, and number of cell types with the overall sign."""
    f, regs, A, C = load(d, adm)
    D, pat, ct = unit_deltas(f, C)
    rows = []
    for j, r in enumerate(regs):
        df = pd.DataFrame({"p": pat, "c": ct, "d": D[:, j]})
        beta = df.groupby("p").d.mean().mean()
        perct = df.groupby("c").d.mean()
        same = int((np.sign(perct.values) == np.sign(beta)).sum()) if beta != 0 else 0
        rows.append((r, beta, same, len(perct)))
    return pd.DataFrame(rows, columns=["regulon", "beta", "same", "nct"]).set_index("regulon")

EST = {(tag, d): e2_and_signs(d, adm) for tag, adm in (("six", SIX), ("seven", SEVEN)) for d in DEFS}

# --------------------------------------------------------------------------- #
# controls: the six-type re-implementation reproduces the stored values
# --------------------------------------------------------------------------- #
e2s = pd.read_csv(S3 / "e2_paired_wilcoxon_reconciliation.csv", index_col=0)
decs = pd.read_csv(S3 / "decomposition_reconciliation.csv", index_col=0)
fun = pd.read_csv(S3 / "nomination_funnel.csv").set_index("regulon")
mod = json.loads((REPO / "modules" / "out" / "corrected_module.json").read_text())
MOD55 = list(mod["regulons"])
assert len(fun) == 201
for d in DEFS:
    r = EST[("six", d)]
    diff = float((r.loc[e2s.index, "beta"] - e2s[f"beta_{d}"]).abs().max())
    assert diff < 1e-12 and len(e2s.index.intersection(r.index)) == 201, f"C1 e2 beta {d}: {diff}"
    ci = decs.index
    mis = int((r.loc[ci, "same"] != decs[f"n_celltypes_same_sign_{d}"]).sum())
    assert mis == 0 and len(ci) == 160, f"C2 same-sign counts {d}: {mis} mismatches"
log("C1/C2 PASS: six-type re-implementation reproduces the stored e2 betas (201) and same-sign counts (160)")

ax2 = ~fun.compartments_agree

def funnel(tag, bar):
    R = lambda d, c: EST[(tag, d)][c].reindex(fun.index)
    ax1 = pd.Series(np.sign(R("Diagnosis", "beta")).values == np.sign(R("Relapse", "beta")).values,
                    index=fun.index)
    ax3 = (R("Diagnosis", "same") >= bar) & (R("Relapse", "same") >= bar)
    return ax1, ax3

a1, a3 = funnel("six", 4)
assert bool((a1 == fun.ax1_directional).all()), "C3 ax1 differs from stored"
f6 = a1 & ax2 & a3
assert int(a1.sum()) == 160 and int((a1 & ax2).sum()) == 73 and int(f6.sum()) == 55, "C3 counts"
assert set(f6[f6].index) == set(MOD55), "C3 membership differs from module_definition_T3_repair"
log("C3 PASS: six-type funnel 201 -> 160 -> 73 -> 55 with the exact module of record")

# --------------------------------------------------------------------------- #
# the seven-type funnel
# --------------------------------------------------------------------------- #
summary = {"admitted_celltypes_seven": SEVEN, "paired_patients_per_celltype": ct_n.to_dict(),
           "bar_sensitivity": args.bar, "primary_bar": "4 of 6",
           "controls_passed": ["C1 e2 betas", "C2 same-sign counts", "C3 funnel and membership"],
           "ax2_note": "stored; already computed on seven malignant types (s3br_13)",
           "funnels": {}}
table = pd.DataFrame(index=fun.index)
table["TF"] = [tf_of(i) for i in fun.index]
for tag in ("six", "seven"):
    for d in DEFS:
        table[f"beta_{tag}_{d}"] = EST[(tag, d)].beta.reindex(fun.index)
        table[f"same_{tag}_{d}"] = EST[(tag, d)].same.reindex(fun.index)
table["ax1_six"], _ = funnel("six", 4)
table["ax1_seven"], _ = funnel("seven", 4)
table["ax2_selective"] = ax2
table["in_primary_55"] = table.index.isin(MOD55)
for bar in (4, 5, 6):
    ax1, ax3 = funnel("seven", bar)
    sel = ax1 & ax2 & ax3
    table[f"in_seven_bar{bar}"] = sel
    s, base = set(sel[sel].index), set(MOD55)
    summary["funnels"][f"{bar}_of_7"] = {
        "pool": 201, "ax1": int(ax1.sum()), "ax1_and_ax2": int((ax1 & ax2).sum()),
        "final": len(s), "TFs": len({tf_of(r) for r in s}),
        "kept_from_primary_55": len(s & base), "new": sorted(s - base), "lost": sorted(base - s)}
    log(f"{bar} of 7: 201 -> {int(ax1.sum())} -> {int((ax1 & ax2).sum())} -> {len(s)} "
        f"(kept {len(s & base)}, new {len(s - base)}, lost {len(base - s)})")
summary["funnels"]["primary_4_of_6"] = {"pool": 201, "ax1": 160, "ax1_and_ax2": 73, "final": 55,
                                         "TFs": len({tf_of(r) for r in MOD55})}
summary["ax1_flips_six_to_seven"] = sorted(table.index[table.ax1_six != table.ax1_seven])
table.to_csv(OUT / "sens_clp7_funnel_by_regulon.csv")

SENS = list(table.index[table[f"in_seven_bar{args.bar}"]])
sens_tfs = sorted({tf_of(r) for r in SENS})
summary["sensitivity_module"] = {"bar": f"{args.bar} of 7", "n_regulons": len(SENS),
                                 "n_TFs": len(sens_tfs), "regulons": SENS, "TFs": sens_tfs}
summary["TF_level_vs_primary"] = {
    "primary_TFs": len({tf_of(r) for r in MOD55}),
    "shared_TFs": len(set(sens_tfs) & {tf_of(r) for r in MOD55}),
    "sensitivity_only_TFs": sorted(set(sens_tfs) - {tf_of(r) for r in MOD55}),
    "primary_only_TFs": sorted({tf_of(r) for r in MOD55} - set(sens_tfs))}

# --------------------------------------------------------------------------- #
# Figure 4d analogue: overlap with the mixed-model regulons and the AP-1/IEG list
# --------------------------------------------------------------------------- #
rec = json.load(open(S3 / "s3_recalibrated.json", encoding="utf-8"))
LMM26 = rec["results"]["E1_random_intercept"]["confirmed_regulons"]
man = json.load(open(REPO / "figures" / "data" / "manifest.json", encoding="utf-8"))
ieg = set(man["s3b.ap1.ieg_definition"]["value"])
t26 = {r.split("_")[0] for r in LMM26}
fig4d = {}
for name, S in (("primary_55", set(MOD55)), (f"sensitivity_{args.bar}of7", set(SENS))):
    tS = {tf_of(r) for r in S}
    fig4d[name] = {"n_regulons": len(S), "intersection_with_LMM26_regulons": len(S & set(LMM26)),
                   "jaccard_regulon": len(S & set(LMM26)) / len(S | set(LMM26)),
                   "intersection_TFs": len(tS & t26), "jaccard_TF": len(tS & t26) / len(tS | t26),
                   "n_regulons_whose_TF_in_IEG_list": int(sum(tf_of(r) in ieg for r in S)),
                   "ieg_list_size": len(ieg)}
summary["fig4d_analogue"] = fig4d
assert fig4d["primary_55"]["intersection_with_LMM26_regulons"] == 3, "control: primary overlap is 3"

# --------------------------------------------------------------------------- #
# pediatric decomposition with seven types (estimators verbatim from f6_07 / s3r_21)
# --------------------------------------------------------------------------- #
def decompose(tab, keep_patient=False):
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
    res = dict(n_patients=len(a), naive_mean=float(naive.mean()), rewiring_mean=float(rewire.mean()),
               composition_mean=float(comp.mean()),
               p_rewiring=float(stats.wilcoxon(rewire, zero_method="wilcox").pvalue)
               if not np.allclose(rewire, 0) else 1.0,
               p_composition=float(stats.wilcoxon(comp, zero_method="wilcox").pvalue)
               if not np.allclose(comp, 0) else 1.0)
    denom = abs(res["rewiring_mean"]) + abs(res["composition_mean"])
    res["composition_share"] = float(abs(res["composition_mean"]) / denom) if denom > 0 else np.nan
    res["same_sign_naive_vs_rewiring"] = bool(np.sign(res["naive_mean"]) == np.sign(res["rewiring_mean"]))
    return res

def heterogeneity(tab):
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
    return dict(p_heterogeneity=float(kw.pvalue), n_celltypes=len(per_ct), n_celltypes_same_sign=same)

ORIENT = mod["orientation"]
dec_out = {}
per_reg_rows = []
for d in DEFS:
    f, regs, A, C = load(d, SEVEN)
    ix = {r: regs.index(r) for r in regs}
    for label, members in (("primary_55_with_CLP", MOD55), (f"sensitivity_{args.bar}of7", SENS)):
        jj = [ix[r] for r in members]
        sgn = np.array([ORIENT[d][r] for r in members], float)
        scales = {"offsetcentred": C[:, jj].mean(axis=1), "raw": A[:, jj].mean(axis=1),
                  "oriented_offsetcentred": (C[:, jj] * sgn).mean(axis=1)}
        for sc, y in scales.items():
            tab = pd.DataFrame({"patient": f.patient.values, "celltype": f.celltype.values,
                                "stage_num": f.stage_num.values, "n_cells": f.n_cells.values, "y": y})
            r = decompose(tab)
            r.update(heterogeneity(tab) or {})
            dec_out.setdefault(label, {}).setdefault(sc, {})[d] = r
            log(f"{label} [{sc}] {d}: naive {r['naive_mean']:+.6f} = within {r['rewiring_mean']:+.6f} "
                f"+ composition {r['composition_mean']:+.6f}; share {r['composition_share']:.4f}")
    # per-regulon components for the union of both sets
    for reg in sorted(set(MOD55) | set(SENS)):
        tab = pd.DataFrame({"patient": f.patient.values, "celltype": f.celltype.values,
                            "stage_num": f.stage_num.values, "n_cells": f.n_cells.values,
                            "y": C[:, ix[reg]]})
        r = decompose(tab) or {}
        r.update(heterogeneity(tab) or {})
        per_reg_rows.append({"vocabulary": d, "regulon": reg, "in_primary_55": reg in MOD55,
                             "in_sensitivity": reg in SENS, **r})

# control: the frozen 55 with seven types equals the earlier whole-script sandbox run, if present
chk = REPO / "analysis" / "out" / "clp_rerun_20261021" / "seven_type" / "f6_decomposition.json"
if chk.exists():
    ref = json.load(open(chk))
    for sc in ("offsetcentred", "raw", "oriented_offsetcentred"):
        for d in DEFS:
            for k in ("naive_mean", "rewiring_mean", "composition_mean", "composition_share",
                      "p_rewiring", "p_composition", "p_heterogeneity"):
                a, b = dec_out["primary_55_with_CLP"][sc][d][k], ref["module_decomposition"][sc][d][k]
                assert abs(a - b) < 1e-9, f"C4 {sc} {d} {k}: {a} vs {b}"
    summary["controls_passed"].append("C4 frozen 55 with CLP equals the earlier sandbox run of f6_07")
    log("C4 PASS: frozen 55 with seven types equals the earlier f6_07 sandbox run")

pr = pd.DataFrame(per_reg_rows)
for d in DEFS:
    for flag in ("in_primary_55", "in_sensitivity"):
        m = (pr.vocabulary == d) & pr[flag]
        for test in ("rewiring", "composition", "heterogeneity"):
            pr.loc[m, f"q_{test}_within_{flag}"] = bh(pr.loc[m, f"p_{test}"].fillna(1.0).values)
pr.to_csv(OUT / "sens_clp7_per_regulon_decomposition.csv", index=False)
summary["decomposition"] = dec_out
summary["decomposition_note"] = ("seven admitted cell types; patient-level Wilcoxon p-values; "
                                 "permutation calibration NOT re-run; orientation frozen from the "
                                 "module definition (carried verbatim for all 201 regulons)")
json.dump(summary, open(OUT / "sens_clp7_summary.json", "w"), indent=2, default=str)
log("wrote sens_clp7_summary.json, sens_clp7_funnel_by_regulon.csv, sens_clp7_per_regulon_decomposition.csv")
