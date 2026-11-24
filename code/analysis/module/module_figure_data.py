#!/usr/bin/env python
"""Support tables for the Fig. 6 module panels: mixture/activity, arm enrichment, ribosome-drop survival."""
import itertools
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "analysis", "out")
REF = os.path.join(ROOT, "analysis", "ref")
R3 = os.path.join(ROOT, "longitudinal/corrected", "out")

VOCABS = ("Diagnosis", "Relapse")
SIZE_MIN, SIZE_MAX, Q_BAR = 10, 500, 0.05
LIBS = {"GOBP": "GO_Biological_Process_2025.gmt",
        "REAC": "Reactome_Pathways_2024.gmt",
        "HALL": "MSigDB_Hallmark_2020.gmt"}

R = {"_status": "SUPPORT FILE for figures/main/make_fig05_module.R. POST-HOC and DESCRIPTIVE, "
                "inheriting the post-hoc status of the analyses it reads. No new inference: this re-reads committed files, aggregates, and "
                "runs the same Fisher/BH enrichment the committed runs use.",
     "_arm_definition": "s1 is the sign of the TF-expression-to-target-gene correlation. '+' arms "
                        "are activating, '-' arms repressing. Pooling different TFs' arms of the "
                        "same sign is legitimate; merging one TF's own two arms is not, and is not "
                        "done here (see the unit-of-analysis correction, f6_06)."}

def log(m):
    print(f"[f6sup] {m}", flush=True)

def clean(gs):
    return {str(g).strip() for g in gs if str(g).strip()}

def bh(p):
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def read_gmt(p):
    s = {}
    for line in open(p, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            g = clean(f[2:])
            if g:
                s[f[0]] = g
    return s

MOD = json.load(open(os.path.join(ROOT, "modules", "out",
                                  "corrected_module.json"), encoding="utf-8"))
M55, POOL = list(MOD["regulons"]), list(MOD["background_regulons"])
assert len(M55) == 55 and len(POOL) == 201
TGT = {v: {r: clean(MOD["targets"][v][r]) for r in POOL} for v in VOCABS}
UNIV = {}
for v in VOCABS:
    u = pd.read_csv(os.path.join(OUT, f"f1_universe_{v}.csv"))
    UNIV[v] = clean(u.symbol[u.symbol.notna()])
GMT = {k: read_gmt(os.path.join(REF, fn)) for k, fn in LIBS.items()}
DEC = json.load(open(os.path.join(OUT, "f6_decomposition.json"), encoding="utf-8"))
ADM = list(DEC["admitted_celltypes"])
RIB = json.load(open(os.path.join(OUT, "f6_ribosome_check.json"), encoding="utf-8"))
RIBO = set(RIB["gene_sets"]["DROP_A"]["genes"])
log(f"module {len(M55)} regulons; pool {len(POOL)}; admitted cell types {ADM}")
log(f"ribosomal gene set (DROP_A, committed): {len(RIBO)} genes")

# ======================================================================== (1) mixture + activity
S55J = json.load(open(os.path.join(OUT, "f6_set55.json"), encoding="utf-8"))
CT = pd.DataFrame(S55J["celltype_summary"])          
rows = []
for v in VOCABS:
    pb = pd.read_csv(os.path.join(R3, f"pseudobulk_def{v}.csv"),
                     usecols=["patient", "celltype", "stage", "n_cells"])
    pb = pb[pb.celltype.isin(ADM)]
    tot = pb.groupby("stage").n_cells.sum()
    w = pb.groupby(["stage", "celltype"]).n_cells.sum() / tot
    for st in ("Diagnosis", "Relapse"):
        for c in ADM:
            act = CT[(CT.stage == st) & (CT.celltype == c)].mean_AUC
            rows.append({"vocabulary": v, "stage": st, "celltype": c,
                         "cell_share": float(w.loc[(st, c)]),
                         "n_cells": int(pb.n_cells[(pb.stage == st) & (pb.celltype == c)].sum()),
                         "module_mean_AUC": float(act.iloc[0]) if len(act) else np.nan})
MIX = pd.DataFrame(rows)
MIX.to_csv(os.path.join(OUT, "f6_fig06_mixture.csv"), index=False)

# the flatness statistic -- the panel's real content
d0 = MIX[(MIX.vocabulary == "Diagnosis") & (MIX.stage == "Diagnosis")]
d1 = MIX[(MIX.vocabulary == "Diagnosis") & (MIX.stage == "Relapse")]
w0 = d0.set_index("celltype").cell_share
w1 = d1.set_index("celltype").cell_share
m0 = d0.set_index("celltype").module_mean_AUC
m1 = d1.set_index("celltype").module_mean_AUC
R["mixture_and_activity"] = {
    "admitted_celltypes": ADM,
    "n_admitted": len(ADM),
    "n_celltypes_in_the_full_profile": int(CT.celltype.nunique()),
    "n_celltypes_in_the_substrate": int(pd.read_csv(
        os.path.join(R3, "pseudobulk_defDiagnosis.csv"), usecols=["celltype"]).celltype.nunique()),
    "_denominator_note": (
        "Three different counts are easy to confuse. The Stage-3 pseudobulk the decomposition "
        "reads carries 19 cell types. The module's own activity profile lists 16. Six are "
        "admitted. The honest denominator for 'how much of the hierarchy is the composition term "
        "blind to' is 19, the set the decomposition could in principle have used. A figure of 20 "
        "has circulated and matches no committed file."),
    "activity_spread_Diagnosis": {
        "min": float(m0.min()), "max": float(m0.max()),
        "range": float(m0.max() - m0.min()), "mean": float(m0.mean()),
        "range_as_pct_of_mean": float(100 * (m0.max() - m0.min()) / m0.mean())},
    "largest_share_change": {
        "celltype": str((w1 - w0).abs().idxmax()),
        "from": float(w0[(w1 - w0).abs().idxmax()]), "to": float(w1[(w1 - w0).abs().idxmax()]),
        "abs_change": float((w1 - w0).abs().max())},
    "reweighting_only_effect": float((w1 * m0).sum() - (w0 * m0).sum()),
    "global_activity_fall": float((m1 - m0).mean()),
    "per_celltype_activity_change": {c: float(m1[c] - m0[c]) for c in ADM},
    "_note": ("reweighting_only_effect freezes activity at Diagnosis and applies the Relapse "
              "mixture, so it isolates the mixture with no regulatory change. It is computed on "
              "cell-type means POOLED ACROSS PATIENTS, so it is not identical in construction to "
              "the per-patient decomposition in f6_decomposition.json and is not a substitute for "
              "it; it establishes the RELATIVE SIZE of the mixture effect only.")}
log(f"activity spread at Diagnosis {m0.min():.4f}-{m0.max():.4f} "
    f"({R['mixture_and_activity']['activity_spread_Diagnosis']['range_as_pct_of_mean']:.0f}% of mean)")
log(f"largest share change: {R['mixture_and_activity']['largest_share_change']['celltype']} "
    f"{R['mixture_and_activity']['largest_share_change']['from']:.4f} -> "
    f"{R['mixture_and_activity']['largest_share_change']['to']:.4f}")
log(f"reweighting alone moves the module {R['mixture_and_activity']['reweighting_only_effect']:+.5f} "
    f"against a global activity fall of {R['mixture_and_activity']['global_activity_fall']:+.5f}")

# ======================================================================== (2) arm enrichment
def s1_of(r):
    return r.split("_direct_")[1][0]

ARMS = {"activating (+)": [r for r in M55 if s1_of(r) == "+"],
        "repressing (-)": [r for r in M55 if s1_of(r) == "-"]}
log(f"arms: activating {len(ARMS['activating (+)'])} regulons, "
    f"repressing {len(ARMS['repressing (-)'])} regulons")
# the two arms of the WHOLE module can overlap (different TFs); only same-TF arms cannot
for v in VOCABS:
    a = set().union(*(TGT[v][r] for r in ARMS["activating (+)"]))
    b = set().union(*(TGT[v][r] for r in ARMS["repressing (-)"]))
    log(f"  {v}: activating union {len(a)} genes, repressing union {len(b)}, shared {len(a & b)}")

erows = []
for v in VOCABS:
    bg_full = set().union(*(TGT[v][r] for r in POOL)) & UNIV[v]
    for gsname, drop in (("full", set()), ("ribosome_depleted", RIBO)):
        bg = bg_full - drop
        for arm, regs in ARMS.items():
            fg = (set().union(*(TGT[v][r] for r in regs)) & UNIV[v]) - drop
            for lib in LIBS:
                rows = []
                for term, gs in GMT[lib].items():
                    t = gs & bg
                    if not (SIZE_MIN <= len(t) <= SIZE_MAX):
                        continue
                    ov = len(t & fg)
                    orr, p = stats.fisher_exact(
                        [[ov, len(fg) - ov],
                         [len(t) - ov, len(bg) - len(fg) - len(t) + ov]], alternative="greater")
                    rows.append({"vocabulary": v, "gene_set": gsname, "arm": arm, "library": lib,
                                 "term": term, "n_overlap": ov, "n_term_in_bg": len(t),
                                 "n_foreground": len(fg), "n_background": len(bg),
                                 "odds_ratio": orr, "p": p})
                d = pd.DataFrame(rows)
                if len(d):
                    d["q"] = bh(d.p.values)          # BH within (library x vocabulary x arm x set)
                    erows.append(d)
ENR = pd.concat(erows, ignore_index=True)
ENR.to_csv(os.path.join(OUT, "f6_fig06_arm_enrichment.csv"), index=False)
log(f"arm enrichment: {len(ENR)} rows over "
    f"{ENR.groupby(['vocabulary','gene_set','arm']).ngroups} (vocabulary x set x arm) cells")

# reproducible-across-vocabularies terms, on the FULL sets, per arm
full = ENR[ENR.gene_set == "full"]
REPRO = {}
for arm in ARMS:
    sig = {v: set(full.term[(full.vocabulary == v) & (full.arm == arm) & (full.q < Q_BAR)])
           for v in VOCABS}
    both = sorted(sig["Diagnosis"] & sig["Relapse"])
    REPRO[arm] = {"n_sig_Diagnosis": len(sig["Diagnosis"]), "n_sig_Relapse": len(sig["Relapse"]),
                  "n_reproducible_both": len(both), "terms": both}
    log(f"  {arm}: {len(sig['Diagnosis'])} sig Dx, {len(sig['Relapse'])} sig Rel, "
        f"{len(both)} reproducible in both")
R["arm_enrichment"] = {
    "arms": {k: {"n_regulons": len(v2)} for k, v2 in ARMS.items()},
    "background": "union of all 201 pool regulons' targets, intersected with the expressed universe",
    "background_saturation_note": ("the expressed-universe background saturates -- 44 of 50 "
                                   "Hallmark terms significant -- so it is not used and no term "
                                   "from it is quoted"),
    "size_window": [SIZE_MIN, SIZE_MAX], "bh_family": "library x vocabulary x arm x gene set",
    "reproducible_across_vocabularies": REPRO}

# survives the ribosome drop? per (arm, term): significant in BOTH vocabularies after depletion
dep = ENR[ENR.gene_set == "ribosome_depleted"]
SURV = {}
for arm in ARMS:
    sig = {v: set(dep.term[(dep.vocabulary == v) & (dep.arm == arm) & (dep.q < Q_BAR)])
           for v in VOCABS}
    SURV[arm] = sorted(sig["Diagnosis"] & sig["Relapse"])
    keep = set(REPRO[arm]["terms"]) & set(SURV[arm])
    log(f"  {arm}: of {len(REPRO[arm]['terms'])} reproducible terms, "
        f"{len(keep)} still reproducible after removing ribosomal protein genes")
R["survives_ribosome_drop"] = {
    "definition": ("a term SURVIVES if it is still significant at q < 0.05 in BOTH vocabularies "
                   "after the 167 ribosomal protein genes of f6_ribosome_check.json DROP_A are "
                   "removed from both foreground and background"),
    "n_ribosomal_genes_removed": len(RIBO),
    **{arm: {"n_surviving": len(set(REPRO[arm]["terms"]) & set(SURV[arm])),
             "surviving_terms": sorted(set(REPRO[arm]["terms"]) & set(SURV[arm]))}
       for arm in ARMS}}

json.dump(R, open(os.path.join(OUT, "f6_fig06_support.json"), "w", encoding="utf-8"),
          indent=1, default=str)
log("wrote f6_fig06_support.json, f6_fig06_mixture.csv, f6_fig06_arm_enrichment.csv")
