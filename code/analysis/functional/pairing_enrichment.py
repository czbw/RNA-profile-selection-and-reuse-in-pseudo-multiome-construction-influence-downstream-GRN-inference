"""FIG1 / S1.c - GO / pathway enrichment of each arm's target-gene union."""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
RNA = ROOT / "paml/malignant" / "data" / "rna"
REF = ROOT / "analysis" / "ref"
OUT = ROOT / "analysis" / "out"

STAGES = {"Diagnosis": "malignant_diagnosis.h5ad", "Relapse": "malignant_relapse.h5ad"}
ARMS = ["link", "pseudomultiome", "repair"]
GRID = [0.001, 0.005, 0.01, 0.05, 0.10]
RETENTION_BAR = 0.995
SIZE_MIN, SIZE_MAX = 10, 500
SIZE_SENS = [(5, 1000)]
Q_BAR = 0.05
LIBS = {"GOBP": "GO_Biological_Process_2025.gmt",
        "REAC": "Reactome_Pathways_2024.gmt",
        "HALL": "MSigDB_Hallmark_2020.gmt"}

def log(m):
    print(f"[f1-enrich] {m}", flush=True)

def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def detection(h5path):
    """per-gene fraction of cells with a nonzero entry, from the CSR indices."""
    with h5py.File(h5path, "r") as h:
        shp = h["X"].attrs["shape"]
        n_cells, n_genes = int(shp[0]), int(shp[1])
        idx = h["X/indices"]
        nnz = np.zeros(n_genes, dtype=np.int64)
        CH = 50_000_000
        for s in range(0, idx.shape[0], CH):
            nnz += np.bincount(idx[s:s + CH], minlength=n_genes).astype(np.int64)
        ens = h["var/_index"][:]
        s = h["var/symbol"]
        if isinstance(s, h5py.Group):          # anndata categorical encoding
            cats = s["categories"][:]
            codes = s["codes"][:]
            dec0 = [x.decode() if isinstance(x, bytes) else str(x) for x in cats]
            sym = [dec0[c] if c >= 0 else "" for c in codes]
        else:
            sym = s[:]
    dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    return pd.DataFrame({"ensembl": dec(ens), "symbol": dec(sym),
                         "n_cells_detected": nnz,
                         "detection": nnz / n_cells}), n_cells

def read_gmt(p):
    sets = {}
    for line in open(p, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        g = {x.strip() for x in f[2:] if x.strip()}
        if g:
            sets[f[0]] = g
    return sets

# ---------------------------------------------------------------- target unions
unions = {}
for stage in STAGES:
    for arm in ARMS:
        d = pd.read_csv(S2 / f"eRegulons_{stage}_{arm}.tsv", sep="\t",
                        usecols=["Gene"], low_memory=False)
        unions[(stage, arm)] = set(d.Gene.astype(str))
        log(f"{stage}/{arm} target union: {len(unions[(stage,arm)])} genes")

# ---------------------------------------------------------------- universe
res = {"status": "DESCRIPTIVE / EXPLORATORY",
       "prespec": "analysis/FIG1_FIG6_PRESPEC.md S1.c",
       "universe_rule": ("largest d in %s at which >= %.1f%% of BOTH arms' target genes at BOTH "
                         "timepoints are retained" % (GRID, RETENTION_BAR * 100)),
       "detection_grid": {}, "universe": {}, "libraries": {}, "counts": {},
       "top_divergent": {}, "top_per_arm": {}}

det, ncell = {}, {}
for stage, fn in STAGES.items():
    det[stage], ncell[stage] = detection(RNA / fn)
    log(f"{stage}: {ncell[stage]} cells, {len(det[stage])} genes in the RNA object")

# symbol-level detection = max over the Ensembl IDs mapping to that symbol
sym_det = {}
for stage in STAGES:
    sym_det[stage] = det[stage].groupby("symbol")["detection"].max()

# evaluate the grid: retention of every arm's target genes at every d, at both timepoints
grid_tab = []
for d in GRID:
    keep = {s: set(sym_det[s][sym_det[s] >= d].index) for s in STAGES}
    worst = 1.0
    for stage in STAGES:
        for arm in ARMS:
            u = unions[(stage, arm)]
            r = len(u & keep[stage]) / len(u)
            grid_tab.append({"d": d, "stage": stage, "arm": arm,
                             "n_target_genes": len(u), "n_retained": len(u & keep[stage]),
                             "retention": r, "universe_size": len(keep[stage])})
            if arm in ("link", "pseudomultiome"):
                worst = min(worst, r)
    res["detection_grid"][str(d)] = {"worst_retention_over_the_two_Fig1_arms": worst,
                                     "meets_bar": bool(worst >= RETENTION_BAR),
                                     "universe_size": {s: len(keep[s]) for s in STAGES}}
GRIDDF = pd.DataFrame(grid_tab)
GRIDDF.to_csv(OUT / "f1_universe_grid.csv", index=False)

ok = [d for d in GRID if res["detection_grid"][str(d)]["meets_bar"]]
D_CHOSEN = max(ok) if ok else min(GRID)
res["d_chosen"] = D_CHOSEN
res["d_chosen_note"] = ("no grid point met the retention bar; the smallest d was used and this is "
                        "reported as a deviation" if not ok else "rule satisfied")
log(f"d chosen = {D_CHOSEN} ({res['d_chosen_note']})")

ann = pd.read_csv(ROOT / "paml/malignant" / "databases" / "hg38_genome_annotation.tsv",
                  sep="\t", low_memory=False)
ann_col = "Gene" if "Gene" in ann.columns else ann.columns[0]
ANN_SYM = set(ann[ann_col].astype(str))
log(f"hg38 annotation: column '{ann_col}', {len(ANN_SYM)} distinct symbols")

U = {}
for stage in STAGES:
    up = set(sym_det[stage][sym_det[stage] >= D_CHOSEN].index)
    us = up & ANN_SYM
    U[stage] = {"primary": up, "sensitivity": us}
    res["universe"][stage] = {"d": D_CHOSEN, "n_primary": len(up), "n_sensitivity": len(us),
                              "n_genes_in_rna_object": int(len(sym_det[stage]))}
    pd.DataFrame({"symbol": sorted(up),
                  "in_hg38_annotation": [g in ANN_SYM for g in sorted(up)]}
                 ).to_csv(OUT / f"f1_universe_{stage}.csv", index=False)
    log(f"{stage}: U_primary={len(up)}  U_sensitivity={len(us)}")

# ---------------------------------------------------------------- enrichment
libs = {k: read_gmt(REF / v) for k, v in LIBS.items()}
for k, v in libs.items():
    res["libraries"][k] = {"file": LIBS[k], "n_sets": len(v)}

div_rows, arm_rows = [], []
for uname in ("primary", "sensitivity"):
    for smin, smax in [(SIZE_MIN, SIZE_MAX)] + SIZE_SENS:
        for stage in STAGES:
            uni = U[stage][uname]
            fg = {a: unions[(stage, a)] & uni for a in ARMS}
            if uname == "primary" and (smin, smax) == (SIZE_MIN, SIZE_MAX):
                res["counts"][stage] = {a: {"target_union": len(unions[(stage, a)]),
                                            "in_universe": len(fg[a]),
                                            "frac_of_universe": len(fg[a]) / len(uni)}
                                        for a in ARMS}
            for lib, sets in libs.items():
                pend_div, pend_arm = [], []
                for term, g in sets.items():
                    gu = g & uni
                    if not (smin <= len(gu) <= smax):
                        continue
                    # per-arm one-sided hypergeometric
                    for a in ARMS:
                        k = len(gu & fg[a])
                        p = stats.hypergeom.sf(k - 1, len(uni), len(gu), len(fg[a]))
                        exp = len(gu) * len(fg[a]) / len(uni)
                        pend_arm.append({"universe": uname, "size_window": f"{smin}-{smax}",
                                         "stage": stage, "library": lib, "arm": a, "term": term,
                                         "n_term_in_universe": len(gu), "n_arm_in_universe": len(fg[a]),
                                         "n_overlap": k, "expected": exp,
                                         "fold_enrichment": (k / exp) if exp > 0 else np.nan,
                                         "p_one_sided": float(p)})
                    # between-arm, two-sided
                    for A, B in [("link", "pseudomultiome"), ("link", "repair"),
                                 ("pseudomultiome", "repair")]:
                        a1 = len(gu & fg[A]); a0 = len(fg[A]) - a1
                        b1 = len(gu & fg[B]); b0 = len(fg[B]) - b1
                        odds, p = stats.fisher_exact([[a1, a0], [b1, b0]])
                        pend_div.append({"universe": uname, "size_window": f"{smin}-{smax}",
                                         "stage": stage, "library": lib,
                                         "arm_A": A, "arm_B": B, "term": term,
                                         "n_term_in_universe": len(gu),
                                         "n_in_A": a1, "n_A_total": len(fg[A]),
                                         "n_in_B": b1, "n_B_total": len(fg[B]),
                                         "frac_A": a1 / len(fg[A]), "frac_B": b1 / len(fg[B]),
                                         "odds_ratio_A_over_B": float(odds),
                                         "p_two_sided": float(p)})
                # BH within (library x timepoint) for divergence, per arm-pair
                dv = pd.DataFrame(pend_div)
                if len(dv):
                    dv["bh_q"] = np.nan
                    for (A, B), sub in dv.groupby(["arm_A", "arm_B"]):
                        dv.loc[sub.index, "bh_q"] = bh(sub.p_two_sided.values)
                    div_rows.append(dv)
                ar = pd.DataFrame(pend_arm)
                if len(ar):
                    ar["bh_q"] = np.nan
                    for a, sub in ar.groupby("arm"):
                        ar.loc[sub.index, "bh_q"] = bh(sub.p_one_sided.values)
                    arm_rows.append(ar)
                log(f"{uname}/{smin}-{smax}/{stage}/{lib}: {len(sets)} sets -> "
                    f"{len(pend_div)//3} testable")

DIV = pd.concat(div_rows, ignore_index=True)
ARM = pd.concat(arm_rows, ignore_index=True)
DIV["significant"] = DIV.bh_q < Q_BAR
ARM["significant"] = ARM.bh_q < Q_BAR
DIV.to_csv(OUT / "f1_enrichment_divergence.csv", index=False)
ARM.to_csv(OUT / "f1_enrichment_per_arm.csv", index=False)

P = (DIV.universe == "primary") & (DIV.size_window == f"{SIZE_MIN}-{SIZE_MAX}")
for stage in STAGES:
    for lib in LIBS:
        s = DIV[P & (DIV.stage == stage) & (DIV.library == lib)
                & (DIV.arm_A == "link") & (DIV.arm_B == "pseudomultiome")]
        res["top_divergent"][f"{stage}_{lib}"] = {
            "n_terms_tested": int(len(s)),
            "n_significant_q<0.05": int(s.significant.sum()),
            "n_sig_higher_in_link": int(((s.bh_q < Q_BAR) & (s.odds_ratio_A_over_B > 1)).sum()),
            "n_sig_higher_in_pseudomultiome": int(((s.bh_q < Q_BAR) & (s.odds_ratio_A_over_B < 1)).sum()),
            "top20": json.loads(s.nsmallest(20, "p_two_sided")[
                ["term", "n_term_in_universe", "n_in_A", "frac_A", "n_in_B", "frac_B",
                 "odds_ratio_A_over_B", "p_two_sided", "bh_q"]].to_json(orient="records")),
        }
PA = (ARM.universe == "primary") & (ARM.size_window == f"{SIZE_MIN}-{SIZE_MAX}")
for stage in STAGES:
    for arm in ARMS:
        s = ARM[PA & (ARM.stage == stage) & (ARM.library == "GOBP") & (ARM.arm == arm)]
        res["top_per_arm"][f"{stage}_{arm}_GOBP"] = {
            "n_terms_tested": int(len(s)),
            "n_significant_q<0.05": int(s.significant.sum()),
            "top15": json.loads(s.nsmallest(15, "p_one_sided")[
                ["term", "n_term_in_universe", "n_overlap", "fold_enrichment",
                 "p_one_sided", "bh_q"]].to_json(orient="records")),
        }

# "significant in one arm, not the other" - reported, and labelled
sig_only = []
for stage in STAGES:
    for lib in LIBS:
        a = ARM[PA & (ARM.stage == stage) & (ARM.library == lib) & (ARM.arm == "link")
                ].set_index("term")
        b = ARM[PA & (ARM.stage == stage) & (ARM.library == lib) & (ARM.arm == "pseudomultiome")
                ].set_index("term")
        common = a.index.intersection(b.index)
        for t in common:
            sa, sb = bool(a.loc[t, "significant"]), bool(b.loc[t, "significant"])
            if sa != sb:
                sig_only.append({"stage": stage, "library": lib, "term": t,
                                 "sig_in_link": sa, "sig_in_pseudomultiome": sb,
                                 "q_link": float(a.loc[t, "bh_q"]),
                                 "q_pseudomultiome": float(b.loc[t, "bh_q"]),
                                 "fold_link": float(a.loc[t, "fold_enrichment"]),
                                 "fold_pseudomultiome": float(b.loc[t, "fold_enrichment"])})
SO = pd.DataFrame(sig_only)
SO["CAVEAT"] = "difference in significance, NOT a significant difference - see the divergence table"
SO.to_csv(OUT / "f1_enrichment_sig_in_one_arm_only.csv", index=False)
res["n_sig_in_one_arm_only"] = {f"{s}_{l}": int(len(SO[(SO.stage == s) & (SO.library == l)]))
                                for s in STAGES for l in LIBS}

with open(OUT / "f1_enrichment.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote f1_enrichment.json + 5 csv companions")
