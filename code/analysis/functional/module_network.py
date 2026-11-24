#!/usr/bin/env python
"""F6 network analysis of the corrected arm's 55-regulon module. Descriptive and exploratory."""
import itertools
import json
import os
from collections import Counter, defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "analysis", "out")
REF = os.path.join(ROOT, "analysis", "ref")
S2 = os.path.join(ROOT, "paml/malignant", "out")
MODJ = os.path.join(ROOT, "modules", "out", "corrected_module.json")

VOCABS = ("Diagnosis", "Relapse")
SEED = 555
N_DRAWS = 1000
SIZE_MIN, SIZE_MAX = 10, 500
SIZE_SENS = (5, 1000)
Q_BAR = 0.05
TOP_N = 15
LIBS = {"GOBP": "GO_Biological_Process_2025.gmt",
        "REAC": "Reactome_Pathways_2024.gmt",
        "HALL": "MSigDB_Hallmark_2020.gmt"}
ALL_MODE = ["IRF8", "KLF10", "SPI1", "ZBTB7A"]

R = {"_status": "DESCRIPTIVE / EXPLORATORY; inherits the post-hoc status of the 55",
     "_prespec": "analysis/F6_NETWORK_PRESPEC.md",
     "_seed": SEED, "_n_draws": N_DRAWS}

def log(m):
    print(f"[f6-net] {m}", flush=True)

def hdr(m):
    print(f"\n{'=' * 78}\n{m}\n{'=' * 78}", flush=True)

def bh(p):
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def clean(gs):
    return {str(g).strip() for g in gs if str(g).strip()}

def jac(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0

def parse_mode(reg):
    """('SPI1', '+', '+') from 'SPI1_direct_+/+'."""
    tf, mode = reg.split("_direct_")
    return tf, mode[0], mode[2]

# ======================================================================== substrate
MOD = json.load(open(MODJ, encoding="utf-8"))
M55 = list(MOD["regulons"])
POOL = list(MOD["background_regulons"])
assert len(M55) == 55 and len(POOL) == 201 and set(M55) <= set(POOL)
TGT = {v: {r: clean(MOD["targets"][v][r]) for r in POOL} for v in VOCABS}
MTF = sorted({parse_mode(r)[0] for r in M55})
assert len(MTF) == 31
log(f"module {len(M55)} regulons / {len(MTF)} TFs; pool {len(POOL)}; source {os.path.relpath(MODJ, ROOT)}")

UNIV = {}
for v in VOCABS:
    u = pd.read_csv(os.path.join(OUT, f"f1_universe_{v}.csv"))
    u = u[u.symbol.notna() & (u.symbol.astype(str).str.strip() != "")]
    UNIV[v] = {"all": clean(u.symbol), "hg38": clean(u.symbol[u.in_hg38_annotation])}
    log(f"universe {v}: {len(UNIV[v]['all'])} genes (hg38 subset {len(UNIV[v]['hg38'])})")
R["_universe"] = {v: {"n_all": len(UNIV[v]["all"]), "n_hg38": len(UNIV[v]["hg38"]),
                      "source": f"analysis/out/f1_universe_{v}.csv"} for v in VOCABS}

# ======================================================================== check 0: sign-mode convention
hdr("CHECK 0 -- the sign-mode convention, verified not assumed")
g0 = {"rule": "zero off-diagonal cells in both cross-tabulations, both timepoints"}
g0_pass = True
for v in VOCABS:
    d = pd.read_csv(os.path.join(S2, f"eRegulons_{v}_repair.tsv"), sep="\t",
                    usecols=["eRegulon_name", "rho_TF2G", "rho_R2G", "TF", "Gene"],
                    low_memory=False)
    m = d.eRegulon_name.str.extract(r"_([+-])/([+-])$")
    s1o = np.where(d.rho_TF2G > 0, "+", "-")
    s2o = np.where(d.rho_R2G > 0, "+", "-")
    off1 = int((m[0].values != s1o).sum())
    off2 = int((m[1].values != s2o).sum())
    g0[v] = {"n_triplets": len(d), "off_diagonal_s1_vs_rho_TF2G": off1,
             "off_diagonal_s2_vs_rho_R2G": off2}
    g0_pass &= (off1 == 0 and off2 == 0)
    log(f"{v}: {len(d):,} triplets | s1 vs sign(rho_TF2G) mismatches {off1} | "
        f"s2 vs sign(rho_R2G) mismatches {off2}")
    if v == "Diagnosis":
        # check 1 needs the (TF, gene) -> distinct s1 count over the whole pool
        d["s1"] = m[0]
        g1c = d.groupby(["TF", "Gene"]).s1.nunique()
        G1 = {"n_TF_gene_pairs": int(len(g1c)),
              "n_pairs_in_2plus_regulons": int((d.groupby(["TF", "Gene"])
                                                .eRegulon_name.nunique() >= 2).sum()),
              "max_distinct_s1_per_TF_gene_pair": int(g1c.max()),
              "n_pairs_with_2_distinct_s1": int((g1c >= 2).sum())}
g0["PASS"] = bool(g0_pass)
R["G0_sign_mode_convention"] = g0
if not g0_pass:
    log("CHECK 0 FAILED -- Part 2 does not run; the failure is the result")
    json.dump(R, open(os.path.join(OUT, "f6_network.json"), "w", encoding="utf-8"), indent=1)
    raise SystemExit(1)
log("CHECK 0 PASS -- s1 = sign(rho_TF2G) [TF->gene], s2 = sign(rho_R2G) [region->gene]")

hdr("CHECK 1 -- is s1 single-valued per (TF, gene)?")
G1["structurally_determined"] = bool(G1["max_distinct_s1_per_TF_gene_pair"] == 1)
R["G1_s1_single_valued_per_TF_gene"] = G1
log(f"(TF,gene) pairs {G1['n_TF_gene_pairs']:,}; in >=2 regulons "
    f"{G1['n_pairs_in_2plus_regulons']:,}; max distinct s1 "
    f"{G1['max_distinct_s1_per_TF_gene_pair']}; pairs with 2 distinct s1 "
    f"{G1['n_pairs_with_2_distinct_s1']}")
log("=> stratum (A) primary result is STRUCTURALLY DETERMINED"
    if G1["structurally_determined"] else "=> stratum (A) is informative under both derivations")

# ======================================================================== Part 1: target-set overlap
def pair_jaccards(regs, v):
    return np.array([jac(TGT[v][a], TGT[v][b]) for a, b in itertools.combinations(regs, 2)])

def jac_graph(regs, v):
    g = nx.Graph()
    g.add_nodes_from(regs)
    for a, b in itertools.combinations(regs, 2):
        w = jac(TGT[v][a], TGT[v][b])
        if w > 0:
            g.add_edge(a, b, weight=w)
    return g

def modularity_of(regs, v, method="greedy"):
    g = jac_graph(regs, v)
    if g.number_of_edges() == 0:
        return 0.0, [{r} for r in regs]
    if method == "greedy":
        comm = list(nx.community.greedy_modularity_communities(g, weight="weight"))
    else:
        comm = list(nx.community.louvain_communities(g, weight="weight", seed=SEED))
    comm = [set(c) for c in comm]
    iso = set(regs) - set().union(*comm) if comm else set(regs)
    comm += [{r} for r in iso]
    return nx.community.modularity(g, comm, weight="weight"), comm

hdr("PART 1 TARGET-SET OVERLAP ACROSS THE 55")
rng = np.random.default_rng(SEED)
DRAWS = [list(rng.choice(POOL, size=55, replace=False)) for _ in range(N_DRAWS)]

# size-matched draws, per vocabulary (bin edges depend on the vocabulary's target sizes)
DRAWS_SM = {}
for v in VOCABS:
    sz = pd.Series({r: len(TGT[v][r]) for r in POOL})
    binid = pd.qcut(sz.rank(method="first"), 5, labels=False)
    want = Counter(binid[M55])
    rng2 = np.random.default_rng(SEED + 1)
    dd = []
    for _ in range(N_DRAWS):
        pick = []
        for b, n in want.items():
            cand = list(binid.index[binid == b])
            pick += list(rng2.choice(cand, size=n, replace=False))
        dd.append(pick)
    DRAWS_SM[v] = dd
    log(f"size-matched null {v}: quintile counts in the module {dict(sorted(want.items()))}")

S1, PJ_ROWS = {}, []
for v in VOCABS:
    pj = pair_jaccards(M55, v)
    for (a, b), w in zip(itertools.combinations(M55, 2), pj):
        PJ_ROWS.append({"vocabulary": v, "regulon_a": a, "regulon_b": b, "jaccard": w,
                        "same_TF": parse_mode(a)[0] == parse_mode(b)[0]})
    union = set().union(*(TGT[v][r] for r in M55))
    hit = Counter()
    for r in M55:
        hit.update(TGT[v][r])
    sizes = [len(TGT[v][r]) for r in M55]
    Qg, Cg = modularity_of(M55, v, "greedy")
    Ql, Cl = modularity_of(M55, v, "louvain")
    nullQ = np.array([modularity_of(d, v, "greedy")[0] for d in DRAWS])
    p_emp = (1 + int((nullQ >= Qg).sum())) / (1 + N_DRAWS)
    big = sum(1 for c in Cg if len(c) >= 3)
    blocked = bool(p_emp < 0.05 and big >= 2)
    nullQl = np.array([modularity_of(d, v, "louvain")[0] for d in DRAWS[:200]])
    p_l = (1 + int((nullQl >= Ql).sum())) / (1 + len(nullQl))
    S1[v] = {
        "n_pairs": int(len(pj)),
        "jaccard": {"min": float(pj.min()), "p05": float(np.percentile(pj, 5)),
                    "p25": float(np.percentile(pj, 25)), "median": float(np.median(pj)),
                    "p75": float(np.percentile(pj, 75)), "p95": float(np.percentile(pj, 95)),
                    "max": float(pj.max()), "mean": float(pj.mean()),
                    "frac_pairs_zero": float((pj == 0).mean())},
        "target_union_size": len(union),
        "targets_per_regulon": {"min": int(min(sizes)), "median": float(np.median(sizes)),
                                "max": int(max(sizes))},
        "genes_hit_by_ge2": int(sum(1 for g, c in hit.items() if c >= 2)),
        "genes_hit_by_ge3": int(sum(1 for g, c in hit.items() if c >= 3)),
        "genes_hit_by_ge5": int(sum(1 for g, c in hit.items() if c >= 5)),
        "modularity_greedy": {
            "Q_observed": float(Qg), "n_communities": len(Cg),
            "n_communities_ge3": big,
            "null_mean": float(nullQ.mean()), "null_sd": float(nullQ.std(ddof=1)),
            "null_p95": float(np.percentile(nullQ, 95)),
            "p_empirical_one_sided": float(p_emp),
            "percentile_of_observed": float((nullQ < Qg).mean() * 100),
            "communities": [sorted(c) for c in sorted(Cg, key=len, reverse=True)]},
        "modularity_louvain_sensitivity": {
            "Q_observed": float(Ql), "n_communities": len(Cl),
            "n_draws": int(len(nullQl)), "p_empirical_one_sided": float(p_l)},
        "VERDICT": "BLOCK-STRUCTURED" if blocked else "FLAT"}
    log(f"{v}: union {len(union)}, median targets {np.median(sizes):.0f}, "
        f"Jaccard median {np.median(pj):.4f} (max {pj.max():.4f}, {100*(pj==0).mean():.1f}% zero); "
        f"Q {Qg:.4f} vs null {nullQ.mean():.4f}+-{nullQ.std(ddof=1):.4f}, p {p_emp:.4f} "
        f"-> {S1[v]['VERDICT']}")
# ---- block detail: within/between Jaccard and family composition (pre-specified) ----------
S55F = pd.read_csv(os.path.join(OUT, "f6_set55_per_regulon.csv")).set_index("regulon")
for v in VOCABS:
    comms = S1[v]["modularity_greedy"]["communities"]
    cid = {r: i for i, c in enumerate(comms) for r in c}
    win = [jac(TGT[v][a], TGT[v][b]) for a, b in itertools.combinations(M55, 2)
           if cid[a] == cid[b]]
    bet = [jac(TGT[v][a], TGT[v][b]) for a, b in itertools.combinations(M55, 2)
           if cid[a] != cid[b]]
    uw, pw = stats.mannwhitneyu(win, bet, alternative="greater")
    S1[v]["block_detail"] = {
        "mean_within_block_jaccard": float(np.mean(win)),
        "mean_between_block_jaccard": float(np.mean(bet)),
        "median_within_block_jaccard": float(np.median(win)),
        "median_between_block_jaccard": float(np.median(bet)),
        "mannwhitney_within_gt_between_p": float(pw),
        "blocks": [{"block": i, "n_regulons": len(c),
                    "TFs": sorted({parse_mode(r)[0] for r in c}),
                    "n_TFs": len({parse_mode(r)[0] for r in c}),
                    "level_B_families": dict(Counter(S55F.loc[sorted(c)].level_B_family)),
                    "regulons": sorted(c)} for i, c in enumerate(comms)],
        # Post hoc, not pre-specified: the plan covered size, members, TF composition, families
        # and within/between Jaccard; the two cross-tabs below were added after seeing that
        # the blocks do not partition by TF. Not used for any verdict.
        "POSTHOC_not_prespecified": {
            "n_TFs_appearing_in_more_than_one_block":
                int(sum(1 for tf in MTF
                        if len({cid[r] for r in M55 if parse_mode(r)[0] == tf}) > 1)),
            "TFs_in_more_than_one_block":
                sorted(tf for tf in MTF
                       if len({cid[r] for r in M55 if parse_mode(r)[0] == tf}) > 1),
            "block_by_s1": {f"block{i}": dict(Counter(parse_mode(r)[1] for r in c))
                            for i, c in enumerate(comms)},
            "block_by_s2": {f"block{i}": dict(Counter(parse_mode(r)[2] for r in c))
                            for i, c in enumerate(comms)},
            "block_by_mode": {f"block{i}": dict(Counter(f"{parse_mode(r)[1]}/{parse_mode(r)[2]}"
                                                        for r in c))
                              for i, c in enumerate(comms)}}}
    log(f"{v}: within-block J {np.mean(win):.4f} vs between {np.mean(bet):.4f} "
        f"(MW greater p {pw:.3g}); TFs in >1 block "
        f"{S1[v]['block_detail']['POSTHOC_not_prespecified']['n_TFs_appearing_in_more_than_one_block']}"
        f" of {len(MTF)}")
    log(f"    [POST-HOC] block x s2: {S1[v]['block_detail']['POSTHOC_not_prespecified']['block_by_s2']}")

# ---- Post hoc, not pre-specified: does the block partition reproduce across the two
# ---- vocabularies, and where does each TF sit? Added after seeing that the Part 1 modularity
# ---- verdict clears its null only narrowly while per-block enrichment separates by
# ---- orders of magnitude. Not used for any verdict.
CID = {v: {r: i for i, c in enumerate(S1[v]["modularity_greedy"]["communities"]) for r in c}
       for v in VOCABS}
bx = {}
for i, ci in enumerate(S1["Diagnosis"]["modularity_greedy"]["communities"]):
    for j, cj in enumerate(S1["Relapse"]["modularity_greedy"]["communities"]):
        bx[f"Dx_block{i}|Rel_block{j}"] = {"n_shared": len(set(ci) & set(cj)),
                                           "jaccard": float(jac(set(ci), set(cj)))}
agree = sum(1 for r in M55 if CID["Diagnosis"][r] == CID["Relapse"][r])
R["S1_POSTHOC_block_reproducibility"] = {
    "_flag": "POST-HOC, NOT IN THE PRE-SPEC. Reported, used for no verdict.",
    "cross_vocabulary_block_overlap": bx,
    "n_regulons_same_block_index_in_both": int(agree),
    "frac_regulons_same_block_index": float(agree / 55),
    "per_TF_block_assignment": {
        tf: {v: sorted({CID[v][r] for r in M55 if parse_mode(r)[0] == tf}) for v in VOCABS}
        for tf in MTF}}
log(f"[POST-HOC] block membership agrees across vocabularies for {agree}/55 regulons "
    f"({agree / 55:.1%}); Dx0|Rel0 J "
    f"{bx['Dx_block0|Rel_block0']['jaccard']:.4f}, Dx1|Rel1 J "
    f"{bx['Dx_block1|Rel_block1']['jaccard']:.4f}")

R["S1_target_overlap"] = S1
pd.DataFrame(PJ_ROWS).to_csv(os.path.join(OUT, "f6_network_pairwise_jaccard.csv"), index=False)

# ======================================================================== Part 2: regulatory concordance
def net_effects(regs, v, deriv):
    """gene -> list of net effects, and gene -> list of (regulon, TF, s1, s2)."""
    eff = defaultdict(list)
    prov = defaultdict(list)
    for r in regs:
        tf, s1, s2 = parse_mode(r)
        n = s1 if deriv == "primary" else ("+" if s1 == s2 else "-")
        for g in TGT[v][r]:
            eff[g].append(n)
            prov[g].append((r, tf, s1, s2))
    return eff, prov

def concordance(regs, v, deriv, stratum="all"):
    eff, prov = net_effects(regs, v, deriv)
    kept = {}
    for g, e in eff.items():
        p = prov[g]
        if stratum == "sameTF":
            tfs = Counter(t for _, t, _, _ in p)
            idx = [i for i, (_, t, _, _) in enumerate(p) if tfs[t] >= 2]
        elif stratum == "diffTF":
            seen = {}
            idx = [i for i, (_, t, _, _) in enumerate(p) if seen.setdefault(t, i) == i]
            if len({p[i][1] for i in idx}) < 2:
                idx = []
        else:
            idx = list(range(len(p)))
        if len(idx) >= 2:
            kept[g] = [e[i] for i in idx]
    if not kept:
        return {"n_genes": 0, "frac_concordant": float("nan"), "frac_opposed": float("nan")}
    conc = sum(1 for e in kept.values() if len(set(e)) == 1)
    byk = {}
    for lab, lo, hi in (("k2", 2, 2), ("k3", 3, 3), ("k4", 4, 4), ("k5", 5, 5),
                        ("k6_10", 6, 10), ("k_gt10", 11, 10 ** 9)):
        s = [e for e in kept.values() if lo <= len(e) <= hi]
        byk[lab] = {"n_genes": len(s),
                    "frac_concordant": (sum(1 for e in s if len(set(e)) == 1) / len(s))
                    if s else float("nan")}
    return {"n_genes": len(kept), "frac_concordant": conc / len(kept),
            "frac_opposed": 1 - conc / len(kept), "by_k": byk}

hdr("PART 2 REGULATORY CONCORDANCE ON SHARED TARGETS -- the primary analysis")
S2R = {}
for v in VOCABS:
    S2R[v] = {}
    for deriv in ("primary", "secondary"):
        S2R[v][deriv] = {}
        for stratum in ("all", "sameTF", "diffTF"):
            obs = concordance(M55, v, deriv, stratum)
            nulls = {}
            for nm, draws in (("unmatched", DRAWS), ("size_matched", DRAWS_SM[v])):
                nv = np.array([concordance(d, v, deriv, stratum)["frac_concordant"]
                               for d in draws], float)
                nv = nv[~np.isnan(nv)]
                o = obs["frac_concordant"]
                nulls[nm] = {
                    "n_draws_usable": int(nv.size),
                    "mean": float(nv.mean()), "sd": float(nv.std(ddof=1)),
                    "p025": float(np.percentile(nv, 2.5)),
                    "p975": float(np.percentile(nv, 97.5)),
                    "percentile_of_observed": float((nv < o).mean() * 100),
                    "p_empirical_one_sided_greater":
                        float((1 + int((nv >= o).sum())) / (1 + nv.size)),
                    "observed_exceeds_p975": bool(o > np.percentile(nv, 97.5))}
            S2R[v][deriv][stratum] = {"observed": obs, "nulls": nulls}
            log(f"{v:9s} {deriv:9s} {stratum:7s}: {obs['n_genes']:5d} genes, "
                f"concordant {obs['frac_concordant']:.4f} | unmatched null "
                f"{nulls['unmatched']['mean']:.4f} (p {nulls['unmatched']['p_empirical_one_sided_greater']:.4f}) | "
                f"size-matched {nulls['size_matched']['mean']:.4f} "
                f"(p {nulls['size_matched']['p_empirical_one_sided_greater']:.4f})")
    if G1["structurally_determined"]:
        S2R[v]["primary"]["sameTF"]["NOTE"] = (
            "STRUCTURALLY DETERMINED: check 1 shows s1 is single-valued per (TF, gene), so same-TF "
            "co-targets cannot disagree under the primary derivation. Not evidence.")

# the pre-specified verdict rule
VERD = all(S2R[v]["primary"]["all"]["nulls"][n]["observed_exceeds_p975"]
           for v in VOCABS for n in ("unmatched", "size_matched"))
R["S2_concordance"] = S2R
R["S2_VERDICT"] = {"rule": "pre-specified: observed > 97.5th pct of BOTH nulls, PRIMARY derivation, "
                           "BOTH vocabularies",
                   "VERDICT": "ABOVE CHANCE" if VERD else "AT CHANCE"}
log(f"\nPART 2 VERDICT (pre-specified rule): {R['S2_VERDICT']['VERDICT']}")

# raw mode-pair cross-tabulation, so the derivation can be checked
XT = {}
for v in VOCABS:
    c = Counter()
    _, prov = net_effects(M55, v, "primary")
    for g, p in prov.items():
        if len(p) >= 2:
            for (ra, _, a1, a2), (rb, _, b1, b2) in itertools.combinations(p, 2):
                c[(f"{a1}/{a2}", f"{b1}/{b2}")] += 1
    modes = ["+/+", "+/-", "-/+", "-/-"]
    XT[v] = {f"{a}|{b}": int(c[(a, b)] + (c[(b, a)] if a != b else 0))
             for i, a in enumerate(modes) for b in modes[i:]}
R["S2_mode_pair_crosstab"] = XT

ST_ROWS = []
for v in VOCABS:
    e1, prov = net_effects(M55, v, "primary")
    e2, _ = net_effects(M55, v, "secondary")
    for g, p in prov.items():
        if len(p) >= 2:
            ST_ROWS.append({
                "vocabulary": v, "gene": g, "k_regulons": len(p),
                "n_distinct_TFs": len({t for _, t, _, _ in p}),
                "regulons": ";".join(r for r, _, _, _ in p),
                "net_primary": ";".join(e1[g]), "net_secondary": ";".join(e2[g]),
                "concordant_primary": len(set(e1[g])) == 1,
                "concordant_secondary": len(set(e2[g])) == 1})
pd.DataFrame(ST_ROWS).to_csv(os.path.join(OUT, "f6_network_shared_targets.csv"), index=False)

# ======================================================================== Part 3: regulatory network
hdr("PART 3 THE REGULATORY NETWORK")
S3 = {}
TF_TGT = {v: {tf: set().union(*(TGT[v][r] for r in M55 if parse_mode(r)[0] == tf))
              for tf in MTF} for v in VOCABS}
for v in VOCABS:
    B = nx.Graph()
    for r in M55:
        for g in TGT[v][r]:
            B.add_edge(("R", r), ("G", g))
    comps = sorted((len(c) for c in nx.connected_components(B)), reverse=True)
    frac_lcc = comps[0] / B.number_of_nodes()
    gj = jac_graph(M55, v)
    T = nx.Graph()
    T.add_nodes_from(MTF)
    for a, b in itertools.combinations(MTF, 2):
        w = jac(TF_TGT[v][a], TF_TGT[v][b])
        if w > 0:
            T.add_edge(a, b, weight=w)
    deg = pd.Series({tf: len(TF_TGT[v][tf]) for tf in MTF}).sort_values(ascending=False)
    cut = float(np.percentile(deg.values, 90))
    hubs = sorted(deg.index[deg.values >= cut])
    gdeg = Counter()
    for r in M55:
        gdeg.update(TGT[v][r])
    gd = np.array(list(gdeg.values()))
    # same-TF modes: share or partition?
    pj = pd.DataFrame(PJ_ROWS)
    pj = pj[pj.vocabulary == v]
    same = pj.jaccard[pj.same_TF].values
    cross = pj.jaccard[~pj.same_TF].values
    u, pu = stats.mannwhitneyu(same, cross, alternative="two-sided")
    S3[v] = {
        "bipartite": {"n_nodes": B.number_of_nodes(), "n_edges": B.number_of_edges(),
                      "n_components": len(comps), "component_sizes_top5": comps[:5],
                      "frac_nodes_in_largest": float(frac_lcc),
                      "VERDICT": "ONE CONNECTED PROGRAMME" if frac_lcc >= 0.90
                                 else "NOT ONE COMPONENT"},
        "regulon_projection": {"n_edges": gj.number_of_edges(),
                               "n_components": nx.number_connected_components(gj),
                               "density": float(nx.density(gj))},
        "TF_projection": {"n_nodes": T.number_of_nodes(), "n_edges": T.number_of_edges(),
                          "n_components": nx.number_connected_components(T),
                          "density": float(nx.density(T))},
        "TF_degree": {"definition": "distinct module target genes",
                      "min": int(deg.min()), "median": float(deg.median()),
                      "max": int(deg.max()),
                      "hub_cut_p90": cut, "hubs": hubs,
                      "per_TF": {k: int(x) for k, x in deg.items()}},
        "gene_degree": {"min": int(gd.min()), "median": float(np.median(gd)),
                        "max": int(gd.max()), "skew": float(stats.skew(gd))},
        "all_mode_survivors": {tf: {"degree": int(deg[tf]),
                                    "rank_of_31": int(deg.rank(ascending=False)[tf]),
                                    "percentile": float(100 * (deg < deg[tf]).mean()),
                                    "is_hub": tf in hubs} for tf in ALL_MODE},
        "same_TF_modes_share_or_partition": {
            "n_same_TF_pairs": int(same.size), "n_cross_TF_pairs": int(cross.size),
            "median_jaccard_same_TF": float(np.median(same)),
            "median_jaccard_cross_TF": float(np.median(cross)),
            "mannwhitney_U": float(u), "p_two_sided": float(pu),
            "VERDICT": ("SHARE" if pu < 0.05 and np.median(same) > np.median(cross)
                        else "PARTITION" if pu < 0.05 else "NEITHER"),
            # pre-specified: reported for all such TFs and separately for the four all-mode
            # survivors. Per-TF pairwise Jaccard between that TF's own modes in the module.
            "per_TF": {tf: {"n_modes": int(len(md)),
                            "pairwise_jaccard": [float(jac(TGT[v][a], TGT[v][b]))
                                                 for a, b in itertools.combinations(md, 2)],
                            "median_jaccard": float(np.median(
                                [jac(TGT[v][a], TGT[v][b])
                                 for a, b in itertools.combinations(md, 2)])),
                            "max_jaccard": float(max(
                                jac(TGT[v][a], TGT[v][b])
                                for a, b in itertools.combinations(md, 2))),
                            "is_all_mode_survivor": tf in ALL_MODE}
                       for tf in MTF
                       for md in [sorted(r for r in M55 if parse_mode(r)[0] == tf)]
                       if len(md) >= 2}}}
    log(f"{v}: bipartite {B.number_of_nodes()} nodes / {len(comps)} components, "
        f"LCC {frac_lcc:.4f} -> {S3[v]['bipartite']['VERDICT']}; TF degree median "
        f"{deg.median():.0f}, hubs (p90 >= {cut:.0f}) {hubs}")
    log(f"    all-mode survivors: " + ", ".join(
        f"{t} deg {deg[t]} rank {int(deg.rank(ascending=False)[t])}/31"
        f"{' HUB' if t in hubs else ''}" for t in ALL_MODE))
    log(f"    same-TF modes: median J {np.median(same):.4f} vs cross-TF "
        f"{np.median(cross):.4f}, MW p {pu:.3g} -> "
        f"{S3[v]['same_TF_modes_share_or_partition']['VERDICT']}")
    _pt = S3[v]["same_TF_modes_share_or_partition"]["per_TF"]
    log("    all-mode survivors, own-mode Jaccard: " + ", ".join(
        f"{tf} n={_pt[tf]['n_modes']} med {_pt[tf]['median_jaccard']:.4f} "
        f"max {_pt[tf]['max_jaccard']:.4f}" for tf in ALL_MODE if tf in _pt))
R["S3_network"] = S3

# ---- emitted graph. No figure is built.
er, et, ett, nodes = [], [], [], []
S55T = pd.read_csv(os.path.join(OUT, "f6_set55_per_regulon.csv")).set_index("regulon")
for v in VOCABS:
    comm = {r: i for i, c in enumerate(S1[v]["modularity_greedy"]["communities"]) for r in c}
    for r in M55:
        tf, s1, s2 = parse_mode(r)
        for g in TGT[v][r]:
            er.append({"vocabulary": v, "regulon": r, "TF": tf, "gene": g, "s1": s1, "s2": s2,
                       "net_primary": s1, "net_secondary": "+" if s1 == s2 else "-"})
        row = S55T.loc[r]
        nodes.append({"vocabulary": v, "node": r, "node_type": "regulon", "TF": tf,
                      "sign_mode": f"{s1}/{s2}", "degree": len(TGT[v][r]),
                      "community_S1": comm.get(r), "level_B_family": row.level_B_family,
                      "ax4_would_remove": bool(row.ax4_would_remove),
                      "abs_beta_max": float(row.abs_beta_max), "direction": row.direction,
                      "external_delta": float(row[f"external_delta_{v}"])})
    for tf in MTF:
        for g in TF_TGT[v][tf]:
            et.append({"vocabulary": v, "TF": tf, "gene": g})
        nodes.append({"vocabulary": v, "node": tf, "node_type": "TF", "TF": tf,
                      "sign_mode": "", "degree": len(TF_TGT[v][tf]),
                      "community_S1": None,
                      "level_B_family": S55T[S55T.TF == tf].level_B_family.iloc[0],
                      "ax4_would_remove": None, "abs_beta_max": None, "direction": "",
                      "external_delta": None})
    for a, b in itertools.combinations(MTF, 2):
        w = jac(TF_TGT[v][a], TF_TGT[v][b])
        if w > 0:
            ett.append({"vocabulary": v, "TF_a": a, "TF_b": b, "jaccard": w,
                        "n_shared": len(TF_TGT[v][a] & TF_TGT[v][b])})
for nm, rows in (("edges_regulon", er), ("edges_TF", et), ("edges_TFTF", ett), ("nodes", nodes)):
    pd.DataFrame(rows).to_csv(os.path.join(OUT, f"f6_network_{nm}.csv"), index=False)
log(f"emitted graph: {len(er):,} regulon-gene, {len(et):,} TF-gene, {len(ett):,} TF-TF edges, "
    f"{len(nodes):,} node rows -- NO FIGURE BUILT")

# ======================================================================== Part 4: enrichment
hdr("PART 4 GO / PATHWAY ENRICHMENT")

def read_gmt(p):
    s = {}
    for line in open(p, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            g = clean(f[2:])
            if g:
                s[f[0]] = g
    return s

GMT = {k: read_gmt(os.path.join(REF, fn)) for k, fn in LIBS.items()}
for k, s in GMT.items():
    log(f"library {k}: {len(s)} sets from {LIBS[k]}")
R["_libraries"] = {k: {"file": LIBS[k], "n_sets": len(v)} for k, v in GMT.items()}

ENR, S4 = [], {}

def run_enrich(fg, bg, stage, lib, test, label, window):
    lo, hi = window
    rows = []
    fg, bg = fg & bg, set(bg)
    for term, gs in GMT[lib].items():
        t = gs & bg
        if not (lo <= len(t) <= hi):
            continue
        a = len(t & fg)
        tab = [[a, len(fg) - a], [len(t) - a, len(bg) - len(fg) - len(t) + a]]
        orr, p = stats.fisher_exact(tab, alternative="greater")
        rows.append({"stage": stage, "library": lib, "test": test, "label": label,
                     "term": term, "n_term_in_bg": len(t), "n_overlap": a,
                     "n_foreground": len(fg), "n_background": len(bg),
                     "odds_ratio": orr, "p": p, "window": f"{lo}-{hi}"})
    if rows:
        d = pd.DataFrame(rows)
        d["q"] = bh(d.p.values)
        return d
    return pd.DataFrame(columns=["stage", "library", "test", "label", "term", "q"])

for v in VOCABS:
    mod_u = set().union(*(TGT[v][r] for r in M55))
    pool_u = set().union(*(TGT[v][r] for r in POOL))
    U = UNIV[v]["all"]
    fg = mod_u & U
    S4[v] = {"module_union": len(mod_u), "module_union_in_universe": len(fg),
             "module_union_dropped_not_in_universe": len(mod_u - U),
             "pool_union": len(pool_u), "pool_union_in_universe": len(pool_u & U),
             "tests": {}}
    for test, bg in (("primary_vs_universe", U), ("secondary_vs_pool201", pool_u & U)):
        sat = len(fg & bg) / len(bg)
        tt = {"n_foreground": len(fg & bg), "n_background": len(bg),
              "saturation_fg_over_bg": float(sat), "by_library": {}}
        for lib in LIBS:
            d = run_enrich(fg, bg, v, lib, test, "module", (SIZE_MIN, SIZE_MAX))
            ds = run_enrich(fg, bg, v, lib, test, "module", SIZE_SENS)
            ENR.append(d)
            nsig = int((d.q < Q_BAR).sum()) if len(d) else 0
            tt["by_library"][lib] = {
                "n_terms_tested": int(len(d)), "n_sig_q05": nsig,
                "frac_sig": float(nsig / len(d)) if len(d) else float("nan"),
                "min_raw_p": float(d.p.min()) if len(d) else float("nan"),
                "min_q": float(d.q.min()) if len(d) else float("nan"),
                "n_sig_q05_sensitivity_5_1000": int((ds.q < Q_BAR).sum()) if len(ds) else 0,
                "top_terms": (d.nsmallest(TOP_N, "q")[["term", "n_overlap", "n_term_in_bg",
                                                       "odds_ratio", "p", "q"]]
                              .to_dict("records") if len(d) else [])}
            log(f"{v:9s} {test:22s} {lib:5s}: {len(d):5d} terms, {nsig:5d} sig q<0.05 "
                f"({100 * nsig / max(len(d), 1):.1f}%), min raw p "
                f"{d.p.min() if len(d) else float('nan'):.3g}")
        tt["SATURATED"] = bool(sat > 0.40)
        S4[v]["tests"][test] = tt
        log(f"{v:9s} {test:22s}: foreground {len(fg & bg):,} / background {len(bg):,} "
            f"= {sat:.3f}{'  [SATURATED, >0.40]' if sat > 0.40 else ''}")

# ---- conditional: per block, only if Part 1 found block structure
COND = {"per_block": {}, "per_hub_TF": {}}
for v in VOCABS:
    U = UNIV[v]["all"]
    bg = (set().union(*(TGT[v][r] for r in POOL))) & U
    if S1[v]["VERDICT"] == "BLOCK-STRUCTURED":
        for i, c in enumerate(S1[v]["modularity_greedy"]["communities"]):
            if len(c) < 3:
                continue
            f = set().union(*(TGT[v][r] for r in c)) & U
            for lib in LIBS:
                d = run_enrich(f, bg, v, lib, "block", f"block{i}", (SIZE_MIN, SIZE_MAX))
                ENR.append(d)
                COND["per_block"].setdefault(f"{v}|block{i}", {})[lib] = {
                    "n_regulons": len(c), "n_terms": int(len(d)),
                    "n_sig_q05": int((d.q < Q_BAR).sum()) if len(d) else 0,
                    "top_terms": (d.nsmallest(5, "q")[["term", "odds_ratio", "q"]]
                                  .to_dict("records") if len(d) else [])}
    else:
        COND["per_block"][v] = "NOT RUN -- Part 1 returned FLAT, and the per-block analysis is conditional on blocks"
    for tf in S3[v]["TF_degree"]["hubs"]:
        f = TF_TGT[v][tf] & U
        for lib in LIBS:
            d = run_enrich(f, bg, v, lib, "hub", tf, (SIZE_MIN, SIZE_MAX))
            ENR.append(d)
            COND["per_hub_TF"].setdefault(f"{v}|{tf}", {})[lib] = {
                "n_targets": len(f), "n_terms": int(len(d)),
                "n_sig_q05": int((d.q < Q_BAR).sum()) if len(d) else 0,
                "top_terms": (d.nsmallest(5, "q")[["term", "odds_ratio", "q"]]
                              .to_dict("records") if len(d) else [])}
R["S4_enrichment"] = S4
R["S4_conditional"] = COND
pd.concat([e for e in ENR if len(e)], ignore_index=True).to_csv(
    os.path.join(OUT, "f6_network_enrichment.csv"), index=False)

json.dump(R, open(os.path.join(OUT, "f6_network.json"), "w", encoding="utf-8"), indent=1,
          default=str)
hdr("DONE")
log(f"PART 2 VERDICT: {R['S2_VERDICT']['VERDICT']}")
log(f"PART 1 VERDICT: " + " | ".join(f"{v} {S1[v]['VERDICT']}" for v in VOCABS))
log(f"PART 3 VERDICT: " + " | ".join(f"{v} {S3[v]['bipartite']['VERDICT']}" for v in VOCABS))
log("wrote analysis/out/f6_network.json and six companions")
