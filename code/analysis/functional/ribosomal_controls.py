#!/usr/bin/env python
"""Is the module's translation block a ribosomal-gene detection artifact?"""
import itertools
import json
import os
import re
from collections import Counter

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "analysis", "out")
REF = os.path.join(ROOT, "analysis", "ref")
GTW = os.path.join(ROOT, "analysis/out/gene_detection")
MODJ = os.path.join(ROOT, "modules", "out", "corrected_module.json")
NETJ = os.path.join(OUT, "f6_network.json")

VOCABS = ("Diagnosis", "Relapse")
SEED, N_DRAWS = 555, 1000                     # same as the network analysis (f6_04)
SIZE_MIN, SIZE_MAX, Q_BAR = 10, 500, 0.05     # same as the network analysis (f6_04)
LIBS = {"GOBP": "GO_Biological_Process_2025.gmt",
        "REAC": "Reactome_Pathways_2024.gmt",
        "HALL": "MSigDB_Hallmark_2020.gmt"}

R = {"_status": "POST-HOC ROBUSTNESS CHECK, requested after seeing the enrichment",
     "_prespec_rule_reused": "the pre-specified Part 1 rule of f6_04, unchanged",
     "_block_numbering": "1-indexed here; f6_network.json is 0-indexed. block N here = "
                         "communities[N-1] there",
     "_seed": SEED, "_n_draws": N_DRAWS}

def log(m):
    print(f"[ribo] {m}", flush=True)

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

def parse_mode(r):
    tf, m = r.split("_direct_")
    return tf, m[0], m[2]

def read_gmt(p):
    s = {}
    for line in open(p, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            g = clean(f[2:])
            if g:
                s[f[0]] = g
    return s

# ======================================================================== substrate
MOD = json.load(open(MODJ, encoding="utf-8"))
NET = json.load(open(NETJ, encoding="utf-8"))
M55, POOL = list(MOD["regulons"]), list(MOD["background_regulons"])
TGT0 = {v: {r: clean(MOD["targets"][v][r]) for r in POOL} for v in VOCABS}
MTF = sorted({parse_mode(r)[0] for r in M55})
UNIV = {}
for v in VOCABS:
    u = pd.read_csv(os.path.join(OUT, f"f1_universe_{v}.csv"))
    UNIV[v] = clean(u.symbol[u.symbol.notna()])
GMT = {k: read_gmt(os.path.join(REF, fn)) for k, fn in LIBS.items()}
log(f"module {len(M55)} / pool {len(POOL)}; universe "
    f"{ {v: len(UNIV[v]) for v in VOCABS} }")

COMM0 = {v: [sorted(c) for c in NET["S1_target_overlap"][v]["modularity_greedy"]["communities"]]
         for v in VOCABS}
BLK0 = {v: {r: i + 1 for i, c in enumerate(COMM0[v]) for r in c} for v in VOCABS}
R["_block_map_to_f6_network_json"] = {
    v: {f"block {i+1} (here)": f"communities[{i}] (there), n={len(c)}"
        for i, c in enumerate(COMM0[v])} for v in VOCABS}

# ======================================================================== the gene sets
hdr("RIBOSOMAL GENE DEFINITIONS")
ALLSYM = sorted(UNIV["Diagnosis"] | UNIV["Relapse"])
CYTO = {g for g in ALLSYM if re.match(r"^RP[LS]", g)}
KIN = {g for g in CYTO if re.match(r"^RPS6K", g)}          # S6 kinases, NOT ribosomal proteins
NOTRP = {"RPS19BP1"}                                       # RPS19 *binding protein*, not a RP
CYTO = CYTO - KIN - NOTRP
MITO = {g for g in ALLSYM if re.match(r"^MRP[LS]", g)}
DROP_A = CYTO | MITO
# stricter sensitivity: the three ribosomal proteins whose symbols do not start RP[LS],
# plus every gene in the two translation terms that actually drove the enrichment
EXTRA = {g for g in ("UBA52", "FAU", "RACK1") if g in set(ALLSYM)}
TERMS = {}
for lib, term in (("REAC", "Peptide Chain Elongation"),
                  ("GOBP", "Cytoplasmic Translation (GO:0002181)")):
    hit = [t for t in GMT[lib] if t == term]
    if hit:
        TERMS[f"{lib}:{term}"] = GMT[lib][hit[0]] & set(ALLSYM)
DROP_B = DROP_A | EXTRA | set().union(*TERMS.values()) if TERMS else DROP_A | EXTRA
log(f"DROP-A (primary): ^RP[LS] minus ^RPS6K* ({len(KIN)}) minus RPS19BP1, plus ^MRP[LS] "
    f"-> cyto {len(CYTO)} + mito {len(MITO)} = {len(DROP_A)} genes")
log(f"DROP-B (sensitivity): DROP-A + {sorted(EXTRA)} + "
    f"{ {k: len(v) for k, v in TERMS.items()} } -> {len(DROP_B)} genes")
R["gene_sets"] = {
    "DROP_A": {"rule": "symbol =~ ^RP[LS] excluding ^RPS6K* and RPS19BP1, union symbol =~ ^MRP[LS]",
               "source": "symbols of analysis/out/f1_universe_{vocab}.csv",
               "n_cytoplasmic": len(CYTO), "n_mitochondrial": len(MITO), "n_total": len(DROP_A),
               "excluded_RPS6K": sorted(KIN), "excluded_other": sorted(NOTRP),
               "genes": sorted(DROP_A)},
    "DROP_B": {"rule": "DROP_A + non-RP-named ribosomal proteins + the two enriched translation "
                       "term gene sets (the maximally aggressive drop)",
               "extra_symbols": sorted(EXTRA),
               "terms_included": {k: len(v) for k, v in TERMS.items()},
               "n_total": len(DROP_B), "genes": sorted(DROP_B)}}
pd.DataFrame([{"symbol": g, "in_DROP_A": g in DROP_A, "in_DROP_B": g in DROP_B,
               "class": ("cytoplasmic_RP" if g in CYTO else
                         "mitochondrial_RP" if g in MITO else "other")}
              for g in sorted(DROP_B)]).to_csv(
    os.path.join(OUT, "f6_ribosome_gene_list.csv"), index=False)

# ======================================================================== Part 1 composition
hdr("PART 1 COMPOSITION -- how much of each block, and of its EDGES, is ribosomal")
COMP = {}
for v in VOCABS:
    COMP[v] = {}
    for bi, c in enumerate(COMM0[v], start=1):
        un = set().union(*(TGT0[v][r] for r in c))
        # the overlaps that define the edges: pooled pairwise intersections within the block
        inter = Counter()
        for a, b in itertools.combinations(c, 2):
            inter.update(TGT0[v][a] & TGT0[v][b])
        ninter = sum(inter.values())
        nribo = sum(n for g, n in inter.items() if g in DROP_A)
        COMP[v][f"block{bi}"] = {
            "n_regulons": len(c), "n_target_union": len(un),
            "n_ribosomal_in_union": len(un & DROP_A),
            "frac_ribosomal_of_union": len(un & DROP_A) / len(un) if un else float("nan"),
            "n_edge_overlap_gene_instances": ninter,
            "n_ribosomal_edge_overlap_instances": nribo,
            "frac_ribosomal_of_edge_overlaps": nribo / ninter if ninter else float("nan"),
            "n_distinct_overlap_genes": len(inter),
            "frac_ribosomal_of_distinct_overlap_genes":
                (len([g for g in inter if g in DROP_A]) / len(inter)) if inter else float("nan")}
        b = COMP[v][f"block{bi}"]
        log(f"{v:9s} block{bi} (n={len(c):2d}): union {len(un):5d} genes, "
            f"{b['frac_ribosomal_of_union']:.4f} ribosomal | edge overlaps "
            f"{ninter:6d} instances, {b['frac_ribosomal_of_edge_overlaps']:.4f} ribosomal "
            f"({b['frac_ribosomal_of_distinct_overlap_genes']:.4f} of distinct)")
    mu = set().union(*(TGT0[v][r] for r in M55))
    COMP[v]["module_union"] = {"n": len(mu), "n_ribosomal": len(mu & DROP_A),
                               "frac_ribosomal": len(mu & DROP_A) / len(mu)}
R["S1_composition"] = COMP

def jac_graph(regs, tgt):
    g = nx.Graph()
    g.add_nodes_from(regs)
    for a, b in itertools.combinations(regs, 2):
        w = jac(tgt[a], tgt[b])
        if w > 0:
            g.add_edge(a, b, weight=w)
    return g

def modularity_of(regs, tgt):
    g = jac_graph(regs, tgt)
    if g.number_of_edges() == 0:
        return 0.0, [{r} for r in regs]
    comm = [set(c) for c in nx.community.greedy_modularity_communities(g, weight="weight")]
    iso = set(regs) - set().union(*comm) if comm else set(regs)
    comm += [{r} for r in iso]
    return nx.community.modularity(g, comm, weight="weight"), comm

# ---- PARITY CHECK ------------------------------------------------------------------
hdr("PARITY CHECK -- reproduce the committed Q before touching anything")
for v in VOCABS:
    q, _ = modularity_of(M55, TGT0[v])
    qc = NET["S1_target_overlap"][v]["modularity_greedy"]["Q_observed"]
    if abs(q - qc) > 1e-12:
        raise SystemExit(f"PARITY FAIL {v}: recomputed Q {q!r} != committed {qc!r}")
    log(f"{v}: Q {q:.6f} == committed {qc:.6f}  PASS")
R["_parity_gate"] = {"rule": "recomputed undepleted Q must equal f6_network.json to 1e-12",
                     "PASS": True,
                     "Q": {v: NET["S1_target_overlap"][v]["modularity_greedy"]["Q_observed"]
                           for v in VOCABS}}

# ======================================================================== Part 2 drop test
hdr("PART 2 THE DROP TEST -- same rule, same null, ribosomal genes removed")
rng = np.random.default_rng(SEED)
DRAWS = [list(rng.choice(POOL, size=55, replace=False)) for _ in range(N_DRAWS)]
DROPS = {"DROP_A": DROP_A, "DROP_B": DROP_B}
S2 = {}
for dname, dset in DROPS.items():
    S2[dname] = {}
    for v in VOCABS:
        tgt = {r: (TGT0[v][r] - dset) for r in POOL}
        empt = [r for r in M55 if not tgt[r]]
        q, comm = modularity_of(M55, tgt)
        nullq = np.array([modularity_of(d, tgt)[0] for d in DRAWS])
        p = (1 + int((nullq >= q).sum())) / (1 + N_DRAWS)
        comm = sorted(comm, key=len, reverse=True)
        big = sum(1 for c in comm if len(c) >= 3)
        verdict = "BLOCK-STRUCTURED" if (p < 0.05 and big >= 2) else "FLAT"
        blk = {r: i + 1 for i, c in enumerate(comm) for r in c}
        # membership change: best-matching new block for each old block, then count movers
        match = {}
        for i, c0 in enumerate(COMM0[v], start=1):
            j = max(range(len(comm)), key=lambda k: len(set(c0) & comm[k]))
            match[i] = j + 1
        moved = [r for r in M55 if match[BLK0[v][r]] != blk[r]]
        S2[dname][v] = {
            "n_targets_removed_module_union": len(
                set().union(*(TGT0[v][r] for r in M55)) & dset),
            "regulons_left_with_no_targets": empt,
            "Q_observed": float(q),
            "Q_committed_undepleted": NET["S1_target_overlap"][v]["modularity_greedy"]["Q_observed"],
            "null_mean": float(nullq.mean()), "null_sd": float(nullq.std(ddof=1)),
            "null_p95": float(np.percentile(nullq, 95)),
            "p_empirical_one_sided": float(p),
            "percentile_of_observed": float((nullq < q).mean() * 100),
            "n_communities": len(comm), "n_communities_ge3": big,
            "community_sizes": [len(c) for c in comm],
            "communities": [sorted(c) for c in comm],
            "old_to_new_block_match": match,
            "n_regulons_changing_block": len(moved), "regulons_changing_block": sorted(moved),
            "VERDICT": verdict}
        log(f"{dname} {v:9s}: removed {S2[dname][v]['n_targets_removed_module_union']:4d} genes "
            f"from the union | Q {q:.4f} (was {S2[dname][v]['Q_committed_undepleted']:.4f}) vs "
            f"null {nullq.mean():.4f}+-{nullq.std(ddof=1):.4f}, p {p:.4f}, "
            f"{len(comm)} communities {[len(c) for c in comm]}, "
            f"{len(moved)}/55 moved -> {verdict}")
R["S2_drop_test"] = S2
pd.DataFrame([{"vocabulary": v, "regulon": r, "TF": parse_mode(r)[0],
               "block_before": BLK0[v][r],
               "block_after_DROP_A": {rr: i + 1 for i, c in
                                      enumerate(S2["DROP_A"][v]["communities"]) for rr in c}[r],
               "block_after_DROP_B": {rr: i + 1 for i, c in
                                      enumerate(S2["DROP_B"][v]["communities"]) for rr in c}[r]}
              for v in VOCABS for r in M55]).to_csv(
    os.path.join(OUT, "f6_ribosome_block_membership.csv"), index=False)

# ======================================================================== Part 2b enrichment
hdr("PART 2b PER-BLOCK ENRICHMENT ON THE DEPLETED SETS (background = depleted 201-pool)")
ENR = []

def run_enrich(fg, bg, lib, label, stage, tag):
    fg, bg = fg & bg, set(bg)
    rows = []
    for term, gs in GMT[lib].items():
        t = gs & bg
        if not (SIZE_MIN <= len(t) <= SIZE_MAX):
            continue
        a = len(t & fg)
        orr, p = stats.fisher_exact(
            [[a, len(fg) - a], [len(t) - a, len(bg) - len(fg) - len(t) + a]],
            alternative="greater")
        rows.append({"stage": stage, "drop_set": tag, "library": lib, "block": label,
                     "term": term, "n_term_in_bg": len(t), "n_overlap": a,
                     "n_foreground": len(fg), "n_background": len(bg),
                     "odds_ratio": orr, "p": p})
    d = pd.DataFrame(rows)
    if len(d):
        d["q"] = bh(d.p.values)
    return d

S2B = {}
for dname, dset in DROPS.items():
    for v in VOCABS:
        tgt = {r: (TGT0[v][r] - dset) for r in POOL}
        bg = (set().union(*(tgt[r] for r in POOL))) & UNIV[v]
        for bi, c in enumerate(S2[dname][v]["communities"], start=1):
            if len(c) < 3:
                continue
            fg = set().union(*(tgt[r] for r in c)) & UNIV[v]
            for lib in LIBS:
                d = run_enrich(fg, bg, lib, f"block{bi}", v, dname)
                ENR.append(d)
                k = f"{dname}|{v}|block{bi}|{lib}"
                S2B[k] = {"n_regulons": len(c), "n_foreground": len(fg), "n_background": len(bg),
                          "n_terms": int(len(d)),
                          "n_sig_q05": int((d.q < Q_BAR).sum()) if len(d) else 0,
                          "top_terms": (d.nsmallest(8, "q")[["term", "n_overlap", "n_term_in_bg",
                                                             "odds_ratio", "q"]].to_dict("records")
                                        if len(d) else [])}
                if dname == "DROP_A":
                    top = "; ".join(f"{t['term'][:44]} q{t['q']:.2g}"
                                    for t in S2B[k]["top_terms"][:3])
                    log(f"{v:9s} block{bi} {lib:5s}: {S2B[k]['n_sig_q05']:4d} sig | {top}")
R["S2b_enrichment_after_drop"] = S2B
pd.concat([e for e in ENR if len(e)], ignore_index=True).to_csv(
    os.path.join(OUT, "f6_ribosome_enrichment.csv"), index=False)

# ---- Myc Targets V1 specifically, before and after (both drops) ------------------
hdr("MYC TARGETS V1 -- does it survive the drop independently of the translation terms?")
MYC = {}
allenr = pd.concat([e for e in ENR if len(e)], ignore_index=True)
for v in VOCABS:
    for dname in DROPS:
        s = allenr[(allenr.stage == v) & (allenr.drop_set == dname) &
                   (allenr.library == "HALL") & (allenr.term == "Myc Targets V1")]
        for _, row in s.iterrows():
            MYC[f"{dname}|{v}|{row.block}"] = {
                "q": float(row.q), "odds_ratio": float(row.odds_ratio),
                "n_overlap": int(row.n_overlap), "n_term_in_bg": int(row.n_term_in_bg),
                "significant_q05": bool(row.q < Q_BAR)}
    myc_genes = GMT["HALL"].get("Myc Targets V1", set())
    MYC[f"_overlap_with_DROP_A|{v}"] = {
        "n_myc_genes": len(myc_genes),
        "n_myc_genes_ribosomal": len(myc_genes & DROP_A),
        "frac_myc_genes_ribosomal": len(myc_genes & DROP_A) / len(myc_genes) if myc_genes else 0.0}
for k, x in sorted(MYC.items()):
    if k.startswith("_"):
        log(f"{k}: {x}")
    elif k.startswith("DROP_A"):
        log(f"{k}: q {x['q']:.3g} OR {x['odds_ratio']:.2f} "
            f"{x['n_overlap']}/{x['n_term_in_bg']} sig={x['significant_q05']}")
R["MYC_targets_v1"] = MYC

# ======================================================================== Part 3 detection
hdr("PART 3 DETECTION CONFOUND -- are block 1's targets systematically higher expressed?")
S3 = {}
for v in VOCABS:
    det = pd.read_csv(os.path.join(GTW, f"gene_detection_malignant_{v}.csv"))
    det = det[det.symbol.notna()].drop_duplicates("symbol").set_index("symbol")
    covered = set(det.index)
    grp = {}
    for bi, c in enumerate(COMM0[v], start=1):
        grp[f"block{bi}"] = (set().union(*(TGT0[v][r] for r in c))) & covered
    grp["module_union"] = (set().union(*(TGT0[v][r] for r in M55))) & covered
    grp["universe"] = UNIV[v] & covered
    grp["ribosomal_DROP_A"] = DROP_A & covered
    st = {}
    for k, gs in grp.items():
        s = det.loc[sorted(gs)]
        st[k] = {"n_genes": len(gs),
                 "mean_detection_rate": float(s.detection_rate.mean()),
                 "median_detection_rate": float(s.detection_rate.median()),
                 "mean_mean_count": float(s.mean_count.mean()),
                 "median_mean_count": float(s.mean_count.median())}
    b1 = det.loc[sorted(grp["block1"])]
    b2 = det.loc[sorted(grp["block2"])]
    uni = det.loc[sorted(grp["universe"])]
    tests = {}
    for nm, a, b in (("block1_vs_block2", b1, b2), ("block1_vs_universe", b1, uni),
                     ("block2_vs_universe", b2, uni)):
        u1, p1 = stats.mannwhitneyu(a.detection_rate, b.detection_rate, alternative="two-sided")
        u2, p2 = stats.mannwhitneyu(a.mean_count, b.mean_count, alternative="two-sided")
        tests[nm] = {"detection_rate_MW_p": float(p1),
                     "detection_rate_median_a": float(a.detection_rate.median()),
                     "detection_rate_median_b": float(b.detection_rate.median()),
                     "mean_count_MW_p": float(p2),
                     "mean_count_median_a": float(a.mean_count.median()),
                     "mean_count_median_b": float(b.mean_count.median())}
    # the same comparison with ribosomal genes removed -- is the shift only the RP genes?
    b1n = det.loc[sorted(grp["block1"] - DROP_A)]
    b2n = det.loc[sorted(grp["block2"] - DROP_A)]
    un, pn = stats.mannwhitneyu(b1n.detection_rate, b2n.detection_rate, alternative="two-sided")
    tests["block1_vs_block2_ribosome_removed"] = {
        "n_a": len(b1n), "n_b": len(b2n),
        "detection_rate_median_a": float(b1n.detection_rate.median()),
        "detection_rate_median_b": float(b2n.detection_rate.median()),
        "detection_rate_MW_p": float(pn),
        "mean_count_median_a": float(b1n.mean_count.median()),
        "mean_count_median_b": float(b2n.mean_count.median())}
    S3[v] = {"_source": f"analysis/out/gene_detection/gene_detection_malignant_{v}.csv",
             "_panel_note": "15,615-gene externally fixed panel (Arm B); coverage of the "
                            "universe is stated per group",
             "n_universe_covered": len(grp["universe"]), "n_universe_total": len(UNIV[v]),
             "groups": st, "tests": tests}
    log(f"{v}: median detection block1 {st['block1']['median_detection_rate']:.4f} | "
        f"block2 {st['block2']['median_detection_rate']:.4f} | universe "
        f"{st['universe']['median_detection_rate']:.4f} | ribosomal "
        f"{st['ribosomal_DROP_A']['median_detection_rate']:.4f}")
    log(f"    b1 vs b2 MW p {tests['block1_vs_block2']['detection_rate_MW_p']:.3g}; "
        f"with ribosomal genes removed p {pn:.3g} "
        f"(medians {b1n.detection_rate.median():.4f} vs {b2n.detection_rate.median():.4f})")
R["S3_detection"] = S3

# ======================================================================== Part 4 generalisation
hdr("PART 4 DOES THE CONCERN GENERALISE -- ribosomal share, module union vs the 201-pool union")
S4 = {}
for v in VOCABS:
    mu = set().union(*(TGT0[v][r] for r in M55))
    pu = set().union(*(TGT0[v][r] for r in POOL))
    a, b = len(mu & DROP_A), len(mu - DROP_A)
    c, d = len(pu & DROP_A) - a, len(pu - DROP_A) - b
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    S4[v] = {"module_union": len(mu), "module_ribosomal": a,
             "module_frac_ribosomal": a / len(mu),
             "pool_union": len(pu), "pool_ribosomal": len(pu & DROP_A),
             "pool_frac_ribosomal": len(pu & DROP_A) / len(pu),
             "pool_minus_module_ribosomal": c, "pool_minus_module_other": d,
             "fisher_odds_ratio": float(orr), "fisher_p_greater": float(p),
             "ENRICHED_IN_MODULE": bool(p < 0.05)}
    log(f"{v}: module {a}/{len(mu)} = {a/len(mu):.4f} ribosomal | pool "
        f"{len(pu & DROP_A)}/{len(pu)} = {len(pu & DROP_A)/len(pu):.4f} | OR {orr:.3f} "
        f"p {p:.4g} -> {'ENRICHED' if p < 0.05 else 'NOT enriched'}")
R["S4_generalisation"] = S4

json.dump(R, open(os.path.join(OUT, "f6_ribosome_check.json"), "w", encoding="utf-8"),
          indent=1, default=str)
hdr("DONE")
for dname in DROPS:
    log(f"{dname}: " + " | ".join(f"{v} {S2[dname][v]['VERDICT']} (Q {S2[dname][v]['Q_observed']:.4f}, "
                                  f"p {S2[dname][v]['p_empirical_one_sided']:.4f})" for v in VOCABS))
log("wrote analysis/out/f6_ribosome_check.json and three companions")
