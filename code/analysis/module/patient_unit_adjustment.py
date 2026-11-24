#!/usr/bin/env python
"""Unit-of-analysis correction: rerun the Part 3 degree and Part 4 per-hub enrichment at finer units (see f6_04)."""
import itertools
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "analysis", "out")
REF = os.path.join(ROOT, "analysis", "ref")
MODJ = os.path.join(ROOT, "modules", "out", "corrected_module.json")
NETJ = os.path.join(OUT, "f6_network.json")

VOCABS = ("Diagnosis", "Relapse")
SIZE_MIN, SIZE_MAX, Q_BAR = 10, 500, 0.05          # identical to the pre-specified run
LIBS = {"GOBP": "GO_Biological_Process_2025.gmt",
        "REAC": "Reactome_Pathways_2024.gmt",
        "HALL": "MSigDB_Hallmark_2020.gmt"}
HUBS_PREV = ["SPI1", "IRF8", "KLF4", "ZBTB7A", "MAFB"]   # union of the two Part 3 hub lists (f6_04)
ALL_MODE = ["IRF8", "KLF10", "SPI1", "ZBTB7A"]

R = {"_status": "POST-HOC CORRECTION of a unit-of-analysis error, requested 2026-08-31",
     "_what_changed": "only the unit of aggregation; question, background, libraries, window, "
                      "BH scheme and seed are unchanged",
     "_superseded": "per_TF results merge structurally disjoint target sets and MUST NOT BE QUOTED"}

def log(m):
    print(f"[unit] {m}", flush=True)

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

MOD = json.load(open(MODJ, encoding="utf-8"))
NET = json.load(open(NETJ, encoding="utf-8"))
M55, POOL = list(MOD["regulons"]), list(MOD["background_regulons"])
TGT = {v: {r: clean(MOD["targets"][v][r]) for r in POOL} for v in VOCABS}
MTF = sorted({parse_mode(r)[0] for r in M55})
UNIV = {}
for v in VOCABS:
    u = pd.read_csv(os.path.join(OUT, f"f1_universe_{v}.csv"))
    UNIV[v] = clean(u.symbol[u.symbol.notna()])
GMT = {k: read_gmt(os.path.join(REF, fn)) for k, fn in LIBS.items()}

# ---- the three unit definitions -------------------------------------------------------
def units_of(kind):
    """name -> list of regulons composing it."""
    d = defaultdict(list)
    for r in M55:
        tf, s1, s2 = parse_mode(r)
        key = {"per_TF": tf, "per_TF_s1": f"{tf}_{s1}", "per_regulon": r}[kind]
        d[key].append(r)
    return dict(d)

UNITS = {k: units_of(k) for k in ("per_TF", "per_TF_s1", "per_regulon")}
log(f"units: per_TF {len(UNITS['per_TF'])}, per_TF_s1 {len(UNITS['per_TF_s1'])}, "
    f"per_regulon {len(UNITS['per_regulon'])}")

# ======================================================================== degree
hdr("DEGREE AT ALL THREE UNITS")
DEG, DROWS = {}, []
for v in VOCABS:
    DEG[v] = {}
    for kind, u in UNITS.items():
        dd = {name: len(set().union(*(TGT[v][r] for r in regs))) for name, regs in u.items()}
        DEG[v][kind] = dd
        cut = float(np.percentile(list(dd.values()), 90))
        hubs = sorted([k for k, x in dd.items() if x >= cut], key=lambda k: -dd[k])
        s = pd.Series(dd)
        R.setdefault("degree", {}).setdefault(v, {})[kind] = {
            "n_units": len(dd), "min": int(s.min()), "median": float(s.median()),
            "max": int(s.max()), "hub_cut_p90": cut, "hubs": hubs,
            "per_unit": {k: int(x) for k, x in sorted(dd.items(), key=lambda t: -t[1])}}
        for name, regs in u.items():
            DROWS.append({"vocabulary": v, "unit_kind": kind, "unit": name,
                          "TF": parse_mode(regs[0])[0], "n_regulons": len(regs),
                          "degree": dd[name], "is_hub_p90": dd[name] >= cut})
        log(f"{v:9s} {kind:12s}: n={len(dd):3d} median {s.median():6.1f} max {s.max():5d} "
            f"p90 cut {cut:6.1f} | top5 {[f'{k}:{dd[k]}' for k in hubs[:5]]}")
pd.DataFrame(DROWS).to_csv(os.path.join(OUT, "f6_unit_correction_degree.csv"), index=False)

# ---- the exact decomposition: cross-s1 sets are disjoint, so degrees must sum -----------
hdr("DECOMPOSITION CHECK -- per-TF degree must equal the SUM of its s1 components")
DEC = {}
for v in VOCABS:
    DEC[v] = {}
    bad = 0
    for tf in MTF:
        comps = {k: DEG[v]["per_TF_s1"][k] for k in DEG[v]["per_TF_s1"] if k.rsplit("_", 1)[0] == tf}
        tot, ssum = DEG[v]["per_TF"][tf], sum(comps.values())
        ok = (tot == ssum)
        bad += (not ok)
        DEC[v][tf] = {"per_TF_degree": tot, "s1_components": comps,
                      "sum_of_components": ssum, "sums_exactly": bool(ok),
                      "largest_component_share": max(comps.values()) / tot if tot else float("nan"),
                      "n_components": len(comps)}
    R.setdefault("decomposition", {})[v] = DEC[v]
    log(f"{v}: {len(MTF)} TFs, exact-sum failures {bad}  "
        f"{'PASS -- cross-s1 target sets are disjoint, as SCENICPLUS_SIGN_MODE_SEMANTICS.md predicts' if bad == 0 else 'FAIL'}")
    for tf in HUBS_PREV:
        d = DEC[v][tf]
        log(f"    {tf:7s} {d['per_TF_degree']:5d} = " +
            " + ".join(f"{k.rsplit('_',1)[1]}:{x}" for k, x in sorted(d['s1_components'].items())) +
            f"   (largest component {100*d['largest_component_share']:.0f}%)")

# ======================================================================== enrichment
hdr("ENRICHMENT AT ALL THREE UNITS (background = 201-pool union n universe, window 10-500)")
EROWS = []

def run_enrich(fg, bg, lib, stage, kind, unit, n_regs):
    fg2 = fg & bg
    rows = []
    for term, gs in GMT[lib].items():
        t = gs & bg
        if not (SIZE_MIN <= len(t) <= SIZE_MAX):
            continue
        a = len(t & fg2)
        orr, p = stats.fisher_exact(
            [[a, len(fg2) - a], [len(t) - a, len(bg) - len(fg2) - len(t) + a]],
            alternative="greater")
        rows.append({"stage": stage, "unit_kind": kind, "unit": unit, "n_regulons": n_regs,
                     "library": lib, "term": term, "n_term_in_bg": len(t), "n_overlap": a,
                     "n_foreground": len(fg2), "n_background": len(bg),
                     "odds_ratio": orr, "p": p})
    d = pd.DataFrame(rows)
    if len(d):
        d["q"] = bh(d.p.values)          # BH within (library x stage x unit), as pre-specified
    return d

ENR = {}
for v in VOCABS:
    bg = (set().union(*(TGT[v][r] for r in POOL))) & UNIV[v]
    for kind, u in UNITS.items():
        for unit, regs in u.items():
            fg = set().union(*(TGT[v][r] for r in regs)) & UNIV[v]
            for lib in LIBS:
                d = run_enrich(fg, bg, lib, v, kind, unit, len(regs))
                if len(d):
                    EROWS.append(d)
                    ENR[(v, kind, unit, lib)] = d
    log(f"{v}: background {len(bg)} genes; "
        f"{sum(len(u) for u in UNITS.values())} units x 3 libraries done")
ALLE = pd.concat(EROWS, ignore_index=True)
ALLE.to_csv(os.path.join(OUT, "f6_unit_correction_enrichment.csv"), index=False)

def top(v, kind, unit, lib, n=5):
    d = ENR.get((v, kind, unit, lib))
    if d is None or not len(d):
        return []
    return d.nsmallest(n, "q")[["term", "n_overlap", "n_term_in_bg", "odds_ratio", "q"]] \
            .to_dict("records")

def nsig(v, kind, unit, lib):
    d = ENR.get((v, kind, unit, lib))
    return int((d.q < Q_BAR).sum()) if d is not None and len(d) else 0

# ---- the five previous hubs, decomposed -----------------------------------------------
hdr("THE FIVE PREVIOUS HUBS -- does each signature survive, split, or disappear?")
HUB = {}
for v in VOCABS:
    for tf in HUBS_PREV:
        comps = sorted(k for k in UNITS["per_TF_s1"] if k.rsplit("_", 1)[0] == tf)
        rec = {"per_TF_SUPERSEDED": {lib: {"n_sig": nsig(v, "per_TF", tf, lib),
                                           "top": top(v, "per_TF", tf, lib, 3)} for lib in LIBS},
               "per_TF_s1": {c: {"n_regulons": len(UNITS["per_TF_s1"][c]),
                                 "degree": DEG[v]["per_TF_s1"][c],
                                 **{lib: {"n_sig": nsig(v, "per_TF_s1", c, lib),
                                          "top": top(v, "per_TF_s1", c, lib, 3)} for lib in LIBS}}
                             for c in comps},
               "per_regulon": {r: {"degree": DEG[v]["per_regulon"][r],
                                   **{lib: {"n_sig": nsig(v, "per_regulon", r, lib),
                                            "top": top(v, "per_regulon", r, lib, 2)}
                                      for lib in LIBS}}
                               for r in sorted(UNITS["per_TF"][tf])}}
        HUB[f"{v}|{tf}"] = rec
        log(f"\n{v} / {tf}: per-TF HALL {nsig(v,'per_TF',tf,'HALL')} sig, "
            f"REAC {nsig(v,'per_TF',tf,'REAC')} sig  [SUPERSEDED]")
        for c in comps:
            t3 = "; ".join(f"{x['term'][:38]} q{x['q']:.1g}" for x in top(v, "per_TF_s1", c, "REAC", 2))
            log(f"    {c:12s} deg {DEG[v]['per_TF_s1'][c]:5d} n={len(UNITS['per_TF_s1'][c])} | "
                f"REAC {nsig(v,'per_TF_s1',c,'REAC'):3d} sig | {t3}")
R["hub_decomposition"] = HUB

# ---- per-regulon enrichment power, reported as a limitation ----------------------------
hdr("PER-REGULON POWER -- how many of the 55 can enrich at all")
POW = {}
for v in VOCABS:
    sizes = {r: len(TGT[v][r] & UNIV[v]) for r in M55}
    any_sig = {r: any(nsig(v, "per_regulon", r, lib) > 0 for lib in LIBS) for r in M55}
    POW[v] = {"n_regulons": 55,
              "foreground_min": int(min(sizes.values())),
              "foreground_median": float(np.median(list(sizes.values()))),
              "foreground_max": int(max(sizes.values())),
              "n_with_any_sig_term": int(sum(any_sig.values())),
              "n_with_no_sig_term": int(55 - sum(any_sig.values())),
              "regulons_with_no_sig_term": sorted(r for r in M55 if not any_sig[r]),
              "median_foreground_of_silent": float(np.median(
                  [sizes[r] for r in M55 if not any_sig[r]]) if any(not x for x in any_sig.values())
                  else float("nan")),
              "median_foreground_of_enriching": float(np.median(
                  [sizes[r] for r in M55 if any_sig[r]]))}
    log(f"{v}: {POW[v]['n_with_any_sig_term']}/55 regulons yield >=1 term at q<0.05; "
        f"silent median foreground {POW[v]['median_foreground_of_silent']:.0f} vs enriching "
        f"{POW[v]['median_foreground_of_enriching']:.0f}")
R["per_regulon_power"] = POW

# ---- confirm the finer units do not touch Part 1 blocks or module-level enrichment ----------
hdr("CONFIRMATION -- Part 1 blocks and module-level enrichment were already per-regulon")
CONF = {}
for v in VOCABS:
    mu = set().union(*(TGT[v][r] for r in M55)) & UNIV[v]
    bg = (set().union(*(TGT[v][r] for r in POOL))) & UNIV[v]
    prev = NET["S4_enrichment"][v]["tests"]["secondary_vs_pool201"]
    now = {}
    for lib in LIBS:
        d = run_enrich(mu, bg, lib, v, "module", "module_union", 55)
        now[lib] = {"n_terms": int(len(d)), "n_sig": int((d.q < Q_BAR).sum())}
    CONF[v] = {
        "module_union_enrichment_unchanged": {
            lib: {"now_n_sig": now[lib]["n_sig"],
                  "committed_n_sig": prev["by_library"][lib]["n_sig_q05"],
                  "match": bool(now[lib]["n_sig"] == prev["by_library"][lib]["n_sig_q05"])}
            for lib in LIBS},
        "S1_blocks_unit": "per regulon in the committed run; the Jaccard graph's nodes are the 55 "
                          "regulons, so no aggregation occurred and no block can change",
        "S1_communities_committed": [len(c) for c in
                                     NET["S1_target_overlap"][v]["modularity_greedy"]["communities"]]}
    ok = all(x["match"] for x in CONF[v]["module_union_enrichment_unchanged"].values())
    log(f"{v}: module-level enrichment reproduces committed counts exactly: {ok} "
        f"({ {lib: (now[lib]['n_sig'], prev['by_library'][lib]['n_sig_q05']) for lib in LIBS} })")
R["confirmation"] = CONF

json.dump(R, open(os.path.join(OUT, "f6_unit_correction.json"), "w", encoding="utf-8"),
          indent=1, default=str)
hdr("DONE")
log("wrote analysis/out/f6_unit_correction.json + degree.csv + enrichment.csv")
