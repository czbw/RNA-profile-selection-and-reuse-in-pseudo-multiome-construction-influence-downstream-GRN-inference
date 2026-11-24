#!/usr/bin/env python
"""The reproducibility ceiling: replicate 1 vs replicate 2 of the SAME arm."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

TOPN = 100_000
regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[ceiling] {m}", flush=True)

def die(m):
    raise SystemExit(f"[ceiling] ABORT — {m}")

def jac(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else float("nan")

def need(p: Path) -> Path:
    if not p.is_file() or p.stat().st_size == 0:
        die(f"missing/empty {p}")
    return p

def top_edges(path: Path, cols: list[str], score: str, n: int = TOPN) -> set:
    """Top-n rows by `score`, as a set of tuples over `cols`. Ties at the boundary are INCLUDED
    on both sides rather than cut arbitrarily, so the comparison cannot be made to look better or
    worse by where a tie happens to fall."""
    df = pd.read_csv(path, sep="\t")
    for c in cols + [score]:
        if c not in df.columns:
            die(f"{path.name} has no column {c!r}; it has {list(df.columns)[:12]}")
    df = df.sort_values(score, ascending=False)
    if len(df) > n:
        cut = df[score].iloc[n - 1]
        df = df[df[score] >= cut]
    return set(map(tuple, df[cols].astype(str).values))

def regulon_targets(path: Path) -> dict[str, set]:
    df = pd.read_csv(path, sep="\t", usecols=["eRegulon_name", "Gene"])
    df["key"] = [regkey(v) for v in df.eRegulon_name]
    return {k: set(map(str, g.unique())) for k, g in df.groupby("key").Gene}

def tfs_of(path: Path) -> set:
    df = pd.read_csv(path, sep="\t", usecols=["TF"])
    return set(map(str, df.TF.unique()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="benchmarks/pbmc/out")
    ap.add_argument("--rep", required=True, help="this bundle's out/")
    ap.add_argument("--stage", default="pbmc10k")
    ap.add_argument("--arm", default="truth")
    ap.add_argument("--reptag", default="rep2")
    a = ap.parse_args()

    S, R = Path(a.src), Path(a.rep)
    t1, t2 = f"{a.stage}_{a.arm}", f"{a.stage}_{a.arm}_{a.reptag}"
    log(f"replicate 1 = {t1}  (in {S})")
    log(f"replicate 2 = {t2}  (in {R})")

    out = {"replicate_1": t1, "replicate_2": t2,
           "what_this_is": "run-to-run reproducibility ceiling of ONE arm on identical inputs; "
                           "scenicplus 1.0a2 exposes no seed, so this is a RE-RUN replicate",
           "top_n_edges": TOPN, "levels": {}}

    # ---- L1 region-gene ----------------------------------------------------------------------
    p1, p2 = need(S / f"r2g_{t1}.tsv"), need(R / f"r2g_{t2}.tsv")
    A = top_edges(p1, ["region", "target"], "importance")
    B = top_edges(p2, ["region", "target"], "importance")
    out["levels"]["L1_region_gene_links"] = {"n_1": len(A), "n_2": len(B), "jaccard": jac(A, B)}
    log(f"L1 region-gene links   : {jac(A,B):.4f}   ({len(A):,} vs {len(B):,})")

    # ---- L2 TF-gene --------------------------------------------------------------------------
    p1, p2 = need(S / f"tf2g_{t1}.tsv"), need(R / f"tf2g_{t2}.tsv")
    A = top_edges(p1, ["TF", "target"], "importance")
    B = top_edges(p2, ["TF", "target"], "importance")
    out["levels"]["L2_tf_gene_edges"] = {"n_1": len(A), "n_2": len(B), "jaccard": jac(A, B)}
    log(f"L2 TF-gene edges       : {jac(A,B):.4f}   ({len(A):,} vs {len(B):,})")

    # ---- L3/L4/L5 from the eRegulon tables ------------------------------------------------------
    e1, e2 = need(S / f"eRegulons_{t1}.tsv"), need(R / f"eRegulons_{t2}.tsv")
    T1, T2 = regulon_targets(e1), regulon_targets(e2)
    shared = sorted(set(T1) & set(T2))
    per = [jac(T1[k], T2[k]) for k in shared]
    med = float(np.median(per)) if per else float("nan")
    out["levels"]["L3_target_sets_median"] = {
        "n_shared_regulons": len(shared), "median_jaccard": med,
        "q25": float(np.percentile(per, 25)) if per else None,
        "q75": float(np.percentile(per, 75)) if per else None}
    log(f"L3 target sets (median): {med:.4f}   over {len(shared):,} shared-name regulons")

    j4 = jac(set(T1), set(T2))
    out["levels"]["L4_eregulons"] = {"n_1": len(T1), "n_2": len(T2), "jaccard": j4}
    log(f"L4 eRegulons           : {j4:.4f}   ({len(T1):,} vs {len(T2):,})")

    F1, F2 = tfs_of(e1), tfs_of(e2)
    j5 = jac(F1, F2)
    out["levels"]["L5_transcription_factors"] = {"n_1": len(F1), "n_2": len(F2), "jaccard": j5}
    log(f"L5 TFs                 : {j5:.4f}   ({len(F1)} vs {len(F2)})")

    # ---- the verdict ---------------------------------------------------------------------------
    vals = {k: v.get("jaccard", v.get("median_jaccard")) for k, v in out["levels"].items()}
    worst = min(v for v in vals.values() if v == v)
    out["ceiling_min_across_levels"] = worst
    if worst >= 0.999:
        verdict = ("DETERMINISTIC — every level reproduces exactly. The ceiling is 1.0 and every "
                   "arm-vs-arm difference in the ladder is attributable to the pairing.")
    elif worst >= 0.90:
        verdict = ("NEAR-DETERMINISTIC — run-to-run noise is small relative to the arm-vs-arm gaps "
                   "the ladder reports. Report the ceiling beside them.")
    else:
        verdict = ("STOCHASTIC — run-to-run variability is material. EVERY topology claim in the "
                   "ladder must be restated against this ceiling before publication. Any arm-vs-arm "
                   "gap smaller than the gap between these two replicates is not interpretable.")
    out["verdict"] = verdict
    log("")
    log(f"CEILING (worst level)  : {worst:.4f}")
    log(f"VERDICT: {verdict}")

    dest = R / "seedrep_ceiling.json"
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    log(f"wrote {dest}")

if __name__ == "__main__":
    main()
