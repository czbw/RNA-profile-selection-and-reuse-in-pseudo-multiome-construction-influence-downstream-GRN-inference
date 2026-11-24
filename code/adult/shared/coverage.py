"""Target-coverage checks from feasibility.py, recomputed for a given module."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

import scoring_common as C

COVERAGE_BAR = 0.90          # unchanged from feasibility.py

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    log = lambda m: print(f"[cov:{a.tag}] {m}", flush=True)

    M = C.resolve_module(a.module, a.arm)
    OUT = C.out_dir(a.arm)
    log(f"{M['name']}  {M['n_regulons']} regulons / {M['n_TFs']} TFs | arm vocabulary "
        f"'{M['arm']}' | background {len(M['background'])}")

    X = pd.read_csv(C.MAIN, sep="\t", index_col=0, compression="gzip")
    X = X[~X.index.duplicated(keep="first")]
    X.index = X.index.astype(str)
    genes = set(X.index)
    log(f"matrix {X.shape[0]:,} genes x {X.shape[1]} samples")

    res = {"provenance": C.provenance(M), "coverage_bar": COVERAGE_BAR,
           "n_genes_in_matrix": int(X.shape[0]), "by_vocabulary": {}}
    rows = []
    for v in C.DEFS:
        tmap = M["targets"][v]
        for r in M["regulons"]:
            tg = set(tmap[r])
            hit = len(tg & genes)
            rows.append(dict(vocabulary=v, regulon=r, n_targets=len(tg), n_present=hit,
                             coverage=hit / len(tg)))
        u = set().union(*[set(tmap[r]) for r in M["regulons"]])
        up = u & genes
        g = pd.DataFrame([x for x in rows if x["vocabulary"] == v])
        sub = X.loc[sorted(up)]
        res["by_vocabulary"][v] = {
            "n_regulons": int(len(g)),
            "median_coverage": float(g.coverage.median()),
            "min_coverage": float(g.coverage.min()),
            "frac_regulons_below_bar": float((g.coverage < COVERAGE_BAR).mean()),
            "passes_coverage_bar": bool(g.coverage.median() >= COVERAGE_BAR),
            "median_n_targets": float(g.n_targets.median()),
            "target_union_size": int(len(u)),
            "target_union_present_in_matrix": int(len(up)),
            "target_union_frac_of_transcriptome": float(len(up) / X.shape[0]),
            "module_genes_zero_everywhere": int((sub.sum(axis=1) == 0).sum()),
            "module_genes_median_detection_rate": float(((sub > 0).mean(axis=1)).median()),
        }
        x = res["by_vocabulary"][v]
        log(f"{v}: median per-regulon coverage {x['median_coverage']:.4f} (bar {COVERAGE_BAR}) -> "
            f"{'PASS' if x['passes_coverage_bar'] else 'FAIL'}; worst {x['min_coverage']:.3f}")
        log(f"{v}: TARGET UNION {x['target_union_size']:,} genes, {x['target_union_present_in_matrix']:,} "
            f"present = {x['target_union_frac_of_transcriptome']:.1%} of the {X.shape[0]:,}-gene "
            f"matrix")

    pd.DataFrame(rows).to_csv(OUT / f"target_coverage_{a.tag}.csv", index=False)
    (OUT / f"coverage_{a.tag}.json").write_text(json.dumps(res, indent=2))
    log(f"wrote {OUT.name}/target_coverage_{a.tag}.csv and coverage_{a.tag}.json")

if __name__ == "__main__":
    main()
