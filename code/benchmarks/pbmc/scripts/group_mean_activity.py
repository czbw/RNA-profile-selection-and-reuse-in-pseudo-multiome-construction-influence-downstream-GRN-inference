'Cell-type mean-profile correlation of regulon activity, each arm vs native reference.'
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')
SEED_BOOT, B = 555, 2000

ARMS = ["link_1to1", "rand_subcluster_matched", "rand_subcluster_1to1", "rand_lineage_matched",
        "rand_compartment_matched", "rand_global"]

def log(m):
    print(f"[p09] {m}", flush=True)

def read_truth(p):
    with h5py.File(p, "r") as f:
        g = f["mod"]["Gene_based"]
        X = g["X"][:]
        cells = [c.decode() if isinstance(c, bytes) else str(c) for c in g["obs"]["Cell"][:]]
        regs = [r.decode() if isinstance(r, bytes) else str(r) for r in g["var"]["_index"][:]]
    return X, cells, regs

def profiles(M, labels, order):
    """(cells x regulons) -> (regulons x n_celltypes) matrix of within-cell-type MEAN activity."""
    return np.vstack([M[labels == c].mean(axis=0) for c in order]).T

def corr_matrix(Pa, Pt):
    'C[x, y] = Pearson(arm profile of regulon x, native reference profile of regulon y). Row-centred.'
    A = Pa - Pa.mean(axis=1, keepdims=True)
    T = Pt - Pt.mean(axis=1, keepdims=True)
    A /= np.linalg.norm(A, axis=1, keepdims=True)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    return A @ T.T

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="10k")
    ap.add_argument("--seeds", type=int, nargs="+", default=[555, 556, 557])
    a = ap.parse_args()
    OUT = ROOT / "out"

    X, cells, regs = read_truth(OUT / f"eRegulon_AUC_pbmc{a.dataset}_truth.h5mu")
    row_of = {c.split("___")[0]: i for i, c in enumerate(cells)}
    man = pd.read_csv(ROOT / "work" / f"cell_data_pbmc{a.dataset}.csv")
    if "cistopic_id" not in man.columns:
        man = man.rename(columns={"cell_id": "cistopic_id"})
    ct = dict(zip(man["true_barcode"], man["Classified_Celltype"]))
    log(f"truth matrix {X.shape[0]:,} x {X.shape[1]}")

    rng = np.random.default_rng(SEED_BOOT)
    results, rows = {}, []
    for seed in a.seeds:
        for arm in ARMS + ["link"]:
            if arm == "link":
                df = man[["cistopic_id", "Linked_scRNA_barcode", "true_barcode"]].copy()
                ref = pd.read_csv(OUT / f"arm_rand_global_{a.dataset}_seed{seed}.csv")
                df = df[df["cistopic_id"].isin(set(ref["cistopic_id"]))]
                df = df[df["Linked_scRNA_barcode"].isin(row_of)].reset_index(drop=True)
            else:
                df = pd.read_csv(OUT / f"arm_{arm}_{a.dataset}_seed{seed}.csv")

            own = np.array([row_of[b] for b in df["true_barcode"]])
            got = np.array([row_of[b] for b in df["Linked_scRNA_barcode"]])
            lab = np.array([ct[b] for b in df["true_barcode"]])
            order = sorted(c for c in np.unique(lab) if (lab == c).sum() >= 20)

            Pt = profiles(X[own], lab, order)          # native reference's cell-type profiles
            Pa = profiles(X[got], lab, order)          # the arm's
            keep = (Pt.std(axis=1) > 0) & (Pa.std(axis=1) > 0)
            C = corr_matrix(Pa[keep], Pt[keep])
            n = C.shape[0]
            matched = np.diag(C)
            off = C[~np.eye(n, dtype=bool)]

            # identifiability: percentile of the true match within its own row (chance = 50)
            pct = np.array([stats.percentileofscore(C[i], C[i, i], kind="rank") for i in range(n)])
            # regulon-clustered bootstrap on the matched median
            bs = np.array([np.median(matched[rng.integers(0, n, n)]) for _ in range(B)])

            r = {"arm": arm, "seed": seed, "n_regulons_used": int(n),
                 "n_regulons_dropped_constant": int((~keep).sum()),
                 "n_celltypes": len(order), "n_cells": int(len(df)),
                 "matched_median": float(np.median(matched)),
                 "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                 "mismatched_null_median": float(np.median(off)),
                 "null_p95": float(np.percentile(off, 95)),
                 "identifiability_median_pct": float(np.median(pct)),
                 "own_is_rank1_frac": float(np.mean(pct >= 100 * (n - 1) / n)),
                 "frac_matched_above_0.5": float(np.mean(matched > 0.5))}
            results[f"{arm}|{seed}"] = r
            rows.append(pd.DataFrame({"arm": arm, "seed": seed,
                                      "regulon": np.array(regs)[keep], "r_profile": matched}))
            log(f"seed {seed} {arm:26s} matched {r['matched_median']:+.4f} "
                f"[{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}] | null {r['mismatched_null_median']:+.4f} "
                f"| ident {r['identifiability_median_pct']:.1f}% | >0.5 {r['frac_matched_above_0.5']:.1%}")

    pd.concat(rows).to_csv(OUT / f"stat3_per_regulon_{a.dataset}.csv", index=False)
    with open(OUT / f"stat3_celltype_profile_{a.dataset}.json", "w") as fh:
        json.dump(results, fh, indent=2)

    log("")
    log("=" * 86)
    log("Cell-type mean-profile correlation vs truth (same resolution as the AML analysis)")
    log("=" * 86)
    log(f"{'arm':28s}{'matched':>10}{'95% CI':>20}{'null':>9}{'ident%':>9}{'>0.5':>8}")
    for arm in ["link"] + ARMS:
        v = [results[f"{arm}|{s}"] for s in a.seeds]
        m = np.mean([x["matched_median"] for x in v])
        lo = np.mean([x["ci95"][0] for x in v]); hi = np.mean([x["ci95"][1] for x in v])
        nu = np.mean([x["mismatched_null_median"] for x in v])
        idp = np.mean([x["identifiability_median_pct"] for x in v])
        f5 = np.mean([x["frac_matched_above_0.5"] for x in v])
        log(f"{arm:28s}{m:>+10.4f}{f'[{lo:+.4f},{hi:+.4f}]':>20}{nu:>+9.4f}{idp:>8.1f}%{f5:>8.1%}")
    log("")
    log("AML comparator (Option A vs Option B, no ground truth): healthy 0.844/0.900/0.862, "
        "malignant 0.820/0.710")

if __name__ == "__main__":
    main()
