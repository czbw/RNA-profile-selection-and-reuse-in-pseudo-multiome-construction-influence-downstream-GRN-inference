"Score native reference's regulon vocabulary on every arm's paired data (exact gather from the native reference AUCell matrix)."
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')

ARMS = ["link_1to1", "rand_subcluster_matched", "rand_subcluster_1to1", "rand_lineage_matched",
        "rand_compartment_matched", "rand_global"]

def log(m):
    print(f"[p08] {m}", flush=True)

def read_truth(h5mu_path):
    with h5py.File(h5mu_path, "r") as f:
        g = f["mod"]["Gene_based"]
        X = g["X"][:]
        cells = [c.decode() if isinstance(c, bytes) else str(c) for c in g["obs"]["Cell"][:]]
        regs = [r.decode() if isinstance(r, bytes) else str(r) for r in g["var"]["_index"][:]]
    return X, cells, regs

def spearman_cols(A, B):
    """Per-column Spearman between two (n x p) matrices. NaN where a column is constant."""
    out = np.full(A.shape[1], np.nan)
    for j in range(A.shape[1]):
        a, b = A[:, j], B[:, j]
        if np.ptp(a) == 0 or np.ptp(b) == 0:
            continue
        out[j] = stats.spearmanr(a, b).statistic
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="10k")
    ap.add_argument("--seeds", type=int, nargs="+", default=[555, 556, 557])
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    OUT = Path(a.outdir) if a.outdir else ROOT / "out"

    X, cells, regs = read_truth(ROOT / "out" / f"eRegulon_AUC_pbmc{a.dataset}_truth.h5mu")
    log(f"truth activity matrix: {X.shape[0]:,} cells x {X.shape[1]} regulons")

    # native reference's rows are keyed by cistopic_id (barcode___sample); arms name donors by bare barcode
    bare = [c.split("___")[0] for c in cells]
    row_of = {b: i for i, b in enumerate(bare)}
    assert len(row_of) == len(bare), "duplicate barcodes in truth matrix"

    manifest = pd.read_csv(ROOT / "work" / f"cell_data_pbmc{a.dataset}.csv")
    if "cistopic_id" not in manifest.columns:
        manifest = manifest.rename(columns={"cell_id": "cistopic_id"})
    ct = dict(zip(manifest["true_barcode"], manifest["Classified_Celltype"]))

    results, per_reg_dump = {}, []
    for seed in a.seeds:
        for arm in ARMS + ["link"]:
            if arm == "link":
                df = manifest[["cistopic_id", "Linked_scRNA_barcode", "true_barcode"]].copy()
                # restrict link to the same analysis set the arms use
                ref = pd.read_csv(OUT / f"arm_rand_global_{a.dataset}_seed{seed}.csv")
                df = df[df["cistopic_id"].isin(set(ref["cistopic_id"]))].reset_index(drop=True)
            else:
                df = pd.read_csv(OUT / f"arm_{arm}_{a.dataset}_seed{seed}.csv")

            # `link` drew partners from the full RNA object, so some donors have no row in native reference's
            # scored matrix and are unscoreable; drop those cells and report the count.
            # The rand_* arms draw from the analysis set and must have no missing donors.
            missing = set(df["Linked_scRNA_barcode"]) - set(row_of)
            n_before = len(df)
            if missing:
                assert arm == "link", f"{arm}: {len(missing)} donors absent -- rand_* must be clean"
                df = df[~df["Linked_scRNA_barcode"].isin(missing)].reset_index(drop=True)
                log(f"    [{arm}] {len(missing)} donors absent from truth's matrix -> dropped "
                    f"{n_before - len(df):,} of {n_before:,} cells ({1 - len(df)/n_before:.2%})")
            n_dropped = n_before - len(df)

            own = np.array([row_of[b] for b in df["true_barcode"]])
            got = np.array([row_of[b] for b in df["Linked_scRNA_barcode"]])
            T, A = X[own], X[got]                      # native reference vs arm, same cells, same vocabulary
            labels = np.array([ct[b] for b in df["true_barcode"]])

            rho_all = spearman_cols(A, T)
            within = {}
            for c in np.unique(labels):
                m = labels == c
                if m.sum() < 20:
                    continue
                within[c] = spearman_cols(A[m], T[m])
            W = np.vstack([within[c] for c in sorted(within)])   # celltypes x regulons

            res = {
                "arm": arm, "seed": seed, "n_cells": int(len(df)), "n_cells_dropped_unscoreable": int(n_dropped),
                "n_regulons": int(X.shape[1]), "n_celltypes": int(len(within)),
                "rho_all_median": float(np.nanmedian(rho_all)),
                "rho_within_median": float(np.nanmedian(W)),
                "rho_within_q25": float(np.nanquantile(W, .25)),
                "rho_within_q75": float(np.nanquantile(W, .75)),
                "frac_within_above_0.5": float(np.nanmean(W > .5)),
                "frac_within_below_0.2": float(np.nanmean(W < .2)),
                "attenuation": float(np.nanmedian(rho_all) - np.nanmedian(W)),
                "self_hit_rate": float((df["Linked_scRNA_barcode"] == df["true_barcode"]).mean()),
            }
            results[f"{arm}|{seed}"] = res
            per_reg_dump.append(pd.DataFrame({
                "arm": arm, "seed": seed, "regulon": regs,
                "rho_all": rho_all, "rho_within_median": np.nanmedian(W, axis=0)}))
            log(f"seed {seed}  {arm:26s} rho_all {res['rho_all_median']:+.3f} | "
                f"within {res['rho_within_median']:+.3f} "
                f"[{res['rho_within_q25']:+.3f},{res['rho_within_q75']:+.3f}] | "
                f"attenuation {res['attenuation']:.3f} | selfhit {res['self_hit_rate']:.4%}")

    pd.concat(per_reg_dump).to_csv(OUT / f"trackA_L6_per_regulon_{a.dataset}.csv", index=False)
    with open(OUT / f"trackA_L6_{a.dataset}.json", "w") as fh:
        json.dump(results, fh, indent=2)

    # ---- pre-specified primary comparison: link vs rand_subcluster_matched, within cell type ----
    log("")
    log("=" * 78)
    log("Primary prediction: link ~= rand_subcluster_matched (within cell type)")
    log("=" * 78)
    lk = [results[f"link|{s}"]["rho_within_median"] for s in a.seeds]
    or_ = [results[f"rand_subcluster_matched|{s}"]["rho_within_median"] for s in a.seeds]
    log(f"  link                     within-celltype rho = {np.mean(lk):+.4f}  (seeds: "
        f"{', '.join(f'{v:+.4f}' for v in lk)})")
    log(f"  rand_subcluster_matched  within-celltype rho = {np.mean(or_):+.4f}  (seeds: "
        f"{', '.join(f'{v:+.4f}' for v in or_)})")
    log(f"  margin (link - oracle)                       = {np.mean(lk) - np.mean(or_):+.4f}")
    log(f"  oracle spread across seeds (sd)              = {np.std(or_, ddof=1):.4f}")
    log("")
    log("  ladder, within-cell-type rho by rung:")
    for arm in ["link"] + ARMS:
        v = [results[f"{arm}|{s}"]["rho_within_median"] for s in a.seeds]
        log(f"    {arm:26s} {np.mean(v):+.4f}")
    log(f"\nwrote {OUT}/trackA_L6_{a.dataset}.json")

if __name__ == "__main__":
    main()
