#!/usr/bin/env python
"""L6 ceiling: is the weak within-cell-type activity agreement a bad pairing, or a limit of the task?"""
import json
import os

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
BM = os.path.join(ROOT, "benchmarks/bmmc")
PB = os.path.join(ROOT, "benchmarks/pbmc")
OUT = os.path.join(ROOT, "analysis", "out", "l6_ceiling_test.json")
N_HVG, N_PC, SEED = 2000, 50, 0

def log(m):
    print(f"[l6ceil] {m}", flush=True)

# ------------------------------------------------------------------ shared machinery ----
def pcs_from_counts(X):
    """log-normalise -> variance-ranked HVGs -> z-score -> PCA. No arm is read here, so the
    space cannot favour one arm over another."""
    lib = np.asarray(X.sum(1)).ravel()
    lib[lib == 0] = 1
    X = X.tocsr().astype(np.float32)
    X.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    mean = np.asarray(X.mean(0)).ravel()
    var = np.maximum(np.asarray(X.multiply(X).mean(0)).ravel() - mean ** 2, 0)
    k = min(N_HVG, X.shape[1])
    hvg = np.sort(np.argsort(-var)[:k])
    D = np.asarray(X[:, hvg].todense(), dtype=np.float32)
    D -= D.mean(0)
    sd = D.std(0)
    sd[sd == 0] = 1
    D /= sd
    np.clip(D, -10, 10, out=D)
    p = PCA(min(N_PC, D.shape[1]), svd_solver="randomized", random_state=SEED)
    return p.fit_transform(D).astype(np.float32)

def within_median(A, T, labels, order):
    """The L6 estimator used for the manifest values: per-regulon Spearman inside each group, median of the grid.
    NaN where a column is constant (0 there would be a spurious measured value)."""
    W = []
    for c in order:
        m = labels == c
        a, t = A[m], T[m]
        row = np.full(a.shape[1], np.nan)
        for j in range(a.shape[1]):
            if np.ptp(a[:, j]) and np.ptp(t[:, j]):
                row[j] = stats.spearmanr(a[:, j], t[:, j]).statistic
        W.append(row)
    return float(np.nanmedian(np.vstack(W)))

def oracle_partner(P, strat):
    """Nearest other cell in RNA PC space, inside the same stratum."""
    partner = np.arange(len(strat))
    for s in np.unique(strat):
        idx = np.flatnonzero(strat == s)
        if idx.size < 2:
            continue                      # nothing else in the block; keeps its own row
        _, ind = NearestNeighbors(n_neighbors=2).fit(P[idx]).kneighbors(P[idx])
        partner[idx] = idx[ind[:, 1]]     # column 0 is the cell itself
    return partner

# ------------------------------------------------------------------------ datasets ----
def load_bmmc():
    with h5py.File(os.path.join(BM, "out", "eRegulon_AUC_bmmc_truth.h5mu"), "r") as f:
        ACT = f["mod/Gene_based/X"][:]
        acells = np.char.replace(f["obs/Cell"][:].astype(str), "___bmmc", "")
    with h5py.File(os.path.join(BM, "out", "bmmc_rna_counts.h5ad"), "r") as f:
        shp = tuple(f["X"].attrs["shape"])
        RNA = sparse.csr_matrix((f["X"]["data"][:], f["X"]["indices"][:], f["X"]["indptr"][:]),
                                shape=shp)
        rbc = f["obs"]["_index"][:].astype(str)
    meta = pd.read_csv(os.path.join(BM, "out", "b00_cell_metadata.csv"))
    ct = dict(zip(meta["h5ad_barcode"], meta["cell_type"]))
    dn = dict(zip(meta["h5ad_barcode"], meta["donor"]))
    rrow = {b: i for i, b in enumerate(rbc)}
    RNA = RNA[[rrow[b] for b in acells]]
    group = np.array([ct[b] for b in acells])
    strat = np.array([f"{ct[b]}|{dn[b]}" for b in acells])
    arms = {}
    for role, fn in (("floor", "arm_rand_celltype_donor_1to1_seed555.csv"),
                     ("observed", "arm_link_x_celltype_donor_1to1_seed555.csv")):
        df = pd.read_csv(os.path.join(BM, "out", fn),
                         usecols=["true_rna_barcode", "assigned_rna_barcode"])
        arms[role] = (df["true_rna_barcode"].to_numpy(), df["assigned_rna_barcode"].to_numpy())
    return dict(name="BMMC", ACT=ACT, cells=acells, RNA=RNA, group=group, strat=strat,
                arms=arms, stratum_def="cell type x donor", group_def="cell type",
                floor_key="bmmc.L6.within.rand_celltype_x_donor_1to1",
                observed_key="bmmc.L6.within.link_x_celltype_donor_1to1")

def load_pbmc():
    with h5py.File(os.path.join(PB, "out", "eRegulon_AUC_pbmc10k_truth.h5mu"), "r") as f:
        ACT = f["mod/Gene_based/X"][:]
        acells = np.char.replace(f["obs/Cell"][:].astype(str), "___pbmc_10k", "")
    # the 10x multiome filtered matrix; the Gene Expression modality only
    p = os.path.join(ROOT, "data/pbmc",
                     "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
    with h5py.File(p, "r") as f:
        m = f["matrix"]
        shp = tuple(m["shape"][:])                     # (features, cells) -- CSC by cell
        M = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shp)
        bcs = m["barcodes"][:].astype(str)
        ftype = m["features/feature_type"][:].astype(str)
    M = M[ftype == "Gene Expression", :].T.tocsr()     # cells x genes
    rrow = {b: i for i, b in enumerate(bcs)}
    miss = [b for b in acells if b not in rrow]
    if miss:
        raise SystemExit(f"PBMC: {len(miss)} activity barcodes absent from the 10x matrix")
    RNA = M[[rrow[b] for b in acells]]
    cd = pd.read_csv(os.path.join(PB, "work", "cell_data_pbmc10k.csv"))
    sub = dict(zip(cd["barcode"], cd["Classified_Celltype"]))
    group = np.array([sub[b] for b in acells])
    strat = group.copy()          # PBMC: Classified_Celltype is the subcluster
    arms = {}
    for role, fn in (("floor", "arm_rand_subcluster_1to1_10k_seed555.csv"),
                     ("observed", "arm_link_subcluster_1to1_10k_seed555.csv")):
        df = pd.read_csv(os.path.join(PB, "out", fn),
                         usecols=["true_barcode", "Linked_scRNA_barcode"])
        arms[role] = (df["true_barcode"].to_numpy(), df["Linked_scRNA_barcode"].to_numpy())
    return dict(name="PBMC", ACT=ACT, cells=acells, RNA=RNA, group=group, strat=strat,
                arms=arms, stratum_def="subcluster", group_def="subcluster (Classified_Celltype)",
                floor_key="pbmcx.L6.within.rand_subcluster_1to1",
                observed_key="pbmcx.L6.within.link_subcluster_1to1")

def run(cfg, manifest):
    n = cfg["name"]
    log(f"--- {n}: {cfg['ACT'].shape[0]:,} cells x {cfg['ACT'].shape[1]} regulons, "
        f"{len(np.unique(cfg['strat']))} strata ({cfg['stratum_def']})")
    P = pcs_from_counts(cfg["RNA"])
    log(f"{n}: RNA PCA built")
    row = {b: i for i, b in enumerate(cfg["cells"])}
    ACT = cfg["ACT"]

    # One cell set and one grouping for all three bars, so they are comparable. The arms do
    # not always cover the whole activity matrix: PBMC's arms are solved on the 10,881 cells
    # that clear the 50-cell subcluster floor, against 10,925 rows of activity. Scoring the
    # oracle over all 10,925 would use a different denominator (and a group present in the
    # activity matrix can be empty in the arm), so the oracle is evaluated on the arms' own rows.
    eval_rows = np.array([row[b] for b in cfg["arms"]["observed"][0]])
    other = np.array([row[b] for b in cfg["arms"]["floor"][0]])
    if not np.array_equal(np.sort(eval_rows), np.sort(other)):
        raise SystemExit(f"{n}: the two arms are scored on different cell sets; "
                         f"observed {len(eval_rows)}, floor {len(other)}")
    glab = cfg["group"][eval_rows]
    order = sorted(c for c in np.unique(glab) if (glab == c).sum() >= 20)
    log(f"{n}: evaluation set {len(eval_rows):,} cells, {len(order)} groups scored "
        f"(activity matrix holds {ACT.shape[0]:,})")

    partner = oracle_partner(P, cfg["strat"])
    self_hits = int((partner == np.arange(len(partner))).sum())
    n_distinct = int(len(np.unique(partner)))
    res = {
        "stratum": cfg["stratum_def"], "l6_grouping": cfg["group_def"],
        "n_cells": int(ACT.shape[0]), "n_regulons": int(ACT.shape[1]),
        "n_groups_scored": len(order), "n_strata": int(len(np.unique(cfg["strat"]))),
        "oracle_self_assignments": self_hits,
        "oracle_n_distinct_rna_cells": n_distinct,
        "oracle_reuse_factor": float(ACT.shape[0] / n_distinct),
        "n_cells_evaluated": int(len(eval_rows)),
        "ceiling_rna_nn_within_stratum": within_median(ACT[partner[eval_rows]], ACT[eval_rows],
                                                       glab, order),
    }
    log(f"{n}: CEILING {res['ceiling_rna_nn_within_stratum']:+.4f}  "
        f"(oracle reuse {res['oracle_reuse_factor']:.4f}x, {self_hits} self-assignments)")

    for role in ("floor", "observed"):
        own_bc, got_bc = cfg["arms"][role]
        own = np.array([row[b] for b in own_bc])
        got = np.array([row[b] for b in got_bc])
        res[role] = within_median(ACT[got], ACT[own], cfg["group"][own], order)
        k = cfg[f"{role}_key"]
        mv = manifest.get(k, {}).get("value")
        res[f"{role}_manifest_key"] = k
        res[f"{role}_manifest_value"] = mv
        # the recomputation is checked against the manifest value to confirm it is the same estimator
        agree = mv is not None and abs(res[role] - mv) < 5e-4
        res[f"{role}_reproduces_manifest"] = bool(agree)
        log(f"{n}: {role:8s} {res[role]:+.4f}  vs manifest {mv}  "
            f"{'OK' if agree else '*** DISAGREES ***'}")

    c, f_, o = res["ceiling_rna_nn_within_stratum"], res["floor"], res["observed"]
    res.update({"unreachable_share_of_0_1": 1 - c,
                "observed_share_of_ceiling": o / c,
                "observed_on_floor_to_ceiling_scale": (o - f_) / (c - f_),
                "floor_share_of_ceiling": f_ / c,
                "matched_contrast_observed_minus_floor": o - f_})
    return res

def main():
    manifest = json.load(open(os.path.join(ROOT, "figures", "data", "manifest.json"),
                              encoding="utf-8"))
    out = {"_provenance": {
        "estimator": "b_06/p_08 L6: per-regulon Spearman within group, median of the grid",
        "rna_space": f"RNA-only PCA, {N_HVG} HVG, {N_PC} PCs, seed {SEED}, built before any arm",
        "oracle": "nearest OTHER cell in RNA PC space within the arm's own stratum",
        "bmmc_rna": "benchmarks/bmmc/out/bmmc_rna_counts.h5ad",
        "pbmc_rna": "data/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
                    " (Gene Expression modality)",
        "bmmc_activity": "benchmarks/bmmc/out/eRegulon_AUC_bmmc_truth.h5mu (Gene_based)",
        "pbmc_activity": "benchmarks/pbmc/out/eRegulon_AUC_pbmc10k_truth.h5mu (Gene_based)"}}
    for loader in (load_bmmc, load_pbmc):
        cfg = loader()
        out[cfg["name"]] = run(cfg, manifest)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w", encoding="utf-8"), indent=1)
    log(f"wrote {os.path.relpath(OUT, ROOT)}")
    print("\n  dataset   floor     observed  ceiling   oracle reuse")
    for d in ("BMMC", "PBMC"):
        r = out[d]
        print(f"  {d:9s} {r['floor']:+.4f}  {r['observed']:+.4f}  "
              f"{r['ceiling_rna_nn_within_stratum']:+.4f}   {r['oracle_reuse_factor']:.4f}x")

if __name__ == "__main__":
    main()
