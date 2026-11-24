#!/usr/bin/env python
"""The `link` and `link_1to1` arms for BMMC: the construction the field publishes."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"

def log(m): print(f"[b04] {m}", flush=True)

def read_h5ad(p: Path):
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        M = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
        obs = np.array([x.decode() for x in f["obs/_index"][:]], dtype=object)
        var = np.array([x.decode() for x in f["var/_index"][:]], dtype=object)
    return M, obs, var

def norm_log(M: sp.csr_matrix, target=1e4) -> sp.csr_matrix:
    """Library-size normalise to `target` counts per cell, then log1p. Seurat's LogNormalize."""
    tot = np.asarray(M.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    M = M.multiply((target / tot)[:, None]).tocsr()
    M.data = np.log1p(M.data)
    return M

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--knn", type=int, default=200, help="forward kNN width for the bijection graph")
    ap.add_argument("--rev-knn", type=int, default=15, help="reverse kNN width")
    ap.add_argument("--npc", type=int, default=50)
    ap.add_argument("--nhvg", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=555)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    t0 = time.time()

    RNA, obs_r, var_r = read_h5ad(OUT / "bmmc_rna_counts.h5ad")
    GA, obs_a, var_a = read_h5ad(OUT / "bmmc_atac_gene_activity.h5ad")
    if not np.array_equal(obs_r, obs_a):
        raise SystemExit("RNA and gene-activity cell order differ — the truth pairing would be wrong")
    n = len(obs_r)
    md = pd.read_csv(OUT / "b00_cell_metadata.csv")
    md = md.set_index("h5ad_barcode").loc[obs_r].reset_index()
    log(f"{n:,} cells | RNA {RNA.shape[1]:,} genes | gene-activity {GA.shape[1]:,} genes")

    shared = np.intersect1d(var_r, var_a)
    log(f"shared vocabulary: {len(shared):,} genes")
    RNA = RNA[:, pd.Index(var_r).get_indexer(shared)]
    GA = GA[:, pd.Index(var_a).get_indexer(shared)]

    RNA = norm_log(RNA); GA = norm_log(GA)

    # variable genes chosen on the REFERENCE only, as Seurat does
    mu = np.asarray(RNA.mean(axis=0)).ravel()
    sq = np.asarray(RNA.multiply(RNA).mean(axis=0)).ravel()
    vr = np.maximum(sq - mu ** 2, 0.0)
    hv = np.argsort(vr)[::-1][:a.nhvg]
    log(f"variable genes selected on the reference: {len(hv):,}")
    RNA = RNA[:, hv]; GA = GA[:, hv]

    # Scale each dataset by its own feature statistics (Seurat's CCA path, not the
    # pcaproject path where the query is scaled with reference statistics).
    def scale(M):
        X = M.toarray().astype(np.float32)
        m = X.mean(0); s = X.std(0); s[s < 1e-9] = 1.0
        X = (X - m) / s
        np.clip(X, -10, 10, out=X)
        return X
    A = scale(RNA); B = scale(GA); del RNA, GA
    log(f"scaled reference {A.shape} and query {B.shape} | {time.time()-t0:.0f}s")

    # CCA as the top singular vectors of A @ B.T, never materialised (see docstring).
    op = LinearOperator((n, n), matvec=lambda v: A @ (B.T @ v),
                        rmatvec=lambda u: B @ (A.T @ u), dtype=np.float32)
    U, S, Vt = svds(op, k=a.npc, random_state=a.seed)
    o = np.argsort(S)[::-1]; U, S, Vt = U[:, o], S[o], Vt[o]
    Pr = np.ascontiguousarray(U.astype(np.float32))
    Pa = np.ascontiguousarray(Vt.T.astype(np.float32))
    del A, B
    log(f"CCA: {a.npc} components via matrix-free SVD, top singular value {S[0]:.3g} "
        f"| {time.time()-t0:.0f}s")

    # L2-normalise before neighbour search (Seurat L2Dim)
    Pr /= np.maximum(np.linalg.norm(Pr, axis=1, keepdims=True), 1e-12)
    Pa /= np.maximum(np.linalg.norm(Pa, axis=1, keepdims=True), 1e-12)

    # ---- link: greedy, independent per ATAC cell ------------------------------------
    nn = NearestNeighbors(n_neighbors=a.knn, algorithm="auto").fit(Pr)
    dist, idx = nn.kneighbors(Pa)
    greedy = idx[:, 0]
    log(f"greedy nearest-neighbour assignment done | {time.time()-t0:.0f}s")

    nn_r = NearestNeighbors(n_neighbors=a.rev_knn).fit(Pa)
    dr, ir = nn_r.kneighbors(Pr)

    def build(kf: int, fallback: bool):
        d2, i2 = (dist, idx) if kf == a.knn else nn.kneighbors(Pa, n_neighbors=kf)
        rows = np.repeat(np.arange(n), kf); cols = i2.ravel(); vals = d2.ravel().astype(np.float64)
        rows = np.concatenate([rows, ir.ravel()])
        cols = np.concatenate([cols, np.repeat(np.arange(n), a.rev_knn)])
        vals = np.concatenate([vals, dr.ravel().astype(np.float64)])
        n_fb = 0
        if fallback:
            perm = rng.permutation(n)
            big = float(vals.max()) * 1e3
            rows = np.concatenate([rows, np.arange(n)])
            cols = np.concatenate([cols, perm])
            vals = np.concatenate([vals, np.full(n, big)])
            n_fb = n
        key = rows.astype(np.int64) * n + cols
        _, first = np.unique(key, return_index=True)
        rows, cols, vals = rows[first], cols[first], vals[first]
        return sp.csr_matrix((vals + 1e-12, (rows, cols)), shape=(n, n)), n_fb

    bij, C, used_k, used_fb, n_fallback_used = None, None, None, False, 0
    # Order matters: try every k on the true-cost graph before using the permutation fallback,
    # so a denser real graph is preferred to the fallback.
    for fb in (False, True):
        for kf in [a.knn, a.knn * 2, a.knn * 4, a.knn * 8]:
            if kf > n: continue
            C, _ = build(kf, fb)
            log(f"bipartite graph: k_fwd={kf} rev={a.rev_knn} fallback={fb} "
                f"edges={C.nnz:,} — solving | {time.time()-t0:.0f}s")
            try:
                r_i, c_i = min_weight_full_bipartite_matching(C)
            except ValueError as e:
                log(f"  infeasible ({e}); escalating")
                continue
            bij = np.empty(n, dtype=np.int64); bij[r_i] = c_i
            used_k, used_fb = kf, fb
            if fb:
                thr = float(np.asarray(C[r_i, c_i]).ravel().max()) * 0.5
                n_fallback_used = int((np.asarray(C[r_i, c_i]).ravel() > thr).sum())
            break
        if bij is not None: break
    if bij is None:
        raise SystemExit("no perfect matching found even with a permutation fallback — "
                         "this should be impossible; investigate the cost graph")
    log(f"bijection solved with k_fwd={used_k}, fallback={used_fb}"
        + (f", {n_fallback_used:,} prohibitive-cost edges used ({n_fallback_used/n:.2%})" if used_fb else "")
        + f" | {time.time()-t0:.0f}s")

    # ---- measure --------------------------------------------------------------------
    donor = md["donor"].to_numpy(); site = md["site"].to_numpy(); ct = md["cell_type"].to_numpy()
    rep = {"n_cells": int(n), "n_shared_genes": int(len(shared)), "n_hvg": int(a.nhvg),
           "n_pcs": int(a.npc), "top_singular_value": float(S[0]), "reduction": "cca (matrix-free svds)",
           "graph": {"forward_knn_requested": a.knn, "forward_knn_used": int(used_k),
                     "reverse_knn": a.rev_knn, "n_edges": int(C.nnz),
                     "permutation_fallback_needed": bool(used_fb),
                     "fallback_edges_in_solution": int(n_fallback_used),
                     "fallback_note": "edges added at 1000x max cost purely to guarantee a "
                                      "perfect matching exists; a large count means the arm is "
                                      "not a nearest-neighbour bijection"},
           "seed": a.seed, "method": "Seurat pcaproject reimplemented; see module docstring",
           "arms": {}}
    for name, asg in (("link", greedy), ("link_1to1", bij)):
        cnt = np.bincount(asg, minlength=n)
        used = int((cnt > 0).sum())
        rep["arms"][name] = {
            "reuse_mean_per_used_rna_cell": float(n / used),
            "reuse_max": int(cnt.max()),
            "distinct_rna_cells": used,
            "pct_pool_used": float(used / n),
            "exact_cell_accuracy": float((asg == np.arange(n)).mean()),
            "cross_donor_rate": float((donor != donor[asg]).mean()),
            "cross_site_rate": float((site != site[asg]).mean()),
            "cross_cell_type_rate": float((ct != ct[asg]).mean()),
        }
        s = rep["arms"][name]
        log(f"{name:10} reuse {s['reuse_mean_per_used_rna_cell']:.3f}x (max {s['reuse_max']}) | "
            f"distinct {used:,} ({s['pct_pool_used']:.1%}) | exact-cell {s['exact_cell_accuracy']:.4%} | "
            f"cross-donor {s['cross_donor_rate']:.1%} | cross-site {s['cross_site_rate']:.1%} | "
            f"cross-celltype {s['cross_cell_type_rate']:.1%}")
        pd.DataFrame({"atac_barcode": obs_r, "rna_barcode": obs_r[asg]}).to_csv(
            OUT / f"arm_{name}_seed{a.seed}.csv", index=False)

    # a bijection that is not a permutation is a bug, so check it rather than trust it
    b = rep["arms"]["link_1to1"]
    if b["distinct_rna_cells"] != n or b["reuse_max"] != 1:
        raise SystemExit(f"link_1to1 is not a bijection: {b['distinct_rna_cells']} distinct, "
                         f"max reuse {b['reuse_max']}")
    rep["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "b04_link_report.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote arm_link_seed{a.seed}.csv, arm_link_1to1_seed{a.seed}.csv, b04_link_report.json")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
