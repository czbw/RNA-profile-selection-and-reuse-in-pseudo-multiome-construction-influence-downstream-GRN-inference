#!/usr/bin/env python
"""Stratum-restricted CCA linkage arms for BMMC."""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import queue
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"

# (arm token, human description, columns of b00_cell_metadata.csv whose product is the block)
BLOCK_DEFS = [
    ("link_x_donor", "donor", ["donor"]),
    ("link_x_celltype_donor", "cell type x donor", ["cell_type", "donor"]),
]

COST_FORMS = ("raw", "rescaled", "quantized")
QUANT_GRID = 1e9

def log(m): print(f"[b04b] {m}", flush=True)

def read_h5ad(p: Path):
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        M = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
        obs = np.array([x.decode() for x in f["obs/_index"][:]], dtype=object)
        var = np.array([x.decode() for x in f["var/_index"][:]], dtype=object)
    return M, obs, var

def read_obs(p: Path) -> np.ndarray:
    with h5py.File(p, "r") as f:
        return np.array([x.decode() for x in f["obs/_index"][:]], dtype=object)

def norm_log(M: sp.csr_matrix, target=1e4) -> sp.csr_matrix:
    """Library-size normalise to `target` counts per cell, then log1p. Seurat's LogNormalize."""
    tot = np.asarray(M.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    M = M.multiply((target / tot)[:, None]).tocsr()
    M.data = np.log1p(M.data)
    return M

def block_vector(md_indexed: pd.DataFrame, order: np.ndarray, cols: list) -> np.ndarray:
    """The block label for each cell of `order`, built by joining that cell order to the
    metadata. Called once per modality so the two sides are derived independently."""
    sub = md_indexed.loc[order, cols]
    v = sub[cols[0]].astype(str)
    for c in cols[1:]:
        v = v + "|" + sub[c].astype(str)
    return v.to_numpy()

def build_cca_space(npc: int, nhvg: int, seed: int):
    """The infer_pairings.py shared space, unchanged. -> (Pa, Pr, obs, S, shared_genes).

    Factored out so check_assignment_solvers.py builds the same space from the same code
    rather than from a copy.
    """
    t0 = time.time()
    RNA, obs_r, var_r = read_h5ad(OUT / "bmmc_rna_counts.h5ad")
    GA, obs_a, var_a = read_h5ad(OUT / "bmmc_atac_gene_activity.h5ad")
    if not np.array_equal(obs_r, obs_a):
        raise SystemExit("RNA and gene-activity cell order differ — the truth pairing would "
                         "be wrong")
    n = len(obs_r)
    shared = np.intersect1d(var_r, var_a)
    log(f"{n:,} cells | RNA {RNA.shape[1]:,} genes | gene-activity {GA.shape[1]:,} genes | "
        f"shared vocabulary {len(shared):,}")
    RNA = RNA[:, pd.Index(var_r).get_indexer(shared)]
    GA = GA[:, pd.Index(var_a).get_indexer(shared)]
    RNA = norm_log(RNA); GA = norm_log(GA)

    # variable genes chosen on the REFERENCE only, as Seurat does
    mu = np.asarray(RNA.mean(axis=0)).ravel()
    sq = np.asarray(RNA.multiply(RNA).mean(axis=0)).ravel()
    hv = np.argsort(np.maximum(sq - mu ** 2, 0.0))[::-1][:nhvg]
    log(f"variable genes selected on the reference: {len(hv):,}")
    RNA = RNA[:, hv]; GA = GA[:, hv]

    # each dataset scaled by its own feature statistics — Seurat's CCA path
    def scale(M):
        X = M.toarray().astype(np.float32)
        m = X.mean(0); s = X.std(0); s[s < 1e-9] = 1.0
        X = (X - m) / s
        np.clip(X, -10, 10, out=X)
        return X
    A = scale(RNA); B = scale(GA); del RNA, GA
    log(f"scaled reference {A.shape} and query {B.shape} | {time.time()-t0:.0f}s")

    # CCA as the top singular vectors of A @ B.T, never materialised
    op = LinearOperator((n, n), matvec=lambda v: A @ (B.T @ v),
                        rmatvec=lambda u: B @ (A.T @ u), dtype=np.float32)
    U, S, Vt = svds(op, k=npc, random_state=seed)
    o = np.argsort(S)[::-1]; U, S, Vt = U[:, o], S[o], Vt[o]
    Pr = np.ascontiguousarray(U.astype(np.float32))
    Pa = np.ascontiguousarray(Vt.T.astype(np.float32))
    del A, B
    log(f"CCA: {npc} components via matrix-free SVD, top singular value {S[0]:.6g} "
        f"| {time.time()-t0:.0f}s")
    Pr /= np.maximum(np.linalg.norm(Pr, axis=1, keepdims=True), 1e-12)   # Seurat L2Dim
    Pa /= np.maximum(np.linalg.norm(Pa, axis=1, keepdims=True), 1e-12)
    return Pa, Pr, obs_r, S, shared

def gate_reproduces_committed_link(Pa, Pr, obs_r, seed: int) -> dict:
    ''
    n = len(obs_r)
    link_csv = OUT / f"arm_link_seed{seed}.csv"
    if not link_csv.is_file():
        raise SystemExit(f"MISSING {link_csv} — the reproduction gate cannot run, and without "
                         f"it there is no evidence the shared space matches the committed "
                         f"link arms")
    ref = pd.read_csv(link_csv)
    rcol = "rna_barcode" if "rna_barcode" in ref.columns else "assigned_rna_barcode"
    if not np.array_equal(ref["atac_barcode"].to_numpy().astype(object), obs_r):
        raise SystemExit(f"{link_csv.name} row order differs from the h5ad cell order — the "
                         f"reproduction gate would compare the wrong cells")
    ref_idx = pd.Index(obs_r).get_indexer(ref[rcol].to_numpy().astype(object))
    if (ref_idx < 0).any():
        raise SystemExit(f"{link_csv.name} names RNA barcodes absent from the h5ad cell set")
    greedy_full = NearestNeighbors(n_neighbors=1).fit(Pr).kneighbors(Pa)[1][:, 0]
    n_disagree = int((greedy_full != ref_idx).sum())
    if n_disagree:
        raise SystemExit(
            f"REPRODUCTION CHECK FAILED: the re-derived CCA space disagrees with the committed "
            f"{link_csv.name} on {n_disagree:,} of {n:,} cells ({n_disagree/n:.4%}). The "
            f"restricted arms would then not be one-variable contrasts against link_1to1. "
            f"Investigate before shipping anything.")
    log(f"reproduction gate PASSED: re-derived greedy nearest neighbour equals the committed "
        f"{link_csv.name} on all {n:,} cells")
    return {"checked_against": link_csv.name,
            "statistic": "per-cell equality of the UNRESTRICTED greedy nearest neighbour "
                         "implied by the re-derived space",
            "n_cells_compared": int(n), "n_disagreements": 0, "passed": True}

def _matching_worker(tq, rq):
    """Long-lived child: pull a cost matrix, return the matching. See BudgetedMatcher."""
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching as _mwfbm
    while True:
        C = tq.get()
        if C is None:
            return
        try:
            r_i, c_i = _mwfbm(C)
            rq.put(("ok", (r_i, c_i)))
        except ValueError as e:               # "no full matching exists" — a Hall violation
            rq.put(("infeasible", str(e)))
        except Exception as e:                # never let the worker die silently
            rq.put(("error", f"{type(e).__name__}: {e}"))

class BudgetedMatcher:
    """min_weight_full_bipartite_matching, run in a worker process under a wall-clock budget."""

    def __init__(self, budget_s: float):
        self.budget = budget_s
        self.ctx = mp.get_context("spawn")
        self.p = self.tq = self.rq = None
        self.n_restarts = 0

    def _start(self):
        self.tq = self.ctx.Queue(); self.rq = self.ctx.Queue()
        self.p = self.ctx.Process(target=_matching_worker, args=(self.tq, self.rq), daemon=True)
        self.p.start()

    def solve(self, C):
        """-> ("ok", (rows, cols)) | ("infeasible", msg) | ("error", msg) | ("stalled", None)."""
        if self.p is None or not self.p.is_alive():
            self._start()
        self.tq.put(C)
        try:
            # the child blocks in put() until the pipe is drained, so read before joining
            return self.rq.get(timeout=self.budget)
        except queue.Empty:
            self.p.terminate(); self.p.join(); self.p = None
            self.n_restarts += 1
            return ("stalled", None)

    def close(self):
        if self.p is not None and self.p.is_alive():
            self.tq.put(None); self.p.join(timeout=10)
            if self.p.is_alive():
                self.p.terminate()
        self.p = None

def solve_blocks(Pa, Pr, blk, levels, knn, rev_knn, rng, tag, matcher):
    """Greedy and bijective assignment, per block. Returns (greedy, bij, per_block, totals)."""
    n = len(blk)
    greedy = np.empty(n, dtype=np.int64)
    bij = np.empty(n, dtype=np.int64)
    per_block, tot_edges, any_fb, tot_fb_used = {}, 0, False, 0
    n_stalls, n_conditioned, n_escalated = 0, 0, 0
    t0 = time.time()
    for li, lv in enumerate(levels):
        ix = np.flatnonzero(blk == lv)      # same index set on both sides (square, asserted)
        m = len(ix)
        Pa_b, Pr_b = Pa[ix], Pr[ix]
        kf0 = min(knn, m)
        kr = min(rev_knn, m)
        nn_f = NearestNeighbors(n_neighbors=kf0, algorithm="auto").fit(Pr_b)
        d_f, i_f = nn_f.kneighbors(Pa_b)
        greedy[ix] = ix[i_f[:, 0]]
        nn_r = NearestNeighbors(n_neighbors=kr).fit(Pa_b)
        d_r, i_r = nn_r.kneighbors(Pr_b)

        def build(kf, fallback, form):
            d2, i2 = (d_f, i_f) if kf == kf0 else nn_f.kneighbors(Pa_b, n_neighbors=kf)
            rows = np.repeat(np.arange(m), kf); cols = i2.ravel()
            vals = d2.ravel().astype(np.float64)
            rows = np.concatenate([rows, i_r.ravel()])
            cols = np.concatenate([cols, np.repeat(np.arange(m), kr)])
            vals = np.concatenate([vals, d_r.ravel().astype(np.float64)])
            if fallback:
                perm = rng.permutation(m)
                big = float(vals.max()) * 1e3
                rows = np.concatenate([rows, np.arange(m)])
                cols = np.concatenate([cols, perm])
                vals = np.concatenate([vals, np.full(m, big)])
            key = rows.astype(np.int64) * m + cols
            _, first = np.unique(key, return_index=True)
            rows, cols, vals = rows[first], cols[first], vals[first]
            lo, hi = float(vals.min()), float(vals.max())
            rg = max(hi - lo, 1e-12)
            if form == "rescaled":            # affine, EXACT: argmin provably unchanged
                vals = (vals - lo) / rg
            elif form == "quantized":         # affine + rounding onto a 1e9 grid
                vals = np.round((vals - lo) / rg * QUANT_GRID)
            # costs are shifted off zero: the routine treats an explicit zero as an edge
            return sp.csr_matrix((vals + 1e-12, (rows, cols)), shape=(m, m))

        sol, used_k, used_fb, used_form, n_fb_used, nnz = None, None, False, "raw", 0, 0
        attempts = []
        # Same order as infer_pairings.py: exhaust the true-cost graph at every k before
        # using the prohibitive-cost permutation. Within a k, try the cost forms in
        # increasing order of deviation before escalating k, because a stall is an
        # arithmetic artefact and not a property of the graph.
        tried = []
        for fb in (False, True):
            for kf_raw in (kf0, kf0 * 2, kf0 * 4, kf0 * 8):
                kf = min(kf_raw, m)
                if (kf, fb) in tried:
                    continue
                tried.append((kf, fb))
                for form in COST_FORMS:
                    C = build(kf, fb, form)
                    kind, payload = matcher.solve(C)
                    attempts.append({"k": int(kf), "fallback": bool(fb),
                                     "cost_form": form, "outcome": kind})
                    if kind == "ok":
                        r_i, c_i = payload
                        sol = np.empty(m, dtype=np.int64); sol[r_i] = c_i
                        used_k, used_fb, used_form, nnz = kf, fb, form, int(C.nnz)
                        if fb:
                            w = np.asarray(C[r_i, c_i]).ravel()
                            n_fb_used = int((w > w.max() * 0.5).sum())
                        break
                    if kind == "infeasible":
                        # feasibility is a property of the graph, not of the cost values, so
                        # no other cost form can help; go straight to the next k
                        log(f"  [{tag}/{lv}] k_fwd={kf} fallback={fb}: infeasible; escalating")
                        break
                    if kind == "error":
                        raise SystemExit(f"block {lv}: matching worker error: {payload}")
                    n_stalls += 1
                    nxt = (COST_FORMS[COST_FORMS.index(form) + 1]
                           if form != COST_FORMS[-1] else None)
                    log(f"  [{tag}/{lv}] k_fwd={kf} fallback={fb} costs={form}: solver did not "
                        f"return within {matcher.budget:.0f}s; "
                        + (f"retrying with `{nxt}` costs" if nxt else "escalating k"))
                if sol is not None:
                    break
            if sol is not None:
                break
        if sol is None:
            raise SystemExit(f"block {lv}: no perfect matching even with a permutation fallback "
                             f"— impossible for a square block; investigate the cost graph")
        bij[ix] = ix[sol]
        any_fb |= used_fb
        tot_fb_used += n_fb_used
        tot_edges += nnz
        n_conditioned += int(used_form != "raw")
        n_escalated += int(used_k != kf0)
        per_block[str(lv)] = {"n_cells": int(m), "forward_knn_used": int(used_k),
                              "reverse_knn": int(kr), "n_edges": nnz,
                              "permutation_fallback_needed": bool(used_fb),
                              "fallback_edges_in_solution": int(n_fb_used),
                              "cost_form": used_form,
                              "solver_attempts": attempts}
        if used_fb or used_form != "raw" or used_k != kf0 or m >= 5000 or li % 40 == 0:
            log(f"  [{tag}/{lv}] n={m:,} k_fwd={used_k} edges={nnz:,} fallback={used_fb}"
                + (f" ({n_fb_used:,} prohibitive edges used)" if used_fb else "")
                + (f" [{used_form} costs]" if used_form != "raw" else "")
                + f" | +{time.time()-t0:.0f}s")
    return greedy, bij, per_block, {
        "n_blocks": len(levels),
        "n_edges_total": int(tot_edges),
        "permutation_fallback_needed_any_block": bool(any_fb),
        "fallback_edges_in_solution_total": int(tot_fb_used),
        "n_solver_stalls": int(n_stalls),
        "n_blocks_not_solved_on_raw_costs": int(n_conditioned),
        "n_blocks_needing_k_escalation": int(n_escalated),
        "solver_note": "a stall is scipy LAPJVsp failing to converge, not an infeasible "
                       "graph; see BudgetedMatcher and out/b04b_solver_equivalence.json",
    }

def main() -> int:
    ap = argparse.ArgumentParser()

    ap.add_argument("--knn", type=int, default=100, help="forward kNN width, within block")
    ap.add_argument("--rev-knn", type=int, default=15, help="reverse kNN width, within block")
    ap.add_argument("--npc", type=int, default=50)
    ap.add_argument("--nhvg", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--solver-budget", type=float, default=90.0,
                    help="seconds a single min_weight_full_bipartite_matching call may take "
                         "before it is treated as stalled; see BudgetedMatcher")
    a = ap.parse_args()
    t0 = time.time()

    md_raw = pd.read_csv(OUT / "b00_cell_metadata.csv", dtype=str).set_index("h5ad_barcode")
    obs_r = read_obs(OUT / "bmmc_rna_counts.h5ad")
    obs_a = read_obs(OUT / "bmmc_atac_gene_activity.h5ad")
    n = len(obs_r)
    absent_r = int(np.sum(~pd.Index(obs_r).isin(md_raw.index)))
    absent_a = int(np.sum(~pd.Index(obs_a).isin(md_raw.index)))
    if absent_r or absent_a:
        raise SystemExit(f"b00_cell_metadata.csv does not cover every cell: "
                         f"{absent_r} RNA / {absent_a} ATAC barcodes absent")

    # ---- block squareness, derived from the two cell orders SEPARATELY --------------------
    # The feasibility precondition for a within-block bijection. Asserted, not assumed;
    # counting one modality's cell order twice would prove nothing.
    squareness = {}
    for tok, desc, cols in BLOCK_DEFS:
        b_r = block_vector(md_raw, obs_r, cols)
        b_a = block_vector(md_raw, obs_a, cols)
        c_r = pd.Series(b_r).value_counts()
        c_a = pd.Series(b_a).value_counts()
        lv = sorted(set(c_r.index) | set(c_a.index))
        tab, bad = {}, []
        for x in lv:
            nr, na = int(c_r.get(x, 0)), int(c_a.get(x, 0))
            tab[str(x)] = {"n_atac": na, "n_rna": nr, "square": na == nr}
            if na != nr:
                bad.append((x, na, nr))
        sizes = np.array([tab[str(x)]["n_atac"] for x in lv])
        if bad:
            for x, na, nr in bad:
                log(f"  NON-SQUARE {desc} block {x}: {na} ATAC vs {nr} RNA")
            raise SystemExit(
                f"{len(bad)} of {len(lv)} `{desc}` blocks are not square, so a within-{desc} "
                f"bijection over all {n:,} cells is INFEASIBLE. Stopping rather than dropping "
                f"cells — a silent drop is the failure mode this bundle exists to avoid.")
        if int(sizes.sum()) != n:
            raise SystemExit(f"`{desc}` block sizes sum to {sizes.sum():,}, not {n:,} — "
                             f"the blocks are not a partition of the cell set")
        # Chance level for exact-cell recovery under a uniform within-block bijection: a cell
        # in a block of size m is matched to itself with probability 1/m, so the expectation
        # over cells is sum_b m_b*(1/m_b) / n = n_blocks / n — exactly one expected fixed
        # point per block, for any block size (arms.py makes the same point).
        chance = float(np.sum(sizes * (1.0 / sizes)) / n)
        squareness[tok] = {
            "block_definition": desc, "block_columns": cols, "n_blocks": len(lv),
            "all_blocks_square": True, "sizes_sum_to_n_cells": True,
            "min_block_size": int(sizes.min()), "max_block_size": int(sizes.max()),
            "median_block_size": float(np.median(sizes)),
            "n_blocks_le_10": int((sizes <= 10).sum()),
            "n_blocks_le_25": int((sizes <= 25).sum()),
            "cells_in_blocks_le_10": int(sizes[sizes <= 10].sum()),
            "cells_in_blocks_le_25": int(sizes[sizes <= 25].sum()),
            "random_within_block_expected_exact_accuracy": chance,
            "per_block": tab}
        log(f"squareness `{desc}`: {len(lv)} blocks, all square, sizes {sizes.min()}-"
            f"{sizes.max()} (median {np.median(sizes):.0f}), sum {sizes.sum():,} = {n:,} — "
            f"within-block bijection feasible | chance exact-cell accuracy {chance:.4%}")

    Pa, Pr, obs_chk, S, shared = build_cca_space(a.npc, a.nhvg, a.seed)
    if not np.array_equal(obs_chk, obs_r):
        raise SystemExit("the CCA space's cell order differs from the metadata join order")
    gate = gate_reproduces_committed_link(Pa, Pr, obs_r, a.seed)
    md = md_raw.loc[obs_r].reset_index()

    matcher = BudgetedMatcher(a.solver_budget)
    donor = md["donor"].to_numpy(); site = md["site"].to_numpy()
    ct = md["cell_type"].to_numpy(); batch = md["batch"].to_numpy()
    frag = md["fragment_barcode"].to_numpy()
    idn = np.arange(n)

    # the site structure inside each donor: what a within-donor constraint does and does NOT
    # hold fixed. Measured, because the donor-aware arm's composite caveat rests on it.
    site_structure = {}
    for d in sorted(set(donor.tolist())):
        ix = np.flatnonzero(donor == d)
        site_structure[str(d)] = {"n_cells": int(len(ix)),
                                  "n_sites": int(len(set(site[ix].tolist()))),
                                  "n_batches": int(len(set(batch[ix].tolist())))}
    multi_site_cells = int(sum(v["n_cells"] for v in site_structure.values() if v["n_sites"] > 1))

    rep = {"n_cells": int(n), "n_shared_genes": int(len(shared)), "n_hvg": int(a.nhvg),
           "n_pcs": int(a.npc), "top_singular_value": float(S[0]),
           "reduction": "cca (matrix-free svds) — identical construction to infer_pairings.py",
           "reproduction_gate": gate,
           "block_squareness": squareness,
           "donor_site_structure": site_structure,
           "cells_in_multi_site_donors": multi_site_cells,
           "cells_in_multi_site_donors_pct": float(multi_site_cells / n),
           "graph": {"forward_knn_requested": a.knn, "reverse_knn": a.rev_knn, "per_family": {}},
           "cost_conditioning": {
               "forms_tried_in_order": list(COST_FORMS),
               "rescaled": "(c - min c)/(max c - min c); EXACT — every full matching uses "
                           "exactly m edges, so a positive-slope affine map is strictly "
                           "increasing in the matching cost and cannot reorder matchings",
               "quantized": f"the rescaled costs rounded onto a {QUANT_GRID:.0e} grid; "
                            "approximate, total error over a full matching bounded by "
                            "m * 0.5/GRID of the cost range (~1.2e-5 raw cost units at the "
                            "largest block here, against optimal costs of order 1e4). "
                            "Deterministic; no RNG.",
               "empirical_check": "check_assignment_solvers.py re-solves every block on "
                                  "the raw and the conditioned costs and compares the two "
                                  "matchings cell by cell; verdict in "
                                  "out/b04b_solver_equivalence.json",
               "why_needed": "scipy LAPJVsp does not converge on some of these graphs — "
                             "measured, and not an infeasibility (maximum_bipartite_matching "
                             "finds a perfect matching on the same graph in ~0.02 s). See the "
                             "BudgetedMatcher docstring for the timings."},
           "seed": a.seed,
           "method": "infer_pairings.py's CCA linkage with the candidate edge set restricted to "
                     "within-block pairs; see module docstring",
           "arms": {}}

    for tok, desc, cols in BLOCK_DEFS:
        blk = block_vector(md_raw, obs_r, cols)
        levels = sorted(set(blk.tolist()))
        rng = np.random.default_rng(a.seed)   # per family, so each is reproducible on its own
        log(f"--- {desc}: {len(levels)} blocks ---")
        greedy, bij, per_block, totals = solve_blocks(
            Pa, Pr, blk, levels, a.knn, a.rev_knn, rng, tok, matcher)
        totals["per_block"] = per_block
        rep["graph"]["per_family"][tok] = totals

        # ---- assert what the arm claims, not a proxy ------------------------------------
        if not np.array_equal(np.sort(bij), idn):
            raise SystemExit(f"{tok}_1to1 is not a permutation of the cell set")
        for nm, asg in ((tok, greedy), (f"{tok}_1to1", bij)):
            if (blk[asg] != blk).any():
                raise SystemExit(f"{nm} left its block for {int((blk[asg] != blk).sum()):,} cells")
        # a within-block permutation must be a permutation of its own block, per block:
        # 'is injective' and 'stayed inside the block' established together and exactly
        for x in levels:
            ix = np.flatnonzero(blk == x)
            if not np.array_equal(np.sort(bij[ix]), np.sort(ix)):
                raise SystemExit(f"block {x} is NOT a permutation of its own cells")
        log(f"{desc}: invariants hold — bijection over all {n:,} cells, every pair inside its "
            f"block, every block a permutation of itself")

        for nm, asg in ((tok, greedy), (f"{tok}_1to1", bij)):
            cnt = np.bincount(asg, minlength=n)
            used = int((cnt > 0).sum())
            rep["arms"][nm] = {
                "block_definition": desc,
                "reuse_mean_per_used_rna_cell": float(n / used),
                "reuse_max": int(cnt.max()),
                "distinct_rna_cells": used,
                "pct_pool_used": float(used / n),
                "exact_cell_accuracy": float((asg == idn).mean()),
                "cross_donor_rate": float((donor != donor[asg]).mean()),
                "cross_site_rate": float((site != site[asg]).mean()),
                "cross_cell_type_rate": float((ct != ct[asg]).mean()),
                "random_within_block_expected_exact_accuracy":
                    squareness[tok]["random_within_block_expected_exact_accuracy"],
            }
            s = rep["arms"][nm]
            log(f"{nm:28} reuse {s['reuse_mean_per_used_rna_cell']:.3f}x (max {s['reuse_max']}) | "
                f"distinct {used:,} ({s['pct_pool_used']:.1%}) | exact-cell "
                f"{s['exact_cell_accuracy']:.4%} (chance "
                f"{s['random_within_block_expected_exact_accuracy']:.4%}) | cross-donor "
                f"{s['cross_donor_rate']:.2%} | cross-site {s['cross_site_rate']:.2%} | "
                f"cross-celltype {s['cross_cell_type_rate']:.2%}")

            # 14-column bundle assignment table (the arms.py schema).
            # build_arm_rna.py resolves `atac_barcode` -> `assigned_rna_barcode`; the
            # rest is provenance carried by the rand_* arms.
            # `subject` == `donor` in this bundle (arms.py: md["subject"] =
            # md["donor"]), so cross_subject is the donor test, not a second axis.
            frame = pd.DataFrame({
                "atac_barcode": obs_r,
                "atac_fragment_barcode": frag,
                "assigned_rna_barcode": obs_r[asg],
                "true_rna_barcode": obs_r,
                "stratum": blk.astype(str),
                "atac_donor": donor, "rna_donor": donor[asg],
                "atac_site": site, "rna_site": site[asg],
                "atac_cell_type": ct,
                "is_self": idn == asg,
                "cross_donor": donor != donor[asg],
                "cross_subject": donor != donor[asg],
                "cross_site": site != site[asg],
            })
            dest = OUT / f"arm_{nm}_seed{a.seed}.csv"
            frame.to_csv(dest, index=False, lineterminator="\n")
            rep["arms"][nm]["file"] = dest.name
            rep["arms"][nm]["row_unit"] = "one row per ATAC cell"

        # the definitional checks, stated as the definitions
        b = rep["arms"][f"{tok}_1to1"]
        if b["distinct_rna_cells"] != n or b["reuse_max"] != 1:
            raise SystemExit(f"{tok}_1to1 is not a bijection: {b['distinct_rna_cells']} distinct, "
                             f"max reuse {b['reuse_max']}")
        for nm in (tok, f"{tok}_1to1"):
            if rep["arms"][nm]["cross_donor_rate"] != 0.0:
                raise SystemExit(f"{nm} is not within-donor: cross-donor rate "
                                 f"{rep['arms'][nm]['cross_donor_rate']}")
            if "cell_type" in cols and rep["arms"][nm]["cross_cell_type_rate"] != 0.0:
                raise SystemExit(f"{nm} claims a within-cell-type restriction but its "
                                 f"cross-cell-type rate is "
                                 f"{rep['arms'][nm]['cross_cell_type_rate']} — the restriction "
                                 f"did not take and the arm is wrong")

    matcher.close()
    rep["graph"]["solver_budget_s"] = a.solver_budget
    rep["graph"]["worker_restarts_after_stall"] = int(matcher.n_restarts)

    ref_rep = json.loads((OUT / "b04_link_report.json").read_text(encoding="utf-8"))
    rep["reference_arms_from_b04_link_report"] = ref_rep["arms"]
    rep["reference_graph_from_b04_link_report"] = ref_rep["graph"]

    rep["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "b04b_link_restricted_report.json").write_text(json.dumps(rep, indent=1),
                                                          encoding="utf-8")
    log("wrote " + ", ".join(f"arm_{k}_seed{a.seed}.csv" for k in rep["arms"])
        + ", b04b_link_restricted_report.json")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
