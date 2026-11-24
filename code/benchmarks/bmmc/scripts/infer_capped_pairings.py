#!/usr/bin/env python
"""The capped (cell type x donor) CCA linkage arm for BMMC, at a forced reuse of ~1.4."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

import sys
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from infer_stratified_pairings import (  # noqa: E402
    COST_FORMS, QUANT_GRID, BudgetedMatcher, block_vector, build_cca_space,
    gate_reproduces_committed_link, read_obs)

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"

# The one block family this arm is defined on: prep_04b's BLOCK_DEFS entry, restated here rather
# than indexed by position, which could silently select the wrong stratum.
BLOCK_COLS = ["cell_type", "donor"]
BLOCK_DESC = "cell type x donor"

def log(m): print(f"[b04d] {m}", flush=True)

def capacity_vector(block_label: str, m: int, cap: int, frac: float, seed: int) -> np.ndarray:
    """c_j for each RNA cell of a block: `cap` for the first round(frac*m) cells of a fixed
    per-block permutation, 1 for the rest.

    The permutation depends only on (seed, block label) -- not on `frac`, not on the iteration
    order of the blocks -- so the cap-2 set is nested as frac rises and a sweep moves one number.
    hashlib, never hash(), which is salted per process and would not reproduce across runs."""
    h = int.from_bytes(hashlib.sha256(block_label.encode()).digest()[:8], "big")
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), h]))
    perm = rng.permutation(m)
    n_big = int(round(frac * m))
    c = np.ones(m, dtype=np.int64)
    c[perm[:n_big]] = int(cap)
    return c

def solve_block_capped(Pa_b, Pr_b, cap_vec, knn, rev_knn, matcher, tag, lv):
    """One block: build the kNN cost graph exactly as prep_04b does, replicate each RNA column"""
    m = len(cap_vec)
    kf0 = min(knn, m)
    kr = min(rev_knn, m)
    nn_f = NearestNeighbors(n_neighbors=kf0, algorithm="auto").fit(Pr_b)
    d_f, i_f = nn_f.kneighbors(Pa_b)
    nn_r = NearestNeighbors(n_neighbors=kr).fit(Pa_b)
    d_r, i_r = nn_r.kneighbors(Pr_b)

    # column layout: RNA cell j owns columns [off[j], off[j] + c_j)
    off = np.concatenate([[0], np.cumsum(cap_vec)]).astype(np.int64)
    n_capacity_slots = int(off[-1])
    rna_column_owner = np.repeat(np.arange(m, dtype=np.int64), cap_vec)

    def build(kf, form):
        d2, i2 = (d_f, i_f) if kf == kf0 else nn_f.kneighbors(Pa_b, n_neighbors=kf)
        rows = np.repeat(np.arange(m), kf)
        cols = i2.ravel()
        vals = d2.ravel().astype(np.float64)
        rows = np.concatenate([rows, i_r.ravel()])
        cols = np.concatenate([cols, np.repeat(np.arange(m), kr)])
        vals = np.concatenate([vals, d_r.ravel().astype(np.float64)])
        # dedup on the BASE (atac, rna) pair, before replication: prep_04b dedups the same key on
        # the same graph, so the base edge set here is the base edge set there
        key = rows.astype(np.int64) * m + cols
        _, first = np.unique(key, return_index=True)
        rows, cols, vals = rows[first], cols[first], vals[first]
        lo, hi = float(vals.min()), float(vals.max())
        rg = max(hi - lo, 1e-12)
        if form == "rescaled":            # affine, exact: argmin provably unchanged
            vals = (vals - lo) / rg
        elif form == "quantized":         # affine + rounding onto a 1e9 grid
            vals = np.round((vals - lo) / rg * QUANT_GRID)
        # replicate every edge onto all c_j columns of its RNA cell, at the same cost
        rep = cap_vec[cols]
        e_rows = np.repeat(rows, rep)
        base = np.repeat(off[cols], rep)
        within = np.arange(len(e_rows)) - np.repeat(np.cumsum(rep) - rep, rep)
        e_cols = base + within
        e_vals = np.repeat(vals, rep)
        # costs shifted off zero: the routine treats an explicit zero as a missing edge
        return sp.csr_matrix((e_vals + 1e-12, (e_rows, e_cols)), shape=(m, n_capacity_slots))

    sol, used_k, used_form, nnz, attempts = None, None, "raw", 0, []
    for kf_raw in (kf0, kf0 * 2, kf0 * 4, kf0 * 8):
        kf = min(kf_raw, m)
        if any(att["k"] == int(kf) for att in attempts):
            continue
        for form in COST_FORMS:
            Cm = build(kf, form)
            kind, payload = matcher.solve(Cm)
            attempts.append({"k": int(kf), "cost_form": form, "outcome": kind})
            if kind == "ok":
                r_i, c_i = payload
                sol = np.empty(m, dtype=np.int64)
                sol[r_i] = c_i
                used_k, used_form, nnz = kf, form, int(Cm.nnz)
                break
            if kind == "infeasible":
                # feasibility is a property of the graph, not of the cost values, so no other cost
                # form can help; go straight to the next k
                log(f"  [{tag}/{lv}] k_fwd={kf}: infeasible; escalating k")
                break
            if kind == "error":
                raise SystemExit(f"block {lv}: matching worker error: {payload}")
            nxt = (COST_FORMS[COST_FORMS.index(form) + 1]
                   if form != COST_FORMS[-1] else None)
            log(f"  [{tag}/{lv}] k_fwd={kf} costs={form}: solver did not return within "
                f"{matcher.budget:.0f}s; "
                + (f"retrying with `{nxt}` costs" if nxt else "escalating k"))
        if sol is not None:
            break
    if sol is None:
        # prep_04b's prohibitive-permutation fallback is deliberately absent: c_j >= 1 for every RNA
        # cell, so the expanded graph contains the square graph prep_04b already solved to a perfect
        # matching on this block. A stall at every k and cost form is a solver problem, not a graph
        # problem, and inventing edges would hide it.
        raise SystemExit(
            f"block {lv}: no full row matching returned at any k or cost form. The expanded graph "
            f"contains prep_04b's square within-block graph, which WAS solved, so this is the "
            f"solver failing to converge and not an infeasibility. Investigate; do not relax.")

    rna_local = rna_column_owner[sol]
    used_counts = np.bincount(rna_local, minlength=m)
    if (used_counts > cap_vec).any():
        bad = int((used_counts > cap_vec).sum())
        raise SystemExit(f"block {lv}: {bad} RNA cell(s) exceed their capacity — node duplication "
                         f"did not impose the cap and the arm is wrong")
    rec = {"n_cells": int(m), "n_columns": int(n_capacity_slots),
           "n_cells_with_capacity_2": int((cap_vec > 1).sum()),
           "forward_knn_used": int(used_k), "reverse_knn": int(kr), "n_edges": nnz,
           "cost_form": used_form,
           "distinct_rna_cells_used": int((used_counts > 0).sum()),
           "reuse_max": int(used_counts.max()),
           "solver_attempts": attempts}
    return rna_local, rec

def build_arm(Pa, Pr, blk, levels, cap, frac, a, matcher, tag):
    """The whole family at one capacity scheme. -> (assignment over all n cells, per_block, totals)."""
    n = len(blk)
    asg = np.empty(n, dtype=np.int64)
    per_block, tot_edges, n_cond, n_esc, tot_slots = {}, 0, 0, 0, 0
    t0 = time.time()
    for li, lv in enumerate(levels):
        ix = np.flatnonzero(blk == lv)          # same index set on both sides (square, asserted)
        m = len(ix)
        cap_vec = capacity_vector(str(lv), m, cap, frac, a.seed)
        rna_local, rec = solve_block_capped(Pa[ix], Pr[ix], cap_vec, a.knn, a.rev_knn,
                                            matcher, tag, lv)
        asg[ix] = ix[rna_local]
        per_block[str(lv)] = rec
        tot_edges += rec["n_edges"]
        tot_slots += rec["n_columns"]
        n_cond += int(rec["cost_form"] != "raw")
        n_esc += int(rec["forward_knn_used"] != min(a.knn, m))
        if rec["cost_form"] != "raw" or rec["forward_knn_used"] != min(a.knn, m) \
                or m >= 2000 or li % 40 == 0:
            log(f"  [{tag}/{lv}] n={m:,} cols={rec['n_columns']:,} k_fwd={rec['forward_knn_used']} "
                f"edges={rec['n_edges']:,} distinct={rec['distinct_rna_cells_used']:,} "
                f"max={rec['reuse_max']}"
                + (f" [{rec['cost_form']} costs]" if rec["cost_form"] != "raw" else "")
                + f" | +{time.time()-t0:.0f}s")
    totals = {"n_blocks": len(levels), "n_edges_total": int(tot_edges),
              "n_columns_total": int(tot_slots),
              "n_blocks_not_solved_on_raw_costs": int(n_cond),
              "n_blocks_needing_k_escalation": int(n_esc),
              "solve_seconds": round(time.time() - t0, 1)}
    return asg, per_block, totals

def measure(asg, n, idn, donor, site, ct, chance_exact, n_blocks):
    cnt = np.bincount(asg, minlength=n)
    used = int((cnt > 0).sum())
    mult = cnt[cnt > 0]
    return {
        "n_assignments": int(n),
        "reuse_mean_per_used_rna_cell": float(n / used),
        "reuse_mean_per_assignment": float((mult.astype(float) ** 2).sum() / n),
        "reuse_max": int(cnt.max()),
        "reuse_distribution": {str(k): int(v) for k, v in
                               zip(*np.unique(mult, return_counts=True))},
        "share_of_assignments_on_a_reused_cell": float(mult[mult > 1].sum() / n),
        "distinct_rna_cells": used,
        "pct_pool_used": float(used / n),
        "self_hits": int((asg == idn).sum()),
        "exact_cell_accuracy": float((asg == idn).mean()),
        "random_within_block_expected_self_hits": int(n_blocks),
        "random_within_block_expected_exact_accuracy": float(chance_exact),
        "cross_donor_rate": float((donor != donor[asg]).mean()),
        "cross_site_rate": float((site != site[asg]).mean()),
        "cross_cell_type_rate": float((ct != ct[asg]).mean()),
    }

def load_space(a):
    """The prep_04b CCA space. --space-cache only skips the h5ad read and the SVD; the
    reproduction check below still runs on whatever came back, so a stale or corrupt cache is
    caught by the same per-cell check a fresh build is."""
    cache = Path(a.space_cache) if a.space_cache else None
    if cache is not None and cache.is_file():
        z = np.load(cache, allow_pickle=True)
        log(f"CCA space loaded from cache {cache.name} (the reproduction gate still runs on it)")
        return z["Pa"], z["Pr"], z["obs"], z["S"], z["shared"]
    Pa, Pr, obs, S, shared = build_cca_space(a.npc, a.nhvg, a.seed)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache, Pa=Pa, Pr=Pr, obs=obs, S=S, shared=shared)
        log(f"CCA space cached to {cache}")
    return Pa, Pr, obs, S, shared

def main() -> int:
    ap = argparse.ArgumentParser()

    ap.add_argument("--knn", type=int, default=100, help="forward kNN width, within block")
    ap.add_argument("--rev-knn", type=int, default=15, help="reverse kNN width, within block")
    ap.add_argument("--npc", type=int, default=50)
    ap.add_argument("--nhvg", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--cap", type=int, default=2, help="capacity given to the sampled fraction")
    ap.add_argument("--cap-frac", type=float, default=None,
                    help="fraction of each block's RNA cells given capacity --cap; the rest get 1. "
                         "1.0 is a uniform cap.")
    ap.add_argument("--sweep", default=None,
                    help="comma-separated cap fractions to MEASURE. Writes the sweep report and "
                         "NO arm table; the band gate does not apply to a sweep point.")
    ap.add_argument("--token", default="link_x_celltype_donor_cap14",
                    help="arm token; the output table is arm_<token>_seed<seed>.csv")
    ap.add_argument("--reuse-band", type=float, nargs=2, default=(1.35, 1.45),
                    metavar=("LO", "HI"),
                    help="hard gate on the realised reuse_mean_per_used_rna_cell")
    ap.add_argument("--outdir", default=str(OUT))
    ap.add_argument("--space-cache", default=None,
                    help="npz to load/save the CCA space; the reproduction gate runs either way")
    ap.add_argument("--solver-budget", type=float, default=90.0)
    a = ap.parse_args()
    if (a.cap_frac is None) == (a.sweep is None):
        raise SystemExit("supply exactly one of --cap-frac (build the arm) or --sweep (measure)")
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    md_raw = pd.read_csv(OUT / "b00_cell_metadata.csv", dtype=str).set_index("h5ad_barcode")
    obs_r = read_obs(OUT / "bmmc_rna_counts.h5ad")
    obs_a = read_obs(OUT / "bmmc_atac_gene_activity.h5ad")
    n = len(obs_r)
    if not np.array_equal(obs_r, obs_a):
        raise SystemExit("RNA and gene-activity cell order differ — the truth pairing would be wrong")
    absent = int(np.sum(~pd.Index(obs_r).isin(md_raw.index)))
    if absent:
        raise SystemExit(f"b00_cell_metadata.csv does not cover {absent} of the {n} cells")

    # ---- block squareness, derived from the two cell orders separately --------------------
    # Not needed for feasibility here (the expanded graph has C >= m columns by construction), but
    # it is the precondition for the claim that this arm sits on the same blocks as prep_04b's;
    # counting one modality's cell order twice would prove nothing.
    b_r = block_vector(md_raw, obs_r, BLOCK_COLS)
    b_a = block_vector(md_raw, obs_a, BLOCK_COLS)
    c_r = pd.Series(b_r).value_counts()
    c_a = pd.Series(b_a).value_counts()
    levels = sorted(set(c_r.index) | set(c_a.index))
    bad = [(x, int(c_a.get(x, 0)), int(c_r.get(x, 0))) for x in levels
           if int(c_a.get(x, 0)) != int(c_r.get(x, 0))]
    if bad:
        for x, na, nr in bad:
            log(f"  NON-SQUARE {BLOCK_DESC} block {x}: {na} ATAC vs {nr} RNA")
        raise SystemExit(f"{len(bad)} of {len(levels)} `{BLOCK_DESC}` blocks are not square — this "
                         f"arm would not sit on prep_04b's blocks. Stopping.")
    sizes = np.array([int(c_a[x]) for x in levels])
    if int(sizes.sum()) != n:
        raise SystemExit(f"block sizes sum to {sizes.sum():,}, not {n:,} — not a partition")
    chance = float(np.sum(sizes * (1.0 / sizes)) / n)
    log(f"squareness `{BLOCK_DESC}`: {len(levels)} blocks, all square, sizes {sizes.min()}-"
        f"{sizes.max()} (median {np.median(sizes):.0f}), sum {sizes.sum():,} = {n:,} | chance "
        f"exact-cell accuracy {chance:.4%} ({len(levels)} expected self-hits under a random "
        f"within-block draw)")

    Pa, Pr, obs_chk, S, shared = load_space(a)
    if not np.array_equal(obs_chk, obs_r):
        raise SystemExit("the CCA space's cell order differs from the metadata join order")
    gate = gate_reproduces_committed_link(Pa, Pr, obs_r, a.seed)

    md = md_raw.loc[obs_r].reset_index()
    donor = md["donor"].to_numpy(); site = md["site"].to_numpy()
    ct = md["cell_type"].to_numpy(); frag = md["fragment_barcode"].to_numpy()
    idn = np.arange(n)
    blk = b_r
    matcher = BudgetedMatcher(a.solver_budget)

    common = {
        "n_cells": int(n), "n_shared_genes": int(len(shared)), "n_hvg": int(a.nhvg),
        "n_pcs": int(a.npc), "top_singular_value": float(S[0]), "seed": a.seed,
        "block_definition": BLOCK_DESC, "block_columns": BLOCK_COLS,
        "n_blocks": int(len(levels)),
        "block_sizes": {"min": int(sizes.min()), "max": int(sizes.max()),
                        "median": float(np.median(sizes)), "sum": int(sizes.sum())},
        "random_within_block_expected_exact_accuracy": chance,
        "reproduction_gate": gate,
        "graph": {"forward_knn_requested": a.knn, "reverse_knn": a.rev_knn,
                  "solver_budget_s": a.solver_budget},
        "capacity_mechanism":
            "node duplication: RNA cell j is replicated into c_j columns carrying identical costs, "
            "and the arm is a minimum-cost FULL ROW matching of the m x C matrix computed by the "
            "same scipy min_weight_full_bipartite_matching prep_04b uses. The estimator is the "
            "same routine on a larger graph.",
        "capacity_scheme":
            f"c_j = {a.cap} for the first round(frac*m) cells of a fixed per-block permutation "
            f"seeded by (seed, block label), 1 for the rest; nested in frac, independent of the "
            f"CCA space",
        "method": "infer_stratified_pairings.py's within-block CCA linkage with a per-RNA-cell "
                  "capacity imposed by node duplication; see module docstring",
    }

    # ---------------------------------------------------------------- sweep mode ----
    if a.sweep is not None:
        fracs = [float(x) for x in a.sweep.split(",") if x.strip()]
        points = []
        for f in fracs:
            log(f"--- sweep point: cap={a.cap}, cap_frac={f} ---")
            asg, per_block, totals = build_arm(Pa, Pr, blk, levels, a.cap, f, a,
                                               matcher, f"cap{a.cap}f{f}")
            st = measure(asg, n, idn, donor, site, ct, chance, len(levels))
            pt = {"cap": int(a.cap), "cap_frac": f,
                  "n_cells_with_capacity_2": int(sum(
                      r["n_cells_with_capacity_2"] for r in per_block.values())),
                  "n_columns_total": totals["n_columns_total"],
                  "solve_seconds": totals["solve_seconds"], **st}
            points.append(pt)
            log(f"cap_frac={f}: reuse {st['reuse_mean_per_used_rna_cell']:.4f}x "
                f"(max {st['reuse_max']}) | distinct {st['distinct_rna_cells']:,} "
                f"({st['pct_pool_used']:.2%}) | exact-cell {st['exact_cell_accuracy']:.4%} | "
                f"cross-donor {st['cross_donor_rate']:.4f} cross-celltype "
                f"{st['cross_cell_type_rate']:.4f}")
        matcher.close()
        rep = dict(common)
        rep["sweep_points"] = points
        rep["target_band"] = list(a.reuse_band)
        rep["runtime_s"] = round(time.time() - t0, 1)
        dest = outdir / "b04d_capacity_sweep.json"
        dest.write_text(json.dumps(rep, indent=1), encoding="utf-8")
        log(f"wrote {dest.name} — NO arm table written (sweep mode)")
        return 0

    # ---------------------------------------------------------------- build mode ----
    log(f"--- {BLOCK_DESC}: {len(levels)} blocks, cap={a.cap}, cap_frac={a.cap_frac} ---")
    asg, per_block, totals = build_arm(Pa, Pr, blk, levels, a.cap, a.cap_frac, a,
                                       matcher, a.token)
    matcher.close()

    # ---- assert what the arm claims, on the realised assignment, not a proxy ----------
    if len(asg) != n:
        raise SystemExit(f"{a.token}: {len(asg)} assignments, expected {n}")
    if (blk[asg] != blk).any():
        raise SystemExit(f"{a.token} left its block for {int((blk[asg] != blk).sum()):,} cells")
    cnt_all = np.bincount(asg, minlength=n)
    cap_all = np.empty(n, dtype=np.int64)
    for lv in levels:
        ix = np.flatnonzero(blk == lv)
        cap_all[ix] = capacity_vector(str(lv), len(ix), a.cap, a.cap_frac, a.seed)
    if (cnt_all > cap_all).any():
        raise SystemExit(f"{a.token}: {int((cnt_all > cap_all).sum())} RNA cells used above their "
                         f"capacity — the cap did not take")

    st = measure(asg, n, idn, donor, site, ct, chance, len(levels))
    lo, hi = a.reuse_band
    if st["cross_donor_rate"] != 0.0:
        raise SystemExit(f"{a.token}: cross-donor rate {st['cross_donor_rate']} — must be 0")
    if st["cross_cell_type_rate"] != 0.0:
        raise SystemExit(f"{a.token}: cross-cell-type rate {st['cross_cell_type_rate']} — must be 0")
    if st["reuse_max"] > a.cap:
        raise SystemExit(f"{a.token}: reuse_max {st['reuse_max']} exceeds cap {a.cap}")
    if not (lo <= st["reuse_mean_per_used_rna_cell"] <= hi):
        raise SystemExit(
            f"{a.token}: realised reuse {st['reuse_mean_per_used_rna_cell']:.4f} is outside the "
            f"target band [{lo}, {hi}]. STOP — re-sweep --cap-frac; do not widen the band.")
    log(f"{a.token}: reuse {st['reuse_mean_per_used_rna_cell']:.4f}x (max {st['reuse_max']}) | "
        f"distinct {st['distinct_rna_cells']:,} ({st['pct_pool_used']:.2%}) | exact-cell "
        f"{st['exact_cell_accuracy']:.4%} = {st['self_hits']:,} self-hits (random-within-block "
        f"expectation {len(levels)}) | cross-donor {st['cross_donor_rate']:.4f} | cross-site "
        f"{st['cross_site_rate']:.4f} | cross-celltype {st['cross_cell_type_rate']:.4f}")

    # 14-column bundle assignment table, the arms.py schema. `subject` == `donor` in this
    # bundle, so cross_subject is the donor test, not a second axis.
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
    dest = outdir / f"arm_{a.token}_seed{a.seed}.csv"
    frame.to_csv(dest, index=False, lineterminator="\n")

    rep = dict(common)
    rep["cap"] = int(a.cap)
    rep["cap_frac"] = float(a.cap_frac)
    rep["n_cells_with_capacity_2"] = int(cap_all[cap_all > 1].size)
    rep["reuse_band_required"] = list(a.reuse_band)
    rep["graph"].update(totals)
    rep["graph"]["per_block"] = per_block
    rep["arm"] = {"token": a.token, "file": dest.name,
                  "row_unit": "one row per ATAC cell", **st}
    b04b = json.loads((OUT / "b04b_link_restricted_report.json").read_text(encoding="utf-8"))
    rep["reference_arms_from_b04b"] = {
        k: b04b["arms"][k] for k in ("link_x_celltype_donor", "link_x_celltype_donor_1to1")
        if k in b04b["arms"]}
    rep["runtime_s"] = round(time.time() - t0, 1)
    (outdir / "b04d_link_capped_report.json").write_text(json.dumps(rep, indent=1),
                                                         encoding="utf-8")
    log(f"wrote {dest.name}, b04d_link_capped_report.json -> {outdir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
