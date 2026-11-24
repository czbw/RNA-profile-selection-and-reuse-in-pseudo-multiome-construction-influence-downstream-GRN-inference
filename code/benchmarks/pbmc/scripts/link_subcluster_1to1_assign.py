#!/usr/bin/env python
"""Build the `link_subcluster_1to1` arm: the link's own CCA ranking, forced to a one-to-one matching within each subcluster block.

Same construction as p_17, but candidate edges are restricted to pairs sharing a subcluster label.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from sklearn.neighbors import NearestNeighbors

ROOT = (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')
ARM = "link_subcluster_1to1"
REFERENCE_ARM = "rand_subcluster_1to1"

def log(m):
    print(f"[p19] {m}", flush=True)

def solve_block(Ab, Rb, knn, dense=False):
    """Minimum-total-distance bijection inside one block.

    Returns (r_i, c_i, assigned_d, nearest_in_block_d, n_edges, complete).
    Edge set: symmetric forward(k) + reverse(m) kNN as in p_17, with k and m capped at the block
    size. `dense=True` uses the complete block graph instead.
    """
    n_a, n_r = Ab.shape[0], Rb.shape[0]
    assert n_a == n_r, f"block is not square: {n_a} ATAC vs {n_r} RNA"

    nn_f = NearestNeighbors(n_neighbors=n_r).fit(Rb)
    dist_all, idx_all = nn_f.kneighbors(Ab, n_neighbors=min(max(knn, 1), n_r))
    nearest_in_block = dist_all[:, 0].copy()

    if dense or knn >= n_r:
        # complete block graph
        from sklearn.metrics import pairwise_distances
        D = pairwise_distances(Ab, Rb)
        rows = np.repeat(np.arange(n_a), n_r)
        cols = np.tile(np.arange(n_r), n_a)
        vals = D.ravel().astype(np.float64)
        complete = True
    else:
        rows = np.repeat(np.arange(n_a), idx_all.shape[1])
        cols = idx_all.ravel()
        vals = dist_all.ravel().astype(np.float64)
        m_rev = min(max(5, knn // 20), n_a)
        nn_r = NearestNeighbors(n_neighbors=m_rev).fit(Ab)
        dist_r, idx_r = nn_r.kneighbors(Rb)
        rows = np.concatenate([rows, idx_r.ravel()])
        cols = np.concatenate([cols, np.repeat(np.arange(n_r), idx_r.shape[1])])
        vals = np.concatenate([vals, dist_r.ravel().astype(np.float64)])
        # de-duplicate: csr_matrix sums duplicate (i,j) entries, which would inflate the costs of
        # mutually nearest pairs (as in p_17)
        key = rows.astype(np.int64) * n_r + cols.astype(np.int64)
        _, first = np.unique(key, return_index=True)
        rows, cols, vals = rows[first], cols[first], vals[first]
        complete = False

    # min_weight_full_bipartite_matching treats EXPLICIT zeros as edges; shift so no cost is 0
    C = csr_matrix((vals + 1e-12, (rows, cols)), shape=(n_a, n_r))
    reach = np.unique(cols).size
    if reach < n_a:
        raise RuntimeError(f"only {reach} of {n_r} RNA cells reachable for {n_a} ATAC cells")
    r_i, c_i = min_weight_full_bipartite_matching(C)
    assigned_d = np.asarray(C[r_i, c_i]).ravel() - 1e-12
    return r_i, c_i, assigned_d, nearest_in_block, int(len(vals)), complete

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="10k")
    ap.add_argument("--knn", type=int, default=240,
                    help="candidate RNA cells per ATAC cell WITHIN a block; capped at block size. "
                         "240 is what the committed link_1to1 run used "
                         "(out/link_1to1_report_10k.json -> graph.forward_knn).")
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--stage", default="pbmc10k")
    ap.add_argument("--rna-pool", default="analysis_set", choices=["analysis_set", "all"],
                    help="verbatim from p_17: analysis_set restricts the donor pool to RNA cells "
                         "that are rows of truth's AUCell matrix, so the arm is a PERMUTATION of "
                         "the truth pairing exactly like rand_subcluster_1to1.")
    ap.add_argument("--dense-check", action="store_true",
                    help="re-solve every block on the COMPLETE block graph and report whether the "
                         "kNN sparsification changed any assignment. Diagnostic only; the shipped "
                         "tables always come from the kNN graph, matching p_17.")
    a = ap.parse_args()

    fr = ROOT / "work" / a.dataset / "archr_fullref"
    cca_fp = fr / f"jointcca_{a.dataset}.csv"
    lnk_fp = fr / f"link_fullref_{a.dataset}.csv"
    cd_fp = ROOT / "work" / f"cell_data_{a.stage}.csv"
    rep_fp = ROOT / "out" / f"random_arms_report_{a.dataset}_seed{a.seed}.json"
    ref_a_fp = ROOT / "out" / f"arm_{REFERENCE_ARM}_{a.dataset}_seed{a.seed}.csv"
    ref_b_fp = ROOT / "out" / f"arm_{REFERENCE_ARM}_{a.dataset}_seed{a.seed}_full.csv"
    for fp in (cca_fp, lnk_fp, cd_fp, rep_fp, ref_a_fp, ref_b_fp):
        if not fp.is_file():
            print(f"[p19] MISSING {fp}", file=sys.stderr)
            return 2

    # ---- the shared space (as in p_17) -------------------------------------
    J = pd.read_csv(cca_fp).set_index("cell")
    cc = [c for c in J.columns if c.startswith("CC_")]
    J.index = J.index.str.replace(r"_(reference|query)$", "", regex=True)
    is_rna = J.index.str.startswith("RNACELL")
    is_atac = J.index.str.contains("#", regex=False)
    R = J.loc[is_rna, cc].to_numpy(dtype=np.float64)
    A = J.loc[is_atac, cc].to_numpy(dtype=np.float64)
    rna_ids = J.index[is_rna].to_numpy()
    atac_ids = J.index[is_atac].to_numpy()
    log(f"joint CCA: {A.shape[0]:,} ATAC x {R.shape[0]:,} RNA over {len(cc)} CC dims")

    L = pd.read_csv(lnk_fp)
    m = L.dropna(subset=["predictedCell_internal", "Linked_scRNA_barcode"])
    rna2bc = dict(zip(m["predictedCell_internal"], m["Linked_scRNA_barcode"]))
    bad = [k for k, v in rna2bc.items() if not isinstance(v, str)]
    assert not bad, f"{len(bad)} RNACELL ids map to a non-string barcode"

    # resolve every RNACELL id to a barcode via the h5 order, checked against the link table
    import h5py
    h5 = ROOT.parents[1] / "data/pbmc" / f"pbmc_granulocyte_sorted_{a.dataset}_filtered_feature_bc_matrix.h5"
    if not h5.is_file():
        print(f"[p19] MISSING {h5}", file=sys.stderr)
        return 5
    with h5py.File(h5, "r") as f:
        bcs = np.array([b.decode() for b in f["matrix"]["barcodes"][:]])
    chk = [(rid, rna2bc[rid], bcs[int(rid.replace("RNACELL", "")) - 1]) for rid in list(rna2bc)[:2000]]
    wrong = [c for c in chk if c[1] != c[2]]
    assert not wrong, (f"h5 barcode order disagrees with the link table for {len(wrong)} ids, "
                       f"e.g. {wrong[:3]} -- refusing to guess")
    log(f"h5 barcode order verified against {len(chk):,} link-table mappings")
    rid2bc = {rid: rna2bc.get(rid, bcs[int(rid.replace("RNACELL", "")) - 1]) for rid in rna_ids}

    # ---- RNA pool restriction, before matching (as in p_17) ----------------
    if a.rna_pool == "analysis_set":
        auc_fp = ROOT / "out" / f"eRegulon_AUC_{a.stage}_truth.h5mu"
        if not auc_fp.is_file():
            print(f"[p19] MISSING {auc_fp} -- needed to define the analysis-set RNA pool",
                  file=sys.stderr)
            return 7
        import mudata
        md = mudata.read_h5mu(auc_fp)
        m0 = md.mod[list(md.mod.keys())[0]]
        allowed_bc = {str(x).split("___")[0] for x in m0.obs_names}
        keep = np.array([rid2bc[r] in allowed_bc for r in rna_ids])
        keep_a = np.array([x.split("#")[-1] in allowed_bc for x in atac_ids])
        log(f"analysis set: RNA {keep.sum():,}/{len(rna_ids):,} kept, "
            f"ATAC {keep_a.sum():,}/{len(atac_ids):,} kept")
        R, rna_ids = R[keep], rna_ids[keep]
        A, atac_ids = A[keep_a], atac_ids[keep_a]

    rna_bc = np.array([rid2bc[r] for r in rna_ids])
    atac_bc = np.array([x.split("#")[-1] for x in atac_ids])
    assert len(set(rna_bc)) == len(rna_bc), "duplicate barcodes on the RNA side of the CCA"
    assert len(set(atac_bc)) == len(atac_bc), "duplicate barcodes on the ATAC side of the CCA"

    # ---- the difference from p_17: subcluster labels on both sides ----------
    cd = pd.read_csv(cd_fp, dtype=str)
    if "cistopic_id" not in cd.columns and "cell_id" in cd.columns:
        cd = cd.rename(columns={"cell_id": "cistopic_id"})
    strat_col = next((c for c in ("subcluster", "Classified_Celltype", "annotated_celltype")
                      if c in cd.columns), None)
    if strat_col is None:
        print(f"[p19] no stratum column in {cd_fp.name}: {list(cd.columns)[:12]}", file=sys.stderr)
        return 6
    # In a true multiome the RNA and ATAC profiles of one barcode come from one cell, so the
    # subcluster of RNA barcode b is the Classified_Celltype of the ATAC cell with that barcode.
    # This matches p_07's donor pool: the true barcodes of the ATAC cells in the subcluster.
    bc2sub = dict(zip(cd["true_barcode"], cd[strat_col]))
    assert len(bc2sub) == len(cd), "true_barcode is not unique in the cell_data manifest"

    a_sub = np.array([bc2sub.get(b, "__UNLABELLED__") for b in atac_bc])
    r_sub = np.array([bc2sub.get(b, "__UNLABELLED__") for b in rna_bc])
    n_unlab_a, n_unlab_r = int((a_sub == "__UNLABELLED__").sum()), int((r_sub == "__UNLABELLED__").sum())
    if n_unlab_a or n_unlab_r:
        print(f"[p19] REFUSING: {n_unlab_a} ATAC and {n_unlab_r} RNA cells carry no subcluster "
              f"label; a within-subcluster restriction is undefined for them.", file=sys.stderr)
        return 9
    strata = sorted(set(a_sub) | set(r_sub))
    log(f"{len(strata)} subcluster blocks over {len(a_sub):,} ATAC / {len(r_sub):,} RNA cells")

    # ---- each block must be square (equal ATAC and RNA counts) --------------
    sq = []
    for s in strata:
        na, nr = int((a_sub == s).sum()), int((r_sub == s).sum())
        sq.append({"stratum": s, "n_atac": na, "n_rna": nr, "square": na == nr})
        log(f"   block {s:20s} ATAC {na:5,d}  RNA {nr:5,d}  {'square' if na == nr else 'NOT SQUARE'}")
    notsq = [b for b in sq if not b["square"]]
    if notsq:
        print(f"[p19] REFUSING: {len(notsq)} block(s) are not square, so no within-subcluster "
              f"bijection exists: {notsq}", file=sys.stderr)
        return 10

    # ---- global nearest-neighbour reference (displacement baseline, as in p_17) --
    nn_glob = NearestNeighbors(n_neighbors=1).fit(R)
    d_glob, _ = nn_glob.kneighbors(A)
    nearest_global = d_glob[:, 0]

    # ---- solve each block independently -------------------------------------
    rows_out, per_block, total_edges = [], {}, 0
    disp_glob, disp_blk = np.zeros(len(a_sub)), np.zeros(len(a_sub))
    dense_disagree = 0
    for s in strata:
        ia = np.flatnonzero(a_sub == s)
        ir = np.flatnonzero(r_sub == s)
        r_i, c_i, ad, nb, ne, complete = solve_block(A[ia], R[ir], a.knn)
        total_edges += ne
        order = np.argsort(r_i)
        r_i, c_i, ad = r_i[order], c_i[order], ad[order]
        # `ad` from solve_block is aligned to the pre-sort order; recompute in CCA space after the
        # sort so each distance belongs to the pair it is attached to.
        gi = ia[r_i]                                   # global ATAC index of each pair
        gj = ir[c_i]                                   # global RNA index of each pair
        adist = np.linalg.norm(A[gi] - R[gj], axis=1)
        # the cost the solver optimised must equal the CCA distance of the pair it returned
        assert np.max(np.abs(adist - ad)) < 1e-8, (
            f"block {s}: solver cost and recomputed CCA distance disagree by "
            f"{np.max(np.abs(adist - ad)):.3g}")
        disp_glob[gi] = adist - nearest_global[gi]
        disp_blk[gi] = adist - nb[r_i]
        correct = atac_bc[gi] == rna_bc[gj]
        n_s = len(ia)
        # greedy-within-block: is the true partner the nearest CCA neighbour inside the block?
        nn_b = NearestNeighbors(n_neighbors=1).fit(R[ir])
        _, nb_idx = nn_b.kneighbors(A[ia])
        greedy_correct = atac_bc[ia] == rna_bc[ir[nb_idx[:, 0]]]
        per_block[s] = {
            "n_s": int(n_s),
            "accuracy_exact_cell": float(correct.mean()),
            "n_correct": int(correct.sum()),
            # A uniform permutation inside a block of n hits the true partner with probability 1/n
            # per cell, so exactly one hit is expected per block for any n. Accuracy rates are not
            # comparable across block sizes, so the readout is the lift (observed hits / 1.0).
            "chance_accuracy_1_over_n": float(1.0 / n_s),
            "expected_correct_by_chance": 1.0,
            "lift_vs_chance": float(correct.sum()),
            "greedy_within_block_accuracy": float(greedy_correct.mean()),
            "n_edges": int(ne),
            "graph_complete": bool(complete),
            "median_displacement_vs_block_nearest": float(np.median(adist - nb[r_i])),
        }
        rows_out.append(pd.DataFrame({
            "cistopic_id": [b + f"___pbmc_{a.dataset}" for b in atac_bc[gi]],
            "Linked_scRNA_barcode": rna_bc[gj],
            "true_barcode": atac_bc[gi],
            "stratum": s,
        }))
        if a.dense_check and not complete:
            r2, c2, _, _, _, _ = solve_block(A[ia], R[ir], a.knn, dense=True)
            o2 = np.argsort(r2)
            same = np.array_equal(c2[o2], c_i)
            if not same:
                dense_disagree += int((c2[o2] != c_i).sum())
                dcorr = float((atac_bc[ia[r2[o2]]] == rna_bc[ir[c2[o2]]]).mean())
                per_block[s]["dense_check"] = {
                    "assignments_changed": int((c2[o2] != c_i).sum()),
                    "dense_accuracy_exact_cell": dcorr}
            else:
                per_block[s]["dense_check"] = {"assignments_changed": 0}

    out = pd.concat(rows_out, ignore_index=True)
    out["link_is_correct"] = out["Linked_scRNA_barcode"] == out["true_barcode"]
    out = out[["cistopic_id", "Linked_scRNA_barcode", "true_barcode", "stratum", "link_is_correct"]]

    # ---- invariants, checked on the output itself ---------------------------
    n_distinct = out["Linked_scRNA_barcode"].nunique()
    assert n_distinct == len(out), (f"NOT 1-to-1: {len(out):,} rows, {n_distinct:,} distinct donors")
    assigned_sub = out["Linked_scRNA_barcode"].map(bc2sub)
    assert assigned_sub.notna().all(), "an assigned barcode has no subcluster label"
    cross = float((assigned_sub.to_numpy() != out["stratum"].to_numpy()).mean())
    log(f"cross-subcluster rate = {cross:.10f}  (must be exactly 0)")
    if cross != 0.0:
        print(f"[p19] REFUSING: cross-subcluster rate {cross:.6%} != 0; the restriction did not "
              f"take and the arm is wrong.", file=sys.stderr)
        return 11

    # cell sets (full and floor-subset) must equal rand_subcluster_1to1's, in both directions
    ref_b = set(pd.read_csv(ref_b_fp, usecols=["cistopic_id"])["cistopic_id"])
    mine_b = set(out["cistopic_id"])
    if mine_b != ref_b:
        print(f"[p19] REFUSING: Track B cell set differs from {ref_b_fp.name}: "
              f"{len(mine_b - ref_b)} extra, {len(ref_b - mine_b)} missing", file=sys.stderr)
        return 12
    log(f"Track B cell-set assertion PASSED: {len(mine_b):,} cells, identical to {ref_b_fp.name}")

    rr = json.loads(rep_fp.read_text())
    floor = int(rr["min_stratum_cells_random_arm"])
    n_admit = int(rr["subclusters_admitted"])
    below = rr["subclusters_below_floor"]
    sizes = out["stratum"].value_counts()
    admitted = sorted(sizes[sizes >= floor].index)
    log(f"admission rule read from {rep_fp.name}: floor={floor} -> {len(admitted)} admitted "
        f"(report says {n_admit}); below floor here {dict(sizes[sizes < floor])} "
        f"vs report {below}")
    if len(admitted) != n_admit or {str(k): int(v) for k, v in sizes[sizes < floor].items()} != below:
        print(f"[p19] REFUSING: admission derived here disagrees with the committed report.",
              file=sys.stderr)
        return 13
    sub = out[out["stratum"].isin(admitted)].reset_index(drop=True)
    ref_a = set(pd.read_csv(ref_a_fp, usecols=["cistopic_id"])["cistopic_id"])
    mine_a = set(sub["cistopic_id"])
    if mine_a != ref_a:
        print(f"[p19] REFUSING: Track A cell set differs from {ref_a_fp.name}: "
              f"{len(mine_a - ref_a)} extra, {len(ref_a - mine_a)} missing", file=sys.stderr)
        return 14
    log(f"Track A cell-set assertion PASSED: {len(mine_a):,} cells, identical to {ref_a_fp.name}")
    assert sub["Linked_scRNA_barcode"].nunique() == len(sub), "not 1-to-1 after the floor subset"
    # blocks are independent, so dropping whole blocks cannot change any surviving assignment
    assert set(sub["stratum"]) == set(admitted), "the floor subset cut into an admitted block"

    dest_full = ROOT / "out" / f"arm_{ARM}_{a.dataset}_seed{a.seed}_full.csv"
    dest_a = ROOT / "out" / f"arm_{ARM}_{a.dataset}_seed{a.seed}.csv"
    out.to_csv(dest_full, index=False)
    sub.to_csv(dest_a, index=False)

    def disp(mask):
        d1, d2 = disp_glob[mask], disp_blk[mask]
        return {
            "median": float(np.median(d1)), "p90": float(np.percentile(d1, 90)),
            "frac_got_own_nearest": float((d1 <= 1e-9).mean()),
            "within_block_median": float(np.median(d2)),
            "within_block_p90": float(np.percentile(d2, 90)),
            "within_block_frac_got_own_nearest": float((d2 <= 1e-9).mean()),
        }

    admit_mask = np.isin(a_sub, admitted)

    # ---- chance baseline, stated in counts so that block size does not enter ------
    # Under rand_subcluster_1to1's construction the expectation is 1 hit per block whatever the
    # block size, so E[total] = number of blocks. The Poisson tail is used as the reference
    # distribution for the sum of independent per-block hit counts, each with mean 1.
    from scipy import stats as _st
    n_cor_full = int(out["link_is_correct"].sum())
    n_cor_a = int(sub["link_is_correct"].sum())
    chance = {
        "rule": "uniform permutation inside a block of n hits the true partner w.p. 1/n per cell "
                "=> exactly 1 expected hit per block, independent of n (as in p_07)",
        "full": {"n_blocks": len(strata), "expected_correct": float(len(strata)),
                 "observed_correct": n_cor_full,
                 "lift": float(n_cor_full / len(strata)),
                 "poisson_sf_p": float(_st.poisson.sf(n_cor_full - 1, len(strata))),
                 "observed_correct_rand_subcluster_1to1": int(
                     (pd.read_csv(ref_b_fp, dtype=str)["link_is_correct"] == "True").sum())},
        "track_a": {"n_blocks": len(admitted), "expected_correct": float(len(admitted)),
                    "observed_correct": n_cor_a,
                    "lift": float(n_cor_a / len(admitted)),
                    "poisson_sf_p": float(_st.poisson.sf(n_cor_a - 1, len(admitted))),
                    "observed_correct_rand_subcluster_1to1": int(
                        (pd.read_csv(ref_a_fp, dtype=str)["link_is_correct"] == "True").sum())},
    }
    _ns = np.array([per_block[s]["n_s"] for s in strata], float)
    _lift = np.array([per_block[s]["lift_vs_chance"] for s in strata], float)
    _sp = _st.spearmanr(_ns, _lift)
    size_lift = {
        "question": "does the CCA ranking have room to matter, and where? A ranked choice and a "
                    "random draw converge as n -> 1, so the lift is read against block size.",
        "spearman_rho_blocksize_vs_lift": float(_sp.statistic),
        "spearman_p": float(_sp.pvalue),
        "blocks_with_lift_ge_1": int((_lift >= 1).sum()), "n_blocks": len(strata),
    }
    if a.dense_check:
        _dn = sum((v["dense_check"]["dense_accuracy_exact_cell"] * v["n_s"]
                   if v.get("dense_check", {}).get("assignments_changed", 0) > 0
                   else v["accuracy_exact_cell"] * v["n_s"]) for v in per_block.values())
        _dna = sum((per_block[s]["dense_check"]["dense_accuracy_exact_cell"] * per_block[s]["n_s"]
                    if per_block[s].get("dense_check", {}).get("assignments_changed", 0) > 0
                    else per_block[s]["accuracy_exact_cell"] * per_block[s]["n_s"])
                   for s in admitted)
        chance["exact_within_block_optimum_sensitivity"] = {
            "_what": "the shipped tables use p_17's kNN-sparsified graph. This re-solves every "
                     "incomplete block on the COMPLETE block graph -- the exact minimum-cost "
                     "bijection the sparsification approximates -- and reports what changes.",
            "assignments_changed": int(dense_disagree),
            "full_n_correct": float(_dn), "full_accuracy_exact_cell": float(_dn / len(out)),
            "track_a_n_correct": float(_dna),
            "track_a_accuracy_exact_cell": float(_dna / len(sub)),
        }

    rep = {
        "arm": ARM, "dataset": a.dataset, "seed": a.seed, "rna_pool": a.rna_pool,
        "restriction": "candidate edges restricted to pairs sharing Classified_Celltype",
        "graph": {"forward_knn": int(a.knn), "reverse_knn": max(5, a.knn // 20),
                  "n_edges": int(total_edges), "per_block": True,
                  "blocks_with_complete_graph": int(sum(1 for v in per_block.values()
                                                        if v["graph_complete"]))},
        "cross_subcluster_rate": cross,
        "n_cells": int(len(out)),
        "distinct_rna": int(n_distinct),
        "reuse_factor": float(len(out) / n_distinct),
        "accuracy_exact_cell": float(out["link_is_correct"].mean()),
        "displacement_vs_own_nearest": disp(np.ones(len(a_sub), bool)),
        "track_a": {
            "min_stratum_cells_random_arm": floor,
            "n_cells": int(len(sub)),
            "n_strata": len(admitted),
            "distinct_rna": int(sub["Linked_scRNA_barcode"].nunique()),
            "reuse_factor": float(len(sub) / sub["Linked_scRNA_barcode"].nunique()),
            "accuracy_exact_cell": float(sub["link_is_correct"].mean()),
            "cross_subcluster_rate": float((sub["Linked_scRNA_barcode"].map(bc2sub).to_numpy()
                                            != sub["stratum"].to_numpy()).mean()),
            "displacement_vs_own_nearest": disp(admit_mask),
        },
        "block_squareness": sq,
        "per_stratum": per_block,
        "dense_check_total_assignments_changed": (int(dense_disagree) if a.dense_check
                                                  else "not run (--dense-check off)"),
        "chance_baseline_within_subcluster": chance,
        "size_vs_lift": size_lift,
        "comparators_measured_from_committed_files": {
            "_note": "cross_subcluster_rate recomputed here by mapping each assigned barcode to "
                     "its own Classified_Celltype in work/cell_data_pbmc10k.csv",
        },
    }
    for nm, fp in [("link_1to1_full", ROOT / "out" / f"arm_link_1to1_{a.dataset}_seed{a.seed}_full.csv"),
                   ("link_1to1_trackA", ROOT / "out" / f"arm_link_1to1_{a.dataset}_seed{a.seed}.csv"),
                   ("rand_subcluster_1to1_full", ref_b_fp),
                   ("rand_subcluster_1to1_trackA", ref_a_fp)]:
        if not fp.is_file():
            rep["comparators_measured_from_committed_files"][nm] = "unknown, no record found"
            continue
        d = pd.read_csv(fp, dtype=str)
        asub = d["Linked_scRNA_barcode"].map(bc2sub)
        rep["comparators_measured_from_committed_files"][nm] = {
            "file": fp.name, "n_cells": int(len(d)), "n_strata": int(d["stratum"].nunique()),
            "distinct_rna": int(d["Linked_scRNA_barcode"].nunique()),
            "accuracy_exact_cell": float((d["link_is_correct"] == "True").mean()),
            "cross_subcluster_rate": float((asub.to_numpy() != d["stratum"].to_numpy()).mean()),
        }

    rep_out = ROOT / "out" / f"{ARM}_report_{a.dataset}.json"
    rep_out.write_text(json.dumps(rep, indent=2))

    def _flat(s, v):
        d = {"stratum": s, **{k: x for k, x in v.items() if k != "dense_check"}}
        dc = v.get("dense_check")
        if v["graph_complete"]:
            # k >= n_s, so the block graph already IS the complete graph: shipped == exact.
            d["dense_assignments_changed"] = 0 if a.dense_check else None
            d["dense_accuracy_exact_cell"] = v["accuracy_exact_cell"] if a.dense_check else None
        else:
            d["dense_assignments_changed"] = dc.get("assignments_changed") if dc else None
            d["dense_accuracy_exact_cell"] = (dc.get("dense_accuracy_exact_cell")
                                              if dc and dc.get("assignments_changed", 0) > 0
                                              else (v["accuracy_exact_cell"] if dc else None))
        return d

    tab = pd.DataFrame([_flat(s, v) for s, v in per_block.items()])
    tab = tab.sort_values("n_s", ascending=False)
    tab["admitted_track_a"] = tab["stratum"].isin(admitted)
    tab_out = ROOT / "out" / f"{ARM}_per_stratum_{a.dataset}.csv"
    tab.to_csv(tab_out, index=False)

    log(f"wrote {dest_full.name}  ({len(out):,} cells, 18-stratum Track B set)")
    log(f"wrote {dest_a.name}  ({len(sub):,} cells, {len(admitted)}-stratum Track A set)")
    log(f"wrote {rep_out.name} and {tab_out.name}")
    print(json.dumps({k: v for k, v in rep.items()
                      if k not in ("per_stratum", "block_squareness")}, indent=2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
