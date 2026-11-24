#!/usr/bin/env python
"""Build the `link_1to1` arm: the link's own joint-CCA inference, constrained to a one-to-one matching.

Solves a minimum-cost bipartite matching between ATAC and RNA cells on joint-CCA distances.
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

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="10k")
    ap.add_argument("--knn", type=int, default=60,
                    help="candidate RNA cells per ATAC cell in the sparsified cost graph")
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--stage", default="pbmc10k")
    ap.add_argument("--rna-pool", default="analysis_set", choices=["analysis_set", "all"],
                    help="analysis_set (default): draw only from RNA cells that are rows of truth's "
                         "AUCell matrix. Required by p_08, whose exact gather looks each assigned "
                         "donor up in that matrix, so a donor outside it cannot be scored. On the "
                         "'all' pool 865 donors (7.92% of cells) were absent, leaving link_1to1 on "
                         "a different cell set from the other arms. Restricting also makes "
                         "link_1to1 a permutation of the truth pairing, like rand_subcluster_1to1, "
                         "so the two differ only in inference vs random assignment.")
    a = ap.parse_args()

    fr = ROOT / "work" / a.dataset / "archr_fullref"
    cca_fp = fr / f"jointcca_{a.dataset}.csv"
    lnk_fp = fr / f"link_fullref_{a.dataset}.csv"
    cd_fp = ROOT / "work" / f"cell_data_{a.stage}.csv"
    for fp in (cca_fp, lnk_fp, cd_fp):
        if not fp.is_file():
            print(f"[p17] MISSING {fp}", file=sys.stderr)
            return 2

    # ---- the shared space --------------------------------------------------
    J = pd.read_csv(cca_fp)
    J = J.set_index("cell")
    cc = [c for c in J.columns if c.startswith("CC_")]
    # Seurat's joint CCA suffixes each cell with its anchor role ("_reference" for RNA, "_query"
    # for ATAC); strip it before matching against the link table or the h5 barcode order.
    J.index = (J.index.str.replace(r"_(reference|query)$", "", regex=True))
    is_rna = J.index.str.startswith("RNACELL")
    is_atac = J.index.str.contains("#", regex=False)
    R = J.loc[is_rna, cc].to_numpy(dtype=np.float64)
    A = J.loc[is_atac, cc].to_numpy(dtype=np.float64)
    rna_ids = J.index[is_rna].to_numpy()
    atac_ids = J.index[is_atac].to_numpy()
    print(f"[p17] joint CCA: {A.shape[0]:,} ATAC x {R.shape[0]:,} RNA over {len(cc)} CC dims")
    if R.shape[0] < A.shape[0]:
        print(f"[p17] REFUSING: only {R.shape[0]:,} RNA cells for {A.shape[0]:,} ATAC cells -- "
              f"a bijection does not exist. Re-run p_16 with a full reference.", file=sys.stderr)
        return 3

    # RNACELL<n> -> true barcode, derived from the link table (which carries both) and checked
    # against every row.
    L = pd.read_csv(lnk_fp)
    m = L.dropna(subset=["predictedCell_internal", "Linked_scRNA_barcode"])
    rna2bc = dict(zip(m["predictedCell_internal"], m["Linked_scRNA_barcode"]))
    bad = [k for k, v in rna2bc.items() if not isinstance(v, str)]
    assert not bad, f"{len(bad)} RNACELL ids map to a non-string barcode"
    known = np.array([r in rna2bc for r in rna_ids])
    print(f"[p17] RNACELL->barcode map covers {known.sum():,} of {len(rna_ids):,} RNA cells "
          f"(only cells link actually chose are covered; the rest are resolved below)")

    # ---- restrict the RNA pool before matching, if asked -------------------
    # Done before matching, not by post-hoc filtering: a bijection solved on a larger pool and then
    # trimmed is neither optimal nor a bijection.
    if a.rna_pool == "analysis_set":
        auc_fp = ROOT / "out" / f"eRegulon_AUC_{a.stage}_truth.h5mu"
        if not auc_fp.is_file():
            print(f"[p17] MISSING {auc_fp} -- needed to define the analysis-set RNA pool",
                  file=sys.stderr)
            return 7
        import mudata
        md = mudata.read_h5mu(auc_fp)
        m0 = md.mod[list(md.mod.keys())[0]]
        allowed_bc = {str(x).split("___")[0] for x in m0.obs_names}
        # resolve RNACELL ids to barcodes first so the pool can be applied
        _pre = {}
        import h5py
        h5 = ROOT.parents[1] / "data/pbmc" / f"pbmc_granulocyte_sorted_{a.dataset}_filtered_feature_bc_matrix.h5"
        with h5py.File(h5, "r") as f:
            _bcs = np.array([b.decode() for b in f["matrix"]["barcodes"][:]])
        for rid in rna_ids:
            _pre[rid] = rna2bc.get(rid, _bcs[int(rid.replace("RNACELL", "")) - 1])
        keep = np.array([_pre[r] in allowed_bc for r in rna_ids])
        # The ATAC side must be restricted to the same set: native reference's AUCell matrix has 10,925 rows,
        # while the ArchR project carries 11,118 cells, so no bijection exists on the full set; the
        # 193 extra ATAC cells are not in the analysis set and no other arm scores them.
        keep_a = np.array([x.split("#")[-1] in allowed_bc for x in atac_ids])
        print(f"[p17] restricted to truth's analysis set:")
        print(f"[p17]   RNA  {keep.sum():,} of {len(rna_ids):,} kept "
              f"({(~keep).sum():,} outside the AUCell matrix)")
        print(f"[p17]   ATAC {keep_a.sum():,} of {len(atac_ids):,} kept "
              f"({(~keep_a).sum():,} not in the analysis set)")
        if keep.sum() < keep_a.sum():
            print(f"[p17] REFUSING: pool {keep.sum():,} < {keep_a.sum():,} ATAC cells; "
                  f"no bijection exists", file=sys.stderr)
            return 8
        R = R[keep]; rna_ids = rna_ids[keep]
        A = A[keep_a]; atac_ids = atac_ids[keep_a]
        rna2bc.update({r: _pre[r] for r in rna_ids})

    # ---- sparse symmetric cost graph: forward kNN (ATAC -> RNA) plus reverse kNN (RNA -> ATAC) ----
    nn_f = NearestNeighbors(n_neighbors=min(a.knn, R.shape[0])).fit(R)
    dist_f, idx_f = nn_f.kneighbors(A)
    rows = np.repeat(np.arange(A.shape[0]), idx_f.shape[1])
    cols = idx_f.ravel()
    vals = dist_f.ravel().astype(np.float64)

    m_rev = max(5, a.knn // 20)
    nn_r = NearestNeighbors(n_neighbors=min(m_rev, A.shape[0])).fit(A)
    dist_r, idx_r = nn_r.kneighbors(R)          # for each RNA cell, its nearest ATAC cells
    rows = np.concatenate([rows, idx_r.ravel()])
    cols = np.concatenate([cols, np.repeat(np.arange(R.shape[0]), idx_r.shape[1])])
    vals = np.concatenate([vals, dist_r.ravel().astype(np.float64)])

    # drop duplicate (row, col) edges, keeping the first
    key = rows.astype(np.int64) * R.shape[0] + cols.astype(np.int64)
    _, first = np.unique(key, return_index=True)
    rows, cols, vals = rows[first], cols[first], vals[first]
    # min_weight_full_bipartite_matching treats EXPLICIT zeros as edges; shift so no cost is 0
    C = csr_matrix((vals + 1e-12, (rows, cols)), shape=(A.shape[0], R.shape[0]))
    reach = np.unique(cols).size
    print(f"[p17] symmetric graph: forward k={a.knn} + reverse m={m_rev} -> {len(vals):,} edges; "
          f"{reach:,} of {R.shape[0]:,} RNA cells reachable")
    if reach < A.shape[0]:
        print(f"[p17] REFUSING: only {reach:,} RNA cells reachable for {A.shape[0]:,} ATAC cells",
              file=sys.stderr)
        return 4
    try:
        r_i, c_i = min_weight_full_bipartite_matching(C)
    except ValueError as e:
        print(f"[p17] REFUSING: no perfect matching on the symmetric graph ({e}). Raise --knn.",
              file=sys.stderr)
        return 4
    dist, idx = dist_f, idx_f
    print(f"[p17] bijection solved: {len(r_i):,} pairs")

    assigned_d = np.asarray(C[r_i, c_i]).ravel() - 1e-12
    nearest_d = dist[:, 0][r_i]
    displacement = assigned_d - nearest_d

    order = np.argsort(r_i)
    r_i, c_i = r_i[order], c_i[order]
    asn = pd.DataFrame({
        "cistopic_id": [x.split("#")[-1] + f"___pbmc_{a.dataset}" for x in atac_ids[r_i]],
        "rna_internal": rna_ids[c_i],
    })
    # resolve any RNACELL not present in the link table by its numeric suffix against the h5 order
    unresolved = sorted(set(asn["rna_internal"]) - set(rna2bc))
    if unresolved:
        import h5py
        h5 = ROOT.parents[1] / "data/pbmc" / f"pbmc_granulocyte_sorted_{a.dataset}_filtered_feature_bc_matrix.h5"
        if not h5.is_file():
            print(f"[p17] {len(unresolved):,} RNACELL ids unresolved and {h5} missing", file=sys.stderr)
            return 5
        with h5py.File(h5, "r") as f:
            bcs = np.array([b.decode() for b in f["matrix"]["barcodes"][:]])
        # check that the h5 order reproduces every mapping the link table already knows
        chk = [(rid, rna2bc[rid], bcs[int(rid.replace("RNACELL", "")) - 1])
               for rid in list(rna2bc)[:2000]]
        wrong = [c for c in chk if c[1] != c[2]]
        assert not wrong, (f"h5 barcode order disagrees with the link table for {len(wrong)} ids, "
                           f"e.g. {wrong[:3]} -- refusing to guess")
        print(f"[p17] h5 order verified against {len(chk):,} known mappings; resolving "
              f"{len(unresolved):,} further ids")
        for rid in unresolved:
            rna2bc[rid] = bcs[int(rid.replace("RNACELL", "")) - 1]

    asn["Linked_scRNA_barcode"] = asn["rna_internal"].map(rna2bc)
    assert asn["Linked_scRNA_barcode"].notna().all(), "unmapped RNACELL id survived"

    # ---- restrict to the analysis set and attach strata --------------------
    cd = pd.read_csv(cd_fp, index_col=0, dtype=str)
    idcol = "cell_id" if "cell_id" in cd.columns else None
    cd = cd.reset_index().rename(columns={cd.index.name or "index": "cistopic_id"})
    if idcol and "cistopic_id" not in cd.columns:
        cd = cd.rename(columns={idcol: "cistopic_id"})
    strat_col = next((c for c in ("subcluster", "Classified_Celltype", "annotated_celltype")
                      if c in cd.columns), None)
    if strat_col is None:
        print(f"[p17] no stratum column in {cd_fp.name}: {list(cd.columns)[:12]}", file=sys.stderr)
        return 6

    out = asn.merge(cd[["cistopic_id", strat_col]], on="cistopic_id", how="inner")
    out = out.rename(columns={strat_col: "stratum"})
    out["true_barcode"] = out["cistopic_id"].str.split("___").str[0]
    out["link_is_correct"] = out["Linked_scRNA_barcode"] == out["true_barcode"]
    out = out[["cistopic_id", "Linked_scRNA_barcode", "true_barcode", "stratum", "link_is_correct"]]

    n_distinct = out["Linked_scRNA_barcode"].nunique()
    assert n_distinct == len(out), (f"NOT 1-to-1: {len(out):,} rows but {n_distinct:,} distinct "
                                    f"RNA cells -- the bijection was violated by the merge")

    dest = ROOT / "out" / f"arm_link_1to1_{a.dataset}_seed{a.seed}_full.csv"
    out.to_csv(dest, index=False)

    rep = {
        "arm": "link_1to1", "dataset": a.dataset, "seed": a.seed,
        "rna_pool": a.rna_pool,
        "graph": {"forward_knn": int(a.knn), "reverse_knn": int(m_rev), "n_edges": int(len(vals))},
        "n_cells": int(len(out)),
        "distinct_rna": int(n_distinct),
        "reuse_factor": float(len(out) / n_distinct),
        "accuracy_exact_cell": float(out["link_is_correct"].mean()),
        "displacement_vs_own_nearest": {
            "median": float(np.median(displacement)),
            "p90": float(np.percentile(displacement, 90)),
            "frac_got_own_nearest": float((displacement <= 1e-9).mean()),
        },
        "comparator_link": {"accuracy_exact_cell": 0.004227, "reuse_factor": 6.095395,
                            "distinct_rna_used": 1824},
    }
    (ROOT / "out" / f"link_1to1_report_{a.dataset}.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))
    print(f"[p17] wrote {dest}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
