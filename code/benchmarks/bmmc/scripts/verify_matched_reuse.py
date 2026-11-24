#!/usr/bin/env python
"""Verify the reuse-1.4 PAIR, and measure the mechanism variable behind it."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from activity_common import load_gene_based, strip_suffix  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"
TAG = "b04e"

N_EXPECTED = 69249
LINK = "link_x_celltype_donor_cap14"
RAND = "rand_celltype_donor_matched_cap14"

def log(m): print(f"[{TAG}] {m}", flush=True)

def ks_distance(a, b) -> float:
    """Two-sample KS between two multiplicity vectors. arms.py's function, imported by
    value rather than by reference only because this script must not depend on b02's argv."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    grid = np.union1d(a, b)
    ca = np.searchsorted(a, grid, side="right") / a.size
    cb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))

def load_arm(path: Path, pos: pd.Series) -> np.ndarray:
    "-> assignment vector indexed by native reference's row order: asg[row of ATAC cell] = row of the RNA"
    if not path.is_file():
        raise SystemExit(f"[{TAG}] MISSING {path}")
    df = pd.read_csv(path, dtype=str)
    n = len(pos)
    if len(df) != n:
        raise SystemExit(f"[{TAG}] {path.name}: {len(df):,} rows, expected {n:,}")
    ai = df["atac_barcode"].map(pos)
    ri = df["assigned_rna_barcode"].map(pos)
    if ai.isna().any() or ri.isna().any():
        raise SystemExit(f"[{TAG}] {path.name}: {int(ai.isna().sum())} ATAC / "
                         f"{int(ri.isna().sum())} assigned barcodes are not cells of this dataset")
    ai = ai.to_numpy(dtype=np.int64)
    ri = ri.to_numpy(dtype=np.int64)
    if not np.array_equal(np.sort(ai), np.arange(n)):
        raise SystemExit(f"[{TAG}] {path.name}: atac_barcode is not a permutation of the cell set "
                         f"— the arm does not cover every ATAC cell exactly once")
    asg = np.empty(n, dtype=np.int64)
    asg[ai] = ri
    return asg

def var_blocks(Y: np.ndarray, ct_codes: np.ndarray, n_ct: int):
    """(total, within-cell-type) summed regulon variance of an activity matrix.

    total       sum_r Var(Y[:, r]) over all cells
    within      sum_r (1/n) * sum_c SS_c(Y[:, r]) — the pooled within-cell-type sum of squares,
                divided by n (not by n - n_ct), so that it is on the same scale as `total` and the
                ratio of two arms' values is a ratio of the same quantity."""
    n = Y.shape[0]
    total = float(Y.var(axis=0).sum())
    sums = np.zeros((n_ct, Y.shape[1]))
    np.add.at(sums, ct_codes, Y)
    cnt = np.bincount(ct_codes, minlength=n_ct).astype(float)[:, None]
    means = sums / cnt
    resid = Y - means[ct_codes]
    within = float((resid * resid).sum() / n)
    return total, within

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--h5mu", default=str(OUT / "eRegulon_AUC_bmmc_truth.h5mu"))
    a = ap.parse_args()

    md = pd.read_csv(OUT / "b00_cell_metadata.csv", dtype=str)
    for c in ("h5ad_barcode", "cell_type", "donor", "site"):
        if c not in md.columns:
            raise SystemExit(f"[{TAG}] b00_cell_metadata.csv has no '{c}' column")
    if not md["h5ad_barcode"].is_unique:
        raise SystemExit(f"[{TAG}] h5ad_barcode is not unique — it cannot be the cell key")

    # ---- native reference's AUCell matrix defines the cell order everything else is put into ----------
    X, cells, regs, mods = load_gene_based(TAG, Path(a.h5mu))
    bare = np.array(strip_suffix(list(cells)), dtype=object)
    if len(set(bare.tolist())) != len(bare):
        raise SystemExit(f"[{TAG}] truth's stripped row ids are not unique")
    if set(bare.tolist()) != set(md["h5ad_barcode"].tolist()):
        raise SystemExit(f"[{TAG}] truth's row ids are not the metadata h5ad_barcode set")
    n = len(bare)
    if n != N_EXPECTED:
        raise SystemExit(f"[{TAG}] truth's matrix has {n:,} rows, expected {N_EXPECTED:,}")
    pos = pd.Series(np.arange(n), index=bare)
    mdx = md.set_index("h5ad_barcode").loc[bare].reset_index()
    ct = mdx["cell_type"].to_numpy()
    donor = mdx["donor"].to_numpy()
    site = mdx["site"].to_numpy()
    blk = (mdx["cell_type"] + "|" + mdx["donor"]).to_numpy()
    lv, blk_codes = np.unique(blk, return_inverse=True)
    n_blocks = int(len(lv))
    ct_lv, ct_codes = np.unique(ct, return_inverse=True)
    idn = np.arange(n)
    log(f"truth AUCell Gene_based {X.shape[0]:,} x {X.shape[1]} | {n_blocks} (cell type x donor) "
        f"blocks | {len(ct_lv)} cell types | modalities {mods}")

    tot_t, win_t = var_blocks(X, ct_codes, len(ct_lv))
    log(f"truth's own summed regulon variance: total {tot_t:.6f} | within-cell-type {win_t:.6f}")

    # ---- the arms ---------------------------------------------------------------------------
    wanted = [
        (LINK, f"arm_{LINK}_seed{a.seed}.csv", OUT, "pair", True),
        (RAND, f"arm_{RAND}_seed{a.seed}.csv", OUT, "pair", True),
        ("link_x_celltype_donor_1to1", f"arm_link_x_celltype_donor_1to1_seed{a.seed}.csv", OUT,
         "positive control (bijection: both ratios must be exactly 1)", False),
        ("rand_celltype_donor_1to1", f"arm_rand_celltype_donor_1to1_seed{a.seed}.csv",
         OUT / ".." / "out_matched", "positive control (bijection: both ratios must be exactly 1)",
         False),
        ("link_x_celltype_donor", f"arm_link_x_celltype_donor_seed{a.seed}.csv", OUT,
         "reference: the reuse-2.008 pair", False),
        ("rand_celltype_donor_matched", f"arm_rand_celltype_donor_matched_seed{a.seed}.csv",
         OUT / ".." / "out_matched", "reference: the reuse-2.008 pair", False),
    ]

    rep = {"script": "verify_matched_reuse.py", "seed": a.seed,
           "n_cells": int(n), "n_regulons": int(X.shape[1]),
           "n_blocks_celltype_x_donor": n_blocks,
           "truth_h5mu": str(Path(a.h5mu)),
           "truth_variance": {"total": tot_t, "within_cell_type": win_t},
           "variance_definition":
               "an arm's Gene_based activity is EXACTLY a row-gather from truth's AUCell matrix "
               "(activity_common.py); total = sum_r Var(Y[:,r]); within = sum_r (1/n) sum_c SS_c; "
               "both reported as a fraction of truth's own",
           "self_hit_expectation_note":
               "one expected fixed point per stratum under a uniform within-stratum draw, for any "
               "stratum size and any d_s; 190 blocks => 190. A gate for the random arm, a "
               "reference for the link arm.",
           "arms": {}, "missing": {}}

    asgs = {}
    for name, fname, d, role, required in wanted:
        p = (Path(d) / fname).resolve()
        if not p.is_file():
            if required:
                raise SystemExit(f"[{TAG}] MISSING required arm table {p}")
            rep["missing"][name] = str(p)
            log(f"{name:34s} table absent ({p}) — skipped, recorded")
            continue
        asg = load_arm(p, pos)
        asgs[name] = asg
        cnt = np.bincount(asg, minlength=n)
        used = int((cnt > 0).sum())
        mult = cnt[cnt > 0]
        Y = X[asg]
        tot_a, win_a = var_blocks(Y, ct_codes, len(ct_lv))
        e = {
            "role": role, "file": p.name, "dir": str(p.parent.name),
            "n_assignments": int(len(asg)),
            "reuse_mean_per_used_rna_cell": float(n / used),
            "reuse_mean_per_assignment": float((mult.astype(float) ** 2).sum() / n),
            "reuse_max": int(cnt.max()),
            "distinct_rna_cells": used,
            "pct_pool_used": float(used / n),
            "self_hits": int((asg == idn).sum()),
            "self_hits_expected_uniform_within_block": n_blocks,
            "exact_cell_accuracy": float((asg == idn).mean()),
            "random_within_block_expected_exact_accuracy": float(n_blocks / n),
            "cross_donor_rate": float((donor != donor[asg]).mean()),
            "cross_site_rate": float((site != site[asg]).mean()),
            "cross_cell_type_rate": float((ct != ct[asg]).mean()),
            "retained_variance_total_frac_of_truth": float(tot_a / tot_t),
            "retained_variance_within_cell_type_frac_of_truth": float(win_a / win_t),
        }
        rep["arms"][name] = e
        log(f"{name:34s} reuse {e['reuse_mean_per_used_rna_cell']:.4f}x (max {e['reuse_max']}) | "
            f"distinct {used:,} ({e['pct_pool_used']:.2%}) | self-hits {e['self_hits']:,} "
            f"(uniform expectation {n_blocks}) | x-donor {e['cross_donor_rate']:.4f} | x-ct "
            f"{e['cross_cell_type_rate']:.4f} | retained var within "
            f"{e['retained_variance_within_cell_type_frac_of_truth']:.4f} total "
            f"{e['retained_variance_total_frac_of_truth']:.4f}")

    # ---- the estimator's own positive control -------------------------------------------
    # A bijection inside (cell type x donor) permutes rows WITHIN a cell type, so both ratios are
    # exactly 1 by construction. Asserted, because a variance estimator that cannot reproduce that
    # is measuring something other than what it claims.
    for name in ("link_x_celltype_donor_1to1", "rand_celltype_donor_1to1"):
        if name not in rep["arms"]:
            continue
        e = rep["arms"][name]
        for k in ("retained_variance_total_frac_of_truth",
                  "retained_variance_within_cell_type_frac_of_truth"):
            if abs(e[k] - 1.0) > 1e-9:
                raise SystemExit(f"[{TAG}] POSITIVE CONTROL FAILED: {name} is a within-(cell type "
                                 f"x donor) bijection so {k} must be 1.0, got {e[k]:.10f}. The "
                                 f"variance estimator is wrong; the pair's numbers cannot be read.")
    log("positive control PASSED: both within-block bijections retain exactly 1.0 of truth's "
        "total and within-cell-type regulon variance")

    # ---- the hard checks on the pair ------------------------------------------------------
    for name in (LINK, RAND):
        e = rep["arms"][name]
        if e["n_assignments"] != N_EXPECTED:
            raise SystemExit(f"[{TAG}] {name}: {e['n_assignments']} assignments, want {N_EXPECTED}")
        if e["cross_donor_rate"] != 0.0:
            raise SystemExit(f"[{TAG}] {name}: cross-donor {e['cross_donor_rate']}, must be 0")
        if e["cross_cell_type_rate"] != 0.0:
            raise SystemExit(f"[{TAG}] {name}: cross-cell-type {e['cross_cell_type_rate']}, "
                             f"must be 0")
        if not (1.35 <= e["reuse_mean_per_used_rna_cell"] <= 1.45):
            raise SystemExit(f"[{TAG}] {name}: reuse {e['reuse_mean_per_used_rna_cell']:.4f} is "
                             f"outside the [1.35, 1.45] band this pair is defined by")
    # the random arm's self-hit count is a sampler check: fixed points of a random permutation are
    # asymptotically Poisson(1) per stratum, so 5 sigma on 190 plus a 5-count slack
    band = 5.0 * float(np.sqrt(n_blocks)) + 5.0
    sh = rep["arms"][RAND]["self_hits"]
    if abs(sh - n_blocks) > band:
        raise SystemExit(f"[{TAG}] {RAND}: {sh} self-hits against an analytic {n_blocks} "
                         f"(band +/-{band:.1f}) — sampler bug, not noise")

    # ---- the pair: concentration and picks -----------------------------------------------
    def mult_vec(asg):
        c = np.bincount(asg, minlength=n)
        return c[c > 0]

    ks = ks_distance(mult_vec(asgs[LINK]), mult_vec(asgs[RAND]))
    used_l = set(np.unique(asgs[LINK]).tolist())
    used_r = set(np.unique(asgs[RAND]).tolist())
    inter = len(used_l & used_r)
    # analytic expectation: in block b both arms use d_b distinct cells (matched), and the random
    # arm's d_b are a uniform without-replacement draw from the block's m_b, so E|intersection|
    # in block b is d_b^2 / m_b. Summed over blocks.
    exp_inter = 0.0
    per_block_d = {}
    for b in range(n_blocks):
        ix = np.flatnonzero(blk_codes == b)
        m_b = int(ix.size)
        d_l = int(np.unique(asgs[LINK][ix]).size)
        d_r = int(np.unique(asgs[RAND][ix]).size)
        if d_l != d_r:
            raise SystemExit(f"[{TAG}] block {lv[b]}: link used {d_l} distinct RNA cells, random "
                             f"used {d_r} — the arms are NOT concentration-matched per stratum")
        per_block_d[str(lv[b])] = {"n_cells": m_b, "d_s": d_l}
        exp_inter += d_l * d_r / m_b
    rep["pair"] = {
        "arms": [LINK, RAND],
        "ks_multiplicity_distributions": ks,
        "ks_target": 0.0,
        "per_stratum_d_s_identical": True,
        "used_rna_cells_link": len(used_l),
        "used_rna_cells_random": len(used_r),
        "used_set_overlap": int(inter),
        "used_set_overlap_expected_at_chance": float(exp_inter),
        "used_set_overlap_ratio_to_chance": float(inter / exp_inter),
        "used_set_jaccard": float(inter / len(used_l | used_r)),
        "overlap_note": "expectation is sum_b d_b^2 / m_b — both arms use the same d_b in block b "
                        "and the random arm's picks are a uniform draw from the block's m_b cells",
    }
    if ks != 0.0:
        raise SystemExit(f"[{TAG}] the pair's multiplicity distributions differ (KS={ks:.6f}); the "
                         f"contrast would not be concentration-matched")
    log(f"PAIR: KS between multiplicity distributions {ks:.6f} (target 0.0) | used-set overlap "
        f"{inter:,} of {len(used_l):,} vs chance {exp_inter:,.1f} "
        f"(ratio {inter/exp_inter:.4f}, Jaccard {inter/len(used_l | used_r):.4f})")

    if LINK in rep["arms"] and RAND in rep["arms"]:
        dl = rep["arms"][LINK]["retained_variance_within_cell_type_frac_of_truth"]
        dr = rep["arms"][RAND]["retained_variance_within_cell_type_frac_of_truth"]
        rep["pair"]["within_cell_type_variance_gap_random_minus_link"] = float(dr - dl)
        log(f"MECHANISM: within-cell-type regulon variance retained — link {dl:.4f}, random "
            f"{dr:.4f}, gap {dr-dl:+.4f}")

    dest = OUT / "b04e_cap14_pair_verify.json"
    dest.write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote {dest.name}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
