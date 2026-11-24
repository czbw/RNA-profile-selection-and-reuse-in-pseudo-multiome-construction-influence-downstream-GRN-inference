#!/usr/bin/env python
'Figure 3c candidate: agreement of stratum-mean regulon activity between each BMMC arm and the native pairing.'
from __future__ import annotations

import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
BM = ROOT / "benchmarks/bmmc" / "out"
UPLOADS = ROOT / "benchmarks/bmmc" / "uploads"   # b_06 reads the arm tables from here
OUT = ROOT / "analysis" / "out" / "stratum_activity"

MIN_CELLS = 20      # a stratum enters when it holds at least this many cells
MIN_DONORS = 5      # a cell type enters when at least this many strata qualify
MIN_CELLS_CT = 20   # per-cell readout: same floor as b_06 (cell types with >= 20 cells)
SENSITIVITY = ((10, 5), (50, 5), (20, 3), (20, 7))

# (reuse label, role, arm token used by the manifest, table file)
ARMS = [
    ("1.000", "inferred", "link_x_celltype_donor_1to1", "arm_link_x_celltype_donor_1to1_seed555.csv"),
    ("1.000", "random", "rand_celltype_donor_1to1", "arm_rand_celltype_donor_1to1_seed555.csv"),
    ("1.406", "inferred", "link_x_celltype_donor_cap14", "arm_link_x_celltype_donor_cap14_seed555.csv"),
    ("1.406", "random", "rand_celltype_donor_matched_cap14", "arm_rand_celltype_donor_matched_cap14_seed555.csv"),
    ("2.008", "inferred", "link_x_celltype_donor", "arm_link_x_celltype_donor_seed555.csv"),
    ("2.008", "random", "rand_celltype_donor_matched", "arm_rand_celltype_donor_matched_seed555.csv"),
]
B06_KEY = {"link_x_celltype_donor_1to1": "link_x_celltype_donor_1to1",
           "rand_celltype_donor_1to1": "rand_celltype_x_donor_1to1",
           "link_x_celltype_donor_cap14": "link_x_celltype_donor_cap14",
           "rand_celltype_donor_matched_cap14": "rand_celltype_donor_matched_cap14",
           "link_x_celltype_donor": "link_x_celltype_donor",
           "rand_celltype_donor_matched": "rand_celltype_donor_matched"}

def log(m):
    print(f"[f3-01] {m}", flush=True)

def die(m):
    raise SystemExit(f"[f3-01] GATE FAILED: {m}")

def col_spearman(A, B):
    """Per-column Spearman (average ranks); NaN where a column is constant in either matrix."""
    ra = np.apply_along_axis(rankdata, 0, A)
    rb = np.apply_along_axis(rankdata, 0, B)
    ra -= ra.mean(0)
    rb -= rb.mean(0)
    den = np.sqrt((ra ** 2).sum(0) * (rb ** 2).sum(0))
    out = np.full(A.shape[1], np.nan)
    ok = (np.ptp(A, axis=0) > 0) & (np.ptp(B, axis=0) > 0) & (den > 0)
    out[ok] = (ra * rb).sum(0)[ok] / den[ok]
    return out

def per_cell_median(A, T, labels):
    ''
    order = sorted(c for c in np.unique(labels) if (labels == c).sum() >= MIN_CELLS_CT)
    W = []
    for c in order:
        m = labels == c
        row = np.full(A.shape[1], np.nan)
        for j in range(A.shape[1]):
            a, t = A[m, j], T[m, j]
            if np.ptp(a) and np.ptp(t):
                row[j] = spearmanr(a, t).statistic
        W.append(row)
    return float(np.nanmedian(np.vstack(W))), len(order)

def stratum_grid(Tm, Am, cts, min_donors):
    """Rows of the (cell type x regulon) grids over the strata in Tm/Am (already size-filtered)."""
    rho, nerr, kept = [], [], []
    for c in sorted(set(cts)):
        ix = np.flatnonzero(cts == c)
        if len(ix) < min_donors:
            continue
        T, A = Tm[ix], Am[ix]
        rho.append(col_spearman(A, T))
        sd = T.std(axis=0, ddof=1)
        rmse = np.sqrt(((A - T) ** 2).mean(axis=0))
        ne = np.full(T.shape[1], np.nan)
        ne[sd > 0] = rmse[sd > 0] / sd[sd > 0]
        nerr.append(ne)
        kept.append((c, len(ix)))
    return np.vstack(rho), np.vstack(nerr), kept

def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.load(open(ROOT / "figures" / "data" / "manifest.json", encoding="utf-8"))
    b06 = json.load(open(BM / "b06_trackA_L6_bmmc.json", encoding="utf-8"))["arms"]

    with h5py.File(BM / "eRegulon_AUC_bmmc_truth.h5mu", "r") as f:
        g = f["mod"]["Gene_based"]
        X = np.asarray(g["X"][:], dtype=np.float64)
        cells = np.char.replace(g["obs"]["Cell"][:].astype(str), "___bmmc", "")
        regs = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["var"]["_index"][:]])
    meta = pd.read_csv(BM / "b00_cell_metadata.csv").set_index("h5ad_barcode").loc[cells]
    ct = meta["cell_type"].to_numpy()
    donor = meta["donor"].to_numpy()
    stratum_label = np.array([f"{c}|{d}" for c, d in zip(ct, donor)])
    row_of = {b: i for i, b in enumerate(cells)}
    n = len(cells)
    log(f"activity matrix {X.shape}, {len(np.unique(stratum_label))} strata, {len(np.unique(ct))} cell types")

    strata, s_code = np.unique(stratum_label, return_inverse=True)
    s_n = np.bincount(s_code)
    s_ct = np.array([s.split("|")[0] for s in strata])

    def stratum_means(M, code):
        S = np.zeros((len(strata), M.shape[1]))
        np.add.at(S, code, M)
        return S / s_n[:, None]

    # G5: identity control
    Tm_all = stratum_means(X, s_code)
    keep = s_n >= MIN_CELLS
    rho_id, nerr_id, kept_id = stratum_grid(Tm_all[keep], Tm_all[keep], s_ct[keep], MIN_DONORS)
    if not (np.allclose(rho_id[~np.isnan(rho_id)], 1.0) and np.nanmax(np.abs(nerr_id)) == 0):
        die("identity control: truth against truth is not rho 1 / error 0")
    log(f"G5 identity control passed: {len(kept_id)} cell types, {int(keep.sum())} strata")

    rows, grid_rows = [], []
    for level, role, token, fname in ARMS:
        df = pd.read_csv(BM / fname if (BM / fname).is_file() else UPLOADS / fname, usecols=["true_rna_barcode", "assigned_rna_barcode"])
        own = np.array([row_of[b] for b in df["true_rna_barcode"]])
        got = np.array([row_of[b] for b in df["assigned_rna_barcode"]])
        if len(np.unique(own)) != n:
            die(f"{token}: own is not a permutation of all {n} cells")                                  # G3
        exact = float((own == got).mean())
        mkey = f"bmmc.exact.{token}" if f"bmmc.exact.{token}" in manifest else f"bmmc.exact.{B06_KEY[token]}"
        if abs(exact - manifest[mkey]["value"]) > 1e-12:
            die(f"{token}: exact-cell accuracy {exact} differs from manifest {mkey}")                    # G1
        if (s_code[own] != s_code[got]).any():
            die(f"{token}: cross-stratum assignments present")                                          # G2
        code = s_code[own]
        T = X[own]
        A = X[got]
        Tm, Am = stratum_means(T, code), stratum_means(A, code)

        if level == "1.000":                                                                           # G6
            if not np.allclose(Am, Tm, atol=1e-9, rtol=0):
                die(f"{token}: bijection within strata but stratum means differ from native "
                    f"(max {np.abs(Am - Tm).max():.3g})")

        rho_cells, n_ct_cells = per_cell_median(A, T, ct[own])
        committed = b06[B06_KEY[token]]["rho_within_median"]
        if abs(rho_cells - committed) > 5e-4:
            die(f"{token}: per-cell median {rho_cells:.4f} differs from committed {committed:.4f}")     # G4

        res = {"level": level, "role": role, "arm": token, "arm_series": B06_KEY[token], "table": fname,
               "exact_cell_accuracy": exact, "rho_cells": rho_cells,
               "rho_cells_committed": committed}
        for mc, md in ((MIN_CELLS, MIN_DONORS),) + SENSITIVITY:
            k = s_n >= mc
            rho, nerr, kept = stratum_grid(Tm[k], Am[k], s_ct[k], md)
            tag = "" if (mc, md) == (MIN_CELLS, MIN_DONORS) else f"_c{mc}_d{md}"
            res[f"rho_strata{tag}"] = float(np.nanmedian(rho))
            res[f"nerr_strata{tag}"] = float(np.nanmedian(nerr))
            res[f"n_celltypes{tag}"] = len(kept)
            res[f"n_strata{tag}"] = int(k[np.isin(s_ct, [c for c, _ in kept])].sum())
            if tag == "":
                res["rho_strata_q25"] = float(np.nanquantile(rho, .25))
                res["rho_strata_q75"] = float(np.nanquantile(rho, .75))
                res["nerr_strata_q25"] = float(np.nanquantile(nerr, .25))
                res["nerr_strata_q75"] = float(np.nanquantile(nerr, .75))
                res["n_grid_rho"] = int(np.isfinite(rho).sum())
                for (c, nd), r1, e1 in zip(kept, rho, nerr):
                    for j, rg in enumerate(regs):
                        grid_rows.append((level, role, token, c, nd, rg, r1[j], e1[j]))
        rows.append(res)
        log(f"{level} {role:8s} {token:36s} cells {rho_cells:.4f} | strata rho {res['rho_strata']:.4f} "
            f"| nerr {res['nerr_strata']:.4f} | {res['n_celltypes']} cell types")

    pd.DataFrame(grid_rows, columns=["level", "role", "arm", "cell_type", "n_donors", "regulon",
                                     "rho_strata", "nerr_strata"]).to_csv(OUT / "f3_stratum_mean_grid.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "f3_stratum_mean_summary.csv", index=False)
    json.dump({"thresholds": {"min_cells_per_stratum": MIN_CELLS, "min_strata_per_cell_type": MIN_DONORS},
               "sensitivity": [list(s) for s in SENSITIVITY], "arms": rows,
               "runtime_s": round(time.time() - t0, 1)},
              open(OUT / "f3_stratum_mean_summary.json", "w", encoding="utf-8"), indent=1)
    log(f"done in {time.time() - t0:.0f} s")

if __name__ == "__main__":
    main()
