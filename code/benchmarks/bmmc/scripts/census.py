#!/usr/bin/env python
"""BMMC feasibility census. Run first; later steps branch on it."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"
OUT.mkdir(parents=True, exist_ok=True)

# floors used for the power statements; declared here, not tuned after seeing the answer
FLOORS = (25, 50, 100, 200)

def log(m): print(f"[b00] {m}", flush=True)

def cat(f, col):
    """Read an AnnData categorical from the older __categories layout used by this file."""
    codes = f[f"obs/{col}"][:]
    cats = [x.decode() if isinstance(x, bytes) else str(x)
            for x in f[f"obs/__categories/{col}"][:]]
    return np.array([cats[c] for c in codes])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="decompressed multiome h5ad")
    a = ap.parse_args()

    f = h5py.File(a.h5ad, "r")
    idx = np.array([x.decode() for x in f["obs/_index"][:]])
    n = len(idx)
    batch = cat(f, "batch")
    ctype = cat(f, "cell_type")
    site = cat(f, "Site")
    donor = cat(f, "DonorNumber")

    # feature split — the GEX/ATAC boundary the export step depends on
    ft_codes = f["var/feature_types"][:]
    ft_cats = [x.decode() for x in f["var/__categories/feature_types"][:]]
    ft = np.array([ft_cats[c] for c in ft_codes])
    n_gex, n_atac = int((ft == "GEX").sum()), int((ft == "ATAC").sum())
    f.close()

    log(f"{n:,} cells | {n_gex:,} GEX features | {n_atac:,} ATAC peaks")
    log(f"{len(set(batch))} batches | {len(set(donor))} donors | {len(set(site))} sites "
        f"| {len(set(ctype))} cell types")

    # ---- the barcode remap, verified at 100% on batch s1d1 --------------
    # h5ad  AAACAGCCAATTAAGG-1-s1d1   ->   fragment  AAACAGCCAATTAAGG-1
    frag_bc = np.array([b.rsplit("-" + s, 1)[0] for b, s in zip(idx, batch)])
    bad = int(sum(1 for b, fb in zip(idx, frag_bc) if b == fb))
    if bad:
        raise SystemExit(f"barcode remap failed for {bad} cells — suffix not found")
    log(f"barcode remap OK for all {n:,} cells (strip '-<batch>')")

    # ---- stratum definitions --------------------------------------------------------
    # The contrast BMMC exists to provide is donor-blind vs donor-aware.
    defs = {
        "celltype":              [f"{c}" for c in ctype],
        "celltype_x_donor":      [f"{c}|{d}" for c, d in zip(ctype, donor)],
        "celltype_x_batch":      [f"{c}|{b}" for c, b in zip(ctype, batch)],
        "celltype_x_site":       [f"{c}|{s}" for c, s in zip(ctype, site)],
    }

    rows, summary = [], {}
    for name, keys in defs.items():
        cnt = Counter(keys)
        sizes = np.array(sorted(cnt.values()))
        summary[name] = {
            "n_strata": len(cnt),
            "median_size": float(np.median(sizes)),
            **{f"n_strata_ge_{fl}": int((sizes >= fl).sum()) for fl in FLOORS},
            **{f"pct_cells_in_strata_ge_{fl}":
               float(sizes[sizes >= fl].sum() / n) for fl in FLOORS},
        }
        for k, v in sorted(cnt.items()):
            rows.append({"definition": name, "stratum": k, "n_cells": v})
        log(f"  {name:20} {len(cnt):5,} strata | median {np.median(sizes):7.1f} | "
            f">=50: {(sizes>=50).sum():4,} | >=200: {(sizes>=200).sum():4,}")

    # ---- the cross-donor ceiling -----------------------------------------------------
    # A pairing built on cell-type label alone draws uniformly from that cell type. The
    # expected share of assignments landing on a DIFFERENT donor is 1 - sum_d p_d^2
    # within each cell type, averaged over cells. This is the quantity PBMC cannot have.
    xdonor, xsite = 0.0, 0.0
    for c in set(ctype):
        m = ctype == c
        w = m.sum() / n
        pd_ = np.array(list(Counter(donor[m]).values())) / m.sum()
        ps_ = np.array(list(Counter(site[m]).values())) / m.sum()
        xdonor += w * (1 - float((pd_ ** 2).sum()))
        xsite += w * (1 - float((ps_ ** 2).sum()))
    log(f"expected CROSS-DONOR assignment rate, cell-type-only pairing: {xdonor:.1%}")
    log(f"expected CROSS-SITE  assignment rate, cell-type-only pairing: {xsite:.1%}")

    # ---- write ------------------------------------------------------------------------
    import csv
    with open(OUT / "b00_cell_metadata.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["h5ad_barcode", "fragment_barcode", "batch", "donor", "site", "cell_type"])
        for i in range(n):
            w.writerow([idx[i], frag_bc[i], batch[i], donor[i], site[i], ctype[i]])
    with open(OUT / "b00_stratum_census.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["definition", "stratum", "n_cells"])
        w.writeheader(); w.writerows(rows)

    rep = {
        "dataset": "GSE194122 BMMC 10x Multiome (NeurIPS 2021)",
        "n_cells": n, "n_gex_features": n_gex, "n_atac_peaks": n_atac,
        "n_batches": len(set(batch)), "n_donors": len(set(donor)),
        "n_sites": len(set(site)), "n_cell_types": len(set(ctype)),
        "bijection_feasible": True,
        "bijection_note": "true multiome: rebuilding the pairing is a permutation within "
                          "stratum, so n_ATAC == n_RNA identically and there is no "
                          "scarcity floor (contrast: AML 1.39x)",
        "expected_cross_donor_rate_celltype_only": xdonor,
        "expected_cross_site_rate_celltype_only": xsite,
        "strata": summary,
    }
    (OUT / "b00_feasibility.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote b00_cell_metadata.csv, b00_stratum_census.csv, b00_feasibility.json -> {OUT}")

if __name__ == "__main__":
    main()
