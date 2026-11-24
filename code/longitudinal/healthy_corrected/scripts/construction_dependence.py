#!/usr/bin/env python
"""How much does the healthy-cell activity readout depend on how the regulon activities were constructed (link vs repair)?"""
from __future__ import annotations

import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/healthy_corrected')
REPO = ROOT.parents[1]
LINK = REPO / "scoring/healthy"          # searched in out/ then logs/
REPAIR = REPO / "paml/healthy_corrected" / "out"
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[s3br-dep] {m}", flush=True)

def die(m):
    raise SystemExit(f"[s3br-dep] ABORT — {m}")

def find_link(name):
    for d in (LINK / "out", LINK / "logs"):
        p = d / name
        if p.exists() and p.stat().st_size > 0:
            return p
    die(f"missing link input {name} (looked in out/ and logs/)")

def read_auc(path):
    """Read the Gene_based modality only; Region_based is invalid in a cross-scored file (no signature matches)."""
    with h5py.File(path, "r") as f:
        g = f["mod"]["Gene_based"]
        oi, vi = g["obs"].attrs["_index"], g["var"].attrs["_index"]
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
        return (np.asarray(g["X"][:], dtype=np.float64),
                dec(g["obs"][oi][:]), np.array([regkey(v) for v in dec(g["var"][vi][:])]))

def main():
    rep = {"what": "healthy link vs repair activity agreement, per 2x2 cell",
           "caveat": "pairing AND vocabulary both change; intersecting on regulon name does "
                     "not neutralise the vocabulary change (shared name != shared target set)",
           "no_jaccard_because": "the healthy arm infers no network; it is scored under the "
                                 "malignant vocabulary, so there is no healthy topology",
           "cells": {}}

    for c in ("Diagnosis", "Relapse"):
        for d in ("Diagnosis", "Relapse"):
            tag = f"H{c}_defM{d}"
            Xl, ol, vl = read_auc(find_link(f"eRegulon_AUC_H{c}_link_defM{d}.h5mu"))
            Xr, orr, vr = read_auc(REPAIR / f"eRegulon_AUC_H{c}_repair_defM{d}.h5mu")

            cells = np.intersect1d(ol, orr)
            regs = np.intersect1d(vl, vr)
            if len(cells) == 0 or len(regs) == 0:
                die(f"{tag}: empty intersection (cells {len(cells)}, regulons {len(regs)})")

            li = {x: i for i, x in enumerate(ol)}
            ri = {x: i for i, x in enumerate(orr)}
            lv = {x: i for i, x in enumerate(vl)}
            rv = {x: i for i, x in enumerate(vr)}
            A = Xl[np.array([li[x] for x in cells])][:, np.array([lv[x] for x in regs])]
            B = Xr[np.array([ri[x] for x in cells])][:, np.array([rv[x] for x in regs])]
            del Xl, Xr

            # per-cell: mean activity across the shared regulons
            a_cell, b_cell = A.mean(axis=1), B.mean(axis=1)
            # per-regulon: mean activity across the shared cells
            a_reg, b_reg = A.mean(axis=0), B.mean(axis=0)

            rec = {
                "n_cells_link": int(len(ol)), "n_cells_repair": int(len(orr)),
                "n_cells_shared": int(len(cells)),
                "frac_of_repair_cells_shared": float(len(cells) / len(orr)),
                "n_regulons_link": int(len(vl)), "n_regulons_repair": int(len(vr)),
                "n_regulons_shared": int(len(regs)),
                "percell_pearson": float(np.corrcoef(a_cell, b_cell)[0, 1]),
                "percell_spearman": float(stats.spearmanr(a_cell, b_cell).statistic),
                "perregulon_pearson": float(np.corrcoef(a_reg, b_reg)[0, 1]),
                "perregulon_spearman": float(stats.spearmanr(a_reg, b_reg).statistic),
                "mean_activity_link": float(a_cell.mean()),
                "mean_activity_repair": float(b_cell.mean()),
            }
            rep["cells"][tag] = rec
            log(f"{tag}: cells {len(ol):,}/{len(orr):,} -> {len(cells):,} shared "
                f"({rec['frac_of_repair_cells_shared']:.1%} of repair) | "
                f"regulons {len(vl)}/{len(vr)} -> {len(regs)} shared")
            log(f"    per-cell   Pearson {rec['percell_pearson']:+.4f}  "
                f"Spearman {rec['percell_spearman']:+.4f}")
            log(f"    per-regulon Pearson {rec['perregulon_pearson']:+.4f}  "
                f"Spearman {rec['perregulon_spearman']:+.4f}")
            log(f"    mean activity  link {rec['mean_activity_link']:.5f}  "
                f"repair {rec['mean_activity_repair']:.5f}")
            del A, B

    (OUT / "s3br_construction_dependence.json").write_text(
        json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote {OUT / 's3br_construction_dependence.json'}")

if __name__ == "__main__":
    main()
