#!/usr/bin/env python
"""Decide from the region-to-gene adjacencies whether an eGRN failure is a defect or a genuine result."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GRID = (0.01, 0.03, 0.05, 0.10)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--r2g", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--rho-threshold", type=float, default=0.05,
        help="threshold ASSUMED to mirror the scenicplus CLI default. The verdict is cross-checked "
             "against max/min rho, which are reported so it can be re-derived under any value.",
    )
    a = ap.parse_args()

    p = Path(a.r2g)
    if not p.is_file() or p.stat().st_size == 0:
        print(f"[egrn_guard] {p} missing or empty -- cannot confirm; treating failure as FATAL",
              file=sys.stderr)
        return 1

    rho = pd.read_csv(p, sep="\t", usecols=["rho"])["rho"].to_numpy(dtype=float)
    finite = rho[np.isfinite(rho)]
    if finite.size == 0:
        print("[egrn_guard] every rho is non-finite -- unexpected; treating failure as FATAL",
              file=sys.stderr)
        return 1

    t = a.rho_threshold
    n_pos = int((finite > t).sum())
    n_neg = int((finite < -t).sum())
    confirmed = (n_pos == 0) or (n_neg == 0)

    rec = {
        "tag": a.tag,
        "verdict": "eGRN_STRUCTURALLY_IMPOSSIBLE" if confirmed else "NOT_CONFIRMED",
        "reason": (
            "SCENIC+ create_emodules dichotomises region-to-gene on rho; an empty positive or "
            "negative set makes Groupby raise 'max() arg is an empty sequence'. The pairing "
            "destroys RNA-ATAC covariation, so no pair clears the threshold."
        ),
        "r2g_file": str(p),
        "n_pairs": int(rho.size),
        "n_rho_nonfinite": int(rho.size - finite.size),
        "max_rho": float(finite.max()),
        "min_rho": float(finite.min()),
        "mean_abs_rho": float(np.abs(finite).mean()),
        "sd_rho": float(finite.std(ddof=1)),
        "assumed_rho_threshold": t,
        "n_pos_at_assumed": n_pos,
        "n_neg_at_assumed": n_neg,
        # reported so the verdict can be re-derived under whatever default was actually in force
        "counts_by_threshold": {
            f"{g}": {"n_pos": int((finite > g).sum()), "n_neg": int((finite < -g).sum())}
            for g in GRID
        },
        # the noise-scale prediction this outcome was expected to match
        "expected_sd_if_pure_noise": float(1.0 / np.sqrt(10925)),
    }

    Path(a.out).write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2))

    if not confirmed:
        print(f"[egrn_guard] NOT CONFIRMED for {a.tag}: {n_pos} positive and {n_neg} negative pairs "
              f"clear rho={t}, so the empty-adjacency signature does not explain this failure. "
              f"Treating it as FATAL.", file=sys.stderr)
        return 1

    print(f"[egrn_guard] CONFIRMED for {a.tag}: max(rho)={finite.max():.4g}, "
          f"min(rho)={finite.min():.4g}; {n_pos} positive / {n_neg} negative pairs clear "
          f"rho={t}. No eGRN is constructible. Recorded as a RESULT, not a failure.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
