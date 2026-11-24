#!/usr/bin/env python
"""Does the cost conditioning in infer_stratified_pairings.py change the matching?"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import queue
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = Path(__file__).resolve().parent
OUT = ROOT / "benchmarks/bmmc" / "out"

import sys
sys.path.insert(0, str(SCRIPTS))
from infer_stratified_pairings import (  # noqa: E402
    BLOCK_DEFS, COST_FORMS, QUANT_GRID, block_vector, build_cca_space, log)

def _worker(tq, rq):
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching as _mwfbm
    while True:
        C = tq.get()
        if C is None:
            return
        try:
            _, c_i = _mwfbm(C)
            rq.put(("ok", c_i))
        except ValueError as e:
            rq.put(("infeasible", str(e)))
        except Exception as e:
            rq.put(("error", f"{type(e).__name__}: {e}"))

class Runner:
    def __init__(self, budget):
        self.budget = budget; self.ctx = mp.get_context("spawn"); self.p = None

    def solve(self, C):
        if self.p is None or not self.p.is_alive():
            self.tq = self.ctx.Queue(); self.rq = self.ctx.Queue()
            self.p = self.ctx.Process(target=_worker, args=(self.tq, self.rq), daemon=True)
            self.p.start()
        self.tq.put(C)
        try:
            kind, payload = self.rq.get(timeout=self.budget)   # drain before join, or deadlock
        except queue.Empty:
            self.p.terminate(); self.p.join(); self.p = None
            return None
        if kind != "ok":
            raise SystemExit(f"solver returned {kind}: {payload} — this script assumes the "
                             f"graphs prep_04b actually used, which were all feasible")
        return payload

    def close(self):
        if self.p is not None and self.p.is_alive():
            self.tq.put(None); self.p.join(timeout=10)
            if self.p.is_alive():
                self.p.terminate()

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--knn", type=int, default=100)
    ap.add_argument("--rev-knn", type=int, default=15)
    ap.add_argument("--npc", type=int, default=50)
    ap.add_argument("--nhvg", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--budget", type=float, default=30.0,
                    help="seconds one solve may take before it is recorded as stalled")
    a = ap.parse_args()
    t0 = time.time()

    Pa, Pr, obs, S, shared = build_cca_space(a.npc, a.nhvg, a.seed)
    md_raw = pd.read_csv(OUT / "b00_cell_metadata.csv", dtype=str).set_index("h5ad_barcode")
    runner = Runner(a.budget)

    rep = {"seed": a.seed, "knn": a.knn, "rev_knn": a.rev_knn, "budget_s": a.budget,
           "quant_grid": QUANT_GRID, "cost_forms": list(COST_FORMS),
           "statistic": "per-cell equality of the two matchings on the SAME graph",
           "families": {}}

    for tok, desc, cols in BLOCK_DEFS:
        blk = block_vector(md_raw, obs, cols)
        levels = sorted(set(blk.tolist()))
        fam = {"block_definition": desc, "n_blocks": len(levels),
               "rescaled": {"identical": 0, "different": [], "raw_stalled": [],
                            "conditioned_stalled": []},
               "quantized": {"identical": 0, "different": [], "raw_stalled": [],
                             "conditioned_stalled": []}}
        for lv in levels:
            ix = np.flatnonzero(blk == lv)
            m = len(ix)
            kf, kr = min(a.knn, m), min(a.rev_knn, m)
            nn_f = NearestNeighbors(n_neighbors=kf).fit(Pr[ix]); d_f, i_f = nn_f.kneighbors(Pa[ix])
            nn_r = NearestNeighbors(n_neighbors=kr).fit(Pa[ix]); d_r, i_r = nn_r.kneighbors(Pr[ix])
            rows = np.concatenate([np.repeat(np.arange(m), kf), i_r.ravel()])
            cols_ = np.concatenate([i_f.ravel(), np.repeat(np.arange(m), kr)])
            vals = np.concatenate([d_f.ravel(), d_r.ravel()]).astype(np.float64)
            key = rows.astype(np.int64) * m + cols_
            _, first = np.unique(key, return_index=True)
            rows, cols_, vals = rows[first], cols_[first], vals[first]
            lo, hi = float(vals.min()), float(vals.max()); rg = max(hi - lo, 1e-12)
            forms = {"raw": vals,
                     "rescaled": (vals - lo) / rg,
                     "quantized": np.round((vals - lo) / rg * QUANT_GRID)}
            sols = {}
            for nm, v in forms.items():
                sols[nm] = runner.solve(sp.csr_matrix((v + 1e-12, (rows, cols_)), shape=(m, m)))
            for nm in ("rescaled", "quantized"):
                d = fam[nm]
                if sols["raw"] is None and sols[nm] is None:
                    d["raw_stalled"].append(str(lv)); d["conditioned_stalled"].append(str(lv))
                elif sols["raw"] is None:
                    d["raw_stalled"].append(str(lv))
                elif sols[nm] is None:
                    d["conditioned_stalled"].append(str(lv))
                elif np.array_equal(sols["raw"], sols[nm]):
                    d["identical"] += 1
                else:
                    d["different"].append({"block": str(lv), "n_cells": int(m),
                                           "n_cells_differing": int((sols["raw"] != sols[nm]).sum())})
        for nm in ("rescaled", "quantized"):
            d = fam[nm]
            d["n_comparable"] = d["identical"] + len(d["different"])
            d["verdict"] = ("identical on every comparable block"
                            if not d["different"] else
                            f"DIFFERS on {len(d['different'])} block(s)")
            log(f"{desc:20s} {nm:10s}: {d['identical']}/{d['n_comparable']} comparable blocks "
                f"identical, {len(d['raw_stalled'])} not comparable (raw stalled), "
                f"{len(d['conditioned_stalled'])} conditioned stalled — {d['verdict']}")
        rep["families"][tok] = fam
        runner.close()

    rep["runtime_s"] = round(time.time() - t0, 1)
    bad = [f"{t}/{nm}" for t, f in rep["families"].items()
           for nm in ("rescaled", "quantized") if f[nm]["different"]]
    rep["overall_verdict"] = ("cost conditioning did not change any comparable matching"
                              if not bad else f"cost conditioning CHANGED the matching: {bad}")
    (OUT / "b04b_solver_equivalence.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(rep["overall_verdict"])
    log("wrote b04b_solver_equivalence.json")
    return 0 if not bad else 1

if __name__ == "__main__":
    raise SystemExit(main())
