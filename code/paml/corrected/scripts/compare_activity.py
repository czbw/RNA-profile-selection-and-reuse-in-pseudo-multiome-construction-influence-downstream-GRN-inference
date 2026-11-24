#!/usr/bin/env python
"""Test whether the activity signal of the frozen 26-regulon module survives the corrected pairing.

Scores each cell by the mean AUCell value of the module regulons under both pairings and compares
per-patient Diagnosis-to-Relapse change. Pre-specified thresholds are not scored here.
"""
from __future__ import annotations

import json
from pathlib import Path

import re

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
S2 = ROOT / "paml/malignant" / "out"
OUT = ROOT / "paml/corrected" / "out"
MODDEF = ROOT / "adult/published" / "out" / "module_definition.json"

T4 = set(json.loads(MODDEF.read_text(encoding="utf-8"))["regulons"])

def log(m): print(f"[a04] {m}", flush=True)

def read_auc(stage: str, arm: str):
    """Gene-based AUC matrix, cell ids, regulon names."""
    p = S2 / f"eRegulon_AUC_{stage}_{arm}.h5mu"
    with h5py.File(p, "r") as f:
        X = f["mod/Gene_based/X"][:]
        cells = [c.decode() if isinstance(c, bytes) else str(c) for c in f["obs/Cell"][:]]
        var = f["mod/Gene_based/var/_index"][:]
        regs = [v.decode() if isinstance(v, bytes) else str(v) for v in var]
    # AUCell names carry a signature-size suffix (e.g. `SPI1_direct_-/+_(41g)`) that the frozen
    # list lacks. Strip it, and assert it changed something so a format change fails loudly.
    bare = [re.sub(r"_\((\d+)[gr]\)$", "", r) for r in regs]
    if sum(a != b for a, b in zip(regs, bare)) == 0:
        raise SystemExit("no AUCell name carried a _(Ng) suffix — check the naming convention")
    return X, np.array(cells, dtype=object), np.array(bare, dtype=object)

def t4_score(X, regs, cells):
    idx = [i for i, r in enumerate(regs) if r in T4]
    if not idx:
        raise SystemExit("no T4 regulon present in this arm — cannot score")
    return pd.Series(X[:, idx].mean(axis=1), index=cells), len(idx)

def main():
    rep = {"t4_n_regulons_frozen": len(T4),
           "original_primary_LOPO_for_orientation": {
               "Diagnosis": {"mean": 0.006886194656647156, "frac_up": 0.8, "n": 20},
               "Relapse": {"mean": 0.004126727831817664, "frac_up": 0.8, "n": 20}},
           "note": "pre-specified activity thresholds not scored; see module docstring",
           "stages": {}, "transition": {}}

    per_patient = {}
    for stage in ("Diagnosis", "Relapse"):
        log(f"===== {stage} =====")
        meta = pd.read_csv(S2 / f"cell_data_{stage}.csv", low_memory=False)
        # AUCell cell ids are `<barcode>___<sample_id>`, stored as two columns in the metadata.
        # Barcode alone is not unique across samples, so build the full key and assert the join.
        meta["_key"] = meta["barcode"].astype(str) + "___" + meta["sample_id"].astype(str)
        pid_col = "Patient_ID"
        meta = meta.drop_duplicates("_key").set_index("_key")
        log(f"  metadata: {len(meta):,} cells | key <barcode>___<sample_id> | "
            f"{meta[pid_col].nunique()} patients")

        s = {}
        scores = {}
        for arm in ("link", "repair"):
            X, cells, regs = read_auc(stage, arm)
            sc, n_t4 = t4_score(X, regs, cells)
            scores[arm] = sc
            s[arm] = {"n_cells": int(len(cells)), "n_regulons": int(len(regs)),
                      "n_t4_regulons_present": int(n_t4),
                      "t4_score_mean": float(sc.mean()), "t4_score_sd": float(sc.std())}
            log(f"  {arm:7} {len(cells):7,} cells | {len(regs):4} regulons | "
                f"T4 present {n_t4:2}/26 | mean T4 score {sc.mean():.5f}")

        common = scores["link"].index.intersection(scores["repair"].index)
        a, b = scores["link"].loc[common], scores["repair"].loc[common]
        pear = float(np.corrcoef(a.values, b.values)[0, 1])
        spear = float(stats.spearmanr(a.values, b.values).statistic)
        s["intersection"] = {
            "n_cells": int(len(common)),
            "pearson_link_vs_repair": pear, "spearman_link_vs_repair": spear}
        log(f"  intersection {len(common):,} cells | per-cell T4 score agreement: "
            f"Pearson {pear:.4f}  Spearman {spear:.4f}")

        pid = meta.reindex(common)[pid_col]          # NaN check before any str coercion
        ok = pid.notna().values
        matched = float(ok.mean())
        if matched < 0.90:
            raise SystemExit(f"{stage}: only {matched:.1%} of AUC cells joined to metadata — "
                             f"the id convention has changed; refusing to aggregate")
        pp = pd.DataFrame({"patient": pid.values[ok].astype(str),
                           "link": a.values[ok], "repair": b.values[ok]}
                          ).groupby("patient").mean()
        if len(pp) < 5:
            raise SystemExit(f"{stage}: {len(pp)} patient groups — the join produced a "
                             f"degenerate grouping; refusing to report a transition")
        log(f"  metadata join: {matched:.1%} of cells matched")
        per_patient[stage] = pp
        s["n_patients"] = int(len(pp))
        s["per_patient_corr_pearson"] = float(np.corrcoef(pp["link"], pp["repair"])[0, 1]) \
            if len(pp) > 2 else float("nan")
        log(f"  {len(pp)} patients | per-patient T4 mean agreement r = "
            f"{s['per_patient_corr_pearson']:.4f}")
        rep["stages"][stage] = s

    # ---- the transition, within patient, under each arm ------------------------------
    log("===== diagnosis -> relapse, within patient =====")
    dx, rl = per_patient["Diagnosis"], per_patient["Relapse"]
    shared = dx.index.intersection(rl.index)
    log(f"  {len(shared)} patients with both stages")
    tr = {"n_paired_patients": int(len(shared))}
    for arm in ("link", "repair"):
        d = (rl.loc[shared, arm] - dx.loc[shared, arm]).values
        w = stats.wilcoxon(d) if len(d) >= 6 else None
        tr[arm] = {"mean_delta": float(d.mean()), "median_delta": float(np.median(d)),
                   "frac_up": float((d > 0).mean()),
                   "p_wilcoxon": float(w.pvalue) if w else None,
                   "direction": "UP" if np.median(d) > 0 else "DOWN"}
        log(f"  {arm:7} delta(Rel-Dx): mean {d.mean():+.5f} | median {np.median(d):+.5f} | "
            f"frac up {(d>0).mean():.2f} | "
            f"p {w.pvalue:.4g}" if w else f"  {arm}: n too small for Wilcoxon")
    same = tr["link"]["direction"] == tr["repair"]["direction"]
    tr["direction_preserved"] = bool(same)
    tr["delta_correlation"] = float(np.corrcoef(
        (rl.loc[shared, "link"] - dx.loc[shared, "link"]).values,
        (rl.loc[shared, "repair"] - dx.loc[shared, "repair"]).values)[0, 1]) \
        if len(shared) > 2 else float("nan")
    log(f"  DIRECTION PRESERVED: {same}  ({tr['link']['direction']} under link, "
        f"{tr['repair']['direction']} under repair)")
    log(f"  per-patient delta correlation between arms: r = {tr['delta_correlation']:.4f}")
    rep["transition"] = tr

    (OUT / "a04_t1_activity_report.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote {OUT / 'a04_t1_activity_report.json'}")

if __name__ == "__main__":
    main()
