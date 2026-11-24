"""Freeze the module definition (regulons with survives_4 in the nomination table) as a portable, pre-specified object."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'adult/corrected')
REPO = ROOT.parents[1]
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)
S2 = REPO / "paml/malignant" / "out"
S3 = REPO / "longitudinal/corrected" / "out"

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[build] {m}", flush=True)

fun = pd.read_csv(S3 / "nomination_funnel.csv", index_col=0)
e2 = pd.read_csv(S3 / "e2_paired_wilcoxon_reconciliation.csv", index_col=0)

T4 = list(fun.index[fun.survives_4])
T5 = list(fun.index[fun.survives_all5])
BACKGROUND = list(e2.index)                      # the 431 shared-name regulons; null-distribution pool
log(f"T4 = {len(T4)} regulons / {fun.loc[T4,'TF'].nunique()} TFs; "
    f"T5 = {len(T5)}; background = {len(BACKGROUND)}")

# Target gene sets, per vocabulary, as built in the discovery cohort

targets = {}
for stage in ("Diagnosis", "Relapse"):
    e = pd.read_csv(S2 / f"eRegulons_{stage}_repair.tsv", sep="\t",
                    usecols=["eRegulon_name", "Gene"])
    e["key"] = [regkey(v) for v in e.eRegulon_name]
    targets[stage] = {k: sorted(set(map(str, g.unique())))
                      for k, g in e.groupby("key").Gene}
    log(f"{stage} vocabulary: {len(targets[stage])} regulons with target sets")

# Orientation (sign of the discovery-cohort paired effect)

orient = {}
for stage in ("Diagnosis", "Relapse"):
    col = f"beta_{stage}"
    o = {}
    for r in BACKGROUND:
        b = e2.loc[r, col] if r in e2.index else np.nan
        if not np.isfinite(b) or b == 0:
            continue
        o[r] = int(np.sign(b))
    orient[stage] = o
    n_t4 = sum(1 for r in T4 if r in o)
    log(f"{stage} orientation: {len(o)} regulons signed | T4 covered {n_t4}/{len(T4)} | "
        f"T4 up {sum(o[r] > 0 for r in T4 if r in o)} / down "
        f"{sum(o[r] < 0 for r in T4 if r in o)}")

# Assemble + sanity checks

missing = {s: [r for r in T4 if r not in targets[s]] for s in targets}
for s, m in missing.items():
    if m:
        raise SystemExit(f"[build] ABORT — {len(m)} T4 regulons have no target set under {s}: {m[:5]}")
for s in orient:
    m = [r for r in T4 if r not in orient[s]]
    if m:
        raise SystemExit(f"[build] ABORT — {len(m)} T4 regulons unsigned under {s}: {m[:5]}")

mod = {
    "name": "T4_tier1_core",
    "provenance": "longitudinal/published/out/nomination_funnel.csv (survives_4); orientation from "
                  "e2_paired_wilcoxon_reconciliation.csv; target sets from "
                  "paml/malignant/out/eRegulons_{stage}_repair.tsv",
    "frozen": "2026-08-06, before any GSE327386 score was computed",
    "n_regulons": len(T4), "n_TFs": int(fun.loc[T4, "TF"].nunique()),
    "TFs": sorted(fun.loc[T4, "TF"].unique().tolist()),
    "regulons": T4,
    "T5_subset_reference_only": T5,
    "T5_note": "FAILED the internal strictest null (null C) under both vocabularies; carried for "
               "reference, NOT for nomination",
    "background_regulons": BACKGROUND,
    "orientation": orient,
    "targets": {s: {r: targets[s][r] for r in BACKGROUND if r in targets[s]} for s in targets},
}

# target-set sizes
mod["target_set_sizes"] = {
    s: {"T4_median": float(np.median([len(targets[s][r]) for r in T4])),
        "T4_min": int(min(len(targets[s][r]) for r in T4)),
        "T4_max": int(max(len(targets[s][r]) for r in T4))} for s in targets}
for s, v in mod["target_set_sizes"].items():
    log(f"{s} T4 target-set size: median {v['T4_median']:.0f} (range {v['T4_min']}-{v['T4_max']})")

with open(OUT / "module_definition.json", "w") as fh:
    json.dump(mod, fh, indent=1)
log(f"wrote out/module_definition.json — {len(T4)} regulons, "
    f"{sum(len(v) for v in mod['targets'].values())} background target sets across 2 vocabularies")
