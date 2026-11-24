"""Post-run verification (checks V1-V10) of the two cross-scored AUC files."""
import json
import os
from pathlib import Path
import re
import sys
from glob import glob

import h5py
import numpy as np
import pandas as pd

ARM = os.environ.get("AMLGRN_S3_PAIRING", "repair")

if len(sys.argv) > 1:                      # local mode: argv[1] is the repo root
    REPO = sys.argv[1].rstrip("/\\")
    S2 = os.path.join(REPO, "paml/malignant", "out")
    OUT = os.path.join(REPO, "scoring/corrected", "out")
else:                                      # CUBIC mode
    S2 = os.environ.get("AMLGRN_S2", str(Path(__file__).resolve().parents[3] / 'paml/malignant/out'))
    OUT = os.environ.get("AMLGRN_S3_OUT",
                         str(Path(__file__).resolve().parents[3] / 'scoring/corrected/out'))

STAGES = ["Diagnosis", "Relapse"]
OTHER = {"Diagnosis": "Relapse", "Relapse": "Diagnosis"}

# ---- bars ------------------------------------------------------------------------------
MIN_CELLS_PER_UNIT = 20        # V10: minimum cells per (patient x cell type) unit, as in compare_networks.py / prep.py
DEAD_COL_CONST_FRAC = 0.99     # V8: same as DEAD_COL_CONST_FRAC in prep.py
DEAD_COL_EXCESS_MAX = 0.20     # V8: same as DEAD_COL_EXCESS_MAX in prep.py
DEAD_COL_ABS_MAX = 0.20        # V8, absolute companion to the excess bar
FRAC_ZERO_MAX = 0.50           # V8
MEAN_AUC_LO, MEAN_AUC_HI = 0.05, 0.50   # V9
MEAN_RATIO_LO, MEAN_RATIO_HI = 0.5, 2.0  # V9
MIN_PATIENTS = 2               # V7

# V3: expected shapes, measured locally on 2026-08-18 from the two same-stage repair AUC files.
MEASURED = {"Diagnosis": {"cells": 112_556, "defs": 332},
            "Relapse":   {"cells":  70_876, "defs": 352}}

# (cells, defs) -> path. cells == defs are the Stage-2 same-stage files (the references);
# cells != defs are the two files this bundle produces.
FILES = {
    ("Diagnosis", "Diagnosis"): os.path.join(S2, f"eRegulon_AUC_Diagnosis_{ARM}.h5mu"),
    ("Relapse",   "Relapse"):   os.path.join(S2, f"eRegulon_AUC_Relapse_{ARM}.h5mu"),
    ("Diagnosis", "Relapse"):   os.path.join(OUT, f"eRegulon_AUC_Diagnosis_{ARM}_defRelapse.h5mu"),
    ("Relapse",   "Diagnosis"): os.path.join(OUT, f"eRegulon_AUC_Relapse_{ARM}_defDiagnosis.h5mu"),
}
CROSSED = [("Diagnosis", "Relapse"), ("Relapse", "Diagnosis")]

FAILURES = []
report = {"arm": ARM, "S2": S2, "OUT": OUT, "bars": {
    "min_cells_per_unit": MIN_CELLS_PER_UNIT, "dead_col_const_frac": DEAD_COL_CONST_FRAC,
    "dead_col_excess_max": DEAD_COL_EXCESS_MAX, "dead_col_abs_max": DEAD_COL_ABS_MAX,
    "frac_zero_max": FRAC_ZERO_MAX, "mean_auc_range": [MEAN_AUC_LO, MEAN_AUC_HI],
    "mean_ratio_range": [MEAN_RATIO_LO, MEAN_RATIO_HI], "min_patients": MIN_PATIENTS,
    "measured_shapes": MEASURED}, "files": {}, "checks": {}}

strip_sig = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[s3r-verify] {m}", flush=True)

def fail(bar, msg):
    FAILURES.append(f"{bar}  {msg}")
    log(f"FAIL  {bar}  {msg}")

def need(cond, bar, msg):
    if not cond:
        fail(bar, msg)
    return bool(cond)

def die_missing(msg):
    log(f"ABORT — {msg}")
    sys.exit(2)

def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])

def read_gene_based(path):
    """Gene_based modality -> (X float64, obs, var). Region_based is never touched."""
    with h5py.File(path, "r") as f:
        if "mod" not in f or "Gene_based" not in f["mod"]:
            return None, None, None
        g = f["mod"]["Gene_based"]
        oi, vi = g["obs"].attrs["_index"], g["var"].attrs["_index"]
        return (np.asarray(g["X"][:], dtype=np.float64),
                dec(g["obs"][oi][:]), dec(g["var"][vi][:]))

# V1 — presence

log(f"arm={ARM}")
log(f"S2  = {S2}")
log(f"OUT = {OUT}")
for (c, d), p in FILES.items():
    if not (os.path.exists(p) and os.path.getsize(p) > 0):
        if (c, d) in CROSSED:
            fail("V1", f"missing or empty output {p}")
        else:
            die_missing(f"reference file {p} not found — cannot evaluate any bar without it")

strays = sorted(glob(os.path.join(OUT, "*.part.*")))
need(not strays, "V1", f"stray .part sidecar(s) present (a killed run): {strays}")

for s in STAGES:
    for p in (os.path.join(S2, f"eRegulons_{s}_{ARM}.tsv"), os.path.join(S2, f"cell_data_{s}.csv")):
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            die_missing(f"{p} not found — needed for V5/V7/V10")

if FAILURES:                      # nothing downstream is meaningful without the two files
    log("")
    log("=" * 78)
    for f_ in FAILURES:
        log(f_)
    log("=" * 78)
    log("VERDICT: FAILED")
    sys.exit(5)

# V2/V3/V8/V9 — load each matrix once, reduce to scalars, keep only the indices

obs_of, var_of, stats = {}, {}, {}
for (c, d), p in FILES.items():
    tag = f"{c}_cells_x_{d}_defs"
    X, obs, var = read_gene_based(p)
    if X is None:
        fail("V2", f"{tag}: no Gene_based modality in {p}")
        continue
    if not need(X.ndim == 2, "V2", f"{tag}: Gene_based X is {X.ndim}-D, expected 2-D"):
        continue

    n, m = X.shape
    finite = np.isfinite(X)
    n_bad = int(finite.size - finite.sum())
    zero_frac = float(np.mean(X == 0.0))
    colzero = (X == 0.0).sum(axis=0) / n
    colstd = X.std(axis=0)
    dead_frac = float(np.mean((colzero > DEAD_COL_CONST_FRAC) | (colstd == 0.0)))
    st = {"path": p, "bytes": os.path.getsize(p), "n_cells": int(n), "n_regulons": int(m),
          "n_nonfinite": n_bad, "auc_min": float(np.nanmin(X)), "auc_max": float(np.nanmax(X)),
          "auc_mean": float(np.nanmean(X)), "frac_exact_zero": zero_frac,
          "dead_col_frac": dead_frac}
    del X, finite, colzero, colstd
    obs_of[(c, d)], var_of[(c, d)], stats[(c, d)] = obs, var, st
    report["files"][tag] = st
    log(f"{tag}: {n:,} x {m}  mean={st['auc_mean']:.5f}  dead={dead_frac:.4f}  "
        f"zero={zero_frac:.4f}  ({os.path.basename(p)})")

if len(stats) != 4:
    log("ABORT — could not read all four matrices; see V2 failures above")
    sys.exit(5)

# V3 — shape, derived from the same-stage files then cross-checked against measured constants
for s in STAGES:
    n_meas = stats[(s, s)]["n_cells"]
    m_meas = stats[(s, s)]["n_regulons"]
    need(n_meas == MEASURED[s]["cells"], "V3",
         f"same-stage {s}: {n_meas:,} cells, but {MEASURED[s]['cells']:,} was measured locally "
         f"on 2026-08-18 — the reference file changed under us")
    need(m_meas == MEASURED[s]["defs"], "V3",
         f"same-stage {s}: {m_meas} regulons, but {MEASURED[s]['defs']} was measured locally")

for (c, d) in CROSSED:
    exp = (stats[(c, c)]["n_cells"], stats[(d, d)]["n_regulons"])
    got = (stats[(c, d)]["n_cells"], stats[(c, d)]["n_regulons"])
    need(got == exp, "V3",
         f"{c} cells x {d} defs: shape {got} != expected {exp} "
         f"(= n_obs of the {c} same-stage file, n_var of the {d} same-stage file)")

# V4 — cell identity (content AND order)

for (c, d) in CROSSED:
    ref, got = obs_of[(c, c)], obs_of[(c, d)]
    same = len(ref) == len(got) and bool(np.array_equal(ref, got))
    need(same, "V4",
         f"{c} cells x {d} defs: obs differs from eRegulon_AUC_{c}_{ARM}.h5mu "
         f"(content or order). s3_10_prep would abort on this.")
    if same:
        log(f"OK    V4  {c} cells x {d} defs: obs identical to the same-stage file "
            f"({len(ref):,} ids, order preserved)")

# V5 — vocabulary, exact and suffix-aware

for (c, d) in CROSSED:
    ref, got = var_of[(d, d)], var_of[(c, d)]
    same = len(ref) == len(got) and bool(np.array_equal(ref, got))
    need(same, "V5",
         f"{c} cells x {d} defs: var differs from eRegulon_AUC_{d}_{ARM}.h5mu "
         f"(content or order). s3_10_prep would abort on this.")
    if same:
        log(f"OK    V5  {c} cells x {d} defs: var identical to the {d} same-stage file "
            f"({len(ref)} regulons)")

for s in STAGES:
    tsv = os.path.join(S2, f"eRegulons_{s}_{ARM}.tsv")
    t = pd.read_csv(tsv, sep="\t", usecols=["eRegulon_name", "Gene_signature_name"])
    names = set(t.eRegulon_name.astype(str))
    sigs = set(t.Gene_signature_name.astype(str))
    for (c, d) in CROSSED:
        if d != s:
            continue
        raw = set(map(str, var_of[(c, d)]))
        stripped = {strip_sig(v) for v in raw}
        need(stripped == names, "V5",
             f"{c} cells x {s} defs: suffix-stripped var != eRegulon_name in {os.path.basename(tsv)} "
             f"(var-only {len(stripped - names)}, tsv-only {len(names - stripped)})")
        need(raw == sigs, "V5",
             f"{c} cells x {s} defs: raw var != Gene_signature_name in {os.path.basename(tsv)} "
             f"(var-only {len(raw - sigs)}, tsv-only {len(sigs - raw)})")
        need(len(raw) == len(var_of[(c, d)]), "V5",
             f"{c} cells x {s} defs: duplicate regulon names in var "
             f"({len(var_of[(c, d)]) - len(raw)} duplicates)")
    report["checks"][f"defs_{s}"] = {"tsv_rows": int(len(t)), "n_eRegulon_name": len(names),
                                     "n_Gene_signature_name": len(sigs)}

# V6 / V7 — cell-id convention and the manifest join

manifest = {}
for s in STAGES:
    manifest[s] = pd.read_csv(os.path.join(S2, f"cell_data_{s}.csv"), index_col=0, dtype=str)
    log(f"cell_data_{s}.csv: {len(manifest[s]):,} rows, "
        f"{manifest[s].Patient_ID.nunique()} patients, "
        f"{manifest[s].Classified_Celltype.nunique()} cell types")

pat_count_same = {}
for (c, d), obs in obs_of.items():
    tag = f"{c}_cells_x_{d}_defs"
    obs_s = [str(o) for o in obs]

    # V6
    n_sep = sum(1 for o in obs_s if "___" in o)
    need(n_sep == len(obs_s), "V6",
         f"{tag}: {len(obs_s) - n_sep:,} of {len(obs_s):,} obs ids lack the '___' separator "
         f"(expected <barcode>___<sample_id>)")
    need(len(set(obs_s)) == len(obs_s), "V6",
         f"{tag}: {len(obs_s) - len(set(obs_s)):,} duplicate obs ids")

    # V7 — the join must be exact, and an all-NaN join must not pass
    idx = set(manifest[c].index.astype(str))
    hits = sum(1 for o in obs_s if o in idx)
    rate = hits / len(obs_s)
    need(rate == 1.0, "V7",
         f"{tag}: manifest join rate {rate:.6f} ({len(obs_s) - hits:,} of {len(obs_s):,} obs ids "
         f"absent from cell_data_{c}.csv). Bar is exactly 1.000000.")

    sub = manifest[c].reindex(pd.Index(obs_s))
    n_nan_p = int(sub.Patient_ID.isna().sum())
    n_nan_c = int(sub.Classified_Celltype.isna().sum())
    need(n_nan_p == 0, "V7", f"{tag}: {n_nan_p:,} cells with NaN Patient_ID after reindex")
    need(n_nan_c == 0, "V7", f"{tag}: {n_nan_c:,} cells with NaN Classified_Celltype after reindex")
    pats = set(sub.Patient_ID.dropna().astype(str))
    need("nan" not in {p.lower() for p in pats}, "V7",
         f"{tag}: a patient id is the literal string 'nan' — a failed join stringified into a "
         f"fake patient. This is the exact failure this bar exists for.")
    need(len(pats) >= MIN_PATIENTS, "V7",
         f"{tag}: only {len(pats)} distinct patient(s), bar is >= {MIN_PATIENTS}")
    if c == d:
        pat_count_same[c] = len(pats)
    else:
        need(len(pats) == pat_count_same.get(c), "V7",
             f"{tag}: {len(pats)} distinct patients but the {c} same-stage file has "
             f"{pat_count_same.get(c)}")
    report["files"][tag].update({"manifest_join_rate": rate, "n_patients": len(pats),
                                 "n_celltypes": int(sub.Classified_Celltype.nunique())})
    log(f"OK    V6/V7  {tag}: join {rate:.6f}, {len(pats)} patients, "
        f"{sub.Classified_Celltype.nunique()} cell types, 0 NaN")

# V8 — numeric sanity and dead columns

for (c, d) in CROSSED:
    st = stats[(c, d)]
    tag = f"{c}_cells_x_{d}_defs"
    need(st["n_nonfinite"] == 0, "V8", f"{tag}: {st['n_nonfinite']:,} non-finite (NaN/Inf) entries")
    need(st["auc_min"] >= 0.0, "V8", f"{tag}: min AUC {st['auc_min']:.6g} < 0")
    need(st["auc_max"] <= 1.0, "V8", f"{tag}: max AUC {st['auc_max']:.6g} > 1")
    need(st["frac_exact_zero"] <= FRAC_ZERO_MAX, "V8",
         f"{tag}: {st['frac_exact_zero']:.4f} of entries are exactly zero, bar <= {FRAC_ZERO_MAX}")
    need(st["dead_col_frac"] <= DEAD_COL_ABS_MAX, "V8",
         f"{tag}: dead-column fraction {st['dead_col_frac']:.4f}, bar <= {DEAD_COL_ABS_MAX}")
    base = stats[(d, d)]["dead_col_frac"]
    need(st["dead_col_frac"] - base <= DEAD_COL_EXCESS_MAX, "V8",
         f"{tag}: dead-column fraction {st['dead_col_frac']:.4f} vs same-vocabulary baseline "
         f"{base:.4f} (excess {st['dead_col_frac'] - base:.4f} > {DEAD_COL_EXCESS_MAX}). "
         f"s3_10_prep would abort on this.")

# V9 — sane mean AUC

for (c, d) in CROSSED:
    st = stats[(c, d)]
    tag = f"{c}_cells_x_{d}_defs"
    mu = st["auc_mean"]
    need(MEAN_AUC_LO <= mu <= MEAN_AUC_HI, "V9",
         f"{tag}: mean AUC {mu:.5f} outside [{MEAN_AUC_LO}, {MEAN_AUC_HI}]")
    ref = stats[(d, d)]["auc_mean"]
    ratio = mu / ref if ref > 0 else float("inf")
    need(MEAN_RATIO_LO <= ratio <= MEAN_RATIO_HI, "V9",
         f"{tag}: mean AUC {mu:.5f} is {ratio:.3f}x the same-vocabulary same-stage mean "
         f"{ref:.5f}; bar is [{MEAN_RATIO_LO}, {MEAN_RATIO_HI}]. On the link arm this ratio "
         f"was 0.998 and 0.970.")
    report["files"][tag]["mean_ratio_vs_same_vocabulary"] = ratio
    log(f"      V9  {tag}: mean {mu:.5f} = {ratio:.3f}x the {d}-vocabulary same-stage mean {ref:.5f}")

# V10 — pseudobulk feasibility (does a within-patient contrast survive the >= 20-cell floor?)

units = {}
for s in STAGES:
    obs_s = [str(o) for o in obs_of[(s, s)]]
    sub = manifest[s].reindex(pd.Index(obs_s))
    g = sub.groupby(["Patient_ID", "Classified_Celltype"]).size()
    keep = g[g >= MIN_CELLS_PER_UNIT]
    units[s] = set(keep.index)
    need(len(keep) >= 1, "V10",
         f"{s}: no (patient x cell-type) unit reaches {MIN_CELLS_PER_UNIT} cells")
    report["checks"][f"units_{s}"] = {"n_units_all": int(len(g)),
                                      "n_units_ge_floor": int(len(keep))}
    log(f"      V10 {s}: {len(g)} units, {len(keep)} at >= {MIN_CELLS_PER_UNIT} cells")

paired = units["Diagnosis"] & units["Relapse"]
paired_patients = sorted({p for p, _ in paired})
need(len(paired) >= 1, "V10",
     "no PAIRED (patient x cell-type) unit exists at the >= 20-cell floor — there is no "
     "within-patient Dx->Rel contrast to fit, so these files cannot serve Stage 3")
report["checks"]["design"] = {"paired_units_ge_floor": int(len(paired)),
                              "paired_patients": len(paired_patients)}
log(f"      V10 paired: {len(paired)} (patient x cell-type) units across "
    f"{len(paired_patients)} patients")

report["failures"] = FAILURES
report["verdict"] = "FAILED" if FAILURES else "ALL CHECKS PASSED"
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, f"s3r_verify_{ARM}.json"), "w") as fh:
    json.dump(report, fh, indent=2)

log("")
log("=" * 78)
if FAILURES:
    for f_ in FAILURES:
        log(f_)
    log("=" * 78)
    log(f"VERDICT: FAILED — {len(FAILURES)} bar(s) violated. Do NOT download these files and do "
        f"NOT run s3_10_prep against them.")
    sys.exit(5)
log("ALL CHECKS PASSED — V1..V10. The two cross-scored files are what s3_10_prep expects.")
log("Gene_based ONLY: never read Region_based out of a cross-scored file.")
log("=" * 78)
sys.exit(0)
