"""Provenance checks and pseudobulk build for the cross-scored 2x2 (cell stage x regulon-definition stage)."""
import json
import os
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/corrected')
REPO = ROOT.parents[1]
S2 = REPO / "paml/malignant" / "out"
S3 = REPO / "scoring/published" / "out"

# Optional environment switches: cell filter, pinned unit list, output-directory suffix.
CELL_FILTER = os.environ.get("AMLGRN_S3_CELL_FILTER", "none")
if CELL_FILTER not in ("none", "own_timepoint"):
    raise SystemExit(f"[s3-prep] ABORT — AMLGRN_S3_CELL_FILTER={CELL_FILTER!r} is not one of "
                     "none | own_timepoint")
UNIT_SET = os.environ.get("AMLGRN_S3_UNIT_SET", "")
OUT = ROOT / ("out" + os.environ.get("AMLGRN_S3_OUT_SUFFIX", ""))
OUT.mkdir(parents=True, exist_ok=True)

MIN_CELLS_PER_UNIT = 20
DEAD_COL_EXCESS_MAX = 0.20
DEAD_COL_CONST_FRAC = 0.99

EXPECTED = {
    "Diagnosis": {"cells": 112_556, "defs": 332},
    "Relapse":   {"cells":  70_876, "defs": 352},
}

S3R = REPO / "scoring/corrected" / "out"

# (cells stage, defs stage) -> file
FILES = {
    ("Diagnosis", "Diagnosis"): S2 / "eRegulon_AUC_Diagnosis_repair.h5mu",
    ("Relapse",   "Relapse"):   S2 / "eRegulon_AUC_Relapse_repair.h5mu",
    ("Diagnosis", "Relapse"):   S3R / "eRegulon_AUC_Diagnosis_repair_defRelapse.h5mu",
    ("Relapse",   "Diagnosis"): S3R / "eRegulon_AUC_Relapse_repair_defDiagnosis.h5mu",
}

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[s3-prep] {m}", flush=True)

def die(m):
    raise SystemExit(f"[s3-prep] ABORT — {m}")

def read_auc(path):
    """Gene_based modality -> (X float64, obs_names, var_names). Region_based is never touched."""
    with h5py.File(path, "r") as f:
        g = f["mod"]["Gene_based"]
        oi, vi = g["obs"].attrs["_index"], g["var"].attrs["_index"]
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
        return np.asarray(g["X"][:], dtype=np.float64), dec(g["obs"][oi][:]), dec(g["var"][vi][:])

def dead_frac(X):
    """Dead frac."""
    n = X.shape[0]
    dead = 0
    for j in range(X.shape[1]):
        col = X[:, j]
        # modal-value share, computed cheaply: zeros dominate the degenerate case
        z = np.count_nonzero(col == 0.0)
        if z / n > DEAD_COL_CONST_FRAC or col.std() == 0.0:
            dead += 1
    return dead / X.shape[1]

prov = {"bars": {"min_cells_per_unit": MIN_CELLS_PER_UNIT,
                 "dead_col_excess_max": DEAD_COL_EXCESS_MAX,
                 "dead_col_const_frac": DEAD_COL_CONST_FRAC,
                 "expected_shapes": EXPECTED},
        "files": {}, "checks": {}}
# Recorded only when actually set, so a default run's s3_provenance.json is unchanged.
if CELL_FILTER != "none":
    prov["bars"]["cell_filter"] = CELL_FILTER
if UNIT_SET:
    prov["bars"]["unit_set"] = str(UNIT_SET)

# Load manifests

cell_data = {}
for s in ("Diagnosis", "Relapse"):
    d = pd.read_csv(S2 / f"cell_data_{s}.csv", index_col=0, dtype=str)
    cell_data[s] = d
    log(f"cell_data {s}: {len(d):,} cells | {d.Classified_Celltype.nunique()} cell types | "
        f"{d.Patient_ID.nunique()} patients | {d.Driver.nunique()} driver labels")

obs_by_cells = {}
var_by_defs = {}
dead_same = {}
pb_parts = {"Diagnosis": [], "Relapse": []}   # keyed by DEFS stage
n_parts = {"Diagnosis": [], "Relapse": []}

for (cells, defs), path in FILES.items():
    if not path.exists():
        die(f"missing input {path}")
    X, obs, var = read_auc(path)
    tag = f"{cells}_cells_x_{defs}_defs"
    log(f"{tag}: {X.shape[0]:,} x {X.shape[1]}  ({path.name})")

    if X.shape[0] != EXPECTED[cells]["cells"] or X.shape[1] != EXPECTED[defs]["defs"]:
        die(f"shape mismatch for {tag}: got {X.shape}, expected "
            f"({EXPECTED[cells]['cells']}, {EXPECTED[defs]['defs']})")

    if cells in obs_by_cells:
        ref = obs_by_cells[cells]
        if len(ref) != len(obs) or not np.array_equal(ref, obs):
            die(f"cell identity broken for {cells}: crossed file's obs differs from the "
                f"same-stage file's obs (content or order)")
    else:
        obs_by_cells[cells] = obs

    if defs in var_by_defs:
        ref = var_by_defs[defs]
        if len(ref) != len(var) or not np.array_equal(ref, var):
            die(f"vocabulary broken for defs={defs}: {tag} does not carry the same var index "
                f"as the {defs} same-stage file")
    else:
        var_by_defs[defs] = var

    df_dead = dead_frac(X)
    if cells == defs:
        dead_same[defs] = df_dead
    else:
        base = dead_same.get(defs)
        if base is None:
            die(f"dead-column check cannot evaluate {tag}: same-stage baseline for defs={defs} not loaded yet")
        if df_dead - base > DEAD_COL_EXCESS_MAX:
            die(f"dead-column rate for {tag} is {df_dead:.3f} vs same-stage baseline "
                f"{base:.3f} (excess {df_dead - base:.3f} > {DEAD_COL_EXCESS_MAX})")

    prov["files"][tag] = {
        "path": str(path), "bytes": path.stat().st_size,
        "n_cells": int(X.shape[0]), "n_regulons": int(X.shape[1]),
        "dead_col_frac": float(df_dead),
        "auc_min": float(X.min()), "auc_max": float(X.max()),
        "auc_mean": float(X.mean()), "frac_exact_zero": float(np.mean(X == 0.0)),
    }

    # ---- pseudobulk: mean activity per (patient x cell-type), this stage, this vocabulary ----
    cd = cell_data[cells].reindex(obs)
    keep = cd.Classified_Celltype.notna().values & cd.Patient_ID.notna().values
    if keep.sum() != len(obs):
        log(f"  {tag}: {len(obs) - keep.sum():,} cells dropped (no manifest row)")

    # own_timepoint: keep only cells whose linked RNA sample comes from the same time point
    if CELL_FILTER == "own_timepoint":
        own = (cd.sample_id.str.split("_").str[2] ==
               cd.Linked_scRNA_Sample.str.split("_").str[1]).values
        keep = keep & own
        log(f"  {tag}: cell filter own_timepoint -> {int(keep.sum()):,}/{len(obs):,} cells kept, "
            f"{int((~own).sum()):,} imported-label cells dropped")
        prov["files"][tag]["cells_kept_own_timepoint"] = int(keep.sum())

    cols = [regkey(v) for v in var]
    df = pd.DataFrame(X[keep], columns=cols)
    del X
    df["patient"] = cd.Patient_ID.values[keep]
    df["celltype"] = cd.Classified_Celltype.values[keep]
    grp = df.groupby(["patient", "celltype"], sort=True)
    mean = grp.mean()
    n = grp.size().rename("n_cells")
    mean.insert(0, "n_cells", n)
    mean.insert(0, "stage", cells)
    pb_parts[defs].append(mean.reset_index())
    del df, grp, mean

    # driver map, recorded once per stage
    prov["files"][tag]["n_units_all"] = int(len(pb_parts[defs][-1]))

# Write the two definition-set tables

driver_map = {}
for s in ("Diagnosis", "Relapse"):
    for pid, drv in cell_data[s].groupby("Patient_ID").Driver.agg(
            lambda v: v.dropna().unique()[0] if v.dropna().nunique() == 1 else "AMBIGUOUS").items():
        driver_map.setdefault(pid, drv)

for defs in ("Diagnosis", "Relapse"):
    tab = pd.concat(pb_parts[defs], ignore_index=True)
    tab.insert(3, "driver", tab.patient.map(driver_map).fillna("UNLABELLED"))
    tab.to_csv(OUT / f"pseudobulk_def{defs}.csv", index=False)
    log(f"pseudobulk_def{defs}.csv: {len(tab):,} units x "
        f"{tab.shape[1] - 5} regulons  (>= {MIN_CELLS_PER_UNIT} cells applied downstream)")
    prov["checks"][f"pseudobulk_def{defs}"] = {
        "n_units_total": int(len(tab)),
        "n_units_ge_floor": int((tab.n_cells >= MIN_CELLS_PER_UNIT).sum()),
        "n_regulons": int(tab.shape[1] - 5),
    }

tab = pd.read_csv(OUT / "pseudobulk_defDiagnosis.csv", usecols=
                  ["stage", "patient", "celltype", "driver", "n_cells"])
ok = tab[tab.n_cells >= MIN_CELLS_PER_UNIT]
dx = ok[ok.stage == "Diagnosis"].set_index(["patient", "celltype"]).index
rel = ok[ok.stage == "Relapse"].set_index(["patient", "celltype"]).index
paired_units = dx.intersection(rel)

if UNIT_SET:
    pinned = pd.MultiIndex.from_frame(
        pd.read_csv(UNIT_SET, dtype=str)[["patient", "celltype"]])
    log(f"design: PINNED to {len(pinned)} units from {UNIT_SET} "
        f"(re-derived here: {len(paired_units)}; "
        f"{len(pinned.difference(paired_units))} pinned units absent from it)")
    prov["checks"]["design_rederived_paired_units"] = int(len(paired_units))
    prov["checks"]["design_pinned_units_not_rederived"] = int(len(pinned.difference(paired_units)))
    paired_units = pinned

paired_patients = sorted({p for p, _ in paired_units})

pat_all_dx = set(tab[tab.stage == "Diagnosis"].patient)
pat_all_rel = set(tab[tab.stage == "Relapse"].patient)

ct_counts = pd.Series([c for _, c in paired_units]).value_counts()
prov["checks"]["design"] = {
    "patients_dx": len(pat_all_dx), "patients_rel": len(pat_all_rel),
    "patients_both": len(pat_all_dx & pat_all_rel),
    "paired_units_ge_floor": int(len(paired_units)),
    "paired_patients": len(paired_patients),
    "paired_celltypes": int(ct_counts.size),
    "celltype_patient_counts": {k: int(v) for k, v in ct_counts.items()},
    "driver_of_paired_patients": {
        k: int(v) for k, v in
        pd.Series([driver_map.get(p, "UNLABELLED") for p in paired_patients]).value_counts().items()},
}
log(f"design: {len(paired_units)} paired (patient x cell-type) units >= {MIN_CELLS_PER_UNIT} cells, "
    f"{len(paired_patients)} patients, {ct_counts.size} cell types")
log(f"design: cell types carrying the contrast — "
    + ", ".join(f"{k} ({v})" for k, v in ct_counts.head(8).items()))

pd.DataFrame({"patient": paired_patients,
              "driver": [driver_map.get(p, "UNLABELLED") for p in paired_patients]}
             ).to_csv(OUT / "paired_patients.csv", index=False)
pd.DataFrame(list(paired_units), columns=["patient", "celltype"]).to_csv(
    OUT / "paired_units.csv", index=False)

with open(OUT / "s3_provenance.json", "w") as fh:
    json.dump(prov, fh, indent=2)
log("Provenance checks passed. Wrote s3_provenance.json.")
log("NOTE: the dead-column check is a local proxy for the pre-flight check; s3_preflight_report.csv has not yet been produced.")
