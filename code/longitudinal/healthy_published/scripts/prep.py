"""Provenance checks and pseudobulk construction for the healthy-compartment 2x2 (healthy cells x malignant-defined regulon vocabulary)."""
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/healthy_published')
REPO = ROOT.parents[1]
S3B = REPO / "scoring/healthy"
S2 = REPO / "paml/malignant" / "out"
S3 = REPO / "scoring/published" / "out"
MANIFEST = REPO / "preprocessing/out"
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)

MIN_CELLS_PER_UNIT = 20           
BARCODE_MATCH_MIN = 0.99          
DEAD_COL_CONST_FRAC = 0.99        
EXPECTED_DEFS = {"Diagnosis": 686, "Relapse": 793}   

# the 2x2: (healthy cells stage, malignant defs stage) -> file; inputs may be in out/ or logs/.
def find(name):
    for d in (S3B / "out", S3B / "logs"):
        p = d / name
        if p.exists() and p.stat().st_size > 0:
            return p
    raise SystemExit(f"[s3b-prep] ABORT — missing input {name} (looked in out/ and logs/)")

FILES = {(c, d): find(f"eRegulon_AUC_H{c}_link_defM{d}.h5mu")
         for c in ("Diagnosis", "Relapse") for d in ("Diagnosis", "Relapse")}

MALIG_REF = {"Diagnosis": S2 / "eRegulon_AUC_Diagnosis_link.h5mu",
             "Relapse":   S2 / "eRegulon_AUC_Relapse_link.h5mu"}

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[s3b-prep] {m}", flush=True)

def die(m):
    raise SystemExit(f"[s3b-prep] ABORT — {m}")

def read_auc(path, values=True):
    """Gene_based modality -> (X, obs_names, var_names). Region_based is never touched."""
    with h5py.File(path, "r") as f:
        g = f["mod"]["Gene_based"]
        oi, vi = g["obs"].attrs["_index"], g["var"].attrs["_index"]
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
        X = np.asarray(g["X"][:], dtype=np.float64) if values else None
        return X, dec(g["obs"][oi][:]), dec(g["var"][vi][:])

def dead_frac(X):
    n = X.shape[0]
    dead = 0
    for j in range(X.shape[1]):
        col = X[:, j]
        if np.count_nonzero(col == 0.0) / n > DEAD_COL_CONST_FRAC or col.std() == 0.0:
            dead += 1
    return dead / X.shape[1]

prov = {"bars": {"min_cells_per_unit": MIN_CELLS_PER_UNIT,
                 "barcode_match_min": BARCODE_MATCH_MIN,
                 "expected_defs": EXPECTED_DEFS},
        "files": {}, "checks": {}}

cell_data = {}
for s in ("Diagnosis", "Relapse"):
    d = pd.read_csv(MANIFEST / f"cell_data_{s}.csv", index_col=0, dtype=str)
    # rebuild the cisTopic-style index from the columns; never trust the file's own index
    d = d.assign(_cistopic=d.barcode.astype(str) + "___" + d.sample_id.astype(str))
    if d._cistopic.duplicated().any():
        die(f"manifest {s}: duplicate cisTopic barcodes after separator fix")
    d = d.set_index("_cistopic")
    cell_data[s] = d
    log(f"healthy manifest {s}: {len(d):,} cells | {d.Classified_Celltype.nunique()} cell types | "
        f"{d.Patient_ID.nunique()} patients")

malig_var = {}
for d, p in MALIG_REF.items():
    if not p.exists():
        die(f"missing malignant reference {p} — cannot be evaluated")
    _, _, v = read_auc(p, values=False)
    malig_var[d] = v
    log(f"malignant {d} vocabulary: {len(v)} regulons  ({p.name})")

malig_dead = {}
mp = REPO / "longitudinal/published" / "out" / "s3_provenance.json"
if mp.exists():
    j = json.load(open(mp))
    for tag, rec in j.get("files", {}).items():
        malig_dead[tag] = rec.get("dead_col_frac")

obs_by_cells = {}
pb_parts = {"Diagnosis": [], "Relapse": []}     # keyed by DEFS stage

for (cells, defs), path in FILES.items():
    X, obs, var = read_auc(path)
    tag = f"healthy_{cells}_cells_x_malignant_{defs}_defs"
    log(f"{tag}: {X.shape[0]:,} x {X.shape[1]}  ({path.name})")

    if X.shape[1] != EXPECTED_DEFS[defs]:
        die(f"column count for {tag}: got {X.shape[1]}, expected {EXPECTED_DEFS[defs]}")

    if len(var) != len(malig_var[defs]) or not np.array_equal(var, malig_var[defs]):
        die(f"vocabulary broken for {tag}: var does not match the malignant {defs} file")

    if cells in obs_by_cells:
        ref = obs_by_cells[cells]
        if len(ref) != len(obs) or not np.array_equal(ref, obs):
            die(f"cell identity broken for healthy {cells}: the two vocabularies' files "
                f"carry different cells (content or order)")
    else:
        obs_by_cells[cells] = obs

    cd = cell_data[cells]
    matched = np.isin(obs, cd.index.values)
    rate = float(matched.mean())
    log(f"  barcode match {rate:.4f} ({matched.sum():,}/{len(obs):,})")
    if rate < BARCODE_MATCH_MIN:
        die(f"barcode join for {tag}: {rate:.4f} < {BARCODE_MATCH_MIN} — separator or "
            f"manifest mismatch")

    df_dead = dead_frac(X)
    prov["files"][tag] = {
        "path": str(path), "bytes": path.stat().st_size,
        "n_cells": int(X.shape[0]), "n_regulons": int(X.shape[1]),
        "barcode_match_rate": rate,
        "dead_col_frac": float(df_dead),
        "auc_min": float(X.min()), "auc_max": float(X.max()),
        "auc_mean": float(X.mean()), "frac_exact_zero": float(np.mean(X == 0.0)),
    }
    log(f"  dead-column rate {df_dead:.4f}  |  malignant same-stage baseline "
        f"{malig_dead.get(f'{defs}_cells_x_{defs}_defs', float('nan'))}")

    # ---- pseudobulk: mean activity per (patient x cell-type), this stage, this vocabulary ----
    sub = cd.reindex(obs)
    keep = sub.Classified_Celltype.notna().values & sub.Patient_ID.notna().values
    cols = [regkey(v) for v in var]
    df = pd.DataFrame(X[keep], columns=cols)
    del X
    df["patient"] = sub.Patient_ID.values[keep]
    df["celltype"] = sub.Classified_Celltype.values[keep]
    grp = df.groupby(["patient", "celltype"], sort=True)
    mean = grp.mean()
    mean.insert(0, "n_cells", grp.size().rename("n_cells"))
    mean.insert(0, "stage", cells)
    pb_parts[defs].append(mean.reset_index())
    del df, grp, mean

# Driver map + the two definition-set tables

driver_map = {}
for s in ("Diagnosis", "Relapse"):
    agg = cell_data[s].groupby("Patient_ID").Driver.agg(
        lambda v: v.dropna().unique()[0] if v.dropna().nunique() == 1 else "AMBIGUOUS")
    for pid, drv in agg.items():
        driver_map.setdefault(pid, drv)

for defs in ("Diagnosis", "Relapse"):
    tab = pd.concat(pb_parts[defs], ignore_index=True)
    tab.insert(3, "driver", tab.patient.map(driver_map).fillna("UNLABELLED"))
    tab.to_csv(OUT / f"pseudobulk_healthy_def{defs}.csv", index=False)
    log(f"pseudobulk_healthy_def{defs}.csv: {len(tab):,} units x {tab.shape[1] - 5} regulons")
    prov["checks"][f"pseudobulk_healthy_def{defs}"] = {
        "n_units_total": int(len(tab)),
        "n_units_ge_floor": int((tab.n_cells >= MIN_CELLS_PER_UNIT).sum()),
        "n_regulons": int(tab.shape[1] - 5)}

# Paired-design audit + the cell-type overlap with the malignant admitted set

tab = pd.read_csv(OUT / "pseudobulk_healthy_defDiagnosis.csv",
                  usecols=["stage", "patient", "celltype", "driver", "n_cells"])
ok = tab[tab.n_cells >= MIN_CELLS_PER_UNIT]
dx = ok[ok.stage == "Diagnosis"].set_index(["patient", "celltype"]).index
rel = ok[ok.stage == "Relapse"].set_index(["patient", "celltype"]).index
paired_units = dx.intersection(rel)
paired_patients = sorted({p for p, _ in paired_units})
ct_counts = pd.Series([c for _, c in paired_units]).value_counts()

MALIG_ADMITTED = ["Early.Basophil", "Early.Erythrocyte", "GMP", "HSC", "Monocytes", "Progenitor"]
healthy_admitted = sorted(ct_counts[ct_counts >= 5].index)
overlap = sorted(set(healthy_admitted) & set(MALIG_ADMITTED))

prov["checks"]["design"] = {
    "paired_units_ge_floor": int(len(paired_units)),
    "paired_patients": len(paired_patients),
    "paired_celltypes": int(ct_counts.size),
    "celltype_patient_counts": {k: int(v) for k, v in ct_counts.items()},
    "healthy_admitted_ge5": healthy_admitted,
    "malignant_admitted": MALIG_ADMITTED,
    "admitted_overlap": overlap,
}
log(f"design: {len(paired_units)} paired (patient x cell-type) units >= {MIN_CELLS_PER_UNIT} "
    f"cells, {len(paired_patients)} patients, {ct_counts.size} cell types")
log("design: " + ", ".join(f"{k} ({v})" for k, v in ct_counts.items()))
log(f"design: healthy admitted (>=5 paired patients) = {healthy_admitted}")
log(f"design: OVERLAP with malignant admitted = {overlap}")

pd.DataFrame({"patient": paired_patients,
              "driver": [driver_map.get(p, "UNLABELLED") for p in paired_patients]}
             ).to_csv(OUT / "paired_patients_healthy.csv", index=False)
pd.DataFrame(list(paired_units), columns=["patient", "celltype"]).to_csv(
    OUT / "paired_units_healthy.csv", index=False)

with open(OUT / "s3b_provenance.json", "w") as fh:
    json.dump(prov, fh, indent=2, default=str)
log("wrote s3b_provenance.json")
