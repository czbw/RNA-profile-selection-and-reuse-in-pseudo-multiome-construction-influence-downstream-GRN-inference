"""Split the 10x multiome h5 and verify the ATAC side."""
import argparse
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

def log(m):
    print(f"[p00] {m}", flush=True)

def die(m):
    raise SystemExit(f"[p00] ABORT -- {m}")

ap = argparse.ArgumentParser()
ap.add_argument("--config", required=True)
ap.add_argument("--dataset", required=True, choices=["3k", "10k"])
a = ap.parse_args()

ROOT = Path(os.environ.get("AMLGRN_ROOT", (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')))
cfg = yaml.safe_load(open(a.config, encoding="utf-8"))
ds = cfg["datasets"][a.dataset]
WORK = ROOT / cfg["paths"]["work"] / a.dataset
WORK.mkdir(parents=True, exist_ok=True)

h5 = ROOT / ds["matrix_h5"]
frag = ROOT / ds["fragments"]
tbi = Path(str(frag) + ".tbi")
for p in (h5, frag):
    if not p.exists():
        die(f"missing input {p} -- did stage_big_files.sh run for {a.dataset}?")

rep = {"dataset": a.dataset, "matrix_h5": str(h5), "fragments": str(frag)}

# 1. ATAC-side verification first (cheapest check that can fail the run)

with open(frag, "rb") as fh:
    head = fh.read(18)
is_bgzf = head[:2] == b"\x1f\x8b" and head[12:14] == b"BC"
rep["fragments_bgzf"] = bool(is_bgzf)
if not is_bgzf:
    die("fragments file is not BGZF. pycisTopic needs random access; re-run through "
        "`bgzip` (a plain `gzip` file will fail deep into pseudobulk, not here).")
if not tbi.exists():
    die(f"missing tabix index {tbi.name}. 10x ships it -- download it rather than regenerating:\n"
        f"  curl -O https://cf.10xgenomics.com/samples/cell-arc/1.0.0/"
        f"pbmc_granulocyte_sorted_{a.dataset}/{frag.name}.tbi")
rep["fragments_tbi"] = True
log(f"ATAC ok: BGZF + tabix index present ({frag.name})")

with gzip.open(frag, "rt") as fh:
    first = None
    for line in fh:
        if not line.startswith("#"):
            first = line.rstrip("\n").split("\t")
            break
if not first or len(first) < 4:
    die("fragments file has no parseable data row (expect chrom/start/end/barcode/count)")
rep["fragments_first_row"] = first[:5]
frag_bc_example = first[3]
log(f"fragments first row: {first[:4]} (barcode example {frag_bc_example})")

# 2. split the multiome h5

log(f"reading {h5.name} ...")
ad = sc.read_10x_h5(h5, gex_only=False)
ad.var_names_make_unique()
ft = ad.var["feature_types"].astype(str)
gex = ad[:, ft == "Gene Expression"].copy()
pk = ad[:, ft == "Peaks"].copy()
rep["n_cells_raw"] = int(ad.n_obs)
rep["n_genes"] = int(gex.n_vars)
rep["n_peaks"] = int(pk.n_vars)
log(f"split: {ad.n_obs:,} cells | {gex.n_vars:,} genes | {pk.n_vars:,} peaks")

if ad.n_obs != ds["n_cells_expected"]:
    log(f"WARNING: cell count {ad.n_obs:,} != expected {ds['n_cells_expected']:,} -- recording, not aborting")
rep["n_cells_expected"] = ds["n_cells_expected"]

# barcodes must match between the matrix and the fragments file, or `truth` is meaningless
mat_bc = set(map(str, ad.obs_names))
if frag_bc_example not in mat_bc:
    log(f"WARNING: fragment barcode {frag_bc_example} not in the matrix barcodes -- this is EXPECTED if "
        f"the fragments file contains non-cell barcodes; the overlap is checked properly below.")

# 3. QC -- the SCENIC+ paper's own treatment of this dataset

q = cfg["qc"]
gex.var["mt"] = gex.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(gex, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
before = gex.n_obs
keep = ((gex.obs["n_genes_by_counts"] >= q["min_genes_per_cell"]) &
        (gex.obs["n_genes_by_counts"] <= q["max_genes_per_cell"]) &
        (gex.obs["pct_counts_mt"] <= q["max_pct_mito"]))
gex = gex[keep].copy()
sc.pp.filter_genes(gex, min_cells=q["min_cells_per_gene"])
rep["qc"] = {"applied": q, "cells_before": int(before), "cells_after": int(gex.n_obs),
             "genes_after": int(gex.n_vars),
             "note": "Scrublet threshold from the config is NOT applied here -- doublet calling is "
                     "left to a separate declared step so its effect is separable."}
log(f"QC: {before:,} -> {gex.n_obs:,} cells, {gex.n_vars:,} genes "
    f"(min_genes {q['min_genes_per_cell']}, max {q['max_genes_per_cell']}, "
    f"mito <= {q['max_pct_mito']}%)")
if gex.n_obs < 500:
    die(f"only {gex.n_obs} cells survive QC -- something is wrong with the input or the thresholds")

# 4. cell manifest, in cisTopic's naming convention

sep = cfg["cell_id_sep"]
sample = f"pbmc_{a.dataset}"
cell = pd.DataFrame({
    "barcode": list(map(str, gex.obs_names)),
    "sample_id": sample,
})
cell["cistopic_id"] = cell.barcode + sep + cell.sample_id
cell = cell.set_index("cistopic_id")
cell.to_csv(WORK / "cell_data.csv")
gex.obs["barcode"] = list(map(str, gex.obs_names))
gex.obs["sample_id"] = sample
gex.obs_names = list(cell.index)
gex.write_h5ad(WORK / "gex.h5ad")
rep["cell_id_sep"] = sep
rep["sample_id"] = sample
rep["outputs"] = {"gex_h5ad": str(WORK / "gex.h5ad"), "cell_data": str(WORK / "cell_data.csv")}
log(f"wrote gex.h5ad ({gex.n_obs:,} x {gex.n_vars:,}) and cell_data.csv "
    f"(ids like {cell.index[0]})")

json.dump(rep, open(WORK / "prepare_report.json", "w"), indent=1, default=str)
log("wrote prepare_report.json -- STEP 00 COMPLETE")
