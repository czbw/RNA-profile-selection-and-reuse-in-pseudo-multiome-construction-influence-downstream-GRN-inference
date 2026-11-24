"""Build per-timepoint MALIGNANT-cell .h5ad RNA inputs for SCENIC+ (Arm A)."""
import glob
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.dirname(HERE)
RAW = os.environ.get("AMLGRN_RNA_RAW", str(Path(__file__).resolve().parents[3] / 'data/paml/data/paml/GSE235063_RAW'))
OUT = os.environ.get("AMLGRN_RNA_OUT", os.path.join(BUNDLE, "data", "rna"))
os.makedirs(OUT, exist_ok=True)

STAGES = ["Diagnosis", "Relapse"]            # malignant spine; no Remission
TOKEN2TP = {"DX": "Diagnosis", "REM": "Remission", "REL": "Relapse"}

# --- symbol -> Ensembl map from the shared raw reference (ensembl<TAB>symbol) ---
ref_genes = os.path.join(RAW, "GSM7494257_AML16_DX_raw_genes.tsv.gz")
ref = pd.read_csv(ref_genes, sep="\t", header=None, names=["ensembl", "symbol"])
sym2ens = dict(zip(ref["symbol"], ref["ensembl"]))
print(f"symbol->Ensembl map: {len(sym2ens):,} entries")

mtx_files = sorted(glob.glob(os.path.join(RAW, "*_processed_matrix.mtx.gz")))
all_stems = [os.path.basename(m).split("_processed_matrix")[0] for m in mtx_files]
if not all_stems:
    raise SystemExit(f"No *_processed_matrix.mtx.gz in {RAW} — check AMLGRN_RNA_RAW")

def load_sample(stem):
    """Return cells x genes AnnData of MALIGNANT cells, var=Ensembl, X=raw counts."""
    p = lambda suf: os.path.join(RAW, f"{stem}_processed_{suf}")
    M = sio.mmread(p("matrix.mtx.gz")).tocsr()            # genes x cells
    genes = pd.read_csv(p("genes.tsv.gz"), header=None)[0].values
    barcodes = pd.read_csv(p("barcodes.tsv.gz"), header=None)[0].values
    meta = pd.read_csv(p("metadata.tsv.gz"), sep="\t").set_index("Cell_Barcode")

    A = ad.AnnData(X=M.T.tocsr().astype(np.float32))      # cells x genes
    A.obs_names = barcodes
    A.var_names = genes
    A.obs = meta.reindex(A.obs_names)

    A = A[A.obs["Malignant"].astype(str).str.lower() == "malignant"].copy()
    if A.n_obs == 0:
        return None

    ens = pd.Series(A.var_names).map(sym2ens).values
    keep = pd.notna(ens)
    A = A[:, keep].copy()
    A.var["symbol"] = A.var_names
    A.var_names = ens[keep]
    A = A[:, ~A.var_names.duplicated()].copy()
    A.obs_names = [f"{stem}|{b}" for b in A.obs_names]
    A.obs["__token_tp"] = TOKEN2TP[stem.split("_")[-1]]
    return A

by_tp = {s: [] for s in STAGES}
skipped_rem = 0
for stem in all_stems:
    a = load_sample(stem)
    if a is None:
        continue
    tps = a.obs["Patient_Sample"].unique().tolist()
    assert len(tps) == 1, f"{stem}: mixed Patient_Sample {tps}"
    tp = tps[0]
    if tp != a.obs["__token_tp"].iloc[0]:
        print(f"  [token!=Patient_Sample] {stem}: token={a.obs['__token_tp'].iloc[0]} "
              f"-> using Patient_Sample={tp}")
    if tp not in by_tp:                          # malignant cells curated as Remission -> not a Stage-2 stage
        skipped_rem += a.n_obs
        continue
    by_tp[tp].append(a)

if skipped_rem:
    print(f"[note] {skipped_rem:,} malignant cells curated as Remission -> not written (blast-empty spine)")

for tp, parts in by_tp.items():
    print(f"\n=== {tp} ===")
    if not parts:
        print(f"  no malignant cells at {tp} — skipped")
        continue
    for a in parts:
        stem = a.obs_names[0].split("|")[0]
        print(f"  {stem}: {a.n_obs} malignant cells, {a.n_vars} mapped genes")
    adata = ad.concat(parts, join="outer", merge="first", fill_value=0)
    adata.var["symbol"] = adata.var.get("symbol")

    adata.layers["counts"] = adata.X.copy()       # raw counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)                            # X = log1p(CP10k)
    # sc.pp.log1p writes uns['log1p']={'base': None}. Recent anndata encodes that None as an
    # encoding_type='null' element that older anndata (as in the scenicplus environment) cannot
    # read (IORegistryError). Drop it; nothing downstream needs it.
    adata.uns.pop("log1p", None)

    adata.obs["timepoint"] = tp
    out = os.path.join(OUT, f"malignant_{tp.lower()}.h5ad")
    adata.write_h5ad(out)
    print(f"  -> {out}")
    print(f"     cells={adata.n_obs:,}  genes={adata.n_vars:,}  "
          f"patients={adata.obs['Patient_ID'].nunique()}  "
          f"drivers={sorted(adata.obs.get('Expected_Driving_Aberration', pd.Series(dtype=str)).dropna().unique())}")

print("\nnext: scp malignant_diagnosis.h5ad + malignant_relapse.h5ad to CUBIC data/rna (see MANIFEST).")
