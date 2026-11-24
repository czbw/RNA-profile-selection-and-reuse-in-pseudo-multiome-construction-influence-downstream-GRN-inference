'Assemble the three pairing arms (native reference, link, pseudomultiome) into the inputs the ATAC/GRN chain expects.'
import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

def log(m):
    print(f"[p03] {m}", flush=True)

def die(m):
    raise SystemExit(f"[p03] ABORT -- {m}")

ap = argparse.ArgumentParser()
ap.add_argument("--config", required=True)
ap.add_argument("--dataset", required=True, choices=["3k", "10k"])
a = ap.parse_args()

ROOT = Path(os.environ.get("AMLGRN_ROOT", (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')))
cfg = yaml.safe_load(open(a.config, encoding="utf-8"))
STAGE = cfg["stages"][0]
WORK = ROOT / cfg["paths"]["work"]
OUT = ROOT / cfg["paths"]["out"]
DSWORK = WORK / a.dataset
FRAGDIR = ROOT / cfg["paths"]["atac_fragments"]
MINCONF = float(cfg["pairing_params"]["link"]["min_confidence"])
for d in (WORK, OUT, FRAGDIR):
    d.mkdir(parents=True, exist_ok=True)

# 1. inputs

gex_path = DSWORK / "gex_annotated.h5ad"
if not gex_path.exists():
    die(f"missing {gex_path} -- run p_00 then p_01 for {a.dataset} first")
link_csv = WORK / "link" / f"link_{a.dataset}.csv"
if not link_csv.exists():
    die(f"missing {link_csv} -- this is archr_link.R's output, produced locally in WSL. "
        f"It is not computed on CUBIC by design.")

G = sc.read_h5ad(gex_path)
log(f"GEX {G.n_obs:,} x {G.n_vars:,} (obs_names like {G.obs_names[0]})")
L = pd.read_csv(link_csv)
log(f"link table {len(L):,} rows from {link_csv.name}")

# cell ids must round-trip as <barcode><sep><sample_id>
SAMPLE = str(G.obs["sample_id"].iloc[0])
SEP = cfg["cell_id_sep"]
_expect = str(G.obs["barcode"].iloc[0]) + SEP + SAMPLE
if _expect != str(G.obs_names[0]):
    die(f"cell-id convention broken: rebuilt {_expect!r} != actual {G.obs_names[0]!r}")
log(f"sample_id {SAMPLE!r} | sep {SEP!r} | id round-trip verified")

bc2id = dict(zip(G.obs["barcode"].astype(str), G.obs_names.astype(str)))

# 2. confidence filter -- the same >=0.5 threshold as the AML analysis

conf = pd.to_numeric(L["Linked_scRNA_confidence_score"], errors="coerce")
n0 = len(L)
L = L[conf >= MINCONF].copy()
log(f"confidence >= {MINCONF}: {len(L):,} of {n0:,} ({100*len(L)/n0:.1f}%) "
    f"[median {conf.median():.3f} -- NOTE p_02 found this score saturates near 1.0 and does NOT "
    f"track correctness, so this filter is near-inert; kept only for parity with the AML arm]")

# keep rows whose ATAC cell AND linked RNA cell both survived GEX QC
L["atac_id"] = L["barcode"].map(bc2id)
L["rna_id"] = L["Linked_scRNA_barcode"].map(bc2id)
keep = L.atac_id.notna() & L.rna_id.notna()
log(f"both cells present after QC: {keep.sum():,} of {len(L):,}")
L = L[keep].copy()
if len(L) < 500:
    die(f"only {len(L)} usable linked cells -- something is wrong upstream")

# 3. cell_data, with the same columns as the AML cell_data tables

cd = pd.DataFrame({
    "barcode": L["barcode"].values,
    "sample_id": SAMPLE,
    "Classified_Celltype": G.obs["annotated_celltype"].reindex(L.atac_id).values,
    "Linked_scRNA_barcode": L["Linked_scRNA_barcode"].values,
    "Linked_scRNA_Sample": SAMPLE,
    "Linked_scRNA_confidence_score": L["Linked_scRNA_confidence_score"].values,
    "true_barcode": L["true_barcode"].values,
    "link_is_correct": L["link_is_correct"].values,
}, index=pd.Index(L.atac_id.values, name="cell_id"))
# the cell type of the RNA cell the link chose, alongside the ATAC cell's own type (p_02's CSV has
# no cell-type column, so this join is needed to compute cell-type agreement).
cd["Linked_Celltype"] = G.obs["annotated_celltype"].reindex(L.rna_id).values
cd.to_csv(WORK / f"cell_data_{STAGE}.csv")
log(f"cell_data_{STAGE}.csv: {len(cd):,} cells, {cd.Classified_Celltype.nunique()} cell types")

# Does the link land on the right cell type even when it misses the cell? This separates "wrong
# cell but right cell type" (tolerable for a cell-type-level GRN) from simply wrong.
ct_agree = float((cd.Classified_Celltype.astype(str) == cd.Linked_Celltype.astype(str)).mean())
exact = float(cd.link_is_correct.mean())
log(f"link EXACT cell   : {100*exact:6.2f}%  ({int(cd.link_is_correct.sum()):,}/{len(cd):,})")
log(f"link CELL TYPE    : {100*ct_agree:6.2f}%  <- the quantity that decides harmless-vs-wrong")

thin = cd.Classified_Celltype.value_counts()
floor = cfg["pairing_params"]["pseudomultiome"]["min_cells_per_celltype"]
if (thin < floor).any():
    log(f"WARNING: below the {floor}-cell floor and will be dropped before prepare_GEX_ACC: "
        f"{thin[thin < floor].to_dict()}")

# 4. fragments symlink, pycisTopic convention

# atac_chain.py reads only the keys of fragments_dict_<stage>.json and rebuilds each path as
# FRAG/f"{key}.fragments.tsv.gz", so the key must be the sample_id and the file must carry that name.
src = ROOT / cfg["datasets"][a.dataset]["fragments"]
for suf in ("", ".tbi"):
    d = FRAGDIR / f"{SAMPLE}.fragments.tsv.gz{suf}"
    s = Path(str(src) + suf)
    if not s.exists():
        die(f"missing {s}")
    if d.exists() or d.is_symlink():
        d.unlink()
    try:
        d.symlink_to(s)
    except OSError:
        import shutil; shutil.copy2(s, d)
frag_dict = {SAMPLE: str(FRAGDIR / f"{SAMPLE}.fragments.tsv.gz")}
(WORK / f"fragments_dict_{STAGE}.json").write_text(json.dumps(frag_dict, indent=2))
log(f"fragments -> {FRAGDIR}/{SAMPLE}.fragments.tsv.gz (+ .tbi)")
log(f"fragments_dict_{STAGE}.json -> keys {list(frag_dict)}")

# 5. the three RNA objects

def write_paired(rna_ids, tag):
    """GEX rows re-indexed to the ATAC cell ids -- multiome mode for run_grn.sh."""
    sub = G[rna_ids].copy()                 # duplicates are intentional: reuse is the phenomenon
    sub.obs_names = cd.index.values
    sub.obs = cd.copy()
    fp = OUT / f"paired_rna_{STAGE}.{tag}.h5ad"
    sub.write_h5ad(fp)
    nuniq = len(set(rna_ids))
    log(f"  {tag:6s}: {sub.n_obs:,} x {sub.n_vars:,} | distinct RNA {nuniq:,} "
        f"| reuse {sub.n_obs/nuniq:.2f}x -> {fp.name}")
    return {"n_cells": int(sub.n_obs), "distinct_rna": int(nuniq),
            "reuse": float(sub.n_obs / nuniq)}

# the native reference arm: each ATAC cell paired with its own RNA profile
truth_ids = cd["true_barcode"].map(bc2id)
if truth_ids.isna().any():
    die(f"{int(truth_ids.isna().sum())} true_barcodes do not resolve to a QC-surviving cell")
if not (truth_ids.values == cd.index.values).all():
    die("true_barcode does not resolve to the ATAC cell itself -- this is not a real multiome, "
        "or p_02 wrote true_barcode wrong. The truth arm would be silently mispaired.")
log(f"truth arm verified: all {len(cd):,} true_barcodes resolve to their own ATAC cell")

log("building the RNA objects:")
rep = {"dataset": a.dataset, "stage": STAGE, "min_confidence": MINCONF,
       "n_cells": int(len(cd)),
       "truth": write_paired(truth_ids.values, "truth"),
       "link": write_paired(L.rna_id.values, "link")}
rep["link"]["exact_cell_accuracy"] = exact
rep["link"]["celltype_agreement"] = ct_agree
rep["link"]["celltype_confusion"] = (
    pd.crosstab(cd.Classified_Celltype.astype(str), cd.Linked_Celltype.astype(str))
      .to_dict())

# pseudomultiome consumes the un-paired GEX; SCENIC+ builds metacells by group_by
nat = G.copy()
nat.obs["Classified_Celltype"] = nat.obs["annotated_celltype"]
nat.obs["sample_id"] = SAMPLE
nat.write_h5ad(OUT / f"native_rna_{STAGE}.h5ad")
log(f"  pseudo: {nat.n_obs:,} x {nat.n_vars:,} -> native_rna_{STAGE}.h5ad "
    f"(SCENIC+ draws metacells by {cfg['pairing_params']['pseudomultiome']['group_by']})")
rep["pseudomultiome"] = {"n_cells": int(nat.n_obs)}

json.dump(rep, open(WORK / f"pairings_report_{a.dataset}.json", "w"), indent=1)
log(f"link exact {100*exact:.2f}% | celltype {100*ct_agree:.2f}% | "
    f"truth reuse {rep['truth']['reuse']:.2f}x vs link reuse {rep['link']['reuse']:.2f}x")
log("STEP 03 COMPLETE")
