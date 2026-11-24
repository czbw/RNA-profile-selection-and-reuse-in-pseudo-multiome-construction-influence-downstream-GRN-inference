"""Build per-stage HEALTHY cell_data + fragments_dict."""
import glob
import json
import os

import pandas as pd

from common import HEALTHY_LABEL, STAGES, cell_id, load_config, log, parse_sample_stem

cfg = load_config()
META = cfg["P"]["atac_meta"]
FRAG = cfg["P"]["atac_fragments"]
WORK = cfg["P"]["work"]
SEP = cfg["cell_id_sep"]
WORK.mkdir(parents=True, exist_ok=True)

KEEP = ["Classified_Celltype", "Patient_ID", "Subgroup", "Expected_Driving_Aberration",
        "Linked_scRNA_barcode", "Linked_scRNA_Sample", "Linked_scRNA_confidence_score"]

files = sorted(glob.glob(str(META / "*_filtered_metadata.tsv.gz")))
if not files:
    raise SystemExit(f"No metadata in {META} — run extract_atac.sh")

buckets = {s: [] for s in STAGES}
for fp in files:
    base = os.path.basename(fp).replace("_filtered_metadata.tsv.gz", "")   # = sample_id
    df = pd.read_csv(fp, sep="\t", dtype=str)
    df = df[df["Malignant"].str.lower() == HEALTHY_LABEL.lower()].copy()
    if df.empty:
        continue
    tp = df["Patient_Sample"].iloc[0]
    assert df["Patient_Sample"].nunique() == 1, f"{base}: mixed timepoint"
    df["sample_id"] = base
    df["barcode"] = df["Cell_Barcode"]
    df["timepoint"] = tp
    df.index = [cell_id(bc, base, SEP) for bc in df["barcode"]]
    cols = ["barcode", "sample_id", "timepoint"] + KEEP
    buckets[tp].append(df[cols].rename(columns={"Expected_Driving_Aberration": "Driver"}))

for stage in STAGES:
    parts = buckets[stage]
    if not parts:
        log(f"{stage}: no healthy cells — skipped")
        continue
    cd = pd.concat(parts)
    cd.to_csv(WORK / f"cell_data_{stage}.csv")

    # fragments_dict: only samples present at this stage, only if the converted file exists
    fd = {}
    for sid in sorted(cd["sample_id"].unique()):
        path = FRAG / f"{sid}.fragments.tsv.gz"
        if path.exists():
            fd[sid] = str(path)
        else:
            log(f"  [WARN] {stage}: converted fragments missing for {sid} ({path}) "
                f"— run prepare_fragments.sh")
    (WORK / f"fragments_dict_{stage}.json").write_text(json.dumps(fd, indent=2))

    log(f"{stage}: {cd.shape[0]:,} healthy cells | {cd['sample_id'].nunique()} samples | "
        f"{cd['Patient_ID'].nunique()} patients | fragments resolved {len(fd)}/{cd['sample_id'].nunique()}")
    log(f"  -> {WORK / f'cell_data_{stage}.csv'}")

log("next: python 04_pycistopic_consensus_peaks.py")
