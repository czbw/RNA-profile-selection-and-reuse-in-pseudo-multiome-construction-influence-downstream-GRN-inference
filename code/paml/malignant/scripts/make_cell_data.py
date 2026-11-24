"""Build per-stage MALIGNANT cell_data + fragments_dict."""
import glob
import json
import os

import pandas as pd

from common import (cell_id, compartment_label, config_stages, load_config, log,
                    parse_sample_stem)

cfg = load_config()
META = cfg["P"]["atac_meta"]
FRAG = cfg["P"]["atac_fragments"]
WORK = cfg["P"]["work"]
SEP = cfg["cell_id_sep"]
STAGES = config_stages(cfg)
LABEL = compartment_label(cfg)          # "Malignant"
WORK.mkdir(parents=True, exist_ok=True)

KEEP = ["Classified_Celltype", "Patient_ID", "Subgroup", "Expected_Driving_Aberration",
        "Linked_scRNA_barcode", "Linked_scRNA_Sample", "Linked_scRNA_confidence_score"]

files = sorted(glob.glob(str(META / "*_filtered_metadata.tsv.gz")))
if not files:
    raise SystemExit(f"No metadata in {META} — stage the shared data/ (link_shared_assets.sh)")

buckets = {s: [] for s in STAGES}
dropped_stage = {}
for fp in files:
    base = os.path.basename(fp).replace("_filtered_metadata.tsv.gz", "")   # = sample_id
    df = pd.read_csv(fp, sep="\t", dtype=str)
    df = df[df["Malignant"].str.lower() == LABEL.lower()].copy()
    if df.empty:
        continue
    tp = df["Patient_Sample"].iloc[0]
    assert df["Patient_Sample"].nunique() == 1, f"{base}: mixed timepoint"
    if tp not in buckets:
        # A malignant sample whose curated timepoint is Remission (blast-empty compartment,
        # not a Stage-2 stage): count it, don't silently swallow it.
        dropped_stage[tp] = dropped_stage.get(tp, 0) + len(df)
        continue
    df["sample_id"] = base
    df["barcode"] = df["Cell_Barcode"]
    df["timepoint"] = tp
    df.index = [cell_id(bc, base, SEP) for bc in df["barcode"]]
    cols = ["barcode", "sample_id", "timepoint"] + KEEP
    buckets[tp].append(df[cols].rename(columns={"Expected_Driving_Aberration": "Driver"}))

if dropped_stage:
    log(f"malignant cells in non-Stage-2 timepoints (not written): {dropped_stage}")

for stage in STAGES:
    parts = buckets[stage]
    if not parts:
        log(f"{stage}: no malignant cells — skipped")
        continue
    cd = pd.concat(parts)
    cd.to_csv(WORK / f"cell_data_{stage}.csv")

    fd = {}
    for sid in sorted(cd["sample_id"].unique()):
        path = FRAG / f"{sid}.fragments.tsv.gz"
        if path.exists():
            fd[sid] = str(path)
        else:
            log(f"  [WARN] {stage}: converted fragments missing for {sid} ({path}) "
                f"— run prepare_fragments.sh / check the shared data/fragments symlink")
    (WORK / f"fragments_dict_{stage}.json").write_text(json.dumps(fd, indent=2))

    log(f"{stage}: {cd.shape[0]:,} malignant cells | {cd['sample_id'].nunique()} samples | "
        f"{cd['Patient_ID'].nunique()} patients | drivers {sorted(cd['Driver'].dropna().unique())} | "
        f"fragments resolved {len(fd)}/{cd['sample_id'].nunique()}")
    log(f"  -> {WORK / f'cell_data_{stage}.csv'}")

log("next: python build_paired_rna.py ; python build_native_rna.py ; then the ATAC chain")
