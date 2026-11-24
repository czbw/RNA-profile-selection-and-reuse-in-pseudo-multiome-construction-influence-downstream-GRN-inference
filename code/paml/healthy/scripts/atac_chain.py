"""ATAC half, per stage: pseudobulk, MACS3 peaks, consensus peaks, cisTopic object, topic model."""
import os
os.environ.setdefault("AMLGRN_CONFIG", "config.yaml")

import json
import pickle
import time

import pandas as pd

from common import load_config, log

cfg = load_config()
WORK = cfg["P"]["work"]; OUTP = cfg["P"]["out"]; DB = cfg["P"]["databases"]
FRAG = cfg["P"]["atac_fragments"]
STAGE = os.environ.get("AMLGRN_STAGE")
if not STAGE:
    raise SystemExit("Set AMLGRN_STAGE (Diagnosis|Remission|Relapse).")
SEP = "___"                                   # pycisTopic 2.x default; barcodes contain '-1'
MACS = os.environ.get("AMLGRN_MACS", cfg.get("macs_path", "macs3"))
TMP = WORK / "tmp"; TMP.mkdir(parents=True, exist_ok=True)
OUTP.mkdir(parents=True, exist_ok=True)

out_pkl = OUTP / f"cistopic_obj_{STAGE}.pkl"
if out_pkl.exists():
    log(f"{STAGE}: {out_pkl.name} exists — skip ATAC chain.")
    raise SystemExit(0)

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.lda_models import run_cgs_models, evaluate_models

clock = {}
def tic(): return time.time()
def toc(name, t0): clock[name] = round(time.time() - t0, 1); log(f"[{STAGE}/{name}] {clock[name]}s")

cd = pd.read_csv(WORK / f"cell_data_{STAGE}.csv", index_col=0)
cd.index = cd["barcode"].astype(str) + SEP + cd["sample_id"].astype(str)
# Drop the 'Unknown' (unclassified) label: it is not a real cell state, and at sparse stages it
# can yield zero MACS3 peaks. Removing the cells here excludes them from all downstream steps.
_n0 = len(cd)
cd = cd[cd["Classified_Celltype"] != "Unknown"].copy()
if len(cd) < _n0:
    log(f"{STAGE}: dropped {_n0 - len(cd):,} 'Unknown' cells -> {len(cd):,} remain")
fragments_dict = {k: str(FRAG / f"{k}.fragments.tsv.gz")
                  for k in json.loads((WORK / f"fragments_dict_{STAGE}.json").read_text())
                  if (FRAG / f"{k}.fragments.tsv.gz").exists()}
log(f"{STAGE}: {len(cd):,} healthy cells | {len(fragments_dict)} samples | "
    f"celltypes {sorted(cd['Classified_Celltype'].unique())}")

chromsizes = pd.read_csv(DB / "hg38.chrom.sizes", sep="\t", header=None, names=["Chromosome", "End"])
chromsizes["Start"] = 0
chromsizes = chromsizes[["Chromosome", "Start", "End"]]

GROUP = cfg["pseudobulk"]["group_by"]
(OUTP / f"pb_bed_{STAGE}").mkdir(parents=True, exist_ok=True); (OUTP / f"pb_bw_{STAGE}").mkdir(parents=True, exist_ok=True)

celltypes = sorted(cd[GROUP].unique())
existing = {ct: (OUTP / f"pb_bed_{STAGE}" / f"{ct}.fragments.tsv.gz") for ct in celltypes}
if all(p.exists() for p in existing.values()):
    bed_paths = {ct: str(p) for ct, p in existing.items()}
    log(f"{STAGE}/pseudobulk: reusing {len(bed_paths)} BEDs"); clock["pseudobulk"] = 0.0
else:
    t0 = tic()
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cd, variable=GROUP, chromsizes=chromsizes,
        bed_path=str(OUTP / f"pb_bed_{STAGE}"), bigwig_path=str(OUTP / f"pb_bw_{STAGE}"),
        path_to_fragments=fragments_dict, sample_id_col="sample_id",
        n_cpu=cfg.get("pseudobulk_n_cpu", 4),   # memory guard: deliberately not cfg["n_cpu"]
        normalize_bigwig=True, split_pattern=SEP, temp_dir=str(TMP))
    toc("pseudobulk", t0)
    # export_pseudobulk returns (bw_paths, bed_paths); its BED files are '{ct}.fragments.tsv.gz'.
    # Swap only if bed_paths actually point at the .bw bigwigs (guards against a reversed return).
    if all(str(p).endswith(".bw") for p in bed_paths.values()):
        bw_paths, bed_paths = bed_paths, bw_paths
json.dump(bed_paths, open(OUTP / f"bed_paths_{STAGE}.json", "w"), indent=2)

# --- 2. MACS3 peak calling per group (separate env; abs path via AMLGRN_MACS) ---
t0 = tic()
narrow = peak_calling(macs_path=MACS, bed_paths=bed_paths, outdir=str(OUTP / f"macs_{STAGE}"),
                      genome_size="hs", n_cpu=cfg["n_cpu"], input_format="BEDPE",
                      shift=73, ext_size=146, keep_dup="all", q_value=0.05,
                      skip_empty_peaks=True)   # a cell type with no peaks is skipped, not fatal
toc("peak_calling", t0)

# --- 3. consensus peaks ---
t0 = tic()
consensus = get_consensus_peaks(narrow_peaks_dict=narrow, peak_half_width=250, chromsizes=chromsizes)
consensus_bed = OUTP / f"consensus_peaks_{STAGE}.bed"
consensus.to_bed(path=str(consensus_bed), keep=False, compression="infer")
toc("consensus_peaks", t0)
log(f"{STAGE}: consensus peaks {len(consensus):,} -> {consensus_bed}")

# --- 4. cisTopic object (per sample -> merge) ---
t0 = tic()
objs = []
for sid, frag in fragments_dict.items():
    bcs = [c.split(SEP)[0] for c in cd.index[cd["sample_id"] == sid]]
    objs.append(create_cistopic_object_from_fragments(
        path_to_fragments=frag, path_to_regions=str(consensus_bed),
        valid_bc=bcs, n_cpu=cfg.get("pseudobulk_n_cpu", 4), project=sid, split_pattern=SEP))
obj = objs[0] if len(objs) == 1 else merge(objs)
overlap = len(set(obj.cell_names) & set(cd.index)) / max(1, len(obj.cell_names))
log(f"{STAGE}: cisTopic {len(obj.cell_names):,} cells x {len(obj.region_names):,} regions | overlap {overlap:.1%}")
assert overlap > 0.9, f"cell-name overlap {overlap:.1%} < 90% (split_pattern mismatch)"
obj.add_cell_data(cd.loc[obj.cell_names])
toc("cistopic_object", t0)

# --- 5. topic model (lda / collapsed Gibbs) ---
# This is the slowest step. If too slow, narrow n_topics / n_iter in config.
# n_cpu here = number of models run in parallel, which drives RAM use.
t0 = tic()
models = run_cgs_models(obj, n_topics=cfg["topic_model"]["n_topics"],
                        n_cpu=cfg.get("topic_n_cpu", 3), n_iter=cfg["topic_model"]["n_iter"],
                        random_state=cfg["seed"], save_path=str(OUTP / f"models_{STAGE}"))
obj.add_LDA_model(evaluate_models(models, return_model=True))
toc("topic_model", t0)

with open(out_pkl, "wb") as fh:
    pickle.dump(obj, fh)
log(f"=== {STAGE} ATAC timings (s): {clock} ===")
log(f"{STAGE}: cisTopic object -> {out_pkl}")
