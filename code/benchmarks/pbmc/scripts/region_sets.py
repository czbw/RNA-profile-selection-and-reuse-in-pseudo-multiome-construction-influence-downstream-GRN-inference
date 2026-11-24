'Stage-1 Arm A — candidate enhancer region sets per stage (NO DB step; DB used only to guarantee overlap). Writes region_sets_{STAGE}/Topics_otsu/*.bed and .../DARs/*.bed for motif_enrichment_cistarget (run_grn.sh).'
import os
os.environ.setdefault("AMLGRN_CONFIG", "config.yaml")
import pickle
import re
import shutil

import numpy as np
import pandas as pd
import pyranges as pr
from common import load_config, log

cfg = load_config()
OUTP = cfg["P"]["out"]
DBP = cfg["P"]["databases"]
STAGE = os.environ.get("AMLGRN_STAGE")
if not STAGE:
    raise SystemExit("Set AMLGRN_STAGE.")
RS = OUTP / f"region_sets_{STAGE}"

from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import impute_accessibility, find_diff_features

# Only primary-assembly chromosomes carry cisTarget/SCREEN cCRE rankings.
PRIMARY_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
MIN_TOPIC_REGIONS = 20      # a topic must keep >= this many DB-overlapping regions to be written
MIN_DAR_REGIONS = 5         # same floor for DARs (they overlap heavily; this is just a guard)
REGION_RE = re.compile(r"(.+):(\d+)-(\d+)")
# The cisTarget rankings DB used by run_grn.sh ($DB/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather).
RANK = DBP / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"

GROUP_BY = cfg["pseudobulk"]["group_by"]
_rs = cfg.get("region_sets") or {}
# Env wins over config so ONE config serves both partitions: submit_bigmem.sh exports a budget
# large enough to keep every cell; submit_container.sh ('long', 240 G) leaves the safe default.
IMPUTE_BUDGET_GB = float(os.environ.get("AMLGRN_IMPUTE_BUDGET_GB", _rs.get("impute_mem_budget_gb", 64)))
MIN_CELLS_PER_GROUP = int(_rs.get("min_cells_per_group", 100))
SEED = int(cfg.get("seed", 555))

def _stratified_cap(cell_data, group_by, cap, min_per_group=MIN_CELLS_PER_GROUP, seed=SEED):
    """Pick <= ~cap cells, allocated proportionally across the `group_by` levels, with a floor
    so that rare cell types survive (DARs are called PER LEVEL — losing a level loses its DAR)."""
    rng = np.random.default_rng(seed)
    n = len(cell_data)
    if n <= cap:
        return list(cell_data.index)
    frac = cap / n
    keep = []
    for _, idx in cell_data.groupby(group_by, observed=True).groups.items():
        idx = list(idx)
        k = min(len(idx), max(min_per_group, int(round(len(idx) * frac))))
        keep.extend(idx if k >= len(idx) else list(rng.choice(idx, size=k, replace=False)))
    return keep

def _db_region_names(feather):
    """Region universe of the cisTarget DB. Cheap: reads only the Arrow schema (column names),
    not the 33 GB of data. Handles both DB orientations (regions-as-columns or regions-as-rows)."""
    import pyarrow as pa
    try:
        with pa.memory_map(str(feather), "r") as src:
            names = pa.ipc.open_file(src).schema.names
    except Exception as e:                                  # not an IPC file / older feather
        log(f"schema-only read failed ({e}); falling back to a full feather read")
        import pyarrow.feather as pf
        return [str(n) for n in pf.read_feather(feather).columns if REGION_RE.match(str(n))]
    regs = [n for n in names if REGION_RE.match(n)]          # regions-as-columns (aertslab default)
    if len(regs) < 10000:                                    # regions-as-rows: read the index column
        import pyarrow.feather as pf
        col = names[0]
        regs = [str(x) for x in pf.read_feather(feather, columns=[col])[col].tolist()
                if REGION_RE.match(str(x))]
    return regs

def _to_pr(region_names):
    """Parse 'chr:start-end' strings (primary chroms only) into a PyRanges."""
    rows = []
    for r in region_names:
        m = REGION_RE.match(str(r))
        if m and m.group(1) in PRIMARY_CHROMS:
            rows.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return pr.PyRanges(pd.DataFrame(rows, columns=["Chromosome", "Start", "End"]))

# ---- 1) DB region universe -> PyRanges (built once) ----
log(f"{STAGE}: loading cisTarget DB region universe from {RANK.name}")
db_regions = _db_region_names(RANK)
if len(db_regions) < 10000:
    raise SystemExit(f"DB region load looks wrong ({len(db_regions)} regions parsed) — check {RANK}")
DB_PR = _to_pr(db_regions)
log(f"{STAGE}: DB universe = {len(db_regions):,} regions ({len(DB_PR):,} on primary chroms)")

# ---- 2) expensive DAR/topic computation, cached so a resubmit skips it ----
CACHE = OUTP / f"region_raw_{STAGE}.pkl"
if CACHE.exists():
    log(f"{STAGE}: reusing cached DAR/topic computation ({CACHE.name}) — skipping impute/find_diff")
    with open(CACHE, "rb") as fh:
        region_bin, dars = pickle.load(fh)
else:
    with open(OUTP / f"cistopic_obj_{STAGE}.pkl", "rb") as fh:
        obj = pickle.load(fh)
    n_cells, n_regions = len(obj.cell_names), len(obj.region_names)
    log(f"{STAGE}: cisTopic {n_cells} cells x {n_regions} regions")

    budget_cells = max(1, int((IMPUTE_BUDGET_GB * 1e9) / (n_regions * 4)))
    log(f"{STAGE}: dense impute = {n_cells * n_regions * 4 / 1e9:.0f} GB for all cells; "
        f"budget {IMPUTE_BUDGET_GB:.0f} GB -> cap {budget_cells:,} cells")
    sel_cells = None
    if n_cells > budget_cells:
        sel_cells = _stratified_cap(obj.cell_data, GROUP_BY, budget_cells)
        log(f"{STAGE}: DAR branch subsampled {n_cells:,} -> {len(sel_cells):,} cells "
            f"(stratified by {GROUP_BY}, floor {MIN_CELLS_PER_GROUP}/type; topics still use ALL cells)")

    region_bin = binarize_topics(obj, target="region", method="otsu")
    imputed = impute_accessibility(obj, selected_cells=sel_cells, scale_factor=10**6)
    if sel_cells is not None:
        # keep find_diff_features' contrast groups aligned with the cells actually imputed
        obj.cell_data = obj.cell_data.loc[sel_cells]
    dars = find_diff_features(obj, imputed, variable=GROUP_BY, split_pattern="___")
    del imputed
    with open(CACHE, "wb") as fh:
        pickle.dump((region_bin, dars), fh)
    log(f"{STAGE}: cached DAR/topic computation -> {CACHE.name}")

# ---- 3) wipe + rewrite the region-set dirs, keeping only DB-overlapping regions ----
for sub in ("Topics_otsu", "DARs"):
    d = RS / sub
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

def write_bed(region_names, path, min_regions):
    """Write only candidate regions that OVERLAP the cisTarget DB (guarantees a non-empty,
    scoreable set). Returns the number of regions written (0 => not written)."""
    cand = _to_pr(region_names)
    if len(cand) == 0:
        return 0
    kept = cand.overlap(DB_PR)
    if len(kept) < min_regions:
        return 0
    kdf = kept.df if hasattr(kept, "df") else kept
    kdf[["Chromosome", "Start", "End"]].to_csv(path, sep="\t", header=False, index=False)
    return len(kept)

nt = nt_skip = 0
for topic, df in region_bin.items():
    if write_bed(df.index, RS / "Topics_otsu" / f"{topic}.bed", MIN_TOPIC_REGIONS):
        nt += 1
    else:
        nt_skip += 1
log(f"{STAGE}: wrote {nt} topic region sets ({nt_skip} dropped: <{MIN_TOPIC_REGIONS} DB-overlapping regions)")

nd = nd_skip = 0
for grp, df in dars.items():
    if df is not None and len(df):
        if write_bed(df.index, RS / "DARs" / f"{grp}.bed", MIN_DAR_REGIONS):
            nd += 1
        else:
            nd_skip += 1
log(f"{STAGE}: wrote {nd} DAR region sets ({nd_skip} dropped) -> {RS}")

if nt + nd == 0:
    raise SystemExit(f"{STAGE}: no region sets survived DB-overlap filtering — aborting (check DB build/path).")
