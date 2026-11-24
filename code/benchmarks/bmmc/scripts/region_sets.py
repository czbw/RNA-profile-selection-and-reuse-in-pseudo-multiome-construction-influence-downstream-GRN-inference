#!/usr/bin/env python
"""The shared chromatin chain, step 2: candidate region sets (binarised topics + DARs)."""
from __future__ import annotations

import argparse
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

# Carried from the AML region-set step at its values. Not tuned here: only primary-assembly chromosomes carry cisTarget/SCREEN cCRE
# rankings, and a region set below the floor is dropped rather than risking the atomic abort.
PRIMARY_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
MIN_TOPIC_REGIONS = 20
MIN_DAR_REGIONS = 5
BINARIZE_METHOD = "otsu"          # the reference's method; config.yaml declares no knob
REGION_RE = re.compile(r"(.+):(\d+)-(\d+)")

def log(m):
    print(f"[b_02] {m}", flush=True)

def die(msg):
    print(f"[b_02] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def load_config():
    """See atac_chain.py: this bundle ships no common.py, so load_config()'s contract
    (read AMLGRN_CONFIG from scripts/, resolve `paths:` against AMLGRN_ROOT) is inlined."""
    name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    p = HERE / name
    if not p.is_file() or p.stat().st_size == 0:
        die(f"config {p} is missing or empty (AMLGRN_CONFIG={name})")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or "paths" not in cfg:
        die(f"{p} did not parse into a mapping with a `paths:` block")
    cfg["P"] = {k: (ROOT / v) for k, v in cfg["paths"].items()}
    cfg["CFG_PATH"] = p
    return cfg

def partof(final: Path) -> Path:
    return final.with_suffix(".part" + final.suffix)

def commit_bytes(final: Path, payload: bytes):
    tmp = partof(final)
    with open(tmp, "wb") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)

def commit_pickle(final: Path, obj):
    tmp = partof(final)
    try:
        with open(tmp, "wb") as fh:
            pickle.dump(obj, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, final)
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise

def nonempty_beds(d: Path):
    """The beds that actually carry regions. Existence of the directory proves nothing."""
    if not d.is_dir():
        return []
    return [p for p in sorted(d.rglob("*.bed")) if p.stat().st_size > 0]

def safe_name(group: str) -> str:
    """Bed filename for a group label. Only path separators are touched: `G/M prog` would
    otherwise become a directory. Everything else (spaces, '+', '-') is legal in a filename and
    is left alone, because this name reappears as the region-set name downstream."""
    return str(group).replace("/", "_").replace(os.sep, "_").strip()

def db_region_names(feather: Path):
    """Region universe of the cisTarget DB. Cheap: reads only the Arrow schema (column names),"""
    import pyarrow as pa
    try:
        with pa.memory_map(str(feather), "r") as src:
            names = pa.ipc.open_file(src).schema.names
    except Exception as e:
        log(f"schema-only read failed ({e}); falling back to a full feather read")
        import pyarrow.feather as pf
        return [str(n) for n in pf.read_feather(feather).columns if REGION_RE.match(str(n))]
    regs = [n for n in names if REGION_RE.match(n)]          # regions-as-columns (aertslab)
    if len(regs) < 10000:                                    # regions-as-rows: read the index
        import pyarrow.feather as pf
        col = names[0]
        regs = [str(x) for x in pf.read_feather(feather, columns=[col])[col].tolist()
                if REGION_RE.match(str(x))]
    return regs

def to_pyranges(region_names):
    """Parse 'chr:start-end' strings (primary chroms only) into a PyRanges."""
    import pyranges as pr
    rows = []
    for r in region_names:
        m = REGION_RE.match(str(r))
        if m and m.group(1) in PRIMARY_CHROMS:
            rows.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return pr.PyRanges(pd.DataFrame(rows, columns=["Chromosome", "Start", "End"]))

def stratified_cap(cell_data, group_by, cap, min_per_group, seed):
    """Pick <= ~cap cells, allocated proportionally across the `group_by` levels, with a floor
    so rare cell types survive (DARs are called per level -- losing a level loses its DAR).
    At BMMC's 69,249 x 116,490 the dense impute is ~32 GB against a 200 GB budget, so this
    should never fire; it is kept so the script behaves identically if it is ever run smaller."""
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

def run_motif_enrichment(cfg, stage: str, rgs: Path, db: Path, tmp: Path):
    """Standalone-only. Byte-for-byte the invocation run_bmmc.sh phase [2c] uses, including
    --n_cpu 1 and --annotation_version v10nr_clust, written through a .part sidecar."""
    menr = cfg["P"]["out"] / f"menr_cistarget_{stage}.hdf5"     
    if menr.is_file() and menr.stat().st_size > 0:
        log(f"motif enrichment present ({menr.name}) -- skip")
        return menr
    rank = db / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    tbl = db / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    for f in (rank, tbl):
        if not f.is_file() or f.stat().st_size == 0:
            die(f"missing cisTarget database file {f}")
    tmp.mkdir(parents=True, exist_ok=True)
    part = partof(menr)
    cmd = ["scenicplus", "grn_inference", "motif_enrichment_cistarget",
           "--region_set_folder", str(rgs),
           "--cistarget_db_fname", str(rank),
           "--output_fname_cistarget_result", str(part),
           "--temp_dir", str(tmp), "--species", "homo_sapiens",
           "--path_to_motif_annotations", str(tbl),
           "--annotation_version", "v10nr_clust", "--n_cpu", "1"]
    log("motif_enrichment_cistarget (standalone; under run_bmmc.sh this is phase [2c])")
    log("  " + " ".join(cmd))
    rc = subprocess.call(cmd)
    if rc != 0 or not part.is_file() or part.stat().st_size == 0:
        if part.exists():
            part.unlink()
        die(f"motif_enrichment_cistarget failed (rc={rc}); no output was left behind")
    os.replace(part, menr)
    log(f"committed {menr.name} ({menr.stat().st_size / 1e6:.0f} MB)")
    return menr

def main():
    ap = argparse.ArgumentParser(description="BMMC shared chromatin chain: region sets")
    ap.add_argument("--stage", default=os.environ.get("AMLGRN_STAGE", "bmmc"))
    ap.add_argument("--cistopic", default=None,
                    help="override out/cistopic_obj_<STAGE>.pkl")
    ap.add_argument("--motif-enrichment", action="store_true", dest="menr",
                    help="ALSO run motif_enrichment_cistarget. Do NOT pass this under "
                         "run_bmmc.sh: phase [2c] owns that artifact.")
    ap.add_argument("--force", action="store_true",
                    help="rebuild the region sets even if a complete set already exists")
    a = ap.parse_args()

    cfg = load_config()
    STAGE = a.stage
    OUTP = cfg["P"]["out"]
    WORK = cfg["P"]["work"]
    OUTP.mkdir(parents=True, exist_ok=True)
    DB = Path(os.environ.get("AMLGRN_DB", str(cfg["P"]["databases"])))
    SEP = str(cfg["cell_id_sep"])
    SEED = int(os.environ.get("AMLGRN_SEED", cfg["seed"]))
    GROUP_BY = cfg["pseudobulk"]["group_by"]
    _rs = cfg.get("region_sets") or {}
    # env wins over config so one config serves both partitions (a large-memory submit script can
    # set a budget large enough to keep every cell; a smaller node leaves the config default).
    BUDGET_GB = float(os.environ.get("AMLGRN_IMPUTE_BUDGET_GB",
                                     _rs.get("impute_mem_budget_gb", 64)))
    MIN_CELLS_PER_GROUP = int(_rs.get("min_cells_per_group", 100))

    RS = OUTP / f"region_sets_{STAGE}"
    # `.part.d` (not a bare `.part`) so it matches run_bmmc.sh's `out/*.part.*`
    # glob and is reported as a transient sidecar that is always safe to delete.
    RS_PART = OUTP / f"region_sets_{STAGE}.part.d"
    CACHE = OUTP / f"region_raw_{STAGE}.pkl"
    REP = OUTP / f"b_02_region_sets_report_{STAGE}.json"
    CIS = Path(a.cistopic) if a.cistopic else OUTP / f"cistopic_obj_{STAGE}.pkl"
    TMP = WORK / f"grn_tmp_{STAGE}_shared"

    log(f"stage {STAGE} | seed {SEED} | group_by {GROUP_BY} | config {cfg['CFG_PATH'].name}")
    log(f"impute budget {BUDGET_GB:.0f} GB (env AMLGRN_IMPUTE_BUDGET_GB wins over config), "
        f"floor {MIN_CELLS_PER_GROUP} cells/group")

    # ---- already done? content, not existence ----------------------------------------------
    have = nonempty_beds(RS)
    if have and not a.force:
        log(f"SKIP: {RS.name}/ already holds {len(have)} non-empty bed file(s) "
            f"({sum(p.stat().st_size for p in have) / 1e6:.1f} MB total); rebuilding it would "
            f"move the frozen chain.")
        if a.menr:
            run_motif_enrichment(cfg, STAGE, RS, DB, TMP)
        return
    if RS.exists() and not have:
        die(f"{RS} exists but contains NO non-empty .bed file. run_bmmc.sh phase [2b] guards on "
            f"`[ -d ]`, so this directory would be accepted as complete and an empty region set "
            f"would be frozen into the chain. Delete it and resubmit (or pass --force).")
    if RS_PART.exists():
        log(f"removing stale {RS_PART.name}/ from a killed run")
        shutil.rmtree(RS_PART)

    # ---- 1) cisTarget DB region universe ----------------------------------------------------
    RANK = DB / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    if not RANK.is_file() or RANK.stat().st_size == 0:
        die(f"missing cisTarget rankings {RANK} (AMLGRN_DB={DB})")
    log(f"loading the cisTarget DB region universe from {RANK.name}")
    db_regions = db_region_names(RANK)
    if len(db_regions) < 10000:
        die(f"DB region load looks wrong ({len(db_regions)} regions parsed) -- check {RANK}")
    DB_PR = to_pyranges(db_regions)
    if len(DB_PR) == 0:
        die("the DB region universe parsed to ZERO primary-chromosome intervals")
    log(f"DB universe = {len(db_regions):,} regions ({len(DB_PR):,} on primary chroms)")

    # ---- 2) the expensive half, cached so a resubmit skips it -------------------------------
    clock = {}
    if CACHE.is_file() and CACHE.stat().st_size > 0:
        log(f"reusing cached topic/DAR computation ({CACHE.name}) -- skipping impute/find_diff")
        try:
            with open(CACHE, "rb") as fh:
                region_bin, dars = pickle.load(fh)
        except Exception as e:
            die(f"{CACHE} exists but does NOT unpickle ({type(e).__name__}: {e}) -- that is a "
                f"TRUNCATED cache. Delete it and resubmit.")
        if not isinstance(region_bin, dict) or not region_bin:
            die(f"{CACHE.name} holds no binarised topics -- delete it and resubmit")
    else:
        if not CIS.is_file() or CIS.stat().st_size == 0:
            die(f"missing {CIS} -- run_bmmc.sh phase [2a] (atac_chain.py) builds it; this "
                f"script never does")
        from pycisTopic.topic_binarization import binarize_topics
        from pycisTopic.diff_features import impute_accessibility, find_diff_features

        log(f"loading {CIS.name} ({CIS.stat().st_size / 1e9:.2f} GB)")
        try:
            with open(CIS, "rb") as fh:
                obj = pickle.load(fh)
        except Exception as e:
            die(f"{CIS} exists but does NOT unpickle ({type(e).__name__}: {e}) -- a TRUNCATED "
                f"cisTopic object. Delete it, let phase [2a] rebuild it, and do not let "
                f"chain_guard freeze the truncation.")
        if getattr(obj, "selected_model", None) is None:
            die(f"{CIS.name} carries no fitted LDA model (selected_model is None); "
                f"binarize_topics has nothing to binarise")
        n_cells, n_regions = len(obj.cell_names), len(obj.region_names)
        ds = cfg["dataset"]
        if (n_cells, n_regions) != (int(ds["n_cells"]), int(ds["n_atac_peaks"])):
            die(f"{CIS.name} is {n_cells:,} x {n_regions:,}, config declares "
                f"{int(ds['n_cells']):,} x {int(ds['n_atac_peaks']):,}")
        if GROUP_BY not in obj.cell_data.columns:
            die(f"{CIS.name} cell_data has no `{GROUP_BY}` column -- find_diff_features groups "
                f"on it. Columns: {list(obj.cell_data.columns)}")
        log(f"cisTopic {n_cells:,} cells x {n_regions:,} regions | "
            f"{obj.cell_data[GROUP_BY].nunique()} {GROUP_BY} levels")

        # memory guard: impute_accessibility() materialises a dense cells x regions matrix and
        # find_diff_features() then copies it, which can exhaust memory on large datasets. The cap
        # bounds the dense matrix by capping the cells fed to the DAR branch, stratified so no cell
        # type is lost; the topics always use every cell.
        budget_cells = max(1, int((BUDGET_GB * 1e9) / (n_regions * 4)))
        log(f"dense impute = {n_cells * n_regions * 4 / 1e9:.0f} GB for all cells; "
            f"budget {BUDGET_GB:.0f} GB -> cap {budget_cells:,} cells")
        sel_cells = None
        if n_cells > budget_cells:
            sel_cells = stratified_cap(obj.cell_data, GROUP_BY, budget_cells,
                                       MIN_CELLS_PER_GROUP, SEED)
            log(f"DAR branch subsampled {n_cells:,} -> {len(sel_cells):,} cells (stratified by "
                f"{GROUP_BY}, floor {MIN_CELLS_PER_GROUP}/type; topics still use ALL cells)")
        else:
            log("every cell reaches the DAR step (no stratified cap needed)")

        t0 = time.time()
        region_bin = binarize_topics(obj, target="region", method=BINARIZE_METHOD)
        clock["binarize_topics"] = round(time.time() - t0, 1)
        if not region_bin:
            die("binarize_topics returned no topics")
        log(f"binarize_topics({BINARIZE_METHOD}): {len(region_bin)} topics "
            f"({clock['binarize_topics']}s)")

        t0 = time.time()
        imputed = impute_accessibility(obj, selected_cells=sel_cells, scale_factor=10 ** 6)
        if sel_cells is not None:
            # keep find_diff_features' contrast groups aligned with the cells actually imputed
            obj.cell_data = obj.cell_data.loc[sel_cells]
        dars = find_diff_features(obj, imputed, variable=GROUP_BY, split_pattern=SEP)
        del imputed
        clock["impute_and_dars"] = round(time.time() - t0, 1)
        if not dars:
            die(f"find_diff_features returned no groups for `{GROUP_BY}`")
        log(f"find_diff_features: {len(dars)} groups ({clock['impute_and_dars']}s)")
        commit_pickle(CACHE, (region_bin, dars))
        log(f"cached topic/DAR computation -> {CACHE.name}")

    # ---- 3) write the beds into the .part.d staging tree ------------------------------------
    for sub in ("Topics_otsu", "DARs"):
        (RS_PART / sub).mkdir(parents=True, exist_ok=True)

    def write_bed(region_names, path: Path, min_regions: int) -> int:
        """Write only candidate regions that overlap the cisTarget DB, guaranteeing a non-empty,
        scoreable set. Returns the number of regions written (0 => nothing written)."""
        cand = to_pyranges(region_names)
        if len(cand) == 0:
            return 0
        kept = cand.overlap(DB_PR)
        if len(kept) < min_regions:
            return 0
        kdf = kept.df if hasattr(kept, "df") else kept
        kdf[["Chromosome", "Start", "End"]].to_csv(path, sep="\t", header=False, index=False)
        if not path.is_file() or path.stat().st_size == 0:
            die(f"wrote {path.name} but it is empty -- refusing to stage a degenerate bed")
        return len(kept)

    name_map, dropped = {}, {"topics": [], "dars": []}

    nt = 0
    for topic, df in sorted(region_bin.items()):
        fn = safe_name(topic) + ".bed"
        if fn in name_map:
            die(f"bed filename collision on {fn!r} (from {topic!r} and {name_map[fn]!r})")
        n = write_bed(df.index, RS_PART / "Topics_otsu" / fn, MIN_TOPIC_REGIONS)
        if n:
            name_map[fn] = str(topic)
            nt += 1
        else:
            dropped["topics"].append(str(topic))
    log(f"wrote {nt} topic region sets ({len(dropped['topics'])} dropped: "
        f"<{MIN_TOPIC_REGIONS} DB-overlapping regions)")

    nd = 0
    for grp, df in sorted(dars.items(), key=lambda kv: str(kv[0])):
        if df is None or not len(df):
            dropped["dars"].append(str(grp))
            continue
        fn = safe_name(grp) + ".bed"
        if fn in name_map:
            die(f"bed filename collision on {fn!r} (from {grp!r} and {name_map[fn]!r}) -- two "
                f"group labels differing only in a path separator")
        n = write_bed(df.index, RS_PART / "DARs" / fn, MIN_DAR_REGIONS)
        if n:
            name_map[fn] = str(grp)
            nd += 1
        else:
            dropped["dars"].append(str(grp))
    log(f"wrote {nd} DAR region sets ({len(dropped['dars'])} dropped)")
    renamed = {k: v for k, v in name_map.items() if k[:-4] != v}
    if renamed:
        log(f"{len(renamed)} bed name(s) differ from their group label (path separators are not "
            f"legal in a filename): " + ", ".join(f"{v!r} -> {k}" for k, v in renamed.items()))

    if nt + nd == 0:
        shutil.rmtree(RS_PART)
        die("no region set survived DB-overlap filtering -- aborting (check the DB build/path). "
            "Nothing was committed.")
    staged = nonempty_beds(RS_PART)
    if len(staged) != nt + nd:
        die(f"staged {len(staged)} non-empty beds but counted {nt + nd} -- refusing to commit")

    # ---- 4) commit the whole tree with one rename -------------------------------------------
    if RS.exists():
        shutil.rmtree(RS)
    os.replace(RS_PART, RS)
    final = nonempty_beds(RS)
    if len(final) != len(staged):
        die(f"post-commit check: {len(final)} non-empty beds in {RS.name}, expected {len(staged)}")
    log(f"committed {len(final)} region sets -> {RS}")

    rep = {
        "stage": STAGE, "seed": SEED, "group_by": GROUP_BY,
        "cistopic_obj": str(CIS),
        "binarize_method": BINARIZE_METHOD,
        "min_topic_regions": MIN_TOPIC_REGIONS, "min_dar_regions": MIN_DAR_REGIONS,
        "impute_mem_budget_gb": BUDGET_GB, "min_cells_per_group": MIN_CELLS_PER_GROUP,
        "n_topic_sets": nt, "n_dar_sets": nd,
        "dropped": dropped,
        "bed_name_to_group": name_map,
        "timings_s": clock,
    }
    commit_bytes(REP, (json.dumps(rep, indent=1) + "\n").encode("utf-8"))

    if a.menr:
        run_motif_enrichment(cfg, STAGE, RS, DB, TMP)
    log("b_02 COMPLETE")

if __name__ == "__main__":
    main()
