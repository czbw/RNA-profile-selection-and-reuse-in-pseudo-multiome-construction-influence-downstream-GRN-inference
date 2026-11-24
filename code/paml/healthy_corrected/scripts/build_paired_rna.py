#!/usr/bin/env python
"""Build paired healthy-cell RNA objects from the repaired healthy cell_data tables.

Writes paired_rna_<stage>.healthy_repair.h5ad and paired_cells_<stage>.healthy_repair.txt.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

# Input/output locations come from the environment. Every write below is suffixed
# `.healthy_repair` so it cannot overwrite an existing output.
OUT = Path(os.environ["AMLGRN_HEALTHY_OUT"])        # CUBIC: $ROOT/out
WORK = Path(os.environ["AMLGRN_HEALTHY_WORK"])      # CUBIC: $ROOT/work
RNA_DIR = Path(os.environ["AMLGRN_HEALTHY_RNA"])    # CUBIC: $ROOT/data/rna
STAGES = os.environ.get("AMLGRN_STAGES", "Diagnosis Relapse").split()
SEP = os.environ.get("AMLGRN_CELL_ID_SEP", "___")
MINCONF = float(os.environ.get("AMLGRN_MIN_CONF", "0.5"))
ALL_RNA_STAGES = ["Diagnosis", "Remission", "Relapse"]

KEEP_OBS = ["sample_id", "Classified_Celltype", "Patient_ID", "Subgroup", "Driver",
            "timepoint", "Linked_scRNA_confidence_score"]

def log(m: str) -> None:
    print(f"[healthy-repair-gex] {m}", flush=True)

def die(m: str) -> None:
    raise SystemExit(f"[healthy-repair-gex] ABORT - {m}")

def rna_lookup() -> pd.DataFrame:
    """Map (Library_ID, Cell_Barcode) -> (stage, obs_name) across the three healthy RNA h5ads.

    Backed mode: reads .obs only, never loads X. All three are searched even when only two stages
    are built, because a library may be named for one timepoint while the deposit assigns it
    another (the AML1_DX/AML1_REM swap; see p_h00). Pairing is within-timepoint by Patient_Sample.
    """
    rows = []
    for s in ALL_RNA_STAGES:
        fp = RNA_DIR / f"healthy_{s.lower()}.h5ad"
        if not fp.exists():
            die(f"missing RNA h5ad {fp}")
        A = ad.read_h5ad(fp, backed="r")
        obs = A.obs
        barcode = pd.Index(obs.index).str.split("|").str[-1]     # obs_names are "<stem>|<barcode>"
        rows.append(pd.DataFrame({"lib": obs["Library_ID"].astype(str).values,
                                  "barcode": barcode,
                                  "obs_name": obs.index.values,
                                  "rna_stage": s}))
        log(f"  {fp.name}: {obs.shape[0]:,} cells")
        A.file.close()
    allr = pd.concat(rows, ignore_index=True)
    allr["key"] = allr["lib"] + "||" + allr["barcode"]
    dup = int(allr["key"].duplicated().sum())
    if dup:
        die(f"{dup} duplicate (library, barcode) keys in the healthy RNA -- the join would be "
            f"ambiguous and the original builder's 'keep first' is not acceptable here")
    return allr.set_index("key")

def cistopic_names(stage: str) -> set | None:
    """Authoritative cell names, read from the existing link mudata if it is present."""
    fp = OUT / f"mudata_{stage}_link.h5mu"
    if not fp.is_file():
        return None
    try:
        with h5py.File(fp, "r") as f:
            g = f["mod"]["scATAC"] if "scATAC" in f["mod"] else f["mod"][list(f["mod"])[0]]
            idx = g["obs"].attrs["_index"]
            arr = g["obs"][idx][:]
        return {x.decode() if isinstance(x, bytes) else str(x) for x in arr}
    except Exception as e:                                        # pragma: no cover
        log(f"  [warn] could not read cell names from {fp.name}: {e}")
        return None

def build_stage(stage: str, lut: pd.DataFrame) -> None:
    src = WORK / f"cell_data_{stage}.healthy_repaired.csv"
    if not src.is_file():
        die(f"missing {src} -- scp it from paml/healthy_corrected/out/")
    cd = pd.read_csv(src, index_col=0, dtype=str)
    n0 = len(cd)

    # the confidence gate must be a no-op under the placeholder
    conf = pd.to_numeric(cd["Linked_scRNA_confidence_score"], errors="coerce")
    if conf.isna().any():
        die(f"{stage}: non-numeric Linked_scRNA_confidence_score")
    n_gated = int((conf < MINCONF).sum())
    if n_gated:
        die(f"{stage}: {n_gated:,} rows below the {MINCONF} confidence gate. p_h01 writes a 1.0 "
            f"placeholder, so this file is NOT the repaired table -- check what was staged.")

    # rebuild cisTopic-style cell ids from barcode and sample
    ids = cd["barcode"].astype(str).str.strip() + SEP + cd["sample_id"].astype(str).str.strip()
    if ids.duplicated().any():
        die(f"{stage}: duplicate cell ids after rebuilding with separator '{SEP}'")
    cd = cd.assign(_cid=ids.values).set_index("_cid")

    ref = cistopic_names(stage)
    if ref is None:
        log(f"  [warn] {stage}: mudata_{stage}_link.h5mu absent -- cannot verify the cell-id "
            f"convention here. prepare_GEX_ACC will fail loudly if '{SEP}' is wrong.")
    else:
        inter = [c for c in cd.index if c in ref]
        hit = len(inter) / max(1, len(cd))       # repaired ids known to the old arm
        kept = len(inter) / max(1, len(ref))     # old-arm cells surviving into the repaired set
        n_new = len(cd) - len(inter)
        log(f"  {stage}: cell-id convention vs mudata_{stage}_link.h5mu -- "
            f"{hit:.2%} of {len(cd):,} repaired ids are known cisTopic cells; "
            f"{kept:.2%} of the {len(ref):,} original cells survive; "
            f"{n_new:,} repaired cells are new to this arm")
        # A wrong separator gives ~0 overlap; any substantial overlap shows the ids parse correctly.
        if hit < 0.50:
            die(f"{stage}: only {hit:.2%} of rebuilt cell ids parse into the cisTopic naming. "
                f"A correct separator gives PARTIAL overlap (the cohort legitimately changes); a "
                f"wrong one gives ~0. This is at the wrong end of that range -- check "
                f"AMLGRN_CELL_ID_SEP='{SEP}'.")
        # The repair only removes excluded-unit and no-donor cells, so most original cells must
        # remain; a large loss means the wrong table was staged.
        if kept < 0.90:
            die(f"{stage}: only {kept:.2%} of the {len(ref):,} original operative cells survive "
                f"into the repaired set. The repair drops ONLY declared-exclusion and no-donor "
                f"cells (p_h02 measured 5.72% at Diagnosis, 1.69% at Relapse), so losing more "
                f"than 10% means the staged table is not the one p_h02 verified. Refusing to "
                f"build a GEX matrix over a cohort nobody checked.")

    # every donor key must resolve
    key = (cd["Linked_scRNA_Sample"].astype(str).str.strip() + "||"
           + cd["Linked_scRNA_barcode"].astype(str).str.strip())
    cd["rna_obs"] = key.map(lut["obs_name"]).values
    cd["rna_stage"] = key.map(lut["rna_stage"]).values
    matched = cd[cd["rna_obs"].notna()].copy()
    rate = len(matched) / max(1, n0)
    log(f"  {stage}: {n0:,} repaired ATAC cells -> {len(matched):,} RNA-matched ({rate:.4%}); "
        f"linked-RNA-stage mix {matched['rna_stage'].value_counts().to_dict()}")
    if len(matched) != n0:
        die(f"{stage}: {n0-len(matched):,} donors did not resolve in this cluster's healthy RNA. "
            f"p_h02 verified every key resolves locally, so the uploaded table and the RNA files "
            f"here disagree. Refusing to silently shrink the cohort.")

    # sanity: the corrected pairing is within-timepoint, so every donor must come from this stage
    off = matched.loc[matched["rna_stage"] != stage]
    if len(off):
        die(f"{stage}: {len(off):,} donors resolved to another timepoint's h5ad "
            f"({off['rna_stage'].value_counts().to_dict()}) -- the within-timepoint constraint is "
            f"broken in the staged table")

    parts = []
    for src_stage, grp in matched.groupby("rna_stage"):
        A = ad.read_h5ad(RNA_DIR / f"healthy_{str(src_stage).lower()}.h5ad", backed="r")
        sub = A[grp["rna_obs"].values].to_memory()
        A.file.close()
        sub.obs_names = grp.index.values          # the ATAC cell id, 1:1 with grp order
        sub.obs = grp[KEEP_OBS].copy()
        sub.var["ensembl"] = sub.var_names
        sym = sub.var["symbol"].astype(str)
        sub = sub[:, sym.notna() & (sym != "nan")].copy()
        sub.var_names = sub.var["symbol"].astype(str)
        sub.var_names_make_unique()
        parts.append(sub)
    paired = parts[0] if len(parts) == 1 else ad.concat(parts, join="outer", fill_value=0,
                                                        merge="first")
    paired.var = paired.var[[c for c in ["ensembl", "symbol"] if c in paired.var]].astype(str)
    for c in paired.obs.columns:
        if paired.obs[c].dtype == object:
            paired.obs[c] = paired.obs[c].astype(str)

    dest = OUT / f"paired_rna_{stage}.healthy_repair.h5ad"
    tmp = dest.with_suffix(".part.h5ad")
    paired.write_h5ad(tmp)
    tmp.replace(dest)
    pd.Series(paired.obs_names).to_csv(OUT / f"paired_cells_{stage}.healthy_repair.txt",
                                       index=False, header=False)
    log(f"  -> {dest.name}  ({paired.n_obs:,} x {paired.n_vars:,})")

def main() -> int:
    log(f"stages: {STAGES} | cell-id separator '{SEP}' | out {OUT}")
    for p in (OUT, WORK, RNA_DIR):
        if not p.is_dir():
            die(f"missing directory {p}")
    lut = rna_lookup()
    log(f"healthy RNA lookup: {len(lut):,} cells across {ALL_RNA_STAGES}")
    for s in STAGES:
        dest = OUT / f"paired_rna_{s}.healthy_repair.h5ad"
        if dest.is_file() and dest.stat().st_size > 0:
            log(f"  {s}: {dest.name} already present -- skip")
            continue
        build_stage(s, lut)
    log("done")
    return 0

if __name__ == "__main__":
    sys.exit(main())
