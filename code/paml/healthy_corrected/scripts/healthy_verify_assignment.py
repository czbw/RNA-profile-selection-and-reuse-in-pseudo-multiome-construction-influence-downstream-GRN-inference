#!/usr/bin/env python
"""Verify the repaired healthy-cell assignment independently of the script that produced it."""
from __future__ import annotations

import argparse
import glob
import math
import sys
from pathlib import Path

import pandas as pd

ATAC_MANIFEST = ("preprocessing", "out")
ALL_STAGES = ["Diagnosis", "Remission", "Relapse"]
DID_STAGES = ["Diagnosis", "Relapse"]
HEALTHY_LABEL = "Normal"
CONF_GATE = 0.5
BUILDER_RNA_POOL = {"Diagnosis": 18_910, "Remission": 78_432, "Relapse": 35_170}

def log(m: str) -> None:
    print(f"[h02] {m}", flush=True)

def atac_unit(sample_id: pd.Series) -> pd.Series:
    return sample_id.astype(str).str.strip().str.split("_", n=1).str[1]

def consumer_keys(root: Path) -> set:
    """The donor keys build_paired_rna.py can resolve, read from the h5ads themselves."""
    import anndata as ad
    keys = set()
    for tp in ALL_STAGES:
        fp = root / "data/paml" / f"healthy_{tp.lower()}.h5ad"
        if not fp.is_file():
            raise SystemExit(f"[h02] MISSING consumer RNA file {fp}")
        A = ad.read_h5ad(fp, backed="r")
        obs = A.obs
        if int(obs.shape[0]) != BUILDER_RNA_POOL[tp]:
            raise SystemExit(f"[h02] {fp.name} holds {obs.shape[0]:,}, expected "
                             f"{BUILDER_RNA_POOL[tp]:,}")
        bc = pd.Index(obs.index).str.split("|").str[-1]
        keys |= set(obs["Library_ID"].astype(str).values + "||" + bc.values)
        A.file.close()
    log(f"consumer (healthy h5ads): {len(keys):,} resolvable (library, barcode) keys")
    return keys

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument("--include-remission", action="store_true")
    ap.add_argument("--include-aml7-dx", action="store_true")
    a = ap.parse_args()
    root = Path(a.root)
    stages = ALL_STAGES if a.include_remission else DID_STAGES
    exclude = [] if a.include_aml7_dx else ["AML7_DX"]
    outdir = root / "paml/healthy_corrected" / "out"
    fails: list[str] = []

    # ---- RNA pool, re-read from source, independent of p_h01 ------------------------------
    files = sorted(glob.glob(str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")))
    if not files:
        raise SystemExit("[h02] no RNA metadata found")
    R = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in files], ignore_index=True)
    R = R[R["Malignant"].astype(str).str.strip() == HEALTHY_LABEL].copy()
    R["_pid"] = R["Patient_ID"].str.strip()
    R["_tp"] = R["Patient_Sample"].str.strip()
    R["_ct"] = R["Classified_Celltype"].str.strip()
    R["_unit"] = R["Library_ID"].str.strip()
    counts = R["_tp"].value_counts().to_dict()
    for tp in ALL_STAGES:
        if int(counts.get(tp, 0)) != BUILDER_RNA_POOL[tp]:
            raise SystemExit(f"[h02] healthy RNA at {tp} is {counts.get(tp,0):,}, expected "
                             f"{BUILDER_RNA_POOL[tp]:,}")
    if exclude:
        R = R[~R["_unit"].isin(exclude)].copy()
    R["_key"] = R["_unit"] + "||" + R["Cell_Barcode"].str.strip()
    meta = R.set_index("_key")[["_pid", "_tp", "_ct"]].to_dict("index")
    pool_n = R.groupby(["_pid", "_tp", "_ct"]).size()
    log(f"RNA deposit ({HEALTHY_LABEL}): {len(R):,} cells, {R['_key'].nunique():,} unique keys")

    ckeys = consumer_keys(root)
    # the deposit-derived pool must be a subset of what the builder can resolve, or p_h01 could
    # pick a donor that build_paired_rna.py cannot supply
    stray = set(R["_key"]) - ckeys
    if stray:
        fails.append(f"{len(stray):,} deposit keys are absent from the healthy h5ads -- the "
                     f"producer's pool is larger than the consumer's")
    else:
        log("producer pool is a subset of the consumer's key set: OK")

    for stage in stages:
        rep = outdir / f"cell_data_{stage}.healthy_repaired.csv"
        orig = root.joinpath(*ATAC_MANIFEST) / f"cell_data_{stage}.csv"
        if not rep.is_file():
            fails.append(f"{stage}: MISSING {rep.name}")
            continue
        D = pd.read_csv(rep, index_col=0, dtype=str)
        O = pd.read_csv(orig, index_col=0, dtype=str)
        log("")
        log(f"--- {stage}: {len(D):,} repaired rows (original {len(O):,}) ---")

        # (6) drop-in shape
        if list(D.columns) != list(O.columns):
            fails.append(f"{stage}: column set/order differs from the manifest -- not a drop-in "
                         f"(got {list(D.columns)})")
        else:
            log(f"  columns identical to the manifest: {len(D.columns)} cols OK")

        D["_pid"] = D["Patient_ID"].str.strip()
        D["_tp"] = D["timepoint"].str.strip()
        D["_ct"] = D["Classified_Celltype"].str.strip()
        D["_unit"] = atac_unit(D["sample_id"])
        D["_key"] = (D["Linked_scRNA_Sample"].str.strip() + "||"
                     + D["Linked_scRNA_barcode"].str.strip())

        # (1) complete, single-valued
        if D["Linked_scRNA_barcode"].isna().any():
            fails.append(f"{stage}: unassigned cells present")
        if D.index.duplicated().any():
            fails.append(f"{stage}: duplicate ATAC cell ids")
        if exclude and D["_unit"].isin(exclude).any():
            fails.append(f"{stage}: excluded unit {exclude} survived into the output")

        # (7) the placeholder, and what it newly admits
        conf = pd.to_numeric(D["Linked_scRNA_confidence_score"], errors="coerce")
        if not (conf == 1.0).all():
            fails.append(f"{stage}: confidence placeholder is not uniformly 1.0")
        oc = pd.to_numeric(O["Linked_scRNA_confidence_score"], errors="coerce")
        newly = int((oc.reindex(D.index) < CONF_GATE).sum())
        log(f"  placeholder 1.0 on all rows; {newly:,} of them were below the original "
            f"{CONF_GATE} gate and are newly admitted")

        # (3) the join key must resolve in the builder's input files
        unresolved = ~D["_key"].isin(ckeys)
        if unresolved.any():
            fails.append(f"{stage}: {int(unresolved.sum()):,} donor keys absent from the healthy "
                         f"h5ads -- build_paired_rna.py would silently drop them")
        else:
            log(f"  join keys resolve in the consumer: {len(D):,}/{len(D):,} OK")

        # (2) no crossing, re-derived from the deposit rather than from p_h01's maps
        ok = D["_key"].isin(meta.keys())
        if not ok.all():
            fails.append(f"{stage}: {int((~ok).sum()):,} donor keys absent from the deposit pool")
        got = pd.DataFrame([meta[k] for k in D.loc[ok, "_key"]], index=D.index[ok])
        for field, label in (("_pid", "patient"), ("_tp", "timepoint"), ("_ct", "cell type")):
            bad = int((got[field].to_numpy() != D.loc[ok, field].to_numpy()).sum())
            if bad:
                fails.append(f"{stage}: {bad:,} assignments CROSS {label}")
            else:
                log(f"  no {label} crossing: {int(ok.sum()):,} checked")

        # (5) per-stratum reuse never exceeds the independently computed floor
        worst = []
        for key, grp in D.groupby(["_pid", "_tp", "_ct"]):
            n_rna = int(pool_n.get(key, 0))
            if n_rna == 0:
                fails.append(f"{stage}: stratum {key} kept but has no donors")
                continue
            floor = math.ceil(len(grp) / n_rna)
            used = grp["Linked_scRNA_barcode"].value_counts().max()
            if used > floor:
                worst.append((key, int(used), floor))
        if worst:
            fails.append(f"{stage}: {len(worst)} strata exceed their reuse floor, e.g. {worst[:3]}")
        else:
            log("  every stratum at or below its reuse floor")

        # (4) dropped cells are exactly the no-donor ones, or the declared excluded unit
        dropped = O.index.difference(D.index)
        Od = O.loc[dropped].copy()
        Od["_pid"] = Od["Patient_ID"].str.strip()
        Od["_tp"] = Od["timepoint"].str.strip()
        Od["_ct"] = Od["Classified_Celltype"].str.strip()
        Od["_unit"] = atac_unit(Od["sample_id"])
        n_excl = int(Od["_unit"].isin(exclude).sum()) if exclude else 0
        wrongly = 0
        for key, grp in Od[~Od["_unit"].isin(exclude)].groupby(["_pid", "_tp", "_ct"]):
            if int(pool_n.get(key, 0)) > 0:
                wrongly += len(grp)
        if wrongly:
            fails.append(f"{stage}: {wrongly:,} cells dropped despite having donors available")
        else:
            log(f"  {len(dropped):,} dropped ({len(dropped)/len(O):.2%}): "
                f"{n_excl:,} declared-exclusion, {len(dropped)-n_excl:,} no-donor. None silent.")

        donors = D["Linked_scRNA_barcode"].nunique()
        log(f"  achieved reuse {len(D)/donors:.4f}x on {donors:,} donors "
            f"({donors/len(R):.2%} of the healthy pool)")

    log("")
    if fails:
        print("[h02] ===== VERIFICATION FAILED =====", file=sys.stderr)
        for f in fails:
            print(f"[h02]   {f}", file=sys.stderr)
        return 1
    log("===== ALL CHECKS PASSED - safe to upload =====")
    return 0

if __name__ == "__main__":
    sys.exit(main())
