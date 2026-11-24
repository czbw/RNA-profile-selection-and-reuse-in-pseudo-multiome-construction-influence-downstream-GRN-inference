#!/usr/bin/env python
"""Feasibility census for the corrected pairing of healthy cells: how much donor reuse is forced."""
from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# The healthy ATAC cell manifest that was fed to the ATAC chain.
ATAC_MANIFEST = ("preprocessing", "out")

ALL_STAGES = ["Diagnosis", "Remission", "Relapse"]
DID_STAGES = ["Diagnosis", "Relapse"]          # stages used by the downstream analysis
HEALTHY_LABEL = "Normal"                       # deposit's Malignant value for healthy cells

# Healthy RNA cells per stage in data/paml/healthy_*.h5ad (obs counts). Asserted against the deposit,
# since a mismatch would make donors silently fail to join.
BUILDER_RNA_POOL = {"Diagnosis": 18_910, "Remission": 78_432, "Relapse": 35_170}

ATAC_REQ = ["barcode", "sample_id", "timepoint", "Classified_Celltype", "Patient_ID",
            "Linked_scRNA_barcode", "Linked_scRNA_Sample", "Linked_scRNA_confidence_score"]
RNA_REQ = ["Cell_Barcode", "Library_ID", "Patient_Sample", "Classified_Celltype",
           "Malignant", "Patient_ID"]

# A unit is a deposit library. The ATAC side carries a GSM prefix the RNA side lacks
# ("GSM7498708_AML16_DX" vs "AML16_DX"), so the key is normalised before comparison.
def atac_unit(sample_id: pd.Series) -> pd.Series:
    return sample_id.astype(str).str.strip().str.split("_", n=1).str[1]

def log(m: str) -> None:
    print(f"[h00] {m}", flush=True)

def load_atac(root: Path, stages: list[str]) -> pd.DataFrame:
    frames = []
    for tp in stages:
        fp = root.joinpath(*ATAC_MANIFEST) / f"cell_data_{tp}.csv"
        if not fp.is_file():
            raise SystemExit(f"[h00] MISSING {fp}")
        d = pd.read_csv(fp, index_col=0, dtype=str)
        missing = [c for c in ATAC_REQ if c not in d.columns]
        if missing:
            raise SystemExit(f"[h00] {fp.name}: missing columns {missing}; found {list(d.columns)}")
        seen = set(d["timepoint"].astype(str).str.strip())
        if seen != {tp}:
            raise SystemExit(f"[h00] {fp.name}: expected timepoint {{{tp}}}, found {sorted(seen)}")
        d["_src"] = fp.name
        frames.append(d)
        log(f"ATAC {fp.name}: {len(d):,} cells")
    return pd.concat(frames, ignore_index=True)

def load_rna(root: Path, stages: list[str]) -> pd.DataFrame:
    """Healthy RNA pool from the deposit, checked against the downstream h5ad counts."""
    pat = str(root / "data/paml/GSE235063_RAW" / "*processed_metadata.tsv.gz")
    files = sorted(glob.glob(pat))
    if not files:
        raise SystemExit(f"[h00] no RNA metadata matched {pat}")
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t", dtype=str, compression="gzip")
        missing = [c for c in RNA_REQ if c not in d.columns]
        if missing:
            raise SystemExit(f"[h00] {Path(f).name}: missing columns {missing}")
        frames.append(d)
    R = pd.concat(frames, ignore_index=True)
    lab = R["Malignant"].astype(str).str.strip()
    unknown = sorted(set(lab) - {"Malignant", HEALTHY_LABEL})
    if unknown:
        raise SystemExit(f"[h00] unexpected Malignant labels {unknown} -- refusing to guess")
    H = R[lab == HEALTHY_LABEL].copy()
    log(f"RNA deposit: {len(R):,} cells across {len(files)} libraries "
        f"({len(H):,} {HEALTHY_LABEL}, {len(H)/len(R):.1%})")

    H["_tp"] = H["Patient_Sample"].astype(str).str.strip()
    seen = set(H["_tp"])
    if not seen <= set(ALL_STAGES):
        raise SystemExit(f"[h00] unexpected RNA timepoints {sorted(seen - set(ALL_STAGES))}")

    # ---- consistency with the downstream h5ads -------------------------------------------
    # The healthy h5ads are the pool build_paired_rna.py can supply; assert the deposit-derived
    # pool has the same counts per timepoint.
    counts = H["_tp"].value_counts().to_dict()
    for tp in ALL_STAGES:
        got, want = int(counts.get(tp, 0)), BUILDER_RNA_POOL[tp]
        if got != want:
            raise SystemExit(
                f"[h00] {HEALTHY_LABEL} RNA at {tp} is {got:,} but the consumer's "
                f"data/paml/healthy_{tp.lower()}.h5ad holds {want:,}. These MUST match or donors will "
                f"silently fail to join. Investigate -- do not adjust the constant to make it pass.")
    log(f"pool locked to the consumer: " +
        ", ".join(f"{t} {BUILDER_RNA_POOL[t]:,}" for t in ALL_STAGES))
    # Returns the full healthy pool; the caller restricts to the scope's timepoints after the
    # cross-timepoint diagnostics, which need every timepoint.
    return H

def verify_h5ad_pool(root: Path) -> dict:
    """Optional: re-read the h5ad obs and compare key sets (counts can agree while membership differs).

    Skipped with a warning if anndata is unavailable.
    """
    try:
        import anndata as ad
    except Exception as e:                                    # pragma: no cover
        log(f"WARNING: anndata unavailable ({e}); h5ad key-set check SKIPPED (counts still asserted)")
        return {"checked": False, "reason": str(e)}
    out = {"checked": True, "per_stage": {}}
    for tp in ALL_STAGES:
        fp = root / "data/paml" / f"healthy_{tp.lower()}.h5ad"
        if not fp.is_file():
            raise SystemExit(f"[h00] MISSING consumer RNA file {fp}")
        A = ad.read_h5ad(fp, backed="r")
        obs = A.obs
        # obs_names are "<stem>|<barcode>"; the join key is (Library_ID, Cell_Barcode)
        bc = pd.Index(obs.index).str.split("|").str[-1]
        keys = set(obs["Library_ID"].astype(str).values + "||" + bc.values)
        A.file.close()
        out["per_stage"][tp] = {"n_obs": int(obs.shape[0]), "n_keys": len(keys)}
        if int(obs.shape[0]) != BUILDER_RNA_POOL[tp]:
            raise SystemExit(f"[h00] {fp.name} holds {obs.shape[0]:,}, expected "
                             f"{BUILDER_RNA_POOL[tp]:,}")
        out["per_stage"][tp]["keys"] = keys
    return out

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument("--include-remission", action="store_true",
                    help="census Remission too. It is a real healthy timepoint but is not used "
                         "by the downstream difference-in-differences analysis.")
    ap.add_argument("--include-aml7-dx", action="store_true",
                    help="override the AML7_DX exclusion. The exclusion is re-derived below from "
                         "the healthy data, not copied from the malignant run.")
    a = ap.parse_args()
    root = Path(a.root)
    stages = ALL_STAGES if a.include_remission else DID_STAGES
    log(f"scope: healthy {' + '.join(stages)}"
        + ("" if a.include_remission else "   (Remission censused only with --include-remission)"))

    A = load_atac(root, stages)
    RALL = load_rna(root, stages)
    h5 = verify_h5ad_pool(root)

    # normalise join keys
    A["_tp"] = A["timepoint"].astype(str).str.strip()
    for df in (A, RALL):
        df["_pid"] = df["Patient_ID"].astype(str).str.strip()
        df["_ct"] = df["Classified_Celltype"].astype(str).str.strip()
    RALL["_unit"] = RALL["Library_ID"].astype(str).str.strip()
    A["_unit"] = atac_unit(A["sample_id"])
    # `R` is the donor pool (healthy RNA at the scope's timepoints); `RALL` keeps every healthy
    # timepoint, for diagnostics that must resolve cross-timepoint links.
    R = RALL[RALL["_tp"].isin(stages)].copy()

    if A["_unit"].isna().any():
        raise SystemExit("[h00] ATAC sample_id did not parse to a library unit -- check the "
                         "'GSM<n>_<LIB>' convention")
    orphan_units = sorted(set(A["_unit"]) - set(RALL["_unit"]))
    if orphan_units:
        raise SystemExit(f"[h00] ATAC units absent from the {HEALTHY_LABEL} RNA deposit entirely: "
                         f"{orphan_units} -- a patient present in one modality only. STOP.")
    orphan_pid = sorted(set(A["_pid"]) - set(RALL["_pid"]))
    if orphan_pid:
        raise SystemExit(f"[h00] ATAC patients absent from the {HEALTHY_LABEL} RNA: {orphan_pid}")

    # Detect libraries whose name token (DX/REM/REL) disagrees with Patient_Sample.
    tok2tp = {"DX": "Diagnosis", "REM": "Remission", "REL": "Relapse"}
    RALL["_tok_tp"] = RALL["_unit"].str.rsplit("_", n=1).str[1].str.upper().map(tok2tp)
    swapped = RALL[RALL["_tok_tp"].notna() & (RALL["_tok_tp"] != RALL["_tp"])]
    swap_units = {}
    if len(swapped):
        for u, g in swapped.groupby("_unit"):
            swap_units[str(u)] = {"library_token_says": str(g["_tok_tp"].iloc[0]),
                                  "Patient_Sample_says": sorted(set(g["_tp"])),
                                  "n_healthy_rna": int(len(g)),
                                  "patients": sorted(set(g["_pid"]))}
        log("")
        log("DEPOSIT LABEL SWAP detected (library name vs Patient_Sample disagree):")
        for u, v in swap_units.items():
            log(f"   {u:<12} token says {v['library_token_says']:<10} "
                f"Patient_Sample says {','.join(v['Patient_Sample_says']):<10} "
                f"{v['n_healthy_rna']:,} healthy RNA  patients {v['patients']}")
        log("   -> Patient_Sample is treated as authoritative, matching what the malignant arm")
        log("      already did. No library token is ever used as a timepoint below.")

    # Every (patient x timepoint) with ATAC cells must have some healthy RNA at that same
    # (patient x timepoint), or the whole unit would be dropped.
    a_keys = set(map(tuple, A[["_pid", "_tp"]].drop_duplicates().to_numpy()))
    r_keys = set(map(tuple, R[["_pid", "_tp"]].drop_duplicates().to_numpy()))
    missing_pt = sorted(a_keys - r_keys)
    if missing_pt:
        raise SystemExit(
            f"[h00] {len(missing_pt)} (patient x timepoint) unit(s) have healthy ATAC cells but "
            f"ZERO healthy RNA at the same timepoint: {missing_pt}. Every cell in them would be "
            f"dropped. STOP and decide explicitly rather than letting the assignment silently "
            f"discard a unit.")
    log(f"every (patient x timepoint) ATAC unit has healthy RNA at its own timepoint: "
        f"{len(a_keys)} units checked")

    # ---- deposit scan: derive the exclusion for healthy cells independently ----------------
    # p_a01 excluded AML7_DX on the malignant side (2,591 ATAC vs a 36-cell RNA library); whether
    # the healthy compartment has a comparable outlier is checked separately here.
    ca = A.groupby("_unit").size().rename("n_atac")
    cr = R.groupby("_unit").size().rename("n_rna")
    U = pd.concat([ca, cr], axis=1).fillna(0).astype(int)
    U = U[U["n_atac"] > 0].copy()
    U["unit_ratio"] = U["n_atac"] / U["n_rna"].replace(0, np.nan)
    U = U.sort_values("unit_ratio", ascending=False)
    log("")
    log("unit-level deposit scan (healthy ATAC cells per available healthy RNA cell):")
    for unit, row in U.head(8).iterrows():
        log(f"   {unit:<12} ATAC {row['n_atac']:>6,}  RNA {row['n_rna']:>6,}  "
            f"{row['unit_ratio']:>7.2f}x")
    top = U.iloc[0]
    runner_up = U.iloc[1]["unit_ratio"] if len(U) > 1 else float("nan")
    log(f"   worst unit is {U.index[0]} at {top['unit_ratio']:.2f}x, "
        f"{top['unit_ratio']/runner_up:.1f}x worse than the next ({runner_up:.2f}x)")

    excl = [] if a.include_aml7_dx else ["AML7_DX"]
    if excl:
        if excl[0] not in U.index:
            raise SystemExit(f"[h00] cannot evaluate the {excl[0]} exclusion: unit not in this "
                             f"scope. Re-check before excluding blindly.")
        log(f"   -> EXCLUDING {excl}: an independent healthy-side outlier, not a copied decision")
        nA0, nR0 = len(A), len(R)
        A = A[~A["_unit"].isin(excl)].copy()
        R = R[~R["_unit"].isin(excl)].copy()
        RALL = RALL[~RALL["_unit"].isin(excl)].copy()
        log(f"   dropped {nA0-len(A):,} ATAC, {nR0-len(R):,} RNA")

    # ---- per-stratum census --------------------------------------------------------------
    na = A.groupby(["_pid", "_tp", "_ct"]).size().rename("n_atac")
    nr = R.groupby(["_pid", "_tp", "_ct"]).size().rename("n_rna")
    C = pd.concat([na, nr], axis=1).fillna(0).astype(int)
    C = C[C["n_atac"] > 0].copy()
    C["min_reuse"] = [(math.inf if r == 0 else max(1.0, math.ceil(x / r)))
                      for x, r in zip(C["n_atac"], C["n_rna"])]
    C["short"] = C["n_rna"] < C["n_atac"]
    C["unassignable"] = C["n_rna"] == 0

    assignable = C[~C["unassignable"]]
    donors_needed = int(np.minimum(assignable["n_atac"], assignable["n_rna"]).sum())
    if donors_needed == 0:
        raise SystemExit("[h00] no assignable stratum at all -- the pool definition is wrong")
    global_floor = float(assignable["n_atac"].sum() / donors_needed)

    atac_total = int(C["n_atac"].sum())
    atac_short = int(C.loc[C["short"], "n_atac"].sum())
    atac_unassignable = int(C.loc[C["unassignable"], "n_atac"].sum())

    log("")
    log("================= HEALTHY FEASIBILITY CENSUS =================")
    log(f"scope                             : healthy {' + '.join(stages)}")
    log(f"strata                            : {len(C):,}")
    log(f"strata with n_RNA < n_ATAC        : {int(C['short'].sum()):,} "
        f"({C['short'].mean():.1%})")
    log(f"ATAC cells in short strata        : {atac_short:,} / {atac_total:,} "
        f"({atac_short/atac_total:.2%})")
    log(f"ATAC cells UNASSIGNABLE (n_RNA=0) : {atac_unassignable:,} "
        f"({atac_unassignable/atac_total:.2%})")
    log(f"MINIMUM FORCED GLOBAL REUSE       : {global_floor:.3f}x   "
        f"(healthy-specific; the malignant floor is NOT reused)")

    # ---- what the original healthy link actually did --------------------------------------
    # Use RALL, not R: a Diagnosis ATAC cell may link into a Remission transcriptome, and that link
    # must resolve to be counted as a cross-timepoint link.
    tp_of_unit = dict(zip(RALL["_unit"], RALL["_tp"]))
    link_donors = A["Linked_scRNA_barcode"].astype(str).nunique()
    link_reuse = len(A) / link_donors if link_donors else float("nan")
    log("")
    log(f"original healthy link, same cells : {link_reuse:.2f}x reuse")
    log(f"  distinct RNA donors used        : {link_donors:,} of {len(R):,} healthy RNA "
        f"({link_donors/len(R):.2%})")
    log(f"  attributable to SCARCITY        : {global_floor:.3f}x")
    if link_reuse > 1 and global_floor > 1:
        log(f"  attributable to the RULE        : {link_reuse/global_floor:.2f}x on top "
            f"({1 - (global_floor-1)/(link_reuse-1):.1%} of the excess)")

    # Cross-timepoint links: build_paired_rna.py searches all three healthy RNA h5ads, so a
    # Diagnosis ATAC cell linking to a Remission transcriptome is possible.
    linked_tp = A["Linked_scRNA_Sample"].astype(str).str.strip().map(tp_of_unit)
    ok = linked_tp.notna()
    per_stage_link = {}
    for tp in stages:
        m = (A["_tp"] == tp) & ok
        sub = A.loc[m]
        if not len(sub):
            raise SystemExit(f"[h00] no resolvable links at {tp} -- unexpected")
        d = sub["Linked_scRNA_barcode"].nunique()
        per_stage_link[tp] = {
            "n_atac": int((A["_tp"] == tp).sum()),
            "link_reuse": float((A["_tp"] == tp).sum() / d),
            "link_distinct_donors": int(d),
            "pct_links_within_timepoint": float((linked_tp[m] == tp).mean()),
            "pct_links_resolvable_in_healthy_pool": float(m.sum() / max(1, (A["_tp"] == tp).sum())),
        }
    same_tp = float((linked_tp[ok] == A.loc[ok, "_tp"]).mean())
    log("")
    log(f"links staying WITHIN timepoint    : {same_tp:.2%} "
        f"(resolvable for {int(ok.sum()):,} of {len(A):,})")
    for tp in stages:
        v = per_stage_link[tp]
        log(f"   {tp:<10} within-timepoint {v['pct_links_within_timepoint']:.2%} | "
            f"reuse {v['link_reuse']:.2f}x on {v['link_distinct_donors']:,} donors | "
            f"resolvable {v['pct_links_resolvable_in_healthy_pool']:.2%}")
    log("   -> a hard within-timepoint constraint is the single largest correction, and it is a")
    log("      LARGER correction here than on the malignant arm (which was 83.95% within).")

    # ---- the confidence-score gate: an effective filter here, unlike the malignant arm ----
    conf = pd.to_numeric(A["Linked_scRNA_confidence_score"], errors="coerce")
    if conf.isna().any():
        raise SystemExit("[h00] non-numeric Linked_scRNA_confidence_score present")
    below = float((conf < 0.5).mean())
    log("")
    log(f"Linked_scRNA_confidence_score     : median {conf.median():.3f}, min {conf.min():.3f}, "
        f"{below:.2%} below the 0.5 gate")
    log("   -> this was measured as a NO-OP on the malignant arm (median 1.000, 0% below 0.5).")
    log("      On the healthy arm it BITES. Replacing it with the 1.0 placeholder therefore")
    log("      CHANGES the admitted cell set here; p_h01 quantifies the change explicitly.")

    worst = C.sort_values("min_reuse", ascending=False).head(12)
    log("")
    log("worst strata by forced reuse:")
    for (pid, tp, ct), row in worst.iterrows():
        mr = "INF (no RNA)" if not np.isfinite(row["min_reuse"]) else f"{row['min_reuse']:.2f}x"
        log(f"   {pid:>10} {tp:<10} {ct:<20} ATAC {row['n_atac']:>6,}  "
            f"RNA {row['n_rna']:>6,}  floor {mr}")

    out = root / "paml/healthy_corrected" / "out"
    out.mkdir(parents=True, exist_ok=True)
    # The filename carries the scope, so a three-stage run (--include-remission) does not
    # overwrite the two-stage census.
    sfx = "" if not a.include_remission else ".withRemission"
    C.reset_index().rename(columns={"_pid": "patient", "_tp": "timepoint", "_ct": "cell_type"}) \
     .to_csv(out / f"h00_stratum_census{sfx}.csv", index=False)
    U.reset_index().rename(columns={"index": "unit"}).to_csv(
        out / f"h00_unit_census{sfx}.csv", index=False)

    rep = {
        "scope": f"healthy {' + '.join(stages)}",
        "stages": stages,
        "rna_pool_policy": f"{HEALTHY_LABEL} ONLY -- forced by the consumer, not chosen",
        "builder_rna_pool_per_stage": {t: BUILDER_RNA_POOL[t] for t in stages},
        "h5ad_key_set_checked": bool(h5.get("checked")),
        "excluded_units": excl,
        "exclusion_evidence": {
            "worst_unit": str(U.index[0]),
            "worst_unit_ratio": float(U.iloc[0]["unit_ratio"]),
            "worst_unit_n_atac": int(U.iloc[0]["n_atac"]),
            "worst_unit_n_rna": int(U.iloc[0]["n_rna"]),
            "runner_up_ratio": float(runner_up),
        },
        "n_strata": int(len(C)), "n_strata_short": int(C["short"].sum()),
        "atac_total": atac_total, "atac_in_short_strata": atac_short,
        "pct_atac_in_short_strata": atac_short / atac_total,
        "atac_unassignable": atac_unassignable,
        "pct_atac_unassignable": atac_unassignable / atac_total,
        "minimum_forced_global_reuse": global_floor,
        "original_link_reuse": link_reuse,
        "original_link_distinct_donors": int(link_donors),
        "healthy_rna_pool": int(len(R)),
        "pct_rna_pool_used_by_link": link_donors / len(R),
        "pct_links_within_timepoint": same_tp,
        "per_stage_link": per_stage_link,
        "confidence_median": float(conf.median()),
        "confidence_min": float(conf.min()),
        "pct_confidence_below_0p5": below,
        "n_patients": int(A["_pid"].nunique()),
        "n_cell_types": int(A["_ct"].nunique()),
        "atac_per_stage": {t: int((A["_tp"] == t).sum()) for t in stages},
    }
    (out / f"h00_healthy_feasibility_report{sfx}.json").write_text(json.dumps(rep, indent=2))
    log("")
    log(f"wrote h00_stratum_census{sfx}.csv, h00_unit_census{sfx}.csv and "
        f"h00_healthy_feasibility_report{sfx}.json  ->  {out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
