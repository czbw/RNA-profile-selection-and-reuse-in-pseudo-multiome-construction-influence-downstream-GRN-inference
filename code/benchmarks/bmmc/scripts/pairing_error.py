#!/usr/bin/env python
"""The donor axis: what the BMMC pairings actually did, and what that did to the eGRN."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent

STAGE_DEFAULT = "bmmc"
SEED_DEFAULT = 555
TRUTH = "truth"
ARMS_DEFAULT = ["truth", "rand_celltype_1to1", "rand_celltype_x_donor_1to1", "rand_global"]

# arm token -> extra file stems to try, in order, after the canonical one (the two donor-arm spellings are interchangeable).
ARM_TABLE_ALIASES = {
    "rand_celltype_x_donor_1to1": ["rand_celltype_donor_1to1"],
    "rand_celltype_donor_1to1": ["rand_celltype_x_donor_1to1"],
}

DONOR_BLIND = "rand_celltype_1to1"
DONOR_AWARE = "rand_celltype_x_donor_1to1"

EXPECTED = {
    "truth":            {"self_hit": (1.000000, 5e-7)},
    "rand_celltype_1to1": {"same_cell_type": (1.0000, 5e-5),
                           "same_donor": (0.1974, 5e-5),
                           "self_hit": (0.000332, 5e-7)},
    "rand_celltype_x_donor_1to1": {"same_cell_type": (1.0000, 5e-5),
                                   "same_donor": (1.0000, 5e-5),
                                   "self_hit": (0.002700, 5e-7)},
    "rand_global":      {"same_cell_type": (0.0902, 5e-5),
                         "same_donor": (0.1424, 5e-5),
                         "self_hit": (0.000029, 5e-7)},
}
EXPECTED_N_CELLS = 69249          # every arm is an exact bijection over exactly this cell set

# The rungs pulled out of b_04's report for the donor contrast. (level, label, path-in-arm-entry).
# Chosen to span the whole ladder: raw edge evidence (L1) up to "which TFs survive at all" (L5).
GRN_STATS = [
    ("L1", "L1_region_to_gene", "spearman_importance"),
    ("L1", "L1_region_to_gene", "spearman_rho"),
    ("L1", "L1_region_to_gene", "pct_rho_undefined"),
    ("L1", "L1_region_to_gene", "topN_jaccard_vs_truth"),
    ("L1", "L1_region_to_gene", "median_abs_distance_topN"),
    ("L2", "L2_tf_to_gene", "n_TF"),
    ("L2", "L2_tf_to_gene", "full_edge_set.recall"),
    ("L2", "L2_tf_to_gene", "full_edge_set.precision"),
    ("L2", "L2_tf_to_gene", "full_edge_set.jaccard"),
    ("L2", "L2_tf_to_gene", "topN.jaccard"),
    ("L2", "L2_tf_to_gene", "topN_fold_enrichment"),
    ("L2", "L2_tf_to_gene", "spearman_importance_shared"),
    ("L2", "L2_tf_to_gene", "spearman_rho_shared"),
    ("L3", "L3_L4_L5_eRegulons", "L3_target_sets.median_jaccard_genes"),
    ("L3", "L3_L4_L5_eRegulons", "L3_region_sets.median_jaccard_regions"),
    ("L4", "L3_L4_L5_eRegulons", "L4_eRegulon_recovery.recall"),
    ("L4", "L3_L4_L5_eRegulons", "L4_eRegulon_recovery.precision"),
    ("L4", "L3_L4_L5_eRegulons", "L4_eRegulon_recovery.jaccard"),
    ("L5", "L3_L4_L5_eRegulons", "L5_TF_recovery.recall"),
    ("L5", "L3_L4_L5_eRegulons", "L5_TF_recovery.precision"),
    ("L5", "L3_L4_L5_eRegulons", "L5_TF_recovery.jaccard"),
]

def log(m):
    print(f"[b_05] {m}", flush=True)

def die(msg):
    print(f"[b_05] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def die_missing(p: Path, why: str):
    die(f"required input missing: {p}\n       reason: {why}")

def dig(d, dotted: str):
    """Walk a dotted path into nested dicts; None if any hop is absent or not a dict."""
    cur = d
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur

def resolve_table(uploads: Path, stage: str, seed: int, arm: str) -> Path:
    """Canonical name first, then the declared aliases. Fails with every path it tried."""
    stems = [arm] + ARM_TABLE_ALIASES.get(arm, [])
    tried = []
    for stem in stems:
        for name in (f"arm_{stem}_{stage}_seed{seed}_full.csv", f"arm_{stem}_seed{seed}.csv"):
            p = uploads / name
            tried.append(p)
            if p.is_file() and p.stat().st_size > 0:
                if stem != arm:
                    log(f"arm {arm!r}: resolved to ALIAS file {p.name} "
                        f"(builder name {stem!r}; see the arm-token header block)")
                return p
    die("no assignment table for arm {!r}. Tried:\n       {}".format(
        arm, "\n       ".join(str(t) for t in tried)))
    raise AssertionError("unreachable")

def score_arm(name: str, tab: pd.DataFrame, meta: pd.DataFrame) -> dict:
    """Every rate here comes from the join, never from the table's own bookkeeping columns."""
    need = {"atac_barcode", "assigned_rna_barcode"}
    missing = need - set(tab.columns)
    if missing:
        die(f"arm {name!r}: assignment table lacks column(s) {sorted(missing)}; "
            f"the join key is `atac_barcode` -> `assigned_rna_barcode` on the h5ad barcode")

    n = len(tab)
    atac = tab["atac_barcode"].to_numpy()
    rna = tab["assigned_rna_barcode"].to_numpy()

    # --- bijection, asserted as its three separate facts. "It is a permutation" is not the same
    # --- claim as "every ATAC cell appears once" and must not be checked by a proxy for it.
    n_atac_u, n_rna_u = int(pd.unique(atac).size), int(pd.unique(rna).size)
    atac_set, rna_set = set(atac.tolist()), set(rna.tolist())
    bij = {
        "n_rows": n,
        "atac_barcodes_unique": n_atac_u == n,
        "assigned_rna_barcodes_unique": n_rna_u == n,
        "image_equals_domain": atac_set == rna_set,
        "n_atac_unique": n_atac_u,
        "n_assigned_rna_unique": n_rna_u,
        "n_rna_not_in_atac_set": len(rna_set - atac_set),
        "cells_covered_equals_expected": n == EXPECTED_N_CELLS,
        "expected_n_cells": EXPECTED_N_CELLS,
    }
    bij["is_exact_bijection"] = bool(
        bij["atac_barcodes_unique"] and bij["assigned_rna_barcodes_unique"]
        and bij["image_equals_domain"]
    )

    # --- both sides joined to the manifest --------------------------------------------------
    a_side = meta.reindex(atac)
    r_side = meta.reindex(rna)
    n_unjoined_atac = int(a_side["cell_type"].isna().sum())
    n_unjoined_rna = int(r_side["cell_type"].isna().sum())
    if n_unjoined_atac or n_unjoined_rna:
        die(f"arm {name!r}: {n_unjoined_atac} ATAC and {n_unjoined_rna} assigned-RNA barcodes did "
            f"not join to b00_cell_metadata.csv on `h5ad_barcode`. A partial join silently changes "
            f"the denominator of every rate below, so this is fatal rather than reported")

    rates = {"self_hit": float((atac == rna).mean())}
    for key, col in (("same_cell_type", "cell_type"), ("same_donor", "donor"),
                     ("same_site", "site"), ("same_batch", "batch")):
        rates[key] = float((a_side[col].to_numpy() == r_side[col].to_numpy()).mean())

    # --- corroboration only: what the table said about itself ---------------------------------
    declared = {}
    for col, key, same in (("is_self", "self_hit", True),
                           ("cross_donor", "same_donor", False),
                           ("cross_site", "same_site", False)):
        if col not in tab.columns:
            declared[col] = "column absent from this table"
            continue
        v = tab[col]
        b = v.astype(str).str.strip().str.lower().map({"true": True, "false": False})
        if b.isna().any():
            declared[col] = f"unparseable values (not True/False) in {int(b.isna().sum())} rows"
            continue
        rate = float(b.mean()) if same else float((~b).mean())
        declared[col] = {
            "declared_rate": rate,
            "join_derived_rate": rates[key],
            "agrees": bool(abs(rate - rates[key]) < 1e-12),
        }

    return {"bijection": bij, "rates_from_join": rates, "declared_columns_corroboration": declared}

def check_expected(name: str, rates: dict, bij: dict) -> list[dict]:
    out = []
    for stat, (want, tol) in EXPECTED.get(name, {}).items():
        got = rates[stat]
        out.append({
            "arm": name, "statistic": stat, "expected": want, "tolerance": tol,
            "measured": got, "abs_diff": abs(got - want),
            "verdict": "PASS" if abs(got - want) <= tol else "FAIL",
        })
    out.append({
        "arm": name, "statistic": "exact_bijection_over_69249_cells",
        "expected": True, "tolerance": 0,
        "measured": bool(bij["is_exact_bijection"] and bij["cells_covered_equals_expected"]),
        "abs_diff": None,
        "verdict": "PASS" if (bij["is_exact_bijection"]
                              and bij["cells_covered_equals_expected"]) else "FAIL",
    })
    return out

def cohort_block(meta: pd.DataFrame) -> dict:
    per_donor_sites = meta.groupby("donor")["site"].nunique()
    per_donor_batches = meta.groupby("donor")["batch"].nunique()
    per_batch_donors = meta.groupby("batch")["donor"].nunique()
    return {
        "n_cells": int(len(meta)),
        "n_batches": int(meta["batch"].nunique()),
        "n_donors": int(meta["donor"].nunique()),
        "n_sites": int(meta["site"].nunique()),
        "n_cell_types": int(meta["cell_type"].nunique()),
        # The nesting is why the donor contrast is not a pure-donor contrast. Measured, not assumed.
        "donors_spanning_more_than_one_site": sorted(
            per_donor_sites[per_donor_sites > 1].index.tolist()),
        "n_donors_spanning_more_than_one_site": int((per_donor_sites > 1).sum()),
        "donors_spanning_more_than_one_batch": sorted(
            per_donor_batches[per_donor_batches > 1].index.tolist()),
        "every_batch_has_exactly_one_donor": bool((per_batch_donors == 1).all()),
        "cells_per_donor": {k: int(v) for k, v in meta["donor"].value_counts().sort_index().items()},
        "cells_per_site": {k: int(v) for k, v in meta["site"].value_counts().sort_index().items()},
    }

def donor_axis_block(arms: dict, trackb: dict | None, trackb_path: Path) -> dict:
    blind = arms.get(DONOR_BLIND)
    aware = arms.get(DONOR_AWARE)
    if blind is None or aware is None:
        return {"status": "NOT COMPUTABLE -- both donor-axis arms must be scored",
                "donor_blind_arm": DONOR_BLIND, "donor_aware_arm": DONOR_AWARE,
                "present": sorted(arms)}

    rb, ra = blind["rates_from_join"], aware["rates_from_join"]

    # --- is this actually a one-variable contrast? Measured, then used as a check ---------------
    checks = {
        "both_exact_bijections": bool(blind["bijection"]["is_exact_bijection"]
                                      and aware["bijection"]["is_exact_bijection"]),
        "same_cell_set_size": blind["bijection"]["n_rows"] == aware["bijection"]["n_rows"],
        "both_same_cell_type_is_1": bool(abs(rb["same_cell_type"] - 1.0) < 5e-5
                                         and abs(ra["same_cell_type"] - 1.0) < 5e-5),
        "donor_concordance_differs": bool(abs(ra["same_donor"] - rb["same_donor"]) > 0.1),
        "aware_arm_is_fully_within_donor": bool(abs(ra["same_donor"] - 1.0) < 5e-5),
    }
    valid = all(checks.values())

    # --- the confound, stated up front rather than discovered later ---------------------------
    confound = {
        "why": ("site and batch are nested under donor for most BMMC donors, so constraining the "
                "bijection to stay within donor also raises SITE and BATCH concordance. The "
                "contrast is therefore donor+site-composite, not pure-donor."),
        "same_site": {"donor_blind": rb["same_site"], "donor_aware": ra["same_site"],
                      "delta_aware_minus_blind": ra["same_site"] - rb["same_site"]},
        "same_batch": {"donor_blind": rb["same_batch"], "donor_aware": ra["same_batch"],
                       "delta_aware_minus_blind": ra["same_batch"] - rb["same_batch"]},
    }

    assignment_side = {
        k: {"donor_blind": rb[k], "donor_aware": ra[k],
            "delta_aware_minus_blind": ra[k] - rb[k]}
        for k in ("self_hit", "same_cell_type", "same_donor", "same_site", "same_batch")
    }

    # --- the GRN side ---------------------------------------------------------------------------
    if trackb is None:
        grn = {
            "status": "PENDING -- b_04 has not been run (or its report is not at the path given)",
            "expected_report": str(trackb_path),
            "produce_it_with": "python network_agreement.py --out-dir ../out --stage bmmc",
        }
    else:
        rows = []
        for level, section, path in GRN_STATS:
            eb = dig(trackb, f"{section}.{DONOR_BLIND}")
            ea = dig(trackb, f"{section}.{DONOR_AWARE}")
            row = {"level": level, "section": section, "statistic": path}
            for side, e in (("donor_blind", eb), ("donor_aware", ea)):
                if e is None:
                    row[side] = None
                    row[f"{side}_state"] = "ARM ABSENT FROM REPORT"
                elif "status" in e:
                    row[side] = None
                    row[f"{side}_state"] = e["status"]
                else:
                    row[side] = dig(e, path)
                    row[f"{side}_state"] = "scored"
            a, b = row["donor_aware"], row["donor_blind"]
            if isinstance(a, (int, float)) and isinstance(b, (int, float)) \
                    and not isinstance(a, bool) and not isinstance(b, bool):
                row["delta_aware_minus_blind"] = a - b
            else:
                row["delta_aware_minus_blind"] = None
                row["delta_reason"] = (
                    "not both numeric -- see the *_state fields. A delta against a missing or "
                    "structurally-absent arm would be an invented number")
            rows.append(row)
        grn = {"status": "scored", "source_report": str(trackb_path),
               "reference_arm": trackb.get("reference_arm"), "rungs": rows}

    return {
        "donor_blind_arm": DONOR_BLIND,
        "donor_aware_arm": DONOR_AWARE,
        "contrast_is_one_variable": valid,
        "contrast_checks": checks,
        "contrast_note": (
            "Both arms are exact bijections within cell type over the same cells; they differ in "
            "whether the bijection is confined to the donor. PBMC 10k CANNOT carry this contrast: "
            "it is a single donor, so every constructed pairing there is trivially within-donor."
            if valid else
            "DELTAS BELOW ARE NOT INTERPRETABLE AS A DONOR EFFECT -- contrast_checks failed; the "
            "two arms differ in more than donor-awareness."),
        "confound_site_and_batch_co_move": confound,
        "assignment_side": assignment_side,
        "grn_side": grn,
    }

def main() -> None:
    ap = argparse.ArgumentParser(
        description="BMMC pairing-error descriptives and the donor-awareness contrast (b_05)")
    ap.add_argument("--stage", default=STAGE_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED_DEFAULT)
    ap.add_argument("--uploads", default=str(BUNDLE / "uploads"),
                    help="directory holding arm_<ARM>_<STAGE>_seed<SEED>_full.csv")
    ap.add_argument("--metadata", default=None,
                    help="b00_cell_metadata.csv (default: uploads/, then out/)")
    ap.add_argument("--arms", nargs="+", default=None,
                    help=f"arms to describe (default: {' '.join(ARMS_DEFAULT)})")
    ap.add_argument("--trackb", default=None,
                    help="b_04's report json (default: out/trackb_rungs_<STAGE>.json). "
                         "Absent -> the GRN half of the donor axis is recorded PENDING")
    ap.add_argument("--out-dir", default=str(BUNDLE / "out"))
    ap.add_argument("--report-name", default=None)
    a = ap.parse_args()

    uploads = Path(a.uploads)
    out = Path(a.out_dir)
    if not uploads.is_dir():
        die(f"--uploads is not a directory: {uploads}")
    out.mkdir(parents=True, exist_ok=True)

    mp = Path(a.metadata) if a.metadata else None
    if mp is None:
        for cand in (uploads / "b00_cell_metadata.csv", out / "b00_cell_metadata.csv"):
            if cand.is_file():
                mp = cand
                break
    if mp is None or not mp.is_file():
        die_missing(uploads / "b00_cell_metadata.csv",
                    "the verified barcode/donor/site/cell_type manifest for all 69,249 cells; "
                    "every rate in this report is computed by joining to it")
    log(f"metadata: {mp}")

    meta = pd.read_csv(mp, dtype=str)
    for col in ("h5ad_barcode", "batch", "donor", "site", "cell_type"):
        if col not in meta.columns:
            die(f"{mp} lacks required column {col!r}")
    if meta["h5ad_barcode"].duplicated().any():
        die(f"{mp}: `h5ad_barcode` is not unique -- it is the join key and a duplicate would "
            f"silently multiply rows on the join")
    meta = meta.set_index("h5ad_barcode")
    log(f"manifest: {len(meta)} cells, {meta['batch'].nunique()} batches, "
        f"{meta['donor'].nunique()} donors, {meta['site'].nunique()} sites, "
        f"{meta['cell_type'].nunique()} cell types")

    arm_names = list(a.arms) if a.arms else list(ARMS_DEFAULT)
    arm_files, arms = {}, {}
    invariants = []
    for name in arm_names:
        p = resolve_table(uploads, a.stage, a.seed, name)
        arm_files[name] = str(p)
        tab = pd.read_csv(p, dtype=str)
        ent = score_arm(name, tab, meta)
        arms[name] = ent
        invariants += check_expected(name, ent["rates_from_join"], ent["bijection"])
        r = ent["rates_from_join"]
        log(f"{name:28s} n={ent['bijection']['n_rows']:6d} "
            f"bijection={ent['bijection']['is_exact_bijection']!s:5s} "
            f"self={r['self_hit']:.6f} cell_type={r['same_cell_type']:.6f} "
            f"donor={r['same_donor']:.6f} site={r['same_site']:.6f} batch={r['same_batch']:.6f}")
        del tab

    failed = [c for c in invariants if c["verdict"] == "FAIL"]
    for c in failed:
        print(f"[b_05] INVARIANT FAILED: arm={c['arm']} {c['statistic']} "
              f"expected {c['expected']} +/- {c['tolerance']}, measured {c['measured']} "
              f"(abs diff {c['abs_diff']})", file=sys.stderr, flush=True)
    if not failed:
        log(f"all {len(invariants)} expected invariants PASS")

    # --- declared-column corroboration --------------------------------------------------------
    disagree = [(nm, col) for nm, e in arms.items()
                for col, d in e["declared_columns_corroboration"].items()
                if isinstance(d, dict) and not d["agrees"]]
    for nm, col in disagree:
        print(f"[b_05] DECLARED-COLUMN DISAGREEMENT: arm={nm} column={col} does not match the "
              f"join-derived rate", file=sys.stderr, flush=True)

    tb_path = Path(a.trackb) if a.trackb else out / f"trackb_rungs_{a.stage}.json"
    trackb = json.loads(tb_path.read_text()) if tb_path.is_file() else None
    if trackb is None:
        log(f"b_04 report not found at {tb_path} -- GRN half of the donor axis recorded PENDING")

    report = {
        "stage": a.stage,
        "seed": a.seed,
        "reference_arm": TRUTH,
        "cohort": cohort_block(meta.reset_index()),
        "arms": arms,
        "expected_invariants": {
            "note": ("independently verified by the investigator on the uploaded tables before this "
                     "script existed; encoded here as assertions with a tolerance of half a unit "
                     "in the last quoted decimal place. A FAIL is a finding, not rounding."),
            "n_checks": len(invariants),
            "n_failed": len(failed),
            "all_pass": not failed,
            "checks": invariants,
        },
        "declared_column_disagreements": [{"arm": n, "column": c} for n, c in disagree],
        "donor_axis": donor_axis_block(arms, trackb, tb_path),
        "provenance": {
            "script": "pairing_error.py",
            "metadata": str(mp),
            "uploads_dir": str(uploads),
            "arm_tables": arm_files,
            "trackb_report": str(tb_path) if trackb is not None else None,
            "join_key": "h5ad_barcode (atac_barcode and assigned_rna_barcode both join to it)",
            "library_versions": {"python": sys.version.split()[0],
                                 "numpy": np.__version__, "pandas": pd.__version__},
        },
    }

    dest = out / (a.report_name or f"b05_pairing_error_{a.stage}.json")
    dest.write_text(json.dumps(report, indent=2))
    log(f"wrote {dest}")

    if failed or disagree:
        die(f"{len(failed)} invariant failure(s), {len(disagree)} declared-column "
            f"disagreement(s) -- report written, exiting non-zero")
    log("b_05 COMPLETE")

if __name__ == "__main__":
    main()
