#!/usr/bin/env python
"""Score the topology rungs L1-L5 of the BMMC pairing ladder."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent

# --- defaults for this ladder -----------------------------------------------------------------
# The four BMMC arms, in ladder order. rand_global is the floor and is expected to exit 45.
STAGE_DEFAULT = "bmmc"
TRUTH_DEFAULT = "truth"
ARMS_DEFAULT = ["rand_celltype_1to1", "rand_celltype_x_donor_1to1", "rand_global"]
# No metacell arm in the default BMMC set (the metacell arms resolve but cannot run here:
# native_rna_bmmc.metacell.h5ad has no builder in this bundle). Empty, not absent: --unpaired
# is the option that would re-enable distribution-only handling if one is ever added.
UNPAIRED_DEFAULT: list[str] = []

# The PBMC ladder, as the parity target. Kept here rather than typed on the command line so the
# parity check cannot be weakened by a typo in an arm name.
PARITY = {
    "stage": "pbmc10k",
    "truth": "truth",
    "arms": ["link", "link_1to1", "pseudomultiome", "pseudomultiome_k10",
             "rand_subcluster_matched", "rand_subcluster_1to1", "rand_global"],

    "unpaired": ["pseudomultiome", "pseudomultiome_k10"],
    "topn": 100_000,
    "report": "trackb_rungs_pbmc10k.json",
}

# Leaf paths this port adds to p_13's report shape. Every one of these is structural (provenance /
# run configuration), carries no statistic, and is skipped by the parity diff. Anything else that
# appears on one side and not the other fails the parity check.
PARITY_ALLOWED_EXTRA = (
    "provenance",          # whole subtree: script name, inputs, arm roles, library versions
    "arms_scored",         # echo of --arms, so a report can never be read against the wrong list
    "unpaired_arms",       # echo of --unpaired; [] for the BMMC default
    "pending_arms",        # arms with neither eRegulons nor a sentinel
    "impossible_arms",     # arms with an eGRN_IMPOSSIBLE sentinel
    # The donor-awareness contrast lives in pairing_error.py, not here. b_04 is a pure
    # producer and b_05 consumes this report; merging b_05's output back in would make the two
    # scripts mutually dependent and force a two-pass run.
)

def log(m):
    print(f"[b_04] {m}", flush=True)

def die(msg):
    print(f"[b_04] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def die_missing(p: Path, why: str):
    """A missing input is named, with the path and the reason it was required. No silent skip and
    no zero-filled placeholder -- a zero here would be indistinguishable from a measured zero."""
    die(f"required input missing: {p}\n       reason: {why}")

# helpers -- identical to p_13's. Do not change them: the parity check compares against p_13.

def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return float(len(a & b)) / u if u else float("nan")

def recall_precision(ref: set, arm: set) -> dict:
    "recall = of native reference's items, how many did this arm recover; precision = of this arm's, how\n    many are native reference's. Deliberately asymmetric -- an arm can inflate its way to high recall."
    return {
        "n_truth": len(ref),
        "n_arm": len(arm),
        "n_shared": len(ref & arm),
        "recall": (len(ref & arm) / len(ref)) if ref else float("nan"),
        "precision": (len(ref & arm) / len(arm)) if arm else float("nan"),
        "jaccard": jaccard(ref, arm),
    }

def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    """Returns (rho, n_effective). n_effective is reported alongside every correlation because it
    is not constant across arms: donor reuse makes genes constant, which makes SCENIC+'s per-pair
    Spearman undefined, so a heavily reused arm silently drops its own worst pairs. A correlation
    quoted without its support would flatter such arms."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 100:
        return float("nan"), int(ok.sum())
    return float(stats.spearmanr(x[ok], y[ok]).statistic), int(ok.sum())

def parse_distance(s: pd.Series) -> np.ndarray:
    """SCENIC+ writes Distance as a stringified 1-element list, e.g. '[-126061]'."""
    return pd.to_numeric(
        s.astype(str).str.strip("[]").str.split(",").str[0], errors="coerce"
    ).to_numpy(dtype=float)

# L1 -- region -> gene

def load_r2g(out: Path, stage: str, arm: str) -> pd.DataFrame:
    return pd.read_csv(
        out / f"r2g_{stage}_{arm}.tsv",
        sep="\t",
        usecols=["target", "region", "importance", "rho", "Distance"],
        dtype={"target": str, "region": str, "importance": np.float64,
               "rho": np.float64, "Distance": str},
    )

def score_L1(out: Path, stage: str, truth: str, arms: list[str], unpaired: set[str],
             topn: int) -> dict:
    ref_p = out / f"r2g_{stage}_{truth}.tsv"
    if not ref_p.exists():
        die_missing(ref_p, f"the reference arm '{truth}' anchors every L1 comparison; without its "
                           f"region->gene table there is nothing to score against")
    ref = load_r2g(out, stage, truth)
    ref_dist = parse_distance(ref["Distance"])
    res: dict[str, dict] = {}

    # top-N by importance for native reference, as the reference selection
    ref_ord = np.argsort(-ref["importance"].to_numpy(), kind="stable")[:topn]
    ref_top = set(zip(ref["region"].to_numpy()[ref_ord], ref["target"].to_numpy()[ref_ord]))

    res["_truth"] = {
        "n_pairs": int(len(ref)),
        "top_n": topn,
        "median_abs_distance_topN": float(np.nanmedian(np.abs(ref_dist[ref_ord]))),
        "median_abs_distance_all": float(np.nanmedian(np.abs(ref_dist))),
    }

    for arm in arms:
        p = out / f"r2g_{stage}_{arm}.tsv"
        if not p.exists():
            res[arm] = {"status": "PENDING"}
            continue
        df = load_r2g(out, stage, arm)
        paired = arm not in unpaired

        entry: dict = {"n_pairs": int(len(df)), "paired_vs_truth": paired}

        if paired:
            # assert the invariant (it is an invariant, not a result).
            same_order = (
                len(df) == len(ref)
                and np.array_equal(df["region"].to_numpy(), ref["region"].to_numpy())
                and np.array_equal(df["target"].to_numpy(), ref["target"].to_numpy())
            )
            if not same_order:
                # fall back to an explicit merge rather than trusting row order
                m = ref.merge(df, on=["region", "target"], suffixes=("_t", "_a"))
                entry["pair_set_identical_in_order"] = False
                entry["n_shared_pairs"] = int(len(m))
                xi, yi = m["importance_t"].to_numpy(), m["importance_a"].to_numpy()
                xr, yr = m["rho_t"].to_numpy(), m["rho_a"].to_numpy()
            else:
                entry["pair_set_identical_in_order"] = True
                entry["n_shared_pairs"] = int(len(df))
                xi, yi = ref["importance"].to_numpy(), df["importance"].to_numpy()
                xr, yr = ref["rho"].to_numpy(), df["rho"].to_numpy()

            entry["spearman_importance"], entry["n_eff_importance"] = spearman(xi, yi)  # all of L1
            entry["spearman_rho"], entry["n_eff_rho"] = spearman(xr, yr)
        else:
            # metacell arm: distributions only, never paired
            entry["pair_set_identical_in_order"] = None
            entry["note"] = "metacell arm -- distribution comparison only, never paired"
            for lbl, col in (("importance", "importance"), ("rho", "rho")):
                a_ = ref[col].to_numpy()
                b_ = df[col].to_numpy()
                a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
                entry[f"ks_{lbl}"] = float(stats.ks_2samp(a_, b_).statistic)

        # Undefined rho is a result, not a nuisance: SCENIC+'s per-pair Spearman is undefined when
        # a gene is constant across cells, and donor reuse manufactures constancy. The count is
        # therefore a direct readout of how far a construction collapses its own RNA support.
        entry["n_rho_undefined"] = int(df["rho"].isna().sum())
        entry["pct_rho_undefined"] = float(df["rho"].isna().mean())

        # top-N selection -- non-degenerate at L1, unlike the full pair set
        ordr = np.argsort(-df["importance"].to_numpy(), kind="stable")[:topn]
        arm_top = set(zip(df["region"].to_numpy()[ordr], df["target"].to_numpy()[ordr]))
        entry["topN_jaccard_vs_truth"] = jaccard(ref_top, arm_top)
        d = parse_distance(df["Distance"])
        entry["median_abs_distance_topN"] = float(np.nanmedian(np.abs(d[ordr])))
        entry["median_abs_distance_all"] = float(np.nanmedian(np.abs(d)))

        res[arm] = entry
        del df
    return res

# L2 -- TF -> gene

def score_L2(out: Path, stage: str, truth: str, arms: list[str], unpaired: set[str],
             topn: int) -> dict:
    def load(arm):
        return pd.read_csv(
            out / f"tf2g_{stage}_{arm}.tsv", sep="\t",
            usecols=["TF", "target", "importance", "rho"],
            dtype={"TF": str, "target": str, "importance": np.float64, "rho": np.float64},
        )

    ref_p = out / f"tf2g_{stage}_{truth}.tsv"
    if not ref_p.exists():
        die_missing(ref_p, f"the reference arm '{truth}' anchors every L2 comparison; without its "
                           f"TF->gene table there is nothing to score against")
    ref = load(truth)
    ref_pairs = set(zip(ref["TF"], ref["target"]))
    ref_ord = np.argsort(-ref["importance"].to_numpy(), kind="stable")[:topn]
    ref_top = set(zip(ref["TF"].to_numpy()[ref_ord], ref["target"].to_numpy()[ref_ord]))

    res = {"_truth": {"n_edges": int(len(ref)), "top_n": topn, "n_TF": int(ref["TF"].nunique())}}
    for arm in arms:
        p = out / f"tf2g_{stage}_{arm}.tsv"
        if not p.exists():
            res[arm] = {"status": "PENDING"}
            continue
        df = load(arm)
        arm_pairs = set(zip(df["TF"], df["target"]))
        ordr = np.argsort(-df["importance"].to_numpy(), kind="stable")[:topn]
        arm_top = set(zip(df["TF"].to_numpy()[ordr], df["target"].to_numpy()[ordr]))

        entry = {
            "n_edges": int(len(df)),
            "n_TF": int(df["TF"].nunique()),
            "full_edge_set": recall_precision(ref_pairs, arm_pairs),
            "topN": recall_precision(ref_top, arm_top),
        }
        # hypergeometric fold-enrichment of the top-N overlap, over the shared candidate universe
        universe = len(ref_pairs | arm_pairs)
        exp = len(ref_top) * len(arm_top) / universe if universe else float("nan")
        entry["topN_fold_enrichment"] = (
            float(len(ref_top & arm_top) / exp) if exp else float("nan")
        )
        entry["topN_expected_overlap_by_chance"] = float(exp)

        # paired importance concordance on the shared edge set
        if arm not in unpaired:
            m = ref.merge(df, on=["TF", "target"], suffixes=("_t", "_a"))
            entry["n_shared_edges"] = int(len(m))
            entry["spearman_importance_shared"], entry["n_eff_importance"] = spearman(
                m["importance_t"].to_numpy(), m["importance_a"].to_numpy()
            )
            entry["spearman_rho_shared"], entry["n_eff_rho"] = spearman(
                m["rho_t"].to_numpy(), m["rho_a"].to_numpy()
            )
            del m
        else:
            entry["note"] = "metacell arm -- shared-edge Spearman not reported as paired"
        res[arm] = entry
        del df
    return res

# L3/L4/L5 -- from the eRegulon table

def load_eregulons(out: Path, stage: str, arm: str) -> pd.DataFrame:
    return pd.read_csv(
        out / f"eRegulons_{stage}_{arm}.tsv", sep="\t",
        usecols=["Region", "Gene", "TF", "is_extended", "eRegulon_name"],
        dtype={"Region": str, "Gene": str, "TF": str, "eRegulon_name": str},
    )

def target_sets(df: pd.DataFrame) -> dict[str, set]:
    return {k: set(v) for k, v in df.groupby("eRegulon_name", observed=True)["Gene"]}

def region_sets(df: pd.DataFrame) -> dict[str, set]:
    return {k: set(v) for k, v in df.groupby("eRegulon_name", observed=True)["Region"]}

def score_L345(out: Path, stage: str, truth: str, arms: list[str]) -> dict:
    ref_p = out / f"eRegulons_{stage}_{truth}.tsv"
    if not ref_p.exists():
        # The reference arm cannot be "structurally impossible" and still be the reference: if
        # native reference supports no eGRN there is no ladder to score, and that is a different report.
        die_missing(ref_p, f"the reference arm '{truth}' anchors L3-L5. If the reference itself "
                           f"exited 45 there is no eGRN to recover and the ladder is void")
    ref = load_eregulons(out, stage, truth)
    ref_tg, ref_rg = target_sets(ref), region_sets(ref)
    ref_names, ref_tfs = set(ref_tg), set(ref["TF"])

    res = {
        "_truth": {
            "n_triplets": int(len(ref)),
            "n_eRegulons": len(ref_names),
            "n_TFs": len(ref_tfs),
            "n_eRegulons_direct": int(ref.loc[~ref["is_extended"], "eRegulon_name"].nunique()),
            # In the PBMC run is_extended was False for every row -- no extended eRegulons were
            # called at all. The direct-only breakout is then vacuous and is suppressed below
            # rather than reported as an independent analysis. Whether that holds for BMMC is
            # measured here, not assumed: the flag below is read off the actual table.
            "any_extended_eRegulons": bool(ref["is_extended"].any()),
            "median_target_set_size": float(np.median([len(v) for v in ref_tg.values()])),
        }
    }

    for arm in arms:
        p = out / f"eRegulons_{stage}_{arm}.tsv"
        if not p.exists():

            sentinel = out / f"eGRN_IMPOSSIBLE_{stage}_{arm}.json"
            if sentinel.exists():
                res[arm] = {
                    "status": "eGRN_STRUCTURALLY_IMPOSSIBLE",
                    "evidence": json.loads(sentinel.read_text()),
                    "L3_target_sets": "does not exist -- no eRegulons were constructible",
                    "L4_eRegulon_recovery": {
                        "n_truth": len(ref_names), "n_arm": 0, "n_shared": 0,
                        "recall": 0.0, "precision": None, "jaccard": 0.0,
                        "note": "precision is UNDEFINED (empty prediction set), not 0.0",
                    },
                    "L5_TF_recovery": {
                        "n_truth": len(ref_tfs), "n_arm": 0, "n_shared": 0,
                        "recall": 0.0, "precision": None, "jaccard": 0.0,
                        "note": "precision is UNDEFINED (empty prediction set), not 0.0",
                    },
                }
            else:
                res[arm] = {"status": "PENDING"}
            continue
        df = load_eregulons(out, stage, arm)
        tg, rg = target_sets(df), region_sets(df)
        names, tfs = set(tg), set(df["TF"])

        # --- L3: target-set agreement, only over eRegulons BOTH arms recovered ---------------
        shared = sorted(ref_names & names)
        jt = np.array([jaccard(ref_tg[k], tg[k]) for k in shared], dtype=float)
        jr = np.array([jaccard(ref_rg[k], rg[k]) for k in shared], dtype=float)

        def summarise(j, label):
            if j.size == 0:
                return {"n": 0}
            return {
                "n": int(j.size),
                f"median_jaccard_{label}": float(np.median(j)),
                f"iqr_{label}": [float(np.percentile(j, 25)), float(np.percentile(j, 75))],
                f"pct_above_0.5_{label}": float((j > 0.5).mean()),
                f"pct_below_0.1_{label}": float((j < 0.1).mean()),
            }

        res[arm] = {
            "n_triplets": int(len(df)),
            "L3_target_sets": summarise(jt, "genes"),
            "L3_region_sets": summarise(jr, "regions"),
            # --- L4: which eRegulons exist at all -------------------------------------------
            "L4_eRegulon_recovery": recall_precision(ref_names, names),
            "L4_eRegulon_recovery_direct_only": (
                recall_precision(
                    set(ref.loc[~ref["is_extended"], "eRegulon_name"]),
                    set(df.loc[~df["is_extended"], "eRegulon_name"]),
                )
                if (ref["is_extended"].any() or df["is_extended"].any())
                else "NOT APPLICABLE -- no extended eRegulons were called in either arm; "
                     "identical to L4_eRegulon_recovery by construction"
            ),
            # --- L5: which TFs survive at all ------------------------------------------------
            "L5_TF_recovery": recall_precision(ref_tfs, tfs),
        }
        del df
    return res

# arm triage -- pending vs structurally impossible, decided once and reused

def triage(out: Path, stage: str, arms: list[str]) -> tuple[list[str], list[str]]:
    pending, impossible = [], []
    for x in arms:
        if (out / f"eRegulons_{stage}_{x}.tsv").exists():
            continue
        if (out / f"eGRN_IMPOSSIBLE_{stage}_{x}.json").exists():
            impossible.append(x)

            for kind in ("r2g", "tf2g"):
                q = out / f"{kind}_{stage}_{x}.tsv"
                if not q.exists():
                    die_missing(
                        q,
                        f"arm '{x}' carries eGRN_IMPOSSIBLE_{stage}_{x}.json, which egrn_guard.py "
                        f"can only have written after reading the r2g table. A sentinel without "
                        f"its {kind} table means the out-dir is incomplete, not that the arm is "
                        f"un-scorable at L1/L2",
                    )
        else:
            pending.append(x)
    return pending, impossible

# the parity check

def flatten(obj, prefix="") -> dict:
    """Leaf-path -> value. Lists are indexed, so a reordered iqr pair is a difference, not a tie."""
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out.update(flatten(v, f"{prefix}[{i}]"))
    else:
        out[prefix] = obj
    return out

def values_equal(a, b, tol: float) -> bool:
    "Exact by default (tol=0.0). NaN == NaN here: both mean 'the estimator declined', and a\n    report where native reference's NaN became the arm's NaN is a match, not a mismatch."
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return True
        if a == b:
            return True
        return tol > 0.0 and abs(float(a) - float(b)) <= tol
    return a == b

def parity_check(ref_json: Path, report: dict, tol: float) -> int:
    """Parity check."""
    if not ref_json.exists():
        die_missing(ref_json, "the parity gate compares against the committed PBMC report")
    ref = json.loads(ref_json.read_text())
    A, B = flatten(ref), flatten(report)

    def allowed(path: str) -> bool:
        head = path.split(".", 1)[0].split("[", 1)[0]
        return head in PARITY_ALLOWED_EXTRA

    only_ref = sorted(k for k in A if k not in B)
    only_new = sorted(k for k in B if k not in A)
    extra_declared = [k for k in only_new if allowed(k)]
    extra_undeclared = [k for k in only_new if not allowed(k)]
    mism = [(k, A[k], B[k]) for k in sorted(A) if k in B and not values_equal(A[k], B[k], tol)]

    log("=" * 96)
    log(f"PARITY CHECK  vs {ref_json}")
    log(f"  leaves in reference report : {len(A)}")
    log(f"  leaves in this port        : {len(B)}")
    log(f"  compared (present in both) : {len(set(A) & set(B))}")
    log(f"  identical                  : {len(set(A) & set(B)) - len(mism)}")
    log(f"  structural extras (declared, skipped) : {len(extra_declared)}")
    for k in extra_declared:
        log(f"      + {k}")
    if only_ref:
        log(f"  MISSING from this port ({len(only_ref)}):")
        for k in only_ref:
            log(f"      - {k}  (reference value: {A[k]!r})")
    if extra_undeclared:
        log(f"  UNDECLARED extras in this port ({len(extra_undeclared)}):")
        for k in extra_undeclared:
            log(f"      + {k} = {B[k]!r}")
    if mism:
        log(f"  VALUE MISMATCHES ({len(mism)}):")
        for k, a, b in mism:
            d = ""
            if isinstance(a, (int, float)) and isinstance(b, (int, float)) \
                    and not isinstance(a, bool) and not isinstance(b, bool):
                d = f"   (abs diff {abs(float(a) - float(b)):.3e})"
            log(f"      ! {k}\n          reference: {a!r}\n          this port: {b!r}{d}")

    ok = not (only_ref or extra_undeclared or mism)
    log("PARITY CHECK: " + ("PASS -- every compared leaf is identical" if ok else "FAIL"))
    log("=" * 96)
    return 0 if ok else 1

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Track B (L1-L5) topology scoring for a pairing ladder; b_04 of the BMMC bundle"
    )
    ap.add_argument("--out-dir", default=str(BUNDLE / "out"),
                    help="directory holding r2g_/tf2g_/eRegulons_/eGRN_IMPOSSIBLE_ files")
    ap.add_argument("--stage", default=STAGE_DEFAULT,
                    help="the STAGE token in every filename (TAG = <STAGE>_<ARM>)")
    ap.add_argument("--reference-arm", default=TRUTH_DEFAULT,
                    help="the arm every other arm is scored against (barcode identity)")
    ap.add_argument("--arms", nargs="+", default=None,
                    help="arms to score, excluding the reference (default: the BMMC ladder)")
    ap.add_argument("--unpaired", nargs="*", default=None,
                    help="metacell arms: distribution-only at L1/L2, never paired. "
                         "Default for BMMC is EMPTY -- no metacell arm is in the default set")
    ap.add_argument("--topn", type=int, default=100_000)
    ap.add_argument("--report-name", default=None,
                    help="output filename (default: trackb_rungs_<STAGE>.json)")
    ap.add_argument("--parity-check", nargs="?", const="AUTO", default=None, metavar="REF_JSON",
                    help="run against the PBMC out-dir and diff every leaf against the committed "
                         "trackb_rungs_pbmc10k.json. With no value, the PBMC stage/arms/unpaired/"
                         "topn are forced so the gate cannot be weakened by a wrong argument")
    ap.add_argument("--parity-tol", type=float, default=0.0,
                    help="absolute tolerance for the parity diff (default 0.0 = bit-exact)")
    ap.add_argument("--no-write", action="store_true",
                    help="compute and diff but do not write a report (parity runs use this)")
    a = ap.parse_args()

    out = Path(a.out_dir)
    if not out.is_dir():
        die(f"--out-dir is not a directory: {out}")

    parity_ref = None
    if a.parity_check is not None:

        if a.parity_check == "AUTO":
            parity_ref = out / PARITY["report"]
        else:
            parity_ref = Path(a.parity_check)
        for flag, forced in (("stage", PARITY["stage"]), ("reference_arm", PARITY["truth"]),
                             ("topn", PARITY["topn"])):
            cur = getattr(a, flag)
            default = {"stage": STAGE_DEFAULT, "reference_arm": TRUTH_DEFAULT,
                       "topn": 100_000}[flag]
            if cur != forced and cur != default:
                die(f"--parity-check forces --{flag.replace('_', '-')}={forced!r}, "
                    f"but {cur!r} was passed")
            setattr(a, flag, forced)
        if a.arms is not None and list(a.arms) != PARITY["arms"]:
            die(f"--parity-check forces --arms {' '.join(PARITY['arms'])}; got {a.arms}")
        a.arms = list(PARITY["arms"])
        if a.unpaired is not None and sorted(a.unpaired) != sorted(PARITY["unpaired"]):
            die(f"--parity-check forces --unpaired {' '.join(PARITY['unpaired'])}; got {a.unpaired}")
        a.unpaired = list(PARITY["unpaired"])
        a.no_write = True

    stage = a.stage
    truth = a.reference_arm
    arms = list(a.arms) if a.arms is not None else list(ARMS_DEFAULT)
    unpaired = set(a.unpaired if a.unpaired is not None else UNPAIRED_DEFAULT)

    if truth in arms:
        # Scoring the reference against itself yields 1.0 everywhere and would be read as a result.
        arms = [x for x in arms if x != truth]
        log(f"reference arm {truth!r} was in --arms; removed (it is the yardstick, not a rung)")
    unknown = unpaired - set(arms)
    if unknown:
        die(f"--unpaired names arm(s) not in --arms: {sorted(unknown)}")

    log(f"stage={stage}  reference={truth}  topn={a.topn}")
    log(f"arms={arms}")
    log(f"unpaired (metacell, distribution-only)={sorted(unpaired) or '[] (none in this set)'}")
    log(f"out-dir={out}")

    pending, impossible = triage(out, stage, arms)
    if pending:
        print(f"[b_04] NOT YET AVAILABLE (scored as PENDING): {', '.join(pending)}")
    if impossible:
        print(f"[b_04] eGRN STRUCTURALLY IMPOSSIBLE, L1/L2 still scored: {', '.join(impossible)}")

    print("[b_04] L1 region->gene ...")
    L1 = score_L1(out, stage, truth, arms, unpaired, a.topn)
    print("[b_04] L2 TF->gene ...")
    L2 = score_L2(out, stage, truth, arms, unpaired, a.topn)
    print("[b_04] L3/L4/L5 eRegulons ...")
    L345 = score_L345(out, stage, truth, arms)

    # Report shape: the first seven keys are p_13's, byte-for-byte, so the parity check compares
    # like with like. Everything this port adds is listed in PARITY_ALLOWED_EXTRA.
    report = {
        "stage": stage,
        "reference_arm": truth,
        "top_n": a.topn,
        "declared_deviation": (
            "L1 'Jaccard of retained pairs' is degenerate by construction (r2g retains 100% of an "
            "identical search space) and is NOT quoted."
        ),
        "L1_region_to_gene": L1,
        "L2_tf_to_gene": L2,
        "L3_L4_L5_eRegulons": L345,
        # ---- additions over p_13 (all structural; see PARITY_ALLOWED_EXTRA) ------------------
        "arms_scored": arms,
        "unpaired_arms": sorted(unpaired),
        "pending_arms": pending,
        "impossible_arms": impossible,
        "provenance": {
            "script": "network_agreement.py",
            "ported_from": "benchmarks/pbmc/scripts/network_agreement.py",
            "out_dir": str(out),
            "inputs": {
                arm: {
                    kind: (str(out / f"{kind}_{stage}_{arm}.tsv")
                           if (out / f"{kind}_{stage}_{arm}.tsv").exists() else None)
                    for kind in ("r2g", "tf2g", "eRegulons")
                } | {
                    "eGRN_IMPOSSIBLE": (
                        str(out / f"eGRN_IMPOSSIBLE_{stage}_{arm}.json")
                        if (out / f"eGRN_IMPOSSIBLE_{stage}_{arm}.json").exists() else None
                    )
                }
                for arm in [truth] + arms
            },
            "library_versions": {
                "python": sys.version.split()[0],
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": __import__("scipy").__version__,
            },
        },
    }

    if not a.no_write:
        dest = out / (a.report_name or f"trackb_rungs_{stage}.json")
        dest.write_text(json.dumps(report, indent=2))
        print(f"[b_04] wrote {dest}")

    if parity_ref is not None:
        raise SystemExit(parity_check(parity_ref, report, a.parity_tol))

    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
