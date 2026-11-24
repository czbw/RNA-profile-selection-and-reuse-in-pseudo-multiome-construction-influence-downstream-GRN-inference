#!/usr/bin/env python
'Score the topology levels L1-L5 of the resolution ladder against the native reference arm.\n\nL1 region->gene, L2 TF->gene, L3 eRegulon target/region sets, L4 eRegulon recovery, L5 TF recovery.\nFor paired arms the region-gene pair set is identical to native reference\'s by construction (asserted, not a\nresult), so the L1 "Jaccard of retained pairs" is degenerate and is not quoted.\n'
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

STAGE = "pbmc10k"
TRUTH = "truth"
ARMS_DEFAULT = [
    "link", "link_1to1", "pseudomultiome", "pseudomultiome_k10",
    "rand_subcluster_matched", "rand_subcluster_1to1", "rand_global",
]

# metacell (non-multiome) arms are not paired to native reference cell by cell
UNPAIRED = {"pseudomultiome", "pseudomultiome_k10"}

# helpers

def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return float(len(a & b)) / u if u else float("nan")

def recall_precision(ref: set, arm: set) -> dict:
    "recall = of native reference's items, how many did this arm recover; precision = of this arm's, how\n    many are native reference's. Reported both ways because an arm can inflate recall by emitting more."
    return {
        "n_truth": len(ref),
        "n_arm": len(arm),
        "n_shared": len(ref & arm),
        "recall": (len(ref & arm) / len(ref)) if ref else float("nan"),
        "precision": (len(ref & arm) / len(arm)) if arm else float("nan"),
        "jaccard": jaccard(ref, arm),
    }

def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    """Returns (rho, n_effective). n_effective is reported with every correlation because it differs
    across arms: donor reuse makes genes constant, which makes SCENIC+'s per-pair Spearman
    undefined, so a heavily reused arm drops its own worst pairs."""
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

def load_r2g(out: Path, arm: str) -> pd.DataFrame:
    return pd.read_csv(
        out / f"r2g_{STAGE}_{arm}.tsv",
        sep="\t",
        usecols=["target", "region", "importance", "rho", "Distance"],
        dtype={"target": str, "region": str, "importance": np.float64,
               "rho": np.float64, "Distance": str},
    )

def score_L1(out: Path, arms: list[str], topn: int) -> dict:
    ref = load_r2g(out, TRUTH)
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
        p = out / f"r2g_{STAGE}_{arm}.tsv"
        if not p.exists():
            res[arm] = {"status": "PENDING"}
            continue
        df = load_r2g(out, arm)
        paired = arm not in UNPAIRED

        entry: dict = {"n_pairs": int(len(df)), "paired_vs_truth": paired}

        if paired:
            # check the invariant (identical pair set and order); it is not a result (see module docstring).
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
            # metacell arm: compare distributions only, never pair
            entry["pair_set_identical_in_order"] = None
            entry["note"] = "metacell arm -- distribution comparison only, never paired"
            for lbl, col in (("importance", "importance"), ("rho", "rho")):
                a_ = ref[col].to_numpy()
                b_ = df[col].to_numpy()
                a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
                entry[f"ks_{lbl}"] = float(stats.ks_2samp(a_, b_).statistic)

        # Undefined rho is itself informative: SCENIC+'s per-pair Spearman is undefined when a gene
        # is constant across cells, and donor reuse creates constant genes, so the count reads out
        # how far a construction collapses its own RNA support.
        entry["n_rho_undefined"] = int(df["rho"].isna().sum())
        entry["pct_rho_undefined"] = float(df["rho"].isna().mean())

        # top-N selection -- informative at L1, unlike the full pair set
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

def score_L2(out: Path, arms: list[str], topn: int) -> dict:
    def load(arm):
        return pd.read_csv(
            out / f"tf2g_{STAGE}_{arm}.tsv", sep="\t",
            usecols=["TF", "target", "importance", "rho"],
            dtype={"TF": str, "target": str, "importance": np.float64, "rho": np.float64},
        )

    ref = load(TRUTH)
    ref_pairs = set(zip(ref["TF"], ref["target"]))
    ref_ord = np.argsort(-ref["importance"].to_numpy(), kind="stable")[:topn]
    ref_top = set(zip(ref["TF"].to_numpy()[ref_ord], ref["target"].to_numpy()[ref_ord]))

    res = {"_truth": {"n_edges": int(len(ref)), "top_n": topn, "n_TF": int(ref["TF"].nunique())}}
    for arm in arms:
        p = out / f"tf2g_{STAGE}_{arm}.tsv"
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
        if arm not in UNPAIRED:
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

def load_eregulons(out: Path, arm: str) -> pd.DataFrame:
    return pd.read_csv(
        out / f"eRegulons_{STAGE}_{arm}.tsv", sep="\t",
        usecols=["Region", "Gene", "TF", "is_extended", "eRegulon_name"],
        dtype={"Region": str, "Gene": str, "TF": str, "eRegulon_name": str},
    )

def target_sets(df: pd.DataFrame) -> dict[str, set]:
    return {k: set(v) for k, v in df.groupby("eRegulon_name", observed=True)["Gene"]}

def region_sets(df: pd.DataFrame) -> dict[str, set]:
    return {k: set(v) for k, v in df.groupby("eRegulon_name", observed=True)["Region"]}

def score_L345(out: Path, arms: list[str]) -> dict:
    ref = load_eregulons(out, TRUTH)
    ref_tg, ref_rg = target_sets(ref), region_sets(ref)
    ref_names, ref_tfs = set(ref_tg), set(ref["TF"])

    res = {
        "_truth": {
            "n_triplets": int(len(ref)),
            "n_eRegulons": len(ref_names),
            "n_TFs": len(ref_tfs),
            "n_eRegulons_direct": int(ref.loc[~ref["is_extended"], "eRegulon_name"].nunique()),
            # In this run is_extended is False for every row (no extended eRegulons were called),
            # so a direct-only breakout would be identical and is not reported as separate.
            "any_extended_eRegulons": bool(ref["is_extended"].any()),
            "median_target_set_size": float(np.median([len(v) for v in ref_tg.values()])),
        }
    }

    for arm in arms:
        p = out / f"eRegulons_{STAGE}_{arm}.tsv"
        if not p.exists():
            # a sentinel file marks arms for which no eRegulons could be built
            sentinel = out / f"eGRN_IMPOSSIBLE_{STAGE}_{arm}.json"
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
        df = load_eregulons(out, arm)
        tg, rg = target_sets(df), region_sets(df)
        names, tfs = set(tg), set(df["TF"])

        # --- L3: target-set agreement, only over eRegulons both arms recovered -------------
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
            # --- L4: which eRegulons exist at all --------------------------------------------
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

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str((Path(__file__).resolve().parents[3] / 'benchmarks/pbmc') / "out"))
    ap.add_argument("--topn", type=int, default=100_000)
    a = ap.parse_args()
    out = Path(a.out_dir)

    arms = [x for x in ARMS_DEFAULT]
    pending, impossible = [], []
    for x in arms:
        if (out / f"eRegulons_{STAGE}_{x}.tsv").exists():
            continue
        (impossible if (out / f"eGRN_IMPOSSIBLE_{STAGE}_{x}.json").exists() else pending).append(x)
    if pending:
        print(f"[p13] NOT YET AVAILABLE (scored as PENDING): {', '.join(pending)}")
    if impossible:
        print(f"[p13] eGRN STRUCTURALLY IMPOSSIBLE, L1/L2 still scored: {', '.join(impossible)}")

    print("[p13] L1 region->gene ...")
    L1 = score_L1(out, arms, a.topn)
    print("[p13] L2 TF->gene ...")
    L2 = score_L2(out, arms, a.topn)
    print("[p13] L3/L4/L5 eRegulons ...")
    L345 = score_L345(out, arms)

    report = {
        "stage": STAGE,
        "reference_arm": TRUTH,
        "top_n": a.topn,
        "declared_deviation": (
            "L1 'Jaccard of retained pairs' is degenerate by construction (r2g retains 100% of an "
            "identical search space) and is NOT quoted; see module docstring."
        ),
        "L1_region_to_gene": L1,
        "L2_tf_to_gene": L2,
        "L3_L4_L5_eRegulons": L345,
    }
    dest = out / f"trackb_rungs_{STAGE}.json"
    dest.write_text(json.dumps(report, indent=2))
    print(f"[p13] wrote {dest}")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
