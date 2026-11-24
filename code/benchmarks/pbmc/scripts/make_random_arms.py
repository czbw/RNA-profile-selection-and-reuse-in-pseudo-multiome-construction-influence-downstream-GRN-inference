"""The five random pairing arms, with reuse matching."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')

# Minimum cells for a stratum to host a random arm. Distinct from config.yaml's
# min_cells (50, pycisTopic), min_cells_per_group (100) and min_cells_per_celltype (25, metacells).
MIN_STRATUM_CELLS_RANDOM_ARM = 50

ARMS = ["rand_subcluster_matched", "rand_subcluster_1to1", "rand_lineage_matched",
        "rand_compartment_matched", "rand_global"]

def log(m):
    print(f"[p07] {m}", flush=True)

# --------------------------------------------------------------- label ladder ----
def to_lineage(sub: str) -> str:
    """subcluster -> lineage. Splits on the first dot: CD4.T.1 -> CD4, Monocyte.CD14.1 -> Monocyte."""
    head = sub.split(".")[0]
    return {"CD4": "CD4", "NK": "NK", "B": "B", "Monocyte": "Monocyte",
            "DC": "DC", "pDC": "pDC", "HSPC": "HSPC"}.get(head, head)

def to_compartment(lin: str) -> str:
    """Lineage -> compartment. pDC is placed in myeloid with the other dendritic populations; this is
    a modelling choice, and a sensitivity analysis with pDC in lymphoid is run separately. HSPC is
    excluded upstream by the 50-cell floor, so no progenitor compartment exists."""
    if lin in ("CD4", "NK", "B"):
        return "lymphoid"
    if lin in ("Monocyte", "DC", "pDC"):
        return "myeloid"
    return lin        # HSPC and anything unexpected keep their own name; floor will drop them

# ------------------------------------------------------------------- sampler ----
def match_stratum(rng, atac_ids, true_bc, link_bc, pool, one_to_one=False):
    """Return {atac_id -> assigned_rna_barcode} for ONE stratum, plus its diagnostics."""
    n_s = len(atac_ids)
    if one_to_one:
        d_s = n_s
        M_s = np.ones(n_s, dtype=int)
    else:
        counts = pd.Series(link_bc).value_counts()
        d_s = int(counts.size)
        M_s = counts.to_numpy()

    # d_s must not exceed n_s; failure means labels are inconsistent between modalities.
    assert d_s <= n_s, f"d_s={d_s} > n_s={n_s}: labels inconsistent between modalities"
    assert M_s.sum() == n_s, f"multiplicity sums to {M_s.sum()}, expected n_s={n_s}"
    # the pool must be able to supply d_s distinct donors without replacement
    assert d_s <= len(pool), f"d_s={d_s} exceeds donor pool {len(pool)} for this stratum"

    donors = rng.choice(np.asarray(pool), size=d_s, replace=False)
    M_shuf = M_s.copy()
    rng.shuffle(M_shuf)
    cells = np.asarray(atac_ids).copy()
    rng.shuffle(cells)

    assign, k = {}, 0
    for donor, m in zip(donors, M_shuf):
        for c in cells[k:k + m]:
            assign[c] = donor
        k += m
    assert k == n_s

    achieved = pd.Series(list(assign.values())).value_counts().to_numpy()
    diag = {
        "n_s": int(n_s), "d_s_target": int(d_s), "d_s_achieved": int(len(achieved)),
        "reuse_target": float(n_s / d_s), "reuse_achieved": float(n_s / len(achieved)),
        "ks_multiplicity": float(ks_distance(M_s, achieved)),
        "self_hits": int(sum(assign[c] == t for c, t in zip(atac_ids, true_bc))),
    }
    return assign, diag

def ks_distance(a, b):
    """Two-sample KS statistic between two multiplicity vectors. Should be 0.0 by construction --
    the achieved vector is M_s permuted. Anything else means the draw or the assignment is wrong."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    grid = np.union1d(a, b)
    ca = np.searchsorted(a, grid, side="right") / a.size
    cb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))

def build_arm(df, key, rng, one_to_one=False):
    assign, per_stratum = {}, {}
    for s, g in df.groupby(key, sort=True):
        pool = g["true_barcode"].unique()          # donor pool: the analysis set for this stratum
        a, d = match_stratum(rng, g["cistopic_id"].to_numpy(), g["true_barcode"].to_numpy(),
                             g["Linked_scRNA_barcode"].to_numpy(), pool, one_to_one)
        assign.update(a)
        per_stratum[str(s)] = d
    return assign, per_stratum

# ---------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["3k", "10k"])
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--min-stratum", type=int, default=MIN_STRATUM_CELLS_RANDOM_ARM,
                    help="minimum cells per stratum for the ladder statistics (default 50). Pass 0 for "
                         "the GRN arms, where the arm object must cover the same cells as truth/link "
                         "(all 10,925) or the arms would differ in cells as well as pairing. "
                         "run_grn.sh applies the 25-cell metacell floor to pseudomultiome only, "
                         "so truth and link were inferred on the full manifest.")
    a = ap.parse_args()
    floor = a.min_stratum

    WORK = ROOT / "work"
    OUT = Path(a.outdir) if a.outdir else ROOT / "out"
    OUT.mkdir(parents=True, exist_ok=True)

    # --- assemble: ATAC cell -> true label, true barcode, link's assignment ---------------
    if a.dataset == "10k":
        df = pd.read_csv(WORK / "cell_data_pbmc10k.csv")
        # p_03 names the ID column `cell_id` in the 10k manifest, whereas the 3k manifest from
        # p_00/p_01 uses `cistopic_id`; normalise to `cistopic_id`.
        if "cistopic_id" not in df.columns and "cell_id" in df.columns:
            df = df.rename(columns={"cell_id": "cistopic_id"})
        label_col = "Classified_Celltype"
    else:
        lab = pd.read_csv(WORK / "cell_data.csv")                       # 3k dry-run annotation
        lnk = pd.read_csv(WORK / "3k" / "archr" / "link_3k.csv")
        df = lab.merge(lnk[["cistopic_id", "Linked_scRNA_barcode", "true_barcode", "link_is_correct"]],
                       on="cistopic_id", how="inner")
        label_col = "annotated_celltype"
    log(f"dataset {a.dataset}: {len(df):,} ATAC cells with a label AND a link")

    df = df.rename(columns={label_col: "subcluster"})
    df["lineage"] = df["subcluster"].map(to_lineage)
    df["compartment"] = df["lineage"].map(to_compartment)
    df["global"] = "all"

    # --- the floor, applied once at the subcluster rung --------------------------------
    sizes = df["subcluster"].value_counts()
    admitted = sorted(sizes[sizes >= floor].index)
    dropped = {str(k): int(v) for k, v in sizes[sizes < floor].items()}
    log(f"subclusters: {len(sizes)} present, {len(admitted)} admitted at "
        f"min_stratum_cells_random_arm={floor}"
        + ("   [floor disabled, full manifest]" if floor == 0 else ""))
    if dropped:
        log(f"  BELOW FLOOR, excluded and reported: {dropped}")
    df = df[df["subcluster"].isin(admitted)].reset_index(drop=True)
    N = len(df)
    log(f"ANALYSIS SET = {N:,} cells, constant across every rung")

    tagsuf = "_full" if floor == 0 else ""
    rng = np.random.default_rng(a.seed)
    report = {"dataset": a.dataset, "seed": a.seed,
              "min_stratum_cells_random_arm": floor,
              "analysis_set_N": N, "subclusters_admitted": len(admitted),
              "subclusters_below_floor": dropped,
              "strata_per_rung": {r: int(df[r].nunique()) for r in
                                  ("subcluster", "lineage", "compartment", "global")},
              "arms": {}}

    for rung in ("subcluster", "lineage", "compartment", "global"):
        log(f"  rung {rung:12s}: {df[rung].nunique():2d} strata  "
            f"[{', '.join(f'{k}={v}' for k, v in df[rung].value_counts().items())}]")

    SPEC = [("rand_subcluster_matched", "subcluster", False),
            ("rand_subcluster_1to1",    "subcluster", True),
            ("rand_lineage_matched",    "lineage",    False),
            ("rand_compartment_matched", "compartment", False),
            ("rand_global",             "global",     False)]

    for arm, key, one2one in SPEC:
        assign, per_stratum = build_arm(df, key, rng, one2one)
        n_strata = df[key].nunique()
        hits = sum(per_stratum[s]["self_hits"] for s in per_stratum)
        exp_hits = n_strata                      # expected self-hits: exactly 1 per stratum
        ks_max = max(d["ks_multiplicity"] for d in per_stratum.values())
        reuse_t = N / sum(d["d_s_target"] for d in per_stratum.values())
        reuse_a = N / sum(d["d_s_achieved"] for d in per_stratum.values())

        out = pd.DataFrame({"cistopic_id": df["cistopic_id"],
                            "Linked_scRNA_barcode": df["cistopic_id"].map(assign),
                            "true_barcode": df["true_barcode"],
                            "stratum": df[key]})
        out["link_is_correct"] = out["Linked_scRNA_barcode"] == out["true_barcode"]
        out.to_csv(OUT / f"arm_{arm}_{a.dataset}_seed{a.seed}{tagsuf}.csv", index=False)

        report["arms"][arm] = {
            "rung": key, "n_strata": int(n_strata),
            "reuse_target": float(reuse_t), "reuse_achieved": float(reuse_a),
            "max_ks_multiplicity": float(ks_max),
            "self_hits_observed": int(hits), "self_hits_expected": int(exp_hits),
            "self_hit_rate_observed": float(hits / N),
            "self_hit_rate_expected": float(exp_hits / N),
            "per_stratum": per_stratum}

        log(f"{arm:28s} rung={key:12s} strata={n_strata:2d} "
            f"reuse {reuse_a:5.2f}x (target {reuse_t:5.2f}x) maxKS={ks_max:.3f} "
            f"self-hits {hits:3d} (expected {exp_hits:3d}, rate {hits/N:.4%} vs {exp_hits/N:.4%})")

    # link's own self-hit rate on the same analysis set
    link_hits = int((df["Linked_scRNA_barcode"] == df["true_barcode"]).sum())
    report["link_on_analysis_set"] = {
        "self_hits": link_hits, "self_hit_rate": float(link_hits / N),
        "vs_subcluster_oracle_expectation": float(link_hits / df["subcluster"].nunique())}
    log(f"\nlink on the same analysis set: {link_hits} self-hits = {link_hits/N:.4%}, "
        f"{link_hits/df['subcluster'].nunique():.2f}x the oracle-subcluster expectation")

    with open(OUT / f"random_arms_report_{a.dataset}_seed{a.seed}{tagsuf}.json", "w") as fh:
        json.dump(report, fh, indent=2)
    log(f"wrote {OUT}/random_arms_report_{a.dataset}_seed{a.seed}.json")

if __name__ == "__main__":
    main()
