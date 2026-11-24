#!/usr/bin/env python
"""Build every pairing arm of the BMMC ladder that does not need Seurat, and measure it."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "benchmarks/bmmc" / "out"
OUT.mkdir(parents=True, exist_ok=True)

# SCENIC+'s rule inside prepare_GEX_ACC: nr_metacells = round(n_cells / k) * 2 per label. The
# doubling makes a metacell arm reuse each cell ~2x rather than partition the cells. Exposed so
# the partition variant (1.0) is reachable as a sensitivity analysis.
METACELL_OVERSAMPLE = 2.0

# The crossing axes every arm is scored on. `donor` is the raw site-local DonorNumber, kept only
# because the census used it and the comparison must be like-for-like; `subject` is the person;
# `cell_type` is included so a label-respecting arm can be shown to respect it.
CROSS_AXES = ("donor", "subject", "site", "cell_type")

ARMS = ["truth", "rand_celltype_1to1", "rand_celltype_donor_1to1",
        "rand_celltype_batch_1to1", "rand_celltype_matched", "rand_celltype_donor_matched",
        "rand_celltype_donor_matched_cap14",
        "rand_global", "pseudomultiome_k5", "pseudomultiome_k10"]

# accepted column spellings in a --reuse-from table, so the Seurat-side builder can use its
# own column names without editing this script
ATAC_ALIASES = ("atac_barcode", "h5ad_barcode", "cistopic_id", "cell_id", "atac_id")
RNA_ALIASES = ("assigned_rna_barcode", "Linked_scRNA_barcode", "rna_barcode", "linked_rna_barcode")

def log(m): print(f"[b02] {m}", flush=True)

def require(cond, msg):
    """assert, but immune to `python -O`. Every use of this is a thing that should be impossible;
    none of them are input validation, which raises SystemExit instead."""
    if not cond:
        raise AssertionError(msg)

def rng_for(seed: int, arm: str) -> np.random.Generator:
    """One independent stream per arm. hashlib, not hash(): hash() is salted per process, so a
    hash()-derived stream would not reproduce across runs even at a fixed seed."""
    h = int.from_bytes(hashlib.sha256(arm.encode()).digest()[:8], "big")
    return np.random.default_rng(np.random.SeedSequence([int(seed), h]))

def groups(stratum) -> list:
    """[(label, positions), ...] in sorted label order. Built with one factorize rather than a
    boolean mask per stratum (O(n_strata * n), the dominant cost at 247 strata x 69k cells). Sorted
    so iteration order, and therefore every draw inside it, is reproducible."""
    s = pd.Series(np.asarray(stratum, dtype=object))
    return [(k, np.asarray(v, dtype=np.int64)) for k, v in sorted(s.groupby(s).indices.items())]

def ks_distance(a, b) -> float:
    """Two-sample KS statistic between two multiplicity vectors. For the reuse-matched arm this is
    0.0 by construction -- the achieved vector is the target vector permuted. Anything else means the
    draw or the assignment is wrong, so it is asserted, not merely logged."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    grid = np.union1d(a, b)
    ca = np.searchsorted(a, grid, side="right") / a.size
    cb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))

# ------------------------------------------------------------------ labels ----
def check_labels(md: pd.DataFrame) -> dict:
    """Establish, from the data, the two label facts the arm definitions depend on:
    (1) the cell key is unique: fragment_barcode repeats across batches, so h5ad_barcode is the
        only safe identity (the fragment barcode as a key would silently merge cells);
    (2) whether DonorNumber is a global subject id. It is: the same person is sampled at more
        than one site, confirmed from donor-constant demographics identical across a donor's
        sites. That replication makes site separable from donor.
    Both are measured and returned for the report."""
    n = len(md)
    require(md["h5ad_barcode"].is_unique, "h5ad_barcode is not unique -- it cannot be the cell key")
    frag_unique = bool(md["fragment_barcode"].is_unique)

    sd = md.groupby(["site", "donor"], sort=True).size()
    donor_sites = md.groupby("donor")["site"].nunique()
    site_local = bool((donor_sites > 1).any())
    # subject == batch is a claim, so check it in both directions rather than assuming the naming
    per_batch = md.groupby("batch")[["site", "donor"]].nunique()
    batch_is_subject = bool((per_batch["site"] == 1).all() and (per_batch["donor"] == 1).all()
                            and md.groupby(["site", "donor"])["batch"].nunique().eq(1).all())
    log(f"labels: {n:,} cells | h5ad_barcode unique=True | fragment_barcode unique={frag_unique}"
        + ("" if frag_unique else " (collides across batches -- NOT usable as the cell key; "
                                 "h5ad_barcode is)"))
    log(f"labels: {md['donor'].nunique()} DonorNumbers but {len(sd)} (site,donor) subjects; "
        f"DonorNumber is site-local={site_local}; subject==batch={batch_is_subject}")
    if site_local:
        multi = {str(d): int(v) for d, v in donor_sites.items() if v > 1}
        log(f"labels:   donors sampled at >1 site: {multi}  <- the replication that makes the "
            f"pairing x site effect identifiable separately from donor")
    require(batch_is_subject,
            "batch is not 1:1 with (site,donor) -- the sample-aware arm's key does not hold")
    require(site_local,
            "no donor appears at more than one site: the site axis is NOT identifiable on this "
            "metadata, and rand_celltype_donor_1to1 would collapse onto the sample-aware arm")
    return {"n_cells": int(n), "fragment_barcode_unique": frag_unique,
            "n_donor_labels": int(md["donor"].nunique()), "n_subjects": int(len(sd)),
            "n_sites": int(md["site"].nunique()), "n_batches": int(md["batch"].nunique()),
            "n_cell_types": int(md["cell_type"].nunique()),
            "donor_label_is_site_local": site_local,
            "subject_equals_batch": batch_is_subject,
            "subject_key": "(site, donor)"}

def analytic_cross(md: pd.DataFrame, stratum: np.ndarray, axis: str,
                   weights: np.ndarray | None = None) -> float:
    """E[share of assignments landing on a different value of `axis`] when the RNA partner is drawn
    uniformly from the ATAC cell's own stratum: mean over cells of (1 - sum_a p_a^2) inside the
    stratum (the census formula, recomputed here). `weights` lets the metacell arms weight by slots."""
    val = md[axis].to_numpy()
    w = np.ones(len(md), float) if weights is None else np.asarray(weights, float)
    tot = w.sum()
    if tot == 0:
        return float("nan")
    acc = 0.0
    for _, idx in groups(stratum):
        ws = w[idx].sum()
        if ws == 0:
            continue
        p = np.array(list(Counter(val[idx].tolist()).values()), float) / idx.size
        acc += (ws / tot) * (1.0 - float((p ** 2).sum()))
    return float(acc)

# ----------------------------------------------------------------- samplers ----
def bijection_arm(md: pd.DataFrame, stratum: np.ndarray, rng: np.random.Generator):
    """Uniform random permutation of RNA cells over ATAC cells within each stratum.

    Feasible unconditionally because BMMC is a true multiome: the RNA pool of a stratum is exactly
    the set of ATAC cells in it, so n_RNA == n_ATAC and there is nothing to ration. Returns
    positional index arrays (atac_idx, rna_idx) into `md`."""
    n = len(md)
    rna_idx = np.empty(n, dtype=np.int64)
    for _, idx in groups(stratum):
        rna_idx[idx] = idx[rng.permutation(idx.size)]
    return np.arange(n, dtype=np.int64), rna_idx

def matched_arm(md: pd.DataFrame, stratum: np.ndarray, ref_rna: np.ndarray,
                rng: np.random.Generator):
    """Random within stratum, reproducing the reference arm's per-stratum reuse concentration.

    Per stratum s, read from the reference only d_s (distinct RNA cells used) and M_s (multiplicity
    vector). Discard its choices. Draw d_s fresh cells without replacement from the stratum's own
    pool, shuffle M_s onto them, shuffle the ATAC cells. Returns (atac_idx, rna_idx, per_stratum)."""
    n = len(md)
    rna_idx = np.empty(n, dtype=np.int64)
    per_stratum = {}
    for s, idx in groups(stratum):
        n_s = idx.size
        counts = pd.Series(ref_rna[idx]).value_counts()
        d_s, M_s = int(counts.size), counts.to_numpy()
        # d_s > n_s means the reference used more distinct RNA cells than there are ATAC cells in
        # the stratum, i.e. the two sides disagree about the labels. Not recoverable here.
        require(d_s <= n_s, f"stratum {s}: reference used d_s={d_s} distinct RNA cells for "
                            f"n_s={n_s} ATAC cells -- labels inconsistent between modalities")
        require(int(M_s.sum()) == n_s, f"stratum {s}: multiplicities sum to {M_s.sum()}, want {n_s}")
        picks = rng.choice(n_s, size=d_s, replace=False)          # positions within the pool
        donors = idx[picks]
        M_shuf = M_s.copy()
        rng.shuffle(M_shuf)
        cells = idx[rng.permutation(n_s)]
        k = 0
        for donor, m in zip(donors, M_shuf):
            rna_idx[cells[k:k + m]] = donor
            k += int(m)
        require(k == n_s, f"stratum {s}: assigned {k} of {n_s} cells")
        achieved = np.unique(rna_idx[idx], return_counts=True)[1]
        ks = ks_distance(M_s, achieved)
        require(ks == 0.0, f"stratum {s}: achieved multiplicity is not a permutation of the "
                           f"target (KS={ks:.4f}) -- the sampler is wrong")
        per_stratum[str(s)] = {"n_s": int(n_s), "d_s_target": d_s,
                               "d_s_achieved": int(achieved.size),
                               "reuse_target": float(n_s / d_s),
                               "reuse_achieved": float(n_s / achieved.size),
                               "ks_multiplicity": ks}
    return np.arange(n, dtype=np.int64), rna_idx, per_stratum

def metacell_arm(md: pd.DataFrame, k: int, rng: np.random.Generator, oversample: float):
    """SCENIC+-style metacells: per cell-type label draw round(n/k)*oversample groups, each a
    without-replacement sample of k cells, drawn independently for the two modalities.

    A label holding fewer than k cells raises `Sample larger than population` inside SCENIC+ and
    kills the whole arm, so such labels are dropped here, named in the log and recorded in the
    report, never silently sampled with replacement."""
    a_idx, r_idx, mc_id, mc_ct = [], [], [], []
    contained, n_slot_by_ct = [], Counter()
    skipped = {}
    for c, idx in groups(md["cell_type"].to_numpy()):
        n_c = idx.size
        if n_c < k:
            skipped[str(c)] = int(n_c)
            continue
        n_meta = max(1, int(round(round(n_c / k) * oversample)))
        for m in range(n_meta):
            a = idx[rng.choice(n_c, size=k, replace=False)]
            r = idx[rng.choice(n_c, size=k, replace=False)]
            require(len(set(a.tolist())) == k and len(set(r.tolist())) == k,
                    f"{c} metacell {m}: members repeat -- sampling was not without replacement")
            rset = set(r.tolist())
            a_idx.append(a); r_idx.append(r)
            mc_id.extend([f"{c}|mc{m:05d}"] * k)
            mc_ct.extend([c] * k)
            contained.extend([int(x) in rset for x in a])
        n_slot_by_ct[c] += n_meta * k
    if skipped:
        log(f"  {len(skipped)} cell type(s) below the k={k} metacell floor, DROPPED from this arm "
            f"(they would raise 'Sample larger than population' inside SCENIC+): {skipped}")
    return (np.concatenate(a_idx), np.concatenate(r_idx), np.array(mc_id), np.array(mc_ct),
            np.array(contained, bool), skipped, n_slot_by_ct)

# --------------------------------------------------------------- measurement ----
def measure(md: pd.DataFrame, atac_idx, rna_idx, expected: dict) -> dict:
    """The central measurement, identical for every arm so the arms are actually comparable."""
    # a missing expectation would come back as NaN, which downstream checks skip, so a check
    # could pass while measuring nothing. Require the keys instead of defaulting them.
    for key in [f"cross_{ax}" for ax in CROSS_AXES] + ["self_hits"]:
        require(key in expected, f"no analytic expectation supplied for '{key}'; a measurement "
                                 f"with nothing to compare against is not a check")
    n = int(len(atac_idx))
    used, counts = np.unique(rna_idx, return_counts=True)
    d = int(used.size)
    hist = Counter(counts.tolist())
    self_hits = int((atac_idx == rna_idx).sum())
    out = {
        "n_assignments": n,
        "reuse_mean_per_used_rna_cell": float(n / d),
        "reuse_mean_per_assignment": float((counts.astype(float) ** 2).sum() / n),
        "reuse_max": int(counts.max()),
        "reuse_distribution": {str(m): int(c) for m, c in sorted(hist.items())},
        "share_of_assignments_on_a_reused_cell": float(counts[counts > 1].sum() / n),
        "distinct_rna_cells": d,
        "distinct_rna_share_of_pool": float(d / len(md)),
        "self_hits": self_hits,
        "self_hit_rate": float(self_hits / n),
        "self_hits_expected": float(expected.get("self_hits", float("nan"))),
        "self_hit_rate_expected": float(expected.get("self_hits", float("nan")) / n),
    }
    for ax in CROSS_AXES:
        v = md[ax].to_numpy()
        out[f"cross_{ax}_rate"] = float((v[atac_idx] != v[rna_idx]).mean())
        out[f"cross_{ax}_rate_expected"] = float(expected.get(f"cross_{ax}", float("nan")))
    return out

def tol(p: float, n: int) -> float:
    """Six standard errors of a rate, floored at 2pp so a small smoke-test frame cannot fail
    spuriously. Assignments inside a stratum are negatively correlated under a permutation, so the
    binomial SE overstates the true spread -- this band is conservative in the safe direction."""
    p = min(max(p, 0.0), 1.0)
    return max(0.02, 6.0 * float(np.sqrt(max(p * (1 - p), 1e-12) / max(n, 1))))

def check_arm(name: str, kind: str, atac_idx, rna_idx, stratum, stats: dict):
    """Every check below exercises the property it claims, on the realised assignment.

    In particular the bijection check compares the multiset of assigned RNA cells against the
    multiset of the stratum's own cells, per stratum -- which establishes 'is a permutation' and
    'stayed inside the stratum' simultaneously and exactly. Checking reuse == 1.0 alone would not:
    an assignment can be injective and still have walked out of its stratum."""
    n = stats["n_assignments"]

    if kind == "bijection":
        for s, idx in groups(stratum):
            require(np.array_equal(np.sort(rna_idx[idx]), np.sort(atac_idx[idx])),
                    f"{name}: stratum {s} is NOT a permutation of its own cells -- the arm claims a "
                    f"within-stratum bijection and does not have one")
        require(stats["reuse_mean_per_used_rna_cell"] == 1.0,
                f"{name}: claims a bijection but reuse is "
                f"{stats['reuse_mean_per_used_rna_cell']}")
        require(stats["distinct_rna_cells"] == n,
                f"{name}: claims a bijection but used {stats['distinct_rna_cells']} of {n} cells")

    if name == "truth":
        require(stats["self_hit_rate"] == 1.0, "truth: self-hit rate is not 1.0")
        for ax in CROSS_AXES:
            require(stats[f"cross_{ax}_rate"] == 0.0, f"truth: crosses {ax} -- identity is broken")

    if name == "rand_celltype_donor_1to1":
        # person-aware: same donor, free to cross site. cross_site == 0 is not asserted: donor1 spans
        # three sites and donor3 two, so a correct person-aware assignment must cross sites for
        # those donors. That site crossing is the quantity this dataset is used to measure.
        for ax in ("donor", "subject"):
            require(stats[f"cross_{ax}_rate"] == 0.0,
                    f"{name}: PERSON-AWARE arm has cross_{ax}_rate="
                    f"{stats[f'cross_{ax}_rate']:.6f}, must be exactly 0")
        require(stats["cross_site_rate"] > 0.0,
                f"{name}: cross_site_rate is 0 -- with donor1 at 3 sites and donor3 at 2, a "
                f"person-aware arm that never crosses a site has silently collapsed onto the "
                f"sample-aware arm and the site axis is unmeasurable")

    if name == "rand_celltype_batch_1to1":
        # sample-aware: the strictest control; crosses nothing on any axis.
        for ax in ("donor", "subject", "site"):
            require(stats[f"cross_{ax}_rate"] == 0.0,
                    f"{name}: SAMPLE-AWARE arm has cross_{ax}_rate="
                    f"{stats[f'cross_{ax}_rate']:.6f}, must be exactly 0")

    if name != "rand_global" and name != "truth":
        require(stats["cross_cell_type_rate"] == 0.0,
                f"{name}: crosses cell type -- every arm but rand_global is label-respecting")

    # self-hits are checked on the count, not the rate: the expected count is small (one per
    # stratum: 22 for the cell-type arms, 1 for rand_global) and a rate tolerance floored at 2pp
    # would pass a thousand spurious fixed points. The number of fixed points of a random
    # permutation is asymptotically Poisson(1) per stratum and sums over independent strata, so the
    # band is 5 sigma plus a 5-count slack for small means.
    exp_hits = stats["self_hits_expected"]
    if np.isfinite(exp_hits):
        band = 5.0 * float(np.sqrt(max(exp_hits, 1.0))) + 5.0
        require(abs(stats["self_hits"] - exp_hits) <= band,
                f"{name}: {stats['self_hits']} self-hits but the arm's own randomisation predicts "
                f"{exp_hits:.1f} (band +/-{band:.1f}) -- the true partner is being drawn at the "
                f"wrong rate, which is a sampler bug, not noise")

    # crossing rates, measured vs the arm's own analytic expectation
    for ax in CROSS_AXES:
        key = f"cross_{ax}_rate"
        exp = stats.get(f"{key}_expected", float("nan"))
        if not np.isfinite(exp):
            continue
        t = tol(exp, n)
        require(abs(stats[key] - exp) <= t,
                f"{name}: {key}={stats[key]:.6f} but analytic expectation is {exp:.6f} "
                f"(|diff| {abs(stats[key]-exp):.6f} > tolerance {t:.6f}) -- sampler bug, not noise")

# ---------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", default=str(OUT / "b00_cell_metadata.csv"),
                    help="b_00's per-cell table (default: the one b_00 writes)")
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--outdir", default=str(OUT))
    ap.add_argument("--reuse-from", default=None,
                    help="assignment table whose per-stratum reuse distribution rand_celltype_"
                         "matched should reproduce -- normally the `link` arm, built with Seurat "
                         "elsewhere. Needs one ATAC-cell column and one assigned-RNA-cell column, "
                         "both holding h5ad_barcode values. Absent => that ONE arm is skipped.")
    ap.add_argument("--metacell-oversample", type=float, default=METACELL_OVERSAMPLE,
                    help="metacells per label = round(n/k) * this. Default 2.0 is SCENIC+'s own "
                         "rule; 1.0 gives the partition variant.")
    a = ap.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    md_fp = Path(a.metadata)
    if not md_fp.is_file():
        raise SystemExit(f"[b02] MISSING {md_fp} -- run census.py first")
    md = pd.read_csv(md_fp, dtype=str)
    for col in ("h5ad_barcode", "fragment_barcode", "batch", "donor", "site", "cell_type"):
        if col not in md.columns:
            raise SystemExit(f"[b02] {md_fp.name} has no '{col}' column; got {list(md.columns)}")
    ## subject = the person. The sample is `batch`; keeping them distinct makes the site axis
    ## measurable (see the docstring).
    md["subject"] = md["donor"]
    md["sample"]  = md["batch"]
    N = len(md)
    log(f"source: {md_fp} -> {N:,} cells")
    labels = check_labels(md)

    # ---- census cross-check ---------------------------------------------------------
    # The census computed the cell-type-only crossing ceiling on the raw DonorNumber label. Recompute
    # it here from the metadata and assert the two agree, so the comparisons below rest on a number
    # this script derived.
    st_ct = md["cell_type"].to_numpy()
    st_sub = (md["cell_type"] + "|" + md["subject"]).to_numpy()   # person-aware
    st_bat = (md["cell_type"] + "|" + md["sample"]).to_numpy()    # sample-aware
    st_all = np.full(N, "all", dtype=object)
    n_by_ct = {c: int(idx.size) for c, idx in groups(st_ct)}
    ceiling = {ax: analytic_cross(md, st_ct, ax) for ax in CROSS_AXES}
    log(f"analytic ceiling (recomputed here, cell-type-only pairing): "
        f"cross-donor {ceiling['donor']:.2%} | cross-subject {ceiling['subject']:.2%} | "
        f"cross-site {ceiling['site']:.2%}")
    census_fp = outdir / "b00_feasibility.json"
    census = None
    if census_fp.is_file():
        c = json.loads(census_fp.read_text(encoding="utf-8"))
        census = {"cross_donor": c.get("expected_cross_donor_rate_celltype_only"),
                  "cross_site": c.get("expected_cross_site_rate_celltype_only"),
                  "n_cells": c.get("n_cells"), "file": census_fp.name}
        log(f"census reference ({census_fp.name}): cross-donor {census['cross_donor']:.2%} | "
            f"cross-site {census['cross_site']:.2%}")
        require(abs(ceiling["donor"] - census["cross_donor"]) < 1e-9
                and abs(ceiling["site"] - census["cross_site"]) < 1e-9,
                f"recomputed ceiling ({ceiling['donor']:.6f}/{ceiling['site']:.6f}) disagrees with "
                f"b_00's census ({census['cross_donor']:.6f}/{census['cross_site']:.6f}) -- the two "
                f"scripts are not looking at the same metadata")
        log("census cross-check: recomputed ceiling matches b_00 exactly")
    else:
        log(f"census reference {census_fp.name} ABSENT -- ceiling used is this script's own "
            f"recomputation, uncorroborated")

    # ---- the reuse reference, if it exists -------------------------------------------
    ref_rna, ref_meta = None, {"supplied": False, "file": a.reuse_from}
    if a.reuse_from:
        rf = Path(a.reuse_from)
        if not rf.is_file():
            raise SystemExit(f"[b02] --reuse-from {rf} does not exist")
        R = pd.read_csv(rf, dtype=str)
        ac = next((c for c in ATAC_ALIASES if c in R.columns), None)
        rc = next((c for c in RNA_ALIASES if c in R.columns), None)
        if ac is None or rc is None:
            raise SystemExit(f"[b02] {rf.name} has no recognisable ATAC/RNA columns; got "
                             f"{list(R.columns)}; expected one of {ATAC_ALIASES} and {RNA_ALIASES}")

        if len(R) != N or set(R[ac]) != set(md["h5ad_barcode"]):
            raise SystemExit(
                f"[b02] {rf.name} covers {len(R):,} rows / {R[ac].nunique():,} ATAC cells but the "
                f"analysis set is {N:,}. The reuse-matched arm must cover the SAME cells as every "
                f"other arm. Rebuild the reference on the full manifest -- do not rescale here.")
        unknown = set(R[rc]) - set(md["h5ad_barcode"])
        if unknown:
            raise SystemExit(f"[b02] {rf.name} assigns {len(unknown):,} RNA barcodes that are not "
                             f"cells of this dataset, e.g. {sorted(unknown)[:3]}")
        pos = pd.Series(np.arange(N), index=md["h5ad_barcode"])
        ref_rna = R.set_index(ac).loc[md["h5ad_barcode"], rc].map(pos).to_numpy(dtype=np.int64)
        ref_d = int(np.unique(ref_rna).size)
        ref_meta = {"supplied": True, "file": str(rf), "arm": rf.stem,
                    "reference_reuse": float(N / ref_d), "reference_distinct_rna_cells": ref_d,
                    "reference_self_hit_rate": float((ref_rna == np.arange(N)).mean())}
        log(f"reuse reference: {rf.name} | reuse {N/ref_d:.3f}x | {ref_d:,} distinct RNA cells "
            f"({ref_d/N:.1%} of pool) | its own self-hit rate {ref_meta['reference_self_hit_rate']:.4%}")
    else:
        log("reuse reference: NOT SUPPLIED (--reuse-from). The reuse-matched arms will be SKIPPED; "
            "rerun with --reuse-from once the link arm exists. No other arm depends on them.")

    # ---- build ------------------------------------------------------------------------
    report = {
        "script": "arms.py",
        "dataset": "GSE194122 BMMC 10x Multiome (NeurIPS 2021)",
        "seed": a.seed, "metadata_source": str(md_fp), "metacell_oversample": a.metacell_oversample,
        "labels": labels,
        "true_pairing_note": "true multiome: the RNA partner of ATAC cell i IS cell i, so every arm "
                             "is a permutation (or, for the metacell arms, a resampling) of the "
                             "identity and the self-hit rate needs no join to compute",
        "analytic_ceiling_celltype_only": ceiling,
        "census_reference": census,
        "reuse_reference": ref_meta,
        "arms_built": [], "arms_skipped": {}, "arms": {},
    }

    plan = [
        ("truth", "identity", st_ct, "none (identity)"),
        ("rand_celltype_1to1", "bijection", st_ct, "cell_type"),
        ("rand_celltype_donor_1to1", "bijection", st_sub, "cell_type x donor (person; may cross site)"),
        ("rand_celltype_batch_1to1", "bijection", st_bat, "cell_type x batch (sample; crosses nothing)"),
        ("rand_celltype_matched", "matched", st_ct, "cell_type"),
        ("rand_celltype_donor_matched", "matched", st_sub, "cell_type x donor"),
        ("rand_celltype_donor_matched_cap14", "matched", st_sub, "cell_type x donor", (1.35, 1.45)),
        ("rand_global", "bijection", st_all, "none (global)"),
    ]

    for row in plan:
        name, kind, stratum, sdesc = row[:4]
        ref_band = row[4] if len(row) > 4 else None
        # keyed on kind, not on a name list: every matched arm reproduces a supplied reuse
        # distribution, so none can be built without --reuse-from
        if kind == "matched" and ref_rna is None:
            report["arms_skipped"][name] = ("no --reuse-from supplied; this arm exists only to "
                                            "reproduce the link arm's measured reuse distribution")
            log(f"{name:26s} SKIPPED -- no reuse distribution to match")
            continue
        if ref_band is not None:
            ref_reuse = ref_meta["reference_reuse"]
            if not (ref_band[0] <= ref_reuse <= ref_band[1]):
                report["arms_skipped"][name] = (
                    f"--reuse-from reference has measured reuse {ref_reuse:.4f}x, outside this "
                    f"arm's required band {list(ref_band)}; this arm exists to match a reference "
                    f"inside that band and would otherwise copy the wrong concentration")
                log(f"{name:26s} SKIPPED -- reference reuse {ref_reuse:.4f}x is outside the band "
                    f"{list(ref_band)} this arm is defined by")
                continue
        rng = rng_for(a.seed, name)
        per_stratum = None
        if kind == "identity":
            atac_idx = np.arange(N, dtype=np.int64)
            rna_idx = atac_idx.copy()
        elif kind == "bijection":
            atac_idx, rna_idx = bijection_arm(md, stratum, rng)
        else:
            atac_idx, rna_idx, per_stratum = matched_arm(md, stratum, ref_rna, rng)

        n_strata = int(len(set(stratum.tolist())))
        # keys must be the same strings measure() reads, or the comparison silently degrades to
        # NaN and every crossing assertion below becomes a no-op
        exp = {f"cross_{ax}": (0.0 if kind == "identity" else analytic_cross(md, stratum, ax))
               for ax in CROSS_AXES}
        # exactly one expected fixed point per stratum, for any stratum size and any d_s
        exp["self_hits"] = float(N) if kind == "identity" else float(n_strata)
        stats = measure(md, atac_idx, rna_idx, exp)
        check_arm(name, kind, atac_idx, rna_idx, stratum, stats)

        frame = pd.DataFrame({
            "atac_barcode": md["h5ad_barcode"].to_numpy(),
            "atac_fragment_barcode": md["fragment_barcode"].to_numpy(),
            "assigned_rna_barcode": md["h5ad_barcode"].to_numpy()[rna_idx],
            "true_rna_barcode": md["h5ad_barcode"].to_numpy(),
            "stratum": stratum.astype(str),
            "atac_donor": md["donor"].to_numpy(), "rna_donor": md["donor"].to_numpy()[rna_idx],
            "atac_site": md["site"].to_numpy(), "rna_site": md["site"].to_numpy()[rna_idx],
            "atac_cell_type": md["cell_type"].to_numpy(),
            "is_self": atac_idx == rna_idx,
            "cross_donor": md["donor"].to_numpy() != md["donor"].to_numpy()[rna_idx],
            "cross_subject": md["subject"].to_numpy() != md["subject"].to_numpy()[rna_idx],
            "cross_site": md["site"].to_numpy() != md["site"].to_numpy()[rna_idx],
        })
        dest = outdir / f"arm_{name}_seed{a.seed}.csv"
        frame.to_csv(dest, index=False, lineterminator="\n")

        entry = {"kind": kind, "stratum_definition": sdesc, "n_strata": n_strata,
                 "row_unit": "one row per ATAC cell", "file": dest.name, **stats}
        if per_stratum is not None:
            entry["per_stratum"] = per_stratum
        report["arms"][name] = entry
        report["arms_built"].append(name)
        log(f"{name:26s} strata={n_strata:4d} | reuse {stats['reuse_mean_per_used_rna_cell']:.3f}x "
            f"(per-assignment {stats['reuse_mean_per_assignment']:.3f}x, max "
            f"{stats['reuse_max']}) | distinct {stats['distinct_rna_cells']:,} "
            f"({stats['distinct_rna_share_of_pool']:.1%} of pool)")
        log(f"{'':26s} self-hits {stats['self_hits']:,} = {stats['self_hit_rate']:.4%} "
            f"(analytic {stats['self_hits_expected']:.0f} = "
            f"{stats['self_hit_rate_expected']:.4%}) | CROSS-DONOR "
            f"{stats['cross_donor_rate']:.2%} (analytic {stats['cross_donor_rate_expected']:.2%}) "
            f"| CROSS-SUBJECT {stats['cross_subject_rate']:.2%} | CROSS-SITE "
            f"{stats['cross_site_rate']:.2%} (analytic {stats['cross_site_rate_expected']:.2%})")

    # ---- metacell arms ------------------------------------------------------------------
    for k in (5, 10):
        name = f"pseudomultiome_k{k}"
        rng = rng_for(a.seed, name)
        a_idx, r_idx, mc_id, mc_ct, contained, skipped, slots_by_ct = metacell_arm(
            md, k, rng, a.metacell_oversample)
        # slot weights: the crossing expectation is per slot, and a cell type's slot count is
        # round(n_c/k)*oversample*k, which is not proportional to n_c (rounding, the k-floor and the
        # doubling bend it). Weighting by cells instead of slots would answer a different question.
        w = np.array([slots_by_ct.get(c, 0) / n_by_ct[c] for c in md["cell_type"]], dtype=float)
        exp = {f"cross_{ax}": analytic_cross(md, st_ct, ax, weights=w) for ax in CROSS_AXES}
        n_slots = int(a_idx.size)
        # per slot the two modalities draw independently and uniformly over the label, so
        # P(same cell) = 1/n_c; summed over the label's slots that is slots_c / n_c
        exp["self_hits"] = float(sum(s / n_by_ct[c] for c, s in slots_by_ct.items()))
        stats = measure(md, a_idx, r_idx, exp)
        check_arm(name, "metacell", a_idx, r_idx, mc_ct, stats)

        n_meta = int(pd.unique(mc_id).size)
        cont_rate = float(contained.mean())
        cont_exp = float(sum(s * k / n_by_ct[c] for c, s in slots_by_ct.items()) / n_slots)
        require(abs(cont_rate - cont_exp) <= tol(cont_exp, n_slots),
                f"{name}: containment rate {cont_rate:.6f} vs analytic {cont_exp:.6f} -- the two "
                f"modalities are not being drawn independently")

        frame = pd.DataFrame({
            "metacell_id": mc_id, "cell_type": mc_ct,
            "slot": np.concatenate([np.arange(k)] * n_meta),
            "atac_barcode": md["h5ad_barcode"].to_numpy()[a_idx],
            "atac_fragment_barcode": md["fragment_barcode"].to_numpy()[a_idx],
            "rna_barcode": md["h5ad_barcode"].to_numpy()[r_idx],
            "atac_donor": md["donor"].to_numpy()[a_idx], "rna_donor": md["donor"].to_numpy()[r_idx],
            "atac_site": md["site"].to_numpy()[a_idx], "rna_site": md["site"].to_numpy()[r_idx],
            "is_self": a_idx == r_idx,
            "own_rna_in_own_metacell": contained,
            "cross_donor": md["donor"].to_numpy()[a_idx] != md["donor"].to_numpy()[r_idx],
            "cross_subject": md["subject"].to_numpy()[a_idx] != md["subject"].to_numpy()[r_idx],
            "cross_site": md["site"].to_numpy()[a_idx] != md["site"].to_numpy()[r_idx],
        })
        dest = outdir / f"arm_{name}_seed{a.seed}.csv"
        frame.to_csv(dest, index=False, lineterminator="\n")

        report["arms"][name] = {
            "kind": "metacell", "stratum_definition": "cell_type", "k": k,
            "n_metacells": n_meta, "n_strata": int(len(slots_by_ct)),
            "row_unit": "one row per (metacell, slot); the metacell is the unit and the ATAC/RNA "
                        "sides are drawn independently -- the slot pairing is bookkeeping only",
            "cell_types_below_k_floor": skipped,
            "own_rna_in_own_metacell_rate": cont_rate,
            "own_rna_in_own_metacell_rate_expected": cont_exp,
            "file": dest.name, **stats}
        report["arms_built"].append(name)
        log(f"{name:26s} k={k} metacells={n_meta:,} slots={n_slots:,} | reuse "
            f"{stats['reuse_mean_per_used_rna_cell']:.3f}x (per-assignment "
            f"{stats['reuse_mean_per_assignment']:.3f}x, max {stats['reuse_max']}) | distinct "
            f"{stats['distinct_rna_cells']:,} ({stats['distinct_rna_share_of_pool']:.1%} of pool)")
        log(f"{'':26s} slot self-hits {stats['self_hits']:,} = {stats['self_hit_rate']:.4%} "
            f"(analytic {stats['self_hit_rate_expected']:.4%}); own RNA anywhere in own metacell "
            f"{cont_rate:.4%} (analytic {cont_exp:.4%}) | CROSS-DONOR {stats['cross_donor_rate']:.2%}"
            f" (analytic {stats['cross_donor_rate_expected']:.2%}) | CROSS-SITE "
            f"{stats['cross_site_rate']:.2%} (analytic {stats['cross_site_rate_expected']:.2%})")

    # ---- the headline, echoed with its source ---------------------------------------------
    blind = report["arms"].get("rand_celltype_1to1")
    aware = report["arms"].get("rand_celltype_donor_1to1")
    if blind and aware:
        log("")
        log(f"HEADLINE (b02_arm_report.json -> arms.rand_celltype_1to1): the DONOR-BLIND bijection "
            f"gives {blind['cross_donor_rate']:.1%} of cells an RNA profile from a different "
            f"DonorNumber and {blind['cross_site_rate']:.1%} from a different SITE"
            + (f"; b_00's census predicted {census['cross_donor']:.1%} / {census['cross_site']:.1%}"
               if census else "; no census file to compare against"))
        log(f"HEADLINE (b02_arm_report.json -> arms.rand_celltype_donor_1to1): the PERSON-AWARE "
            f"bijection gives {aware['cross_donor_rate']:.1%} cross-donor (asserted exactly 0) "
            f"but {aware['cross_site_rate']:.1%} cross-SITE -- it is free to cross sites because "
            f"donor1 spans three and donor3 two, and that residual is the site effect priced "
            f"separately from donor identity")
        report["headline"] = {
            "donor_blind_arm": "rand_celltype_1to1",
            "donor_blind_cross_donor_rate": blind["cross_donor_rate"],
            "donor_blind_cross_subject_rate": blind["cross_subject_rate"],
            "donor_blind_cross_site_rate": blind["cross_site_rate"],
            "census_predicted_cross_donor_rate": (census or {}).get("cross_donor"),
            "census_predicted_cross_site_rate": (census or {}).get("cross_site"),
            "donor_aware_arm": "rand_celltype_donor_1to1",
            "donor_aware_cross_donor_rate": aware["cross_donor_rate"],
            "donor_aware_cross_site_rate": aware["cross_site_rate"],
            "interpretation": "the donor-blind arm is the pairing a cell-type-only method builds; "
                              "the donor-aware arm is the same construction with the batch "
                              "constraint added. PBMC 10k (one donor) cannot separate these two.",
        }

    rep_fp = outdir / "b02_arm_report.json"
    rep_fp.write_text(json.dumps(report, indent=1), encoding="utf-8")
    log(f"wrote {len(report['arms_built'])} arm tables + {rep_fp.name} -> {outdir}")
    if report["arms_skipped"]:
        log(f"SKIPPED: {sorted(report['arms_skipped'])} (see arms_skipped in the report)")

if __name__ == "__main__":
    main()
