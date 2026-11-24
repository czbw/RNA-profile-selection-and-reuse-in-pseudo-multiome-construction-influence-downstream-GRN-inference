#!/usr/bin/env python
"""Shared loaders, the premise check, and the two fixed declarations for the regulon-activity
scoring scripts (b_06 to b_09)."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent

STAGE_DEFAULT = "bmmc"
SEED_DEFAULT = 555
EXPECTED_N_CELLS = 69249          # the cell set every BMMC arm is a bijection over (asserted in b_05)
TRUTH = "truth"

LADDER_ARMS = ["truth", "rand_celltype_1to1", "rand_celltype_x_donor_1to1", "rand_global"]
# Built for BMMC but never carried into the GRN half. Scoring here needs no eGRN, so they are
# scoreable; reported only with ladder=false beside them. See declarations().
NONLADDER_ARMS = ["link_1to1", "link",

                  "link_x_celltype_donor_1to1", "link_x_celltype_donor",
                  "rand_celltype_donor_matched",
                  "link_x_celltype_donor_cap14", "rand_celltype_donor_matched_cap14"]

DONOR_BLIND = "rand_celltype_1to1"
DONOR_AWARE = "rand_celltype_x_donor_1to1"

# arm token -> extra file stems to try after the canonical one. Same list b_05 uses: the builder
# (`arms.py` line 155) names the donor-aware arm `rand_celltype_donor_1to1` while
# `submit_bmmc.sh` / `config.yaml` name it `rand_celltype_x_donor_1to1`.
ARM_TABLE_ALIASES = {
    "rand_celltype_x_donor_1to1": ["rand_celltype_donor_1to1"],
    "rand_celltype_donor_1to1": ["rand_celltype_x_donor_1to1"],
}

MIN_CELLS_PER_TYPE = 20           # p_08/p_09/p_10 all use this; kept identical

def make_log(tag):
    def log(m):
        print(f"[{tag}] {m}", flush=True)
    return log

def die(tag, msg):
    print(f"[{tag}] FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)

def sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()

def peak_rss_bytes() -> int | None:
    """Peak resident set size of this process, or None if it cannot be determined."""
    try:
        import psutil
        mi = psutil.Process().memory_info()
        return int(getattr(mi, "peak_wset", mi.rss))     # peak_wset exists on Windows
    except Exception:
        return None

def library_versions() -> dict:
    from scipy import __version__ as scipy_version
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy_version, "h5py": h5py.__version__}

# loaders

def resolve_table(tag, uploads: Path, stage: str, seed: int, arm: str) -> Path:
    """Canonical name first, then the declared aliases. Fails naming every path it tried."""
    stems = [arm] + ARM_TABLE_ALIASES.get(arm, [])
    tried = []
    for stem in stems:
        for name in (f"arm_{stem}_{stage}_seed{seed}_full.csv", f"arm_{stem}_seed{seed}.csv"):
            p = uploads / name
            tried.append(p)
            if p.is_file() and p.stat().st_size > 0:
                return p
    die(tag, "no assignment table for arm {!r}. Tried:\n       {}".format(
        arm, "\n       ".join(str(t) for t in tried)))
    raise AssertionError("unreachable")

def _decode(arr):
    return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]

def load_gene_based(tag, h5mu_path: Path):
    """Return (X, cell_ids, regulon_names, modalities_present). Gene_based only.

    X is read whole. It is (69,249 x 221) float64 = 116.8 MiB, so a chunked read would buy nothing
    and would make the code harder to check.
    """
    if not h5mu_path.is_file():
        die(tag, f"required input missing: {h5mu_path}\n"
                 f"       reason: Track A scores truth's frozen vocabulary; this IS that matrix")
    with h5py.File(h5mu_path, "r") as f:
        if "mod" not in f:
            die(tag, f"{h5mu_path} has no /mod group -- not a MuData file")
        mods = sorted(f["mod"].keys())
        if "Gene_based" not in mods:
            die(tag, f"{h5mu_path} has no Gene_based modality (found {mods}). Track A is defined "
                     f"on Gene_based only; see the module header.")
        g = f["mod"]["Gene_based"]
        X = np.asarray(g["X"][:], dtype=np.float64)
        cells = _decode(g["obs"]["Cell"][:])
        regs = _decode(g["var"]["_index"][:])
    return X, cells, regs, mods

def load_metadata(tag, meta_path: Path) -> pd.DataFrame:
    if not meta_path.is_file():
        die(tag, f"required input missing: {meta_path}\n"
                 f"       reason: supplies the 22 cell-type labels Track A stratifies by")
    md = pd.read_csv(meta_path)
    need = {"h5ad_barcode", "cell_type", "donor", "site", "batch"}
    miss = need - set(md.columns)
    if miss:
        die(tag, f"{meta_path} lacks column(s) {sorted(miss)}")
    return md

def load_assignment(tag, path: Path, arm: str) -> pd.DataFrame:
    """Normalise every arm table to columns (atac_barcode, assigned_rna_barcode, true_rna_barcode)."""
    df = pd.read_csv(path)
    cols = set(df.columns)
    if {"atac_barcode", "assigned_rna_barcode"} <= cols:
        out = df[["atac_barcode", "assigned_rna_barcode"]].copy()
        if "true_rna_barcode" in cols:
            out["true_rna_barcode"] = df["true_rna_barcode"]
        else:
            out["true_rna_barcode"] = df["atac_barcode"]
    elif {"atac_barcode", "rna_barcode"} <= cols:
        out = df[["atac_barcode"]].copy()
        out["assigned_rna_barcode"] = df["rna_barcode"]
        out["true_rna_barcode"] = df["atac_barcode"]
    else:
        die(tag, f"arm {arm!r}: {path} has none of the two known column layouts "
                 f"(got {sorted(cols)})")
    return out.reset_index(drop=True)

# the premise check

def strip_suffix(cells):
    """cisTopic names BMMC cells `<h5ad_barcode>___bmmc`. Strip on the last '___' only."""
    return [c.rsplit("___", 1)[0] if "___" in c else c for c in cells]

def verify_premise(tag, h5mu_path, X, cells, regs, mods, meta, arm_tables) -> dict:
    """Check the premises the scoring relies on: row counts, id joins and per-arm bijection."""
    log = make_log(tag)
    block = {
        "h5mu": str(h5mu_path),
        "h5mu_sha256": sha256(h5mu_path),
        "P1_gene_based_modality_present": True,          # load_gene_based already died otherwise
        "modalities_in_h5mu": mods,
        "region_based_present": "Region_based" in mods,
        "region_based_note": ("present but deliberately unused: Region_based depends on the ATAC "
                              "side, which is byte-identical across arms, and the fixed "
                              "analysis plan uses Gene_based throughout. Not the object of "
                              "this analysis."),
        "n_rows_gene_based": int(X.shape[0]),
        "n_rows_expected": EXPECTED_N_CELLS,
        "n_regulons_gene_based": int(X.shape[1]),
    }
    if X.shape[0] != EXPECTED_N_CELLS or len(cells) != EXPECTED_N_CELLS:
        die(tag, f"row-count check FAILED: Gene_based has {X.shape[0]} X rows / {len(cells)} obs ids, "
                 f"expected {EXPECTED_N_CELLS}")
    block["P2_rows_equal_expected_cell_count"] = True
    if X.shape[1] != len(regs):
        die(tag, f"row-count check FAILED: X has {X.shape[1]} columns but var/_index has {len(regs)}")

    bare = strip_suffix(cells)
    suffixes = sorted({c.rsplit("___", 1)[1] for c in cells if "___" in c})
    block["cell_id_example"] = cells[0]
    block["cell_id_suffixes_observed"] = suffixes
    block["cell_ids_without_suffix"] = int(sum("___" not in c for c in cells))
    if len(set(bare)) != len(bare):
        die(tag, f"unique-id check FAILED: stripped row ids are not unique "
                 f"({len(bare) - len(set(bare))} duplicates)")
    block["P3_stripped_row_ids_unique"] = True

    meta_set = set(meta["h5ad_barcode"])
    if len(meta_set) != len(meta):
        die(tag, f"metadata-join check FAILED: b00_cell_metadata.csv h5ad_barcode is not unique")
    if set(bare) != meta_set:
        only_h5mu = sorted(set(bare) - meta_set)[:5]
        only_meta = sorted(meta_set - set(bare))[:5]
        die(tag, f"metadata-join check FAILED: stripped h5mu row ids != metadata h5ad_barcode set. "
                 f"only-in-h5mu e.g. {only_h5mu}; only-in-metadata e.g. {only_meta}. "
                 f"A convention that drifts by one separator matches zero cells.")
    block["P4_row_ids_equal_metadata_barcode_set"] = True
    log(f"premise checks 1-4 PASS: Gene_based {X.shape[0]:,} x {X.shape[1]}, row ids join to "
        f"b00_cell_metadata.csv exactly (suffix {suffixes})")

    row_of = {b: i for i, b in enumerate(bare)}
    per_arm = {}
    for arm, path in arm_tables.items():
        df = load_assignment(tag, path, arm)
        n = len(df)
        assigned = df["assigned_rna_barcode"]
        missing = int((~assigned.isin(row_of)).sum())
        n_missing_unique = int(len(set(assigned) - set(row_of)))
        
        if missing != 0:
            die(tag, f"assignment-coverage check FAILED for arm {arm!r}: {missing:,} of {n:,} assigned RNA barcodes "
                     f"({n_missing_unique:,} distinct) have no row in truth's AUCell matrix. "
                     f"PBMC's p_08 absorbed this as a declared deviation for its `link` arm; the "
                     f"BMMC design says it cannot happen, so a non-zero count is a FINDING and "
                     f"this script refuses to proceed.")
        true_eq_atac = bool((df["true_rna_barcode"] == df["atac_barcode"]).all())
        dom_unique = bool(df["atac_barcode"].is_unique)
        img_unique = bool(assigned.is_unique)
        img_eq_dom = bool(set(assigned) == set(df["atac_barcode"]))
        is_ladder = arm in LADDER_ARMS
        
        if is_ladder and not true_eq_atac:
            die(tag, f"truth-identity check FAILED for arm {arm!r}: true_rna_barcode != atac_barcode. BMMC is a "
                     f"true multiome; Track A derives every truth row index from that identity.")

        if is_ladder and not (dom_unique and img_unique and img_eq_dom):
            die(tag, f"bijection check FAILED for arm {arm!r}: not an exact bijection "
                     f"(domain unique {dom_unique}, image unique {img_unique}, "
                     f"image==domain {img_eq_dom}). Track A's gather is only a row permutation "
                     f"of truth's matrix if this holds.")
        per_arm[arm] = {
            "table": str(path),
            "table_sha256": sha256(path),
            "ladder_arm": is_ladder,
            "n_rows": n,
            "n_cells_dropped_unscoreable": missing,
            "n_distinct_rna_cells_used": int(assigned.nunique()),
            "P6_true_rna_equals_atac": true_eq_atac,
            "P7_domain_unique": dom_unique,
            "P7_image_unique": img_unique,
            "P7_image_equals_domain": img_eq_dom,
            "is_exact_bijection": bool(dom_unique and img_unique and img_eq_dom),
        }
        log(f"premise checks 5-7 {arm:28s} n={n:,} dropped={missing} "
            f"bijection={per_arm[arm]['is_exact_bijection']} "
            f"distinct_rna={per_arm[arm]['n_distinct_rna_cells_used']:,}"
            + ("" if is_ladder else "   [non-ladder]"))
    block["per_arm"] = per_arm
    block["P5_total_cells_dropped_unscoreable"] = int(
        sum(v["n_cells_dropped_unscoreable"] for v in per_arm.values()))
    return block, row_of

def celltype_labels(tag, meta) -> dict:
    return dict(zip(meta["h5ad_barcode"], meta["cell_type"]))

def celltype_order(labels: np.ndarray, min_cells: int = MIN_CELLS_PER_TYPE):
    """Cell types with at least `min_cells` cells, sorted. Same rule p_08/p_09/p_10 use."""
    return sorted(c for c in np.unique(labels) if (labels == c).sum() >= min_cells)

def declarations() -> dict:
    """The two fixed declarations, carried into every scoring report so they travel with the
    numbers rather than living only in a docstring."""
    return {
        "decision_A_icc_construction_axis": {
            "question": "p_11's PRIMARY axis is {truth, link}; BMMC's ladder has no link arm.",
            "declared_before_any_result": True,
            "primary": ["truth", DONOR_AWARE],
            "reason": ("within the pre-specified four-arm ladder this is the best-informed "
                       "CONSTRUCTIBLE pairing -- cell-type labels and donor identity are both "
                       "known at pairing time -- so it occupies link's role as 'what you use when "
                       "you do not have truth'. rand_celltype_1to1 is the same construction with "
                       "that knowledge deliberately discarded (a degraded variant, not the "
                       "default choice); rand_global carries no cell-type signal and is excluded "
                       "from every ICC axis, as PBMC excluded it, and is shown separately as "
                       "floor_check so the exclusion is visible."),
            "also_reported": {
                "donor_contrast": ["truth", DONOR_BLIND, DONOR_AWARE],
                "primary_link_nonladder": ["truth", "link"],
                "sensitivity_ladder": ["truth", DONOR_AWARE, DONOR_BLIND, "link_1to1", "link"],
                "floor_check": ["truth", "rand_global"],
            },
            "link_arms_are_scoreable_note": (
                "link and link_1to1 were built for BMMC (infer_pairings.py, out/b04_link_report.json) "
                "but never carried into the GRN half. Track A needs only an assignment table, so "
                "they are scoreable here; both cover all 69,249 ATAC cells with zero barcodes "
                "absent from truth's matrix. Flagged ladder=false wherever they appear."),
        },
        "decision_B_ksweep_k_list": {
            "question": "p_10's k-list and k=26 bar were described as PBMC-specific.",
            "declared_before_any_result": True,
            "k_list": [1, 2, 5, 10, 26, 50],
            "bar_k": 26,
            "bar": {"STRONG": ">= 0.50", "WEAK": "0.30 - 0.50", "NOT_RESCUED": "< 0.30"},
            "verdict": "L7 IS PORTED; the PBMC k-list and bar are carried over UNCHANGED.",
            "reason": ("The pre-specified bar is at k = 26, the size of the AML reference module. "
                       "26 is an AML-project constant -- that module is 26 regulons / 14 TFs -- "
                       "not a PBMC-derived one, so it is wrong to call it PBMC-specific. "
                       "Re-deriving a different k would destroy the only comparison the bar "
                       "exists to support. Retention is checked against BMMC's own structure: "
                       "BMMC's truth vocabulary is 221 Gene_based eRegulons (measured), so every "
                       "k is drawable without replacement and the top rung is 50/221 = 22.6% of "
                       "the vocabulary against PBMC's 50/245 = 20.4%. k=1 is not a free choice "
                       "either: it is p_10's built-in consistency anchor against L6."),
        },
        "declared_deviation_shared_draws": (
            "b_08 generates each module draw-set ONCE per k and scores every arm on the identical "
            "modules, where p_10 consumed one running RNG and gave each arm different modules. "
            "The BMMC headline contrast is a one-variable arm pair; scoring the two arms on "
            "different random modules would confound the arm difference with a draw difference. "
            "No estimator changes."),
        "not_ported": (
            "The pre-specified third L6 variant -- activity residualised on cell-type "
            "one-hot + 5 diffusion components -- is not computed. PBMC's committed p_08 never "
            "computed it either, so there is no counterpart to compare against. Outstanding."),
    }

def provenance(script: str, source_script: str, inputs: dict, extra: dict | None = None) -> dict:
    p = {
        "script": script,
        "ported_from": source_script,
        "bundle": str(BUNDLE),
        "inputs": inputs,
        "library_versions": library_versions(),
        "peak_rss_bytes": peak_rss_bytes(),
    }
    if extra:
        p.update(extra)
    return p

CELLTYPE_MEAN_DEGENERACY = (
    "an exact bijection whose same-cell-type rate is 1.000 permutes truth's own RNA rows WITHIN "
    "each cell type, and a mean is invariant under permutation, so its cell-type mean profile is "
    "identical to truth's. Any statistic computed on cell-type mean profiles therefore returns its "
    "maximum by construction and has zero power to separate such arms. The value is an IDENTITY, "
    "not a result. PBMC's committed rand_subcluster_1to1 hits the same degeneracy at "
    "matched_median = 1.000000 and was not flagged there.")

def within_celltype_bijection(df, ct, is_exact_bijection):
    """(same_cell_type_rate, degenerate_at_celltype_mean_resolution) for one arm.

    Both sides of the pair are joined to the manifest's own `cell_type`; the arm table's `stratum`
    column is not used, because that column was written by the same builder that wrote the
    assignment (b_05's argument, kept)."""
    lab_true = np.array([ct[b] for b in df["true_rna_barcode"]])
    lab_got = np.array([ct[b] for b in df["assigned_rna_barcode"]])
    rate = float((lab_true == lab_got).mean())
    return rate, bool(is_exact_bijection and rate == 1.0)

PBMC_OUT = BUNDLE.parents[1] / "benchmarks/pbmc" / "out"

def read_pbmc_committed(name: str):
    """Returns the parsed object, or an explicit not-found record with the command that would get it."""
    p = PBMC_OUT / name
    if not p.is_file():
        return {"status": "unknown, no record found", "expected_path": str(p),
                "command": f"ls -l {p}"}
    return json.loads(p.read_text())

def write_json(tag, path: Path, obj):
    path.write_text(json.dumps(obj, indent=2))
    make_log(tag)(f"wrote {path}  ({path.stat().st_size:,} B)")
