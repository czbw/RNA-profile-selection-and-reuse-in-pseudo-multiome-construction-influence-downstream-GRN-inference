"""Feasibility checks for replicating the module in the GSE327386 cohort (design, identifiers, target coverage, depth, controls)."""
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'adult/published')
REPO = ROOT.parents[1]
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)
S2 = REPO / "paml/malignant" / "out"
S3 = REPO / "longitudinal/published" / "out"

MAIN = ROOT / "GSE327386_raw_count_matrix.txt.gz"
SERIAL = ROOT / "GSE327386_raw_count_matrix_serial_samples.txt.gz"
COVERAGE_BAR = 0.90          # F3: minimum median target coverage
MIN_PAIRED = 5               # F1: minimum paired patients for the paired estimator

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))

def log(m):
    print(f"[feas] {m}", flush=True)

res = {"substrate": "GSE327386", "bars": {"coverage": COVERAGE_BAR, "min_paired": MIN_PAIRED}}

# F1 — design

with gzip.open(MAIN, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
samples = header[1:]
nbm = [s for s in samples if s.upper().startswith("NBM")]
aml = [s for s in samples if not s.upper().startswith("NBM")]

def split_sample(s):
    m = re.match(r"^(AML\d+)_(Dx|Rel\d*)$", s)
    return (m.group(1), m.group(2)) if m else (None, None)

parsed = [split_sample(s) for s in aml]
unparsed = [s for s, (p, t) in zip(aml, parsed) if p is None]
pat = {}
for s, (p, t) in zip(aml, parsed):
    if p:
        pat.setdefault(p, {})[t] = s

paired = sorted([p for p, d in pat.items() if "Dx" in d and any(k.startswith("Rel") for k in d)])
dx_only = sorted([p for p, d in pat.items() if "Dx" in d and not any(k.startswith("Rel") for k in d)])
rel_only = sorted([p for p, d in pat.items() if "Dx" not in d])

res["F1_design"] = {
    "n_samples_total": len(samples), "n_NBM": len(nbm), "n_AML_samples": len(aml),
    "n_patients": len(pat), "n_paired_Dx_Rel": len(paired),
    "n_Dx_only": len(dx_only), "n_Rel_only": len(rel_only),
    "unparsed_sample_names": unparsed,
    "passes_min_paired": len(paired) >= MIN_PAIRED,
    "our_internal_cohort_n_paired": 20}
log(f"F1 samples={len(samples)} (NBM {len(nbm)}, AML {len(aml)}) | patients={len(pat)} | "
    f"PAIRED Dx+Rel = {len(paired)}  (our internal cohort: 20)")
if unparsed:
    log(f"F1 unparsed sample names ({len(unparsed)}): {unparsed[:8]}")

# Load the matrix (genes x samples)

log("loading count matrix ...")
X = pd.read_csv(MAIN, sep="\t", index_col=0, compression="gzip")
log(f"matrix: {X.shape[0]:,} genes x {X.shape[1]} samples")
genes = set(map(str, X.index))

# F2 — identifier space

ens = sum(1 for g in list(genes)[:5000] if str(g).startswith("ENSG"))
res["F2_identifier"] = {
    "n_genes": int(X.shape[0]),
    "looks_like_ensembl": ens > 100,
    "space": "ensembl" if ens > 100 else "symbol",
    "examples": list(map(str, X.index[:5]))}
log(f"F2 identifier space = {res['F2_identifier']['space']} "
    f"(examples: {', '.join(res['F2_identifier']['examples'][:4])})")

# F3 — per-regulon target coverage for the module regulons

fun = pd.read_csv(S3 / "nomination_funnel.csv", index_col=0)
T4 = list(fun.index[fun.survives_4])
T5 = list(fun.index[fun.survives_all5])

targets = {}
for stage in ("Diagnosis", "Relapse"):
    tsv = S2 / f"eRegulons_{stage}_link.tsv"
    if not tsv.exists():
        continue
    e = pd.read_csv(tsv, sep="\t", usecols=["eRegulon_name", "Gene"])
    e["key"] = [regkey(v) for v in e.eRegulon_name]
    for k, g in e.groupby("key").Gene:
        targets.setdefault(stage, {})[k] = set(map(str, g.unique()))

cov_rows = []
for stage, tmap in targets.items():
    for reg in T4:
        tg = tmap.get(reg)
        if not tg:
            continue
        hit = len(tg & genes)
        cov_rows.append(dict(vocabulary=stage, regulon=reg, TF=fun.loc[reg, "TF"],
                             n_targets=len(tg), n_present=hit, coverage=hit / len(tg),
                             in_T5=reg in T5))
cov = pd.DataFrame(cov_rows)
cov.to_csv(OUT / "t4_target_coverage.csv", index=False)

res["F3_coverage"] = {}
for stage, g in cov.groupby("vocabulary"):
    med = float(g.coverage.median())
    res["F3_coverage"][stage] = {
        "n_regulons_found": int(len(g)),
        "n_regulons_expected": len(T4),
        "median_coverage": med,
        "min_coverage": float(g.coverage.min()),
        "frac_regulons_below_bar": float((g.coverage < COVERAGE_BAR).mean()),
        "median_n_targets": float(g.n_targets.median()),
        "passes": med >= COVERAGE_BAR}
    x = res["F3_coverage"][stage]
    log(f"F3 {stage} vocab: {x['n_regulons_found']}/{len(T4)} T4 regulons found | "
        f"median target coverage {med:.4f} (bar {COVERAGE_BAR}) -> "
        f"{'PASS' if x['passes'] else 'FAIL'} | worst {x['min_coverage']:.3f} | "
        f"{x['frac_regulons_below_bar']:.1%} of regulons below bar")

# union gene space actually needed by the module
needed = set().union(*[t for tm in targets.values() for k, t in tm.items() if k in T4])
res["F3_coverage"]["module_gene_union"] = {
    "n_needed": len(needed), "n_present": len(needed & genes),
    "frac_present": len(needed & genes) / max(len(needed), 1)}
log(f"F3 module gene union: {len(needed & genes):,}/{len(needed):,} present "
    f"({len(needed & genes)/max(len(needed),1):.3f})")

# F4 — library depth / detection

tot = X.sum(axis=0)
det = (X > 0).sum(axis=0)
res["F4_depth"] = {
    "median_library_size": float(tot.median()),
    "min_library_size": float(tot.min()),
    "median_genes_detected": float(det.median()),
    "min_genes_detected": float(det.min()),
    "n_samples_below_1M_counts": int((tot < 1e6).sum())}
log(f"F4 library size median {tot.median():,.0f} (min {tot.min():,.0f}); "
    f"genes detected median {det.median():,.0f}; "
    f"{int((tot < 1e6).sum())} samples under 1M counts")

# detection of the module's own genes specifically
mg = sorted(needed & genes)
sub = X.loc[mg]
res["F4_depth"]["module_genes_zero_everywhere"] = int((sub.sum(axis=1) == 0).sum())
res["F4_depth"]["module_genes_median_detection_rate"] = float(((sub > 0).mean(axis=1)).median())
log(f"F4 module genes: {res['F4_depth']['module_genes_zero_everywhere']} of {len(mg):,} are zero in "
    f"every sample; median per-gene detection rate "
    f"{res['F4_depth']['module_genes_median_detection_rate']:.3f}")

# F5 — controls

res["F5_controls"] = {
    "n_NBM": len(nbm), "names": nbm,
    "supports_malignant_vs_healthy": len(nbm) >= 3,
    "supports_within_patient_selectivity": False,
    "note": "NBM are SEPARATE normal-bone-marrow donors, not each patient's OWN healthy "
            "compartment. They can support a malignant-vs-healthy contrast but CANNOT replicate "
            "the within-patient selectivity design, which stays unreplicated."}
log(f"F5 {len(nbm)} NBM controls — malignant-vs-healthy: yes; within-patient selectivity: NO")

# Serial-sample matrix (reported, not used)

with gzip.open(SERIAL, "rt") as fh:
    sh = fh.readline().rstrip("\n").split("\t")[1:]
res["serial_matrix"] = {"n_samples": len(sh), "samples": sh,
                        "note": "multi-relapse timepoints from 2 patients; a separate matrix with a "
                                "different gene annotation — do NOT merge with the main matrix "
                                "without re-checking the identifier space"}
log(f"serial matrix: {len(sh)} samples from a separate annotation — kept aside")

# Verdict

ok_design = res["F1_design"]["passes_min_paired"]
ok_ident = res["F2_identifier"]["space"] == "symbol"
ok_cov = all(v["passes"] for k, v in res["F3_coverage"].items() if isinstance(v, dict)
             and "passes" in v)
res["verdict"] = {
    "design": ok_design, "identifier": ok_ident, "coverage": ok_cov,
    "FEASIBLE": bool(ok_design and ok_ident and ok_cov)}
log("")
log(f"VERDICT — design {ok_design} | identifier {ok_ident} | coverage {ok_cov} => "
    f"{'FEASIBLE' if res['verdict']['FEASIBLE'] else 'NOT FEASIBLE AS-IS'}")

with open(OUT / "feasibility.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote out/feasibility.json + out/t4_target_coverage.csv")
