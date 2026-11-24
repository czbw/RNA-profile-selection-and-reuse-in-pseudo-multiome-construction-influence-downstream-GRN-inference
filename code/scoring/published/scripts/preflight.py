"""Pre-flight check for the cross-stage re-scoring: are the Diagnosis and Relapse gene universes
compatible enough to score each stage's cells with the other stage's regulons?"""
import os
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

S2 = os.environ.get("AMLGRN_S2", str(Path(__file__).resolve().parents[3] / 'paml/malignant/out'))
OUT = os.environ.get("AMLGRN_S3_OUT", str(Path(__file__).resolve().parents[3] / 'scoring/published/out'))
STAGES = ["Diagnosis", "Relapse"]
PAIRING = "link"                       # Option B carries no patient identity, so Stage 3 uses Option A only
JACCARD_MIN = 0.95
COVERAGE_MIN = 0.90
os.makedirs(OUT, exist_ok=True)

def log(m):
    print(f"[s3-preflight] {m}", flush=True)

def idx(h5, mod, axis):
    g = h5["mod"][mod][axis]
    key = g.attrs["_index"]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[key][:]])

report, ok = [], True
uni_gene, uni_region, n_cells = {}, {}, {}

for s in STAGES:
    p = os.path.join(S2, f"mudata_{s}_{PAIRING}.h5mu")
    if not os.path.exists(p):
        log(f"FATAL: missing {p}")
        sys.exit(3)
    with h5py.File(p, "r") as f:
        mods = list(f["mod"].keys())
        rna = next((m for m in mods if "rna" in m.lower() or "gex" in m.lower()), None)
        atac = next((m for m in mods if "atac" in m.lower() or "acc" in m.lower()), None)
        if rna is None or atac is None:
            log(f"FATAL: could not identify RNA/ATAC modalities in {p} (found {mods})")
            sys.exit(3)
        uni_gene[s] = set(idx(f, rna, "var"))
        uni_region[s] = set(idx(f, atac, "var"))
        n_cells[s] = len(idx(f, rna, "obs"))
    log(f"{s}: {n_cells[s]:,} cells | {len(uni_gene[s]):,} genes | {len(uni_region[s]):,} regions")

# ---- CHECK 1: GEX gene-universe agreement -------------------------------------------
gA, gB = uni_gene["Diagnosis"], uni_gene["Relapse"]
jac = len(gA & gB) / len(gA | gB)
c1 = jac >= JACCARD_MIN
ok &= c1
log(f"CHECK 1  GEX gene-universe Jaccard = {jac:.4f}  (bar {JACCARD_MIN})  "
    f"{'PASS' if c1 else 'FAIL'}   shared {len(gA & gB):,} | Dx-only {len(gA - gB):,} | "
    f"Rel-only {len(gB - gA):,}")
report.append(dict(check="gex_gene_universe_jaccard", value=jac, bar=JACCARD_MIN, passed=bool(c1)))

# ---- CHECK 2: per-regulon target coverage in the OTHER stage's GEX -------------------
targets = {}
for s in STAGES:
    d = pd.read_csv(os.path.join(S2, f"eRegulons_{s}_{PAIRING}.tsv"), sep="\t",
                    usecols=["eRegulon_name", "Gene", "Region"])
    d["regulon"] = d.eRegulon_name.str.replace(r"_\(\d+[gr]\)$", "", regex=True)
    targets[s] = {r: set(g.Gene.astype(str)) for r, g in d.groupby("regulon")}
    globals().setdefault("regions_of", {})[s] = set(d.Region.astype(str))
    log(f"{s}: {len(targets[s]):,} regulons, {d.Gene.nunique():,} distinct target genes")

for src, dst in [("Relapse", "Diagnosis"), ("Diagnosis", "Relapse")]:
    cov = np.array([len(t & uni_gene[dst]) / len(t) for t in targets[src].values() if t])
    med = float(np.median(cov))
    c2 = med >= COVERAGE_MIN
    ok &= c2
    log(f"CHECK 2  {src}-defined regulons scored on {dst} cells: median target coverage "
        f"= {med:.4f} (bar {COVERAGE_MIN})  {'PASS' if c2 else 'FAIL'}   "
        f"min {cov.min():.3f} | frac<0.5 {(cov < 0.5).mean():.3f}")
    report.append(dict(check=f"target_coverage_{src}_defs_on_{dst}_cells",
                       value=med, bar=COVERAGE_MIN, passed=bool(c2)))

# ---- CHECK 3: region overlap (informational, expected to be low) ---------------------------
rA, rB = uni_region["Diagnosis"], uni_region["Relapse"]
rjac = len(rA & rB) / len(rA | rB)
log(f"CHECK 3  ACC region-universe Jaccard = {rjac:.4f}  (INFORMATIONAL — no bar)")
for src, dst in [("Relapse", "Diagnosis"), ("Diagnosis", "Relapse")]:
    hit = len(regions_of[src] & uni_region[dst]) / max(len(regions_of[src]), 1)
    log(f"         {src}-defined regulon regions present in {dst} ACC: {hit:.4f}")
    report.append(dict(check=f"region_hit_{src}_defs_on_{dst}_cells", value=hit,
                       bar=None, passed=None))
report.append(dict(check="acc_region_universe_jaccard", value=rjac, bar=None, passed=None))
log("         ^ THIS IS WHY STAGE 3 USES Gene_based ONLY. Region_based in a cross-scored")
log("           file is meaningless by construction — do not read it.")

pd.DataFrame(report).to_csv(os.path.join(OUT, "s3_preflight_report.csv"), index=False)

verdict = ("PASS — gene universes are compatible; submit submit_crossscore.sh"
           if ok else
           "FAIL — gene universes are NOT compatible; do NOT spend bigmem. The cross-stage\n"
           "       contrast would not be on a common scale. Reconsider: rebuild both GEX\n"
           "       matrices on a shared gene set, or fall back to pooled joint inference.")
log("")
log("=" * 76)
log(f"VERDICT: {verdict}")
log("=" * 76)
with open(os.path.join(OUT, "s3_preflight_verdict.txt"), "w") as fh:
    fh.write(verdict + "\n")
sys.exit(0 if ok else 3)
