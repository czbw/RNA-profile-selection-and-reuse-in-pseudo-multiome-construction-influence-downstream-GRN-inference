"""Pre-flight check for scoring healthy cells against malignant regulons: are the gene
universes compatible enough to put the scores on a common scale?"""
import os
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

HEALTHY = os.environ.get("AMLGRN_HEALTHY",
                         str(Path(__file__).resolve().parents[3] / 'out'))
MALIG = os.environ.get("AMLGRN_MALIG",
                       str(Path(__file__).resolve().parents[3] / 'paml/malignant/out'))
OUT = os.environ.get("AMLGRN_S3B_OUT",
                     str(Path(__file__).resolve().parents[3] / 'scoring/healthy/out'))
STAGES = ["Diagnosis", "Relapse"]
PAIRING = "link"                 # Option B metacells carry no patient identity
JACCARD_MIN = 0.95
COVERAGE_MIN = 0.90
EXPECTED_HEALTHY_CELLS = {"Diagnosis": 18125}     
os.makedirs(OUT, exist_ok=True)

def log(m):
    print(f"[s3b-preflight] {m}", flush=True)

def idx(h5, mod, axis):
    g = h5["mod"][mod][axis]
    key = g.attrs["_index"]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[key][:]])

def id_space(names):
    """'ensembl' | 'symbol' | 'mixed' — which identifier space a gene index is in.

    The h5ads are keyed by Ensembl ID with a separate `symbol` column, while eRegulon `Gene`
    entries are symbols. Comparing the two spaces gives an overlap of exactly zero, which would
    look like a coverage failure but only means the wrong columns were compared. The space is
    therefore detected and reported explicitly rather than inferred from a zero overlap.
    """
    n = list(names[:2000])
    if not n:
        return "empty"
    ens = sum(1 for x in n if str(x).startswith("ENSG"))
    if ens > 0.9 * len(n):
        return "ensembl"
    if ens < 0.1 * len(n):
        return "symbol"
    return "mixed"

def assert_same_space(spaces, what):
    """spaces: {label: space}. Abort with a distinct exit code if they disagree."""
    distinct = set(spaces.values())
    for lab, sp in spaces.items():
        log(f"ID-SPACE  {lab}: {sp}")
    if len(distinct) > 1:
        log("")
        log("FATAL: identifier-space mismatch — " + what)
        log("       This is NOT a data incompatibility. The gene indices are in different")
        log("       namespaces (Ensembl IDs vs gene symbols), so every overlap statistic below")
        log("       would read as zero and the gate would 'fail' for the wrong reason.")
        log("       Fix: remap to one namespace (the h5ads carry a `symbol` column) and")
        log("       re-run. Do not relax the bars.")
        sys.exit(5)
    return distinct.pop()

def read_mudata(path, tag):
    if not os.path.exists(path):
        log(f"FATAL: missing {path}")
        sys.exit(3)
    with h5py.File(path, "r") as f:
        mods = list(f["mod"].keys())
        rna = next((m for m in mods if "rna" in m.lower() or "gex" in m.lower()), None)
        atac = next((m for m in mods if "atac" in m.lower() or "acc" in m.lower()), None)
        if rna is None or atac is None:
            log(f"FATAL: could not identify RNA/ATAC modalities in {path} (found {mods})")
            sys.exit(3)
        genes = set(idx(f, rna, "var"))
        regions = set(idx(f, atac, "var"))
        ncell = len(idx(f, rna, "obs"))
    log(f"{tag}: {ncell:,} cells | {len(genes):,} genes | {len(regions):,} regions")
    return genes, regions, ncell

report, ok = [], True
h_gene, h_region, h_cells = {}, {}, {}
m_gene = {}

# ---- load ---------------------------------------------------------------------------
for s in STAGES:
    h_gene[s], h_region[s], h_cells[s] = read_mudata(
        os.path.join(HEALTHY, f"mudata_{s}_{PAIRING}.h5mu"), f"HEALTHY {s}")
for s in STAGES:
    g, _, n = read_mudata(os.path.join(MALIG, f"mudata_{s}_{PAIRING}.h5mu"), f"MALIGNANT {s}")
    m_gene[s] = g

for s, exp in EXPECTED_HEALTHY_CELLS.items():
    got = h_cells.get(s)
    if got != exp:
        log(f"FATAL: healthy {s} mudata holds {got:,} cells, expected {exp:,}. This is not the "
            f"file Stage 1 produced — refusing to proceed.")
        sys.exit(3)
    log(f"PROVENANCE  healthy {s} cell count {got:,} matches the Stage-1 record. OK")

# ---- CHECK 0: identifier space — must precede every overlap statistic -----------------
log("")
gex_space = assert_same_space(
    {**{f"healthy {s} GEX": id_space(sorted(h_gene[s])) for s in STAGES},
     **{f"malignant {s} GEX": id_space(sorted(m_gene[s])) for s in STAGES}},
    "the four GEX matrices are not all in the same gene namespace")
report.append(dict(check="gex_identifier_space", value=gex_space, bar=None, passed=True))

# ---- CHECK 1: healthy vs malignant gene universes ------------------------------------
for hs in STAGES:
    for ms in STAGES:
        a, b = h_gene[hs], m_gene[ms]
        jac = len(a & b) / len(a | b)
        c = jac >= JACCARD_MIN
        ok &= c
        log(f"CHECK 1  healthy {hs} GEX vs malignant {ms} GEX Jaccard = {jac:.4f} "
            f"(bar {JACCARD_MIN})  {'PASS' if c else 'FAIL'}")
        report.append(dict(check=f"gene_universe_healthy{hs}_vs_malignant{ms}",
                           value=jac, bar=JACCARD_MIN, passed=bool(c)))

# ---- CHECK 2: the two healthy stages against each other ------------------------------
a, b = h_gene["Diagnosis"], h_gene["Relapse"]
jac = len(a & b) / len(a | b)
c = jac >= JACCARD_MIN
ok &= c
log(f"CHECK 2  healthy Diagnosis vs healthy Relapse GEX Jaccard = {jac:.4f} "
    f"(bar {JACCARD_MIN})  {'PASS' if c else 'FAIL'}")
log("         ^ without this the healthy Dx->Rel contrast is not on a common scale either.")
report.append(dict(check="gene_universe_healthyDx_vs_healthyRel",
                   value=jac, bar=JACCARD_MIN, passed=bool(c)))

# ---- CHECK 3: malignant regulon target coverage in the healthy GEX --------------------
targets, regions_of = {}, {}
for ms in STAGES:
    d = pd.read_csv(os.path.join(MALIG, f"eRegulons_{ms}_{PAIRING}.tsv"), sep="\t",
                    usecols=["eRegulon_name", "Gene", "Region"])
    d["regulon"] = d.eRegulon_name.str.replace(r"_\(\d+[gr]\)$", "", regex=True)
    targets[ms] = {r: set(g.Gene.astype(str)) for r, g in d.groupby("regulon")}
    regions_of[ms] = set(d.Region.astype(str))
    log(f"malignant {ms}: {len(targets[ms]):,} regulons, {d.Gene.nunique():,} target genes")
    
    tsp = id_space(sorted(set(d.Gene.astype(str))))
    if tsp != gex_space:
        log("")
        log(f"FATAL: eRegulons_{ms}_link.tsv Gene column is '{tsp}' but the GEX matrices are "
            f"'{gex_space}'.")
        log("       Coverage would compute as zero for every regulon and the gate would 'fail'")
        log("       for the wrong reason. Remap to one namespace and re-run.")
        sys.exit(5)

for ms in STAGES:
    for hs in STAGES:
        cov = np.array([len(t & h_gene[hs]) / len(t) for t in targets[ms].values() if t])
        med = float(np.median(cov))
        c = med >= COVERAGE_MIN
        ok &= c
        log(f"CHECK 3  malignant {ms} regulons scored on healthy {hs} cells: median coverage "
            f"= {med:.4f} (bar {COVERAGE_MIN})  {'PASS' if c else 'FAIL'}   "
            f"min {cov.min():.3f} | frac<0.5 {(cov < 0.5).mean():.3f}")
        report.append(dict(check=f"target_coverage_malignant{ms}_on_healthy{hs}",
                           value=med, bar=COVERAGE_MIN, passed=bool(c)))

# ---- CHECK 4: region overlap (informational, expected to be low) ----------------------------
for ms in STAGES:
    for hs in STAGES:
        hit = len(regions_of[ms] & h_region[hs]) / max(len(regions_of[ms]), 1)
        log(f"CHECK 4  malignant {ms} regulon regions present in healthy {hs} ACC: {hit:.4f} "
            f"(INFORMATIONAL)")
        report.append(dict(check=f"region_hit_malignant{ms}_on_healthy{hs}",
                           value=hit, bar=None, passed=None))
log("         ^ THIS IS WHY STAGE 3B USES Gene_based ONLY. Region_based in a cross-scored file")
log("           is meaningless by construction — do not read it.")

# ---- memory guidance -----------------------------------------------------------------
log("")
log("MEMORY GUIDANCE (AUCell peak ~ dense cells x regions, twice — matrix + rank matrix):")
rec = 0
for hs in STAGES:
    gb = h_cells[hs] * len(h_region[hs]) * 8 / 1e9
    need = gb * 2 * 1.6                      # matrix + rank matrix + 60% headroom
    rec = max(rec, need)
    log(f"  healthy {hs}: {h_cells[hs]:,} cells x {len(h_region[hs]):,} regions "
        f"-> ~{gb:.0f} GB dense, peak ~{need:.0f} GB")
    report.append(dict(check=f"estimated_peak_GB_healthy{hs}", value=need, bar=None, passed=None))
log(f"  -> recommend --mem={int(np.ceil(rec / 50.0) * 50)}G  "
    f"(bigmem MaxMemPerCPU=64430M, so --cpus-per-task >= {int(np.ceil(rec * 1024 / 64430))})")

pd.DataFrame(report).to_csv(os.path.join(OUT, "s3b_preflight_report.csv"), index=False)

verdict = ("PASS — healthy and malignant gene universes are compatible; submit "
           "submit_crossscore.sh"
           if ok else
           "FAIL — gene universes are NOT compatible. Scoring malignant regulons on healthy\n"
           "       cells would not be on a common scale, and the bystander comparison would be\n"
           "       invalid. Do NOT submit the cross-scoring. Read s3b_preflight_report.csv.")
log("")
log("=" * 78)
log(f"VERDICT: {verdict}")
log("=" * 78)
with open(os.path.join(OUT, "s3b_preflight_verdict.txt"), "w") as fh:
    fh.write(verdict + "\n")
sys.exit(0 if ok else 3)
