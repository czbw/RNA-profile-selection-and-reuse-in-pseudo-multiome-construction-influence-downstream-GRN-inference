"""Score the external AML count matrix with both target vocabularies.

All background regulons are scored for offset centering and null comparisons.
Module and background sizes are read from the supplied definition."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ctxcore.recovery import aucs

ROOT = (Path(__file__).resolve().parents[3] / 'adult/corrected')
OUT = ROOT / "out"

MAIN = ROOT.parents[1] / "adult/published" / "GSE327386_raw_count_matrix.txt.gz"

SEED = 555
AUC_THRESHOLD = 0.05          # SCENIC+ default, as in the discovery-cohort AUCell runs
N_POOL_PER_SLOT = 200         # random-gene-set null pool depth per size slot
MIN_DETECT_FRAC = 0.10        # a gene is "expressed" if detected in >=10% of samples (NULL B pool)
rng = np.random.default_rng(SEED)

def log(m):
    print(f"[score] {m}", flush=True)

mod = json.loads((OUT / "module_definition.json").read_text())
module_regulons = mod["regulons"]
background_regulons = mod["background_regulons"]

# Expression matrix + sample table

log("loading counts ...")
X = pd.read_csv(MAIN, sep="\t", index_col=0, compression="gzip")
X = X[~X.index.duplicated(keep="first")]
genes = X.index.astype(str)
X.index = genes
log(f"matrix {X.shape[0]:,} genes x {X.shape[1]} samples")

samp = []
for s in X.columns:
    if s.upper().startswith("NBM"):
        samp.append(dict(sample=s, group="NBM", patient=None, stage=None))
    else:
        p, t = s.rsplit("_", 1)
        samp.append(dict(sample=s, group="AML", patient=p, stage=t))
st = pd.DataFrame(samp)
st.to_csv(OUT / "sample_table.csv", index=False)
log(f"samples: {(st.group=='NBM').sum()} NBM | "
    f"{(st.group=='AML').sum()} AML over {st.patient.nunique()} patients")

# Rankings — pyscenic convention (shuffle genes, rank descending, ties by first)

log("building within-sample rankings (ties broken randomly, seed 555) ...")
ex = X.T                                                    # samples x genes
ex = ex.sample(frac=1.0, replace=False, axis=1, random_state=np.random.RandomState(SEED))
gene_ranks = ex.rank(axis=1, ascending=False, method="first").astype(np.int32) - 1   # 0-based
TOTAL_GENES = gene_ranks.shape[1]
log(f"rankings: {gene_ranks.shape[0]} samples x {TOTAL_GENES:,} genes; "
    f"auc_threshold={AUC_THRESHOLD} -> top {int(AUC_THRESHOLD*TOTAL_GENES):,} ranks")

gene_set = set(gene_ranks.columns)

# rank-ascending matrix reused by singscore (higher expression -> higher rank)
ascending_gene_ranks = (TOTAL_GENES - 1) - gene_ranks
# z-scores across samples, reused by mean-z
Z = ((X.T - X.T.mean(axis=0)) / X.T.std(axis=0).replace(0, np.nan))

def score_aucell(sig_genes):
    g = [x for x in sig_genes if x in gene_set]
    if not g:
        return np.full(gene_ranks.shape[0], np.nan)
    sub = gene_ranks[g]
    w = np.ones(len(g), dtype=np.float64)
    return aucs(sub, TOTAL_GENES, w, AUC_THRESHOLD)

def score_singscore(sig_genes):
    """Foroutan et al. 2018: normalised mean rank, mapped to [0,1]."""
    g = [x for x in sig_genes if x in gene_set]
    k = len(g)
    if k == 0:
        return np.full(gene_ranks.shape[0], np.nan)
    mean_rank = ascending_gene_ranks[g].mean(axis=1).values + 1.0        # ranks 1..n
    lo = (k + 1) / 2.0
    hi = TOTAL_GENES - (k - 1) / 2.0
    return (mean_rank - lo) / (hi - lo)

def score_meanz(sig_genes):
    g = [x for x in sig_genes if x in gene_set]
    if not g:
        return np.full(gene_ranks.shape[0], np.nan)
    return Z[g].mean(axis=1).values

METHODS = {"aucell": score_aucell, "singscore": score_singscore, "meanz": score_meanz}

# Score the 431 background regulons under both vocabularies, all three methods

for vocab in ("Diagnosis", "Relapse"):
    tmap = mod["targets"][vocab]
    regs = [r for r in background_regulons if r in tmap]
    log(f"{vocab}: scoring {len(regs)} background regulons x {len(METHODS)} methods ...")
    for mname, fn in METHODS.items():
        M = np.column_stack([fn(tmap[r]) for r in regs])
        df = pd.DataFrame(M, index=gene_ranks.index, columns=regs)
        df.to_csv(OUT / f"scores_{mname}_{vocab}.csv")
        t4_present = [r for r in module_regulons if r in df.columns]
        log(f"  {mname}: {df.shape} | T4 columns present {len(t4_present)}/{len(module_regulons)} | "
            f"module mean {df[t4_present].values.mean():.4f}")

# Random-gene-set null pool (null B) — size-matched random gene sets (primary estimator only)

detected = X.index[(X > 0).mean(axis=1) >= MIN_DETECT_FRAC].astype(str)
pool_genes = [g for g in detected if g in gene_set]
log(f"NULL B universe: {len(pool_genes):,} genes detected in >={MIN_DETECT_FRAC:.0%} of samples")

for vocab in ("Diagnosis", "Relapse"):
    tmap = mod["targets"][vocab]
    sizes = [len([x for x in tmap[r] if x in gene_set]) for r in module_regulons]
    log(f"NULL B pool {vocab}: {len(sizes)} slots, sizes {min(sizes)}-{max(sizes)}, "
        f"{N_POOL_PER_SLOT} sets/slot = {len(sizes)*N_POOL_PER_SLOT:,} random sets")
    pool = np.empty((gene_ranks.shape[0], len(sizes), N_POOL_PER_SLOT), dtype=np.float32)
    for si, k in enumerate(sizes):
        for j in range(N_POOL_PER_SLOT):
            pick = rng.choice(len(pool_genes), size=min(k, len(pool_genes)), replace=False)
            pool[:, si, j] = score_aucell([pool_genes[i] for i in pick])
        if (si + 1) % 10 == 0:
            log(f"  slot {si+1}/{len(sizes)} done")
    np.savez_compressed(OUT / f"nullB_pool_{vocab}.npz", pool=pool,
                        sizes=np.array(sizes), samples=np.array(gene_ranks.index, dtype=object))
    log(f"  saved nullB_pool_{vocab}.npz {pool.shape}")

json.dump({"seed": SEED, "auc_threshold": AUC_THRESHOLD, "total_genes": int(TOTAL_GENES),
           "n_pool_per_slot": N_POOL_PER_SLOT, "min_detect_frac": MIN_DETECT_FRAC,
           "nullB_universe_size": len(pool_genes),
           "ranking_convention": "pyscenic: shuffle gene order, rank descending method='first', "
                                 "0-based; ties among undetected genes broken by the shuffle"},
          open(OUT / "scoring_params.json", "w"), indent=1)
log("wrote scoring_params.json — scoring complete")
