'Figure 1d control: size- and expression-matched random gene sets.'
import argparse
import json
import random
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from numba import njit, prange
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant"
OUT = ROOT / "analysis" / "out"

# ----------------------------------------------------------------------------- constants
STAGES = ["Diagnosis", "Relapse"]
SEED = 555                    # project-wide seed
AUC_THRESHOLD = 0.05          # SCENIC+ default (scenicplus.eregulon_enrichment.score_eRegulons)
MIN_UNITS = 20                # cell types with fewer scored units are excluded (as in f1_04)
N_EXPR_BINS = 20              # mean-expression bins used for matching
NR_CELLS_PER_METACELL = 5     # run_grn.sh --nr_cells_per_metacells
METACELL_FLOOR = 25           # config.yaml min_cells_per_celltype (both modalities)
CHUNK = 1500                  # RNA profiles ranked per block
GATE_PROFILE_MEDIAN, GATE_PROFILE_P5, GATE_UNIT_SPEARMAN = 0.99, 0.95, 0.95
N_UNIT_PARITY = 2000

def to_builtin(o):
    """json default: numpy scalars and arrays to native Python types."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serialisable: {type(o)}")

def log(msg):
    print(f"[f1-null] {msg}", flush=True)

def strip_size(name):
    """Drop the '_(<n>g)' / '_(<n>r)' size suffix so identities match across constructions."""
    return name.rsplit("_(", 1)[0]

def decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]

# ----------------------------------------------------------------------------- stored data
def read_stored_auc(stage, arm):
    """Gene-based AUCell matrix of one stored construction -> (unit names, regulons, matrix)."""
    path = S2 / "out" / f"eRegulon_AUC_{stage}_{arm}.h5mu"
    with h5py.File(path, "r") as h:
        units = decode(h["obs/Cell"][:])
        regs = [strip_size(v) for v in decode(h["mod/Gene_based/var/_index"][:])]
        mat = h["mod/Gene_based/X"][:].astype(np.float64)
    return units, regs, mat

def unit_celltypes(arm, units, cell_data):
    """Cell type of each stored unit (link: ATAC-cell label; metacell: parsed from its name)."""
    if arm == "pseudomultiome":
        return np.array([u.rsplit("_", 1)[0] for u in units])
    ct = cell_data["Classified_Celltype"].reindex(units).to_numpy()
    assert not pd.isna(ct).any(), "link units without a cell-type label"
    return ct.astype(str)

def stored_profiles(mat, regs, cts):
    """Cell-type mean AUCell score per regulon, restricted to cell types with >= MIN_UNITS units."""
    types, counts = np.unique(cts, return_counts=True)
    keep = [t for t, n in zip(types, counts) if n >= MIN_UNITS]
    prof = pd.DataFrame({t: mat[cts == t].mean(0, dtype=np.float64) for t in keep}, index=regs).T
    return prof, dict(zip(types, counts.astype(int)))

def load_gene_sets(stage, arm):
    """Target-gene sets of every eRegulon of one construction, keyed by size-stripped name."""
    tsv = pd.read_csv(S2 / "out" / f"eRegulons_{stage}_{arm}.tsv", sep="\t",
                      usecols=["Gene", "Gene_signature_name"])
    sets = {strip_size(k): sorted(set(g)) for k, g in tsv.groupby("Gene_signature_name").Gene}
    return sets

# ----------------------------------------------------------------------------- RNA objects
def rna_path(stage):
    return S2 / "data" / "rna" / f"malignant_{stage.lower()}.h5ad"

def gene_universe(path):
    """Symbol-indexed gene list exactly as build_paired_rna.py / build_native_rna.py derive it."""
    A = ad.read_h5ad(path, backed="r")
    sym = A.var["symbol"].astype(str)
    keep = np.flatnonzero((sym.notna() & (sym != "nan")).to_numpy())
    names = ad.utils.make_index_unique(pd.Index(sym.to_numpy()[keep]))
    obs = A.obs[["Library_ID", "Classified_Celltype"]].copy()
    obs["barcode"] = pd.Index(A.obs.index).str.split("|").str[-1]
    A.file.close()
    return keep, list(names), obs

def read_csr_rows(path, rows, keep_cols):
    """Rows of the X matrix of an h5ad (CSR) as a float32 CSR restricted to keep_cols."""
    rows = np.asarray(rows)
    with h5py.File(path, "r") as h:
        g = h["X"]
        indptr = g["indptr"][:]
        n_genes = int(g.attrs["shape"][1])
        data, indices = g["data"], g["indices"]
        blocks_d, blocks_i, ptr = [], [], [0]
        for r in rows:
            a, b = indptr[r], indptr[r + 1]
            blocks_d.append(data[a:b])
            blocks_i.append(indices[a:b])
            ptr.append(ptr[-1] + (b - a))
    M = sp.csr_matrix((np.concatenate(blocks_d).astype(np.float32),
                       np.concatenate(blocks_i), np.asarray(ptr)), shape=(len(rows), n_genes))
    return M[:, keep_cols]

def build_link_arm(stage, units, cts):
    """Distinct linked RNA profiles of the link construction and unit counts per cell type."""
    cell_data = pd.read_csv(S2 / "out" / f"cell_data_{stage}.csv", index_col=0)
    lut, sources = {}, {}
    for s in STAGES:                       # lookup order of build_paired_rna.py: first key wins
        keep, names, obs = gene_universe(rna_path(s))
        sources[s] = (keep, names, obs)
        for i, (lib, bc) in enumerate(zip(obs["Library_ID"].astype(str), obs["barcode"])):
            lut.setdefault(f"{lib}||{bc}", (s, i))
    keep0, names0, _ = sources[STAGES[0]]
    for s in STAGES[1:]:                   # the paired object is an outer join over one symbol set
        assert sources[s][1] == names0, "gene universes differ between stage RNA objects"
    cd = cell_data.reindex(units)
    key = cd["Linked_scRNA_Sample"].astype(str) + "||" + cd["Linked_scRNA_barcode"].astype(str)
    hits = [lut.get(k) for k in key]
    missing = sum(h is None for h in hits)
    assert missing == 0, f"{missing} stored link units have no linked RNA profile"
    uniq = sorted(set(hits))
    col = {u: j for j, u in enumerate(uniq)}
    types = sorted(set(cts))
    C = np.zeros((len(types), len(uniq)), np.float32)
    tix = {t: i for i, t in enumerate(types)}
    for h, t in zip(hits, cts):
        C[tix[t], col[h]] += 1.0
    blocks, pos, n_before = [], {}, 0
    for s in STAGES:
        rows = [i for (src, i) in uniq if src == s]
        if rows:
            blocks.append(read_csr_rows(rna_path(s), rows, sources[s][0]))
            for j, i in enumerate(rows):
                pos[(s, i)] = n_before + j
            n_before += len(rows)
    X = sp.vstack(blocks).tocsr()[[pos[u] for u in uniq]]
    log(f"{stage}/link: {len(units):,} units -> {len(uniq):,} distinct RNA profiles")
    return {"X": X.tocsr(), "genes": names0, "types": types, "C": C}

def build_metacell_arm(stage, cell_data):
    """Reproduce the SCENIC+ pseudo-multiome metacells of one stage (mean of 5 sampled cells)."""
    keep, names, obs = gene_universe(rna_path(stage))
    ct_rna = obs["Classified_Celltype"].astype(str).to_numpy()
    ct_atac = cell_data["Classified_Celltype"].astype(str).to_numpy()
    n_rna = pd.Series(ct_rna).value_counts()
    n_atac = pd.Series(ct_atac).value_counts()
    common = sorted(set(n_rna.index) & set(n_atac.index))
    kept = [c for c in common if min(n_rna[c], n_atac[c]) >= max(METACELL_FLOOR, NR_CELLS_PER_METACELL)]
    n_meta = {c: int(round(min(n_rna[c], n_atac[c]) / NR_CELLS_PER_METACELL)) * 2 for c in kept}
    sel = np.flatnonzero(np.isin(ct_rna, kept))        # GEX_anndata after annotation filtering
    sel_ct = ct_rna[sel]
    with h5py.File(rna_path(stage), "r") as h:
        g = h["X"]
        X = sp.csr_matrix((g["data"][:].astype(np.float32), g["indices"][:], g["indptr"][:]),
                          shape=tuple(g.attrs["shape"]))
    X = X[sel][:, keep].tocsr()
    groups = {c: np.flatnonzero(sel_ct == c) for c in sorted(kept)}   # Groupby: sorted keys
    names_units, members = [], []
    for c, idx in groups.items():
        for x in range(n_meta[c]):
            random.seed(x)
            members.append(random.sample(list(idx), NR_CELLS_PER_METACELL))
            names_units.append(f"{c}_{x}")
    log(f"{stage}/metacell: {len(names_units):,} metacells over {len(groups)} cell types")
    return {"X": X, "genes": names, "members": members, "units": names_units,
            "types": sorted(groups), "n_meta": n_meta}

# ----------------------------------------------------------------------------- scoring kernel
@njit(parallel=True)
def _rank_weights(X, cut, seed0, W):
    """W[c, g] = max(cut - rank, 0) for cell c; ties broken by a seeded random permutation."""
    n, g = X.shape
    for c in prange(n):
        np.random.seed(seed0 + c)
        perm = np.random.permutation(g)
        order = perm[np.argsort(-X[c][perm], kind="mergesort")]
        for j in range(g):
            W[c, j] = 0.0
        for i in range(cut):
            W[c, order[i]] = cut - i

def accumulate(blocks, n_genes, type_names, seed0, capture=()):
    """Per-cell-type mean of w_g, unit counts per type and the w rows of selected profiles.

    blocks yields (dense float32 [m, G], count matrix [T, m]) where count[t, j] is the number of
    scored units of type t that carry profile j. capture lists global profile indices whose w rows
    are returned (used for the unit-level parity gate).
    """
    cut = int(round(AUC_THRESHOLD * n_genes)) - 1
    tot = np.zeros((len(type_names), n_genes), np.float64)
    n_units = np.zeros(len(type_names))
    capture = np.asarray(sorted(capture), dtype=np.int64)
    kept, offset = [], 0
    for dense, cnt in blocks:
        W = np.empty(dense.shape, np.float32)
        _rank_weights(dense, cut, seed0 + offset, W)
        tot += cnt.astype(np.float64) @ W
        n_units += cnt.sum(1)
        inside = capture[(capture >= offset) & (capture < offset + dense.shape[0])] - offset
        if len(inside):
            kept.append(W[inside])
        offset += dense.shape[0]
    captured = np.vstack(kept) if kept else np.empty((0, n_genes), np.float32)
    return tot / np.maximum(n_units, 1)[:, None], n_units, cut, captured

def link_blocks(arm):
    X, C = arm["X"], arm["C"]
    for a in range(0, X.shape[0], CHUNK):
        yield X[a:a + CHUNK].toarray(), C[:, a:a + CHUNK]

def metacell_blocks(arm):
    X, mem, units = arm["X"], arm["members"], arm["units"]
    tix = {t: i for i, t in enumerate(arm["types"])}
    for a in range(0, len(mem), CHUNK):
        rows = mem[a:a + CHUNK]
        dense = np.vstack([np.asarray(X[r].mean(axis=0)).ravel() for r in rows]).astype(np.float32)
        cnt = np.zeros((len(tix), len(rows)), np.float32)
        for j, u in enumerate(units[a:a + CHUNK]):
            cnt[tix[u.rsplit("_", 1)[0]], j] = 1.0
        yield dense, cnt

def set_profile(wbar, idx, cut):
    """Cell-type mean AUCell score of a gene set: columns of wbar summed over idx."""
    return wbar[:, idx].sum(1) / ((cut + 1) * len(idx))

def present(gene_index, genes):
    return np.array([gene_index[g] for g in genes if g in gene_index], dtype=np.int64)

# ----------------------------------------------------------------------------- null draws
def expression_bins(mean_expr):
    """Bin id per gene: genes with zero mean form bin 0, the rest N_EXPR_BINS quantile bins."""
    bins = np.zeros(len(mean_expr), np.int64)
    nz = mean_expr > 0
    edges = np.quantile(mean_expr[nz], np.linspace(0, 1, N_EXPR_BINS + 1)[1:-1])
    bins[nz] = 1 + np.searchsorted(edges, mean_expr[nz], side="right")
    return bins

def draw_matched(rng, set_idx, bins, universe_by_bin, k):
    """k random gene sets with the same number of genes per expression bin as set_idx."""
    own = set(set_idx.tolist())
    out, short = [], 0
    for b, c in zip(*np.unique(bins[set_idx], return_counts=True)):
        pool = universe_by_bin.get(b)
        if pool is None:
            pool = np.array([], np.int64)
        pool = pool[~np.isin(pool, list(own))] if len(pool) else pool
        take = min(c, len(pool))
        short += c - take
        if take:
            keys = rng.random((k, len(pool)))
            out.append(pool[np.argpartition(keys, take - 1, axis=1)[:, :take]])
    return np.concatenate(out, axis=1), short

def pearson_rows(profile_a, profiles_b):
    """Pearson correlation of one vector with each column of profiles_b (types x draws)."""
    a = profile_a - profile_a.mean()
    b = profiles_b - profiles_b.mean(0)
    den = np.sqrt((a ** 2).sum() * (b ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a @ b) / den

def pearson_cols(pa, pb):
    """Row-wise Pearson correlation of two (draws x types) arrays."""
    a = pa - pa.mean(1, keepdims=True)
    b = pb - pb.mean(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))

def bh(p):
    """Benjamini-Hochberg adjusted p values."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    adj[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(adj, 1.0)

def emp_p(obs, null):
    null = null[np.isfinite(null)]
    return (1 + (null >= obs).sum()) / (len(null) + 1)

# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stages", nargs="+", default=STAGES)
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--max-regulons", type=int, default=0, help="0 = all shared eRegulons")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    t0 = time.time()
    suffix = f"_{args.tag}" if args.tag else ""

    meta_log = json.load(open(OUT / "aml_n_metacells.json"))
    ref_csv = pd.read_csv(OUT / "f1_celltype_agreement_allshared.csv")
    ref_csv = ref_csv[ref_csv.level == "regulon"]

    result = {"status": "DESCRIPTIVE / EXPLORATORY", "seed": SEED, "draws": args.draws,
              "auc_threshold": AUC_THRESHOLD, "n_expression_bins": N_EXPR_BINS,
              "gates": {}, "stages": {},
              "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                           "anndata": __import__("importlib.metadata").metadata.version("anndata"), "python": sys.version.split()[0]}}
    rows, draws_store, failures = [], {}, []

    for si, stage in enumerate(args.stages):
        cell_data = pd.read_csv(S2 / "out" / f"cell_data_{stage}.csv", index_col=0)
        stored, sets, prof = {}, {}, {}
        for arm in ("link", "pseudomultiome"):
            units, regs, mat = read_stored_auc(stage, arm)
            cts = unit_celltypes(arm, units, cell_data)
            prof[arm], counts = stored_profiles(mat, regs, cts)
            stored[arm] = {"units": units, "regs": regs, "mat": mat, "cts": cts}
            sets[arm] = load_gene_sets(stage, arm)
        shared_t = sorted(set(prof["link"].index) & set(prof["pseudomultiome"].index))
        shared_r = sorted(set(prof["link"].columns) & set(prof["pseudomultiome"].columns))
        if args.max_regulons:
            shared_r = shared_r[:args.max_regulons]
        log(f"{stage}: {len(shared_t)} shared cell types, {len(shared_r)} shared eRegulons")

        obs_stored = {r: np.corrcoef(prof["link"].loc[shared_t, r],
                                     prof["pseudomultiome"].loc[shared_t, r])[0, 1] for r in shared_r}
        ref = ref_csv[ref_csv.stage == stage].set_index("unit")["pearson_r"]
        g4 = max(abs(obs_stored[r] - ref[r]) for r in shared_r)
        result["gates"][f"{stage}_G4_max_abs_diff_vs_committed_table"] = float(g4)
        if g4 >= 1e-9:
            failures.append(f"{stage} G4 {g4:.3g}")

        # ---- rebuild both constructions' RNA and accumulate per-gene cell-type means
        arms = {"link": build_link_arm(stage, stored["link"]["units"], stored["link"]["cts"]),
                "pseudomultiome": build_metacell_arm(stage, cell_data)}
        got = {c: sum(1 for u in arms["pseudomultiome"]["units"] if u.rsplit("_", 1)[0] == c)
               for c in arms["pseudomultiome"]["types"]}
        logged = meta_log["stages"][stage]["per_cell_type"]
        g1 = {c: (got[c], logged.get(c)) for c in got if got[c] != logged.get(c)}
        extra = sorted(set(logged) - set(got))
        result["gates"][f"{stage}_G1_metacell_count_mismatches"] = {
            "mismatched": {k: list(v) for k, v in g1.items()}, "logged_but_not_rebuilt": extra}
        if g1 or extra:
            failures.append(f"{stage} G1 {g1} {extra}")

        wbar, ncut, nunits = {}, {}, {}
        parity_idx = np.unique(np.linspace(0, len(arms["pseudomultiome"]["units"]) - 1,
                                           N_UNIT_PARITY).astype(int))
        unit_w_meta = None
        for ai, arm in enumerate(("link", "pseudomultiome")):
            A = arms[arm]
            G = A["X"].shape[1]
            blocks = link_blocks(A) if arm == "link" else metacell_blocks(A)
            wb, nu, cut, uw = accumulate(blocks, G, A["types"], SEED + 10_000_000 * (si * 2 + ai),
                                         capture=parity_idx if arm == "pseudomultiome" else ())
            wbar[arm], ncut[arm], nunits[arm] = (wb, A["types"]), cut, nu
            if arm == "pseudomultiome":
                unit_w_meta = uw
            log(f"{stage}/{arm}: ranked {int(nu.sum()):,} units, {G:,} genes, rank cut {cut}")

        # ---- real eRegulons from the rebuilt matrices, gates G2 and G3
        gidx = {arm: {g: i for i, g in enumerate(arms[arm]["genes"])} for arm in arms}
        recon = {}
        for arm in arms:
            wb, types = wbar[arm]
            pr = {}
            for r in sets[arm]:
                idx = present(gidx[arm], sets[arm][r])
                ok = len(idx) / len(sets[arm][r]) >= 0.8
                pr[r] = set_profile(wb, idx, ncut[arm]) if ok and len(idx) else np.zeros(len(types))
            recon[arm] = pd.DataFrame(pr, index=types)
        r_prof = {}
        for arm in arms:
            keep_t = [t for t in recon[arm].index if t in prof[arm].index]
            vals = [np.corrcoef(recon[arm].loc[keep_t, r], prof[arm].loc[keep_t, r])[0, 1]
                    for r in shared_r]
            vals = np.array(vals)
            r_prof[arm] = vals
            result["gates"][f"{stage}_G2_{arm}"] = {
                "median": float(np.nanmedian(vals)), "p5": float(np.nanpercentile(vals, 5)),
                "n_regulons": int(np.isfinite(vals).sum())}
            if np.nanmedian(vals) < GATE_PROFILE_MEDIAN or np.nanpercentile(vals, 5) < GATE_PROFILE_P5:
                failures.append(f"{stage} G2 {arm} {np.nanmedian(vals):.4f}/{np.nanpercentile(vals, 5):.4f}")
        # G3: unit-level parity for the first N_UNIT_PARITY metacells
        Am = arms["pseudomultiome"]
        names_u = [Am["units"][i] for i in parity_idx]
        col = {u: i for i, u in enumerate(stored["pseudomultiome"]["units"])}
        sidx = [col[u] for u in names_u]
        smat = stored["pseudomultiome"]["mat"]
        sreg = {r: i for i, r in enumerate(stored["pseudomultiome"]["regs"])}
        sp_vals = []
        for r in shared_r[:60]:
            idx = present(gidx["pseudomultiome"], sets["pseudomultiome"][r])
            auc = unit_w_meta[:, idx].sum(1) / ((ncut["pseudomultiome"] + 1) * len(idx))
            sp_vals.append(stats.spearmanr(auc, smat[sidx, sreg[r]]).statistic)
        result["gates"][f"{stage}_G3_metacell_unit_spearman_median"] = float(np.nanmedian(sp_vals))
        if np.nanmedian(sp_vals) < GATE_UNIT_SPEARMAN:
            failures.append(f"{stage} G3 {np.nanmedian(sp_vals):.4f}")
        del unit_w_meta

        if failures:
            result["gates"]["FAILED"] = failures
            json.dump(result, open(OUT / f"f1_random_geneset_control{suffix}.json", "w"), indent=1, default=to_builtin)
            log(f"GATE FAILURE - null analysis not run: {failures}")
            sys.exit(2)
        log(f"{stage}: gates passed ({result['gates'][f'{stage}_G2_link']}, "
            f"{result['gates'][f'{stage}_G2_pseudomultiome']})")

        # ---- observed values and nulls
        tA = [t for t in wbar["link"][1] if t in shared_t]
        tB = [t for t in wbar["pseudomultiome"][1] if t in shared_t]
        assert tA == tB
        obs_recon = {r: np.corrcoef(recon["link"].loc[shared_t, r],
                                    recon["pseudomultiome"].loc[shared_t, r])[0, 1] for r in shared_r}
        # N3 (stored, full cohort): every other shared eRegulon's profile in the second construction
        PA = prof["link"].loc[shared_t, shared_r].to_numpy()
        PB = prof["pseudomultiome"].loc[shared_t, shared_r].to_numpy()
        n3 = {r: pearson_rows(PA[:, i], np.delete(PB, i, axis=1)) for i, r in enumerate(shared_r)}

        rng_base = [SEED, si]
        bin_of, uni = {}, {}
        for arm in arms:
            if arm == "link":     # mean over scored units: distinct profiles weighted by unit count
                w = arms[arm]["C"].sum(0)
                mean_expr = np.asarray(arms[arm]["X"].T @ w).ravel() / w.sum()
            else:
                mean_expr = np.asarray(arms[arm]["X"].mean(axis=0)).ravel()
            bin_of[arm] = expression_bins(mean_expr)
            allg = np.arange(len(mean_expr))
            tgt = present(gidx[arm], sorted({g for s in sets[arm].values() for g in s}))
            uni[arm] = {"N1": allg, "N2": tgt}
        nulls = {k: np.full((len(shared_r), args.draws), np.nan) for k in ("N1", "N2")}
        short = {k: 0 for k in nulls}
        for ri, r in enumerate(shared_r):
            for ni, nk in enumerate(("N1", "N2")):
                profs = {}
                for ai, arm in enumerate(("link", "pseudomultiome")):
                    idx = present(gidx[arm], sets[arm][r])
                    rng = np.random.default_rng(rng_base + [ri, ni, ai])
                    u = uni[arm][nk]
                    by_bin = {b: u[bin_of[arm][u] == b] for b in np.unique(bin_of[arm][u])}
                    D, sh = draw_matched(rng, idx, bin_of[arm], by_bin, args.draws)
                    short[nk] += sh
                    wb = wbar[arm][0]
                    sel = [wbar[arm][1].index(t) for t in shared_t]
                    profs[arm] = wb[sel][:, D].sum(2).T / ((ncut[arm] + 1) * D.shape[1])
                nulls[nk][ri] = pearson_cols(profs["link"], profs["pseudomultiome"])
            if (ri + 1) % 20 == 0:
                log(f"{stage}: null draws for {ri + 1}/{len(shared_r)} eRegulons "
                    f"({time.time() - t0:.0f} s)")
        draws_store[stage] = {"N1": nulls["N1"], "N2": nulls["N2"],
                              "N3": np.array([n3[r] for r in shared_r])}

        stage_rows = []
        for ri, r in enumerate(shared_r):
            row = {"stage": stage, "regulon": r, "TF": r.split("_")[0],
                   "n_genes_link": len(sets["link"][r]), "n_genes_metacell": len(sets["pseudomultiome"][r]),
                   "r_observed_stored": obs_stored[r], "r_observed_rebuilt": obs_recon[r]}
            for nk, vec in (("N1", nulls["N1"][ri]), ("N2", nulls["N2"][ri]), ("N3", n3[r])):
                v = vec[np.isfinite(vec)]
                row.update({f"{nk}_median": float(np.median(v)), f"{nk}_q95": float(np.quantile(v, .95)),
                            f"{nk}_p": emp_p(obs_stored[r], vec)})
            stage_rows.append(row)
        df = pd.DataFrame(stage_rows)
        for nk in ("N1", "N2", "N3"):
            df[f"{nk}_q"] = bh(df[f"{nk}_p"].to_numpy())
        rows.append(df)

        def pooled(a):
            a = a[np.isfinite(a)]
            return {"median": float(np.median(a)), "q95": float(np.quantile(a, .95)),
                    "q99": float(np.quantile(a, .99)), "max": float(a.max()), "n": int(a.size)}

        result["stages"][stage] = {
            "n_shared_celltypes": len(shared_t), "celltypes": shared_t,
            "n_shared_eregulons": len(shared_r),
            "observed_stored": {"median": float(np.median(list(obs_stored.values()))),
                                "q25": float(np.quantile(list(obs_stored.values()), .25)),
                                "q75": float(np.quantile(list(obs_stored.values()), .75))},
            "observed_rebuilt_median": float(np.median(list(obs_recon.values()))),
            "pooled_null": {nk: pooled(draws_store[stage][nk]) for nk in ("N1", "N2", "N3")},
            "n_eregulons_above_own_null_q95": {nk: int((df.r_observed_stored > df[f"{nk}_q95"]).sum())
                                               for nk in ("N1", "N2", "N3")},
            "n_eregulons_BH_q_lt_0.05": {nk: int((df[f"{nk}_q"] < 0.05).sum()) for nk in ("N1", "N2", "N3")},
            "matching_shortfall_genes": short}
        log(f"{stage}: observed median r {result['stages'][stage]['observed_stored']['median']:.3f}; "
            f"pooled null medians "
            f"{ {k: round(v['median'], 3) for k, v in result['stages'][stage]['pooled_null'].items()} }")

    result["runtime_seconds"] = round(time.time() - t0)
    pd.concat(rows).to_csv(OUT / f"f1_random_geneset_control{suffix}.csv", index=False)
    np.savez_compressed(OUT / f"f1_random_geneset_control_draws{suffix}.npz",
                        **{f"{s}_{k}": v for s, d in draws_store.items() for k, v in d.items()})
    json.dump(result, open(OUT / f"f1_random_geneset_control{suffix}.json", "w"), indent=1, default=to_builtin)
    log(f"done in {result['runtime_seconds']} s")

if __name__ == "__main__":
    main()
