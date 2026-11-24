"""Resolve module definitions, score target signatures and align patient pairs.

The optimized AUCell implementation is checked against the ctxcore kernel.
Target vocabularies and discovery-derived orientations remain arm-specific."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from ctxcore.recovery import aucs, auc2d, derive_rank_cutoff

SEED = 555
REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "adult/published" / "GSE327386_raw_count_matrix.txt.gz"
LAMBO_S3 = REPO / "data/reference/lambo_cell_markers.xlsx"
DEFS = ("Diagnosis", "Relapse")

BUNDLE = {"link": REPO / "adult/published", "repair": REPO / "adult/corrected"}

def out_dir(arm):
    return BUNDLE[arm] / "out"

# module resolution

def resolve_module(path, arm):
    """Return a flat view of a module definition for ONE arm.

    Raises rather than guessing. For an arm-keyed module the arm selects the vocabulary block;
    for a single-arm module the requested arm must be the module's own arm.
    """
    mod = json.loads(Path(path).read_text())
    if arm not in BUNDLE:
        raise SystemExit(f"[mod] unknown arm {arm!r}")
    if mod.get("arm_keyed"):
        for key in ("orientation", "targets", "background_regulons"):
            if arm not in mod[key]:
                raise SystemExit(f"[mod] arm-keyed module has no {key}[{arm!r}]")
        orient = mod["orientation"][arm]
        targets = mod["targets"][arm]
        background = mod["background_regulons"][arm]
    else:
        vocs = mod.get("vocabularies") or []
        if list(vocs) != [arm]:
            raise SystemExit(f"[mod] module {mod['name']!r} declares vocabularies={vocs}; "
                             f"it cannot be scored under arm {arm!r}")
        orient, targets = mod["orientation"], mod["targets"]
        background = mod["background_regulons"]
    for v in DEFS:
        if v not in orient or v not in targets:
            raise SystemExit(f"[mod] resolved block is missing vocabulary {v}")
    regs = list(mod["regulons"])
    for v in DEFS:
        miss_t = [r for r in regs if r not in targets[v]]
        miss_o = [r for r in regs if r not in orient[v]]
        if miss_t:
            raise SystemExit(f"[mod] {len(miss_t)} members have no {v} target set: {miss_t[:5]}")
        if miss_o:
            raise SystemExit(f"[mod] {len(miss_o)} members unsigned under {v}: {miss_o[:5]}")
        bad = [r for r in regs if r not in background]
        if bad:
            raise SystemExit(f"[mod] {len(bad)} members outside the background: {bad[:5]}")
    return {"name": mod["name"], "arm": arm, "arm_keyed": bool(mod.get("arm_keyed")),
            "file": str(Path(path)), "regulons": regs,
            "n_regulons": len(regs), "n_TFs": mod["n_TFs"], "TFs": list(mod["TFs"]),
            "background": list(background), "orientation": orient, "targets": targets,
            "frozen": mod.get("frozen", ""), "raw": mod}

def provenance(M, extra=None):
    p = {"module_name": M["name"], "module_file": M["file"],
         "arm_vocabulary_used": M["arm"],
         "arm_keyed_module": M["arm_keyed"],
         "arm_selection_note":
             ("T4_intersection is arm-keyed: orientation/targets/background were taken from the "
              f"['{M['arm']}'] block, i.e. the {M['arm']} arm's SCENIC+ vocabulary, matching the "
              "bundle this output is written into. The other arm's block was never read, and the "
              "two arms are never averaged or merged.") if M["arm_keyed"] else
             (f"single-arm module; scored under its own ({M['arm']}) vocabulary only"),
         "n_regulons": M["n_regulons"], "n_TFs": M["n_TFs"],
         "n_background_regulons": len(M["background"]),
         "module_frozen": M["frozen"]}
    if extra:
        p.update(extra)
    return p

# Expression rankings

class Expr:
    def __init__(self, log=print):
        X = pd.read_csv(MAIN, sep="\t", index_col=0, compression="gzip")
        X = X[~X.index.duplicated(keep="first")]
        X.index = X.index.astype(str)
        self.X = X
        ex = X.T.sample(frac=1.0, replace=False, axis=1,
                        random_state=np.random.RandomState(SEED))
        self.RNK = ex.rank(axis=1, ascending=False, method="first").astype(np.int32) - 1
        self.TOTAL = self.RNK.shape[1]
        self.samples = list(self.RNK.index)
        self.genes = list(self.RNK.columns)
        self.gset = set(self.genes)
        self._M = np.ascontiguousarray(self.RNK.values)
        self._col = {g: i for i, g in enumerate(self.genes)}
        self.DETECT = (ex > 0).mean(axis=0)
        log(f"[expr] {X.shape[0]:,} genes x {X.shape[1]} samples; rankings {self.RNK.shape}")

    def auc_fast(self, genes, auc_threshold):
        g = [x for x in genes if x in self.gset]
        if not g:
            return None
        idx = np.fromiter((self._col[x] for x in g), dtype=np.int64, count=len(g))
        sub = self._M[:, idx]
        w = np.ones(len(g), dtype=np.float64)
        rc = derive_rank_cutoff(auc_threshold, self.TOTAL)
        return auc2d(sub, w, rc, float((rc + 1) * w.sum()))

    def auc(self, genes, auc_threshold):
        v = self.auc_fast(genes, auc_threshold)
        return None if v is None else pd.Series(v, index=self.RNK.index)

    def selftest_auc(self, auc_threshold, n=3, size=250, log=print):
        """Assert the fast path is bit-identical to ctxcore.recovery.aucs on real draws."""
        rng = np.random.default_rng(12345)
        worst = 0.0
        for _ in range(n):
            g = list(rng.choice(self.genes, size, replace=False))
            a = aucs(self.RNK[g], self.TOTAL, np.ones(len(g)), auc_threshold)
            b = self.auc_fast(g, auc_threshold)
            worst = max(worst, float(np.abs(a - b).max()))
        if worst != 0.0:
            raise SystemExit(f"[expr] fast AUCell path differs from ctxcore.aucs by {worst}")
        log(f"[expr] fast-AUCell self-test: {n} draws of {size} genes, max|diff| vs "
            f"ctxcore.recovery.aucs = {worst:.1e} (bit-identical)")

# cohort

def load_pairs(arm):
    st = pd.read_csv(out_dir(arm) / "sample_table.csv")
    aml = st[st.group == "AML"].copy()
    aml["is_rel"] = aml.stage.str.startswith("Rel")
    pairs = {}
    for p, g in aml.groupby("patient"):
        dx, rl = sorted(g[~g.is_rel]["sample"]), sorted(g[g.is_rel]["sample"])
        if dx and rl:
            pairs[p] = (dx[0], rl[0])
    patients = sorted(pairs)
    nbm = list(st[st.group == "NBM"]["sample"])
    return st, pairs, patients, nbm

def load_scores(arm, method, vocab, background):
    """Load the arm's precomputed score matrix restricted to the background regulons."""
    S = pd.read_csv(out_dir(arm) / f"scores_{method}_{vocab}.csv", index_col=0)
    keep = [c for c in background if c in S.columns]
    if len(keep) != len(background):
        raise SystemExit(f"[scores] {method}/{vocab}: {len(background)-len(keep)} background "
                         f"regulons missing from the committed score matrix")
    return S[keep]

def module_score(S, regs, signs):
    C = S.values - S.values.mean(axis=1, keepdims=True)
    ix = [S.columns.get_loc(r) for r in regs]
    sg = np.array([signs[r] for r in regs], dtype=float)
    return pd.Series((C[:, ix] * sg).mean(axis=1), index=S.index)

def delta_of(series, pairs, patients):
    return np.array([series[pairs[p][1]] - series[pairs[p][0]] for p in patients])

def target_union(M, vocab, gset=None):
    u = set().union(*[set(M["targets"][vocab][r]) for r in M["regulons"]])
    return u if gset is None else (u, u & gset)
