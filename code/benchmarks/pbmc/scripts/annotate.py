"""Cluster the RNA and assign marker-based cell-type labels (needed by the `pseudomultiome` arm)."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

MARKERS = {
    "CD4.T":       ["IL7R", "CCR7", "CD3D", "CD3E", "LTB"],
    "CD8.T":       ["CD8A", "CD8B", "GZMK", "CD3D"],
    "NK":          ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B":           ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Monocyte.CD14": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
    "Monocyte.CD16": ["FCGR3A", "MS4A7", "CDKN1C"],
    "DC":          ["FCER1A", "CST3", "CLEC10A"],
    "pDC":         ["LILRA4", "IRF7", "GZMB"],
    "Granulocyte": ["FCGR3B", "CSF3R", "CXCR2", "ELANE", "MPO"],
    "Platelet":    ["PPBP", "PF4"],
    "HSPC":        ["CD34", "KIT", "SPINK2"],
}

def log(m):
    print(f"[p01] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument("--config", required=True)
ap.add_argument("--dataset", required=True)
ap.add_argument("--resolution", type=float, default=1.0)
a = ap.parse_args()

ROOT = Path(os.environ.get("AMLGRN_ROOT", (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')))
cfg = yaml.safe_load(open(a.config, encoding="utf-8"))
WORK = ROOT / cfg["paths"]["work"] / a.dataset
seed = cfg.get("seed", 555)

ad = sc.read_h5ad(WORK / "gex.h5ad")
log(f"loaded {ad.n_obs:,} cells x {ad.n_vars:,} genes")

ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000)
ad.raw = ad
adh = ad[:, ad.var.highly_variable].copy()
sc.pp.scale(adh, max_value=10)
sc.tl.pca(adh, svd_solver="arpack", random_state=seed)
sc.pp.neighbors(adh, n_neighbors=15, random_state=seed)

# clustering: Leiden (leidenalg, else igraph flavour), then Louvain, then KMeans as a last resort
backend = None
for flavor, need in (("leidenalg", "leidenalg"), ("igraph", "igraph")):
    try:
        __import__(need)
    except ImportError:
        continue
    kw = dict(resolution=a.resolution, random_state=seed, key_added="leiden")
    if flavor == "igraph":
        kw.update(flavor="igraph", n_iterations=2, directed=False)
    sc.tl.leiden(adh, **kw)
    backend = f"leiden[{flavor}]"
    break

if backend is None:
    try:
        import louvain  # noqa: F401
        sc.tl.louvain(adh, resolution=a.resolution, random_state=seed, key_added="leiden")
        backend = "louvain"
    except ImportError:
        from sklearn.cluster import KMeans
        k = int(cfg.get("fallback_k", 12))
        log(f"WARNING: no leidenalg/igraph/louvain in this environment -- falling back to "
            f"KMeans(k={k}) on the PCA. This is a DIFFERENT clustering algorithm; the label set "
            f"is still reproducible and declared, but install leidenalg on CUBIC for the real run.")
        km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(adh.obsm["X_pca"][:, :30])
        adh.obs["leiden"] = [str(x) for x in km.labels_]
        backend = f"kmeans[k={k}]"

ad.obs["leiden"] = adh.obs["leiden"].astype(str).values
n_cl = ad.obs.leiden.nunique()
log(f"clustering backend = {backend} @ resolution {a.resolution}: {n_cl} clusters "
    f"(sizes {sorted(ad.obs.leiden.value_counts().tolist(), reverse=True)[:8]} ...)")

# marker scoring -- descriptive naming only

scores = {}
for name, genes in MARKERS.items():
    present = [g for g in genes if g in ad.raw.var_names]
    if not present:
        continue
    sc.tl.score_genes(ad, present, score_name=f"_s_{name}", random_state=seed, use_raw=True)
    scores[name] = f"_s_{name}"
    if len(present) < len(genes):
        log(f"  {name}: {len(present)}/{len(genes)} markers present")

S = ad.obs[list(scores.values())]
S.columns = list(scores.keys())
per_cluster = S.groupby(ad.obs.leiden.values).mean()
best = per_cluster.idxmax(axis=1)
margin = per_cluster.max(axis=1) - per_cluster.apply(lambda r: r.nlargest(2).iloc[-1], axis=1)

# a cluster with no clear winner keeps its numeric id rather than being forced into a label
MARGIN_MIN = 0.05
label = {}
for cl in per_cluster.index:
    label[cl] = f"{best[cl]}" if margin[cl] >= MARGIN_MIN else f"cluster{cl}"
counts = pd.Series(label).value_counts()
if (counts > 1).any():                       # disambiguate duplicate names
    seen = {}
    for cl in per_cluster.index:
        b = label[cl]
        if counts.get(b, 0) > 1:
            seen[b] = seen.get(b, 0) + 1
            label[cl] = f"{b}.{seen[b]}"
ad.obs["annotated_celltype"] = ad.obs.leiden.map(label).astype(str)

vc = ad.obs.annotated_celltype.value_counts()
log("annotation:")
for k, v in vc.items():
    log(f"  {k:18s} {v:6,d}")

mn = cfg["pairing_params"]["pseudomultiome"]["min_cells_per_celltype"]
thin = vc[vc < mn]
if len(thin):
    log(f"WARNING: {len(thin)} label(s) below the {mn}-cell floor and will be dropped before "
        f"prepare_GEX_ACC: {dict(thin)}  (a 3-4 cell type crashes SCENIC+ metacell sampling)")

ad.write_h5ad(WORK / "gex_annotated.h5ad")
cd = pd.read_csv(WORK / "cell_data.csv", index_col=0)
cd["annotated_celltype"] = ad.obs["annotated_celltype"].reindex(cd.index).values
cd["leiden"] = ad.obs["leiden"].reindex(cd.index).values
cd.to_csv(WORK / "cell_data.csv")

json.dump({"resolution": a.resolution, "seed": seed, "n_clusters": int(n_cl),
           "clustering_backend": backend,
           "labels": {k: int(v) for k, v in vc.items()},
           "below_floor": {k: int(v) for k, v in thin.items()},
           "margin_min": MARGIN_MIN,
           "note": "labels feed ONLY the pseudomultiome arm; truth and link do not consult them"},
          open(WORK / "annotate_report.json", "w"), indent=1)
log("wrote gex_annotated.h5ad + annotate_report.json -- STEP 01 COMPLETE")
