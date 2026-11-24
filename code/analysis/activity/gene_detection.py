'Calculate gene detection rates in malignant RNA profiles for the retained diagnostic input tables.'
import os
import sys

import h5py
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.path.join(ROOT, "analysis", "out", "gene_detection")
PROJECT_ROOT = os.environ.get("AMLGRN_PROJECT", ROOT)
PANEL_F = os.path.join(PROJECT_ROOT, "data/rna_panel/panel_genes.txt")
META_F = os.path.join(PROJECT_ROOT, "data/rna_panel/panel_meta.csv")
H5AD = {"Diagnosis": os.path.join(PROJECT_ROOT, "paml/malignant/data/rna/malignant_diagnosis.h5ad"),
        "Relapse": os.path.join(PROJECT_ROOT, "paml/malignant/data/rna/malignant_relapse.h5ad")}
CHUNK = 20_000_000          # nonzero entries per slice
os.makedirs(WORK, exist_ok=True)

ens = open(PANEL_F).read().split()
meta = pd.read_csv(META_F)
assert list(meta["ensembl"]) == ens, "panel_meta.csv order != panel_genes.txt"
print(f"[gt02b] panel: {len(ens):,} genes | TFs flagged: {int(meta.is_tf.sum()):,}")

for stage in sys.argv[1:] or ["Diagnosis", "Relapse"]:
    path = H5AD[stage]
    with h5py.File(path, "r") as f:
        X = f["X"]
        assert X.attrs["encoding-type"] == "csr_matrix", f"unexpected X encoding in {path}"
        n_obs, n_var = (int(v) for v in X.attrs["shape"])
        var_names = np.array([x.decode() if isinstance(x, bytes) else str(x)
                              for x in f["var"][f["var"].attrs["_index"]][:]])
        nnz = X["indices"].shape[0]
        print(f"[gt02b] {stage}: {n_obs:,} cells x {n_var:,} genes | {nnz:,} stored entries")

        # per-gene nonzero cell count, streamed: bincount over column indices, weighted by
        # (data != 0) so an explicitly-stored zero is not counted as a detection.
        counts = np.zeros(n_var, dtype=np.int64)
        sums = np.zeros(n_var, dtype=np.float64)
        for lo in range(0, nnz, CHUNK):
            hi = min(lo + CHUNK, nnz)
            idx = X["indices"][lo:hi]
            dat = X["data"][lo:hi]
            nz = dat != 0
            counts += np.bincount(idx[nz], minlength=n_var)
            sums += np.bincount(idx, weights=dat.astype(np.float64), minlength=n_var)
            print(f"          .. {hi:,}/{nnz:,}", end="\r", flush=True)
    print()

    pos = pd.Series(np.arange(n_var), index=var_names)
    missing = [g for g in ens if g not in pos.index]
    assert not missing, f"{len(missing)} panel genes absent from {path} (e.g. {missing[:3]})"
    take = pos.reindex(ens).to_numpy()

    det = counts[take] / n_obs
    mean = sums[take] / n_obs
    out = pd.DataFrame({"ensembl": ens, "symbol": meta["symbol"].values,
                        "is_tf": meta["is_tf"].values,
                        "detection_rate": det, "mean_count": mean})
    dest = os.path.join(WORK, f"gene_detection_malignant_{stage}.csv")
    out.to_csv(dest, index=False)
    print(f"[gt02b] {stage}: detection median {np.median(det):.4f} "
          f"IQR [{np.percentile(det,25):.4f}, {np.percentile(det,75):.4f}] "
          f"| all-zero genes {int((det==0).sum()):,}")
    print(f"[gt02b] -> {dest}")

# Comparison against the healthy covariate, so the swap is documented rather than silent.
h = os.path.join(WORK, "gene_detection.csv")
if os.path.exists(h):
    a = pd.read_csv(h).set_index("ensembl").detection_rate
    for stage in ["Diagnosis", "Relapse"]:
        p = os.path.join(WORK, f"gene_detection_malignant_{stage}.csv")
        if os.path.exists(p):
            b = pd.read_csv(p).set_index("ensembl").detection_rate
            j = pd.concat([a, b], axis=1, keys=["healthy", "malignant"]).dropna()
            print(f"[gt02b] healthy vs malignant {stage}: Pearson r="
                  f"{j.healthy.corr(j.malignant):.4f}  Spearman="
                  f"{j.healthy.corr(j.malignant, method='spearman'):.4f}")
