#!/usr/bin/env python
"""Split the BMMC multiome h5ad into the two modality objects the ladder needs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
BLOCK = 5000            # rows per streaming block

def log(m): print(f"[b01] {m}", flush=True)

def cat(f, col):
    codes = f[f"obs/{col}"][:]
    cats = [x.decode() if isinstance(x, bytes) else str(x)
            for x in f[f"obs/__categories/{col}"][:]]
    return np.array([cats[c] for c in codes], dtype=object)

def digest(arr) -> str:
    h = hashlib.sha256()
    for x in arr:
        h.update(str(x).encode())
    return h.hexdigest()[:16]

def write_h5ad(path: Path, M: sp.csr_matrix, obs_names, var_names, obs_cols: dict):
    """Minimal AnnData-compatible writer: X as CSR, obs/var as string-index dataframes."""
    import h5py as _h5
    with _h5.File(path, "w") as o:
        o.attrs["encoding-type"] = "anndata"; o.attrs["encoding-version"] = "0.1.0"
        g = o.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"; g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array(M.shape, dtype="int64")
        g.create_dataset("data", data=M.data, compression="gzip", compression_opts=1)
        g.create_dataset("indices", data=M.indices, compression="gzip", compression_opts=1)
        g.create_dataset("indptr", data=M.indptr, compression="gzip", compression_opts=1)
        st = _h5.string_dtype("utf-8")
        for name, idx, cols in (("obs", obs_names, obs_cols), ("var", var_names, {})):
            d = o.create_group(name)
            d.attrs["encoding-type"] = "dataframe"; d.attrs["encoding-version"] = "0.1.0"
            d.attrs["_index"] = "_index"
            ## must go through attrs.create with an explicit variable-length string type:
            ## a plain object-dtype array has no native HDF5 equivalent, and an EMPTY
            ## column list (the var frame) fails the same way
            d.attrs.create("column-order", np.array(list(cols), dtype=object), dtype=st)
            d.create_dataset("_index", data=np.array(idx, dtype=object), dtype=st)
            for c, v in cols.items():
                d.create_dataset(c, data=np.array(v, dtype=object), dtype=st)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default=str(ROOT / "benchmarks/bmmc" / "out"))
    a = ap.parse_args()
    OUT = Path(a.outdir); OUT.mkdir(parents=True, exist_ok=True)

    f = h5py.File(a.h5ad, "r")
    obs_names = np.array([x.decode() for x in f["obs/_index"][:]], dtype=object)
    var_names = np.array([x.decode() for x in f["var/_index"][:]], dtype=object)
    n_cells, n_feat = len(obs_names), len(var_names)

    ftc = f["var/feature_types"][:]
    ftcat = [x.decode() for x in f["var/__categories/feature_types"][:]]
    ft = np.array([ftcat[c] for c in ftc])
    gex_i = np.where(ft == "GEX")[0]
    atac_i = np.where(ft == "ATAC")[0]

    # --- the contiguity assumption is CHECKED, not trusted -------------------------
    if not (gex_i.max() < atac_i.min()
            and np.array_equal(gex_i, np.arange(gex_i.min(), gex_i.max() + 1))
            and np.array_equal(atac_i, np.arange(atac_i.min(), atac_i.max() + 1))):
        raise SystemExit("GEX/ATAC blocks are NOT contiguous — the column slice below "
                         "would scramble features. Refusing to write.")
    SPLIT = int(atac_i.min())
    log(f"{n_cells:,} cells x {n_feat:,} features | GEX 0..{SPLIT-1} | ATAC {SPLIT}..{n_feat-1}")

    # --- peak names to pycisTopic form ---------------------------------------------
    atac_names = np.array([v.replace("-", ":", 1) for v in var_names[atac_i]], dtype=object)
    changed = int(sum(1 for o, n in zip(var_names[atac_i], atac_names) if o != n))
    if changed != len(atac_i):
        raise SystemExit(f"peak rename changed {changed}/{len(atac_i)} names — unexpected format")
    log(f"renamed all {changed:,} peaks  chr1-9776-10668 -> {atac_names[0]}")

    # --- stream the counts layer ----------------------------------------------------
    src = f["layers/counts"] if "counts" in f.get("layers", {}) else f["X"]
    log(f"reading {'layers/counts' if src is not f.get('X') else 'X'} in blocks of {BLOCK}")
    indptr, data, indices = src["indptr"], src["data"], src["indices"]
    rna_blocks, atac_blocks = [], []
    for s in range(0, n_cells, BLOCK):
        e = min(s + BLOCK, n_cells)
        lo, hi = int(indptr[s]), int(indptr[e])
        d = data[lo:hi]; j = indices[lo:hi]
        ip = indptr[s:e + 1] - lo
        blk = sp.csr_matrix((d, j, ip), shape=(e - s, n_feat))
        rna_blocks.append(blk[:, :SPLIT])
        atac_blocks.append(blk[:, SPLIT:])
        if (s // BLOCK) % 5 == 0:
            log(f"  rows {s:,}-{e:,}")
    f_obs = {}
    for c in ("batch", "donor", "site", "cell_type"):
        key = {"donor": "DonorNumber", "site": "Site"}.get(c, c)
        f_obs[c] = cat(f, key)
    f.close()

    RNA = sp.vstack(rna_blocks).tocsr(); del rna_blocks
    ATAC = sp.vstack(atac_blocks).tocsr(); del atac_blocks
    log(f"RNA  {RNA.shape} nnz={RNA.nnz:,}")
    log(f"ATAC {ATAC.shape} nnz={ATAC.nnz:,}")

    write_h5ad(OUT / "bmmc_rna_counts.h5ad", RNA, obs_names, var_names[gex_i], f_obs)
    write_h5ad(OUT / "bmmc_atac_counts.h5ad", ATAC, obs_names, atac_names, f_obs)

    rep = {
        "source_h5ad": str(a.h5ad),
        "n_cells": n_cells, "split_column": SPLIT,
        "rna": {"shape": list(RNA.shape), "nnz": int(RNA.nnz)},
        "atac": {"shape": list(ATAC.shape), "nnz": int(ATAC.nnz)},
        "obs_index_sha256_16": digest(obs_names),
        "atac_var_sha256_16": digest(atac_names),
        "rna_var_sha256_16": digest(var_names[gex_i]),
        "peak_name_format": str(atac_names[0]),
        "counts_layer_used": "layers/counts",
    }
    (OUT / "b01_export_report.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote bmmc_rna_counts.h5ad, bmmc_atac_counts.h5ad, b01_export_report.json -> {OUT}")
    log(f"ATAC var digest {rep['atac_var_sha256_16']} — this is the chain's identity; "
        f"assert it unchanged for every arm")

if __name__ == "__main__":
    main()
