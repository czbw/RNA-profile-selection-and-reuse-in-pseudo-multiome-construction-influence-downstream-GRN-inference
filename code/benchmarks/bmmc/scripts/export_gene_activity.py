#!/usr/bin/env python
"""Export the ATAC gene-activity matrix, the input the `link` arm needs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]

def log(m): print(f"[b03] {m}", flush=True)

def cat(f, col):
    codes = f[f"obs/{col}"][:]
    cats = [x.decode() if isinstance(x, bytes) else str(x)
            for x in f[f"obs/__categories/{col}"][:]]
    return np.array([cats[c] for c in codes], dtype=object)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default=str(ROOT / "benchmarks/bmmc" / "out"))
    a = ap.parse_args()
    OUT = Path(a.outdir); OUT.mkdir(parents=True, exist_ok=True)

    f = h5py.File(a.h5ad, "r")
    if "ATAC_gene_activity" not in f["obsm"]:
        raise SystemExit("obsm/ATAC_gene_activity is absent — the link arm cannot be built "
                         "from this file")
    g = f["obsm/ATAC_gene_activity"]
    shape = tuple(int(x) for x in g.attrs["shape"])
    enc = g.attrs.get("encoding-type", b"")
    enc = enc.decode() if isinstance(enc, bytes) else str(enc)
    log(f"obsm/ATAC_gene_activity: {shape[0]:,} x {shape[1]:,}  encoding={enc}")

    M = (sp.csc_matrix if enc == "csc_matrix" else sp.csr_matrix)(
        (g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
    M = M.tocsr()
    log(f"loaded, nnz={M.nnz:,}  density={M.nnz / (shape[0]*shape[1]):.4f}")

    obs_names = np.array([x.decode() for x in f["obs/_index"][:]], dtype=object)
    obs_cols = {}
    for c in ("batch", "donor", "site", "cell_type"):
        obs_cols[c] = cat(f, {"donor": "DonorNumber", "site": "Site"}.get(c, c))

    # Feature names are supplied, in uns/ATAC_gene_activity_var_names (not in obsm). Read them
    # and assert the width.
    gv = f["uns/ATAC_gene_activity_var_names"][:]
    var_names = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in gv],
                         dtype=object)
    if len(var_names) != shape[1]:
        raise SystemExit(f"uns/ATAC_gene_activity_var_names has {len(var_names)} entries "
                         f"but the matrix has {shape[1]} columns")
    uns_note = ("feature names from uns/ATAC_gene_activity_var_names; these ARE gene symbols "
                "and are joinable to the GEX vocabulary by intersection")
    log(f"gene-activity features named from uns: {len(var_names):,} symbols, "
        f"e.g. {list(var_names[:3])}")
    f.close()

    st = h5py.string_dtype("utf-8")
    p = OUT / "bmmc_atac_gene_activity.h5ad"
    with h5py.File(p, "w") as o:
        o.attrs["encoding-type"] = "anndata"; o.attrs["encoding-version"] = "0.1.0"
        X = o.create_group("X")
        X.attrs["encoding-type"] = "csr_matrix"; X.attrs["encoding-version"] = "0.1.0"
        X.attrs["shape"] = np.array(M.shape, dtype="int64")
        X.create_dataset("data", data=M.data, compression="gzip", compression_opts=1)
        X.create_dataset("indices", data=M.indices, compression="gzip", compression_opts=1)
        X.create_dataset("indptr", data=M.indptr, compression="gzip", compression_opts=1)
        for name, idx, cols in (("obs", obs_names, obs_cols), ("var", var_names, {})):
            d = o.create_group(name)
            d.attrs["encoding-type"] = "dataframe"; d.attrs["encoding-version"] = "0.1.0"
            d.attrs["_index"] = "_index"
            d.attrs.create("column-order", np.array(list(cols), dtype=object), dtype=st)
            d.create_dataset("_index", data=np.array(idx, dtype=object), dtype=st)
            for c, v in cols.items():
                d.create_dataset(c, data=np.array(v, dtype=object), dtype=st)
        u = o.create_group("uns")
        u.create_dataset("gene_activity_note", data=np.array(uns_note, dtype=object), dtype=st)

    (OUT / "b03_gene_activity_report.json").write_text(json.dumps({
        "shape": list(M.shape), "nnz": int(M.nnz), "source_encoding": enc,
        "feature_names": "gene symbols from uns/ATAC_gene_activity_var_names", "note": uns_note,
    }, indent=1), encoding="utf-8")
    log(f"wrote {p.name} -> {OUT}")

if __name__ == "__main__":
    main()
