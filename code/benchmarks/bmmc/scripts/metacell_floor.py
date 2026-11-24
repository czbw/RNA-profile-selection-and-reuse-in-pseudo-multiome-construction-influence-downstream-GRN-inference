#!/usr/bin/env python
"""Metacell preflight check, run before a metacell arm at metacell size k.

Fact 1: every shared cell type holds at least k cells in both modalities (otherwise SCENIC+ raises
"Sample larger than population"). Fact 2: the GEX cell ids are the frozen chain's analysis set, or
that set minus whole dropped cell types. Exit codes: 0 ok, 2 unusable input, 3 floor violated,
4 cell-id or cell-set violation. Adapted from the AML bundle's preflight script.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AMLGRN_ROOT", str(HERE.parent)))

def log(m):
    print(f"[b_12] {m}", flush=True)

def err(m):
    print(f"[b_12] {m}", file=sys.stderr, flush=True)

def load_config():
    name = os.environ.get("AMLGRN_CONFIG", "config.yaml")
    p = HERE / name
    if not p.is_file() or p.stat().st_size == 0:
        err(f"config {p} is missing or empty (AMLGRN_CONFIG={name})")
        raise SystemExit(2)
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or "paths" not in cfg:
        err(f"{p} did not parse into a mapping with a `paths:` block")
        raise SystemExit(2)
    cfg["P"] = {k: (ROOT / v) for k, v in cfg["paths"].items()}
    cfg["CFG_PATH"] = p
    return cfg

def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)

def _index_of(g) -> np.ndarray:
    key = g.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return _decode(g[str(key)][:])

def _obs_col(f, col) -> np.ndarray:
    node = f[f"obs/{col}"]
    if isinstance(node, h5py.Group):
        cats = _decode(node["categories"][:])
        return np.array([cats[c] for c in node["codes"][:]], dtype=object)
    legacy = f"obs/__categories/{col}"
    if legacy in f:
        cats = _decode(f[legacy][:])
        return np.array([cats[c] for c in node[:]], dtype=object)
    return _decode(node[:])

def atac_labels(a, cfg, group_by, sep):
    """ATAC-side per-cell labels, keyed on the bare h5ad barcode."""
    OUTP = cfg["P"]["out"]

    if a.cistopic:
        p = Path(a.cistopic)
        if not p.is_file() or p.stat().st_size == 0:
            return None, f"--cistopic {p} is missing or empty"
        try:
            import pickle
            with open(p, "rb") as fh:
                obj = pickle.load(fh)
            cd = obj.cell_data
            col = next((c for c in (group_by, "cell_type") if c in cd.columns), None)
            if col is None:
                return None, (f"{p.name} cell_data has neither `{group_by}` nor `cell_type` "
                              f"(columns {list(cd.columns)})")
            s = pd.Series(cd[col].astype(str).values,
                          index=[str(i).split(sep)[0] for i in cd.index])
            del obj
            return s, f"{p} (cell_data/{col})"
        except Exception as e:
            return None, f"{p.name} could not be unpickled ({type(e).__name__}: {e})"

    p = Path(a.atac_h5ad) if a.atac_h5ad else OUTP / "bmmc_atac_counts.h5ad"
    if p.is_file() and p.stat().st_size > 0:
        with h5py.File(p, "r") as f:
            ids = _index_of(f["obs"])
            col = next((c for c in (group_by, "cell_type") if f"obs/{c}" in f), None)
            if col is None:
                return None, (f"{p.name} obs has neither `{group_by}` nor `cell_type` "
                              f"(columns {list(f['obs'].keys())})")
            lab = _obs_col(f, col)
        return (pd.Series(lab, index=[str(i).split(sep)[0] for i in ids]),
                f"{p} (obs/{col})")

    if a.cell_data:
        q = Path(a.cell_data)
        if not q.is_file() or q.stat().st_size == 0:
            return None, f"--cell-data {q} is missing or empty"
        cd = pd.read_csv(q, index_col=0, dtype=str)
        col = next((c for c in (group_by, "cell_type") if c in cd.columns), None)
        if col is None:
            return None, f"{q.name} has neither `{group_by}` nor `cell_type`"
        return (pd.Series(cd[col].astype(str).values,
                          index=[str(i).split(sep)[0] for i in cd.index]),
                f"{q} ({col})")

    names_p = Path(a.analysis_set) if a.analysis_set else OUTP / f"cistopic_cell_names_{a.stage}.txt"
    meta_p = Path(a.metadata) if a.metadata else None
    if names_p.is_file() and names_p.stat().st_size > 0 and meta_p and meta_p.is_file():
        names = [ln.strip() for ln in names_p.read_text(encoding="utf-8").splitlines()
                 if ln.strip()]
        base = [n.split(sep)[0] for n in names]
        meta = pd.read_csv(meta_p, dtype=str)
        if "h5ad_barcode" not in meta.columns or "cell_type" not in meta.columns:
            return None, f"{meta_p.name} needs h5ad_barcode and cell_type"
        meta = meta.set_index("h5ad_barcode")
        missing = [b for b in base if b not in meta.index]
        if missing:
            return None, (f"{len(missing)} chain cell id(s) have no row in {meta_p.name} "
                          f"(e.g. {missing[:3]})")
        return (pd.Series(meta.loc[base, "cell_type"].astype(str).values, index=base),
                f"{names_p} joined to {meta_p} (cell_type)")

    return None, ("no ATAC source: give --cistopic, or --atac-h5ad, or --cell-data, or "
                  "--analysis-set together with --metadata")

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gex", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--cistopic", default=None, help="cistopic_obj_<STAGE>.pkl (authoritative, "
                                                     "needs pycisTopic; opt-in)")
    ap.add_argument("--atac-h5ad", default=None, help="default out/bmmc_atac_counts.h5ad")
    ap.add_argument("--cell-data", default=None)
    ap.add_argument("--analysis-set", default=None,
                    help="out/cistopic_cell_names_<STAGE>.txt; enables FACT 2")
    ap.add_argument("--metadata", default=None, help="uploads/b00_cell_metadata.csv")
    ap.add_argument("--require-analysis-set", action="store_true",
                    help="exit 2 if --analysis-set was not supplied or does not exist, instead "
                         "of reporting FACT 2 as not-checked")
    a = ap.parse_args()

    cfg = load_config()
    SEP = str(cfg["cell_id_sep"])
    pm = cfg["pairing_params"]["pseudomultiome"]
    GROUP_BY = str(pm.get("group_by", cfg["pseudobulk"]["group_by"]))
    if a.k < 2:
        err(f"metacell size must be >= 2 (got {a.k})")
        return 2

    gex_fp = Path(a.gex)
    if not gex_fp.is_file() or gex_fp.stat().st_size == 0:
        err(f"MISSING {gex_fp}")
        return 2

    import anndata as ad
    A = ad.read_h5ad(gex_fp)
    if GROUP_BY not in A.obs.columns:
        err(f"{gex_fp.name} has no '{GROUP_BY}' column; columns = {list(A.obs.columns)}")
        return 2
    rna_n = A.obs[GROUP_BY].astype(str).value_counts()
    gex_ids = [str(x) for x in A.obs_names]

    atac_s, atac_src = atac_labels(a, cfg, GROUP_BY, SEP)
    if atac_s is None:
        err(f"WARNING: ATAC side NOT checked -- {atac_src}")
        atac_n = None
    else:
        atac_n = atac_s.astype(str).value_counts()
        log(f"ATAC side from {atac_src}")

    labels = sorted(set(rna_n.index) if atac_n is None else set(rna_n.index) & set(atac_n.index))
    log(f"{a.stage}: k={a.k}  GEX {A.n_obs:,} cells x {A.n_vars:,} genes  "
        f"| {len(labels)} shared cell types"
        f"{' (RNA only -- ATAC UNCHECKED)' if atac_n is None else ''}")

    # ---- FACT 1: the k floor ----------------------------------------------------------------
    bad = []
    for ct in labels:
        r = int(rna_n.get(ct, 0))
        t = int(atac_n.get(ct, 0)) if atac_n is not None else None
        worst = r if t is None else min(r, t)
        flag = "  <-- BELOW FLOOR" if worst < a.k else ""
        print(f"[b_12]   {ct:26s} RNA={r:>6,}"
              + (f"  ATAC={t:>6,}" if t is not None else "") + f"  min={worst:>6,}{flag}",
              flush=True)
        if worst < a.k:
            bad.append((ct, r, t))

    if bad:
        err(f"REFUSING to run at k={a.k}: {len(bad)} cell type(s) hold fewer than {a.k} cells "
            f"and would raise 'Sample larger than population' inside SCENIC+:")
        for ct, r, t in bad:
            err(f"  {ct}: RNA={r}" + (f", ATAC={t}" if t is not None else ""))
        err("Rebuild the metacell GEX with a floor of at least this k "
            "(filter_gex_for_metacells.py --k), or lower k.")
        return 3
    thin = min((int(rna_n.get(ct, 0)) if atac_n is None
                else min(int(rna_n.get(ct, 0)), int(atac_n.get(ct, 0)))) for ct in labels)
    log(f"FACT 1 OK: every shared cell type clears the {a.k}-cell floor "
        f"(thinnest holds {thin:,}).")

    # ---- FACT 2: the cell-id convention -----------------------------------------------------
    as_p = Path(a.analysis_set) if a.analysis_set else None
    if as_p is None or not as_p.is_file() or as_p.stat().st_size == 0:
        msg = (f"FACT 2 NOT CHECKED: no usable --analysis-set"
               f"{f' ({as_p})' if as_p else ''}. The GEX cell ids were not compared against the "
               f"frozen chain. Pass --analysis-set out/cistopic_cell_names_{a.stage}.txt on "
               f"CUBIC, where the chain exists.")
        if a.require_analysis_set:
            err(msg)
            return 2
        err("WARNING: " + msg)
        log("b_12 OK (FACT 1 only).")
        return 0

    names = [ln.strip() for ln in as_p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(set(names)) != len(names):
        err(f"{as_p.name} contains duplicate cell ids")
        return 2
    nameset = set(names)
    foreign = [i for i in gex_ids if i not in nameset]
    n_ok = len(gex_ids) - len(foreign)
    log(f"FACT 2a: {n_ok:,}/{len(gex_ids):,} GEX cell ids are ids of the analysis set "
        f"({as_p.name}, {len(names):,} ids)")
    if foreign:
        err(f"CELL-ID CONVENTION VIOLATED: {len(foreign):,} GEX cell id(s) are NOT in the frozen "
            f"analysis set. Examples GEX={foreign[:3]} vs chain={names[:3]}.")
        if n_ok == 0:
            err("The overlap is ZERO -- a convention that drifts "
                "by one separator matches no cells at all, silently. Rebuild the GEX with "
                "native_rna.py on a host where "
                f"out/cistopic_cell_names_{a.stage}.txt exists, then re-run b_11.")
        return 4

    absent = [n for n in names if n not in set(gex_ids)]
    if absent:
        kept = set(rna_n.index)
        if atac_s is None:
            err(f"FACT 2b NOT CHECKED: {len(absent):,} analysis-set id(s) are absent from the "
                f"GEX and no ATAC label source is available to confirm they belong to dropped "
                f"cell types.")
            return 2
        base_absent = [n.split(SEP)[0] for n in absent]
        unknown = [b for b in base_absent if b not in atac_s.index]
        if unknown:
            err(f"{len(unknown):,} absent analysis-set id(s) have no ATAC label "
                f"(e.g. {unknown[:3]})")
            return 4
        wrongly_absent = sorted({str(atac_s.loc[b]) for b in base_absent} & kept)
        if wrongly_absent:
            err(f"CELL-SET VIOLATED: {len(absent):,} analysis-set cell(s) are missing from the "
                f"GEX although their cell type was KEPT: {wrongly_absent[:5]}. The floor filter "
                f"drops whole cell types, never individual cells.")
            return 4
        log(f"FACT 2b OK: {len(absent):,} analysis-set cell(s) are absent, all of dropped cell "
            f"type(s) {sorted({str(atac_s.loc[b]) for b in base_absent})}")
    else:
        log(f"FACT 2b OK: the GEX covers the analysis set exactly ({len(names):,} cells); "
            f"the floor filter dropped nothing")

    log(f"b_12 OK: every shared cell type clears the {a.k}-cell floor and the GEX cell ids are "
        f"the frozen analysis set.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
