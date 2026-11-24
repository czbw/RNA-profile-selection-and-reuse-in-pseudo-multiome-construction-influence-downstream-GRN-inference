"""Build paired RNA objects from cell_data_<stage>.repaired.csv (one per stage).

Writes paired_rna_<stage>.repair.h5ad and paired_cells_<stage>.repair.txt.
"""
import anndata as ad
import pandas as pd

from common import config_stages, load_config, log, rna_prefix, stage_key

cfg = load_config()
WORK = cfg["P"]["work"]
OUT = cfg["P"]["out"]
RNA_DIR = cfg["P"]["rna_h5ad_dir"]
STAGES = config_stages(cfg)
PREFIX = rna_prefix(cfg)          # "malignant"
OUT.mkdir(parents=True, exist_ok=True)

LIB = cfg["linkage"]["rna_sample_field"]       # Library_ID
MINCONF = float(cfg["linkage"]["min_confidence"])

def rna_lookup():
    """Map (Library_ID, Cell_Barcode) -> (stage, obs_name) across the malignant RNA h5ads.
    Uses backed mode: reads .obs only, never loads X."""
    rows = []
    for s in STAGES:
        fp = RNA_DIR / f"{PREFIX}_{stage_key(s)}.h5ad"
        if not fp.exists():
            log(f"  [WARN] RNA h5ad missing: {fp}")
            continue
        A = ad.read_h5ad(fp, backed="r")
        obs = A.obs
        barcode = pd.Index(obs.index).str.split("|").str[-1]   # obs_names = "{stem}|{barcode}"
        df = pd.DataFrame({"lib": obs[LIB].values, "barcode": barcode,
                           "obs_name": obs.index.values, "rna_stage": s})
        rows.append(df)
        A.file.close()
    if not rows:
        raise SystemExit("No malignant RNA h5ads found — run build_malignant_h5ad.py first.")
    allr = pd.concat(rows, ignore_index=True)
    allr["key"] = allr["lib"].astype(str) + "||" + allr["barcode"].astype(str)
    dup = allr["key"].duplicated().sum()
    if dup:
        log(f"  [note] {dup} duplicate (lib,barcode) keys in RNA — keeping first")
        allr = allr.drop_duplicates("key")
    return allr.set_index("key")

def build_stage(stage, lut):
    cd = pd.read_csv(WORK / f"cell_data_{stage}.repaired.csv", index_col=0)
    n0 = len(cd)
    conf = pd.to_numeric(cd["Linked_scRNA_confidence_score"], errors="coerce")
    cd = cd[conf >= MINCONF].copy()
    cd["key"] = (cd["Linked_scRNA_Sample"].astype(str) + "||"
                 + cd["Linked_scRNA_barcode"].astype(str))
    cd["rna_obs"] = cd["key"].map(lut["obs_name"])
    cd["rna_stage"] = cd["key"].map(lut["rna_stage"])
    matched = cd[cd["rna_obs"].notna()].copy()
    log(f"{stage}: ATAC malignant {n0:,} -> conf>={MINCONF} {len(cd):,} -> RNA-matched "
        f"{len(matched):,} ({len(matched)/max(1,n0):.1%}). "
        f"linked-RNA-stage mix: {matched['rna_stage'].value_counts().to_dict()}")
    if matched.empty:
        return

    keep_obs = ["sample_id", "Classified_Celltype", "Patient_ID", "Subgroup", "Driver",
                "timepoint", "Linked_scRNA_confidence_score"]
    parts = []
    for src_stage, grp in matched.groupby("rna_stage"):
        A = ad.read_h5ad(RNA_DIR / f"{PREFIX}_{stage_key(src_stage)}.h5ad", backed="r")
        sub = A[grp["rna_obs"].values].to_memory()
        A.file.close()
        sub.uns.clear()                   # defensive: drop any uns (e.g. log1p null) before re-write
        sub.obs_names = grp.index.values          # ATAC cell_id (1:1 with grp order)
        sub.obs = grp[keep_obs].copy()
        sub.var["ensembl"] = sub.var_names
        sym = sub.var["symbol"].astype(str)
        sub = sub[:, sym.notna() & (sym != "nan")].copy()
        sub.var_names = sub.var["symbol"].astype(str)
        sub.var_names_make_unique()
        parts.append(sub)
    paired = parts[0] if len(parts) == 1 else ad.concat(parts, join="outer", fill_value=0,
                                                        merge="first")
    paired.var = paired.var[[c for c in ["ensembl", "symbol"] if c in paired.var]].astype(str)
    for c in paired.obs.columns:
        if paired.obs[c].dtype == object:
            paired.obs[c] = paired.obs[c].astype(str)
    paired.write_h5ad(OUT / f"paired_rna_{stage}.repair.h5ad")
    pd.Series(paired.obs_names).to_csv(OUT / f"paired_cells_{stage}.repair.txt", index=False, header=False)
    log(f"  -> {OUT / f'paired_rna_{stage}.repair.h5ad'}  ({paired.n_obs:,} x {paired.n_vars:,})")

if __name__ == "__main__":
    lut = rna_lookup()
    log(f"RNA lookup: {len(lut):,} malignant RNA cells across {STAGES}")
    for s in STAGES:
        if (WORK / f"cell_data_{s}.repaired.csv").exists():
            build_stage(s, lut)
        else:
            log(f"{s}: no cell_data_{s}.repaired.csv in work/ — SKIPPED (stage it first)")
    log("next: the ATAC chain (atac_chain.py) then run_grn.sh")
