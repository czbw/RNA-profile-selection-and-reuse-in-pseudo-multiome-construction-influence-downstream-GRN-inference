"""Native symbol-indexed MALIGNANT RNA per stage (Option B input)."""
import anndata as ad

from common import config_stages, load_config, log, rna_prefix, stage_key

cfg = load_config()
OUT = cfg["P"]["out"]
RNA_DIR = cfg["P"]["rna_h5ad_dir"]
STAGES = config_stages(cfg)
PREFIX = rna_prefix(cfg)          # "malignant"
OUT.mkdir(parents=True, exist_ok=True)

KEEP_OBS = ["Classified_Celltype", "Patient_ID", "Subgroup",
            "Expected_Driving_Aberration", "Patient_Sample", "Library_ID"]

def build(stage):
    fp = RNA_DIR / f"{PREFIX}_{stage_key(stage)}.h5ad"
    if not fp.exists():
        log(f"{stage}: RNA h5ad missing ({fp}) — run build_malignant_h5ad.py & stage it — skip")
        return
    A = ad.read_h5ad(fp)
    A.uns.clear()                     # defensive: drop any uns (e.g. log1p null) before re-write
    A.var["ensembl"] = A.var_names
    sym = A.var["symbol"].astype(str)
    A = A[:, sym.notna() & (sym != "nan")].copy()
    A.var_names = A.var["symbol"].astype(str)
    A.var_names_make_unique()
    A.var = A.var[["ensembl", "symbol"]].astype(str)
    keep = [c for c in KEEP_OBS if c in A.obs.columns]
    A.obs = A.obs[keep].copy()
    A.obs = A.obs.rename(columns={"Expected_Driving_Aberration": "Driver"})
    A.obs["timepoint"] = stage
    for c in A.obs.columns:
        if A.obs[c].dtype == object:
            A.obs[c] = A.obs[c].astype(str)
    A.write_h5ad(OUT / f"native_rna_{stage}.h5ad")
    log(f"{stage}: native malignant RNA {A.n_obs:,} x {A.n_vars:,} -> {OUT / f'native_rna_{stage}.h5ad'}  "
        f"celltypes={A.obs['Classified_Celltype'].nunique()}")

if __name__ == "__main__":
    for s in STAGES:
        build(s)
    log("Option B GEX ready. (Option A GEX = paired_rna_{stage}.h5ad from build_paired_rna.py.)")
