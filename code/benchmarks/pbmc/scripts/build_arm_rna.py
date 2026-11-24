"""Build the paired RNA object for a random arm. RUNS ON CUBIC."""
import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

ROOT = (Path(__file__).resolve().parents[3] / 'benchmarks/pbmc')

def log(m):
    print(f"[p12] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True)
ap.add_argument("--dataset", required=True)
ap.add_argument("--arm", required=True)
ap.add_argument("--assignment", required=True)
a = ap.parse_args()

OUT = ROOT / "out"
WORK = ROOT / "work"

G = sc.read_h5ad(WORK / a.dataset / "gex_annotated.h5ad")
log(f"GEX {G.n_obs:,} x {G.n_vars:,} (obs_names like {G.obs_names[0]})")

cd = pd.read_csv(WORK / f"cell_data_{a.stage}.csv")
idcol = "cistopic_id" if "cistopic_id" in cd.columns else "cell_id"
cd = cd.set_index(idcol)
log(f"cell_data {len(cd):,} rows, indexed by {idcol}")

asg = pd.read_csv(a.assignment)
log(f"assignment {len(asg):,} rows from {Path(a.assignment).name}")

# --- the cell-set assertion: identical to what native reference/link were inferred on --------
missing = set(cd.index) - set(asg["cistopic_id"])
extra = set(asg["cistopic_id"]) - set(cd.index)
if missing or extra:
    raise SystemExit(
        f"CELL SET MISMATCH -- refusing to build.\n"
        f"  in cell_data but not in the assignment : {len(missing)}\n"
        f"  in the assignment but not in cell_data : {len(extra)}\n"
        f"  truth/link were inferred on all {len(cd):,} cells. Rebuild the assignment with\n"
        f"  make_random_arms.py --min-stratum 0 (writes *_full.csv).")
log(f"cell-set assertion PASSED: assignment covers all {len(cd):,} manifest cells")

asg = asg.set_index("cistopic_id").loc[cd.index]          # align to manifest order

# bare barcode -> GEX obs_name, exactly as p_03 does it
bc2id = dict(zip(G.obs["barcode"].astype(str), G.obs_names.astype(str)))
unknown = set(asg["Linked_scRNA_barcode"]) - set(bc2id)
if unknown:
    raise SystemExit(f"{len(unknown)} assigned barcodes absent from the GEX object, e.g. "
                     f"{sorted(unknown)[:5]}")
rna_ids = asg["Linked_scRNA_barcode"].map(bc2id).values

sub = G[rna_ids].copy()                   # duplicates intentional -- reuse is the phenomenon
sub.obs_names = cd.index.values
sub.obs = cd.copy()
fp = OUT / f"paired_rna_{a.stage}.{a.arm}.h5ad"
sub.write_h5ad(fp)

n_distinct = len(set(rna_ids))
log(f"{a.arm}: {sub.n_obs:,} x {sub.n_vars:,} | distinct RNA {n_distinct:,} "
    f"| reuse {sub.n_obs / n_distinct:.2f}x -> {fp.name}")
log("STEP 12 COMPLETE")
