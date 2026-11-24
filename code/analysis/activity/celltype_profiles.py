"""FIG1 / S1.e  +  FIG6 / (c) - cell-type distribution of each arm's regulon activity."""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
L3 = ROOT / "longitudinal/published" / "out"
R3 = ROOT / "longitudinal/corrected" / "out"
OUT = ROOT / "analysis" / "out"

STAGES = ["Diagnosis", "Relapse"]
ARMS = ["link", "pseudomultiome", "repair"]
MIN_UNITS = 20

def log(m):
    print(f"[f1-ct] {m}", flush=True)

def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

def read_gene_based(path):
    with h5py.File(path, "r") as h:
        cells = dec(h["obs/Cell"][:])
        var = dec(h["mod/Gene_based/var/_index"][:])
        X = h["mod/Gene_based/X"][:]

    reg = [v.rsplit("_(", 1)[0] for v in var]
    return cells, reg, X

# module definitions
fl = pd.read_csv(L3 / "nomination_funnel.csv", index_col=0)
fr = pd.read_csv(R3 / "nomination_funnel.csv", index_col=0)
rec = json.load(open(R3 / "s3_recalibrated.json"))
MODULES = {
    "link_T4_26": list(fl.index[fl.survives_4]),
    "repair_core34": list(fr.index[fr.ax4_tier1 & fr.ax1_directional & fr.ax2_selective
                                   & fr.ax3_celltype_robust]),
    "repair_LMM26": rec["results"]["E1_random_intercept"]["confirmed_regulons"],
}

res = {"status": "DESCRIPTIVE / EXPLORATORY", "prespec": "analysis/FIG1_FIG6_PRESPEC.md S1.e",
       "join": {}, "celltype_counts": {}, "excluded_below_floor": {}, "modules": {k: len(v) for k, v in MODULES.items()}}
prof_rows, mod_rows = [], []

for stage in STAGES:
    cd = pd.read_csv(S2 / f"cell_data_{stage}.csv", index_col=0)
    ct_map = cd["Classified_Celltype"].to_dict()
    known_ct = set(cd["Classified_Celltype"].unique())
    for arm in ARMS:
        cells, reg, X = read_gene_based(S2 / f"eRegulon_AUC_{stage}_{arm}.h5mu")
        if arm == "pseudomultiome":
            ct = np.array([c.rsplit("_", 1)[0] for c in cells])
            unmatched = sorted(set(ct) - known_ct)
            res["join"][f"{stage}_{arm}"] = {
                "method": "metacell name parse `<CellType>_<i>`",
                "n_units": len(cells), "n_celltypes_parsed": int(len(set(ct))),
                "parsed_labels_not_in_cell_data_vocabulary": unmatched,
                "match_rate": float(np.mean([c in known_ct for c in ct])),
            }
        else:
            ct = np.array([ct_map.get(c, "__UNMATCHED__") for c in cells])
            mr = float(np.mean(ct != "__UNMATCHED__"))
            res["join"][f"{stage}_{arm}"] = {
                "method": "composite key <barcode>___<GSM>_<sample> vs cell_data index",
                "n_units": len(cells), "n_matched": int((ct != "__UNMATCHED__").sum()),
                "match_rate": mr,
            }
            assert mr == 1.0, f"{stage}/{arm}: composite-key join match rate {mr}, expected 1.000"
        log(f"{stage}/{arm}: {len(cells)} units, {X.shape[1]} regulons, "
            f"{len(set(ct))} cell types, match {res['join'][f'{stage}_{arm}']['match_rate']:.4f}")

        cts, counts = np.unique(ct, return_counts=True)
        res["celltype_counts"][f"{stage}_{arm}"] = {a: int(b) for a, b in zip(cts, counts)}
        keep = [c for c, n in zip(cts, counts) if n >= MIN_UNITS]
        res["excluded_below_floor"][f"{stage}_{arm}"] = {
            a: int(b) for a, b in zip(cts, counts) if b < MIN_UNITS}

        M = pd.DataFrame(X, columns=reg)
        M["__ct__"] = ct
        prof = M[M.__ct__.isin(keep)].groupby("__ct__").mean()
        # within-arm z across cell types, per regulon
        z = (prof - prof.mean(axis=0)) / prof.std(axis=0, ddof=0).replace(0, np.nan)
        for c in prof.index:
            for r in prof.columns:
                prof_rows.append({"stage": stage, "arm": arm, "celltype": c, "regulon": r,
                                  "n_units": int(res["celltype_counts"][f"{stage}_{arm}"][c]),
                                  "mean_AUC_raw": float(prof.loc[c, r]),
                                  "z_within_arm_across_celltypes": float(z.loc[c, r])})
        # module profiles: mean over the module's regulons that this arm actually has
        for mname, mregs in MODULES.items():
            present = [r for r in mregs if r in prof.columns]
            if not present:
                continue
            pr = prof[present].mean(axis=1)
            zr = z[present].mean(axis=1)
            for c in pr.index:
                mod_rows.append({"stage": stage, "arm": arm, "module": mname,
                                 "n_module_regulons_in_arm": len(present),
                                 "n_module_regulons_total": len(mregs),
                                 "celltype": c,
                                 "n_units": int(res["celltype_counts"][f"{stage}_{arm}"][c]),
                                 "module_mean_AUC_raw": float(pr[c]),
                                 "module_mean_z_within_arm": float(zr[c])})
        del X, M, prof, z

P = pd.DataFrame(prof_rows)
MOD = pd.DataFrame(mod_rows)
P.to_csv(OUT / "f1_celltype_profile.csv", index=False)
MOD.to_csv(OUT / "f1_celltype_module_profile.csv", index=False)

# cross-arm profile-shape agreement, on the module z-profiles over shared cell types
agree = {}
for stage in STAGES:
    for mname in MODULES:
        for A, B in [("link", "pseudomultiome"), ("link", "repair"), ("pseudomultiome", "repair")]:
            a = MOD[(MOD.stage == stage) & (MOD.arm == A) & (MOD.module == mname)
                    ].set_index("celltype").module_mean_z_within_arm
            b = MOD[(MOD.stage == stage) & (MOD.arm == B) & (MOD.module == mname)
                    ].set_index("celltype").module_mean_z_within_arm
            ix = a.index.intersection(b.index)
            if len(ix) >= 5:
                agree[f"{stage}|{mname}|{A}_vs_{B}"] = {
                    "n_shared_celltypes": int(len(ix)),
                    "pearson_r": float(np.corrcoef(a[ix], b[ix])[0, 1]),
                    "spearman_r": float(pd.Series(a[ix]).corr(pd.Series(b[ix]), method="spearman")),
                    "n_module_regulons_in_A": int(MOD[(MOD.stage == stage) & (MOD.arm == A)
                                                      & (MOD.module == mname)].n_module_regulons_in_arm.iloc[0]),
                    "n_module_regulons_in_B": int(MOD[(MOD.stage == stage) & (MOD.arm == B)
                                                      & (MOD.module == mname)].n_module_regulons_in_arm.iloc[0]),
                }
res["module_profile_shape_agreement"] = agree
res["NOTE_normalisation"] = ("raw AUC means are NOT comparable across arms - single cells vs 5-cell "
                             "metacells. Compare the z-scored profile SHAPE only; metacell averaging "
                             "compresses within-cell-type variance so spread is not comparable.")

with open(OUT / "f1_celltype.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote f1_celltype_profile.csv, f1_celltype_module_profile.csv, f1_celltype.json")
