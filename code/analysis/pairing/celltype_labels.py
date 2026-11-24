"""FIG1 disclosure - is the AML ATAC cell-type label the link's own output?"""
import json
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
RNA = ROOT / "paml/malignant" / "data" / "rna"
OUT = ROOT / "analysis" / "out"

STAGES = {"Diagnosis": "malignant_diagnosis.h5ad", "Relapse": "malignant_relapse.h5ad"}

def log(m):
    print(f"[f1-labelprov] {m}", flush=True)

res = {
    "status": ("corroboration of an earlier result, computed against the MALIGNANT RNA object, not a "
               "reproduction of that result's own denominator - see the module docstring"),
    "mechanism_source": ("Longitudinal_pAML_paper/Data_processing/Annotate_ATAC.R:34-41 - the ATAC "
                         "cell type is copied from the single best-matching RNA cell "
                         "(`predictedCell`), which is the object deposited as "
                         "`Linked_scRNA_barcode` and is this project's pairing source"),
    "join_key": "(Linked_scRNA_Sample, Linked_scRNA_barcode) -> (Library_ID, barcode)",
    "by_stage": {},
}

for stage, fn in STAGES.items():
    cd = pd.read_csv(S2 / f"cell_data_{stage}.csv", index_col=0)
    a = ad.read_h5ad(RNA / fn, backed="r")
    obs = a.obs[["Library_ID", "Classified_Celltype"]].copy()
    # obs_names are "<GSM>_<Library_ID>|<barcode>"
    obs["barcode"] = [str(x).split("|", 1)[1] for x in a.obs_names]

    # Barcode alone is ambiguous, so it may not be used as the key (asserted below).
    n_bc = obs["barcode"].nunique()
    ambiguous = len(obs) - n_bc
    assert ambiguous > 0, (
        "barcode alone was expected to be ambiguous across RNA libraries; it is not, which "
        "means this object is not what the join assumed and the composite key must be re-derived")

    rna_ct = dict(zip(zip(obs["Library_ID"].astype(str), obs["barcode"].astype(str)),
                      obs["Classified_Celltype"].astype(str)))

    key = list(zip(cd["Linked_scRNA_Sample"].astype(str),
                   cd["Linked_scRNA_barcode"].astype(str)))
    partner = [rna_ct.get(k) for k in key]
    resolvable = [p is not None for p in partner]
    n_atac = len(cd)
    n_res = sum(resolvable)
    agree = sum(1 for p, c in zip(partner, cd["Classified_Celltype"].astype(str))
                if p is not None and p == c)

    # cross-timepoint: the linked partner's library is from the other timepoint
    tp_tok = "_DX" if stage == "Diagnosis" else "_REL"
    cross = sum(1 for s in cd["Linked_scRNA_Sample"].astype(str) if not s.upper().endswith(tp_tok))

    res["by_stage"][stage] = {
        "n_atac_cells": int(n_atac),
        "n_rna_cells_in_malignant_object": int(len(obs)),
        "barcode_only_collisions": int(ambiguous),
        "barcode_only_is_ambiguous": True,
        "n_partner_resolvable_in_malignant_object": int(n_res),
        "frac_partner_resolvable": n_res / n_atac,
        "n_label_agrees": int(agree),
        "celltype_label_agreement_among_resolvable": agree / n_res if n_res else None,
        "n_partner_unresolvable": int(n_atac - n_res),
        "frac_partner_unresolvable": (n_atac - n_res) / n_atac,
        "frac_linked_partner_from_other_timepoint": cross / n_atac,
    }
    s = res["by_stage"][stage]
    log(f"{stage}: {n_atac} ATAC cells | partner resolvable {n_res} "
        f"({s['frac_partner_resolvable']:.4f}) | label agreement among those "
        f"{s['celltype_label_agreement_among_resolvable']:.6f} | "
        f"cross-timepoint {s['frac_linked_partner_from_other_timepoint']:.4f} | "
        f"barcode-only collisions {ambiguous}")

with open(OUT / "f1_label_provenance.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote f1_label_provenance.json")
