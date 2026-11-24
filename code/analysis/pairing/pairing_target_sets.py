"""FIG1 / S1.0 + S1.a - set-level agreement and target-set agreement, link vs pseudo-multiome."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
OUT = ROOT / "analysis" / "out"
OUT.mkdir(parents=True, exist_ok=True)

STAGES = ["Diagnosis", "Relapse"]
ARMS = ["link", "pseudomultiome", "repair"]
NEST_BAR = 0.80          # pre-specified threshold
JACCARD_BAR = 0.50       

USECOLS = ["Region", "Gene", "TF", "is_extended", "eRegulon_name"]

def log(m):
    print(f"[f1-sets] {m}", flush=True)

def load(stage, arm):
    p = S2 / f"eRegulons_{stage}_{arm}.tsv"
    d = pd.read_csv(p, sep="\t", usecols=USECOLS, low_memory=False)
    # every Stage-2 row is `direct` (results.json /structure: extended == 0 in all four units).
    # asserted, because an `extended` row would double-count a regulon name.
    n_ext = int((d.is_extended.astype(str).str.lower() == "true").sum())
    assert n_ext == 0, f"{p.name}: {n_ext} extended rows, unexpected"
    return p, d

def jac(a, b):
    a, b = set(a), set(b)
    u = len(a | b)
    return (len(a & b) / u) if u else np.nan

def summarise(v):
    v = np.asarray(v, dtype=float)
    return {
        "n": int(v.size),
        "median": float(np.median(v)) if v.size else None,
        "q25": float(np.percentile(v, 25)) if v.size else None,
        "q75": float(np.percentile(v, 75)) if v.size else None,
        "mean": float(v.mean()) if v.size else None,
        "frac_gt_0.5": float((v > 0.5).mean()) if v.size else None,
        "frac_lt_0.1": float((v < 0.1).mean()) if v.size else None,
        "n_ge_bar_0.50": int((v >= JACCARD_BAR).sum()) if v.size else 0,
    }

res = {
    "status": "DESCRIPTIVE / EXPLORATORY - inherits the post-hoc label of both AML arms",
    "prespec": "analysis/FIG1_FIG6_PRESPEC.md S1.0, S1.a",
    "bars": {"nesting_TF_containment": NEST_BAR, "target_jaccard": JACCARD_BAR},
    "inputs": {},
    "S1_0_sets": {},
    "S1_a_targets": {},
    "crosscheck_vs_results_json": {},
}

tables = {}
for stage in STAGES:
    for arm in ARMS:
        p, d = load(stage, arm)
        tables[(stage, arm)] = d
        res["inputs"][f"{stage}_{arm}"] = {
            "path": str(p.relative_to(ROOT)).replace("\\", "/"),
            "rows": int(len(d)),
            "n_regulons": int(d.eRegulon_name.nunique()),
            "n_TFs": int(d.TF.nunique()),
            "n_genes_union": int(d.Gene.nunique()),
            "n_regions_union": int(d.Region.nunique()),
        }
        log(f"{stage}/{arm}: {len(d)} rows, {d.eRegulon_name.nunique()} regulons, "
            f"{d.TF.nunique()} TFs, {d.Gene.nunique()} genes")

# ---------------------------------------------------------------- S1.0 set agreement
PAIRS = [("link", "pseudomultiome"), ("link", "repair"), ("pseudomultiome", "repair")]
member_rows = []

for stage in STAGES:
    for A, B in PAIRS:
        dA, dB = tables[(stage, A)], tables[(stage, B)]
        block = {}
        for level, col in (("TF", "TF"), ("regulon", "eRegulon_name")):
            sA, sB = set(dA[col]), set(dB[col])
            inter = sA & sB
            block[level] = {
                "n_A": len(sA), "n_B": len(sB), "n_shared": len(inter),
                "n_union": len(sA | sB),
                "containment_B_in_A": len(inter) / len(sB) if sB else None,
                "containment_A_in_B": len(inter) / len(sA) if sA else None,
                "jaccard": len(inter) / len(sA | sB) if (sA | sB) else None,
                "n_unique_A": len(sA - sB), "n_unique_B": len(sB - sA),
                "unique_A": sorted(sA - sB) if level == "TF" else None,
                "unique_B": sorted(sB - sA) if level == "TF" else None,
            }
        res["S1_0_sets"][f"{stage}__{A}_vs_{B}"] = block

    # long-form membership table across all three arms, both levels
    for level, col in (("TF", "TF"), ("regulon", "eRegulon_name")):
        universe = set()
        for arm in ARMS:
            universe |= set(tables[(stage, arm)][col])
        for m in sorted(universe):
            row = {"stage": stage, "level": level, "member": m}
            for arm in ARMS:
                row[f"in_{arm}"] = m in set(tables[(stage, arm)][col])
            member_rows.append(row)

pd.DataFrame(member_rows).to_csv(OUT / "f1_set_members.csv", index=False)

# verdict
verdict = {}
fired = True
for stage in STAGES:
    c = res["S1_0_sets"][f"{stage}__link_vs_pseudomultiome"]["TF"]["containment_B_in_A"]
    verdict[stage] = c
    fired = fired and (c >= NEST_BAR)
res["S1_0_sets"]["VERDICT"] = {
    "bar": NEST_BAR,
    "containment_pseudomultiome_in_link_TF_level": verdict,
    "NESTED": bool(fired),
    "reading": ("pseudo-multiome is LARGELY NESTED inside link at TF level at both timepoints"
                if fired else
                "pseudo-multiome is NOT largely nested inside link at TF level at both timepoints"),
}
log(f"VERDICT nested={fired}  {verdict}")

# ---------------------------------------------------------------- S1.a target agreement
reg_rows, tf_rows = [], []
for stage in STAGES:
    for A, B in PAIRS:
        dA, dB = tables[(stage, A)], tables[(stage, B)]
        gA = dA.groupby("eRegulon_name")["Gene"].apply(set)
        gB = dB.groupby("eRegulon_name")["Gene"].apply(set)
        rA = dA.groupby("eRegulon_name")["Region"].apply(set)
        rB = dB.groupby("eRegulon_name")["Region"].apply(set)
        shared_reg = sorted(set(gA.index) & set(gB.index))
        for r in shared_reg:
            reg_rows.append({
                "stage": stage, "arm_A": A, "arm_B": B, "regulon": r,
                "TF": r.split("_")[0],
                "sign_mode": r.split("_")[-1],
                "n_genes_A": len(gA[r]), "n_genes_B": len(gB[r]),
                "jaccard_genes": jac(gA[r], gB[r]),
                "n_regions_A": len(rA[r]), "n_regions_B": len(rB[r]),
                "jaccard_regions": jac(rA[r], rB[r]),
            })

        tgA = dA.groupby("TF")["Gene"].apply(set)
        tgB = dB.groupby("TF")["Gene"].apply(set)
        trA = dA.groupby("TF")["Region"].apply(set)
        trB = dB.groupby("TF")["Region"].apply(set)
        for t in sorted(set(tgA.index) & set(tgB.index)):
            tf_rows.append({
                "stage": stage, "arm_A": A, "arm_B": B, "TF": t,
                "n_regulons_A": int((dA.TF == t).groupby(dA.eRegulon_name).any().sum()),
                "n_genes_A": len(tgA[t]), "n_genes_B": len(tgB[t]),
                "jaccard_genes_union": jac(tgA[t], tgB[t]),
                "jaccard_regions_union": jac(trA[t], trB[t]),
            })

REG = pd.DataFrame(reg_rows)
TFJ = pd.DataFrame(tf_rows)
REG.to_csv(OUT / "f1_shared_regulon_jaccard.csv", index=False)
TFJ.to_csv(OUT / "f1_shared_tf_jaccard.csv", index=False)

for stage in STAGES:
    for A, B in PAIRS:
        k = f"{stage}__{A}_vs_{B}"
        r = REG[(REG.stage == stage) & (REG.arm_A == A) & (REG.arm_B == B)]
        t = TFJ[(TFJ.stage == stage) & (TFJ.arm_A == A) & (TFJ.arm_B == B)]
        res["S1_a_targets"][k] = {
            "L3_regulon_genes": summarise(r.jaccard_genes.values),
            "L3_regulon_regions": summarise(r.jaccard_regions.values),
            "TF_level_gene_union": summarise(t.jaccard_genes_union.values),
            "TF_level_region_union": summarise(t.jaccard_regions_union.values),
        }

# ---------------------------------------------------------------- cross-check
rj = json.load(open(ROOT / "analysis/network_comparison" / "out" / "results.json"))
conv = {c["stage"]: c for c in rj["convergence"]}
struct = {(s["stage"], s["pairing"]): s for s in rj["structure"]}
cc = {}
for stage in STAGES:
    mine = res["S1_0_sets"][f"{stage}__link_vs_pseudomultiome"]
    theirs = conv[stage]
    mine_j = res["S1_a_targets"][f"{stage}__link_vs_pseudomultiome"]["L3_regulon_genes"]["median"]
    cc[stage] = {
        "n_regulons_link": [mine["regulon"]["n_A"], struct[(stage, "link")]["regulons"]],
        "n_regulons_pm": [mine["regulon"]["n_B"], struct[(stage, "pseudomultiome")]["regulons"]],
        "n_TFs_link": [mine["TF"]["n_A"], struct[(stage, "link")]["tfs"]],
        "n_TFs_pm": [mine["TF"]["n_B"], struct[(stage, "pseudomultiome")]["tfs"]],
        "shared_regulons": [mine["regulon"]["n_shared"], theirs["shared_regulons"]],
        "shared_TFs": [mine["TF"]["n_shared"], theirs["tf_shared"]],
        "median_target_jaccard": [mine_j, theirs["median_target_jaccard"]],
    }
    cc[stage]["all_match"] = all(
        (abs(a - b) < 1e-9) if isinstance(a, float) else (a == b)
        for a, b in cc[stage].values() if isinstance(a, (int, float))
    )
res["crosscheck_vs_results_json"] = cc

with open(OUT / "f1_sets_and_targets.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote f1_sets_and_targets.json, f1_set_members.csv, f1_shared_regulon_jaccard.csv, "
    "f1_shared_tf_jaccard.csv")
