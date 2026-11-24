"""The evidence-stacked nomination funnel."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/corrected')
REPO = ROOT.parents[1]
OUT = ROOT / "out"
S3B = REPO / "longitudinal/healthy_corrected" / "out"
GT = REPO / "data/reference_support"
S2A = REPO / "analysis/network_comparison" / "out"

CT_BAR = 4          # minimum cell types (of 6) in which the direction must hold
tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[s3-nom] {m}", flush=True)

# Assemble the axes

did = pd.read_csv(S3B / "did_reconciliation.csv", index_col=0)
dec = pd.read_csv(OUT / "decomposition_reconciliation.csv", index_col=0)
e2 = pd.read_csv(OUT / "e2_paired_wilcoxon_reconciliation.csv", index_col=0)
cat = pd.read_csv(S2A / "regulon_catalogue.csv")
tier = {r: str(t).split()[0] for r, t in zip(cat.regulon, cat.tier)}

# The ChIP-seq axis needs a K562 reference for the repair pairing; skip it if absent.
AX5_AVAILABLE = (GT / "gate_malig_Diagnosis_repair_K562_w5000.csv").exists()
if AX5_AVAILABLE:
    k = {s: pd.read_csv(GT / f"gate_malig_{s}_repair_K562_w5000.csv").set_index("tf")
         for s in ("Diagnosis", "Relapse")}
else:
    print("[s3-nom] ax5 UNAVAILABLE on the repair arm (no repair-arm K562 ChIP-seq reference; the "
          "axis is not used). ax5=False for all; read survives_4, not survives_all5.", flush=True)
    k = None

d = did.join(dec[["n_celltypes_same_sign_Diagnosis", "n_celltypes_same_sign_Relapse"]])
d["TF"] = [tf_of(i) for i in d.index]
d["tier"] = [tier.get(i, "none") for i in d.index]
d["beta_mal_Dx"] = d.beta_malignant_Diagnosis
d["beta_mal_Rel"] = d.beta_malignant_Relapse

d["ax1_directional"] = [i in set(e2.index[e2.sign_agree]) for i in d.index]
d["ax2_selective"] = ~d.compartments_agree
d["ax3_celltype_robust"] = ((d.n_celltypes_same_sign_Diagnosis >= CT_BAR) &
                            (d.n_celltypes_same_sign_Relapse >= CT_BAR)).fillna(False)
d["ax4_tier1"] = d.tier == "T1"
d["margin_K562_Dx"] = ([k["Diagnosis"].margin.get(t, np.nan) for t in d.TF]
                       if k is not None else np.nan)
d["margin_K562_Rel"] = ([k["Relapse"].margin.get(t, np.nan) for t in d.TF]
                        if k is not None else np.nan)
d["ax5_chip"] = (((d.margin_K562_Dx > 0) & (d.margin_K562_Rel > 0)).fillna(False)
                 if k is not None else False)
d["n_axes"] = d[[f"ax{i}_" + n for i, n in enumerate(
    ["directional", "selective", "celltype_robust", "tier1", "chip"], 1)]].sum(axis=1)
d["direction"] = np.where(d.did_min_abs > 0, "up at relapse", "down at relapse")

# Regulon activity direction is not TF output direction. A leading `-` in the sign mode means the
# TF's presence is associated with the region closing, so the targets are genes the TF suppresses:
# rising activity of a `-`-mode regulon implies less suppression, i.e. falling TF output.
lead = pd.Series([(re.search(r"_([+-])/[+-]$", str(i)).group(1)
                   if re.search(r"_([+-])/[+-]$", str(i)) else "?") for i in d.index], index=d.index)
act_up = np.sign(d.beta_mal_Dx) > 0
d["sign_mode_lead"] = lead
d["tf_output"] = np.where(lead == "+", np.where(act_up, "up", "down"),
                   np.where(lead == "-", np.where(act_up, "down", "up"), "?"))
d["rank_stat"] = np.minimum(d.beta_mal_Dx.abs(), d.beta_mal_Rel.abs())      # conservative effect size

# The funnel; every depth is reported

steps = [("all shared-name regulons", pd.Series(True, index=d.index)),
         ("1. directional under both vocabularies", d.ax1_directional),
         ("2. + selective vs the patient's own healthy cells", d.ax1_directional & d.ax2_selective),
         ("3. + direction holds in >=4 of 6 cell types",
          d.ax1_directional & d.ax2_selective & d.ax3_celltype_robust),
         ("4. + tier 1 (both pairings, both stages)",
          d.ax1_directional & d.ax2_selective & d.ax3_celltype_robust & d.ax4_tier1),
         ("5. + TF ChIP-seq supported at both stages",
          d.ax1_directional & d.ax2_selective & d.ax3_celltype_robust & d.ax4_tier1 & d.ax5_chip)]

res = {"status": "EXPLORATORY hypothesis-ranking — NOT a significance test, NOT FDR-controlled",
       "criteria": {"cell_type_bar_of_6": CT_BAR,
                    "chip_reference": "ENCODE K562, link pairing, TSS+/-5kb, margin>0 both stages",
                    "rank_statistic": "min(|beta_Dx|, |beta_Rel|) — the conservative reading"},
       "funnel": []}
for name, m in steps:
    res["funnel"].append({"step": name, "n_regulons": int(m.sum()),
                          "n_TFs": int(d.loc[m, "TF"].nunique())})
    log(f"{name}: {int(m.sum())} regulons / {int(d.loc[m,'TF'].nunique())} TFs")

d["survives_all5"] = steps[-1][1]
d["survives_4"] = steps[-2][1]
d.sort_values(["n_axes", "rank_stat"], ascending=[False, False]).to_csv(
    OUT / "nomination_funnel.csv")

# also report the alternative depths a reader might prefer
alt = {"drop_chip_axis (tier-1 core)": steps[-2][1],
       "drop_tier_axis (chip-supported, any tier)":
           d.ax1_directional & d.ax2_selective & d.ax3_celltype_robust & d.ax5_chip}
for name, m in alt.items():
    res.setdefault("alternative_depths", {})[name] = {
        "n_regulons": int(m.sum()), "n_TFs": int(d.loc[m, "TF"].nunique()),
        "TFs": sorted(d.loc[m, "TF"].unique().tolist())}
    log(f"[alt] {name}: {int(m.sum())} regulons / {int(d.loc[m,'TF'].nunique())} TFs")

top = d[d.survives_all5].sort_values("rank_stat", ascending=False)
res["all_five_axes"] = {
    "n_regulons": int(len(top)), "n_TFs": int(top.TF.nunique()),
    "up_at_relapse": int((top.direction == "up at relapse").sum()),
    "down_at_relapse": int((top.direction == "down at relapse").sum()),
    "regulons": [{"regulon": i, "TF": r.TF, "direction": r.direction,
                  "beta_mal_Dx": float(r.beta_mal_Dx), "beta_mal_Rel": float(r.beta_mal_Rel),
                  "beta_healthy_Dx": float(r.beta_healthy_Diagnosis),
                  "did_min_abs": float(abs(r.did_min_abs)),
                  "celltypes_same_sign": [int(r.n_celltypes_same_sign_Diagnosis),
                                          int(r.n_celltypes_same_sign_Relapse)],
                  "K562_margin": [float(r.margin_K562_Dx), float(r.margin_K562_Rel)]}
                 for i, r in top.iterrows()]}
# per-TF consistency: do all of a TF's surviving regulons imply the SAME output direction?
tfsum = []
for tf, g in top.groupby("TF"):
    dirs = set(g.tf_output)
    tfsum.append({"TF": tf, "n_regulons": int(len(g)),
                  "tf_output": (dirs.pop() if len(dirs) == 1 else "INCONSISTENT"),
                  "consistent": len(dirs) == 0,
                  "regulons": list(g.index),
                  "K562_margin_Dx": float(g.margin_K562_Dx.iloc[0]),
                  "K562_margin_Rel": float(g.margin_K562_Rel.iloc[0])})
for r in tfsum:
    r["consistent"] = r["tf_output"] != "INCONSISTENT"
res["all_five_axes"]["per_TF"] = tfsum
for r in tfsum:
    log(f"  {r['TF']}: {r['n_regulons']} regulons -> TF output {r['tf_output'].upper()} at relapse "
        f"({'consistent' if r['consistent'] else 'INCONSISTENT — do not read'})")
log(f"ALL FIVE AXES: {len(top)} regulons / {top.TF.nunique()} TFs — "
    f"{', '.join(sorted(top.TF.unique()))}")

with open(OUT / "s3_nomination.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote s3_nomination.json")
