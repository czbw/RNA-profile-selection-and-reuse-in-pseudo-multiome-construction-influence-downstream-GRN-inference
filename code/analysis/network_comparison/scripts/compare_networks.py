"""Full local analysis of the Dx + Rel SCENIC+ outputs."""
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = (Path(__file__).resolve().parents[3] / 'analysis/network_comparison')
S2 = ROOT.parents[1] / "paml/malignant" / "out"
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)

STAGES = ["Diagnosis", "Relapse"]
PAIRINGS = ["link", "pseudomultiome"]
RNG = np.random.default_rng(555)
N_BOOT = 2000

# Pre-specified before any cross-stage number was computed (module D).
DRIFT_JACCARD_MIN = 0.50
MIN_CELLS_PER_UNIT = 20

regkey = lambda s: re.sub(r"_\(\d+[gr]\)$", "", str(s))
tf_of = lambda s: re.split(r"[_(]", str(s))[0]

def log(m):
    print(f"[s2] {m}", flush=True)

# I/O

def read_auc(path):
    """Gene_based modality of an AUCell .h5mu -> (X float64, obs_names, var_names)."""
    with h5py.File(path, "r") as f:
        g = f["mod"]["Gene_based"]
        X = np.asarray(g["X"][:], dtype=np.float64)
        oi, vi = g["obs"].attrs["_index"], g["var"].attrs["_index"]
        dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
        return X, np.array(dec(g["obs"][oi][:])), np.array(dec(g["var"][vi][:]))

def collapse_dup_cols(df):
    """Average columns that share a regulon key (regkey can collide in principle)."""
    return df.T.groupby(level=0).mean().T if df.columns.duplicated().any() else df

def z(M):
    m = M - M.mean(0, keepdims=True)
    s = M.std(0, keepdims=True)
    s[s == 0] = np.inf
    return m / s

results = {"meta": {"seed": 555, "n_boot": N_BOOT,
                    "drift_jaccard_min": DRIFT_JACCARD_MIN,
                    "min_cells_per_unit": MIN_CELLS_PER_UNIT}}

# Load cell manifests

cell_data = {}
for s in STAGES:
    d = pd.read_csv(S2 / f"cell_data_{s}.csv", index_col=0, dtype=str)
    cell_data[s] = d
    log(f"cell_data {s}: {len(d):,} cells | {d.Classified_Celltype.nunique()} celltypes | "
        f"{d.Patient_ID.nunique()} patients")

# Pass 1 — activity profiles (the only step that touches the big matrices)

prof_ct = {}        # (stage,pairing) -> DataFrame celltype x regulon
prof_pat = {}       # (stage,'link')  -> DataFrame (patient,celltype) x regulon
auc_audit = []

for s in STAGES:
    cd = cell_data[s]
    ct_by_cell = cd.Classified_Celltype
    pt_by_cell = cd.Patient_ID
    for p in PAIRINGS:
        X, obs, var = read_auc(S2 / f"eRegulon_AUC_{s}_{p}.h5mu")
        cols = [regkey(v) for v in var]
        df = pd.DataFrame(X, columns=cols)
        del X
        if p == "link":
            ct = ct_by_cell.reindex(obs).values
            pat = pt_by_cell.reindex(obs).values
        else:
            ct = np.array([re.sub(r"_\d+$", "", n) for n in obs])
            pat = np.full(len(obs), None)
        keep = pd.notna(ct)
        auc_audit.append(dict(stage=s, pairing=p, units=len(obs), mapped=int(keep.sum()),
                              regulons=len(set(cols)),
                              celltypes=int(pd.Series(ct[keep]).nunique())))
        sub = df.loc[keep]
        prof_ct[(s, p)] = collapse_dup_cols(sub.groupby(ct[keep]).mean())
        if p == "link":
            idx = pd.MultiIndex.from_arrays([pat[keep], ct[keep]], names=["patient", "celltype"])
            g = sub.set_index(idx).groupby(level=["patient", "celltype"])
            prof_pat[s] = collapse_dup_cols(g.mean())
            prof_pat[s + "_n"] = g.size()
        del df, sub
        log(f"AUC {s}/{p}: {len(obs):,} units, {int(keep.sum()):,} mapped, "
            f"{len(set(cols))} regulons, {int(pd.Series(ct[keep]).nunique())} celltypes")

results["auc_audit"] = auc_audit

# MODULE B — structure (edge tables)

SIGN_RE = re.compile(r"_(direct|extended)_([+-])/([+-])")
edges, targets, regions, struct = {}, {}, {}, []

for s in STAGES:
    for p in PAIRINGS:
        d = pd.read_csv(S2 / f"eRegulons_{s}_{p}.tsv", sep="\t")
        d["regulon"] = d.eRegulon_name.map(regkey)
        edges[(s, p)] = set(map(tuple, d[["TF", "Gene"]].drop_duplicates().values))
        targets[(s, p)] = {r: set(g.Gene.astype(str)) for r, g in d.groupby("regulon")}
        regions[(s, p)] = {r: set(g.Region.astype(str)) for r, g in d.groupby("regulon")}
        mode = d.eRegulon_name.str.extract(SIGN_RE)
        sizes = pd.Series({r: len(v) for r, v in targets[(s, p)].items()})
        struct.append(dict(
            stage=s, pairing=p, rows=len(d),
            tfs=int(d.TF.nunique()), regulons=int(d.regulon.nunique()),
            edges=len(edges[(s, p)]), genes=int(d.Gene.nunique()),
            regions_n=int(d.Region.nunique()),
            direct=int((mode[0] == "direct").sum()), extended=int((mode[0] == "extended").sum()),
            activating=int(d.eRegulon_name.str.contains(r"_\+/\+", regex=True).sum()),
            repressing=int(d.eRegulon_name.str.contains(r"_-/-", regex=True).sum()),
            median_targets=float(sizes.median()), max_targets=int(sizes.max()),
            median_regions=float(np.median([len(v) for v in regions[(s, p)].values()]))))
        log(f"struct {s}/{p}: {len(d):,} rows, {d.TF.nunique()} TFs, "
            f"{d.regulon.nunique()} regulons, {len(edges[(s,p)]):,} unique TF-gene edges")
        del d

results["structure"] = struct

# MODULE A — pairing convergence, with mismatched-regulon and identifiability controls

conv, per_reg = [], {}

for s in STAGES:
    PA, PB = prof_ct[(s, "link")], prof_ct[(s, "pseudomultiome")]
    cts = sorted(set(PA.index) & set(PB.index))
    regs = sorted(set(PA.columns) & set(PB.columns))
    tA = {tf_of(r) for r in PA.columns}
    tB = {tf_of(r) for r in PB.columns}
    eA, eB = edges[(s, "link")], edges[(s, "pseudomultiome")]

    A, B = z(PA.loc[cts, regs].values), z(PB.loc[cts, regs].values)
    C = (A.T @ B) / len(cts)
    matched = np.diag(C).copy()
    off = C[~np.eye(len(regs), dtype=bool)]

    bm = np.array([np.median(matched[RNG.integers(0, len(regs), len(regs))]) for _ in range(N_BOOT)])
    pct = np.array([(C[i] < C[i, i]).mean() for i in range(len(regs))])
    rank = np.array([(C[i] > C[i, i]).sum() + 1 for i in range(len(regs))])
    bp = np.array([np.median(pct[RNG.integers(0, len(regs), len(regs))]) for _ in range(N_BOOT)])

    # target-set agreement per shared regulon (level 3, per regulon rather than pooled)
    tj = []
    for r in regs:
        a, b = targets[(s, "link")].get(r, set()), targets[(s, "pseudomultiome")].get(r, set())
        tj.append(len(a & b) / len(a | b) if (a | b) else np.nan)

    conv.append(dict(
        stage=s, cell_types=len(cts), shared_regulons=len(regs),
        regulons_A=len(set(PA.columns)), regulons_B=len(set(PB.columns)),
        tf_A=len(tA), tf_B=len(tB), tf_shared=len(tA & tB),
        edge_jaccard=len(eA & eB) / max(len(eA | eB), 1),
        act_median_r=float(np.median(matched)),
        act_lo=float(np.percentile(bm, 2.5)), act_hi=float(np.percentile(bm, 97.5)),
        null_median_r=float(np.median(off)), null_p95=float(np.percentile(off, 95)),
        ident_top1=float((rank == 1).mean()), ident_pct=float(np.median(pct)),
        ident_pct_lo=float(np.percentile(bp, 2.5)), ident_pct_hi=float(np.percentile(bp, 97.5)),
        frac_r_gt_half=float((matched > 0.5).mean()),
        margin_over_null=float(np.median(matched) - np.median(off)),
        median_target_jaccard=float(np.nanmedian(tj))))

    per_reg[s] = pd.DataFrame({"regulon": regs, "TF": [tf_of(r) for r in regs],
                               "r": matched, "rank": rank, "pctile": pct,
                               "target_jaccard": tj}).sort_values("r", ascending=False)
    per_reg[s].to_csv(OUT / f"convergence_per_regulon_{s}.csv", index=False)
    np.save(OUT / f"null_offdiag_{s}.npy", off[RNG.integers(0, len(off), min(len(off), 200000))])
    np.save(OUT / f"matched_diag_{s}.npy", matched)
    np.save(OUT / f"ident_pct_{s}.npy", pct)
    log(f"CONVERGENCE {s}: r={np.median(matched):+.3f} [{np.percentile(bm,2.5):+.3f},"
        f"{np.percentile(bm,97.5):+.3f}] null={np.median(off):+.3f} "
        f"ident={100*np.median(pct):.1f}% top1={100*(rank==1).mean():.1f}% "
        f"edgeJ={len(eA&eB)/max(len(eA|eB),1):.3f}")

results["convergence"] = conv

# MODULE C — composition + Stage-3 power

comp_rows, drivers = [], []
for s in STAGES:
    cd = cell_data[s]
    for ct, n in cd.Classified_Celltype.value_counts().items():
        comp_rows.append(dict(stage=s, celltype=ct, cells=int(n),
                              frac=float(n / len(cd)),
                              patients=int(cd[cd.Classified_Celltype == ct].Patient_ID.nunique())))
    for dv, n in cd.Driver.value_counts().items():
        drivers.append(dict(stage=s, driver=str(dv), cells=int(n),
                            patients=int(cd[cd.Driver == dv].Patient_ID.nunique())))
comp = pd.DataFrame(comp_rows)
comp.to_csv(OUT / "composition_by_celltype.csv", index=False)
pd.DataFrame(drivers).to_csv(OUT / "driver_distribution.csv", index=False)
results["composition"] = comp_rows
results["drivers"] = drivers

pats = {s: set(cell_data[s].Patient_ID.unique()) for s in STAGES}
paired = sorted(pats["Diagnosis"] & pats["Relapse"])
results["patients"] = {"Diagnosis": len(pats["Diagnosis"]), "Relapse": len(pats["Relapse"]),
                       "paired": len(paired), "paired_ids": paired}
log(f"patients: Dx {len(pats['Diagnosis'])}, Rel {len(pats['Relapse'])}, paired {len(paired)}")

# (patient x celltype x stage) units that clear the floor — Stage-3 power
unit_rows = []
for s in STAGES:
    g = cell_data[s].groupby(["Patient_ID", "Classified_Celltype"]).size()
    for (pt, ct), n in g.items():
        unit_rows.append(dict(stage=s, patient=pt, celltype=ct, cells=int(n)))
units = pd.DataFrame(unit_rows)
units.to_csv(OUT / "units_patient_celltype.csv", index=False)

wide = units.pivot_table(index=["patient", "celltype"], columns="stage", values="cells")
both = wide.dropna()
usable = both[(both.Diagnosis >= MIN_CELLS_PER_UNIT) & (both.Relapse >= MIN_CELLS_PER_UNIT)]
results["power"] = {
    "units_total": int(len(units)),
    "units_both_stages": int(len(both)),
    "units_usable": int(len(usable)),
    "patients_with_usable_unit": int(usable.reset_index().patient.nunique()),
    "celltypes_with_usable_unit": int(usable.reset_index().celltype.nunique()),
    "usable_by_celltype": usable.reset_index().celltype.value_counts().to_dict()}
log(f"power: {len(both)} patient x celltype units at both stages, "
    f"{len(usable)} clear the {MIN_CELLS_PER_UNIT}-cell floor")

# MODULE D — cross-stage target-set drift, then the exploratory contrast

drift_rows = []
for p in PAIRINGS:
    shared = set(targets[("Diagnosis", p)]) & set(targets[("Relapse", p)])
    for r in sorted(shared):
        a, b = targets[("Diagnosis", p)][r], targets[("Relapse", p)][r]
        drift_rows.append(dict(pairing=p, regulon=r, TF=tf_of(r),
                               n_dx=len(a), n_rel=len(b),
                               jaccard=len(a & b) / len(a | b) if (a | b) else np.nan))
drift = pd.DataFrame(drift_rows)
drift.to_csv(OUT / "cross_stage_target_drift.csv", index=False)
results["drift"] = {
    p: {"n_shared_regulons": int((drift.pairing == p).sum()),
        "median_jaccard": float(drift[drift.pairing == p].jaccard.median()),
        "frac_above_threshold": float((drift[drift.pairing == p].jaccard >= DRIFT_JACCARD_MIN).mean()),
        "n_above_threshold": int((drift[drift.pairing == p].jaccard >= DRIFT_JACCARD_MIN).sum())}
    for p in PAIRINGS}
for p in PAIRINGS:
    log(f"drift {p}: {results['drift'][p]['n_shared_regulons']} shared regulon names, "
        f"median target Jaccard {results['drift'][p]['median_jaccard']:.3f}, "
        f"{results['drift'][p]['n_above_threshold']} clear {DRIFT_JACCARD_MIN}")

# Exploratory within-patient Dx->Rel contrast (link pairing only — metacells carry no patient ID)
PD, PR = prof_pat["Diagnosis"], prof_pat["Relapse"]
nD, nR = prof_pat["Diagnosis_n"], prof_pat["Relapse_n"]
eligible = set(drift[(drift.pairing == "link") & (drift.jaccard >= DRIFT_JACCARD_MIN)].regulon)
regs_c = sorted(set(PD.columns) & set(PR.columns) & eligible)

keyD = {(a, b) for a, b in PD.index}
keyR = {(a, b) for a, b in PR.index}
shared_units = sorted(
    u for u in keyD & keyR
    if nD.get(u, 0) >= MIN_CELLS_PER_UNIT and nR.get(u, 0) >= MIN_CELLS_PER_UNIT
    and u[0] in paired)

contrast = []
if regs_c and shared_units:
    su = pd.MultiIndex.from_tuples(shared_units, names=["patient", "celltype"])
    D = PD.loc[su, regs_c]
    R = PR.loc[su, regs_c]
    delta = (R - D)
    # average within patient across the cell types that survive at BOTH stages, so a patient
    # with many cell types does not dominate and composition shifts are partly absorbed
    per_patient = delta.groupby(level="patient").mean()
    for r in regs_c:
        v = per_patient[r].dropna().values
        if len(v) < 5:
            continue
        w = stats.wilcoxon(v, alternative="two-sided", zero_method="wilcox")
        contrast.append(dict(regulon=r, TF=tf_of(r), n_patients=len(v),
                             mean_delta=float(np.mean(v)), median_delta=float(np.median(v)),
                             stat=float(w.statistic), p=float(w.pvalue)))
    cdf = pd.DataFrame(contrast).sort_values("p")
    if len(cdf):
        rank_ = np.arange(1, len(cdf) + 1)
        cdf["q_BH"] = np.minimum.accumulate(
            (cdf.p.values * len(cdf) / rank_)[::-1])[::-1].clip(0, 1)
        cdf.to_csv(OUT / "within_patient_dx_rel_contrast.csv", index=False)
    results["contrast"] = {
        "regulons_tested": int(len(cdf)), "patients": len(per_patient),
        "units": len(shared_units),
        "n_q05": int((cdf.q_BH < 0.05).sum()) if len(cdf) else 0,
        "top_up": cdf.nlargest(10, "mean_delta")[
            ["regulon", "n_patients", "mean_delta", "p", "q_BH"]].to_dict("records") if len(cdf) else [],
        "top_down": cdf.nsmallest(10, "mean_delta")[
            ["regulon", "n_patients", "mean_delta", "p", "q_BH"]].to_dict("records") if len(cdf) else []}
    log(f"contrast: {len(cdf)} regulons over {len(per_patient)} paired patients, "
        f"{int((cdf.q_BH < 0.05).sum())} at q<0.05")
else:
    results["contrast"] = {"regulons_tested": 0, "patients": 0, "units": 0, "n_q05": 0,
                           "top_up": [], "top_down": []}
    log("contrast: NOT RUN (no eligible regulon x unit combination)")

# MODULE E — confidence-tiered catalogue + canonical-TF face validity

sets = {(s, p): set(targets[(s, p)]) for s in STAGES for p in PAIRINGS}
allr = sorted(set().union(*sets.values()))
cat = []
for r in allr:
    inn = {f"{s[:2]}_{p[:4]}": (r in sets[(s, p)]) for s in STAGES for p in PAIRINGS}
    both_dx = inn["Di_link"] and inn["Di_pseu"]
    both_rel = inn["Re_link"] and inn["Re_pseu"]
    tier = ("T1 both pairings, both stages" if both_dx and both_rel else
            "T2 both pairings, one stage" if both_dx or both_rel else
            "T3 single pairing")
    rr = {}
    for s in STAGES:
        m = per_reg[s][per_reg[s].regulon == r]
        rr[s] = float(m.r.iloc[0]) if len(m) else np.nan
    cat.append(dict(regulon=r, TF=tf_of(r), tier=tier, **inn,
                    pairing_r_Diagnosis=rr["Diagnosis"], pairing_r_Relapse=rr["Relapse"]))
catalogue = pd.DataFrame(cat)
catalogue.to_csv(OUT / "regulon_catalogue.csv", index=False)
results["catalogue"] = {"total": int(len(catalogue)),
                        "by_tier": catalogue.tier.value_counts().to_dict(),
                        "tier1_tfs": int(catalogue[catalogue.tier.str.startswith("T1")].TF.nunique())}
log(f"catalogue: {len(catalogue)} regulons | " +
    " | ".join(f"{k}={v}" for k, v in catalogue.tier.value_counts().items()))

CANON = {
    "SPI1": "PU.1 — myeloid master regulator",
    "RUNX1": "core-binding factor; RUNX1-RUNX1T1 is a cohort driver",
    "CBFB": "CBF partner of RUNX1",
    "CEBPA": "granulocytic commitment; mutated in AML",
    "GATA1": "erythroid master regulator",
    "GATA2": "HSC/progenitor maintenance",
    "TAL1": "SCL — HSC/erythroid",
    "KLF1": "erythroid",
    "NFE2": "erythroid/megakaryocytic",
    "LMO2": "HSC complex partner of TAL1",
    "MYB": "progenitor proliferation; AML dependency",
    "MYC": "proliferation",
    "IRF8": "monocyte/DC",
    "ERG": "HSC; poor-prognosis AML",
    "FLI1": "megakaryocytic/HSC",
    "MEF2C": "KMT2A-r AML dependency",
    "HOXA9": "KMT2A-r AML hallmark",
    "MEIS1": "HOXA9 cofactor in KMT2A-r AML",
    "PBX3": "HOX cofactor",
    "GFI1": "granulocytic",
    "GFI1B": "erythroid/megakaryocytic",
    "IKZF1": "lymphoid",
    "BCL11A": "lymphoid/HSC",
    "STAT5A": "cytokine signalling; FLT3-ITD axis",
    "LYL1": "HSC",
}
face = []
for tf, role in CANON.items():
    row = dict(TF=tf, role=role)
    for s in STAGES:
        m = per_reg[s][per_reg[s].TF == tf]
        row[f"{s}_regulons"] = int(len(m))
        row[f"{s}_max_r"] = float(m.r.max()) if len(m) else np.nan
    row["tier1"] = bool(((catalogue.TF == tf) & catalogue.tier.str.startswith("T1")).any())
    face.append(row)
faced = pd.DataFrame(face)
faced.to_csv(OUT / "canonical_tf_face_validity.csv", index=False)
results["face_validity"] = faced.replace({np.nan: None}).to_dict("records")
log(f"face validity: {int(faced.tier1.sum())}/{len(faced)} canonical TFs in tier 1")

with open(OUT / "results.json", "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=2, default=str)
log(f"-> {OUT/'results.json'}")
