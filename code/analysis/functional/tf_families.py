"""FIG1 / S1.b - regulator family composition, link vs pseudo-multiome (and the corrected arm)."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
REF = ROOT / "analysis" / "ref"
OUT = ROOT / "analysis" / "out"

STAGES = ["Diagnosis", "Relapse"]
ARMS = ["link", "pseudomultiome", "repair"]
MIN_FAMILY_TFS = 5   # pre-specified
Q_BAR = 0.05         # pre-specified

IEG = {"FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
       "ATF3", "KLF2", "KLF4", "KLF6", "NR4A1", "IER2", "ZFP36", "DUSP1", "NFIL3"}

AP1_CORE = {"FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"}

def log(m):
    print(f"[f1-family] {m}", flush=True)

def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

# ---------------------------------------------------------------- LEVEL A
lam = pd.read_csv(REF / "humantfs_v1.01_DatabaseExtract.csv", low_memory=False)
lam.columns = [c.strip() for c in lam.columns]
lam_tf = lam[lam["Is TF?"] == "Yes"][["HGNC symbol", "DBD"]].dropna(subset=["HGNC symbol"])
DBD = dict(zip(lam_tf["HGNC symbol"].astype(str), lam_tf["DBD"].astype(str)))
# Lambert also carries records it judges not to be TFs; keep them separately so a SCENIC+ regulator
# that Lambert calls a co-factor is labelled as such rather than as "absent".
lam_no = lam[lam["Is TF?"] == "No"][["HGNC symbol"]].dropna()
NOT_TF = set(lam_no["HGNC symbol"].astype(str))
log(f"Lambert: {len(DBD)} TFs with a DBD; {len(NOT_TF)} records marked Is TF? = No")

def level_b(sym, dbd):
    """Refine only bZIP and C2H2 ZF, by symbol regex, within the Lambert class."""
    if dbd == "bZIP":
        if sym in AP1_CORE:
            return "AP-1"
        if re.fullmatch(r"CEBP[A-Z]", sym):
            return "C/EBP"
        if re.match(r"^(MAF|NFE2|BACH)", sym):
            return "MAF/NFE2/BACH"
        if re.match(r"^(ATF|CREB)", sym) or sym == "XBP1":
            return "ATF/CREB/XBP"
        return "bZIP-other"
    if dbd == "C2H2 ZF":
        if re.fullmatch(r"KLF[0-9]+", sym) or re.fullmatch(r"SP[0-9]+", sym):
            return "KLF/SP"
        return "C2H2-other"
    return dbd

def collapse_b_to_a(fam_b, dbd):
    """Assert Level B is a strict refinement of Level A."""
    if fam_b in ("AP-1", "C/EBP", "MAF/NFE2/BACH", "ATF/CREB/XBP", "bZIP-other"):
        return "bZIP"
    if fam_b in ("KLF/SP", "C2H2-other"):
        return "C2H2 ZF"
    return fam_b

# ---------------------------------------------------------------- load arms
tabs = {}
for stage in STAGES:
    for arm in ARMS:
        d = pd.read_csv(S2 / f"eRegulons_{stage}_{arm}.tsv", sep="\t",
                        usecols=["TF", "eRegulon_name"], low_memory=False)
        tabs[(stage, arm)] = d.drop_duplicates()

all_tfs = sorted({t for k in tabs for t in tabs[k].TF.unique()})
rows = []
for t in all_tfs:
    dbd = DBD.get(t)
    if dbd is None:
        a = "not-in-Lambert (co-factor/chromatin per Lambert)" if t in NOT_TF else "not-in-Lambert"
        b = a
    else:
        a, b = dbd, level_b(t, dbd)
        assert collapse_b_to_a(b, dbd) == dbd, f"Level B broke the refinement for {t}"
    rows.append({"TF": t, "level_A_Lambert_DBD": a, "level_B_family": b,
                 "level_C_project_IEG": t in IEG})
ASSIGN = pd.DataFrame(rows)
ASSIGN.to_csv(OUT / "f1_family_assignment.csv", index=False)
log(f"assigned {len(ASSIGN)} distinct TFs across all arms/timepoints; "
    f"in Lambert-with-DBD: {(ASSIGN.level_A_Lambert_DBD.str.startswith('not-in') == False).sum()}")

AMAP = dict(zip(ASSIGN.TF, ASSIGN.level_A_Lambert_DBD))
BMAP = dict(zip(ASSIGN.TF, ASSIGN.level_B_family))

# ---------------------------------------------------------------- counts and tests
count_rows, test_rows = [], []
res = {
    "status": "DESCRIPTIVE / EXPLORATORY",
    "prespec": "analysis/FIG1_FIG6_PRESPEC.md S1.b",
    "bars": {"min_family_TFs_to_test": MIN_FAMILY_TFS, "q": Q_BAR},
    "level_C_definition_verbatim": sorted(IEG),
    "level_C_source": "longitudinal/healthy_published/scripts/shared_module.py:33-34",
    "lambert_coverage": {},
}

for stage in STAGES:
    for arm in ARMS:
        d = tabs[(stage, arm)]
        tfs = sorted(d.TF.unique())
        cov = sum(1 for t in tfs if not AMAP[t].startswith("not-in"))
        res["lambert_coverage"][f"{stage}_{arm}"] = {
            "n_TFs": len(tfs), "in_Lambert_with_DBD": cov,
            "frac": cov / len(tfs) if tfs else None,
            "absent": sorted(t for t in tfs if AMAP[t].startswith("not-in")),
        }
        for level, mapper in (("A_Lambert_DBD", AMAP), ("B_refined", BMAP)):
            for unit in ("TF", "regulon"):
                if unit == "TF":
                    s = pd.Series([mapper[t] for t in tfs])
                else:
                    s = d.eRegulon_name.drop_duplicates().map(lambda r: mapper[r.split("_")[0]])
                for fam, n in s.value_counts().items():
                    count_rows.append({"stage": stage, "arm": arm, "level": level,
                                       "unit": unit, "family": fam, "n": int(n),
                                       "total": int(s.size), "frac": float(n / s.size)})
        # level C
        for unit in ("TF", "regulon"):
            src = tfs if unit == "TF" else [r.split("_")[0] for r in d.eRegulon_name.drop_duplicates()]
            n = sum(1 for t in src if t in IEG)
            count_rows.append({"stage": stage, "arm": arm, "level": "C_project_IEG",
                               "unit": unit, "family": "AP-1/IEG (project list)", "n": int(n),
                               "total": len(src), "frac": n / len(src) if src else None})

COUNTS = pd.DataFrame(count_rows)
COUNTS.to_csv(OUT / "f1_family_counts.csv", index=False)

for stage in STAGES:
    for A, B in [("link", "pseudomultiome"), ("link", "repair"), ("pseudomultiome", "repair")]:
        for level in ("A_Lambert_DBD", "B_refined", "C_project_IEG"):
            for unit in ("TF", "regulon"):
                dA, dB = tabs[(stage, A)], tabs[(stage, B)]
                if unit == "TF":
                    mA = sorted(dA.TF.unique()); mB = sorted(dB.TF.unique())
                else:
                    mA = list(dA.eRegulon_name.drop_duplicates())
                    mB = list(dB.eRegulon_name.drop_duplicates())
                if level == "C_project_IEG":
                    famA = ["AP-1/IEG (project list)" if (m if unit == "TF" else m.split("_")[0]) in IEG
                            else "other" for m in mA]
                    famB = ["AP-1/IEG (project list)" if (m if unit == "TF" else m.split("_")[0]) in IEG
                            else "other" for m in mB]
                    fams = ["AP-1/IEG (project list)"]
                else:
                    mp = AMAP if level == "A_Lambert_DBD" else BMAP
                    famA = [mp[m if unit == "TF" else m.split("_")[0]] for m in mA]
                    famB = [mp[m if unit == "TF" else m.split("_")[0]] for m in mB]
                    fams = sorted(set(famA) | set(famB))
                # the family-size filter is always evaluated on TFs in the two-arm union (pre-specified)
                tfU = set(dA.TF.unique()) | set(dB.TF.unique())
                sizeU = {}
                for f in fams:
                    if level == "C_project_IEG":
                        sizeU[f] = len(tfU & IEG)
                    else:
                        mp = AMAP if level == "A_Lambert_DBD" else BMAP
                        sizeU[f] = sum(1 for t in tfU if mp[t] == f)
                pend = []
                for f in fams:
                    a1 = sum(1 for x in famA if x == f); a0 = len(famA) - a1
                    b1 = sum(1 for x in famB if x == f); b0 = len(famB) - b1
                    tested = sizeU[f] >= MIN_FAMILY_TFS
                    if tested:
                        odds, p = stats.fisher_exact([[a1, a0], [b1, b0]])
                    else:
                        odds, p = np.nan, np.nan
                    pend.append({"stage": stage, "arm_A": A, "arm_B": B, "level": level,
                                 "unit": unit, "family": f,
                                 "n_in_A": a1, "n_A_total": len(famA),
                                 "n_in_B": b1, "n_B_total": len(famB),
                                 "frac_A": a1 / len(famA) if famA else None,
                                 "frac_B": b1 / len(famB) if famB else None,
                                 "family_size_TF_union": sizeU[f],
                                 "tested": tested,
                                 "odds_ratio_A_over_B": float(odds) if tested else None,
                                 "fisher_p": float(p) if tested else None})
                ps = [r["fisher_p"] for r in pend if r["tested"]]
                qs = bh(ps) if ps else []
                it = iter(qs)
                for r in pend:
                    r["bh_q"] = float(next(it)) if r["tested"] else None
                    r["significant_q<0.05"] = bool(r["bh_q"] is not None and r["bh_q"] < Q_BAR)
                    r["note"] = None if r["tested"] else f"not tested (family TF union n<{MIN_FAMILY_TFS})"
                test_rows += pend

TESTS = pd.DataFrame(test_rows)
TESTS.to_csv(OUT / "f1_family_tests.csv", index=False)

prim = TESTS[(TESTS.arm_A == "link") & (TESTS.arm_B == "pseudomultiome")
             & (TESTS.level == "B_refined") & (TESTS.unit == "TF") & TESTS.tested]
res["primary_link_vs_pm_levelB_TF"] = json.loads(
    prim.sort_values("fisher_p").to_json(orient="records"))
res["n_significant_primary"] = int(prim["significant_q<0.05"].sum())

with open(OUT / "f1_family.json", "w") as f:
    json.dump(res, f, indent=1)
log(f"primary table: {len(prim)} families tested, {res['n_significant_primary']} at q<{Q_BAR}")
log("wrote f1_family.json, f1_family_assignment.csv, f1_family_counts.csv, f1_family_tests.csv")
