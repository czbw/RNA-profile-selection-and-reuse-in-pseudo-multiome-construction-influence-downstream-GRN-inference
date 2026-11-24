#!/usr/bin/env python
"""Compare the network from the corrected pairing with the published-pairing network.

Reports the recovery of the frozen 26-regulon / 14-TF module against pre-specified thresholds
(pass 0.50, qualified 0.25) and whole-network overlap between the two pairings.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
S2 = ROOT / "paml/malignant" / "out"
OUT = ROOT / "paml/corrected" / "out"

MODDEF = ROOT / "adult/published" / "out" / "module_definition.json"
_md = json.loads(MODDEF.read_text(encoding="utf-8"))
T4_REGULONS = set(_md["regulons"])
T4_TFS = list(_md["TFs"])
T4_BACKGROUND = set(_md.get("background_regulons", []))
assert len(T4_REGULONS) == _md["n_regulons"] == 26, "frozen T4 is not 26 regulons"
assert len(T4_TFS) == _md["n_TFs"] == 14, "frozen T4 is not 14 TFs"
PASS_T, QUAL_T = 0.50, 0.25
CAL = {"pbmc_link_1to1_vs_truth_L4": 0.670, "pbmc_link_vs_truth_L4": 0.337}

def log(m): print(f"[a03] {m}", flush=True)

def jac(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else float("nan")

def load(stage: str, arm: str):
    p = S2 / f"eRegulons_{stage}_{arm}.tsv"
    if not p.exists():
        raise SystemExit(f"missing {p}")
    df = pd.read_csv(p, sep="\t", usecols=["Gene", "TF", "eRegulon_name"])
    tgt = df.groupby("eRegulon_name")["Gene"].apply(lambda s: set(s)).to_dict()
    return {"regulons": set(df["eRegulon_name"].unique()),
            "tfs": set(df["TF"].unique()),
            "targets": tgt,
            "tf_of": df.drop_duplicates("eRegulon_name").set_index("eRegulon_name")["TF"].to_dict(),
            "n_rows": len(df)}

def main():
    rep = {"framing": "de novo primary; overlap reported as construction dependence",
           "prereg": {"pass": PASS_T, "qualified": QUAL_T, "calibration": CAL},
           "t4_tfs": T4_TFS, "stages": {}}

    for stage in ("Diagnosis", "Relapse"):
        log(f"===== {stage} =====")
        old = load(stage, "link")
        new = load(stage, "repair")
        log(f"  published (link)  {len(old['regulons']):5,} regulons | {len(old['tfs']):4} TFs "
            f"| {old['n_rows']:9,} rows")
        log(f"  corrected (repair){len(new['regulons']):5,} regulons | {len(new['tfs']):4} TFs "
            f"| {new['n_rows']:9,} rows")

        # ---- frozen module, from its frozen definition -------------------------------------------
        t4 = set(T4_REGULONS)
        t4_in_old = t4 & old["regulons"]
        log(f"  T4 (frozen, {MODDEF.name}): {len(t4)} regulons; "
            f"{len(t4_in_old)} present in the published network for this stage")
        t4_exact = True

        t4_recovered = t4 & new["regulons"]
        t4_jac = jac(t4, t4 | (new["regulons"] & t4))  # recovery of the frozen set
        t4_frac = len(t4_recovered) / len(t4) if t4 else float("nan")
        verdict = ("PASS" if t4_frac >= PASS_T else
                   "QUALIFIED" if t4_frac >= QUAL_T else "FAIL")

        # ---- whole-network agreement, old vs new -------------------------------------
        shared = old["regulons"] & new["regulons"]
        tgt_j = [jac(old["targets"][r], new["targets"][r]) for r in sorted(shared)]
        tgt_med = float(pd.Series(tgt_j).median()) if tgt_j else float("nan")

        s = {
            "published_n_regulons": len(old["regulons"]), "published_n_tfs": len(old["tfs"]),
            "corrected_n_regulons": len(new["regulons"]), "corrected_n_tfs": len(new["tfs"]),
            "t4_reconstruction_exact": t4_exact, "t4_n_regulons_reconstructed": len(t4),
            "t4_regulons_recovered": len(t4_recovered),
            "t4_fraction_recovered": t4_frac,
            "t4_jaccard_within_t4_tfs": t4_jac,
            "t4_tfs_recovered": sorted(set(T4_TFS) & new["tfs"]),
            "t4_tfs_lost": sorted(set(T4_TFS) - new["tfs"]),
            "prereg_verdict": verdict,
            "network_regulon_jaccard": jac(old["regulons"], new["regulons"]),
            "network_tf_jaccard": jac(old["tfs"], new["tfs"]),
            "n_shared_regulons": len(shared),
            "median_target_set_jaccard_shared": tgt_med,
            "corrected_only_regulons": len(new["regulons"] - old["regulons"]),
            "published_only_regulons": len(old["regulons"] - new["regulons"]),
        }
        rep["stages"][stage] = s

        log(f"  T4: {len(t4_recovered)}/{len(t4)} regulons recovered "
            f"= {t4_frac:.3f}  ->  PRE-SPECIFIED VERDICT: {verdict}")
        log(f"  T4 TFs recovered {len(s['t4_tfs_recovered'])}/14"
            + (f"; LOST: {', '.join(s['t4_tfs_lost'])}" if s["t4_tfs_lost"] else "; none lost"))
        log(f"  whole network: regulon Jaccard {s['network_regulon_jaccard']:.3f} | "
            f"TF Jaccard {s['network_tf_jaccard']:.3f} | "
            f"median target-set Jaccard over {len(shared)} shared {tgt_med:.3f}")
        log(f"  de novo: {s['corrected_only_regulons']:,} regulons appear ONLY under the "
            f"corrected pairing; {s['published_only_regulons']:,} only under the published one")

    (OUT / "a03_t2_report.json").write_text(json.dumps(rep, indent=1), encoding="utf-8")
    log(f"wrote {OUT / 'a03_t2_report.json'}")

if __name__ == "__main__":
    main()
