#!/usr/bin/env python
'AML reuse recounted with the full RNA-profile identity (library plus barcode).'
import argparse
import csv
import json
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]))
a = ap.parse_args()
ROOT = Path(a.root)
PUB = ROOT / "paml/malignant" / "out"
COR = ROOT / "paml/corrected" / "out"
OUT = ROOT / "analysis" / "out"
STAGES = ("Diagnosis", "Relapse")

def key(d):
    return d["Linked_scRNA_Sample"].astype(str) + "||" + d["Linked_scRNA_barcode"].astype(str)

def counts(d):
    n = len(d)
    comp = int(key(d).nunique())
    bc = int(d["Linked_scRNA_barcode"].nunique())
    return {"n": n, "distinct": comp, "reuse": n / comp, "distinct_barcode_only": bc,
            "reuse_barcode_only": n / bc}

census = list(csv.DictReader(open(COR / "a00_stratum_census.csv", encoding="utf-8")))
cols = ["Linked_scRNA_barcode", "Linked_scRNA_Sample"]
res = {"definition": "distinct RNA profile = Linked_scRNA_Sample plus Linked_scRNA_barcode",
       "stages": {}, "checks": {}}
frames = []
for s in STAGES:
    pub = pd.read_csv(PUB / f"cell_data_{s}.csv", usecols=cols, low_memory=False)
    cor = pd.read_csv(COR / f"cell_data_{s}.repaired.csv", usecols=cols, low_memory=False)
    frames.append(pub)
    p, c = counts(pub), counts(cor)
    rows = [r for r in census if r["timepoint"] == s]
    una = sum(int(r["n_atac"]) for r in rows if str(r.get("unassignable", "")).strip().lower() == "true")
    atac = sum(int(r["n_atac"]) for r in rows)
    smin = sum(min(int(r["n_atac"]), int(r["n_rna"])) for r in rows)
    floor = (atac - una) / smin
    res["stages"][s] = {
        "old_reuse": p["reuse"], "old_distinct_donors": p["distinct"], "n_in": p["n"],
        "achieved_reuse": c["reuse"], "distinct_donors": c["distinct"], "n_out": c["n"],
        "old_reuse_barcode_only": p["reuse_barcode_only"],
        "old_distinct_donors_barcode_only": p["distinct_barcode_only"],
        "achieved_reuse_barcode_only": c["reuse_barcode_only"],
        "distinct_donors_barcode_only": c["distinct_barcode_only"],
        "scarcity_floor": floor}
    res["checks"][f"{s}_corrected_reuse_equals_floor"] = bool(abs(c["reuse"] - floor) < 1e-9)
    assert res["checks"][f"{s}_corrected_reuse_equals_floor"], (s, c["reuse"], floor)
pooled = counts(pd.concat(frames))
res["published_pooled"] = {"lambo_link_reuse": pooled["reuse"], "lambo_link_distinct_donors": pooled["distinct"],
                           "n": pooled["n"], "lambo_link_reuse_barcode_only": pooled["reuse_barcode_only"],
                           "lambo_link_distinct_donors_barcode_only": pooled["distinct_barcode_only"]}
OUT.mkdir(parents=True, exist_ok=True)
json.dump(res, open(OUT / "aml_reuse_composite_key.json", "w"), indent=2)
print(json.dumps({s: {k: round(v, 4) if isinstance(v, float) else v for k, v in d.items()
                      if k in ("old_reuse", "old_distinct_donors", "achieved_reuse", "distinct_donors",
                               "scarcity_floor")} for s, d in res["stages"].items()}, indent=1))
print("pooled", res["published_pooled"])
