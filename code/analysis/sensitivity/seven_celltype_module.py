"""SENSITIVITY (Supplementary Data 6), step 2: write the module definition of the seven-type,
5-of-7 sensitivity module in the schema read by adult/shared/nullB_pool.py and
test.py (which are run unmodified).

Reads sens_clp7_summary.json (step 1). Orientation, target sets and background are carried over
VERBATIM from modules/out/corrected_module.json (they exist for all 201 regulons);
only the membership fields change. The primary module file is not modified.

Run from the project root with --repo . --out analysis/out/seven_celltypes.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
repo, out = Path(a.repo), Path(a.out)

src = json.loads((repo / "modules" / "out" / "corrected_module.json").read_text())
summ = json.loads((out / "sens_clp7_summary.json").read_text())
regs = list(summ["sensitivity_module"]["regulons"])
tfs = sorted({re.split(r"[_(]", r)[0] for r in regs})
assert set(regs) <= set(src["background_regulons"]), "members outside the 201 background"
assert all(r in src["orientation"][v] and r in src["targets"][v] for r in regs for v in ("Diagnosis", "Relapse"))

mod = dict(src)
mod["name"] = "SENS_CLP7_5of7"
mod["frozen"] = ("SENSITIVITY module, not a result of record: seven admitted cell types (CLP added), "
                 "cell-type bar 5 of 7; derived by seven_celltype_nomination.py; "
                 "orientation, targets and background carried verbatim from T3_repair")
mod["regulons"] = regs
mod["n_regulons"] = len(regs)
mod["TFs"] = tfs
mod["n_TFs"] = len(tfs)
def _sizes(v):
    n = [len(src["targets"][v][r]) for r in regs]
    return {"SENS_CLP7_5of7_median": float(np.median(n)), "SENS_CLP7_5of7_min": int(min(n)),
            "SENS_CLP7_5of7_max": int(max(n)), "SENS_CLP7_5of7_n_regulons_measured": len(n)}

mod["target_set_sizes"] = {v: _sizes(v) for v in ("Diagnosis", "Relapse")}
(out / "module_definition_SENS_CLP7_5of7.json").write_text(json.dumps(mod, indent=2))
print(f"wrote module_definition_SENS_CLP7_5of7.json: {len(regs)} regulons / {len(tfs)} TFs")
