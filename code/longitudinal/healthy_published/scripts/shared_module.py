"""Characterise the shared (non-selective) regulons, including enrichment of immediate-early/AP-1 TFs."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = (Path(__file__).resolve().parents[3] / 'longitudinal/healthy_published')
OUT = ROOT / "out"

# canonical immediate-early / AP-1 programme, fixed before inspection of the verdict table
IEG = {"FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
       "ATF3", "KLF2", "KLF4", "KLF6", "NR4A1", "IER2", "ZFP36", "DUSP1", "NFIL3"}

def log(m):
    print(f"[s3b-module] {m}", flush=True)

w = pd.read_csv(OUT / "did_reconciliation.csv", index_col=0)
w["is_ieg"] = w.TF.isin(IEG)
w["is_shared"] = w.verdict == "SHARED"

tab = pd.crosstab(w.is_ieg, w.is_shared)
odds, p = stats.fisher_exact(tab.reindex(index=[False, True], columns=[False, True]).values)

res = {"status": "descriptive / exploratory",
       "ieg_definition": sorted(IEG),
       "n_regulons_tested": int(len(w)),
       "n_ieg_regulons": int(w.is_ieg.sum()),
       "shared_rate_ieg": float(w[w.is_ieg].is_shared.mean()),
       "shared_rate_other": float(w[~w.is_ieg].is_shared.mean()),
       "odds_ratio": float(odds), "fisher_p": float(p),
       "contingency": {"ieg_shared": int(tab.loc[True, True]) if True in tab.index else 0,
                       "ieg_not_shared": int(tab.loc[True, False]) if True in tab.index else 0,
                       "other_shared": int(tab.loc[False, True]),
                       "other_not_shared": int(tab.loc[False, False])},
       "ieg_TFs_present": sorted(w[w.is_ieg].TF.unique().tolist())}

log(f"IEG/AP-1 regulons: {res['n_ieg_regulons']} of {res['n_regulons_tested']} "
    f"({', '.join(res['ieg_TFs_present'])})")
log(f"SHARED rate: IEG {res['shared_rate_ieg']:.3f} vs other {res['shared_rate_other']:.3f} "
    f"-> OR {res['odds_ratio']:.2f}, Fisher p = {res['fisher_p']:.4g}")

# and the mirror image: what does the SELECTIVE set look like?
sel = w[w.verdict == "DIRECTIONAL-SELECTIVE"]
res["selective_top_TFs"] = sel.TF.value_counts().head(15).to_dict()
res["shared_top_TFs"] = w[w.is_shared].TF.value_counts().head(15).to_dict()
log(f"SELECTIVE set: {len(sel)} regulons over {sel.TF.nunique()} TFs")

# healthy effect size inside vs outside the IEG set — is the shared module simply the one the
# healthy compartment moves on at all?
for label, g in (("IEG", w[w.is_ieg]), ("other", w[~w.is_ieg])):
    res[f"healthy_median_abs_beta_{label}"] = float(g.beta_healthy_Diagnosis.abs().median())
log(f"healthy |beta| (Dx vocab): IEG {res['healthy_median_abs_beta_IEG']:.5f} vs "
    f"other {res['healthy_median_abs_beta_other']:.5f}")

with open(OUT / "s3b_shared_module.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
log("wrote s3b_shared_module.json")
