"""FIG1 panel c diagnostic (supplementary) - does the all-shared per-regulon number reproduce"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "paml/malignant" / "out"
OUT = ROOT / "analysis" / "out"

STAGES = ["Diagnosis", "Relapse"]
ARMS = ["link", "pseudomultiome"]
FLOOR = 20

def log(m):
    print(f"[f1-floor] {m}", flush=True)

def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

def read_profile(stage, arm, ct_map):
    """(cell type x regulon) mean Gene_based AUC, and the unit count per cell type."""
    with h5py.File(S2 / f"eRegulon_AUC_{stage}_{arm}.h5mu", "r") as h:
        cells = dec(h["obs/Cell"][:])
        var = dec(h["mod/Gene_based/var/_index"][:])
        X = h["mod/Gene_based/X"][:]
    reg = [v.rsplit("_(", 1)[0] for v in var]
    if arm == "pseudomultiome":
        ct = np.array([c.rsplit("_", 1)[0] for c in cells])
    else:
        ct = np.array([ct_map.get(c, "__UNMATCHED__") for c in cells])
        assert (ct != "__UNMATCHED__").all(), f"{stage}/{arm}: composite-key join incomplete"
    M = pd.DataFrame(X, columns=reg)
    M["__ct__"] = ct
    prof = M.groupby("__ct__").mean()
    counts = pd.Series(ct).value_counts()
    return prof, counts

def median_matched(pa, pb, cts, regs):
    """Median matched."""
    A = pa.loc[cts, regs].values.astype(float)
    B = pb.loc[cts, regs].values.astype(float)
    zA = (A - A.mean(0)) / A.std(0, ddof=0)
    zB = (B - B.mean(0)) / B.std(0, ddof=0)
    C = (zA.T @ zB) / len(cts)
    m = np.diag(C)
    m = m[np.isfinite(m)]
    return float(np.median(m)), m

res = {"status": "DIAGNOSTIC - report only",
       "committed_source": ("analysis/network_comparison/out/results.json :: convergence[].act_median_r, "
                            "produced by paml/healthy/scripts/pairing_convergence.py"),
       "estimator": ("z within arm across cell types; C = (A.T @ B)/n_celltypes; "
                     "act_median_r = median(diag(C)). Identical to f1_07's per-regulon Pearson."),
       "by_stage": {}}

committed = {c["stage"]: c for c in
             json.load(open(ROOT / "analysis/network_comparison" / "out" / "results.json"))["convergence"]}

for stage in STAGES:
    cd = pd.read_csv(S2 / f"cell_data_{stage}.csv", index_col=0)
    ct_map = cd["Classified_Celltype"].to_dict()
    prof, cnt = {}, {}
    for arm in ARMS:
        prof[arm], cnt[arm] = read_profile(stage, arm, ct_map)
    regs = sorted(set(prof["link"].columns) & set(prof["pseudomultiome"].columns))

    cts_all = sorted(set(prof["link"].index) & set(prof["pseudomultiome"].index))
    keep = {a: {c for c in cnt[a].index if cnt[a][c] >= FLOOR} for a in ARMS}
    cts_floor = sorted(set(cts_all) & keep["link"] & keep["pseudomultiome"])

    m_all, v_all = median_matched(prof["link"], prof["pseudomultiome"], cts_all, regs)
    m_flr, v_flr = median_matched(prof["link"], prof["pseudomultiome"], cts_floor, regs)

    com = committed[stage]
    res["by_stage"][stage] = {
        "n_shared_regulons": len(regs),
        "committed_n_celltypes": com["cell_types"],
        "committed_shared_regulons": com["shared_regulons"],
        "committed_act_median_r": com["act_median_r"],
        "committed_ci95": [com["act_lo"], com["act_hi"]],
        "recomputed_NO_floor": {"n_celltypes": len(cts_all), "celltypes": cts_all,
                                "median_r": m_all,
                                "abs_diff_vs_committed": abs(m_all - com["act_median_r"])},
        "recomputed_WITH_20_unit_floor": {
            "n_celltypes": len(cts_floor), "celltypes": cts_floor, "median_r": m_flr,
            "celltypes_dropped_by_floor": sorted(set(cts_all) - set(cts_floor)),
            "unit_counts_of_dropped": {c: {a: int(cnt[a].get(c, 0)) for a in ARMS}
                                       for c in sorted(set(cts_all) - set(cts_floor))}},
        "floor_cost_median_r": m_flr - m_all,
        "inside_committed_ci95": bool(com["act_lo"] <= m_flr <= com["act_hi"]),
        "PARITY_no_floor_reproduces_committed": bool(abs(m_all - com["act_median_r"]) < 1e-9),
    }
    r = res["by_stage"][stage]
    log(f"{stage}: shared regulons {len(regs)} (committed {com['shared_regulons']}) | "
        f"no floor {len(cts_all)} cts -> {m_all:.6f} vs committed {com['act_median_r']:.6f} "
        f"(diff {r['recomputed_NO_floor']['abs_diff_vs_committed']:.2e}) | "
        f"floored {len(cts_floor)} cts -> {m_flr:.6f} | dropped "
        f"{r['recomputed_WITH_20_unit_floor']['celltypes_dropped_by_floor']}")

res["ALL_PARITY_PASS"] = bool(all(v["PARITY_no_floor_reproduces_committed"]
                                  for v in res["by_stage"].values()))
with open(OUT / "f1_celltype_floor_sensitivity.json", "w") as f:
    json.dump(res, f, indent=1)
log(f"parity (no-floor reproduces committed activity-level result): {res['ALL_PARITY_PASS']}")
log("wrote f1_celltype_floor_sensitivity.json")
