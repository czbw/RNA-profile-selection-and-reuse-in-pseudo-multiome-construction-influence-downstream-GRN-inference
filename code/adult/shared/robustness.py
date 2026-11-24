"""Robustness checks (purity confound, leave-one-out, per-regulon agreement, centring) for a given module."""
import argparse
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import scoring_common as C

N_SIG = 200

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    log = lambda m: print(f"[robust:{a.tag}] {m}", flush=True)

    M = C.resolve_module(a.module, a.arm)
    OUT = C.out_dir(a.arm)
    REGS, BG = M["regulons"], M["background"]
    params = json.loads((OUT / "scoring_params.json").read_text())
    st, pairs, patients, nbm = C.load_pairs(a.arm)
    log(f"{M['name']} {len(REGS)} regulons | arm vocabulary '{M['arm']}' | {len(patients)} patients")

    E = C.Expr(log=log)
    E.selftest_auc(params["auc_threshold"], log=log)

    res = {"provenance": C.provenance(M),
           "status": "POST-HOC — not part of the pre-specified analysis; prompted by the NBM ordering and the "
                     "internal-vs-external magnitude gap",
           "n_patients": len(patients)}

    # ---------------- R1 purity ----------------
    cpm = E.X / E.X.sum(axis=0) * 1e6
    lg = np.log2(cpm + 1)
    d = lg[nbm].mean(axis=1) - lg[[s for p in patients for s in pairs[p]]].mean(axis=1)
    norm_sig = list(d.sort_values(ascending=False).head(N_SIG).index)
    res["R1_purity"] = {"signature_size": len(norm_sig), "top_genes": norm_sig[:15],
                        "definition": "top genes by mean log2CPM difference, NBM vs all AML samples",
                        "membership_independent": True}

    proxy = E.auc(norm_sig, params["auc_threshold"])
    dprox = C.delta_of(proxy, pairs, patients)
    res["R1_purity"]["proxy_means"] = {"NBM": float(proxy[nbm].mean()),
                                       "Dx": float(proxy[[pairs[p][0] for p in patients]].mean()),
                                       "Rel": float(proxy[[pairs[p][1] for p in patients]].mean())}
    res["R1_purity"]["delta_proxy_mean"] = float(dprox.mean())
    res["R1_purity"]["delta_proxy_p"] = float(stats.wilcoxon(dprox, zero_method="wilcox").pvalue)
    log(f"R1 proxy: NBM {proxy[nbm].mean():.4f} | Dx {res['R1_purity']['proxy_means']['Dx']:.4f} | "
        f"Rel {res['R1_purity']['proxy_means']['Rel']:.4f}; delta {dprox.mean():+.5f} "
        f"(p={res['R1_purity']['delta_proxy_p']:.3g})")

    res["R1_purity"]["by_vocabulary"] = {}
    for vocab in C.DEFS:
        signs = M["orientation"][vocab]
        S = C.load_scores(a.arm, "aucell", vocab, BG)
        sc = C.module_score(S, REGS, signs)
        dmod = C.delta_of(sc, pairs, patients)
        r_p = stats.pearsonr(dmod, dprox)
        r_s = stats.spearmanr(dmod, dprox)
        fit = sm.OLS(dmod, sm.add_constant(dprox)).fit()
        resid = dmod - (fit.params[1] * dprox)
        w_res = stats.wilcoxon(resid, zero_method="wilcox")
        res["R1_purity"]["by_vocabulary"][vocab] = {
            "corr_delta_module_vs_delta_purity_pearson": float(r_p[0]), "corr_p": float(r_p[1]),
            "spearman": float(r_s.statistic),
            "intercept_adjusted_mean_delta": float(fit.params[0]),
            "intercept_p": float(fit.pvalues[0]),
            "slope_on_purity": float(fit.params[1]), "slope_p": float(fit.pvalues[1]),
            "r_squared": float(fit.rsquared),
            "purity_residualised_mean_delta": float(resid.mean()),
            "purity_residualised_wilcoxon_p": float(w_res.pvalue),
            "survives_purity_adjustment": bool(fit.params[0] > 0 and fit.pvalues[0] < 0.05)}
        x = res["R1_purity"]["by_vocabulary"][vocab]
        log(f"R1 {vocab}: corr {x['corr_delta_module_vs_delta_purity_pearson']:+.3f} "
            f"(p={x['corr_p']:.3g}) R2={x['r_squared']:.4f}; adjusted intercept "
            f"{x['intercept_adjusted_mean_delta']:+.5f} (p={x['intercept_p']:.3g}) -> survives "
            f"{x['survives_purity_adjustment']}")

    # ---------------- R2 / R3 ----------------
    res["R2_leave_one_out"], res["R3_per_regulon"] = {}, {}
    for vocab in C.DEFS:
        signs = M["orientation"][vocab]
        S = C.load_scores(a.arm, "aucell", vocab, BG)
        Cn = pd.DataFrame(S.values - S.values.mean(axis=1, keepdims=True),
                          index=S.index, columns=S.columns)
        a_dx = [pairs[p][0] for p in patients]
        a_rl = [pairs[p][1] for p in patients]

        loo = []
        for drop in REGS:
            keep = [r for r in REGS if r != drop]
            sg = np.array([signs[r] for r in keep], float)
            s = pd.Series((Cn[keep].values * sg).mean(axis=1), index=Cn.index)
            dd = np.array([s[b] - s[aa] for aa, b in zip(a_dx, a_rl)])
            loo.append(dict(dropped=drop, mean_delta=float(dd.mean()),
                            p=float(stats.wilcoxon(dd, zero_method="wilcox").pvalue)))
        L = pd.DataFrame(loo)
        res["R2_leave_one_out"][vocab] = {
            "n_members_swept": int(len(L)),
            "min_mean_delta": float(L.mean_delta.min()), "max_mean_delta": float(L.mean_delta.max()),
            "max_p": float(L.p.max()), "all_still_positive": bool((L.mean_delta > 0).all()),
            "all_still_p05": bool((L.p < 0.05).all()),
            "most_influential": L.loc[L.mean_delta.idxmin(), "dropped"]}
        x = res["R2_leave_one_out"][vocab]
        log(f"R2 {vocab}: sweeping {x['n_members_swept']} members -> delta in "
            f"[{x['min_mean_delta']:+.5f},{x['max_mean_delta']:+.5f}], worst p={x['max_p']:.3g}, "
            f"all positive {x['all_still_positive']}, all p<0.05 {x['all_still_p05']}")

        per = []
        for r in REGS:
            dd = np.array([Cn.loc[b, r] - Cn.loc[aa, r] for aa, b in zip(a_dx, a_rl)]) * signs[r]
            per.append(dict(regulon=r, oriented_mean_delta=float(dd.mean()),
                            frac_positive=float((dd > 0).mean()),
                            p=float(stats.wilcoxon(dd, zero_method="wilcox").pvalue),
                            agrees_with_internal_sign=bool(dd.mean() > 0)))
        P = pd.DataFrame(per).sort_values("oriented_mean_delta", ascending=False)
        P.to_csv(OUT / f"per_regulon_external_{a.tag}_{vocab}.csv", index=False)
        res["R3_per_regulon"][vocab] = {
            "n_agreeing_with_internal_sign": int(P.agrees_with_internal_sign.sum()),
            "n_total": int(len(P)), "n_individually_p05": int((P.p < 0.05).sum()),
            "top3": P.head(3).regulon.tolist(), "bottom3": P.tail(3).regulon.tolist()}
        x = res["R3_per_regulon"][vocab]
        log(f"R3 {vocab}: {x['n_agreeing_with_internal_sign']}/{x['n_total']} move in their "
            f"internal direction; {x['n_individually_p05']} individually p<0.05")

    # ---------------- R4 ----------------
    res["R4_centring"] = {}
    for vocab in C.DEFS:
        S = C.load_scores(a.arm, "aucell", vocab, BG)
        off = pd.Series(S.values.mean(axis=1), index=S.index)
        doff = C.delta_of(off, pairs, patients)
        sc = C.module_score(S, REGS, M["orientation"][vocab])
        dmod = C.delta_of(sc, pairs, patients)
        r = stats.pearsonr(dmod, doff)
        res["R4_centring"][vocab] = {"global_offset_delta_mean": float(doff.mean()),
                                     "corr_module_vs_offset": float(r[0]), "p": float(r[1])}
        log(f"R4 {vocab}: global offset delta {doff.mean():+.5f}; "
            f"corr(dModule,dOffset) = {r[0]:+.3f} (p={r[1]:.3g})")

    (OUT / f"robustness_{a.tag}.json").write_text(json.dumps(res, indent=2, default=str))
    log(f"wrote robustness_{a.tag}.json")

if __name__ == "__main__":
    main()
