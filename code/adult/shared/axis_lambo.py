"""Maturity-axis analysis using Lambo et al. markers for a given module."""
import argparse
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import scoring_common as C

N_DRAWS = 2000
R05_PRIMITIVE = ["CD34", "KIT", "MEIS1", "HOXA9", "PROM1", "FLT3", "MECOM", "ERG", "GATA2", "MLLT3"]
R05_MATURE = ["ITGAM", "FCGR3A", "CD14", "LYZ", "CSF1R", "FCN1", "S100A8", "S100A9", "CEBPE", "MNDA"]
CACHE = C.BUNDLE["link"] / "out" / "axis_null_draws_shared.npz"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    log = lambda m: print(f"[axis:{a.tag}] {m}", flush=True)

    M = C.resolve_module(a.module, a.arm)
    OUT = C.out_dir(a.arm)
    REGS, BG = M["regulons"], M["background"]
    params = json.loads((OUT / "scoring_params.json").read_text())
    AUC_T, MIN_DET = params["auc_threshold"], params["min_detect_frac"]
    st, pairs, patients, nbm = C.load_pairs(a.arm)
    log(f"{M['name']} {len(REGS)} regulons | arm vocabulary '{M['arm']}' | "
        f"{len(patients)} paired patients, {len(nbm)} NBM")

    s3 = pd.read_excel(C.LAMBO_S3, sheet_name="Cell markers", header=0)
    s3.columns = [str(c).strip() for c in s3.columns]
    MARK = {c: [g for g in s3[c].dropna().astype(str).str.strip() if g and g.lower() != "nan"]
            for c in s3.columns}
    log("Lambo S3 cell types: " + ", ".join(f"{k}({len(v)})" for k, v in MARK.items()))

    E = C.Expr(log=log)
    E.selftest_auc(AUC_T, log=log)
    UNIVERSE = np.array([g for g in E.genes if E.DETECT[g] >= MIN_DET])
    log(f"{E.TOTAL} genes; expressed universe (detect>={MIN_DET}) = {len(UNIVERSE)}")

    def auc_of(genes):
        g = [x for x in genes if x in E.gset]
        if not g:
            return None, []
        return E.auc(g, AUC_T), g

    def delta(series):
        return C.delta_of(series, pairs, patients)

    MOD_DELTA = {}
    for vocab in C.DEFS:
        S = C.load_scores(a.arm, "aucell", vocab, BG)
        MOD_DELTA[vocab] = delta(C.module_score(S, REGS, M["orientation"][vocab]))
        log(f"{vocab}: module delta mean {MOD_DELTA[vocab].mean():+.5f}")

    TGT = {v: set().union(*[set(M["targets"][v][r]) for r in REGS]) for v in C.DEFS}
    for v in C.DEFS:
        log(f"{v}: module target union = {len(TGT[v])} genes "
            f"({len(TGT[v] & E.gset)} present, {len(TGT[v] & E.gset)/E.TOTAL:.1%} of the matrix)")

    def fit_axis(d_mod, d_axis):
        r = stats.pearsonr(d_mod, d_axis)
        f = sm.OLS(d_mod, sm.add_constant(d_axis)).fit()
        resid = d_mod - f.params[1] * d_axis
        w = stats.wilcoxon(resid, zero_method="wilcox")
        return {"corr": float(r[0]), "corr_p": float(r[1]), "r_squared": float(f.rsquared),
                "adjusted_intercept": float(f.params[0]),
                "adjusted_intercept_p": float(f.pvalues[0]),
                "slope": float(f.params[1]), "slope_p": float(f.pvalues[1]),
                "residualised_mean_delta": float(resid.mean()),
                "residualised_wilcoxon_p": float(w.pvalue),
                "survives_adjustment": bool(f.params[0] > 0 and f.pvalues[0] < 0.05)}

    AXES = {"A1_lambo_HSC_minus_Monocyte": (MARK["HSC"], MARK["Monocyte"]),
            "A2_lambo_HSCProg_minus_MonoGMP": (MARK["HSC"] + MARK["Progenitor"],
                                               MARK["Monocyte"] + MARK["GMP"]),
            "A4_r05_handbuilt_reference": (R05_PRIMITIVE, R05_MATURE)}

    # shared, module-independent null draws
    
    if CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        DRAWS = {k: z[k] for k in z.files if not k.endswith("_deg")}
        DEG = {k[:-4]: int(z[k]) for k in z.files if k.endswith("_deg")}
        log(f"loaded shared axis null draws from {CACHE.name} "
            f"({', '.join(f'{k}:{v.shape[0]}' for k, v in DRAWS.items())})")
    else:
        rng = np.random.default_rng(C.SEED)
        det_universe = E.DETECT[UNIVERSE].values
        order = np.argsort(det_universe)
        uni_sorted, det_sorted = UNIVERSE[order], det_universe[order]

        def match_by_detection(genes):
            out = []
            for g in genes:
                j = int(np.searchsorted(det_sorted, E.DETECT[g]))
                lo, hi = max(0, j - 25), min(len(uni_sorted), j + 25)
                cand = [x for x in uni_sorted[lo:hi] if x not in out]
                out.append(rng.choice(cand) if cand else rng.choice(uni_sorted))
            return out

        def null_axis_deltas(n_pos, n_neg, pos_genes, neg_genes, mode):
            out, n_degenerate = [], 0
            for _ in range(N_DRAWS):
                if mode == "size":
                    gp = list(rng.choice(UNIVERSE, n_pos, replace=False))
                    gn = list(rng.choice(UNIVERSE, n_neg, replace=False))
                else:
                    gp, gn = match_by_detection(pos_genes), match_by_detection(neg_genes)
                ap_, _ = auc_of(gp)
                an_, _ = auc_of(gn)
                d = delta(ap_ - an_)
                if np.ptp(d) == 0 or not np.all(np.isfinite(d)):
                    n_degenerate += 1
                    continue
                out.append(d)
            return np.array(out), n_degenerate

        DRAWS, DEG = {}, {}
        for name, (pos, neg) in AXES.items():
            _, gp = auc_of(pos)
            _, gn = auc_of(neg)
            for mode, key in (("size", "size"), ("detect", "detect")):
                d, deg = null_axis_deltas(len(gp), len(gn), gp, gn, mode)
                DRAWS[f"{name}__{key}"] = d
                DEG[f"{name}__{key}"] = deg
                log(f"drew {name} {key}-matched null: {len(d)} valid ({deg} degenerate)")
        np.savez_compressed(CACHE, **DRAWS, **{f"{k}_deg": np.array(v) for k, v in DEG.items()})
        log(f"cached shared axis null draws -> {CACHE}")

    def null_summary(d_mod, draws, n_degenerate, observed_r2):
        r2 = np.array([stats.pearsonr(d_mod, d)[0] ** 2 for d in draws])
        r2 = r2[np.isfinite(r2)]
        return {"n_draws_requested": N_DRAWS, "n_draws_valid": int(len(r2)),
                "n_draws_degenerate_discarded": int(n_degenerate),
                "median_r2": float(np.median(r2)), "mean_r2": float(r2.mean()),
                "p95_r2": float(np.percentile(r2, 95)),
                "ci95": [float(np.percentile(r2, 2.5)), float(np.percentile(r2, 97.5))],
                "p_empirical": float((np.sum(r2 >= observed_r2) + 1) / (len(r2) + 1))}

    res = {"provenance": C.provenance(M, {
               "null_draws": "shared, module-independent; drawn once with rng(555) in the r_06 "
                             "order (A1 size, A1 detect, A2 size, A2 detect, A4 size, A4 detect) "
                             f"and cached at {CACHE.name}"}),
           "status": "POST-HOC maturity-axis analysis under the module "
                     "redefinition (ax4_tier1 retired). The gene-disjoint variant is recomputed "
                     "against THIS module's target union.",
           "seed": C.SEED, "n_draws": N_DRAWS, "n_patients": len(patients),
           "auc_threshold": AUC_T, "min_detect_frac": MIN_DET,
           "expressed_universe": int(len(UNIVERSE)),
           "lambo_S3_cell_types": {k: len(v) for k, v in MARK.items()},
           "module_delta": {v: float(MOD_DELTA[v].mean()) for v in C.DEFS},
           "module_target_union": {v: {"n_genes": int(len(TGT[v])),
                                       "n_present_in_matrix": int(len(TGT[v] & E.gset)),
                                       "frac_of_matrix": float(len(TGT[v] & E.gset) / E.TOTAL)}
                                   for v in C.DEFS},
           "axes": {}}

    for name, (pos, neg) in AXES.items():
        ap_, gp = auc_of(pos)
        an_, gn = auc_of(neg)
        ax = ap_ - an_
        d_ax = delta(ax)
        entry = {"pos_genes_requested": len(pos), "pos_genes_used": len(gp),
                 "neg_genes_requested": len(neg), "neg_genes_used": len(gn),
                 "axis_means": {"NBM": float(ax[nbm].mean()),
                                "Dx": float(ax[[pairs[p][0] for p in patients]].mean()),
                                "Rel": float(ax[[pairs[p][1] for p in patients]].mean())},
                 "delta_axis": {"mean": float(d_ax.mean()),
                                "p_wilcoxon": float(stats.wilcoxon(d_ax, zero_method="wilcox").pvalue),
                                "frac_positive": float((d_ax > 0).mean())},
                 "by_vocabulary": {}}
        log("")
        log(f"### {name}: {len(gp)}/{len(pos)} pos, {len(gn)}/{len(neg)} neg genes present; "
            f"delta {d_ax.mean():+.5f}")

        dr1, deg1 = DRAWS[f"{name}__size"], DEG[f"{name}__size"]
        dr2, deg2 = DRAWS[f"{name}__detect"], DEG[f"{name}__detect"]

        for vocab in C.DEFS:
            d_mod = MOD_DELTA[vocab]
            f = fit_axis(d_mod, d_ax)

            tgt = TGT[vocab]
            dp = [g for g in gp if g not in tgt]
            dn = [g for g in gn if g not in tgt]
            if dp and dn:
                adp, _ = auc_of(dp)
                adn, _ = auc_of(dn)
                fd = fit_axis(d_mod, delta(adp - adn))
                fd["pos_genes_kept"], fd["neg_genes_kept"] = len(dp), len(dn)
                fd["pos_genes_dropped"] = sorted(set(gp) - set(dp))
                fd["neg_genes_dropped"] = sorted(set(gn) - set(dn))
                fd["module_target_union_size"] = int(len(tgt))
            else:
                fd = {"unavailable": True,
                      "reason": f"the module's {vocab} target union ({len(tgt)} genes) absorbs a "
                                f"whole pole of the axis: {len(dp)} positive / {len(dn)} negative "
                                f"marker genes survive the disjoint filter, so the variant cannot "
                                f"be formed",
                      "pos_genes_kept": len(dp), "neg_genes_kept": len(dn),
                      "pos_genes_dropped": sorted(set(gp) - set(dp)),
                      "neg_genes_dropped": sorted(set(gn) - set(dn)),
                      "module_target_union_size": int(len(tgt))}

            s1 = null_summary(d_mod, dr1, deg1, f["r_squared"])
            s2 = null_summary(d_mod, dr2, deg2, f["r_squared"])
            log(f"    [{vocab}] full: r={f['corr']:+.3f} R2={f['r_squared']:.4f} "
                f"adj.int={f['adjusted_intercept']:+.5f} (p={f['adjusted_intercept_p']:.3g}) "
                f"survives={f['survives_adjustment']}")
            log(f"      N1 size R2 med {s1['median_r2']:.4f} p={s1['p_empirical']:.4g} | "
                f"N2 detect R2 med {s2['median_r2']:.4f} p={s2['p_empirical']:.4g}")
            if "unavailable" in fd:
                log(f"      gene-disjoint: UNAVAILABLE — {fd['pos_genes_kept']} pos / "
                    f"{fd['neg_genes_kept']} neg genes survive")
            else:
                log(f"      gene-disjoint ({fd['pos_genes_kept']}+{fd['neg_genes_kept']} genes kept "
                    f"of {len(gp)}+{len(gn)}): r={fd['corr']:+.3f} R2={fd['r_squared']:.4f} "
                    f"adj.int={fd['adjusted_intercept']:+.5f} "
                    f"(p={fd['adjusted_intercept_p']:.3g}) survives={fd['survives_adjustment']}")

            entry["by_vocabulary"][vocab] = {"primary": f, "gene_disjoint": fd,
                                             "null_size_matched": s1, "null_detection_matched": s2}
        res["axes"][name] = entry

    hsc, ghsc = auc_of(MARK["HSC"])
    d_hsc = delta(hsc)
    res["A3_lambo_HSC_alone"] = {
        "genes_used": len(ghsc),
        "delta": {"mean": float(d_hsc.mean()),
                  "p_wilcoxon": float(stats.wilcoxon(d_hsc, zero_method="wilcox").pvalue),
                  "frac_positive": float((d_hsc > 0).mean())},
        "by_vocabulary": {v: fit_axis(MOD_DELTA[v], d_hsc) for v in C.DEFS}}
    log("")
    log(f"### A3 Lambo HSC alone ({len(ghsc)} genes): delta {d_hsc.mean():+.5f}")
    for v in C.DEFS:
        x = res["A3_lambo_HSC_alone"]["by_vocabulary"][v]
        log(f"    [{v}] R2={x['r_squared']:.4f} adj.int={x['adjusted_intercept']:+.5f} "
            f"(p={x['adjusted_intercept_p']:.3g}) survives={x['survives_adjustment']}")

    (OUT / f"axis_lambo_{a.tag}.json").write_text(json.dumps(res, indent=2, default=str))
    log(f"wrote axis_lambo_{a.tag}.json")

if __name__ == "__main__":
    main()
