"""Evaluate oriented patient changes and three empirical null comparisons.

Module membership is supplied explicitly. Historical decision-rule fields are
retained for output compatibility; the current 55-member selection is post hoc."""
import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

import scoring_common as C

N_DRAWS = 2000
METHODS = ("aucell", "singscore", "meanz")

historical_decision_rule = ("Pre-specified decision rule: the check passes if, under both vocabularies: "
             "1. the per-patient delta is positive (agrees with the internal direction), and "
             "2. the paired Wilcoxon reaches p < 0.05, and 3. it clears NULL A at p < 0.05. "
             "Clearing NULL C additionally is recorded as the strong form of the result. Failing "
             "null C while passing 1-3 means the module moves but is not distinguishable from an "
             "arbitrary internal regulon set - reported as such, not as a pass.")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--arm", required=True, choices=["link", "repair"])
    ap.add_argument("--tag", required=True)
    ap.add_argument("--nullb-tag", default=None,
                    help="basename of the NULL B pool to use; defaults to --tag. Only set this "
                         "for the frozen-old-module reproduction check.")
    ap.add_argument("--out-suffix", default=None)
    ap.add_argument("--nullb-name-fmt", default="nullB_pool_{tag}_{vocab}.npz",
                    help="filename pattern for the NULL B pool; the default is the re-drawn, "
                         "size-matched pool. Set to 'nullB_pool_{vocab}.npz' ONLY for the "
                         "frozen-old-module reproduction check against the committed pool.")
    a = ap.parse_args()
    tag = a.tag
    suffix = a.out_suffix or tag
    nbtag = a.nullb_tag or tag
    log = lambda m: print(f"[test:{tag}] {m}", flush=True)

    M = C.resolve_module(a.module, a.arm)
    OUT = C.out_dir(a.arm)
    module_regulons, background_regulons = M["regulons"], M["background"]
    st, pairs, patients, nbm = C.load_pairs(a.arm)
    log(f"{M['name']} {len(module_regulons)} regulons / {M['n_TFs']} TFs | arm vocabulary '{M['arm']}' | "
        f"background {len(background_regulons)} | paired patients {len(patients)} | NBM {len(nbm)}")

    res = {"provenance": C.provenance(M, {
               "scores_source": "committed scores_{method}_{vocab}.csv — membership-INDEPENDENT "
                                "(columns are the arm's full background; values depend only on the "
                                "per-regulon target sets, which the new definitions carry over "
                                "verbatim)",
               "nullB_pool": a.nullb_name_fmt.format(tag=nbtag, vocab="<vocab>")}),
           "seed": C.SEED, "n_draws": N_DRAWS, "n_paired_patients": len(patients),
           "preregistered_decision_rule": historical_decision_rule,
           "decision_rule": "delta>0 AND wilcoxon p<0.05 AND nullA p<0.05, under BOTH vocabularies",
           "primary_scale": "offset-centred", "primary_estimator": "aucell",
           "nulls_used": {}, "by_vocabulary": {}}
    delta_rows = []

    for vocab in C.DEFS:
        signs = M["orientation"][vocab]
        S = C.load_scores(a.arm, "aucell", vocab, background_regulons)
        sc = C.module_score(S, module_regulons, signs)
        d = C.delta_of(sc, pairs, patients)
        obs = float(d.mean())
        w = stats.wilcoxon(d, zero_method="wilcox")
        for p, v in zip(patients, d):
            delta_rows.append(dict(vocabulary=vocab, patient=p, delta=v))

        # ---- NULL A ----
        r2 = np.random.default_rng(C.SEED)
        nullA = np.array([np.mean(d * r2.choice([-1.0, 1.0], size=len(d))) for _ in range(N_DRAWS)])
        pA = float((np.sum(np.abs(nullA) >= abs(obs)) + 1) / (N_DRAWS + 1))

        # ---- NULL B (size-matched pool) ----
        z = np.load(OUT / a.nullb_name_fmt.format(tag=nbtag, vocab=vocab), allow_pickle=True)
        pool, pool_samples = z["pool"], list(z["samples"])
        if pool.shape[1] != len(module_regulons):
            raise SystemExit(f"[test] NULL B pool has {pool.shape[1]} slots but the module has "
                             f"{len(module_regulons)} members — refusing to score a module against a null "
                             f"size-matched to a different module")
        off = pd.Series(S.values.mean(axis=1), index=S.index)
        sg = np.array([signs[r] for r in module_regulons], dtype=float)
        r3 = np.random.default_rng(C.SEED + 1)
        pix = {s: i for i, s in enumerate(pool_samples)}
        idx_dx = [pix[pairs[p][0]] for p in patients]
        idx_rl = [pix[pairs[p][1]] for p in patients]
        offd = np.array([off[pairs[p][1]] - off[pairs[p][0]] for p in patients])
        nullB = np.empty(N_DRAWS)
        for b in range(N_DRAWS):
            pick = r3.integers(0, pool.shape[2], size=pool.shape[1])
            sel = pool[:, np.arange(pool.shape[1]), pick]
            m = (sel * sg).mean(axis=1)
            nullB[b] = np.mean((m[idx_rl] - m[idx_dx]) - offd * sg.mean())
        pB = float((np.sum(np.abs(nullB) >= abs(obs)) + 1) / (N_DRAWS + 1))

        # ---- NULL C (random modules of the SAME size from the background) ----
        Cm = S.values - S.values.mean(axis=1, keepdims=True)
        allsg = np.array([signs[r] for r in S.columns], dtype=float)
        six = {s: i for i, s in enumerate(S.index)}
        a_dx = np.array([six[pairs[p][0]] for p in patients])
        a_rl = np.array([six[pairs[p][1]] for p in patients])
        r4 = np.random.default_rng(C.SEED + 2)
        nullC = np.empty(N_DRAWS)
        for b in range(N_DRAWS):
            pick = r4.choice(Cm.shape[1], size=len(module_regulons), replace=False)
            m = (Cm[:, pick] * allsg[pick]).mean(axis=1)
            nullC[b] = np.mean(m[a_rl] - m[a_dx])
        pC = float((np.sum(np.abs(nullC) >= abs(obs)) + 1) / (N_DRAWS + 1))

        res["nulls_used"][vocab] = {
            "nullA_label_permutation": {"n_draws": N_DRAWS, "seed": C.SEED,
                                        "unit": f"sign flip of each of {len(patients)} patient deltas"},
            "nullB_random_gene_sets": {"n_draws": N_DRAWS, "seed": C.SEED + 1,
                                       "n_slots": int(pool.shape[1]),
                                       "sets_per_slot": int(pool.shape[2]),
                                       "n_random_gene_sets_scored": int(pool.shape[1] * pool.shape[2]),
                                       "slot_sizes_min": int(z["sizes"].min()),
                                       "slot_sizes_max": int(z["sizes"].max()),
                                       "slot_sizes_median": float(np.median(z["sizes"])),
                                       "size_matched_to": f"{len(module_regulons)}-regulon module {M['name']}"},
            "nullC_random_regulon_modules": {"n_draws": N_DRAWS, "seed": C.SEED + 2,
                                             "module_size_drawn": len(module_regulons),
                                             "background_pool": len(background_regulons)}}

        sens = {}
        for meth in METHODS[1:]:
            S2 = C.load_scores(a.arm, meth, vocab, background_regulons)
            s2 = C.module_score(S2, module_regulons, signs)
            d2 = C.delta_of(s2, pairs, patients)
            sens[meth] = {"mean_delta": float(d2.mean()),
                          "frac_positive": float((d2 > 0).mean()),
                          "p_wilcoxon": float(stats.wilcoxon(d2, zero_method="wilcox").pvalue),
                          "same_sign_as_primary": bool(np.sign(d2.mean()) == np.sign(obs))}

        ixR = [S.columns.get_loc(r) for r in module_regulons]
        raw = pd.Series((S.values[:, ixR] * sg).mean(axis=1), index=S.index)
        draw = C.delta_of(raw, pairs, patients)

        v = {"mean_delta": obs, "median_delta": float(np.median(d)),
             "frac_patients_positive": float((d > 0).mean()),
             "p_wilcoxon": float(w.pvalue),
             "nullA_label_permutation": {"mean": float(nullA.mean()), "sd": float(nullA.std()),
                                         "p_empirical": pA},
             "nullB_random_gene_sets": {"mean": float(nullB.mean()), "sd": float(nullB.std()),
                                        "p_empirical": pB},
             "nullC_random_regulon_modules": {"mean": float(nullC.mean()), "sd": float(nullC.std()),
                                              "p_empirical": pC,
                                              "percentile_of_observed":
                                                  float(stats.percentileofscore(np.abs(nullC), abs(obs)))},
             "sensitivities": sens,
             "raw_scale": {"mean_delta": float(draw.mean()),
                           "frac_positive": float((draw > 0).mean()),
                           "p_wilcoxon": float(stats.wilcoxon(draw, zero_method="wilcox").pvalue)},
             "passes_primary_criteria": bool(obs > 0 and w.pvalue < 0.05 and pA < 0.05)}
        res["by_vocabulary"][vocab] = v

        log(f"{vocab}: delta {obs:+.5f} ({v['frac_patients_positive']:.0%} of {len(patients)} "
            f"patients positive), Wilcoxon p={w.pvalue:.4g}")
        log(f"   nullA p={pA:.4g} | nullB p={pB:.4g} ({pool.shape[1]} slots x {pool.shape[2]}) | "
            f"nullC p={pC:.4g} ({v['nullC_random_regulon_modules']['percentile_of_observed']:.1f}th pct)")
        log("   sensitivities: " + ", ".join(
            f"{k} {x['mean_delta']:+.5f} (p={x['p_wilcoxon']:.3g}, same sign {x['same_sign_as_primary']})"
            for k, x in sens.items()))
        log(f"   raw scale: {draw.mean():+.5f} (p={v['raw_scale']['p_wilcoxon']:.3g})")

    both = res["by_vocabulary"]
    res["GATE5"] = {
        "both_positive": bool(all(both[v]["mean_delta"] > 0 for v in C.DEFS)),
        "both_wilcoxon_sig": bool(all(both[v]["p_wilcoxon"] < 0.05 for v in C.DEFS)),
        "both_clear_nullA": bool(all(both[v]["nullA_label_permutation"]["p_empirical"] < 0.05
                                     for v in C.DEFS)),
        "both_clear_nullB": bool(all(both[v]["nullB_random_gene_sets"]["p_empirical"] < 0.05
                                     for v in C.DEFS)),
        "both_clear_nullC": bool(all(both[v]["nullC_random_regulon_modules"]["p_empirical"] < 0.05
                                     for v in C.DEFS))}
    res["GATE5"]["PASS"] = bool(res["GATE5"]["both_positive"] and res["GATE5"]["both_wilcoxon_sig"]
                                and res["GATE5"]["both_clear_nullA"])
    res["GATE5"]["PASS_strong_form"] = bool(res["GATE5"]["PASS"] and res["GATE5"]["both_clear_nullC"])

    desc = {}
    for vocab in C.DEFS:
        S = C.load_scores(a.arm, "aucell", vocab, background_regulons)
        sc = C.module_score(S, module_regulons, M["orientation"][vocab])
        g_nbm = sc[[s for s in sc.index if s.upper().startswith("NBM")]]
        g_dx = sc[[pairs[p][0] for p in patients]]
        g_rl = sc[[pairs[p][1] for p in patients]]
        desc[vocab] = {"NBM_mean": float(g_nbm.mean()), "Dx_mean": float(g_dx.mean()),
                       "Rel_mean": float(g_rl.mean()),
                       "p_NBM_vs_Dx": float(stats.mannwhitneyu(g_nbm, g_dx).pvalue),
                       "p_NBM_vs_Rel": float(stats.mannwhitneyu(g_nbm, g_rl).pvalue)}
    res["descriptive_vs_normal_marrow"] = desc
    res["descriptive_note"] = ("NBM are SEPARATE normal-marrow donors, not each patient's own "
                               "healthy compartment — this is NOT a replication of the within-patient "
                               "selectivity analysis and is reported descriptively only.")

    pd.DataFrame(delta_rows).to_csv(OUT / f"per_patient_deltas_{suffix}.csv", index=False)
    (OUT / f"gate5_result_{suffix}.json").write_text(json.dumps(res, indent=2, default=str))

    g = res["GATE5"]
    log("")
    log(f"Acceptance check — positive both: {g['both_positive']} | wilcoxon both: {g['both_wilcoxon_sig']} | "
        f"nullA both: {g['both_clear_nullA']} | nullB both: {g['both_clear_nullB']} | "
        f"nullC both: {g['both_clear_nullC']}")
    log(f"Acceptance check verdict: {'PASS' if g['PASS'] else 'FAIL'}"
        f"{' (STRONG FORM)' if g['PASS_strong_form'] else ''}")
    log(f"wrote gate5_result_{suffix}.json")

if __name__ == "__main__":
    main()
