'Post-hoc sensitivity: exclude maturity markers from regulons, retain intact axes.\n\nUses the existing adult expression/ranking and scoring implementation without editing it.\nAll 201 regulons are rescored so centering cannot reintroduce excluded marker targets.\nModule membership, pediatric signs, patients, ranking universe and axis scoring stay fixed.\nNo re-selection, new null calibration, or interpretation as a preregistered test is made.\nRun from any directory: python analysis/sensitivity/exclude_maturity_markers.py\n'
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'adult/shared'))
import scoring_common as C
import argparse
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def summarize(y, x):
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    ci = fit.conf_int(alpha=.05)
    return dict(n_patients=len(y), mean_delta=float(y.mean()),
                median_delta=float(np.median(y)), n_positive=int((y > 0).sum()),
                wilcoxon_p=float(stats.wilcoxon(y, zero_method='wilcox').pvalue),
                r_squared=float(fit.rsquared), intercept=float(fit.params[0]),
                intercept_ci_low=float(ci[0, 0]), intercept_ci_high=float(ci[0, 1]),
                intercept_p=float(fit.pvalues[0]), slope=float(fit.params[1]),
                slope_ci_low=float(ci[1, 0]), slope_ci_high=float(ci[1, 1]),
                slope_p=float(fit.pvalues[1]))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'analysis/out/marker_exclusion')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    module_file = ROOT / 'modules/out/corrected_module.json'
    axis_file = ROOT / 'adult/corrected/out/axis_lambo_T3_repair.json'
    gate_file = ROOT / 'adult/corrected/out/gate5_result_T3_repair.json'
    params_file = ROOT / 'adult/corrected/out/scoring_params.json'
    old_axis = json.loads(axis_file.read_text())
    old_gate = json.loads(gate_file.read_text())
    params = json.loads(params_file.read_text())
    M = C.resolve_module(module_file, 'repair')
    regs, bg = M['regulons'], M['background']
    assert len(regs) == 55 and len(bg) == 201
    _, pairs, patients, _ = C.load_pairs('repair')
    assert len(patients) == 59
    assert len({s for p in patients for s in pairs[p]}) == 118
    inputs = [Path(__file__), Path(C.__file__), C.MAIN, C.LAMBO_S3, module_file,
              axis_file, gate_file, params_file, C.out_dir('repair') / 'sample_table.csv']
    inputs += [C.out_dir('repair') / f'scores_aucell_{v}.csv' for v in C.DEFS]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    protocol = dict(status='Post-hoc sensitivity requested by author 2026-11-11',
        created_utc=datetime.now(timezone.utc).isoformat(),
        module='T3_repair: 55 nominated regulons; 201-regulon centering background',
        operation='For each axis separately, remove its marker genes from every regulon target set; rescore and recenter, retaining the intact axis.',
        fixed=['55-regulon membership', 'pediatric orientation by vocabulary', '59 patient pairs',
               'full 23250-gene ranking universe', 'AUCell top 5% cutoff', 'original axis marker lists, including their original multiplicities'],
        empty_target_policy='Stop if any background regulon has zero retained detected targets; do not silently drop regulons.',
        small_target_policy='Retain nonempty regulons; report counts below 5 and 10 detected targets and complete size audit.',
        tests='Two-sided paired Wilcoxon on module changes; ordinary least squares of module change on axis change; conventional two-sided t test and 95% CI for intercept and slope. Nominal p values, no multiplicity correction.',
        interpretation='Sensitivity analysis only; the original null distributions are not reused for the modified module.',
        seed=params['seed'], scoring_parameters=params, input_sha256=hashes,
        environment={'python':sys.version,'executable':sys.executable,
                     'packages':{x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','statsmodels','ctxcore','numba','openpyxl']}})
    (out / 'protocol.json').write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    print('Protocol recorded; loading expression and original rankings.', flush=True)
    E = C.Expr(log=lambda s:print(s, flush=True))
    assert E.TOTAL == params['total_genes'] == 23250
    E.selftest_auc(params['auc_threshold'])
    s3 = pd.read_excel(C.LAMBO_S3, sheet_name='Cell markers', header=0)
    s3.columns = [str(c).strip() for c in s3.columns]
    marks = {c:[g for g in s3[c].dropna().astype(str).str.strip() if g and g.lower() != 'nan'] for c in s3}
    axes = {'A1_lambo_HSC_minus_Monocyte':(marks['HSC'], marks['Monocyte']),
            'A2_lambo_HSCProg_minus_MonoGMP':(marks['HSC'] + marks['Progenitor'], marks['Monocyte'] + marks['GMP'])}
    check, scores, original_modules = {}, {}, {}
    for v in C.DEFS:
        for r in bg:
            assert len(M['targets'][v][r]) == len(set(M['targets'][v][r])), (v,r,'duplicate regulon targets')
        fresh = pd.DataFrame({r:E.auc(M['targets'][v][r], params['auc_threshold']) for r in bg})
        stored = C.load_scores('repair', 'aucell', v, bg).loc[fresh.index, bg]
        diff = float(np.max(np.abs(fresh.values - stored.values)))
        assert diff < 1e-12, (v,'baseline score mismatch',diff)
        scores[v] = fresh
        original_modules[v] = C.module_score(fresh, regs, M['orientation'][v])
        d = C.delta_of(original_modules[v], pairs, patients)
        assert abs(d.mean() - old_gate['by_vocabulary'][v]['mean_delta']) < 1e-12
        assert np.isclose(stats.wilcoxon(d, zero_method='wilcox').pvalue,
                          old_gate['by_vocabulary'][v]['p_wilcoxon'], rtol=1e-9, atol=1e-15)
        check[v] = {'max_absolute_baseline_score_difference':diff, 'original_module_mean_and_wilcoxon_reproduced':True}
        print(f'{v}: original 201-regulon scores and module reproduced, max difference {diff:.3g}', flush=True)
    summary, audit, patient_rows, axis_audit, retention, comparisons = [], [], [], [], [], []
    for axis_name, (pos,neg) in axes.items():
        pos_present = [g for g in pos if g in E.gset]
        neg_present = [g for g in neg if g in E.gset]
        exclude = set(pos) | set(neg)
        ax = E.auc(pos, params['auc_threshold']) - E.auc(neg, params['auc_threshold'])
        x = C.delta_of(ax, pairs, patients)
        assert abs(x.mean()-old_axis['axes'][axis_name]['delta_axis']['mean']) < 1e-12
        for gene in sorted(exclude):
            axis_audit.append(dict(axis=axis_name,gene=gene,in_matrix=gene in E.gset,
                                   positive_multiplicity=pos.count(gene),negative_multiplicity=neg.count(gene)))
        for v in C.DEFS:
            pruned = {}
            these = []
            for r in bg:
                targets = M['targets'][v][r]
                present = set(targets) & E.gset
                keep = [g for g in targets if g not in exclude]
                retained = set(keep) & E.gset
                removed = present & exclude
                assert retained and not (retained & exclude), (axis_name,v,r,'empty target set or exclusion failure')
                pruned[r] = E.auc(keep, params['auc_threshold'])
                row = dict(axis=axis_name,vocabulary=v,regulon=r,is_module_member=r in regs,
                           n_targets_total=len(targets),n_targets_detected=len(present),
                           n_removed_detected=len(removed),n_retained_detected=len(retained),
                           fraction_removed=len(removed)/len(present),removed_genes=';'.join(sorted(removed)))
                audit.append(row);these.append(row)
            P = pd.DataFrame(pruned)
            assert list(P.columns)==bg and P.index.equals(scores[v].index)
            new_module = C.module_score(P, regs, M['orientation'][v])
            y = C.delta_of(new_module, pairs, patients)
            orig_y = C.delta_of(original_modules[v], pairs, patients)
            full_fit = summarize(orig_y,x)
            committed = old_axis['axes'][axis_name]['by_vocabulary'][v]['primary']
            for k,old in [('r_squared','r_squared'),('intercept','adjusted_intercept'),('intercept_p','adjusted_intercept_p'),('slope','slope')]:
                assert np.isclose(full_fit[k],committed[old],rtol=1e-9,atol=1e-12),(axis_name,v,k)
            new_fit = summarize(y,x)
            summary.extend([dict(axis=axis_name,vocabulary=v,construction='original_module_intact_axis',**full_fit),
                            dict(axis=axis_name,vocabulary=v,construction='marker_excluded_module_intact_axis',**new_fit)])
            union = set().union(*(set(M['targets'][v][r]) for r in regs))
            loss = pd.DataFrame(these)
            members = loss[loss.is_module_member]
            retention.append(dict(axis=axis_name,vocabulary=v,n_axis_unique_present=len(exclude & E.gset),
                                  n_axis_positive_entries=len(pos_present),n_axis_negative_entries=len(neg_present),
                                  n_module_union_total=len(union),n_module_union_present=len(union & E.gset),
                                  n_overlap_unique_present=len(union & exclude & E.gset),
                                  n_module_regulons_affected=int((members.n_removed_detected>0).sum()),
                                  min_retained_targets_module=int(members.n_retained_detected.min()),
                                  min_retained_targets_background=int(loss.n_retained_detected.min()),
                                  n_module_below5=int((members.n_retained_detected<5).sum()),
                                  n_module_below10=int((members.n_retained_detected<10).sum()),
                                  n_background_below5=int((loss.n_retained_detected<5).sum()),
                                  n_background_below10=int((loss.n_retained_detected<10).sum()),
                                  median_module_target_fraction_removed=float(members.fraction_removed.median()),
                                  max_module_target_fraction_removed=float(members.fraction_removed.max())))
            comparisons.append(dict(axis=axis_name,vocabulary=v,
                                    retained_mean_change_fraction=float(y.mean()/orig_y.mean()),
                                    original_vs_excluded_patient_change_pearson=float(stats.pearsonr(orig_y,y).statistic)))
            P.to_csv(out/f'scores_marker_excluded_{axis_name}_{v}.csv')
            pd.DataFrame({'original_module_score':original_modules[v], 'marker_excluded_module_score':new_module,
                          'intact_axis_score':ax, 'original_background_mean':scores[v].mean(axis=1),
                          'marker_excluded_background_mean':P.mean(axis=1)}).to_csv(out/f'sample_scores_{axis_name}_{v}.csv',index_label='sample')
            for i,p in enumerate(patients):
                patient_rows.append(dict(axis=axis_name,vocabulary=v,patient=p,diagnosis_sample=pairs[p][0],
                                         relapse_sample=pairs[p][1],intact_axis_delta=float(x[i]),
                                         original_module_delta=float(orig_y[i]),marker_excluded_module_delta=float(y[i])))
            print(f'{axis_name} {v}: mean={y.mean():+.8f}, Wilcoxon p={new_fit["wilcoxon_p"]:.6g}, R2={new_fit["r_squared"]:.6f}, intercept={new_fit["intercept"]:+.8f}, p={new_fit["intercept_p"]:.6g}',flush=True)
    for name,rows in [('summary',summary),('target_retention',retention),('per_regulon_target_audit',audit),
                      ('per_patient_changes',patient_rows),('axis_genes',axis_audit),('module_change_comparison',comparisons)]:
        pd.DataFrame(rows).to_csv(out/f'{name}.csv',index=False)
    for p in inputs:
        assert sha(p)==hashes[str(p.relative_to(ROOT))],('input changed',p)
    check.update(original_axis_and_all_four_regressions_reproduced=True,
                 no_empty_regulons=True,all_exclusions_asserted=True,all_input_hashes_unchanged=True,
                 auc_fast_vs_ctxcore_selftest='bit-identical',n_module_members=55,n_background_regulons=201,n_patients=59)
    results=dict(protocol=protocol,checks=check,summary=summary,target_retention=retention,module_change_comparison=comparisons)
    (out/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False),encoding='utf-8')
    print('Completed:',out,flush=True)

if __name__=='__main__':
    main()
