"""Verify saved marker-exclusion scores and regressions independently of the producer."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'analysis/out/marker_exclusion'
module = json.loads((ROOT/'modules/out/corrected_module.json').read_text())
summary = pd.read_csv(OUT/'summary.csv')
patient = pd.read_csv(OUT/'per_patient_changes.csv')
audit = pd.read_csv(OUT/'per_regulon_target_audit.csv')
genes = pd.read_csv(OUT/'axis_genes.csv')
checks = []
for (axis,vocab), d in patient.groupby(['axis','vocabulary']):
    assert len(d)==59 and d.patient.nunique()==59
    S = pd.read_csv(OUT/f'scores_marker_excluded_{axis}_{vocab}.csv',index_col=0)
    samples = pd.read_csv(OUT/f'sample_scores_{axis}_{vocab}.csv',index_col=0)
    regs = module['regulons']; bg = module['background_regulons']
    assert list(S.columns)==bg and len(regs)==55 and len(bg)==201
    signs = pd.Series(module['orientation'][vocab])[regs]
    reconstructed = S[regs].sub(S.mean(axis=1),axis=0).mul(signs,axis=1).mean(axis=1)
    assert np.allclose(reconstructed,samples.marker_excluded_module_score,atol=1e-14,rtol=0)
    excluded = set(genes.loc[(genes.axis==axis)&genes.in_matrix,'gene'])
    subaudit = audit[(audit.axis==axis)&(audit.vocabulary==vocab)]
    assert len(subaudit)==201
    for r in subaudit.itertuples():
        expected = set(module['targets'][vocab][r.regulon]) & excluded
        observed = set(str(r.removed_genes).split(';')) if pd.notna(r.removed_genes) else set()
        assert observed==expected and len(expected)==r.n_removed_detected
        assert r.n_targets_detected-r.n_removed_detected==r.n_retained_detected>0
    for construction,column in [('original_module_intact_axis','original_module_delta'),
                                ('marker_excluded_module_intact_axis','marker_excluded_module_delta')]:
        row = summary[(summary.axis==axis)&(summary.vocabulary==vocab)&(summary.construction==construction)].iloc[0]
        x=d.intact_axis_delta.to_numpy(); y=d[column].to_numpy()
        X=np.column_stack([np.ones(len(x)),x])
        beta=np.linalg.lstsq(X,y,rcond=None)[0]
        residual=y-X@beta
        df=len(y)-2
        variance=(residual@residual)/df
        se=np.sqrt(np.diag(variance*np.linalg.inv(X.T@X)))
        p=2*stats.t.sf(np.abs(beta/se),df)
        ci=beta[0]+np.array([-1,1])*stats.t.ppf(.975,df)*se[0]
        r2=1-(residual@residual)/np.sum((y-y.mean())**2)
        for observed,expected in [(beta[0],row.intercept),(p[0],row.intercept_p),(r2,row.r_squared),
                                  (ci[0],row.intercept_ci_low),(ci[1],row.intercept_ci_high),
                                  (y.mean(),row.mean_delta),
                                  (stats.wilcoxon(y,zero_method='wilcox').pvalue,row.wilcoxon_p)]:
            assert np.isclose(observed,expected,atol=1e-12,rtol=1e-9),(axis,vocab,construction,observed,expected)
        if construction.startswith('marker_excluded'):
            delta=np.array([reconstructed.loc[b]-reconstructed.loc[a] for a,b in zip(d.diagnosis_sample,d.relapse_sample)])
            assert np.allclose(delta,y,atol=1e-14,rtol=0)
        checks.append(dict(axis=axis,vocabulary=vocab,construction=construction,
                           saved_scores_and_59_patient_pairs_verified=True,
                           independent_numpy_ols_and_t_test_verified=True))
result={'passed':True,'n_regressions_checked':len(checks),'all_804_regulon_exclusion_records_verified':True,
        'centering_and_orientation_rebuilt_from_saved_201_regulon_scores':True,'checks':checks}
(OUT/'independent_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
