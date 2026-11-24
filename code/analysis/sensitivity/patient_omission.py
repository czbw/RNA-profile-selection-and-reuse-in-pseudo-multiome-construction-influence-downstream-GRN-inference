"""Patient-omission stability of the corrected AML nomination and module scores.

Uses existing activity matrices; inferred networks and target definitions stay fixed.
Full-cohort parity is required before any omission results are written.
"""
import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

VOCABULARIES = ('Diagnosis', 'Relapse')
META = ['patient', 'celltype', 'stage', 'driver', 'n_cells']
SIX_TYPES = ['Early.Basophil', 'Early.Erythrocyte', 'GMP', 'HSC', 'Monocytes', 'Progenitor']
AXES = ['ax1_directional', 'ax2_selective', 'ax3_celltype_robust']

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def paired_matrix(path, pairs):
    frame = pd.read_csv(path)
    columns = [c for c in frame if c not in META]
    keys = pd.MultiIndex.from_frame(pairs[['patient', 'celltype']])
    row_keys = pd.MultiIndex.from_frame(frame[['patient', 'celltype']])
    frame = frame.loc[row_keys.isin(keys) & (frame.n_cells >= 20)].copy()
    if frame.duplicated(['patient', 'celltype', 'stage']).any():
        raise ValueError('Duplicate patient/cell-type/stage observations: ' + str(path))
    if set(frame.stage) != set(VOCABULARIES):
        raise ValueError('Unexpected disease stages in ' + str(path))
    values = frame[columns].to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite activity values: ' + str(path))
    # Center on the entire stage-specific panel, before selecting the 201-member universe.
    frame.loc[:, columns] = values - values.mean(axis=1, keepdims=True)
    return frame, columns

def unit_deltas(frame, columns, types, omitted=None):
    mask = frame.celltype.isin(types)
    if omitted is not None:
        mask &= frame.patient != omitted
    part = frame.loc[mask]
    by_stage = {s: part.loc[part.stage == s].set_index(['patient', 'celltype'])[columns]
                for s in VOCABULARIES}
    if not by_stage['Diagnosis'].index.equals(by_stage['Relapse'].index):
        if set(by_stage['Diagnosis'].index) != set(by_stage['Relapse'].index):
            raise ValueError('Unpaired units remain after filtering')
    return by_stage['Relapse'].subtract(by_stage['Diagnosis']).sort_index()

def eligible_types(pairs, omitted=None):
    kept = pairs if omitted is None else pairs.loc[pairs.patient != omitted]
    counts = kept.groupby('celltype').patient.nunique()
    return sorted(counts.index[counts >= 5])

def nominate(data, pairs, background, omitted=None, fixed_types=False):
    types = {arm: eligible_types(pairs[arm], None if fixed_types else omitted)
             for arm in ('malignant', 'healthy')}
    result = pd.DataFrame(index=pd.Index(background, name='regulon'))
    info = {'admitted_malignant': types['malignant'], 'admitted_healthy': types['healthy']}
    for vocab in VOCABULARIES:
        malignant, columns = data['malignant'][vocab]
        healthy, healthy_columns = data['healthy'][vocab]
        if columns != healthy_columns:
            raise ValueError('Malignant/healthy vocabulary mismatch')
        units = unit_deltas(malignant, columns, SIX_TYPES, omitted)[background]
        patients = units.groupby(level='patient').mean()
        mean = patients.mean()
        result['e2_beta_' + vocab] = mean
        per_type = units.groupby(level='celltype').mean()
        counts = (np.sign(per_type).eq(np.sign(mean), axis='columns')).sum()
        counts.loc[mean == 0] = 0
        result['n_celltypes_same_sign_' + vocab] = counts
        result['e2_estimable_' + vocab] = ((patients.shape[0] >= 5) &
            ~pd.Series(np.isclose(patients.to_numpy(), 0).all(axis=0), index=background))
        mal = unit_deltas(malignant, columns, types['malignant'], omitted).groupby(level='patient').mean()
        hea = unit_deltas(healthy, columns, types['healthy'], omitted).groupby(level='patient').mean()
        shared = mal.index.intersection(hea.index)
        if len(shared) < 5:
            raise ValueError('Fewer than five patients in compartment comparison')
        result['beta_malignant_' + vocab] = mal.loc[shared, background].mean()
        result['beta_healthy_' + vocab] = hea.loc[shared, background].mean()
        result['did_estimable_' + vocab] = ~pd.Series(np.isclose(
            (mal.loc[shared, background] - hea.loc[shared, background]).to_numpy(), 0).all(axis=0), index=background)
        result['orientation_' + vocab] = np.sign(mean).astype(int)
        info['e2_patients_' + vocab] = list(patients.index)
        info['compartment_patients_' + vocab] = list(shared)
    result[AXES[0]] = ((np.sign(result.e2_beta_Diagnosis) == np.sign(result.e2_beta_Relapse)) &
                       result.e2_estimable_Diagnosis & result.e2_estimable_Relapse)
    same_compartments = ((np.sign(result.beta_malignant_Diagnosis) == np.sign(result.beta_healthy_Diagnosis)) &
                        (np.sign(result.beta_malignant_Relapse) == np.sign(result.beta_healthy_Relapse)))
    result[AXES[1]] = ~same_compartments & result.did_estimable_Diagnosis & result.did_estimable_Relapse
    # The historical decomposition table contains only directionally concordant regulons.
    result[AXES[2]] = (result[AXES[0]] & (result.n_celltypes_same_sign_Diagnosis >= 4) &
                      (result.n_celltypes_same_sign_Relapse >= 4))
    result['selected'] = result[AXES].all(axis=1)
    return result, info

def module_deltas(patient_matrix, members, orientation):
    if not members:
        return pd.Series(np.nan, index=patient_matrix.index)
    return patient_matrix[members].mul(pd.Series(orientation)[members], axis='columns').mean(axis=1)

def effect(values):
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        return {'n': len(values), 'mean': None, 'median': None, 'fraction_positive': None, 'p_wilcoxon': None}
    p = 1.0 if np.all(values == 0) else float(stats.wilcoxon(values, zero_method='wilcox').pvalue)
    return dict(n=len(values), mean=float(values.mean()), median=float(np.median(values)),
                fraction_positive=float((values > 0).mean()), p_wilcoxon=p)

def external_deltas(root, vocab, background):
    directory = root / 'adult/corrected/out'
    samples = pd.read_csv(directory / 'sample_table.csv')
    aml = samples.loc[samples.group == 'AML'].copy()
    aml['is_relapse'] = aml.stage.str.startswith('Rel')
    scores = pd.read_csv(directory / ('scores_aucell_' + vocab + '.csv'), index_col=0)[background]
    if scores.index.duplicated().any():
        raise ValueError('Duplicate external sample identifiers')
    centered = scores.subtract(scores.mean(axis=1), axis='index')
    rows, labels, selections = [], [], []
    for patient, group in aml.groupby('patient'):
        diagnosis = sorted(group.loc[~group.is_relapse, 'sample'])
        relapse = sorted(group.loc[group.is_relapse, 'sample'])
        if diagnosis and relapse:
            labels.append(patient)
            rows.append((centered.loc[relapse[0]] - centered.loc[diagnosis[0]]).to_numpy())
            selections.append({'patient': patient, 'diagnosis': diagnosis[0], 'relapse': relapse[0]})
    return pd.DataFrame(rows, index=pd.Index(labels, name='patient'), columns=background), selections

def decompose_module(frame, columns, members, orientation):
    selected = frame.loc[frame.celltype.isin(SIX_TYPES)].copy()
    selected['module_score'] = selected[members].mul(pd.Series(orientation)[members], axis='columns').mean(axis=1)
    rows = []
    for patient, group in selected.groupby('patient'):
        dx = group.loc[group.stage == 'Diagnosis'].set_index('celltype')
        relapse = group.loc[group.stage == 'Relapse'].set_index('celltype')
        types = dx.index.intersection(relapse.index)
        if len(types) < 2:
            continue
        wd = dx.loc[types, 'n_cells'] / dx.loc[types, 'n_cells'].sum()
        wr = relapse.loc[types, 'n_cells'] / relapse.loc[types, 'n_cells'].sum()
        md, mr = dx.loc[types, 'module_score'], relapse.loc[types, 'module_score']
        total = float((wr * mr).sum() - (wd * md).sum())
        within = float((wd * (mr - md)).sum())
        rows.append(dict(patient=patient, total=total, within_type=within, composition=total-within))
    return pd.DataFrame(rows).set_index('patient')

def parity_checks(root, baseline, module, internal, external, decomposition):
    checks = []
    def numeric(name, actual, expected):
        actual, expected = np.asarray(actual, float), np.asarray(expected, float)
        passed = actual.shape == expected.shape and np.allclose(actual, expected, rtol=1e-9, atol=2e-13)
        checks.append(dict(check=name, passed=bool(passed), max_abs_difference=float(np.max(np.abs(actual-expected)))))
    funnel = pd.read_csv(root/'longitudinal/corrected/out/nomination_funnel.csv', index_col=0).loc[baseline.index]
    for axis in AXES:
        checks.append(dict(check=axis, passed=bool(baseline[axis].equals(funnel[axis]))))
    checks.append(dict(check='55-member identity', passed=set(baseline.index[baseline.selected]) == set(module['regulons'])))
    e2 = pd.read_csv(root/'longitudinal/corrected/out/e2_paired_wilcoxon_reconciliation.csv', index_col=0).loc[baseline.index]
    dec = pd.read_csv(root/'longitudinal/corrected/out/decomposition_reconciliation.csv', index_col=0)
    old_external = pd.read_csv(root/'adult/corrected/out/per_patient_deltas_T3_repair.csv')
    old_dec = pd.read_csv(root/'analysis/out/f6_decomposition_module_by_patient.csv')
    for vocab in VOCABULARIES:
        numeric('E2 beta '+vocab, baseline['e2_beta_'+vocab], e2['beta_'+vocab])
        for kind in ['malignant', 'healthy']:
            key = 'beta_'+kind+'_'+vocab
            numeric(key, baseline[key], funnel[key])
        numeric('cell-type signs '+vocab, baseline.loc[dec.index, 'n_celltypes_same_sign_'+vocab], dec['n_celltypes_same_sign_'+vocab])
        numeric('stored orientation '+vocab, baseline['orientation_'+vocab], [module['orientation'][vocab][r] for r in baseline.index])
        old = old_external.loc[old_external.vocabulary == vocab].set_index('patient')
        numeric('external patient deltas '+vocab, external[vocab].loc[old.index], old.delta)
        old = old_dec.loc[(old_dec.vocabulary == vocab) & (old_dec.scale == 'oriented_offsetcentred')].set_index('patient')
        for current, previous in [('total','naive'), ('within_type','rewiring'), ('composition','composition')]:
            numeric('pediatric '+current+' '+vocab, decomposition[vocab].loc[old.index, current], old[previous])
    return checks

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--out', type=Path, default=Path('analysis/out/patient_omission'))
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.out if args.out.is_absolute() else root/args.out
    output.mkdir(parents=True, exist_ok=True)
    input_paths = []
    def source(relative):
        path = root/relative
        input_paths.append(path)
        return path
    module = json.loads(source('modules/out/corrected_module.json').read_text(encoding='utf-8'))
    members, background = module['regulons'], module['background_regulons']
    assert len(members) == 55 and len(background) == 201
    directories = {'malignant':'longitudinal/corrected/out', 'healthy':'longitudinal/healthy_corrected/out'}
    pairs, data = {}, {}
    for arm, directory in directories.items():
        healthy = '_healthy' if arm == 'healthy' else ''
        pairs[arm] = pd.read_csv(source(directory+'/paired_units'+healthy+'.csv'))
        assert not pairs[arm].duplicated(['patient','celltype']).any()
        data[arm] = {v: paired_matrix(source(directory+'/pseudobulk'+healthy+'_def'+v+'.csv'), pairs[arm])
                     for v in VOCABULARIES}
    baseline, design = nominate(data, pairs, background)
    pediatric_patients = design['e2_patients_Diagnosis']
    internal_matrices, external_matrices, internal, external, decomposition = {}, {}, {}, {}, {}
    for vocab in VOCABULARIES:
        frame, columns = data['malignant'][vocab]
        internal_matrices[vocab] = unit_deltas(frame, columns, SIX_TYPES).groupby(level='patient').mean()[background]
        external_matrices[vocab], sample_pairs = external_deltas(root, vocab, background)
        internal[vocab] = module_deltas(internal_matrices[vocab], members, module['orientation'][vocab])
        external[vocab] = module_deltas(external_matrices[vocab], members, module['orientation'][vocab])
        decomposition[vocab] = decompose_module(frame, columns, members, module['orientation'][vocab])
        source('adult/corrected/out/scores_aucell_'+vocab+'.csv')
    for name in ['longitudinal/corrected/out/nomination_funnel.csv',
                 'longitudinal/corrected/out/e2_paired_wilcoxon_reconciliation.csv',
                 'longitudinal/corrected/out/decomposition_reconciliation.csv',
                 'longitudinal/healthy_corrected/out/did_reconciliation.csv',
                 'adult/corrected/out/sample_table.csv',
                 'adult/corrected/out/per_patient_deltas_T3_repair.csv',
                 'analysis/out/f6_decomposition_module_by_patient.csv']:
        source(name)
    checks = parity_checks(root, baseline, module, internal, external, decomposition)
    (output/'baseline_parity.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    if not all(x['passed'] for x in checks):
        raise RuntimeError('Baseline parity failed; omission analysis stopped. See baseline_parity.json')
    print('Full-cohort parity passed: 55 members and all three axes.', flush=True)
    fold_rows, membership_rows, fold_scores = [], [], []
    full_set = set(members)
    for fixed in [False, True]:
        mode = 'fixed_celltype_eligibility' if fixed else 'recomputed_celltype_eligibility'
        for omitted in pediatric_patients:
            fold, fold_design = nominate(data, pairs, background, omitted, fixed_types=fixed)
            selected = list(fold.index[fold.selected])
            selected_set = set(selected)
            row = dict(mode=mode, omitted_patient=omitted, n_selected=len(selected),
                       n_retained=len(selected_set & full_set), n_lost=len(full_set-selected_set),
                       n_gained=len(selected_set-full_set), jaccard=len(selected_set & full_set)/len(selected_set | full_set),
                       lost=';'.join(sorted(full_set-selected_set)), gained=';'.join(sorted(selected_set-full_set)),
                       admitted_malignant=';'.join(fold_design['admitted_malignant']),
                       admitted_healthy=';'.join(fold_design['admitted_healthy']))
            fold_rows.append(row)
            for regulon in background:
                item = {'mode':mode, 'omitted_patient':omitted, 'regulon':regulon, 'in_original_55':regulon in full_set}
                item.update(fold.loc[regulon].to_dict())
                membership_rows.append(item)
            if not fixed:
                for vocab in VOCABULARIES:
                    signs = fold['orientation_'+vocab].to_dict()
                    held = module_deltas(internal_matrices[vocab].loc[[omitted]], selected, signs).iloc[0]
                    transferred = module_deltas(external_matrices[vocab], selected, signs)
                    fold_scores.append(dict(omitted_patient=omitted, vocabulary=vocab, n_selected=len(selected),
                        heldout_pediatric_delta=float(held), external_mean=float(transferred.mean()),
                        external_p_wilcoxon=effect(transferred)['p_wilcoxon'],
                        n_original_orientation_changes=sum(signs[r] != module['orientation'][vocab][r] for r in members)))
    folds = pd.DataFrame(fold_rows)
    memberships = pd.DataFrame(membership_rows)
    pd.DataFrame(fold_scores).to_csv(output/'reselected_module_scores.csv', index=False)
    folds.to_csv(output/'nomination_by_omitted_patient.csv', index=False)
    memberships.to_csv(output/'all_candidate_fold_decisions.csv', index=False)
    frequency = memberships.groupby(['mode','regulon']).agg(
        in_original_55=('in_original_55','first'), selected_folds=('selected','sum'),
        total_folds=('selected','size'), ax1_pass_folds=(AXES[0],'sum'),
        ax2_pass_folds=(AXES[1],'sum'), ax3_pass_folds=(AXES[2],'sum')).reset_index()
    frequency['selection_fraction'] = frequency.selected_folds / frequency.total_folds
    frequency.to_csv(output/'regulon_selection_frequency.csv', index=False)
    pd.DataFrame(sample_pairs).to_csv(output/'external_sample_pairs.csv', index=False)

    effects, deltas = [], []
    def record_effects(cohort, metric, vocab, values):
        full = effect(values)
        effects.append(dict(cohort=cohort, metric=metric, vocabulary=vocab, omitted_patient='NONE', **full))
        for patient, value in values.items():
            deltas.append(dict(cohort=cohort, metric=metric, vocabulary=vocab, patient=patient, delta=float(value)))
            result = effect(values.drop(patient))
            effects.append(dict(cohort=cohort, metric=metric, vocabulary=vocab, omitted_patient=patient,
                                change_from_full=result['mean']-full['mean'], **result))
    for vocab in VOCABULARIES:
        record_effects('pediatric', 'equal_celltype_weight', vocab, internal[vocab])
        record_effects('external', 'module_delta', vocab, external[vocab])
        for metric in ['total', 'within_type', 'composition']:
            record_effects('pediatric', metric, vocab, decomposition[vocab][metric])
    effects = pd.DataFrame(effects)
    effects.to_csv(output/'fixed_module_patient_influence.csv', index=False)
    pd.DataFrame(deltas).to_csv(output/'fixed_module_patient_deltas.csv', index=False)
    baseline.to_csv(output/'baseline_rebuilt_nomination.csv')

    mode = 'recomputed_celltype_eligibility'
    primary = folds.loc[folds['mode'] == mode]
    original = frequency.loc[(frequency['mode'] == mode) & frequency.in_original_55]
    extra_healthy = sorted(set(pairs['healthy'].patient) - set(pediatric_patients))
    # Healthy-only observations can affect eligibility; report their deletion separately.
    healthy_only_checks = []
    for patient in extra_healthy:
        fold, _ = nominate(data, pairs, background, patient)
        healthy_only_checks.append(dict(patient=patient, membership_identical=bool(fold.selected.equals(baseline.selected))))
    summary = dict(status='Exploratory patient-omission analysis conditional on fixed networks and activity matrices',
        n_pediatric_folds=len(pediatric_patients), n_external_pairs=len(external['Diagnosis']),
        baseline_module_size=len(members), nomination_size_range=[int(primary.n_selected.min()),int(primary.n_selected.max())],
        retained_original_range=[int(primary.n_retained.min()),int(primary.n_retained.max())],
        jaccard_range=[float(primary.jaccard.min()),float(primary.jaccard.max())],
        original_members_selected_in_every_fold=int((original.selected_folds == len(pediatric_patients)).sum()),
        original_members_selected_in_at_least_90pct_folds=int((original.selection_fraction >= .9).sum()),
        original_members_minimum_selected_folds=int(original.selected_folds.min()),
        healthy_only_omission_checks=healthy_only_checks, design=design,
        eligibility_sensitivity_membership_identical=bool(np.array_equal(
            memberships.loc[memberships['mode']==mode,'selected'].to_numpy(),
            memberships.loc[memberships['mode']=='fixed_celltype_eligibility','selected'].to_numpy())),
        fixed_module_effects=[], reselected_external_and_heldout={})
    for (cohort, metric, vocab), group in effects.groupby(['cohort','metric','vocabulary']):
        full = group.loc[group.omitted_patient == 'NONE'].iloc[0]
        omitted = group.loc[group.omitted_patient != 'NONE']
        summary['fixed_module_effects'].append(dict(cohort=cohort, metric=metric, vocabulary=vocab,
            full_mean=float(full['mean']), omitted_mean_min=float(omitted['mean'].min()),
            omitted_mean_max=float(omitted['mean'].max()), max_omitted_p=float(omitted.p_wilcoxon.max()),
            all_omissions_same_mean_sign=bool((np.sign(omitted['mean']) == np.sign(full['mean'])).all())))
    scores = pd.DataFrame(fold_scores)
    for vocab, group in scores.groupby('vocabulary'):
        summary['reselected_external_and_heldout'][vocab] = dict(
            external_mean_range=[float(group.external_mean.min()),float(group.external_mean.max())],
            external_max_p=float(group.external_p_wilcoxon.max()),
            heldout_pediatric=effect(group.heldout_pediatric_delta),
            caveat='Training-fold membership/orientation only; fixed networks still use all discovery patients.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
        packages={'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
        script_sha256=sha256(Path(__file__)),
        inputs=[dict(path=p.relative_to(root).as_posix(),sha256=sha256(p),bytes=p.stat().st_size)
                for p in sorted(set(input_paths))],
        no_network_reinference=True, no_original_output_overwritten=True)
    (output/'run_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['design','fixed_module_effects']}, indent=2))

if __name__ == '__main__':
    main()
