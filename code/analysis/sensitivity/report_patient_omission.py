"""Summarize patient-omission outputs and make descriptive review figures."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'analysis/out/patient_omission'
summary = json.loads((OUT/'summary.json').read_text())
folds = pd.read_csv(OUT/'nomination_by_omitted_patient.csv')
frequency = pd.read_csv(OUT/'regulon_selection_frequency.csv')
decisions = pd.read_csv(OUT/'all_candidate_fold_decisions.csv')
effects = pd.read_csv(OUT/'fixed_module_patient_influence.csv')
mode = 'recomputed_celltype_eligibility'
primary = folds.loc[folds['mode']==mode].sort_values('n_retained')
original = frequency.loc[(frequency['mode']==mode) & frequency.in_original_55].sort_values(['selected_folds','regulon'])
selected = decisions.loc[(decisions['mode']==mode)&decisions.in_original_55].pivot(
    index='regulon',columns='omitted_patient',values='selected').loc[original.regulon, primary.omitted_patient]
labels = {p:f'P{i+1:02d}' for i,p in enumerate(sorted(primary.omitted_patient))}
pd.DataFrame([{'plot_label':v,'patient':k} for k,v in labels.items()]).to_csv(OUT/'plot_patient_labels.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none'})
fig,(ax,bar)=plt.subplots(1,2,figsize=(11,13),gridspec_kw={'width_ratios':[5,1]},sharey=True)
ax.imshow(selected.to_numpy(int),aspect='auto',cmap=ListedColormap(['#ededed','#287ba8']),vmin=0,vmax=1)
ax.set_yticks(np.arange(55), selected.index,fontsize=7)
ax.set_xticks(np.arange(19),[labels[x] for x in selected.columns],rotation=90)
ax.set_xlabel('Omitted pediatric patient (ordered by retained membership)')
ax.set_title('Original 55: blue = selected again; gray = not selected')
bar.barh(np.arange(55),original.selected_folds,color='#287ba8')
bar.set_xlim(0,20);bar.set_xticks([0,10,19]);bar.set_xlabel('Folds selected')
bar.tick_params(axis='y',left=False,labelleft=False)
fig.suptitle('Patient-omission nomination stability\nFixed networks and activity matrices; all three selection rules recomputed',y=.995)
fig.tight_layout()
for extension in ['png','svg']:
    fig.savefig(OUT/('nomination_stability.'+extension),dpi=220,bbox_inches='tight')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,3.8))
for ax,cohort,metric,title in zip(axes,['pediatric','external'],['equal_celltype_weight','module_delta'],
                                 ['Pediatric: fixed 55, 19 patients','External: fixed 55, 59 pairs']):
    for y,vocab in enumerate(['Diagnosis','Relapse']):
        group=effects.loc[(effects.cohort==cohort)&(effects.metric==metric)&(effects.vocabulary==vocab)]
        full=group.loc[group.omitted_patient=='NONE','mean'].iloc[0]
        omission=group.loc[group.omitted_patient!='NONE','mean'].to_numpy()
        ax.scatter(omission,y+np.linspace(-.14,.14,len(omission)),s=14,c='#287ba8',alpha=.7)
        ax.scatter([full],[y],marker='D',s=50,c='#d35c39',zorder=3)
    ax.axvline(0,color='.5',linewidth=.7);ax.set_yticks([0,1],['Diagnosis targets','Relapse targets'])
    ax.set_xlabel('Mean oriented relapse minus diagnosis activity')
    ax.set_title(title);ax.invert_yaxis()
fig.suptitle('Blue: omit one patient; orange diamond: full cohort\nDescriptive influence ranges, not confidence intervals',fontsize=11)
fig.tight_layout()
for extension in ['png','svg']:
    fig.savefig(OUT/('module_patient_influence.'+extension),dpi=220,bbox_inches='tight')
plt.close(fig)

lines=['# Patient-omission results - 29 September 2026','',
'The overall expression-score signal survives every single-patient omission, while exact nomination membership is substantially more sensitive. Keep the current 55 as the recorded retrospective set; do not silently replace it with a selected subset.','',
'## Nomination stability','',
f"- Full-data reconstruction reproduced all three axes and the exact 55-member set. Numerical comparisons of E2 effects, compartment effects, cell-type counts, stored orientations, external patient deltas and pediatric decomposition passed.",
f"- Across {summary['n_pediatric_folds']} omitted pediatric patients, newly selected sets contain {summary['nomination_size_range'][0]}–{summary['nomination_size_range'][1]} regulons; {summary['retained_original_range'][0]}–{summary['retained_original_range'][1]} of the original 55 remain.",
f"- {summary['original_members_selected_in_every_fold']}/55 survive every omission; {summary['original_members_selected_in_at_least_90pct_folds']}/55 survive at least 90% (at least 18/19) of omissions. The least stable original member survives only {summary['original_members_minimum_selected_folds']}/19.",
f"- Jaccard with the original set ranges from {summary['jaccard_range'][0]:.3f} to {summary['jaccard_range'][1]:.3f}.",
'- Recomputing the at-least-five-patient compartment eligibility sets versus holding those sets fixed gives identical membership results. Omitting the healthy-only patient also leaves membership unchanged.','',
'The 26 members present in every omission are a descriptive stability subset. They are not the historical 26-regulon module or the corrected LMM-significant 26, and are not proposed as a replacement module.','',
'### Least stable original members','', '| Regulon | Selected folds | Ax1 passes | Ax2 passes | Ax3 passes |','|---|---:|---:|---:|---:|']
for row in original.head(12).itertuples():
    lines.append(f'| {row.regulon} | {row.selected_folds}/19 | {row.ax1_pass_folds} | {row.ax2_pass_folds} | {row.ax3_pass_folds} |')
lines += ['', 'Axis 3 was historically evaluated only for axis-1 survivors; these pass counts are not independent failure categories. Complete decisions for all 201 candidates and every fold are retained.', '',
'## Fixed-module patient influence','', '| Cohort / metric | Targets | Full mean | Omission mean range | Largest omission Wilcoxon p |','|---|---|---:|---:|---:|']
for item in summary['fixed_module_effects']:
    lines.append(f"| {item['cohort']} / {item['metric']} | {item['vocabulary']} | {item['full_mean']:.6f} | {item['omitted_mean_min']:.6f}–{item['omitted_mean_max']:.6f} | {item['max_omitted_p']:.3g} |")
lines += ['', 'The overall, equal-cell-type-weighted and within-type module effects retain their positive sign after every omission. The smaller composition component can change sign and its sensitivity p-values are not consistently below 0.05. The decomposition uses six admitted types and diagnosis weights; it is not an edge-rewiring measurement or whole-marrow decomposition.','',
'## Re-selection and orientation sensitivity','']
for vocab,item in summary['reselected_external_and_heldout'].items():
    held=item['heldout_pediatric']
    lines.append(f"- {vocab} targets: external mean after pediatric fold-wise re-selection/orientation ranges {item['external_mean_range'][0]:.6f}–{item['external_mean_range'][1]:.6f}; largest Wilcoxon p={item['external_max_p']:.3g}. The omitted-patient score has mean {held['mean']:.6f}, {held['fraction_positive']:.1%} positive and descriptive Wilcoxon p={held['p_wilcoxon']:.3g}.")
lines += ['', 'These omitted-patient scores avoid using that patient in the downstream membership/sign calculation. They are not independent validation: network inference, target definitions and score panels still use the full discovery cohort. Fold scores also share training patients, so their descriptive p-values should not be read as a confirmatory cross-validation test.','',
'## Implementation and interpretation','',
'- `patient_omission.py` uses existing pseudobulk and external AUCell matrices. Pediatric centering uses the complete stage-specific panel (332/352 columns), whereas external centering uses its recorded 201-regulon background. Subsetting before centering would change the analysis.',
'- Six cell types feed directional agreement and the cell-type rule; seven malignant and eight healthy types feed the compartment comparison. Both stages and compartments are removed together for a patient.',
'- All target sets and networks stay fixed; no CUBIC execution is needed. This does not measure uncertainty in network inference or establish generalization to another cohort.',
'- No original results, module definitions, thresholds or main figures were replaced. Full inputs and software hashes are in `run_manifest.json`; full-cohort parity is in `baseline_parity.json`.',
'- The frozen external orientation is the sign of E2 effects: 30 positive and 25 negative module members under both vocabularies. The funnel’s displayed 31-up/24-down split uses a different direction field derived from the compartment difference. These must not be conflated.',
'', 'Suggested Results wording: “Patient-omission analysis supported the aggregate expression-score association but revealed sensitivity of exact membership: 26 of 55 nominated regulons were retained in all 19 omissions, and each fold retained 39–53 members. External module-score changes remained positive under both fixed-set patient deletion and discovery-fold re-selection. These analyses were conditional on the original inferred networks and target definitions.”','']
(OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
print('\n'.join(lines[:13]))
print('Wrote REPORT.md and four PNG/SVG figure files.')
