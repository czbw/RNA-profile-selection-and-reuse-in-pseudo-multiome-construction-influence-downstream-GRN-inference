"""Assemble supplementary tables from preserved study outputs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--workbook-spec', action='store_true', help='Also write the table specification for Excel reading copies.')
args = parser.parse_args()
source = args.source.resolve()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
manifest = json.loads((source/'manifest.json').read_text())
lookup = {old:r['path'] for r in manifest for old in r['previous_supplementary_paths']}
lookup['matched_report_2008']='outputs/benchmarks/bmmc/out_matched/b02_arm_report.json'
for rec in manifest:
    assert hashlib.sha256((source/rec['path']).read_bytes()).hexdigest()==rec['sha256'],rec['path']

titles = {1:'Published link and metacell comparisons',2:'Activity controls and native-pairing benchmarks',3:'Matched-reuse assignments and activity',4:'Network agreement and pAML nomination',5:'Adult expression scores and functional annotation',6:'Pediatric decomposition and sensitivity analyses'}
catalog=[]
columns=[]
books={str(n):{'title':t,'sheets':[]} for n,t in titles.items()}

def readj(key):
    return json.loads((source/lookup[key]).read_text(encoding='utf-8-sig'))

def readcsv(key):
    with (source/lookup[key]).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def typed(v):
    if v is None or v=='': return None
    if not isinstance(v,str): return v
    if v in ('True','False'): return v=='True'
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except ValueError:
        return v

def table(n,name,rows,keys,description,selection,fields=None,sheet=None):
    assert rows, name
    headers=list(dict.fromkeys(k for row in rows for k in row))
    rel=f'Supplementary_Data_{n}/{name}.csv'
    p=output/rel
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(rows)
    catalog.append({'table':rel,'description':description,'source_files':'; '.join(lookup[k] for k in keys),'selection':selection,'rows':len(rows),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    for h in headers:
        columns.append({'table':rel,'column':h,'source_field_or_definition':(fields or {}).get(h,h)})
    if sheet:
        books[str(n)]['sheets'].append({'name':sheet,'table':rel,'headers':headers,'rows':[[typed(r.get(h)) if h not in ('patient','omitted_patient','sample','stratum','regulon','TF','member') else r.get(h) for h in headers] for r in rows]})

def csvtable(n,name,key,description,cols=None,where=None,selection='All records.',rename=None,values=None,sheet=None):
    rows=[]; rename=rename or {}; values=values or {}
    for i,r in enumerate(readcsv(key),2):
        if where and not where(r): continue
        rr={rename.get(k,k):values.get(k,{}).get(r[k],r[k]) for k in (cols or r.keys())}
        rr['source_row']=i
        rows.append(rr)
    fields={rename.get(k,k):k for k in (cols or readcsv(key)[0])}
    fields['source_row']='1-based CSV record number, including the header.'
    table(n,name,rows,[key],description,selection,fields,sheet)

def leaves(x,path=''):
    if isinstance(x,dict):
        for k,v in x.items():
            yield from leaves(v,path+'/'+str(k).replace('~','~0').replace('/','~1'))
    elif isinstance(x,list):
        for i,v in enumerate(x): yield from leaves(v,path+'/'+str(i))
    elif isinstance(x,(int,float,bool)) and (not isinstance(x,float) or math.isfinite(x)):
        yield path,x

def numeric(n,name,key,branches,description):
    j=readj(key);rows=[]
    for branch in branches:
        x=j
        for k in branch: x=x[k]
        for ptr,v in leaves(x,'/'+ '/'.join(branch) if branch else ''):
            rows.append({'source_json_pointer':ptr,'value':v})
    table(n,name,rows,[key],description,'Numeric fields in the listed JSON branches. Text annotations are omitted.')

armnames={'link':'published link','pseudomultiome':'metacell','repair':'corrected'}
d1='Supplementary_Data_1/'
csvtable(1,'membership',d1+'f1_set_members.csv','TF and signed eRegulon membership (Figure 1C).',
    ['stage','level','member','in_link','in_pseudomultiome'],lambda r:r['in_link']=='True' or r['in_pseudomultiome']=='True',
    'Members present in published link or metacell networks.',{'in_link':'in_published_link','in_pseudomultiome':'in_metacell'})
csvtable(1,'shared_regulon_agreement',d1+'f1_activity_vs_target.csv','Activity and target agreement for the same shared eRegulons (Figure 1D).',
    rename={'n_genes_A':'n_genes_published_link','n_genes_B':'n_genes_metacell'},sheet='Shared eRegulons')
csvtable(1,'shared_tf_targets',d1+'f1_shared_tf_jaccard.csv','Target overlap after pooling signed eRegulons within each TF.',
    where=lambda r:r['arm_A']=='link' and r['arm_B']=='pseudomultiome',selection='Published link versus metacell only.',values={'arm_A':armnames,'arm_B':armnames})
j=readj(d1+'f1_sets_and_targets.json'); rows=[]
for key,v in j['inputs'].items():
    stage,arm=key.split('_',1)
    if arm not in ('link','pseudomultiome'):continue
    rows.append({'stage':stage,'construction':armnames[arm],**{k:v[k] for k in ['n_regulons','n_TFs','rows','n_genes_union','n_regions_union']}})
table(1,'network_counts',rows,[d1+'f1_sets_and_targets.json'],'Network sizes underlying Figure 1C.','inputs: published link and metacell at each stage.',{'rows':'inputs/<stage>_<arm>/rows; TF-region-gene triplets.'},'Network counts')
books['1']['sheets'].reverse()

d2='Supplementary_Data_2/'
key=d2+'activity_controls/f1_random_geneset_control.csv'
csvtable(2,'random_gene_controls',key,'Expression- and size-matched random gene controls for Figure 1.',
    ['stage','regulon','TF','n_genes_link','n_genes_metacell','r_observed_stored','r_observed_rebuilt','N1_median','N1_q95','N1_p','N1_q'],
    rename={'n_genes_link':'n_genes_published_link','N1_median':'random_median','N1_q95':'random_95th_percentile','N1_p':'empirical_p','N1_q':'BH_q'})
key=d2+'activity_controls/f1_random_geneset_control.json';j=readj(key)
rows=[{'stage':s,'n_shared_regulons':v['n_shared_eregulons'],'n_cell_types':v['n_shared_celltypes'],'observed_median_r':v['observed_stored']['median'],'random_median_r':v['pooled_null']['N1']['median'],'n_above_own_random_95th_percentile':v['n_eregulons_above_own_null_q95']['N1'],'random_sets_per_regulon':j['draws']} for s,v in j['stages'].items()]
table(2,'random_control_summary',rows,[key],'Summary of the random gene controls cited in Results.','stages/*; N1 uses the full eligible gene universe.','' or None,'Random controls')
for dataset in ['PBMC','BMMC']:
    key=d2+dataset+('/trackb_rungs_pbmc10k_all.json' if dataset=='PBMC' else '/trackb_rungs_bmmc_all.json')
    j=readj(key);rows=[]
    for arm,v in j['L3_L4_L5_eRegulons'].items():
        if arm=='_truth' or 'n_triplets' not in v:continue
        l1=j['L1_region_to_gene'].get(arm,{})
        l2=j['L2_tf_to_gene'].get(arm,{})
        rows.append({'arm':arm,'L1_top_pair_jaccard':l1.get('topN_jaccard_vs_truth'),'L1_rho_spearman':l1.get('spearman_rho'),'L1_effective_pairs':l1.get('n_eff_rho'),'L2_top_edge_jaccard':l2.get('topN',{}).get('jaccard'),'L3_target_gene_jaccard':v.get('L3_target_sets',{}).get('median_jaccard_genes'),'L3_target_region_jaccard':v.get('L3_region_sets',{}).get('median_jaccard_regions'),'L4_regulon_jaccard':v.get('L4_eRegulon_recovery',{}).get('jaccard'),'L5_TF_jaccard':v.get('L5_TF_recovery',{}).get('jaccard'),'n_triplets':v['n_triplets'],'reference_triplets':j['L3_L4_L5_eRegulons']['_truth']['n_triplets'],'triplet_count_ratio':v['n_triplets']/j['L3_L4_L5_eRegulons']['_truth']['n_triplets']})
    table(2,dataset.lower()+'_network_agreement',rows,[key],f'{dataset} native-reference network comparisons (Figure 2A,B).','All completed arms; undefined or uncomputed metrics remain blank.',{'triplet_count_ratio':'n_triplets / reference_triplets.'},dataset+' networks')
    numeric(2,dataset.lower()+'_network_details',key,[[k] for k in ['L1_region_to_gene','L2_tf_to_gene','L3_L4_L5_eRegulons']],f'Counts and additional network metrics for {dataset}.')
for f,name in [('PBMC/trackA_L6_10k.json','pbmc_cell_activity'),('PBMC/trackA_L6_10k_subcluster_arm.json','pbmc_stratified_cell_activity')]:
    numeric(2,name,d2+f,[[]],'Per-cell activity agreement (Figure 2C).')
numeric(2,'bmmc_cell_activity',d2+'BMMC/b06_trackA_L6_bmmc.json',[['arms']],'Per-cell activity agreement (Figure 2C).')
numeric(2,'bmmc_assignment_accuracy',d2+'BMMC/b05_pairing_error_bmmc.json',[['arms']],'Assignment accuracy and label agreement in BMMC.')
numeric(2,'rna_neighbor_comparator',d2+'activity_controls/l6_ceiling_test.json',[['PBMC'],['BMMC']],'RNA nearest-neighbor comparator (Figure 2C).')
for f,name in [('b06_trackA_L6_per_celltype_bmmc.csv','bmmc_cell_type_activity'),('b06_trackA_L6_per_regulon_bmmc.csv','bmmc_regulon_activity'),('b07_stat3_per_regulon_bmmc.csv','bmmc_group_mean_activity')]:
    csvtable(2,name,d2+'BMMC/'+f,'Underlying cell-type or regulon activity comparisons.')

d3='Supplementary_Data_3/'
key=d3+'construction/b04e_cap14_pair_verify.json';j=readj(key)
rows=[]
for arm,v in j['arms'].items():
    rows.append({'arm':arm,'partner_selection':'random' if arm.startswith('rand') else 'inferred','n_assignments':v['n_assignments'],'distinct_RNA_profiles':v['distinct_rna_cells'],'reuse_ratio':v['reuse_mean_per_used_rna_cell'],'maximum_reuse':v['reuse_max'],'exact_cell_recovery':v['exact_cell_accuracy'],'native_cells_recovered':v['self_hits'],'total_AUCell_variance_ratio':v['retained_variance_total_frac_of_truth'],'within_cell_type_AUCell_variance_ratio':v['retained_variance_within_cell_type_frac_of_truth']})
rows.sort(key=lambda r:(r['reuse_ratio'],r['partner_selection']))
table(3,'matched_reuse_summary',rows,[key],'Assignment counts and activity variance for all six matched-reuse arms (Figure 3A,D).','arms: all six entries. Ratios use the native-pairing value as denominator.',{'distinct_RNA_profiles':'distinct_rna_cells','reuse_ratio':'reuse_mean_per_used_rna_cell','total_AUCell_variance_ratio':'retained_variance_total_frac_of_truth','within_cell_type_AUCell_variance_ratio':'retained_variance_within_cell_type_frac_of_truth'},'Assignments and variance')
checks=[]
keys=[d3+'construction/b02_matched_cap14_report.json','matched_report_2008']
for key,arm in zip(keys,['rand_celltype_donor_matched_cap14','rand_celltype_donor_matched']):
    a=readj(key)['arms'][arm]
    for stratum,v in a['per_stratum'].items():
        checks.append({'reuse_ratio':a['reuse_mean_per_used_rna_cell'],'stratum':stratum,'n_ATAC_cells':v['n_s'],'inferred_distinct_RNA_profiles':v['d_s_target'],'random_distinct_RNA_profiles':v['d_s_achieved'],'inferred_stratum_reuse':v['reuse_target'],'random_stratum_reuse':v['reuse_achieved'],'multiplicity_KS_distance':v['ks_multiplicity']})
assert len(checks)==380 and all(r['multiplicity_KS_distance']==0 and r['inferred_distinct_RNA_profiles']==r['random_distinct_RNA_profiles'] for r in checks)
table(3,'multiplicity_checks',checks,keys,'Per-stratum equality checks for the two reuse-matched comparisons.','All 190 strata per reused level. Reuse 1.000 assigns every RNA profile once by construction.',{'multiplicity_KS_distance':'per_stratum/<stratum>/ks_multiplicity; zero denotes equal empirical multiplicity distributions.'})
csvtable(3,'stratum_census',d3+'construction/b00_stratum_census.csv','The 190 donor-cell-type strata used in matched-reuse comparisons.',where=lambda r:r['definition']=='celltype_x_donor',selection='definition == celltype_x_donor.')
csvtable(3,'donor_activity_summary',d3+'stratum_mean_activity/f3_stratum_mean_summary.csv','Donor-stratum mean activity agreement (Figure 3C), with eligibility sensitivity results.',cols=['level','role','arm','rho_strata','n_celltypes','n_strata','rho_strata_q25','rho_strata_q75','n_grid_rho','rho_strata_c10_d5','rho_strata_c50_d5','rho_strata_c20_d3','rho_strata_c20_d7'],rename={'role':'partner_selection'},values={'role':{'pair':'inferred'}},sheet='Donor activity')
csvtable(3,'donor_activity_by_regulon',d3+'stratum_mean_activity/f3_stratum_mean_grid.csv','Per-cell-type-regulon donor-stratum correlations for the primary eligibility rule.',rename={'role':'partner_selection'},values={'role':{'pair':'inferred'}})
key=d3+'construction/b04d_capacity_sweep.json'
rows=[{k:v[k] for k in ['cap','cap_frac','n_assignments','distinct_rna_cells','reuse_mean_per_used_rna_cell','self_hits','exact_cell_accuracy']} for v in readj(key)['sweep_points']]
table(3,'capacity_sweep',rows,[key],'Capacity sweep used to choose the intermediate reuse level.','All sweep points. The source report was assembled from two runs; no new matching was performed.')

d4='Supplementary_Data_4/'
csvtable(4,'matched_reuse_networks',d4+'matched_reuse/figS_reuse_native_agreement.csv','Absolute agreement with the native network (Figure 3B and Supplementary Figure S3).',['metric','rule','level','reuse','arm','agreement'],sheet='Matched-reuse networks')
key=d4+'pAML_construction/aml_reuse_composite_key.json';j=readj(key)
rows=[]
for stage,v in j['stages'].items():
    for arm in ['published link','corrected']:
        old=arm=='published link'
        rows.append({'stage':stage,'construction':arm,'n_assignments':v['n_in' if old else 'n_out'],'distinct_RNA_profiles':v['old_distinct_donors' if old else 'distinct_donors'],'reuse_ratio':v['old_reuse' if old else 'achieved_reuse'],'corrected_scarcity_floor':None if old else v['scarcity_floor']})
table(4,'paml_reuse',rows,[key],'pAML RNA reuse using sample-plus-barcode identities (Figure 4).','stages: composite-key counts only. Stored donor fields count RNA profiles, not patients.',None,'pAML reuse')
network_rows=[]
for name,v in readj(d1+'f1_sets_and_targets.json')['inputs'].items():
    stage,arm=name.split('_',1)
    if arm in ['link','repair']:
        network_rows.append({'stage':stage,'construction':armnames[arm],'n_regulons':v['n_regulons'],'n_TFs':v['n_TFs'],'n_triplets':v['rows']})
table(4,'paml_network_counts',network_rows,[d1+'f1_sets_and_targets.json'],'Published link and corrected network sizes (Figure 4).','inputs: link and repair constructions.',None,'pAML network counts')
funnel_rows=[]
for arm in ['corrected','published_link']:
    key=d4+arm+'/nomination_funnel.csv'
    wanted=['regulon','TF','beta_malignant_Diagnosis','beta_malignant_Relapse','beta_healthy_Diagnosis','beta_healthy_Relapse','did_Diagnosis','did_Relapse','p_did_Diagnosis','p_did_Relapse','q_did_Diagnosis','q_did_Relapse','n_celltypes_same_sign_Diagnosis','n_celltypes_same_sign_Relapse','ax1_directional','ax2_selective','ax3_celltype_robust','direction']
    csvtable(4,arm+'_nomination',key,'All shared candidates and the three further cumulative nomination criteria (Figure 4).',wanted,rename={'ax1_directional':'criterion_2_direction','ax2_selective':'criterion_3_selectivity','ax3_celltype_robust':'criterion_4_cell_type_agreement'})
    selected=readcsv(key)
    for criterion,label in [(None,'Shared signed identity'),('ax1_directional','Consistent direction'),('ax2_selective','Malignant-versus-healthy selectivity'),('ax3_celltype_robust','Cell-type agreement')]:
        if criterion:selected=[r for r in selected if r[criterion]=='True']
        funnel_rows.append({'construction':arm.replace('_',' '),'cumulative_criterion':label,'n_regulons':len(selected),'n_TFs':len({r['TF'] for r in selected})})
    for f,name in [('e2_paired_wilcoxon_fits.csv','paired_tests'),('lmm_fits_by_definition.csv','mixed_models'),('paired_units.csv','paired_units')]:
        csvtable(4,arm+'_'+name,d4+arm+'/'+f,'pAML patient-cell-type units and longitudinal statistical results.')
    numeric(4,arm+'_calibration',d4+arm+'/s3_recalibrated.json',[[]],'Mixed-model calibration diagnostics; numeric results only.')
table(4,'nomination_counts',funnel_rows,[d4+a+'/nomination_funnel.csv' for a in ['corrected','published_link']],'Cumulative candidate counts (Figure 4).','All shared rows, followed by cumulative intersection of ax1_directional, ax2_selective and ax3_celltype_robust. Distinct TFs are counted within each retained set.',None,'Nomination counts')
for f,name in [('a00_stratum_census.csv','paml_stratum_census'),('a01_per_stratum.csv','paml_stratum_assignment')]:
    csvtable(4,name,d4+'pAML_construction/'+f,'pAML eligible RNA and ATAC counts by patient, stage and cell-type.')

d5='Supplementary_Data_5/'
key=d5+'adult_scores/gate5_result_T3_repair.json';j=readj(key)
rows=[];nulls=[]
for vocab,v in j['by_vocabulary'].items():
    rows.append({'vocabulary':vocab,'score':'AUCell','n_patients':j['n_paired_patients'],'mean_change':v['mean_delta'],'median_change':v['median_delta'],'fraction_positive':v['frac_patients_positive'],'paired_Wilcoxon_p':v['p_wilcoxon']})
    for score,w in v['sensitivities'].items():
        rows.append({'vocabulary':vocab,'score':{'singscore':'singscore-based','meanz':'mean standardized expression'}[score],'n_patients':j['n_paired_patients'],'mean_change':w['mean_delta'],'median_change':None,'fraction_positive':w['frac_positive'],'paired_Wilcoxon_p':w['p_wilcoxon']})
    for k in ['nullA_label_permutation','nullB_random_gene_sets','nullC_random_regulon_modules']:
        w=v[k];nulls.append({'vocabulary':vocab,'null':k[:5],'n_draws':j['n_draws'],'mean':w['mean'],'SD':w['sd'],'empirical_p':w['p_empirical']})
table(5,'adult_module_scores',rows,[key],'Oriented, background-centered adult module changes (Figure 5B), including alternative scoring methods.','by_vocabulary: primary AUCell and sensitivities. Missing medians were not stored for alternative scores.',None,'Adult module scores')
table(5,'adult_nulls',nulls,[key],'Null distributions for the primary adult module score (Figure 5B).','by_vocabulary: null A, B and C.',None,'Adult nulls')
for f,name in [('per_patient_deltas_T3_repair.csv','adult_patient_changes'),('per_regulon_external_T3_repair_Diagnosis.csv','adult_regulons_diagnosis'),('per_regulon_external_T3_repair_Relapse.csv','adult_regulons_relapse'),('target_coverage_T3_repair.csv','adult_target_coverage'),('sample_table.csv','adult_samples')]:
    csvtable(5,name,d5+'adult_scores/'+f,'Adult expression-score changes, sample identifiers and target coverage (Figure 5A,B).')
key=d5+'module_definition/module_definition_T3_repair.json';j=readj(key)
members=[];targets=[]
for v in ['Diagnosis','Relapse']:
    for r in j['background_regulons']:
        members.append({'vocabulary':v,'regulon':r,'in_module_55':r in j['regulons'],'pediatric_orientation':j['orientation'][v][r],'n_targets':len(j['targets'][v][r])})
        for gene in j['targets'][v][r]: targets.append({'vocabulary':v,'regulon':r,'in_module_55':r in j['regulons'],'target_gene':gene})
table(5,'module_membership',members,[key],'Module membership and pediatric orientations within the 201-regulon background.','All background regulons under both vocabularies.',None,'Module membership')
table(5,'module_targets',targets,[key],'Target definitions for the module and its centering background.','All target genes in targets/<vocabulary>/<regulon>.')
for f in ['summary.csv','target_retention.csv','axis_genes.csv','per_patient_changes.csv','per_regulon_target_audit.csv','module_change_comparison.csv']:
    csvtable(5,'maturity_marker_exclusion/'+f[:-4],d5+'maturity_marker_exclusion/'+f,'Intact maturity axes compared with original and marker-excluded module scores (Figure 5D).',sheet='Maturity regressions' if f=='summary.csv' else None)
key=d5+'maturity_axes_legacy/axis_lambo_T3_repair.json';j=readj(key)
rows=[];earlier=[];controls=[]
for axis,v in j['axes'].items():
    if not axis.startswith(('A1_','A2_')):continue
    rows.append({'axis':axis,**v['delta_axis'],'positive_markers':v['pos_genes_used'],'negative_markers':v['neg_genes_used']})
    for vocab,w in v['by_vocabulary'].items():
        earlier.append({'axis':axis,'vocabulary':vocab,**{k:x for k,x in w['gene_disjoint'].items() if isinstance(x,(int,float)) and k!='survives_adjustment'}})
        for null in ['null_size_matched','null_detection_matched']:
            controls.append({'axis':axis,'vocabulary':vocab,'control':null,**{k:x for k,x in w[null].items() if isinstance(x,(int,float))}})
table(5,'maturity_axis_changes',rows,[key],'Diagnosis-to-relapse change in the intact axes. Larger scores indicate primitive marker expression.','axes/A1 and A2/delta_axis.',None,'Axis changes')
table(5,'maturity_axes_legacy/axis_target_exclusion',earlier,[key],'Earlier sensitivity excluding module targets from the axes, retained for Supplementary Methods S7 only.','axes/A1 and A2/by_vocabulary/*/gene_disjoint; not the source for Figure 5D.')
table(5,'maturity_axis_random_controls',controls,[key],'Size- and detection-matched random-axis controls (Supplementary Methods S7).','axes/A1 and A2/by_vocabulary/*/null_* fields.')
csvtable(5,'pooled_cell_type_activity',d5+'activity_and_function/f6_fig06_mixture.csv','Secondary pooled raw, uncentered and unoriented activity and six-cell-type shares (Figure 5C).')
csvtable(5,'signed_target_enrichment',d5+'activity_and_function/f6_fig06_arm_enrichment.csv','All tested terms for pooled signed target sets, before and after ribosomal depletion (Figure 5E).',rename={'arm':'TF_gene_correlation_sign'},values={'arm':{'activating (+)':'positive','repressing (-)':'negative'}})
csvtable(5,'tf_signed_target_counts',d5+'activity_and_function/f6_unit_correction_degree.csv','Target counts after separating each TF by its TF-gene correlation sign (Figure 5F).',where=lambda r:r['unit_kind']=='per_TF_s1',selection='unit_kind == per_TF_s1.')
csvtable(5,'tf_signed_enrichment',d5+'activity_and_function/f6_unit_correction_enrichment.csv','All tested terms for signed TF target sets (Figure 5F).',where=lambda r:r['unit_kind']=='per_TF_s1',selection='unit_kind == per_TF_s1; all terms, including nonsignificant results.')
csvtable(5,'ribosomal_genes',d5+'activity_and_function/f6_ribosome_gene_list.csv','Ribosomal protein genes excluded in the reported depletion sensitivity.',['symbol','class'],where=lambda r:r['in_DROP_A']=='True',selection='in_DROP_A == True.')
for f in ['enrichment_counts','gene_membership','go_enrichment_all_terms','go_exclusive_target_contrast','per_regulon_target_overlap','target_coverage']:
    csvtable(5,'stagewise_function/'+f,d5+'stagewise_function/'+f+'.csv','Stage-specific functional enrichment and direct stage-exclusive target comparisons (Supplementary Methods S8).')

d6='Supplementary_Data_6/'
key=d6+'decomposition/f6_decomposition.json';j=readj(key);rows=[]
for scale in ['oriented_offsetcentred','offsetcentred','raw']:
    for vocab,v in j['module_decomposition'][scale].items():
        rows.append({'vocabulary':vocab,'scale':scale,'analysis_role':'primary' if scale=='oriented_offsetcentred' else 'secondary','n_patients':v['n_patients'],'total_change':v['naive_mean'],'within_cell_type_change':v['rewiring_mean'],'composition_component':v['composition_mean'],'within_cell_type_p':v['p_rewiring'],'composition_p':v['p_composition'],'composition_share':v['composition_share']})
table(6,'module_decomposition',rows,[key],'Pediatric module decomposition; oriented centered scores are primary (Figure 5C).','module_decomposition: all three stored scales.',{'within_cell_type_change':'rewiring_mean; activity change, not regulatory edge changes.','total_change':'naive_mean','composition_share':'composition_mean / naive_mean; unstable when total change is close to zero.'},'Module decomposition')
csvtable(6,'patient_decomposition',d6+'decomposition/f6_decomposition_module_by_patient.csv','Patient-level total, within-cell-type and composition components.',rename={'naive':'total_change','rewiring':'within_cell_type_change','composition':'composition_component','abs_composition_gt_abs_rewiring':'abs_composition_gt_abs_within_cell_type','naive_rewiring_same_sign':'total_and_within_cell_type_same_sign'})
csvtable(6,'regulon_decomposition',d6+'decomposition/f6_decomposition_set55_per_regulon.csv','Per-regulon pediatric decomposition and statistical comparisons.',rename={'rewiring_mean':'within_cell_type_change','naive_mean':'total_change','p_rewiring':'within_cell_type_p'})
numeric(6,'clp_sensitivity_summary',d6+'CLP_sensitivity/sens_clp7_summary.json',[[k] for k in ['paired_patients_per_celltype','funnels','TF_level_vs_primary','decomposition']],'Seven-cell-type sensitivity including CLP.')
for f,name in [('sens_clp7_funnel_by_regulon.csv','clp_candidate_criteria'),('sens_clp7_per_regulon_decomposition.csv','clp_regulon_decomposition'),('per_patient_deltas_SENS_CLP7_5of7.csv','clp_adult_patient_changes')]:
    csvtable(6,name,d6+'CLP_sensitivity/'+f,'CLP sensitivity results.',rename={'rewiring_mean':'within_cell_type_change','naive_mean':'total_change','p_rewiring':'within_cell_type_p'})
for f,name in [('gate5_result_SENS_CLP7_5of7.json','clp_adult_scores'),('gate5_result_CTRL55_sens_clp7.json','clp_fixed_module_adult_scores')]:
    numeric(6,name,d6+'CLP_sensitivity/'+f,[['by_vocabulary']],'Adult score and null summaries for the CLP sensitivity.')
for f in ['all_candidate_fold_decisions','baseline_rebuilt_nomination','external_sample_pairs','fixed_module_patient_deltas','fixed_module_patient_influence','nomination_by_omitted_patient','regulon_selection_frequency','reselected_module_scores']:
    csvtable(6,'patient_omission/'+f,d6+'patient_omission/'+f+'.csv','Patient-omission results (Supplementary Figures S1,S2).',sheet='Reselected sets' if f=='nomination_by_omitted_patient' else None)

for n,title in titles.items():
    entries=[r for r in catalog if r['table'].startswith(f'Supplementary_Data_{n}/')]
    text=[f'Supplementary Data {n}. {title}','', 'CSV tables contain the detailed results. The Excel workbook presents selected tables for reading.',
          'The table index identifies source files and selection rules. Source files are preserved under source data/.',
          'source_row refers to the original CSV record, including its header. Blank values are missing or not computed, not zero.',
          'No network inference or statistical tests were rerun during curation. All selections retain nonsignificant results.', '', 'Tables:']
    text.extend(f"{Path(r['table']).relative_to(f'Supplementary_Data_{n}')} - {r['description']}" for r in entries)
    (output/f'Supplementary_Data_{n}'/'README.txt').write_text('\n'.join(text)+'\n',encoding='utf-8')
column_notes = {'matched_reuse_summary.csv': {'n_assignments': 'arms/<arm>/n_assignments',
                               'maximum_reuse': 'arms/<arm>/reuse_max',
                               'native_cells_recovered': 'arms/<arm>/self_hits',
                               'exact_cell_recovery': 'arms/<arm>/exact_cell_accuracy',
                               'partner_selection': 'Arm label: link_* is inferred; rand_* is '
                                                    'random.'},
 'paml_reuse.csv': {'n_assignments': 'stages/<stage>/n_in for published link; n_out for corrected.',
                    'distinct_RNA_profiles': 'stages/<stage>/old_distinct_donors for published '
                                             'link; distinct_donors for corrected. These are RNA '
                                             'profiles.',
                    'reuse_ratio': 'stages/<stage>/old_reuse for published link; achieved_reuse '
                                   'for corrected.',
                    'corrected_scarcity_floor': 'stages/<stage>/scarcity_floor; not assigned to '
                                                'published link.'},
 'adult_module_scores.csv': {'mean_change': 'by_vocabulary/<vocabulary>/mean_delta, or '
                                            'sensitivities/<score>/mean_delta.',
                             'median_change': 'by_vocabulary/<vocabulary>/median_delta for AUCell; '
                                              'not stored for alternative scores.',
                             'paired_Wilcoxon_p': 'by_vocabulary/<vocabulary>/p_wilcoxon, or '
                                                  'sensitivities/<score>/p_wilcoxon.',
                             'fraction_positive': 'frac_patients_positive for AUCell; '
                                                  'sensitivities/<score>/frac_positive for '
                                                  'alternative scores.'},
 'adult_nulls.csv': {'null': 'A: patient sign flips; B: random targets; C: random regulon modules.',
                     'mean': 'by_vocabulary/<vocabulary>/null*/mean',
                     'SD': 'by_vocabulary/<vocabulary>/null*/sd',
                     'empirical_p': 'by_vocabulary/<vocabulary>/null*/p_empirical'},
 'module_membership.csv': {'in_module_55': 'Membership in the regulons list.',
                           'pediatric_orientation': 'orientation/<vocabulary>/<regulon>',
                           'n_targets': 'Length of targets/<vocabulary>/<regulon>.'},
 'module_targets.csv': {'target_gene': 'Each member of targets/<vocabulary>/<regulon>.',
                        'in_module_55': 'Membership in the regulons list.'},
 'module_decomposition.csv': {'scale': 'Key under module_decomposition.',
                              'analysis_role': 'oriented_offsetcentred is primary; other scales '
                                               'are secondary.',
                              'within_cell_type_p': 'module_decomposition/<scale>/<vocabulary>/p_rewiring',
                              'composition_p': 'module_decomposition/<scale>/<vocabulary>/p_composition',
                              'composition_component': 'module_decomposition/<scale>/<vocabulary>/composition_mean'}}
for record in columns:
    definition=column_notes.get(Path(record["table"]).name,{}).get(record["column"])
    if definition:
        record["source_field_or_definition"]=definition

for name,rows in [('table_index',catalog),('column_sources',columns)]:
    with (output/(name+'.csv')).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(output/'curation_manifest.json').write_text(json.dumps(catalog,indent=2)+'\n')
if args.workbook_spec:
    (output/'workbook_tables.json').write_text(json.dumps(books,indent=2,allow_nan=False)+'\n')
print(json.dumps({'tables':len(catalog),'rows':sum(r['rows'] for r in catalog),'workbooks':len(books)}))
