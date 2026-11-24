"""Readable computed findings for the submission planning memo."""
from pathlib import Path
import json
import pandas as pd

root=Path(__file__).resolve().parents[2]
folder=root/'analysis/out/go_stage_comparison'
enr=pd.read_csv(folder/'go_enrichment_all_terms.csv')
con=pd.read_csv(folder/'go_exclusive_target_contrast.csv')
mem=pd.read_csv(folder/'gene_membership.csv')
results={}
for group in ['all55','TF2G_positive','TF2G_negative']:
    d=mem[(mem.background=='common_pool')&(mem.gene_set=='full')&(mem.group==group)]
    a=set(d[d.vocabulary=='Diagnosis'].symbol)
    b=set(d[d.vocabulary=='Relapse'].symbol)
    c=con[(con.group==group)&(con.gene_set=='full')]
    dep=con[(con.group==group)&(con.gene_set=='ribosome_depleted')]
    top=c.sort_values(['q','p','term']).head(8)[['term','relapse_only_hits','relapse_only_total','diagnosis_only_hits','diagnosis_only_total','log2_odds_ratio_rel_vs_dx_Haldane','p','q','q_all_groups']]
    results[group]=dict(diagnosis_targets=len(a),relapse_targets=len(b),shared=len(a&b),
                        diagnosis_only=len(a-b),relapse_only=len(b-a),union_jaccard=len(a&b)/len(a|b),
                        exclusive_annotation_contrast_q_lt_005=int((c.q<.05).sum()),
                        exclusive_annotation_contrast_q_all_groups_lt_005=int((c.q_all_groups<.05).sum()),
                        exclusive_annotation_contrast_depleted_q_lt_005=int((dep.q<.05).sum()),
                        leading_exclusive_annotation_contrasts=top.to_dict(orient='records'))
print(json.dumps(results,indent=2))
(folder/'findings_summary.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
