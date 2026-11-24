#!/usr/bin/env python
"""Compare functional annotations of stage-specific regulon target sets.

Stagewise enrichment and direct exclusive-target contrasts answer different
questions. The analysis is exploratory and does not test patient-level activity."""
import argparse
import hashlib
import json
import platform
import re
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, hypergeom

STAGES = ('Diagnosis', 'Relapse')
GROUPS = ('all55', 'TF2G_positive', 'TF2G_negative')
GMT_SHA = '0fc02458765cfce6c9348dce9f7a9397c6caf7c09ba318f468135d7ac60342ee'

def bh(p):
    p = np.asarray(p, dtype=float)
    if not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError('BH requires finite p-values in [0,1]')
    if not len(p):
        return p
    order = np.argsort(p, kind='stable')
    q = np.empty(len(p))
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def symbols(items):
    return {str(x).strip() for x in items if pd.notna(x) and str(x).strip()}

def read_gmt(path):
    terms = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) >= 3:
            if fields[0] in terms:
                raise ValueError('Duplicate GO term: ' + fields[0])
            terms[fields[0]] = symbols(fields[2:])
    return terms

def union_targets(targets, regs):
    return set().union(*(targets[r] for r in regs))

def ora(fg, bg, terms, minimum=10, maximum=500):
    if not fg <= bg:
        raise ValueError('Foreground outside background')
    if not fg or not bg or fg == bg:
        raise ValueError('Empty or saturated foreground: enrichment is uninformative')
    rows = []
    for term, genes in sorted(terms.items()):
        gs = genes & bg
        if not minimum <= len(gs) <= maximum:
            continue
        hits = gs & fg
        a, b, c = len(hits), len(fg-hits), len(gs-hits)
        d = len(bg)-a-b-c
        if min(a,b,c,d) < 0:
            raise ValueError('Invalid contingency table')
        # Cross product OR: structural zero denominator is kept as infinity, not capped.
        odds = (a*d)/(b*c) if b*c else (np.inf if a*d else np.nan)
        rows.append(dict(term=term, n_overlap=a, n_term_in_bg=len(gs),
                         n_foreground=len(fg), n_background=len(bg),
                         odds_ratio=odds, fold_enrichment=(a/len(fg))/(len(gs)/len(bg)),
                         gene_ratio=a/len(fg), background_ratio=len(gs)/len(bg),
                         overlap_genes=';'.join(sorted(hits)), a=a,b=b,c=c,d=d))
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError('No eligible GO terms')
    frame['p'] = hypergeom.sf(frame.n_overlap-1, frame.n_background,
                              frame.n_term_in_bg, frame.n_foreground)
    frame['q_stage'] = bh(frame.p)
    return frame

def exclusive_contrast(dx, rel, bg, terms, minimum=10, maximum=500):
    """Conditional gene-set annotation contrast; genes, not patients, are units."""
    if not dx <= bg or not rel <= bg:
        raise ValueError('Contrast foreground outside background')
    d_only, r_only = dx-rel, rel-dx
    valid = bool(d_only and r_only)
    rows = []
    for term, genes in sorted(terms.items()):
        gs = genes & bg
        if not minimum <= len(gs) <= maximum:
            continue
        # Row 1 Rel-only, row 2 Dx-only, so positive log2 OR points to Rel-only.
        a,b,c,d = len(r_only & gs),len(r_only-gs),len(d_only & gs),len(d_only-gs)
        if valid:
            odds, p = fisher_exact([[a,b],[c,d]], alternative='two-sided')
            effect = np.log2(((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5)))
        else:
            odds,p,effect = np.nan,np.nan,np.nan
        rows.append(dict(term=term, relapse_only_hits=a, relapse_only_total=a+b,
                         diagnosis_only_hits=c, diagnosis_only_total=c+d,
                         shared_targets=len(dx & rel), odds_ratio_rel_vs_dx=odds,
                         log2_odds_ratio_rel_vs_dx_Haldane=effect,p=p,
                         status='descriptive_gene_set_test' if valid else 'undefined_empty_exclusive_set',
                         relapse_only_hit_genes=';'.join(sorted(r_only & gs)),
                         diagnosis_only_hit_genes=';'.join(sorted(d_only & gs))))
    result = pd.DataFrame(rows)
    result['q'] = bh(result.p) if valid and len(result) else np.nan
    return result

def verify_export(path, expected):
    frame = pd.read_csv(path, sep='\t', usecols=['eRegulon_name','Gene'])
    observed = frame.groupby('eRegulon_name').Gene.agg(symbols).to_dict()
    missing = set(expected)-set(observed)
    mismatch = [r for r in expected if observed.get(r,set()) != expected[r]]
    if missing or mismatch:
        raise ValueError(f'{path}: target mismatch; missing={sorted(missing)}, different={mismatch}')
    return dict(checked_regulons=len(expected), triplets=len(frame), mismatches=0)

def parity_check(computed, previous):
    old = pd.read_csv(previous)
    old = old[old.library == 'GOBP'].copy()
    old['group'] = old.arm.map({'activating (+)':'TF2G_positive','repressing (-)':'TF2G_negative'})
    new = computed[(computed.background == 'stage_pool') & (computed.group != 'all55')].copy()
    keys = ['vocabulary','gene_set','group','term']
    merged = old.merge(new,on=keys,how='outer',suffixes=('_old','_new'),indicator=True,validate='one_to_one')
    if not merged['_merge'].eq('both').all():
        raise ValueError('Existing GOBP term families were not reproduced exactly')
    for col in ['n_overlap','n_term_in_bg','n_foreground','n_background']:
        if not (merged[col+'_old'] == merged[col+'_new']).all():
            raise ValueError('Existing enrichment count mismatch: ' + col)
    if not np.allclose(merged.p_old,merged.p_new,rtol=1e-9,atol=0):
        raise ValueError('Hypergeometric p-values differ from existing Fisher results')
    if not np.allclose(merged.q,merged.q_stage,rtol=1e-9,atol=0):
        raise ValueError('Per-stage BH values differ from existing enrichment')
    return dict(status='PASS', rows=len(merged),
                max_abs_p_difference=float(np.max(np.abs(merged.p_old-merged.p_new))),
                max_abs_q_difference=float(np.max(np.abs(merged.q-merged.q_stage))))

def display_terms(frame, terms, bg, maximum):
    ranked = []
    for stage in STAGES:
        d = frame[(frame.vocabulary == stage) & (frame.q_both_stages < .05)]
        ranked.append(d.sort_values(['q_both_stages','p','term']).term.tolist())
    ordered = []
    for i in range(max(map(len,ranked),default=0)):
        for candidates in ranked:
            if i < len(candidates) and candidates[i] not in ordered:
                ordered.append(candidates[i])
    chosen = []
    for term in ordered:
        gs = terms[term] & bg
        if any(len(gs & (terms[t] & bg))/len(gs | (terms[t] & bg)) > .70 for t in chosen):
            continue
        chosen.append(term)
        if len(chosen) >= maximum:
            break
    return chosen

def draw(frame, selections, figdir, gene_set, all55=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    groups = ['all55'] if all55 else ['TF2G_positive','TF2G_negative']
    # One panel per sign; common q colour scale and overlap-area scale.
    height = max(4.5, max([len(selections[g]) for g in groups]+[1])*.42+1.65)
    width = 7 if all55 else 13
    fig, axes = plt.subplots(1,len(groups),figsize=(width,height),squeeze=False)
    fig.subplots_adjust(left=.32 if all55 else .25,right=.87 if all55 else .89,
                        bottom=.19,top=.89,wspace=2.1)
    if gene_set != 'full':
        fig.suptitle('Ribosomal-protein genes removed from targets and background',
                     fontsize=10, y=.99)
    norm = Normalize(0,20)
    titles = {'all55':'All 55 regulons (target union)',
              'TF2G_positive':'Positive TF-to-gene correlation',
              'TF2G_negative':'Negative TF-to-gene correlation'}
    sizes = [10,50,150]
    for ax,g in zip(axes[0],groups):
        chosen = selections[g]
        d = frame[frame.group == g]
        ax.set_title(titles[g],fontsize=10,pad=17)
        if not chosen:
            ax.text(.5,.5,'No terms pass q < 0.05',ha='center',transform=ax.transAxes)
            ax.set_yticks([])
        else:
            for j,stage in enumerate(STAGES):
                sub = d[d.vocabulary == stage].set_index('term').loc[chosen]
                for i,(_,row) in enumerate(sub.iterrows()):
                    q = float(row.q_both_stages)
                    fill = plt.cm.Blues(norm(min(20,-np.log10(max(q,1e-300)))))
                    if q >= .05:
                        fill = 'white'
                    ax.scatter(j,i,s=18+1.1*row.n_overlap,facecolors=fill,
                               edgecolors='#303030' if q < .05 else '#aaaaaa',linewidths=.7,zorder=3)
            labels = ['\n'.join(textwrap.wrap(re.sub(r' \(GO:\d+\)$','',t),width=35)) for t in chosen]
            ax.set_yticks(np.arange(len(chosen)),labels,fontsize=8)
            ax.set_ylim(len(chosen)-.45,-.6)
            ax.grid(axis='y',color='#e8e8e8',linewidth=.6,zorder=0)
        ax.set_xlim(-.45,1.45)
        labs=[]
        for stage in STAGES:
            n = int(d[d.vocabulary == stage].n_foreground.iloc[0])
            labs.append(f'{stage}\nnetwork\n{n:,} targets')
        ax.set_xticks([0,1],labs,fontsize=9)
        ax.tick_params(axis='y',length=0,pad=5)
        for s in ['top','right','left']:
            ax.spines[s].set_visible(False)
    cax = fig.add_axes([.92 if not all55 else .91,.38,.012,.30])
    sm = plt.cm.ScalarMappable(cmap='Blues',norm=norm)
    cb=fig.colorbar(sm,cax=cax,extend='max')
    cb.set_label('-log10 q (both stages)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    handles=[Line2D([0],[0],marker='o',linestyle='',markerfacecolor='#4292C6',
                    markeredgecolor='#303030',markersize=np.sqrt(18+1.1*n),label=str(n)) for n in sizes]
    handles.append(Line2D([0],[0],marker='o',linestyle='',markerfacecolor='white',
                          markeredgecolor='#aaaaaa',label='q >= 0.05'))
    fig.legend(handles=handles,loc='lower center',ncol=4,frameon=False,
               title='GO overlap (genes); open symbols do not pass the threshold',fontsize=8,title_fontsize=8)
    stem = 'go_stage_all55' if all55 else 'go_stage_sign_split'
    if gene_set != 'full':
        stem += '_ribosome_depleted'
    fig.savefig(figdir/(stem+'.svg'),bbox_inches='tight')
    fig.savefig(figdir/(stem+'.png'),dpi=240,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--out',type=Path)
    ap.add_argument('--figures',type=Path)
    ap.add_argument('--no-plots',action='store_true')
    args=ap.parse_args()
    root=args.root.resolve()
    out=args.out or root/'analysis/out/go_stage_comparison'
    figdir=args.figures or root/'figures/supplementary/functional_enrichment'
    inputs={'module':root/'modules/out/corrected_module.json',
            'GO_BP_2025':root/'analysis/ref/GO_Biological_Process_2025.gmt',
            'ribosomal_genes':root/'analysis/out/f6_ribosome_gene_list.csv',
            'existing_enrichment':root/'analysis/out/f6_fig06_arm_enrichment.csv'}
    for v in STAGES:
        inputs['universe_'+v]=root/f'analysis/out/f1_universe_{v}.csv'
        inputs['network_'+v]=root/f'paml/malignant/out/eRegulons_{v}_repair.tsv'
    absent=[str(p) for p in inputs.values() if not p.is_file()]
    if absent:
        raise FileNotFoundError('Missing required inputs: '+', '.join(absent))
    module=json.loads(inputs['module'].read_text(encoding='utf-8'))
    regs,pool=module['regulons'],module['background_regulons']
    if len(regs)!=55 or len(set(regs))!=55 or len(pool)!=201 or len(set(pool))!=201 or not set(regs)<=set(pool):
        raise ValueError('The module must be the frozen 55 of 201')
    if any(re.fullmatch(r'.+_direct_[+-]/[+-]',r) is None for r in pool):
        raise ValueError('Unexpected regulon sign/name format')
    if len({r.split('_direct_')[0] for r in regs})!=31:
        raise ValueError('Expected 31 TFs')
    groups={'all55':regs,
            'TF2G_positive':[r for r in regs if r.split('_direct_')[1][0]=='+'],
            'TF2G_negative':[r for r in regs if r.split('_direct_')[1][0]=='-']}
    target={v:{r:symbols(module['targets'][v][r]) for r in pool} for v in STAGES}
    checks={v:verify_export(inputs['network_'+v],target[v]) for v in STAGES}
    universe={v:symbols(pd.read_csv(inputs['universe_'+v]).symbol) for v in STAGES}
    poolgenes={v:union_targets(target[v],pool) for v in STAGES}
    common=(poolgenes['Diagnosis']|poolgenes['Relapse']) & (universe['Diagnosis']&universe['Relapse'])
    if not common:
        raise ValueError('No common eligible gene background')
    terms=read_gmt(inputs['GO_BP_2025'])
    if digest(inputs['GO_BP_2025'])!=GMT_SHA:
        raise ValueError('Frozen GO library SHA-256 changed; review before rerunning')
    ribo=pd.read_csv(inputs['ribosomal_genes'],dtype=str)
    drop=symbols(ribo.loc[ribo.in_DROP_A.str.lower().eq('true'),'symbol'])
    if len(drop)!=167:
        raise ValueError('Expected the existing 167-gene DROP_A definition')
    rows=[]
    contrast=[]
    overlaps=[]
    memberships=[]
    regcounts=[]
    for r in sorted(regs):
        a,b=target['Diagnosis'][r],target['Relapse'][r]
        regcounts.append(dict(regulon=r,diagnosis_targets=len(a),relapse_targets=len(b),
                              shared_targets=len(a&b),union_targets=len(a|b),
                              target_jaccard=len(a&b)/len(a|b)))
    for background in ['common_pool','stage_pool']:
        for gsname,deplete in [('full',set()),('ribosome_depleted',drop)]:
            sets={}
            for v in STAGES:
                bg=(common if background=='common_pool' else poolgenes[v]&universe[v])-deplete
                memberships.extend(dict(background=background,gene_set=gsname,vocabulary=v,
                                        group='background',symbol=g) for g in sorted(bg))
                for group,selected in groups.items():
                    raw=union_targets(target[v],selected)
                    fg=raw & bg
                    sets[(v,group)]=fg
                    d=ora(fg,bg,terms)
                    d['background'],d['gene_set'],d['vocabulary'],d['group']=background,gsname,v,group
                    rows.append(d)
                    memberships.extend(dict(background=background,gene_set=gsname,vocabulary=v,
                                            group=group,symbol=g) for g in sorted(fg))
                    overlaps.append(dict(background=background,gene_set=gsname,vocabulary=v,group=group,
                                         n_regulons=len(selected),raw_targets=len(raw),eligible_targets=len(fg),
                                         excluded_targets=len(raw-fg),excluded_genes=';'.join(sorted(raw-fg)),
                                         n_background=len(bg),GO_annotated_targets=len(fg&set().union(*terms.values()))))
            if background=='common_pool':
                for group in groups:
                    c=exclusive_contrast(sets[('Diagnosis',group)],sets[('Relapse',group)],common-deplete,terms)
                    c['group'],c['gene_set']=group,gsname
                    contrast.append(c)
    enrichment=pd.concat(rows,ignore_index=True)
    families=['background','gene_set','group']
    enrichment['q_both_stages']=enrichment.groupby(families).p.transform(lambda s:bh(s.to_numpy()))
    enrichment['q_all_groups_both_stages']=enrichment.groupby(['background','gene_set']).p.transform(lambda s:bh(s.to_numpy()))
    parity=parity_check(enrichment,inputs['existing_enrichment'])
    print('Existing enrichment parity:',parity,flush=True)
    contrast=pd.concat(contrast,ignore_index=True)
    contrast['q_all_groups']=np.nan
    for _,ix in contrast[contrast.p.notna()].groupby('gene_set').groups.items():
        contrast.loc[ix,'q_all_groups']=bh(contrast.loc[ix,'p'])
    out.mkdir(parents=True,exist_ok=True)
    figdir.mkdir(parents=True,exist_ok=True)
    enrichment.to_csv(out/'go_enrichment_all_terms.csv',index=False)
    contrast.to_csv(out/'go_exclusive_target_contrast.csv',index=False)
    pd.DataFrame(overlaps).to_csv(out/'target_coverage.csv',index=False)
    pd.DataFrame(memberships).to_csv(out/'gene_membership.csv',index=False)
    pd.DataFrame(regcounts).to_csv(out/'per_regulon_target_overlap.csv',index=False)
    displayed=[]
    for gsname,deplete in [('full',set()),('ribosome_depleted',drop)]:
        d=enrichment[(enrichment.background=='common_pool')&(enrichment.gene_set==gsname)]
        sel={g:display_terms(d[d.group==g],terms,common-deplete,16 if g=='all55' else 12) for g in groups}
        for g,tt in sel.items():
            for rank,term in enumerate(tt,1):
                displayed.append(dict(gene_set=gsname,group=g,display_order=rank,term=term))
        if not args.no_plots:
            draw(d,sel,figdir,gsname,all55=False)
            draw(d,sel,figdir,gsname,all55=True)
    pd.DataFrame(displayed).to_csv(out/'displayed_terms.csv',index=False)
    counts=[]
    for keys,d in enrichment.groupby(['background','gene_set','group'],sort=False):
        sig={v:set(d[(d.vocabulary==v)&(d.q_both_stages<.05)].term) for v in STAGES}
        counts.append(dict(background=keys[0],gene_set=keys[1],group=keys[2],
                           diagnosis_sig=len(sig['Diagnosis']),relapse_sig=len(sig['Relapse']),
                           both_sig=len(sig['Diagnosis']&sig['Relapse']),
                           diagnosis_only_sig=len(sig['Diagnosis']-sig['Relapse']),
                           relapse_only_sig=len(sig['Relapse']-sig['Diagnosis'])))
    pd.DataFrame(counts).to_csv(out/'enrichment_counts.csv',index=False)
    import matplotlib
    summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),
                 status='POST-HOC DESCRIPTIVE; stage-derived targets, not patient-level pathway change',
                 root=str(root),command=sys.argv,python=sys.version,
                 packages=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
                 platform=platform.platform(),n_regulons=55,n_TFs=31,n_background_regulons=201,
                 regulons_by_group={g:len(r) for g,r in groups.items()},
                 common_background_genes=len(common),ribosomal_definition_genes=len(drop),
                 ribosomal_genes_in_common_background=len(drop&common),
                 target_export_checks=checks,existing_enrichment_parity=parity,
                 enrichment_rows=len(enrichment),exclusive_contrast_rows=len(contrast),enrichment_counts=counts,
                 script_sha256=digest(Path(__file__)),
                 inputs={key:dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=digest(p)) for key,p in inputs.items()},
                 q_both_stages_family='GO terms x two stages within background x depletion x target group',
                 q_all_groups_family='GO terms x two stages x three groups within background x depletion',
                 direct_contrast_family='GO terms per depletion x target group; q_all_groups also supplied',
                 interpretation=['Different significant lists do not establish differential enrichment.',
                                 'Exclusive-target Fisher tests treat genes as exchangeable and do not measure patient/network uncertainty.',
                                 'Positive and negative TF2G correlations do not prove causal activation/repression.',
                                 'The 55 were adopted post hoc; no new confirmatory claims are licensed.'])
    (out/'run_manifest.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    primary=[c for c in counts if c['background']=='common_pool' and c['gene_set']=='full']
    print(json.dumps(dict(common_background_genes=len(common),counts=primary),indent=2))
    caption=('GO Biological Process enrichment of inferred targets of the frozen 55-regulon module, '
             'split by the sign of the TF-to-gene correlation and evaluated separately in networks inferred '
             'at diagnosis and relapse. Both columns use the union of the 201-regulon pool target genes '
             'across networks, restricted to genes in both expressed universes '
             f'(n = {len(common):,}). Dot area represents gene overlap; colour represents -log10 BH q, '
             'adjusted across all eligible GO terms and both network definitions within each sign group. '
             'Open dots do not pass q < 0.05. Terms were selected from either column by a deterministic '
             'rank-based rule with gene-set redundancy filtering; all tests and the display rule are supplied. '
             'The sensitivity figure removes the existing 167-gene ribosomal-protein definition from '
             'foregrounds and backgrounds. Different enrichment patterns describe inferred target annotations; '
             'they do not establish pathway activation, causal regulation or patient-level rewiring. '
             'This analysis is post hoc and descriptive.\n')
    (out/'proposed_caption.txt').write_text(caption,encoding='utf-8')
    print('Wrote',out,'and',figdir,flush=True)

if __name__=='__main__':
    main()
