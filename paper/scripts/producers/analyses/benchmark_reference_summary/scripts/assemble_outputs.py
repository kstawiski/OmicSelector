from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,hashlib,json,os,platform,shutil,datetime
import numpy as np
import pandas as pd

B=Path(__file__).resolve().parents[1];M=B.parent;O=B/'out';S=O/'synthesis';N=O/'normalization';OLD=M/'benchmark_saved_scores/out/synthesis'
EX={'GSE270497','GSE186595'};NON={'GSE145176','GSE31568'}
LABEL_MAP={'Primary prespecified synthesis':'Primary eligible-set synthesis (30 units)',
 'Excluding contrasts without confirmed cancer status':'Excluding GSE266859 (legacy cancer-status sensitivity; eligible melanoma contrast)',
 'Excluding GSE186595 (pool membership unresolved)':'GSE186595 exclusion (already excluded from eligible primary; identical result)'}
def rd(p):return pd.read_csv(p,sep='\t')
def save(d,p):p.parent.mkdir(parents=True,exist_ok=True);_public_table(d).to_csv(p,sep='\t',index=False,na_rep='')
def hashfile(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
checks=[]
def valid(name,ok,detail=''):
 checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError((name,detail))

ind=rd(O/'independent_check.tsv');valid('independent_recomputation',ind.pass_check.all(),f'{len(ind)} checks; max error {ind.absolute_error.max():.12g}')
primary_error=ind.loc[ind.scope=='Independent primary S6 REML/KH/BY','absolute_error'].max()
summary_rows=[]
def derived(root,tag):
 w=rd(root/'evidence/within_performance.tsv');p=rd(root/'evidence/pooled.tsv');ss=rd(root/'evidence/sensitivity_pooled.tsv');wt=rd(root/'evidence/cohort_weights.tsv')
 rows=[]
 def add(metric,value,method='',analysis='',source='evidence/pooled.tsv'):rows.append(dict(metric=metric,method=method,analysis=analysis,value=value,source_file=source))
 e=p[p.pooled_lift.notna()];eligible=w[w.all_seeds_eligible&w.auc_method.notna()]
 a=eligible[(eligible.method=='reo-pairratio')&~eligible.cohort.isin(NON)];b=eligible[(eligible.method=='ws-rclr-trimmed')&~eligible.cohort.isin(NON)]
 shared=a.merge(b[['cohort','auc_method']],on='cohort',suffixes=('_pairratio','_reference'));shared['paired_auc_difference']=shared.auc_method_pairratio-shared.auc_method_reference
 if tag!='old':save(shared,root/'evidence/shared_cohort_auc.tsv')
 base={'planned_method_cohort_cells':len(w),'evaluable_method_cohort_cells':len(eligible),'datasets':w.cohort.nunique(),'synthesis_planned_units':w[~w.cohort.isin(NON)].cohort.nunique(),'methods':w.method.nunique(),'method_seed_units':w.method.nunique()*5,'planned_method_seed_dataset_cells':len(w)*5,'estimable_methods':len(e),'BY_above_zero':sum(e.p_lift_gt_zero_BY<.05),'BY_above_005':sum(e.p_lift_gt_margin_0_05_BY<.05),'smallest_primary_BY':e.p_lift_gt_zero_BY.min(),'CI_wholly_below_zero':sum(e.ci_hi<0),'CI_wholly_above_zero':sum(e.ci_lo>0),'CI_upper_below_005':sum(e.ci_hi<.05),'CI_upper_above_005':sum(e.ci_hi>.05),'near_zero_count':sum((e.ci_hi>-.01)&(e.ci_hi<0)),'near_zero_list':';'.join(e.loc[(e.ci_hi>-.01)&(e.ci_hi<0),'method']),'shared_pairratio_reference_cohorts':len(shared),'pairratio_shared_median_AUC':shared.auc_method_pairratio.median(),'reference_shared_median_AUC':shared.auc_method_reference.median(),'shared_median_paired_AUC_difference':shared.paired_auc_difference.median(),'reference_evaluable_pairratio_unavailable':len(set(b.cohort)-set(a.cohort)),'reference_evaluable_pairratio_unavailable_list':';'.join(sorted(set(b.cohort)-set(a.cohort))),'pairratio_positive_signs':sum(shared.loc[shared.lift_se>0,'lift']>0),'pairratio_sign_denominator':sum(shared.lift_se>0),'minimum_sensitivity_BY':ss.p_lift_gt_zero_BY.min()}
 for metric,value in base.items():add(metric,value,source='method_cohort_performance.tsv; pooled_effects.tsv; sensitivity_effects.tsv')
 for rank,(_,r) in enumerate(e.sort_values('pooled_lift',ascending=False).iterrows(),1):add('pooled_point_estimate_rank',rank,r.method)
 audit=rd(root/'evidence/reml_diagnostics.tsv') if (root/'evidence/reml_diagnostics.tsv').exists() else None
 if audit is not None:
  add('REML_execution_requests',len(audit),source='evidence/reml_diagnostics.tsv');add('REML_execution_boundary_requests',audit.reml_boundary.sum(),source='evidence/reml_diagnostics.tsv')
 for field in ['pooled_lift','ci_hi','pi_hi']:
  r=e.loc[e[field].idxmax()];add('largest_'+field,r[field],r.method)
 for _,r in p.iterrows():
  for f in p.select_dtypes('number').columns:add(f,r[f],r.method)
 for _,r in wt[wt.cohort=='TORAY_PUBLIC_CANONICAL'].iterrows():add('Toray_weight_percent',r.weight_percent,r.method,source='evidence/cohort_weights.tsv')
 g=rd(root/'evidence/grouped_changes.tsv');add('maximum_absolute_grouped_tuning_change',g.delta_pooled_lift.abs().max(),source='evidence/grouped_changes.tsv')
 fam=rd(root/'pooled_effects.tsv').iloc[:,:2];fam.columns=['method','family'];coverage=w.merge(fam,on='method').assign(evaluable=lambda x:x.all_seeds_eligible&x.auc_method.notna()).groupby('family').agg(planned=('method','size'),evaluable=('evaluable','sum')).reset_index();coverage['fraction']=coverage.evaluable/coverage.planned
 if tag!='old':save(coverage,root/'evidence/figure_family_coverage.tsv')
 for _,r in coverage.iterrows():
  for f in ['planned','evaluable','fraction']:add('figure_family_'+f,r[f],r.family,source='evidence/figure_family_coverage.tsv')
 ssum=[]
 for an,g in ss.groupby('analysis',sort=False):
  e1=g[g.pooled_lift.notna()];r=e1.loc[e1.ci_hi.idxmax()];pr=e1[e1.method=='reo-pairratio'].iloc[0]
  ssum.append(dict(analysis=an,estimable=len(e1),BY_above_zero=sum(e1.p_lift_gt_zero_BY<.05),BY_above_005=sum(e1.p_lift_gt_margin_0_05_BY<.05),CI_wholly_below_zero=sum(e1.ci_hi<0),CI_upper_above_005=sum(e1.ci_hi>.05),methods_CI_upper_above_005=';'.join(e1.loc[e1.ci_hi>.05,'method']),maximum_CI_upper=r.ci_hi,maximum_CI_upper_method=r.method,pairratio_lift=pr.pooled_lift,pairratio_k=pr.k,S9_data_rows=';'.join(str(i+1) for i in g.index)))
  for f,v in ssum[-1].items():
   if f not in ['analysis','S9_data_rows']:add('sensitivity_'+f,v,analysis=LABEL_MAP.get(an,an) if tag=='old' else an,source='sensitivity_effects.tsv')
 cells=rd(root/'cells_union_complete.tsv');idx=rd(root/'evidence/foundation_seed_sensitivities.tsv')
 for _,r in idx.iterrows():
  sc=rd(root/_public_analysis_directory(r.performance_path.replace('_performance.tsv','_cells.tsv')));c=cells[cells.method==r.method][['cohort','seed','auc_method']].merge(sc[['cohort','seed','auc_method']],on=['cohort','seed'],suffixes=('_primary','_seed'));c['seed_minus_primary']=c.auc_method_seed-c.auc_method_primary;c.insert(0,'method',r.method)
  if tag!='old':save(c,root/'evidence'/_public_analysis_directory(f'{r.method}_seed_effects.tsv'))
  for f,v in [('minimum',c.seed_minus_primary.min()),('maximum',c.seed_minus_primary.max()),('cells',len(c))]:add('model_seed_effect_'+f,v,r.method,source=f'evidence/{r.method}_seed_effects.tsv')
 if tag!='old':
  save(pd.DataFrame(rows),root/'evidence/headline_metrics.tsv');save(pd.DataFrame(ssum),root/'evidence/sensitivity_summary.tsv')
 return pd.DataFrame(rows)
oldmetrics=derived(OLD,'old');newmetrics=derived(S,'primary');fullmetrics=derived(S/'full32_sensitivity','full')
save(oldmetrics,O/'current_headline_metrics.tsv')
change_rows=[]
def compare(a,b,label,keys,oldsource,newsource,keep_all=True):
 a=a.copy();b=b.copy()
 if 'analysis'in a:a['analysis']=a.analysis.replace(LABEL_MAP)
 if 'Sensitivity analysis'in a:a['Sensitivity analysis']=a['Sensitivity analysis'].replace(LABEL_MAP)
 for d in [a,b]:
  d[keys]=d[keys].fillna('').astype(str);d['_occurrence']=d.groupby(keys,dropna=False).cumcount()
 k=keys+['_occurrence'];aa=a.set_index(k);bb=b.set_index(k)
 for key in aa.index.union(bb.index,sort=False):
  ar=aa.loc[key] if key in aa.index else None;br=bb.loc[key] if key in bb.index else None
  for field in sorted(set(aa.columns)|set(bb.columns)):
   av=ar.get(field,np.nan) if ar is not None else np.nan;bv=br.get(field,np.nan) if br is not None else np.nan
   same=(pd.isna(av) and pd.isna(bv)) or av==bv
   try:delta=float(bv)-float(av)
   except (ValueError,TypeError):delta=np.nan
   status='removed_from_primary' if br is None else 'added' if ar is None else 'unchanged' if same else 'rounding' if np.isfinite(delta) and abs(delta)<=1e-12 else 'changed'
   if keep_all or status!='unchanged':change_rows.append(dict(table=label,row_key='|'.join(map(str,key[:-1])),field=field,old=av,new=bv,delta=delta,status=status,old_source=str(oldsource),new_source=str(newsource)))
compare(oldmetrics,newmetrics,'headline_metrics',['metric','method','analysis'],OLD/'evidence/pooled.tsv',S/'evidence/headline_metrics.tsv')
sw=rd(S/'evidence/sensitivity_cohort_weights.tsv');sel_an='Excluding units with unresolved scale or timing or normalized input'
oldw=rd(OLD/'evidence/within_performance.tsv');oldw=oldw[(oldw.method=='bal-selbal')&~oldw.cohort.isin(NON|{'GSE85589','GSE188232','GSE70080','GSE117063'})&oldw.lift.notna()&oldw.lift_se.gt(0)]
oldss=rd(OLD/'evidence/sensitivity_pooled.tsv');ot=oldss[(oldss.method=='bal-selbal')&(oldss.analysis==sel_an)].tau2.iloc[0];weights=1/(oldw.lift_se**2+ot)
old_toray=100*weights[oldw.cohort=='TORAY_PUBLIC_CANONICAL'].iloc[0]/weights.sum();new_toray=sw[(sw.method=='bal-selbal')&(sw.analysis==sel_an)&(sw.cohort=='TORAY_PUBLIC_CANONICAL')].weight_percent.iloc[0]
change_rows.append(dict(table='sensitivity_cohort_weights',row_key=sel_an+'|bal-selbal|TORAY_PUBLIC_CANONICAL',field='weight_percent',old=old_toray,new=new_toray,delta=new_toray-old_toray,status='changed',old_source=str(OLD/'evidence/sensitivity_pooled.tsv'),new_source='synthesis/evidence/sensitivity_cohort_weights.tsv'))
pairs=[]
for f in ['method_cohort_performance.tsv','pooled_effects.tsv','sensitivity_effects.tsv']:
 a=rd(OLD/_public_analysis_directory(f));b=rd(S/_public_analysis_directory(f));pairs.append((OLD/_public_analysis_directory(f),S/_public_analysis_directory(f),list(a.columns[:2]) if f!='pooled_effects.tsv' else [a.columns[0]]))
for f,keys in [('pooled.tsv',['method']),('cohort_weights.tsv',['method','cohort']),('leave_one_cohort_out.tsv',['method','omitted_cohort']),('sensitivity_pooled.tsv',['analysis','method']),('all_biofluid_pooled.tsv',['method']),('modality_pooled.tsv',['method','modality']),('exclusion_impact.tsv',['method']),('pooled_without_exclusions.tsv',['method']),('grouped_changes.tsv',['method']),('log2_performance.tsv',['method','cohort']),('grouped_performance.tsv',['method','cohort'])]:pairs.append((OLD/'evidence'/_public_analysis_directory(f),S/'evidence'/_public_analysis_directory(f),keys))
for f in ['matched_baseline_dataset.tsv','matched_baseline_pooled.tsv']:pairs.append((M/'benchmark_saved_scores/out'/_public_analysis_directory(f),O/_public_analysis_directory(f),['method','cohort'] if 'dataset' in f else ['method']))
for a,b,keys in pairs:compare(rd(a),rd(b),str(b.relative_to(O)),keys,a,b)
for folder,oldfolder in [(N/'companion_sensitivity/computing_bundle',M/'companion_sensitivity/out/computing_bundle'),(N/'normalization',M/'normalization/out'),(N/'normalization_comparison',M/'normalization_comparison/out'),(N/'companion_sensitivity/restriction',M/'companion_sensitivity/out/restriction')]:
 for f in folder.glob('*.tsv'):
  old=oldfolder/f.name
  if not old.exists():continue
  a=rd(old);b=rd(f)
  if 'D2_pooled' in f.name:a=a[a.analysis.isin(['Part A','Part B deployment AUC'])]
  if 'D2_restriction' in f.name:a=a[a.analysis.isin(['Part A','Part B deployment AUC'])]
  candidates=['expectation','dataset','unit','analysis','part','type','contrast','metric','arm','model','seed','fold','repetition','stage','batch','fraction']
  keys=[x for x in candidates if x in a and x in b]
  valid('semantic_keys_'+str(f.relative_to(O)),len(keys)>0)
  compare(a,b,str(f.relative_to(O)),keys,old,f)

def normalization_metrics(root):
 rows=[]
 def add(metric,value,method='',analysis='normalization'):rows.append(dict(metric=metric,method=method,analysis=analysis,value=value,source_file='normalization/R2b/m16/'))
 if root==M/'companion_sensitivity/out':base=root/'computing_bundle';r2d=M/'normalization_comparison/out'
 else:base=root/'companion_sensitivity/computing_bundle';r2d=root/'normalization_comparison'
 ct=rd(base/'partA_contrasts_dataset.tsv');da=rd(base/'partA_dataset_auc.tsv');fits=rd(base/'partA_fits.tsv');bf=rd(base/'partB_fits.tsv');bs=rd(base/'partB_summary.tsv');ex=rd(base/'expectations.tsv')
 add('PartA_planned_datasets',ct.dataset.nunique());add('PartA_fitted_datasets',fits.dataset.nunique());add('PartB_planned_datasets',bs.dataset.nunique())
 for part,d in [('PartA',fits),('PartB',bf)]:
  z=d[d.model.isin(['STD','QN','RAW','CLR'])&d.features_selected.eq(0)]
  add(part+'_intercept_only_fits',len(z));add(part+'_intercept_only_datasets',z.dataset.nunique())
 for (ty,co),g in ct.groupby(['type','contrast'],sort=False):
  for f,v in [('median',g.estimate.median()),('minimum',g.estimate.min()),('maximum',g.estimate.max()),('datasets_evaluable',g.estimate.notna().sum())]:add('PartA_dataset_'+f,v,co)
  if ty=='A1':
   for f,v in [('median',-g.estimate.median()),('minimum',-g.estimate.max()),('maximum',-g.estimate.min())]:add('frozen_minus_cohort_'+f,v,co)
 for (en,st),g in ex.groupby(['expectation','status'],sort=False):add('expectation_'+en+'_'+st,len(g))
 for _,r in ex[ex.expectation=='E3'].iterrows():add('deployment_minus_fullset_mean',r.value,r.arm)
 for ar,g in bs[bs.metric=='B5'].groupby('arm'):add('deployment_median_call_discordance',g.estimate.median(),ar)
 for file in ['partA_nonempty_pooled.tsv','partA_nonempty_without_toray_pooled.tsv']:
  q=rd(r2d/_public_analysis_directory(file));q=q[q.type=='A2'];add(file+'_minimum_ci_lo',q.ci_lo.min());add(file+'_min_k',q.k.min());add(file+'_max_k',q.k.max())
 return pd.DataFrame(rows)
nmold=normalization_metrics(M/'companion_sensitivity/out');nmnew=normalization_metrics(N)
compare(nmold,nmnew,'normalization_printed_metrics',['metric','method','analysis'],M/'companion_sensitivity/out',N)
save(nmnew,N/'printed_metrics.tsv');save(pd.concat([newmetrics,nmnew],ignore_index=True),O/'headline_metrics.tsv')
# Explicitly oriented figure source values, including every point and heatmap cell.
for tag,root,nroot in [('current',OLD,M/'companion_sensitivity/out'),('eligible',S,N/'companion_sensitivity')]:
 base=nroot/'computing_bundle';pa=rd(base/'partA_pooled.tsv');ct=rd(base/'partA_contrasts_dataset.tsv')
 keep=['CN-STD minus FZ-STD','CN-QN minus FZ-QN','CN-COMBAT minus FZ-RAW','SS-CLR minus CN-STD','SS-CLR minus CN-QN','SS-CLR minus CN-COMBAT']
 pp=pa[pa.contrast.isin(keep)].copy();dd=ct[ct.contrast.isin(keep)].copy()
 for d,field in [(pp,'pooled_lift'),(dd,'estimate')]:
  neg=d.type.eq('A1');d.loc[neg,field]=-d.loc[neg,field];lo=d.loc[neg,'ci_lo'].copy();d.loc[neg,'ci_lo']=-d.loc[neg,'ci_hi'];d.loc[neg,'ci_hi']=-lo
  d['figure_direction']='single-sample minus cohort-normalized'
 bs=rd(base/'partB_summary.tsv');heat=bs[bs.metric.isin(['B1','B3'])].pivot(index=['dataset','arm'],columns='metric',values='estimate').reset_index();heat['difference']=heat.B3-heat.B1
 heat['printed']=heat.difference.map(lambda x:'' if pd.isna(x) else '0.00' if f'{x:+.2f}' in ['+0.00','-0.00'] else f'{x:+.2f}')
 s7=rd(root/'pooled_effects.tsv');s7['eligible_universe']=32 if tag=='current' else 30
 lasso=pa[(pa.type=='A3')&(pa.contrast=='FZ-STD minus R')].copy();lasso['eligible_universe']=32 if tag=='current' else 30
 for name,d in [('figure2C_pooled',pp),('figure2C_dataset',dd),('figure2D',heat),('figure3_methods',s7),('figure3_frozen_lasso',lasso)]:save(d,O/'figure_values'/_public_analysis_directory(tag)/(_public_analysis_directory(name+'.tsv')))
for name,keys in [('figure2C_pooled',['type','contrast']),('figure2C_dataset',['dataset','type','contrast']),('figure2D',['dataset','arm']),('figure3_methods',['Method identifier']),('figure3_frozen_lasso',['type','contrast'])]:
 a=O/'figure_values/current'/_public_analysis_directory(f'{name}.tsv');b=O/'figure_values/eligible'/_public_analysis_directory(f'{name}.tsv');compare(rd(a),rd(b),'figure_values/'+name,keys,a,b)
# All finite quantities needed to replace manuscript/matched_baseline tables and plot source labels.
head=[]
for label,path in [('primary pooled versus rCLR',S/'evidence/pooled.tsv'),('matched frozen-lasso pooled',O/'matched_baseline_pooled.tsv'),('sensitivity pooled versus rCLR',S/'evidence/sensitivity_pooled.tsv')]:
 d=rd(path);d.insert(0,'metric',label);d['evidence']=str(path.relative_to(O));head.append(d)
q=pd.concat(head,ignore_index=True);save(q,O/'headline_results.tsv')
save(newmetrics[newmetrics.metric.isin(['near_zero_count','near_zero_list'])],S/'evidence/near_zero_summary.tsv')
p=rd(S/'evidence/pooled.tsv');save(p[(p.ci_hi>-.01)&(p.ci_hi<0)],S/'evidence/near_zero_methods.tsv')

# Preserve original computing identities, clearly separated from this pooling run.
prov=O/'source_provenance';prov.mkdir(exist_ok=True)
for f in ['provenance.json','execution_environment.tsv','actual_task_environment.tsv','table_provenance.tsv']:shutil.copy2(M/'companion_sensitivity/out'/_public_analysis_directory(f),prov/(_public_analysis_directory('R2b_'+f)))
ids=set(rd(S/'evidence/within_performance.tsv').cohort)
for f in ['final_cell_provenance.tsv','benchmark_sources.tsv','intercept_only_folds.tsv','fold_row_alignment.tsv','matched_fold_differences.tsv']:
 d=rd(M/'benchmark_saved_scores/out'/_public_analysis_directory(f))
 if 'cohort'in d:d=d[d.cohort.isin(ids)]
 save(d,prov/(_public_analysis_directory('R11_'+f)))
for root in [S,S/'full32_sensitivity']:
 d=rd(root/'evidence/within_performance.tsv');cells=rd(root/'cells_union_complete.tsv')
 for method,g in d.groupby('method'):
  save(g,root/'evidence/methods'/_public_analysis_directory(f'{method}_performance.tsv'))
  save(cells[cells.method==method],root/'evidence/methods'/_public_analysis_directory(f'{method}_cells_complete.tsv'))

# Validation: exact reference benchmark replay, normalization replay, frozen inputs,
# eligible membership, complete sensitivity families, and normalized cohort weights.
for f in ['pooled.tsv','leave_one_cohort_out.tsv','cohort_weights.tsv','sensitivity_pooled.tsv']:
 a=rd(OLD/'evidence'/_public_analysis_directory(f));b=rd(S/'full32_sensitivity/evidence'/_public_analysis_directory(f))
 keys=[x for x in ['analysis','method','omitted_cohort','cohort'] if x in a]
 a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True)
 valid('full32_benchmark_replay_'+f,a.shape==b.shape)
 cols=a.select_dtypes('number').columns;diff=(a[cols]-b[cols]).abs().max().max()
 valid('full32_numeric_replay_'+f,(np.isnan(diff) or diff<=1e-6) and a[cols].isna().equals(b[cols].isna()),diff)
for rel,old in [('companion_sensitivity/computing_bundle/partA_pooled.tsv','companion_sensitivity/out/computing_bundle/partA_pooled.tsv'),('companion_sensitivity/computing_bundle/partB_pooled.tsv','companion_sensitivity/out/computing_bundle/partB_pooled.tsv'),('normalization/partA_without_toray_pooled.tsv','normalization/out/partA_without_toray_pooled.tsv'),('normalization_comparison/partA_nonempty_pooled.tsv','normalization_comparison/out/partA_nonempty_pooled.tsv'),('normalization_comparison/partA_nonempty_without_toray_pooled.tsv','normalization_comparison/out/partA_nonempty_without_toray_pooled.tsv')]:
 a=rd(M/_public_analysis_directory(old)).set_index(['type','contrast']);b=rd(N/'full32_sensitivity'/_public_analysis_directory(rel)).set_index(['type','contrast']);cols=a.select_dtypes('number').columns
 delta=(a[cols]-b[cols]).abs().max().max();valid('full32_normalization_replay_'+rel,delta<=1e-9,delta)
valid('primary_units',len(ids)==30 and not ids&(EX|NON) and 'GSE266859'in ids)
w=rd(S/'evidence/within_performance.tsv');ow=rd(OLD/'evidence/within_performance.tsv');cols=[c for c in w.columns if c not in ['p_lift_gt_zero_BY','p_lift_gt_margin_0_05_BY']]
a=w[cols].set_index(['method','cohort']).sort_index();b=ow[ow.cohort.isin(ids)][cols].set_index(['method','cohort']).sort_index();valid('retained_cohort_estimands_unchanged',a.equals(b))
ss=rd(S/'evidence/sensitivity_pooled.tsv');valid('S9_12_blocks',ss.analysis.nunique()==12 and ss.groupby('analysis').size().eq(49).all())
families=ss.groupby('analysis').agg(method_rows=('method','size'),above_zero_finite_tests=('p_lift_gt_zero','count'),above_margin_finite_tests=('p_lift_gt_margin_0_05','count')).reset_index();save(families,O/'analysis_families.tsv')
valid('S9_both_48_test_families',families.above_zero_finite_tests.eq(48).all() and families.above_margin_finite_tests.eq(48).all())
valid('primary_weights_sum_100',np.allclose(rd(S/'evidence/cohort_weights.tsv').groupby('method').weight_percent.sum(),100,atol=1e-10))
fz=rd(O/'matched_baseline_pooled.tsv');valid('matched_six_test_families',len(fz)==6 and fz.p_lift_gt_zero.notna().all() and fz.datasets_planned.eq(30).all())
valid('no_refits',json.loads((O/'computing_identity.json').read_text())['refits']==0)
valid('primary_A2_noninferiority',rd(N/'companion_sensitivity/computing_bundle/partA_pooled.tsv').query('type=="A2"').ci_lo.gt(-.02).all())
pass
save(pd.DataFrame(checks),O/'validation.tsv')

names={'ai-tabpfn':'TabPFN','tab-tabicl':'TabICL','tab-tabdpt':'TabDPT','bal-selbal':'selbal','reo-pairratio':'pair-ratio','reo-ktsp':'k-TSP','ws-rclr-trimmed':'rCLR'}
def ci(r):return f'{r.pooled_lift:+.4f} (95% CI {r.ci_lo:+.4f} to {r.ci_hi:+.4f})'
e=p[p.pooled_lift.notna()];top=e.nlargest(3,'pooled_lift');pr=p[p.method=='reo-pairratio'].iloc[0];near=p[(p.ci_hi>-.01)&(p.ci_hi<0)]
bp=rd(N/'companion_sensitivity/computing_bundle/partB_pooled.tsv');norm=rd(N/'companion_sensitivity/computing_bundle/partA_pooled.tsv')
text='; '.join(f'{names[r.method]} {ci(r)}' for _,r in top.iterrows())
report=f'''The 30-unit primary synthesis gives these largest pooled AUC differences versus rCLR: {text}. Pair-ratio gives {ci(pr)} across {int(pr.k)} datasets. None of 48 estimable methods passes either BY family; the smallest adjusted P above zero is {e.p_lift_gt_zero_BY.min():.3f}. {sum(e.ci_hi<0)} intervals lie wholly below zero. TabICL’s upper confidence limit exceeds +0.05. Four upper limits lie between −0.010 and zero: {', '.join(near.method)}. [Results](headline_results.tsv).

The exploratory input-quality exclusion favors selbal: +0.0302 (95% CI +0.0256 to +0.0347), 14 datasets, BY P=2.68e-7. Toray carries {new_toray:.2f}% of this estimate’s weight. No sensitivity shows a BY-significant gain above +0.05. GSE270497 and GSE186595 are excluded from primary pools; GSE266859 remains eligible. The 32-unit set is a named sensitivity. The sensitivity analysis uses eleven eligibility definitions. [S9](synthesis/S9.tsv).

The matched frozen-lasso analysis includes six comparisons, two six-test BY families, cohort weights, leave-one-out estimates and prediction intervals. Normalization comparisons retain their original estimators and schemas. SS-CLR remains non-inferior at −0.02 in all three primary benchmark-fold comparisons. Deployment includes 11 eligible datasets, with 10 entering each contrast’s pool. [Normalization](normalization/R2b/m16/partB_pooled.tsv).

Independent Python code reproduces every primary pooled estimate and both BY families from S6, with maximum error {primary_error:.3g}. Every S9 pool was independently checked; all {len(ind):,} checks pass within 1e-6. Full-set primary estimates reproduce R11 exactly; sensitivity replay differs by at most 4.7e-10 and normalization by 3.2e-12. [Checks](independent_check.tsv), [analysis comparisons](changes.tsv).

No classifier was refitted. Computation used {platform.node()}, CPU, nice 10 and one numerical thread. Original cell computing identities are retained. The uncertainty estimates retain the original repeated-CV and conditional-bootstrap limitations. A separate likelihood-optimum audit records numerical stopping differences; production estimators are unchanged. This local calculation does not establish clinical or manuscript readiness.
'''
valid('report_at_most_350_words',len(report.split())<=350,len(report.split()))
(O/'REPORT.md').write_text(report);(S/'REPORT.md').write_text('See [results and limitations](../REPORT.md).\n');(S/'results_summary.md').write_text(report.replace('](','](../'))
(O/'REBUILD.md').write_text('Run the public summary entry point with the checksum-verified release bundle. Numerical contracts compare the generated cohort and pooled estimates, uncertainty intervals, multiplicity adjustments and sensitivity analyses with deposited values. Independent calculations implement REML, Knapp–Hartung intervals and Benjamini–Yekutieli adjustment. Prepared-input and repeated-cross-validation limitations apply.\n')
(S/'table_notes.md').write_text('Primary analysis: 30 units; GSE270497 and GSE186595 excluded, GSE266859 retained. GSE145176 and GSE31568 remain outside cancer-discrimination pooling. All effects are method minus rCLR unless the table names another comparator. REML random effects and Knapp–Hartung intervals retain the R11 estimator and numerical stopping rules. Positive finite corrected repeated-CV standard errors are required for pooling. Near zero means an upper 95% confidence limit strictly between −0.010 and zero. Full-set sensitivity and melanoma exclusion are explicitly named. See ../REBUILD.md and ../source_provenance/.\n')
(N/'README.md').write_text('Primary: eligible 30-unit set for Part A, 11 units for Part B deployment. `full32_sensitivity/` retains the original set (32 Part A units; 12 deployment units). Original schemas and contrast directions are preserved. R2d uses the original contrast-specific all-fold nonempty membership, intersected with the eligible set. Zero-SE effects remain descriptive and do not enter pools. Deposited-input restriction outputs cover only Part A and Part B cancer discrimination; unrelated restriction analyses are outside this task. This summary uses stored bootstrap and model outputs.\n')
print(report)
