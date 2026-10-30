from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import pandas as pd
import numpy as np
M=Path('./analyses')
NON={'GSE145176','GSE31568'}
LABEL_MAP={}
def rd(p):return pd.read_csv(p,sep='\t')
def save(d,p):p.parent.mkdir(parents=True,exist_ok=True);_public_table(d).to_csv(p,sep='\t',index=False,na_rep='')
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
