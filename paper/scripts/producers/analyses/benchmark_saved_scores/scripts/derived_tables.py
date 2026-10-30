#!/usr/bin/env python3
"""Recompute the pre-existing headline and sensitivity summaries, preserving schemas."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import os,csv,json,math
import pandas as pd
import numpy as np
B=Path(__file__).resolve().parents[1];O=B/'out/synthesis'
def rd(p):return pd.read_csv(p,sep='\t')
def save(x,p):_public_table(x).to_csv(p,sep='\t',index=False,na_rep='')
w=rd(O/'evidence/within_performance.tsv');p=rd(O/'evidence/pooled.tsv');ss=rd(O/'evidence/sensitivity_pooled.tsv');weights=rd(O/'evidence/cohort_weights.tsv')
a=w[(w.method=='reo-pairratio')&w.all_seeds_eligible&w.auc_method.notna()&~w.cohort.isin(['GSE145176','GSE31568'])]
b=w[(w.method=='ws-rclr-trimmed')&w.all_seeds_eligible&w.auc_method.notna()]
shared=a.merge(b[['cohort','auc_method']],on='cohort',suffixes=('_pairratio','_reference'));shared['paired_auc_difference']=shared.auc_method_pairratio-shared.auc_method_reference
shared['within_source_data_row_pairratio']=shared.cohort.map(w[w.method=='reo-pairratio'].reset_index().set_index('cohort')['index'])+1
shared['within_source_data_row_reference']=shared.cohort.map(w[w.method=='ws-rclr-trimmed'].reset_index().set_index('cohort')['index'])+1
# interval/diagnostic audit fields are supplied in within_performance instead.
old_shared=rd(B/'reference/evidence/shared_cohort_auc.tsv')
shared=shared.reindex(columns=old_shared.columns)
save(shared,O/'evidence/shared_cohort_auc.tsv')
metrics=[]
def metric(name,value,source,rule,method='',cohort='',column=''):
 metrics.append(dict(metric=name,method=method,cohort=cohort,value=value,source_file=source,source_rows_or_filter=rule,source_column=column))
e=p[p.pooled_lift.notna()];top=e.loc[e.pooled_lift.idxmax()];uc=e.loc[e.ci_hi.idxmax()];pi=e.loc[e.pi_hi.idxmax()]
for name,value,source,rule in [
 ('planned_method_cohort_cells',len(w),'synthesis/method_cohort_performance.tsv','all rows'),('evaluable_method_cohort_cells',int((w.all_seeds_eligible&w.auc_method.notna()).sum()),'synthesis/method_cohort_performance.tsv','complete and finite AUC'),('pending_method_seed_cohort_cells',0,'synthesis/evidence/pending_cells.tsv','all rows'),('estimable_methods',len(e),'synthesis/pooled_effects.tsv','finite pooled estimates'),('BY_above_zero',int((e.p_lift_gt_zero_BY<.05).sum()),'synthesis/pooled_effects.tsv','BY-adjusted P < 0.05'),('BY_above_005',int((e.p_lift_gt_margin_0_05_BY<.05).sum()),'synthesis/pooled_effects.tsv','BY-adjusted P above 0.05 margin < 0.05'),('CI_wholly_below_zero',int((e.ci_hi<0).sum()),'synthesis/pooled_effects.tsv','upper CI < 0'),('CI_wholly_below_zero_methods',';'.join(e[e.ci_hi<0].method),'synthesis/pooled_effects.tsv','upper CI < 0'),('CI_wholly_above_zero',int((e.ci_lo>0).sum()),'synthesis/pooled_effects.tsv','lower CI > 0'),('CI_upper_below_005',int((e.ci_hi<.05).sum()),'synthesis/pooled_effects.tsv','upper CI < 0.05'),('shared_pairratio_reference_cohorts',len(shared),'synthesis/evidence/shared_cohort_auc.tsv','all rows'),('pairratio_shared_median_AUC',shared.auc_method_pairratio.median(),'synthesis/evidence/shared_cohort_auc.tsv','all rows'),('reference_shared_median_AUC',shared.auc_method_reference.median(),'synthesis/evidence/shared_cohort_auc.tsv','all rows'),('shared_median_paired_AUC_difference',shared.paired_auc_difference.median(),'synthesis/evidence/shared_cohort_auc.tsv','all rows')]:metric(name,value,source,rule)
for name,row,col in [('largest_pooled_gain',top,'pooled_lift'),('largest_CI_upper',uc,'ci_hi'),('largest_PI_upper',pi,'pi_hi')]:metric(name,row[col],'synthesis/evidence/pooled.tsv',str(row.name+1),row.method,column=col)
for _,r in p.iterrows():
 for col in ['k','pooled_lift','ci_lo','ci_hi','pi_lo','pi_hi','p_lift_gt_zero_BY','p_lift_gt_margin_0_05_BY','I2']:metric(col,r[col],'synthesis/evidence/pooled.tsv',str(r.name+1),r.method,column=col)
for _,r in weights[weights.cohort=='TORAY_PUBLIC_CANONICAL'].iterrows():metric('Toray_weight_percent',r.weight_percent,'synthesis/evidence/cohort_weights.tsv',str(r.name+1),r.method,r.cohort,'weight_percent')
for _,r in ss[ss.analysis=='Excluding the Toray unit'].iterrows():
 for col in ['pooled_lift','ci_lo','ci_hi','pi_lo','pi_hi','I2']:metric('leave_Toray_out_'+col,r[col],'synthesis/evidence/sensitivity_pooled.tsv',str(r.name+1),r.method,column=col)
sourcegrid=rd(O/'chunk_sources.tsv');sourcegrid=sourcegrid[~sourcegrid.chunk_id.str.endswith('__gse83977_corrected')]
meta=w.drop_duplicates('cohort');toray=meta[meta.cohort=='TORAY_PUBLIC_CANONICAL'].iloc[0]
refset=set(b[~b.cohort.isin(['GSE145176','GSE31568'])].cohort)
oldpool=rd(O/'evidence/v6/pooled.tsv');oldmethods=set(oldpool[oldpool.pooled_lift.notna()].method)
groupchange=rd(O/'evidence/grouped_changes.tsv');maxgroup=groupchange.loc[groupchange.delta_pooled_lift.abs().idxmax()]
extras={'available_chunks':int((sourcegrid.status=='available').sum()),'pending_chunks':int((sourcegrid.status=='pending').sum()),'reference_evaluable_pairratio_unavailable':len(refset-set(shared.cohort)),'newly_estimable_methods':len(set(e.method)-oldmethods),'no_longer_estimable_methods':len(oldmethods-set(e.method)),'maximum_absolute_grouped_tuning_change':abs(maxgroup.delta_pooled_lift),'datasets':len(meta),'profiles':int(meta.n_profiles_inventory.sum()),'recorded_groups':int(meta.n_groups_inventory.sum()),'source_accessions':int(meta.n_source_accessions.sum()),'Toray_profiles':int(toray.n_profiles_inventory),'Toray_corpus_percentage':100*toray.n_profiles_inventory/meta.n_profiles_inventory.sum()}
for name,value in extras.items():
 src='synthesis/chunk_sources.tsv' if name in ['available_chunks','pending_chunks'] else 'synthesis/evidence/grouped_changes.tsv' if name=='maximum_absolute_grouped_tuning_change' else 'synthesis/evidence/v6/pooled.tsv; synthesis/evidence/pooled.tsv' if name in ['newly_estimable_methods','no_longer_estimable_methods'] else 'synthesis/evidence/within_performance.tsv'
 metric(name,value,src,'same original metric definition')
old_s6=rd(B/'reference/method_cohort_performance.tsv');new_s6=rd(O/'method_cohort_performance.tsv')
renames={'Corrected repeated-CV standard error':'Nadeau-Bengio corrected repeated-CV standard error (inferential)','Biological-group influence standard error (diagnostic)':'Person-collapsed Nadeau-Bengio standard error (diagnostic)'}
old_s6=old_s6.rename(columns=renames);keys=list(old_s6.columns[:2]);old_s6=old_s6.set_index(keys).sort_index();new_s6=new_s6.set_index(keys).sort_index()
common=list(old_s6.columns);same=(old_s6.eq(new_s6[common])|(old_s6.isna()&new_s6[common].isna())).all(axis=1)
metric('S6_changed_rows',int((~same).sum()),'out/changes.tsv','R1b baseline; changed values in existing columns')
old_s7=rd(B/'reference/pooled_effects.tsv').set_index('Method identifier').sort_index();new_s7=rd(O/'pooled_effects.tsv').set_index('Method identifier').sort_index()
changed_cols=[]
for col in new_s7:
 if col not in old_s7 or not (old_s7[col].eq(new_s7[col])|(old_s7[col].isna()&new_s7[col].isna())).all():changed_cols.append(col)
metric('S7_changed_columns',len(changed_cols),'out/changes.tsv','R1b baseline; changed or added S7 columns')
h=pd.DataFrame(metrics);save(h,B/'out/headline.tsv')
legacy=rd(B/'reference/evidence/headline_metrics.tsv');lookup={r['metric']:r for r in metrics if not r['method']}
for name in ['largest_pooled_gain','largest_CI_upper','largest_PI_upper']:lookup[name]=next(r for r in metrics if r['metric']==name)
for i,r in legacy.iterrows():
 assert r.metric in lookup,r.metric
 q=lookup[r.metric]
 for col in ['value','source_file','source_rows_or_filter']:legacy.loc[i,col]=q[col].replace('synthesis/','').replace('out/','../out/') if col=='source_file' else q[col]
save(legacy,O/'evidence/headline_metrics.tsv')
summary=[]
for label,g in ss.groupby('analysis',sort=False):
 a=g[g.pooled_lift.notna()];mx=a.loc[a.ci_hi.idxmax()];rr=a[a.method=='reo-pairratio']
 summary.append(dict(analysis=label,estimable=len(a),BY_above_zero=int((a.p_lift_gt_zero_BY<.05).sum()),BY_above_005=int((a.p_lift_gt_margin_0_05_BY<.05).sum()),CI_wholly_below_zero=int((a.ci_hi<0).sum()),CI_upper_above_005=int((a.ci_hi>.05).sum()),methods_CI_upper_above_005=';'.join(a[a.ci_hi>.05].method),maximum_CI_upper=mx.ci_hi,maximum_CI_upper_method=mx.method,pairratio_lift=rr.pooled_lift.iloc[0] if len(rr) else np.nan,pairratio_k=rr.k.iloc[0] if len(rr) else np.nan,S9_data_rows=';'.join(str(i+1) for i in g.index)))
summary=pd.DataFrame(summary);assert list(summary)==list(rd(B/'reference/evidence/sensitivity_summary.tsv'));save(summary,O/'evidence/sensitivity_summary.tsv')
# Tag each copied ledger explicitly. Original input/output manifests are retained
active={'within_performance.tsv','pooled.tsv','cohort_weights.tsv','sensitivity_summary.tsv','headline_metrics.tsv','shared_cohort_auc.tsv','validation_checks.tsv','pooled_from_cohort_check.tsv','pending_cells.tsv','within_linear.tsv','log2_performance.tsv','sensitivity_pooled.tsv','leave_one_cohort_out.tsv','within_grouped.tsv','grouped_performance.tsv','grouped_changes.tsv','v6_comparison_summary.tsv','eligibility_transitions.tsv','cells_original_chunks.tsv','all_biofluid_pooled.tsv','modality_pooled.tsv','comparison_method_scope.tsv','within_performance_without_exclusions.tsv','pooled_without_exclusions.tsv','cohort_weights_without_exclusions.tsv','exclusion_impact.tsv','exclusion_checks.tsv','exclusion_fold_counts.tsv','model_warnings.tsv','reml_diagnostics.tsv','foundation_seed_sensitivities.tsv','evidence_status.tsv'}
status=[dict(file=f.name,status='recomputed_R9' if f.name in active else 'preserved_R1b_historical_audit') for f in sorted((O/'evidence').glob('*')) if f.is_file()]
save(pd.DataFrame(status),O/'evidence/evidence_status.tsv')
print('Recomputed',len(legacy),'headline metrics and',len(summary),'sensitivity summaries')
