from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import pandas as pd
import numpy as np
import shutil
B=Path(__file__).resolve().parent.parent;O=B/"out";E=O/"evidence"
spec={"stats2":["item1_scoring_conditions_scoring_conditions.tsv","item2_intercept_only_sensitivity.tsv"]}
def read(p):return pd.read_csv(p,sep="\t")
for f in spec["stats2"]:shutil.copy2(B/"results/statsfinal6"/_public_analysis_directory(f),O/_public_analysis_directory(f))
a=read(O/_public_analysis_directory(spec['stats2'][0]));b=read(O/_public_analysis_directory(spec['stats2'][1]));checks=[]
def check(scope,key,x,y,tol=1e-10):
    ok=(pd.isna(x) and pd.isna(y)) or (np.isfinite(x) and np.isfinite(y) and abs(x-y)<=tol)
    checks.append(dict(scope=scope,key=key,observed=x,expected=y,absolute_difference=abs(x-y),pass_check=ok))
M=B.parent;t=read(M/'companion_sensitivity/out/reassay_design/scoring_conditions_rows.tsv');t=t[t.dataset.eq('GSE268524')]
aa=a[a.scope.eq('person_mean_across_splits')]
for _,r in t.iterrows():
    z=aa[aa.scorer.eq(r.scorer)&aa.condition.eq(aa.primary_table2_condition)].iloc[0]
    for k,k0 in [('estimate','AUC_estimate'),('ci_lo','AUC_ci_low'),('ci_hi','AUC_ci_high')]:check('Table2 primary reconstruction',r.scorer+' '+k,z[k],r[k0])
ref=read(M/'companion_sensitivity/out/computing_bundle/partB_summary.tsv');bb=b[b.row_type.eq('summary')&b.scenario.eq('primary')&b.scope.eq('mean_retained_splits')&~b.metric.eq('any_changed_call')]
z=bb.merge(ref,on=['dataset','arm','metric'],suffixes=('_new','_old'),validate='one_to_one');assert len(z)==len(bb)
for _,r in z.iterrows():
    for k in ['estimate','ci_lo','ci_hi']:check('case mix primary reconstruction',r.dataset+' '+r.arm+' '+r.metric+' '+k,r[k+'_new'],r[k+'_old'])
ref=read(M/'benchmark_sensitivity_summary/out/normalization/companion_sensitivity/computing_bundle/partB_pooled.tsv');bb=b[b.row_type.eq('pooled_deployment_auc_difference')&b.scenario.eq('primary')]
z=bb.merge(ref,on='contrast',suffixes=('_new','_old'),validate='one_to_one')
for _,r in z.iterrows():
    check('deployment primary pool',r.contrast,r.estimate,r.pooled_lift,1e-8)
    for k in ['ci_lo','ci_hi']:check('deployment primary pool',r.contrast+' '+k,r[k+'_new'],r[k+'_old'],1e-8)
for typ in ['summary','dataset_deployment_auc_difference','pooled_deployment_auc_difference']:
    q=b[b.row_type.eq(typ)&b.scenario.isin(['primary','neutral_constant_no_positive'])]
    q=q[q.metric.isin(['B1','B2_10','B2_50','B2_90','B3','B5','any_changed_call','B3_difference'])]
    keys=['dataset','arm','scope','metric','contrast'];q=q.fillna({k:'' for k in keys});left=q[q.scenario.eq('primary')];right=q[q.scenario.eq('neutral_constant_no_positive')];zz=left.merge(right,on=keys,suffixes=('_p','_n'),validate='one_to_one')
    for _,r in zz.iterrows():check('neutral AUC and discordance unchanged',typ+' '+'|'.join(str(r[k]) for k in keys),r.estimate_p,r.estimate_n)
audit=pd.DataFrame(checks);_public_table(audit).to_csv(E/'stats_validation.tsv',sep='\t',index=False);assert audit.pass_check.all(),audit[~audit.pass_check].to_string()

print(audit.groupby("scope").pass_check.agg(["count","all"]).to_string())
