"""Check neutral calls against each fitted model, including pooled deployment."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import pandas as pd
import numpy as np
B=Path(__file__).resolve().parent.parent
d=pd.read_csv(B/'out/item2_intercept_only_sensitivity.tsv',sep='\t',low_memory=False)
fits=pd.read_csv(B/'payload/case_fits.tsv',sep='\t')
arm_model={'FZ-STD':'STD','CN-STD':'STD','FZ-QN':'QN','CN-QN':'QN','FZ-RAW':'RAW','CN-COMBAT':'RAW','SS-CLR':'CLR'}
empty={(r.dataset,r.model,int(r.repetition)) for r in fits.itertuples() if r.features_selected==0 and r.model in arm_model.values()}
q=d[d.row_type.eq('summary')&d.scope.str.startswith('split_',na=False)&d.metric.str.startswith('B4_',na=False)]
keys=['dataset','arm','scope','metric']
z=q[q.scenario.eq('primary')].merge(q[q.scenario.eq('neutral_constant_no_positive')],on=keys,suffixes=('_p','_n'),validate='one_to_one')
checks=[]
for r in z.itertuples():
    constant=(r.dataset,arm_model.get(r.arm),int(r.scope.removeprefix('split_'))) in empty
    for col in ['estimate','ci_lo','ci_hi']:
        expected=float('specificity' in r.metric) if constant else getattr(r,col+'_p')
        observed=getattr(r,col+'_n')
        ok=(pd.isna(observed) and pd.isna(expected)) or (np.isfinite(observed) and np.isfinite(expected) and abs(observed-expected)<1e-12)
        checks.append(dict(dataset=r.dataset,arm=r.arm,scope=r.scope,metric=r.metric,quantity=col,intercept_only=constant,observed=observed,expected=expected,pass_check=ok))
qq=d[d.row_type.eq('summary')&d.metric.str.startswith('B4_',na=False)]
zz=qq[qq.scenario.eq('primary')].merge(qq[qq.scenario.eq('neutral_constant_no_positive')],on=keys,suffixes=('_p','_n'),validate='one_to_one')
for r in zz.itertuples():
    checks.append(dict(dataset=r.dataset,arm=r.arm,scope=r.scope,metric=r.metric,quantity='bootstrap_finite_draws',intercept_only=np.nan,observed=r.bootstrap_finite_draws_n,expected=r.bootstrap_finite_draws_p,pass_check=r.bootstrap_finite_draws_n==r.bootstrap_finite_draws_p))
a=pd.DataFrame(checks);assert len(a)>0
_public_table(a).to_csv(B/'out/evidence/neutral_operating_validation.tsv',sep='\t',index=False)
assert a.pass_check.all(),a[~a.pass_check].to_string(index=False)
print(len(a),'neutral operating-point checks passed')
