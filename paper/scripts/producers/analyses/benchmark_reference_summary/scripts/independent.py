"""Independent REML likelihood and KH/BY calculations from the delivered method_cohort_performance.

No R code, production pooling function or fitted variance estimate is imported.
"""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import t, chi2

B=Path(__file__).resolve().parents[1]; O=B/'out'; S=O/'synthesis'
def read(p): return pd.read_csv(p,sep='\t')
def save(d,p): _public_table(d).to_csv(p,sep='\t',index=False,na_rep='')
def own_meta(y,se,match_stopping=True):
    good=np.isfinite(y)&np.isfinite(se)&(se>0)
    y=np.asarray(y[good],float); v=np.asarray(se[good],float)**2; k=len(y)
    fields='pooled_lift pooled_se ci_lo ci_hi pi_lo pi_hi tau2 I2 Q Q_p df_knha p_lift_gt_zero p_lift_gt_margin_0_05'.split()
    if k<5: return dict(k=k,**{f:np.nan for f in fields})
    scale=max(np.var(y,ddof=1),np.median(v),1e-14)
    v0=v/scale; y0=y/np.sqrt(scale)
    def nll(a):
        w=1/(v0+a); mu=np.dot(w,y0)/sum(w)
        return sum(np.log(v0+a))+np.log(sum(w))+np.dot(w,(y0-mu)**2)
    def score(a):
        w=1/(v0+a); mu=np.dot(w,y0)/sum(w)
        return sum(w)-np.dot(w,w)/sum(w)-np.dot(w*w,(y0-mu)**2)
    hi=max(1.,4*np.var(y0,ddof=1))
    while nll(hi*2)<nll(hi): hi*=2
    grid=np.r_[0,np.geomspace(1e-15,hi*2,1800)]
    sc=np.array([score(a) for a in grid])
    candidates=[0.,hi*2]
    for j in np.where(sc[:-1]*sc[1:]<0)[0]:
        candidates.append(brentq(score,grid[j],grid[j+1],xtol=1e-15,rtol=1e-14))
    a=min(candidates,key=nll); optimal_tau=a*scale; tau=optimal_tau
    if match_stopping:
        # Independently implement the documented metafor Fisher-scoring update
        # and the benchmark_saved_scores evaluate/acceptance thresholds. This reproduces the numerical
        # estimator being checked, without reading reported tau2 or calling R.
        controls=[(1.,1e-5,100),(.5,1e-10,10000),(.25,1e-12,20000),(.125,1e-12,40000)]
        for attempt,(step,tol,maxiter) in enumerate(controls):
            candidate=max(0.,np.var(y,ddof=1)-np.mean(v)); converged=False
            for iteration in range(maxiter+1):
                ww=1/(v+candidate); sw=sum(ww); sw2=np.dot(ww,ww)
                mu=np.dot(ww,y)/sw
                trp=sw-sw2/sw; trpp=sw2-2*sum(ww**3)/sw+(sw2/sw)**2
                adj=step*(np.dot(ww*ww,(y-mu)**2)-trp)/trpp
                while candidate+adj<0: adj/=2
                previous=candidate;candidate+=adj
                if abs(candidate-previous)<=tol:
                    converged=iteration<maxiter;break
            if not converged:continue
            bad_warning=(nll(candidate/scale)-nll(0))/2>np.finfo(float).eps**.25 and candidate>tol
            if bad_warning:candidate=0.
            good=abs(candidate-optimal_tau)<=max(1e-7,1e-3*optimal_tau) and (nll(candidate/scale)-nll(a))/2<=1e-6
            if attempt==0 and (candidate<=1e-12 or bad_warning or not good):continue
            if bad_warning or not good:continue
            tau=candidate;break
    w=1/(v+tau)
    mu=np.dot(w,y)/sum(w); df=k-1
    khse=np.sqrt(np.dot(w,(y-mu)**2)/df/sum(w)); crit=t.ppf(.975,df)
    wf=1/v; m0=np.dot(wf,y)/sum(wf); Q=np.dot(wf,(y-m0)**2)
    vt=df/(sum(wf)-np.dot(wf,wf)/sum(wf))
    return dict(k=k,pooled_lift=mu,pooled_se=khse,ci_lo=mu-crit*khse,ci_hi=mu+crit*khse,
                pi_lo=mu-crit*np.sqrt(khse**2+tau),pi_hi=mu+crit*np.sqrt(khse**2+tau),tau2=tau,
                I2=100*tau/(tau+vt),Q=Q,Q_p=chi2.sf(Q,df),df_knha=df,
                p_lift_gt_zero=t.sf(mu/khse,df),p_lift_gt_margin_0_05=t.sf((mu-.05)/khse,df))
def own_by(p):
    p=np.array(p,float); out=np.full(len(p),np.nan); ix=np.flatnonzero(np.isfinite(p)); n=len(ix)
    order=ix[np.argsort(p[ix],kind='stable')]
    scaled=p[order]*n*sum(1/np.arange(1,n+1))/np.arange(1,n+1)
    out[order]=np.minimum(1,np.minimum.accumulate(scaled[::-1])[::-1]);return out

checks=[]
def check(scope,key,field,a,b,tol=1e-6,required=True):
    both=pd.isna(a) and pd.isna(b); err=abs(a-b) if np.isfinite(a) and np.isfinite(b) else np.nan
    checks.append(dict(scope=scope,method=key,field=field,recomputed=a,reported=b,absolute_error=err,tolerance=tol,required=required,pass_check=both or (np.isfinite(err) and err<=tol)))
s6=read(S/'method_cohort_performance.tsv'); assert s6.iloc[:,1].nunique()==30
assert not s6.iloc[:,1].isin(['GSE270497','GSE186595','GSE145176','GSE31568']).any()
rows=[]
for method,g in s6.groupby(s6.columns[0],sort=True):
    rows.append(dict(method=method,**own_meta(g.iloc[:,7].to_numpy(),g.iloc[:,10].to_numpy())))
ind=pd.DataFrame(rows)
for h in ['p_lift_gt_zero','p_lift_gt_margin_0_05']:ind[h+'_BY']=own_by(ind[h])
p=read(S/'evidence/pooled.tsv').set_index('method')
for _,r in ind.iterrows():
    for f in ind.columns.drop('method'):
        check('Independent primary S6 REML/KH/BY',r.method,f,r[f],p.loc[r.method,f])
save(ind,O/'independent_primary.tsv')
optimal=pd.DataFrame([dict(method=m,**own_meta(g.iloc[:,7].to_numpy(),g.iloc[:,10].to_numpy(),False)) for m,g in s6.groupby(s6.columns[0],sort=True)])
save(optimal,O/'independent_global_optimum.tsv')
audit=[]
for _,r in optimal.iterrows():
    for f in optimal.columns.drop('method'):
        audit.append(dict(method=r.method,field=f,global_optimum=r[f],R11_numerical_estimator=p.loc[r.method,f],difference=r[f]-p.loc[r.method,f]))
save(pd.DataFrame(audit),O/'optimizer_stopping_audit.tsv')
# Independently check both BY families for every sensitivity_effects block and matched six-test family.
ss=read(S/'evidence/sensitivity_pooled.tsv')
base=read(S/'evidence/within_performance.tsv')
exclusions={
 'Excluding four hemolysis-associated units':['GSE266859','GSE117063','GSE227778','GSE188627'],
 'Excluding the two units with repeated profiles per person':['GSE204951','GSE221088'],
 'Excluding units with unresolved scale or timing or normalized input':['GSE85589','GSE188232','GSE70080','GSE117063'],
 'Excluding GSE266859 (legacy cancer-status sensitivity; eligible melanoma contrast)':['GSE266859'],
 'Excluding the Toray unit':['TORAY_PUBLIC_CANONICAL'], 'Excluding GSE85589':['GSE85589']}
sensitivity_weights=[]
for an,g in ss.groupby('analysis',sort=False):
    source=base.copy()
    if an.startswith('Full 32-unit'):source=read(S/'full32_sensitivity/evidence/within_performance.tsv')
    elif an=='GSE85589 entered as deposited log2 signal':source=read(S/'evidence/within_linear.tsv')
    elif an=='Person-grouped inner resampling':source=read(S/'evidence/within_grouped.tsv')
    elif an.startswith('Foundation model seed tied'):
        for _,r in read(S/'evidence/foundation_seed_sensitivities.tsv').iterrows():
            q=read(S/_public_analysis_directory(r.performance_path));source=pd.concat([source[source.method!=r.method],q],ignore_index=True)
    source=source[~source.cohort.isin(exclusions.get(an,[]))]
    independent_sens=pd.DataFrame([dict(method=m,**own_meta(q.lift.to_numpy(),q.lift_se.to_numpy())) for m,q in source.groupby('method') if m!='ws-rclr-trimmed']).set_index('method')
    for _,r in g.iterrows():
        ww=source[(source.method==r.method)&source.lift.notna()&source.lift_se.gt(0)].copy()
        if len(ww)>=5:
            z=1/(ww.lift_se.to_numpy()**2+r.tau2);z=z/z.sum()*100
            for (_,rr),weight in zip(ww.iterrows(),z):sensitivity_weights.append(dict(analysis=an,method=r.method,cohort=rr.cohort,weight_percent=weight,lift=rr.lift,lift_se=rr.lift_se,tau2=r.tau2))
        for f in ['k','pooled_lift','pooled_se','ci_lo','ci_hi','pi_lo','pi_hi','tau2','I2','p_lift_gt_zero','p_lift_gt_margin_0_05']:
            check('Independent S9 pool '+an,r.method,f,independent_sens.loc[r.method,f],r[f])
    for h in ['p_lift_gt_zero','p_lift_gt_margin_0_05']:
        q=own_by(independent_sens.loc[g.method,h])
        for i,(_,r) in enumerate(g.iterrows()):check('S9 '+an,r.method,h+'_BY',q[i],r[h+'_BY'],1e-12)
for branch in [S,S/'full32_sensitivity']:
    d=read(branch/'evidence/matched_baseline_dataset.tsv'); p=read(branch/'evidence/matched_baseline_pooled.tsv').set_index('method')
    q=pd.DataFrame([dict(method=m,**own_meta(g.estimate.to_numpy(),g.corrected_rcv_se.to_numpy())) for m,g in d.groupby('method')])
    for h in ['p_lift_gt_zero','p_lift_gt_margin_0_05']:q[h+'_BY']=own_by(q[h])
    for _,r in q.iterrows():
        for f in ['k','pooled_lift','pooled_se','ci_lo','ci_hi','pi_lo','pi_hi','p_lift_gt_zero_BY','p_lift_gt_margin_0_05_BY']:
            check('Independent matched frozen-lasso '+branch.name,r.method,f,r[f],p.loc[r.method,f])
save(pd.DataFrame(sensitivity_weights),S/'evidence/sensitivity_cohort_weights.tsv')
save(pd.DataFrame(checks),O/'independent_check.tsv')
bad=pd.DataFrame(checks).query('required and not pass_check')
print(bad.to_string(index=False));print('Checks',len(checks),'required failures',len(bad))
assert bad.empty
