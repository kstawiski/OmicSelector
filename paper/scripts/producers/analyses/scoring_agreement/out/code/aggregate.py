from pathlib import Path

import csv,gzip,json,math,hashlib,statistics,sys,itertools,collections

def read(p):
 if not p.exists():return []
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))

def write(p,rows,fields):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def ranks(xs):
 out=[0.]*len(xs);order=sorted(range(len(xs)),key=xs.__getitem__);j=0
 while j<len(xs):
  k=j+1
  while k<len(xs) and xs[order[k]]==xs[order[j]]:k+=1
  for i in order[j:k]:out[i]=(j+1+k)/2
  j=k
 return out

def auc(y,x):
 n1=sum(y);n0=len(y)-n1
 if n1*n0==0:return None
 return (sum(r for r,a in zip(ranks(x),y) if a)-n1*(n1+1)/2)/(n1*n0)

def spearman(a,b):
 a=ranks(a);b=ranks(b);ma=statistics.mean(a);mb=statistics.mean(b);v=sum((x-ma)**2 for x in a)*sum((x-mb)**2 for x in b)
 return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(v) if v else None

def thresholds(y,x):
 uniq=sorted(set(x));ths=[-math.inf]+[(a+b)/2 for a,b in zip(uniq,uniq[1:])]+[math.inf];n1=sum(y);n0=len(y)-n1
 # Use the exact integer numerator of sensitivity + specificity - 1 to retain ties.
 if n1*n0:
  vals=[sum(int(v>=t) for a,v in zip(y,x) if a)*n0+sum(int(v<t) for a,v in zip(y,x) if not a)*n1 for t in ths]
  youden=ths[vals.index(max(vals))]
 else:youden=None
 return dict(fixed=.5,youden=youden,median=statistics.median(x))

def metrics(y,a,b):
 aa=auc(y,a);ab=auc(y,b);r=dict(max_absolute_difference=max(abs(x-z) for x,z in zip(a,b)),spearman=spearman(a,b),auc_benchmark=aa,auc_deployed=ab,auc_difference=ab-aa if aa is not None and ab is not None else None)
 for k,t in thresholds(y,a).items():
  ch=sum((x>=t)!=(z>=t) for x,z in zip(a,b)) if t is not None else None
  r.update({f'threshold_{k}':t,f'changed_{k}':ch,f'share_changed_{k}':ch/len(y) if ch is not None else None})
 return r

summary=[]
for m in METHODS:
 for group in ['all','<=128','>128']:
  rr=[r for r in bykey.values() if r['method']==m and (group=='all' or r['feature_count']==group)]
  def extreme(key,fn=max,absolute=False):
   v=[r[key] for r in rr if r.get(key) is not None];v=[abs(x) for x in v] if absolute else v
   return fn(v) if v else None
  s=dict(method=m,feature_count=group,datasets=len({r['dataset'] for r in rr}),folds=len(rr),specimens=sum(int(r['n']) for r in rr),maximum_difference=extreme('max_absolute_difference'),minimum_spearman=extreme('spearman',min),minimum_auc_difference=extreme('auc_difference',min),maximum_auc_difference=extreme('auc_difference'),maximum_absolute_auc_difference=extreme('auc_difference',absolute=True))
  for t in ['fixed','youden','median']:
   ch=sum(r['changed_'+t] for r in rr if r['changed_'+t] is not None);n=sum(int(r['n']) for r in rr if r['changed_'+t] is not None);s.update({f'changed_{t}':ch,f'denominator_{t}':n,f'share_changed_{t}':ch/n if n else None})
  planned=[r for r in plan if r['method']==m and (group=='all' or r['feature_count']==group)]
  selected_coverage=[r for r in coverage if r['method']==m and (group=='all' or r['feature_count']==group)]
  s.update(planned_datasets=len({r['dataset'] for r in planned}),planned_folds=len(planned),planned_specimens=sum(int(r['n']) for r in planned),scored_folds=sum(int(r['scored_specimens'])>0 for r in selected_coverage),scored_specimens=sum(int(r['scored_specimens']) for r in selected_coverage),pending_folds=sum(r['status'].startswith('pending') for r in selected_coverage),failed_folds=sum(r['status']=='failed' for r in selected_coverage),unavailable_benchmark_folds=sum(r['status']=='benchmark_unavailable' for r in selected_coverage),undefined_spearman_folds=sum(r['spearman'] is None for r in rr))
  summary.append(s)
write(O/'summary.tsv',summary,list(summary[0]))
for c in coverage:
 pass
