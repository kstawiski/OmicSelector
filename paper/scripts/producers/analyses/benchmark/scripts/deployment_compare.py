from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import pandas as pd,numpy as np,ast,statistics,math,csv,gzip,json,hashlib,os
DEPLOYMENT_INPUT_ORDER={'*__speed.tsv': ['deploy__GSE198834__ai-tabpfn__5/cells/deploy__GSE198834__ai-tabpfn__5__speed.tsv', 'deploy__GSE198834__ai-tabpfn__2/cells/deploy__GSE198834__ai-tabpfn__2__speed.tsv', 'deploy__GSE198834__tab-tabdpt__4/cells/deploy__GSE198834__tab-tabdpt__4__speed.tsv', 'deploy__GSE198834__ai-tabpfn__1/cells/deploy__GSE198834__ai-tabpfn__1__speed.tsv', 'deploy__GSE198834__tab-tabdpt__3/cells/deploy__GSE198834__tab-tabdpt__3__speed.tsv', 'deploy__GSE198834__ai-tabpfn__3/cells/deploy__GSE198834__ai-tabpfn__3__speed.tsv', 'deploy__GSE198834__tab-tabdpt__1/cells/deploy__GSE198834__tab-tabdpt__1__speed.tsv', 'deploy__GSE198834__ai-tabpfn__4/cells/deploy__GSE198834__ai-tabpfn__4__speed.tsv', 'deploy__GSE198834__tab-tabdpt__2/cells/deploy__GSE198834__tab-tabdpt__2__speed.tsv', 'deploy__GSE198834__tab-tabdpt__5/cells/deploy__GSE198834__tab-tabdpt__5__speed.tsv'], '*__meta.tsv': ['deploy__GSE198834__ai-tabpfn__5/cells/deploy__GSE198834__ai-tabpfn__5__meta.tsv', 'deploy__GSE198834__ai-tabpfn__2/cells/deploy__GSE198834__ai-tabpfn__2__meta.tsv', 'deploy__GSE198834__tab-tabdpt__4/cells/deploy__GSE198834__tab-tabdpt__4__meta.tsv', 'deploy__GSE198834__ai-tabpfn__1/cells/deploy__GSE198834__ai-tabpfn__1__meta.tsv', 'deploy__GSE198834__tab-tabdpt__3/cells/deploy__GSE198834__tab-tabdpt__3__meta.tsv', 'deploy__GSE198834__ai-tabpfn__3/cells/deploy__GSE198834__ai-tabpfn__3__meta.tsv', 'deploy__GSE198834__tab-tabdpt__1/cells/deploy__GSE198834__tab-tabdpt__1__meta.tsv', 'deploy__GSE198834__ai-tabpfn__4/cells/deploy__GSE198834__ai-tabpfn__4__meta.tsv', 'deploy__GSE198834__tab-tabdpt__2/cells/deploy__GSE198834__tab-tabdpt__2__meta.tsv', 'deploy__GSE198834__tab-tabdpt__5/cells/deploy__GSE198834__tab-tabdpt__5__meta.tsv']}
B=Path(__file__).resolve().parents[1];M=B.parent;O=B/'out/deployment';O.mkdir(exist_ok=True)
source=M/'scoring_agreement/out/code/aggregate.py';tree=ast.parse(source.read_text());scope=dict(Path=Path,csv=csv,gzip=gzip,statistics=statistics,math=math)
for n in tree.body:
 if isinstance(n,ast.FunctionDef):exec(compile(ast.Module(body=[n],type_ignores=[]),str(source),'exec'),scope)
pred=pd.read_csv(os.getenv('D2_PREDICTIONS',str(B/'out/cells/primary_predictions.tsv.gz')),sep='\t');old=pd.read_csv(M/'scoring_agreement/out/deployed_vs_benchmark.tsv',sep='\t');rows=[];checks=[]
for f in sorted((B/'out/retrieved/deploy2/cells').glob('*/cells/*__scores.tsv')):
 d=pd.read_csv(f,sep='\t');method=d.method.iloc[0];fold=int(d.fold.iloc[0]);q=pred.query('method==@method and seed==101 and fold==@fold');assert len(q)==len(d)>0
 q=d.merge(q,on='sample_id',suffixes=('_deployed','_benchmark'),validate='one_to_one');assert (q.y_deployed==q.y_benchmark).all();assert(q.sample_idx_deployed==q.sample_idx_benchmark).all()
 score=[c for c in d.columns if c=='deployed'];assert len(score)==1,score
 a=q.score_m.to_numpy();b=q[score[0]].to_numpy();y=q.y_deployed.to_numpy();assert np.isfinite(a).all() and np.isfinite(b).all()
 r=old[(old.dataset=='GSE198834')&(old.method==method)&(old.fold==fold)].iloc[0].to_dict();md=pd.read_csv(f.with_name(f.name.replace('__scores.tsv','__meta.tsv')),sep='\t').iloc[0].to_dict();r.update({k:v for k,v in md.items() if k in old.columns});r.update(scope['metrics'](y.tolist(),a.tolist(),b.tolist()));r.update(status='complete',n=len(q),error='');rows.append(r)
 from scipy.stats import rankdata,spearmanr
 def ia(v):return (rankdata(v)[y==1].sum()-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(y==0))
 checks.append(dict(method=method,fold=fold,n=len(q),labels_match=True,sample_ids_match=True,auc_benchmark_error=abs(ia(a)-r['auc_benchmark']),auc_deployed_error=abs(ia(b)-r['auc_deployed']),max_abs_difference_error=abs(np.max(np.abs(a-b))-r['max_absolute_difference'])))
 _public_table(q).to_csv(O/_public_analysis_directory(f'{method}__fold_{fold}__paired_scores.tsv'),sep='\t',index=False)
assert len(rows)==10
new=pd.DataFrame(rows);full=pd.concat([old[old.dataset!='GSE198834'],new[old.columns]],ignore_index=True);_public_table(full).to_csv(O/'deployed_vs_benchmark.tsv',sep='\t',index=False);_public_table(new).to_csv(O/'gse198834_replacement.tsv',sep='\t',index=False);_public_table(pd.DataFrame(checks)).to_csv(O/'independent_checks.tsv',sep='\t',index=False)
assert max(r[k] for r in checks for k in r if k.endswith('_error'))<1e-12
(O/'producer_identity.json').write_text(json.dumps(dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),function='metrics; thresholds; ranks; auc; spearman',n_old_eligible=0,n_new_eligible=10),indent=2))
print('Deployment comparisons',len(rows),flush=True)

# Reuse the recorded aggregate summary block on the complete updated comparison grid.
cov=pd.read_csv(M/'scoring_agreement/out/planned_coverage.tsv',sep='\t')
for _,r in new.iterrows():
 mask=(cov.dataset=='GSE198834')&(cov.method==r.method)&(cov.fold==r.fold)
 for c in ['n','n_train','n_features','feature_count']:cov.loc[mask,c]=r[c]
 cov.loc[mask,['planned_specimens','scored_specimens','completed_specimens']]=int(r['n']);cov.loc[mask,'status']='complete';cov.loc[mask,'stored_eligible']=True;cov.loc[mask,'stored_reason']='';cov.loc[mask,'fold_source']='Corrected-label grouped folds, seed 101'
_public_table(cov).to_csv(O/'planned_coverage.tsv',sep='\t',index=False)
def records(d):return d.astype(object).where(d.notna(),None).to_dict('records')
scope.update(O=O,METHODS=['ai-tabpfn','tab-tabdpt'],bykey={(r['dataset'],r['method'],r['fold']):r for r in records(full[full.status=='complete'])},plan=records(cov),coverage=records(cov))
source_text=source.read_text();block=source_text[source_text.index('summary=[]'):source_text.index('for c in coverage:')];exec(compile(block,str(source),'exec'),scope)
# Replace timing, environment, alignment and score evidence from the newly fitted deployment cells.
for name,pattern in [('speed.tsv','*__speed.tsv'),('environment.tsv','*__meta.tsv'),('fold_alignment.tsv','*__meta.tsv')]:
 oldtab=pd.read_csv(M/'scoring_agreement/out'/_public_analysis_directory(name),sep='\t');parts=[pd.read_csv(f,sep='\t') for f in [B/'out/retrieved/deploy2/cells'/_public_analysis_directory(rel) for rel in DEPLOYMENT_INPUT_ORDER[pattern]]];z=pd.concat(parts,ignore_index=True)
 for c in oldtab.columns:
  if c not in z:z[c]=np.nan
 _public_table(pd.concat([oldtab[oldtab.dataset!='GSE198834'],z[oldtab.columns]],ignore_index=True)).to_csv(O/_public_analysis_directory(name),sep='\t',index=False)
oldfail=pd.read_csv(M/'scoring_agreement/out/failures.tsv',sep='\t');_public_table(oldfail[oldfail.dataset!='GSE198834']).to_csv(O/'failures.tsv',sep='\t',index=False)
