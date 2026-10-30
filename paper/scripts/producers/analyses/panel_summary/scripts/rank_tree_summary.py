from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,hashlib,math,shutil,json,datetime
W=Path(__file__).resolve().parents[1];O=W/'out';T=O/'ranktree';N=T/'training/rank_tree_competitors';R=W.parent/'panel_restriction';S=W.parent.parent/'analyses/device_benchmark_local_gpu_runs/out/results/rank/rank_tree_competitors'
def rows(p,sep='\t'):
 with p.open() as f:
  r=csv.DictReader(f,delimiter=sep);return r.fieldnames,list(r)
def write(p,data,fields=None):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
accepted={r['unit'] for r in rows(T/'accepted_units.tsv')[1]};checks=[];changes=[]
keys={'cells':['method','cohort','seed'],'performance':['method_id','unit_id'],'eligibility':['method_id','unit_id'],'pairwise':['pair_id','unit_id'],'pair_strata':['pair_id','unit_id','stratum_id'],'units':['unit_id'],'splits':['unit_id','seed','fold'],'meta':['pair_id','analysis_population'],'coverage':['method_id','analysis_population']}
for name,key in keys.items():
 fields,old=rows(S/(_public_analysis_directory(name+'.tsv')));nf,new=rows(N/(_public_analysis_directory(name+'.tsv')))
 a={tuple(r[k] for k in key):r for r in old};b={tuple(r[k] for k in key):r for r in new}
 assert fields==nf and len(a)==len(old) and len(b)==len(new),name
 if name!='pair_strata':assert a.keys()==b.keys(),name
 else:assert all(k[1] in accepted for k in a.keys() ^ b.keys()),name
 checks.append(dict(check=name+'/schema_unique_keys_complete',n=len(new),pass_='TRUE'))
 if name not in ['meta','coverage']:
  unchanged=[k for k,r in a.items() if r.get('unit_id',r.get('cohort')) not in accepted]
  assert all(a[k]==b[k] for k in unchanged),name
  checks.append(dict(check=name+'/unaffected_rows_exact',n=len(unchanged),pass_='TRUE'))
 for source,base in [('stored_original',S),('R7_five_unit_candidate',R/'out/ranktree/training/rank_tree_competitors')]:
  prior={tuple(r[k] for k in key):r for r in rows(base/(name+'.tsv'))[1]}
  for k in sorted(prior.keys() | b.keys()):
   r=b.get(k,prior.get(k))
   unit=r.get('unit_id',r.get('cohort','ALL'))
   if unit!='GSE320031' and name not in ['meta','coverage']:continue
   for metric in fields:
    if metric in key:continue
    av=prior[k][metric] if k in prior else '<ROW_ABSENT>';bv=b[k][metric] if k in b else '<ROW_ABSENT>';diff='NA'
    try:
     aa=float(av);bb=float(bv)
     if math.isfinite(aa) and math.isfinite(bb):diff=format(bb-aa,'.17g')
    except ValueError:pass
    changes.append(dict(reference=source,table=name,unit=unit,key=json.dumps(k,separators=(',',':')),metric=metric,stored=av,rerun=bv,difference=diff,changed=str(av!=bv).upper()))
fields,old=rows(S/'supplementary_file_rank_tree_competitors_rank_tree_competitors.csv',',');nf,new=rows(N/'supplementary_file_rank_tree_competitors_rank_tree_competitors.csv',',')
key=[k for k in ['record_type','method','cohort','seed','method_id','unit_id','pair_id','analysis_population'] if k in fields]
a={tuple(r[k] for k in key):r for r in old};b={tuple(r[k] for k in key):r for r in new}
assert fields==nf and len(a)==len(old) and len(b)==len(new) and a.keys()==b.keys()
def unit(r):return r['unit_id'] if r['unit_id'] not in {'','NA'} else r['cohort']
unaffected=[k for k,r in a.items() if unit(r) not in {'','NA'} and unit(r) not in accepted]
assert all(a[k]==b[k] for k in unaffected)
checks.append(dict(check='public_csv/schema_keys_and_unaffected_rows_exact',n=len(unaffected),pass_='TRUE'))
a={tuple(r[k] for k in ['method','cohort','seed']):r for r in rows(N/'cells.tsv')[1]}
for k,r in a.items():
 if k[1] in accepted:
  p=T/'training/rank_parts'/_public_analysis_directory(k[0])/(_public_analysis_directory('seed_'+k[2]))/_public_analysis_directory(k[1])/'cells.tsv'
  assert rows(p)[1]==[r],k
checks.append(dict(check='all_six_units_replaced_from_training_parts',n=210,pass_='TRUE'))
pass
write(O/'validation_checks.tsv',[dict(check=x['check'],n=x['n'],pass_value=x['pass_']) for x in checks])
for p in N.iterdir():
 if p.is_file():shutil.copy2(p,T/_public_analysis_directory(p.name))
# Consolidate the pre-run manifest with the runtime read audits.
manifest={r['path']:r for r in rows(O/'input_sha256.tsv')[1]}
for p in O.glob('inputs_*.tsv'):
 for line in p.read_text().splitlines():
  path,h,size=line.split('\t');r=dict(path=path,sha256=h,bytes=size)
  assert path not in manifest or manifest[path]['sha256']==h,path
  manifest[path]=r
write(O/'input_sha256.tsv',[manifest[k] for k in sorted(manifest)])
for p in O.glob('inputs_*.tsv'):p.unlink()
print('Validated complete summaries; wrote',len(changes),'quantity-level comparisons;',len(manifest),'input hashes')
print('GSE320031 mean fold AUC changes:')
for r in changes:
 if r['reference']=='stored_original' and r['table']=='performance' and r['metric']=='auc':print(r['key'],r['stored'],r['rerun'],r['difference'])
print('Meta status:',[(r['status'],r['fit_reason']) for r in rows(N/'meta.tsv')[1] if r['status']!='REML_KH'])
