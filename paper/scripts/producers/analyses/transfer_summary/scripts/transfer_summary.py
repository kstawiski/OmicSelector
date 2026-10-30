from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,subprocess,os,math,shutil,statistics,collections,sys,importlib.util,json,datetime,re
W=Path(__file__).resolve().parents[1];REV=W.parent;D=REV.parent;OLD=D/'analyses';O=W/'out'
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rr,fields=None):
 rr=list(rr);fields=fields or list(dict.fromkeys(k for r in rr for k in r))
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def same(a,b,tol=1e-9):
 if a==b:return True
 if a in ('','NA','NaN') and b in ('','NA','NaN'):return True
 try:return abs(float(a)-float(b))<=tol
 except:return False
def milestone(s):
 with (O/'milestones.log').open('a') as f:f.write(datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s+'\n')
for mode,n,src in [('original',432,OLD/'transfer_transfer_final/out/full_results/merged'),('merged',414,O)]:
 td=W/'tmp'/(_public_analysis_directory('summary_'+mode));td.mkdir(exist_ok=True)
 for dest,origin in [('transfer_cells.tsv',src/(_public_analysis_directory('transfer_cells_r3.tsv' if mode=='original' else 'validated_cells.tsv'))),('transfer_predictions.tsv.gz',src/(_public_analysis_directory('transfer_predictions_r3.tsv.gz' if mode=='original' else 'transfer_predictions.tsv.gz')))]:
  p=td/_public_analysis_directory(dest)
  if not p.exists():p.symlink_to(origin)
 env=os.environ.copy();env.update(R5B_SUMMARY_DIR=str(td),R5B_EXPECTED_SUMMARIES=str(n))
 with (O/'evidence'/(_public_analysis_directory('summary_'+mode+'.log'))).open('w') as log:
  subprocess.run(['bash',str(W/'scripts/run_R.sh'),str(W/'scripts/resummarize.R')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 result=td/'transfer_performance_r3.tsv'
 if mode=='original':
  old=rows(src/'transfer_performance_r3.tsv');new=rows(result);keys=['disease','held_out','method'];aa={tuple(r[k] for k in keys):r for r in old};bb={tuple(r[k] for k in keys):r for r in new};assert aa.keys()==bb.keys()
  checks=[]
  for k in aa:
   for field in aa[k]:checks.append(dict(disease=k[0],held_out=k[1],method=k[2],quantity=field,stored=aa[k][field],replayed=bb[k][field],passed=same(aa[k][field],bb[k][field])))
  write(O/'transfer_summary_gate.tsv',checks);assert all(r['passed'] for r in checks)
  shutil.copyfile(result,O/'evidence/transfer_performance_original_replayed.tsv')
  milestone(f'Original transfer summary gate passed for {len(checks)} fields across 432 method-target summaries.')
 else:
  schema=list(rows(OLD/'transfer_transfer_final/out/full_results/merged/transfer_performance_r3.tsv')[0]);new=rows(result);assert list(new[0])==schema
  shutil.copyfile(result,O/'transfer_performance.tsv');milestone('Rebuilt 414 transfer performance summaries with the original producer.')
# Execute the transfer_cells producer with only input paths and specified grid dimensions changed.
src=OLD/'tables_supplement_package/out/supplement/scripts'
sys.path.insert(0,str(src));import common
common.READ_ROOTS.append(REV)
t=(src/'tables.py').read_text().replace("BASE/'transfer_transfer_final/out/full_results/validated_cells.tsv'",f"Path({str(O/'validated_cells.tsv')!r})")
t=t.replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance_r3.tsv'",f"Path({str(O/'transfer_performance.tsv')!r})")
t=t.replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_summary.tsv'",f"Path({str(O/'absent_summary.tsv')!r})").replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance.tsv'",f"Path({str(O/'absent_performance.tsv')!r})")
t=t.replace('2160','2070').replace('(24,18,5)','(23,18,5)').replace('24 methods 18 targets 5 seeds','23 methods 18 targets 5 seeds').replace('432 method-target','414 method-target').replace('len(out)==432','len(out)==414')
(W/'scripts/tables_adapted.py').write_text(t)
scope={'__name__':'r5b_supplement', '_public_analysis_directory': _public_analysis_directory, '_public_table': _public_table};exec(compile(t,str(W/'scripts/tables_adapted.py'),'exec'),scope)
dest=O/'supplement';dest.mkdir(exist_ok=True);scope['build_transfer'](dest)
write(O/'supplement/producer_checks.tsv',common.CHECKS)
milestone('Original S14 and S14_summary producer completed on the merged 23-method grid.')
# The nine invariant-learning routes plus blood-cell-marker-excl-rclr need no specimen metadata.
methods=['blood-cell-marker-excl-rclr','cvae','dann','dg-fishr','dg-ibirm','dro-group','dro-vrex','icp','moe-gated','reo-metaktsp']
before=rows(OLD/'transfer_transfer_final/out/full_results/validated_cells.tsv');after=rows(O/'validated_cells.tsv');hh=[];ref=[]
for phase,rr in [('before_archival_24',before),('before', [r for r in before if r['method']!='sel-stablemate']),('after',after)]:
 structural={r['held_out'] for r in rr if int(r['n_training_cohorts'])<2};eligible={(r['held_out'],r['seed']) for r in rr if r['held_out'] not in structural};scored=[r for r in rr if r['eligible']=='TRUE']
 for metric,value in [('registered_methods',len({r['method'] for r in rr})),('all_method_target_seed_cells',len(rr)),('structurally_eligible_target_seed_cells',len(eligible)),('structurally_eligible_method_target_seed_cells',sum(r['held_out'] not in structural for r in rr)),('method_target_seed_cells_with_heldout_scores',len(scored)),('target_seed_cells_with_any_heldout_scores',len({(r['held_out'],r['seed']) for r in scored})),('structural_exclusions',sum(r['held_out'] in structural for r in rr))]:hh.append(dict(phase=phase,method='ALL',metric=metric,value=value,denominator='',definition='Cell counts include all refusals. Structurally eligible targets have at least two source cohorts.'))
 for m in methods:
  groups=collections.defaultdict(list)
  for r in rr:
   if r['method']==m:groups[r['held_out']].append(r)
  ok={k:v for k,v in groups.items() if len(v)==5 and all(r['eligible']=='TRUE' for r in v)}
  for metric,value in [('targets_scored',len(ok)),('median_target_auc',statistics.median(statistics.mean(float(r['auc_method']) for r in v) for v in ok.values()) if ok else ''),('median_paired_difference_from_rclr',statistics.median(statistics.mean(float(r['lift']) for r in v) for v in ok.values()) if ok else '')]:hh.append(dict(phase=phase,method=m,metric=metric,value=value,denominator=len({r['held_out'] for r in rr})-len(structural),definition='Target estimates are arithmetic means across all five eligible seeds. Medians give each target equal weight.'))
 codes=collections.Counter((r['method'],r['refusal_code'] or ('STRUCTURAL_INELIGIBLE' if r['held_out'] in structural else ';'.join(sorted(set(re.findall(r'\[([A-Z_]+)\]',r['ineligible_reason'])))) or 'UNCODED_NON_EVALUABLE')) for r in rr if r['eligible']!='TRUE')
 for (m,code),n in sorted(codes.items()):ref.append(dict(phase=phase,method=m,refusal_code=code,n_cells=n))
write(O/'headline.tsv',hh);write(O/'refusals.tsv',ref)
aa={(r['held_out'],r['method'],r['seed']):r for r in before};cc=[]
for r in after:
 a=aa[(r['held_out'],r['method'],r['seed'])]
 for k in ['eligible','refusal_code','outcome_class','auc_method','auc_baseline','lift','lift_se','n_common_features','ineligible_reason']:
  x,y=a[k],r[k]
  try:delta=float(y)-float(x)
  except:delta=''
  cc.append(dict(analysis='transfer',unit=r['held_out'],method=r['method'],seed=r['seed'],quantity=k,before=x,after=y,delta=delta,changed=not same(x,y,0),action='replaced' if r['held_out'] in {z['target'] for z in rows(REV/'transfer/out/cells_to_compute.tsv')} else 'retained'))
pass
milestone('Wrote before/after headline, refusal and cell change tables.')
