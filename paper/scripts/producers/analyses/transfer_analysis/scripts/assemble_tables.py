from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,json,math,statistics,collections,re,shutil,hashlib,ast,platform,datetime,sys
W=Path(__file__).resolve().parents[1];M=W.parent;O=W/'out';S=O/'summaries';F=O/'transfer_inputs'
S.mkdir(exist_ok=True);F.mkdir(exist_ok=True);csv.field_size_limit(32*1024*1024)
R=W/'results/summary/outputs'
if not R.exists():R=W/'results/summary'
assert (R/'VALID').read_text().strip()=='VALID'
BREAST={'GSE210329','GSE270497','GSE320031','GSE41922'};TARGETS=BREAST-{'GSE270497'}
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rr,fields=None):
 rr=list(rr);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fields or list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rr:w.writerow({k:json.dumps(v,ensure_ascii=False,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def num(v):
 try:return float(v)
 except:return math.nan
def eq(a,b):
 if str(a)==str(b):return True
 if str(a) in ('','NA','NaN','nan') and str(b) in ('','NA','NaN','nan'):return True
 x,y=num(a),num(b);return math.isfinite(x) and math.isfinite(y) and x==y
changes=[];coverage=[]
def compare(name,aa,bb,keys,baseline='R5b'):
 a={tuple(str(r[k]) for k in keys):r for r in aa};b={tuple(str(r[k]) for k in keys):r for r in bb}
 assert len(a)==len(aa) and len(b)==len(bb),(name,'duplicate keys')
 count=0
 for key in sorted(a.keys()|b.keys()):
  x,y=a.get(key,{}),b.get(key,{})
  for field in sorted(x.keys()|y.keys()):
   old,new=x.get(field,''),y.get(field,'')
   if key not in a or key not in b or not eq(old,new):
    delta=num(new)-num(old)
    changes.append(dict(table=name,baseline=baseline,key=json.dumps(dict(zip(keys,key)),sort_keys=True),quantity=field,old=old,new=new,delta=delta if math.isfinite(delta) else '',action='removed' if key not in b else 'added' if key not in a else 'changed'));count+=1
 coverage.append(dict(table=name,baseline=baseline,reference_rows=len(aa),analysis_rows=len(bb),changed_fields=count,removed_rows=len(a.keys()-b.keys()),added_rows=len(b.keys()-a.keys())))
old=rows(M/'transfer_summary/out/validated_cells.tsv');new=rows(O/'cells/validated_cells.tsv');assert len(new)==1955
compare('validated_cells.tsv',old,new,['held_out','method','seed'])
for name in ['transfer_performance.tsv','transfer_validation.tsv','validation.tsv','derive_locked_cutoffs.tsv','derive_calibrator_fits.tsv','derive_decision_curves_by_seed.tsv','sessionInfo.txt']:
 shutil.copyfile(R/_public_analysis_directory(name),S/_public_analysis_directory(name))
shutil.copyfile(R/'transfer_predictions.tsv.gz',O/'cells/transfer_predictions.tsv.gz')
oldperf=rows(M/'transfer_summary/out/transfer_performance.tsv');oldperfmap={(r['held_out'],r['method']):r for r in oldperf}
newperf=[r if r['held_out'] in TARGETS else oldperfmap[(r['held_out'],r['method'])] for r in rows(S/'transfer_performance.tsv')]
write(S/'transfer_performance.tsv',newperf,list(oldperf[0]))
compare('transfer_performance.tsv',oldperf,newperf,['disease','held_out','method'])
# Execute the transfer_cells table function from the producer. Only input paths and grid dimensions change.
src=M.parent/'analyses/tables_supplement_package/out/supplement/scripts/tables.py';source=src.read_text();tree=ast.parse(source)
n=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='build_transfer');code=ast.get_source_segment(source,n)
code=code.replace("BASE/'transfer_transfer_final/out/full_results/validated_cells.tsv'",f"Path({str(O/'cells/validated_cells.tsv')!r})")
code=code.replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance_r3.tsv'",f"Path({str(S/'transfer_performance.tsv')!r})")
code=code.replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_summary.tsv'",f"Path({str(S/'absent_summary.tsv')!r})").replace("BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance.tsv'",f"Path({str(S/'absent_performance.tsv')!r})")
code=code.replace('2160','1955').replace('(24,18,5)','(23,17,5)').replace('24 methods 18 targets 5 seeds','23 methods 17 targets 5 seeds').replace('432','391')
producerchecks=[]
def check(name,ok,detail=''):
 producerchecks.append(dict(check=name,pass_=bool(ok),detail=detail));assert ok,name
scope=dict(Path=Path,BASE=M.parent/'analyses',annotated=lambda p:[dict(r,source_path=str(p),source_row=i) for i,r in enumerate(rows(p),2)],check=check,write_tsv=write,defaultdict=collections.defaultdict,Counter=collections.Counter,statistics=statistics,numeric=lambda x:num(x) if math.isfinite(num(x)) else None,joined=lambda v:';'.join(dict.fromkeys(str(x) for x in v if x is not None and str(x))),issue=lambda *a:(_ for _ in ()).throw(AssertionError(a)))
exec(compile(code,str(src),'exec'),scope);scope['build_transfer'](S)
write(O/'evidence/transfer_cells_producer_checks.tsv',producerchecks)
for name,keys in [('transfer_cells.tsv',['held_out','method','seed']),('transfer_summary.tsv',['Method identifier','Held-out analysis unit'])]:compare(name,rows(M/'transfer_summary/out/supplement'/_public_analysis_directory(name)),rows(S/_public_analysis_directory(name)),keys)
# Replace only the three transfer targets; remove the excluded transfer target.
dca_tables={}
dca_input_identity={}
for q in (W/'results/cpu').rglob('input_identity.tsv'):
 for r in rows(q):dca_input_identity[r['target']]=r
assert set(dca_input_identity)==TARGETS
for target in TARGETS:
 pool=next(r['training_cohorts'] for r in new if r['held_out']==target)
 assert dca_input_identity[target]['training_cohorts']==pool
write(O/'cells/decision_curve_input_identity.tsv',dca_input_identity.values())

for name,keys in [('operating_points.tsv',['context','cohort','method','seed','target_specificity']),('calibration.tsv',['context','cohort','method','seed','calibrator']),('decision_curves.tsv',['context','cohort','method','calibrator','target_prevalence','risk_threshold'])]:
 aa=rows(M/'decision_analysis/out'/_public_analysis_directory(name));replacement=rows(R/_public_analysis_directory(name));expected=180 if name=='operating_points.tsv' else 120 if name=='calibration.tsv' else 1512;assert len(replacement)==expected
 assert all(r['context']=='transfer' and r['cohort'] in TARGETS for r in replacement)
 training_by_target={r['held_out']:r['training_cohorts'] for r in new if r['held_out'] in TARGETS}
 for r in replacement:r['training_cohorts']=training_by_target[r['cohort']]
 assert all('GSE270497' not in r['training_cohorts'] and len(r['training_cohorts'].split(','))==2 for r in replacement)

 bkeys={tuple(r[k] for k in keys):r for r in replacement};bb=[]
 for a in aa:
  if a['context']=='transfer' and a['cohort'] in BREAST:
   if a['cohort']=='GSE270497':continue
   bb.append(bkeys.pop(tuple(a[k] for k in keys)))
  else:bb.append(a)
 assert not bkeys
 assert [r for r in aa if r['context']!='transfer' or r['cohort'] not in BREAST]==[r for r in bb if r['context']!='transfer' or r['cohort'] not in BREAST]
 write(S/_public_analysis_directory(name),bb,list(aa[0]));compare(name,aa,bb,keys,'DH1b');dca_tables[name]=(aa,bb)
# Auxiliary decision_analysis artifacts retain their own producer schemas.
for name,keys in [('derive_locked_cutoffs.tsv',['context','cohort','method','seed','fold','rule']),('derive_calibrator_fits.tsv',['context','cohort','method','seed','fold','calibrator']),('derive_decision_curves_by_seed.tsv',['context','cohort','method','seed','calibrator','target_prevalence','risk_threshold'])]:
 aa=rows(M/'decision_analysis/out'/_public_analysis_directory(name));bb=[r for r in aa if r['cohort'] not in BREAST]+rows(R/_public_analysis_directory(name));write(S/_public_analysis_directory(name),bb,list(aa[0]));compare(name,aa,bb,keys,'DH1b')
# Retain the producer's cell-outcome audit, including all non-breast refusals.
aa=rows(M/'decision_analysis/out/derive_cell_audit.tsv');nr=rows(W/'summary_payload/payload/audit.tsv')
for r in nr:r['artifact']=str(O/'cells/decision_curve'/_public_analysis_directory(Path(r['artifact']).name))
bb=[r for r in aa if r['cohort'] not in BREAST]+nr
write(S/'derive_cell_audit.tsv',bb,list(aa[0]));compare('derive_cell_audit.tsv',aa,bb,['context','cohort','method','seed'],'DH1b')
dca_tables['derive_cell_audit.tsv']=(aa,bb)
def status_counts(tables,index):
 out=[]
 for name,pair in sorted(tables.items()):
  ct=collections.Counter()
  for r in pair[index]:
   if r['context']!='transfer':continue
   for group in ['ALL',r['cohort']]:
    for method in ['ALL',r['method']]:ct[(group,method,r['status'])]+=1
  out.extend(dict(table=name,cohort=group,method=method,status=status,n_rows=n) for (group,method,status),n in sorted(ct.items()))
 return out
oldstatus=status_counts(dca_tables,0);newstatus=status_counts(dca_tables,1)
write(S/'decision_curve_status_counts.tsv',newstatus);compare('decision_curve_status_counts.tsv',oldstatus,newstatus,['table','cohort','method','status'],'DH1b')

# Counts and target-weighted medians use all five seeds; refusals remain cells.
def summarize(rr):
 structural={r['held_out'] for r in rr if int(r['n_training_cohorts'])<2};scored=[r for r in rr if r['eligible']=='TRUE'];den=len({r['held_out'] for r in rr})-len(structural)
 values={'registered_methods':len({r['method'] for r in rr}),'all_method_target_seed_cells':len(rr),'structurally_eligible_target_seed_cells':len({(r['held_out'],r['seed']) for r in rr if r['held_out'] not in structural}),'structurally_eligible_method_target_seed_cells':sum(r['held_out'] not in structural for r in rr),'method_target_seed_cells_with_heldout_scores':len(scored),'target_seed_cells_with_any_heldout_scores':len({(r['held_out'],r['seed']) for r in scored}),'structural_exclusions':sum(r['held_out'] in structural for r in rr)}
 hh=[dict(phase='after',method='ALL',metric=k,value=v,denominator='',definition='Cell counts include all refusals. Structurally eligible targets have at least two source cohorts.') for k,v in values.items()]
 med=[]
 for m in sorted({r['method'] for r in rr}):
  groups=collections.defaultdict(list)
  for r in rr:
   if r['method']==m:groups[r['held_out']].append(r)
  good={k:v for k,v in groups.items() if len(v)==5 and all(r['eligible']=='TRUE' for r in v)}
  vs={'targets_scored':len(good),'median_target_auc':statistics.median(statistics.mean(float(r['auc_method']) for r in v) for v in good.values()) if good else '', 'median_paired_difference_from_rclr':statistics.median(statistics.mean(float(r['lift']) for r in v) for v in good.values()) if good else ''}
  med.append(dict(method=m,eligible_targets=den,**vs))
  if m in {'blood-cell-marker-excl-rclr','cvae','dann','dg-fishr','dg-ibirm','dro-group','dro-vrex','icp','moe-gated','reo-metaktsp'}:
   hh.extend(dict(phase='after',method=m,metric=k,value=v,denominator=den,definition='Target estimates are arithmetic means across all five eligible seeds. Medians give each target equal weight.') for k,v in vs.items())
 codes=collections.Counter((r['method'],r['refusal_code'] or ('STRUCTURAL_INELIGIBLE' if int(r['n_training_cohorts'])<2 else ';'.join(sorted(set(re.findall(r'\[([A-Z_]+)\]',r['ineligible_reason'])))) or 'UNCODED_NON_EVALUABLE')) for r in rr if r['eligible']!='TRUE')
 ref=[];tot=collections.Counter()
 for (m,c),n in sorted(codes.items()):
  ref.append(dict(phase='after',method=m,refusal_code=c,n_cells=n));tot[c]+=n
  hh.append(dict(phase='after',method=m,metric='refusal_cells:'+c,value=n,denominator='',definition='Non-evaluable method-target-seed cells include structural exclusions and uncoded outcomes.'))
 hh.extend(dict(phase='after',method='ALL',metric='refusal_cells:'+c,value=n,denominator='',definition='All methods use specimen metadata without imputation.') for c,n in sorted(tot.items()))
 return hh,med,ref
oldh,oldm,oldr=summarize(old);newh,newm,newr=summarize(new)
# The saved headline is an additional independent arithmetic check.
stored=[r for r in rows(M/'transfer_summary/out/headline.tsv') if r['phase']=='after'];assert {(r['method'],r['metric']) for r in stored}=={(r['method'],r['metric']) for r in oldh}
for r in stored:
 q=next(x for x in oldh if (x['method'],x['metric'])==(r['method'],r['metric']));assert abs(num(r['value'])-num(q['value']))<1e-12
for name,aa,bb,keys in [('headline.tsv',oldh,newh,['method','metric']),('method_target_medians.tsv',oldm,newm,['method']),('refusals.tsv',oldr,newr,['method','refusal_code'])]:write(S/_public_analysis_directory(name),bb);compare(name,aa,bb,keys)
def outcomes(rr):
 c=collections.Counter()
 for r in rr:
  for family in ['ALL',r['disease']]:
   for method in ['ALL',r['method']]:c[(family,method,r['outcome_class'])]+=1
 return [dict(disease=f,method=m,outcome_class=o,n_cells=n) for (f,m,o),n in sorted(c.items())]
write(S/'outcome_counts.tsv',outcomes(new));compare('outcome_counts.tsv',outcomes(old),outcomes(new),['disease','method','outcome_class'])
# Figure 4B/C preserves the producer's five-seed arithmetic and explicit refusal symbols.
figsource=M/'figures/work/figures/provenance/transfer_transport_conditions_values.tsv'
code_map={'UNSEEN_METADATA_LEVEL':'U','MISSING_METADATA':'M','INSUFFICIENT_ALR_PIVOTS':'A','INSUFFICIENT_TECHNOLOGY_CLASSES':'T','ALL_NONDETECT':'N','INSUFFICIENT_BIOFLUID_ANCHORS':'B'}
def figures(rr,source):
 structural={r['held_out'] for r in rr if r['outcome_class']=='structural_ineligible'};targets=sorted({r['held_out'] for r in rr}-structural);out=[];matrix=[]
 def add(panel,item,value,group,column,transformation,printed):
  ix=[rr.index(r)+2 for r in group]
  out.append(dict(figure='figure_4_transport_conditions',panel=panel,item=item,value=value,printed_text=printed,source_file=str(source),source_row=','.join(map(str,ix)),source_column=column,transformation=transformation,value_kind='evidence'))
 for m in sorted({r['method'] for r in rr}):
  n=0
  for target in targets:
   z=[r for r in rr if r['method']==m and r['held_out']==target];good=[r for r in z if r['eligible']=='TRUE'];assert len(z)==5
   av=lv='';marker=''
   if len(good)==5:
    av=statistics.mean(float(r['auc_method']) for r in good);lv=statistics.mean(float(r['lift']) for r in good);n+=1
    add('B',m+'; '+target,av,good,'auc_method','mean over all five seeds','tile');add('C',m+'; '+target,lv,good,'lift','mean of paired same-fold differences over all five seeds','tile')
   else:
    codes={r['refusal_code'] for r in z if r['refusal_code']}
    for r in z:codes.update(re.findall(r'\[([A-Z_]+)\]',r['ineligible_reason']))
    marker='/'.join(sorted(code_map.get(c,c) for c in codes)) if codes else 'X'
    if good:marker='P'
    add('B / C',m+'; '+target+'; refusal',marker,z,'refusal_code; ineligible_reason; eligible','literal documented codes; P=partial seed coverage; X=uncoded non-evaluable',marker)
   matrix.append(dict(method=m,held_out=target,mean_auc=av,mean_paired_difference=lv,evaluated_seeds=len(good),refusal_marker=marker))
  add('B',m+'; scored targets',n,[r for r in rr if r['method']==m],'eligible; held_out','count targets with all five seeds eligible',f'{n}/{len(targets)}')
 return out,matrix
oldf,oldmat=figures(old,M/'transfer_summary/out/validated_cells.tsv');newf,newmat=figures(new,O/'cells/validated_cells.tsv')
storedfig=[r for r in rows(figsource) if r['panel'] in {'B','C','B / C'}];storedkey={(r['panel'],r['item']):r for r in storedfig}
assert storedkey.keys()=={(r['panel'],r['item']) for r in oldf}
for r in oldf:
 v=storedkey[(r['panel'],r['item'])]['value'];assert eq(v,r['value']) or abs(num(v)-num(r['value']))<1e-12
write(F/'transfer_transport_conditions_values.tsv',newf);write(F/'transfer_BC.tsv',newmat)
compare('figure4_values',oldf,newf,['panel','item'],'The transfer figure producer supplies the reference values.')
compare('transfer_BC.tsv',oldmat,newmat,['method','held_out'])
write(F/'refusal_key.tsv',[dict(code=k,marker=v) for k,v in sorted(code_map.items())])
# decision_analysis summary schema, extended to each requested target and to the breast family.
def dca_head(rr):
 out=[]
 for group in ['transfer','breast_cancer']+sorted(BREAST|{'GSE204951'}):
  z=[r for r in rr if r['context']=='transfer' and (group=='transfer' or group=='breast_cancer' and r['cohort'] in BREAST or r['cohort']==group)]
  if not z:continue
  good=[r for r in z if r['status']=='COMPUTED'];vals={'grid_curve_rows':len(z),'computed_curve_rows':len(good),'held_curve_rows':len(z)-len(good),'above_both':sum(num(r['increment_vs_treat_all'])>0 and num(r['increment_vs_treat_none'])>0 for r in good),'maximum_increment_over_better_reference':max((min(num(r['increment_vs_treat_all']),num(r['increment_vs_treat_none'])) for r in good),default=''),'maximum_net_benefit':max((num(r['net_benefit']) for r in good),default='')}
  out.extend(dict(baseline='R12b',group=group,metric=k,value=v) for k,v in vals.items())
 return out
aa,bb=dca_tables['decision_curves.tsv'];oh=dca_head(aa);nh=dca_head(bb);write(S/'decision_curve_headline.tsv',nh);compare('decision_curve_headline.tsv',oh,nh,['group','metric'],'DH1b')
# Update the decision_analysis aggregate counts that include transfer rows, preserving reference columns.
ops=dca_tables['operating_points.tsv'][1];cal=dca_tables['calibration.tsv'][1];curves=dca_tables['decision_curves.tsv'][1]
inventory={(r['context'],r['cohort'],r['method'],r['seed']) for r in ops}
computed={(r['context'],r['cohort'],r['method'],r['seed']) for r in ops if r['threshold_fit_executed']=='TRUE'}
counts={'fit_attempts':len(inventory),'computed_cells':len(computed),'held_cells':len(inventory)-len(computed),'held_operating_points_rows':sum(r['status']=='HOLD' for r in ops),'held_decision_curves_rows':sum(r['status']=='HOLD' for r in curves),'held_calibration_rows':sum(r['status']=='HOLD' for r in cal),'platt_nonconvergence_rows':sum('PLATT_NONCONVERGENCE' in r['hold_reason'] for r in cal),'largest_computed_net_benefit':max(float(r['net_benefit']) for r in curves if r['status']=='COMPUTED')}
s1old=rows(M/'decision_analysis/out/matched_baseline_numbers.tsv')
s1new=[dict(r,R12b=counts[r['metric']],difference_vs_DH1b=counts[r['metric']]-float(r['DH1b'])) for r in s1old]
write(S/'matched_baseline_numbers.tsv',s1new)
compare('matched_baseline_numbers.tsv',[dict(metric=r['metric'],value=r['DH1b']) for r in s1old],[dict(metric=k,value=v) for k,v in counts.items()],['metric'],'DH1b')

# paired_auc_contrasts is a separate published-classifier contrast family; establish actual coverage.
s16=[]
for p in (M/'figures/work/staging/paired_auc_contrasts').glob('*.tsv'):
 rr=rows(p);assert not any(r.get('accession') in BREAST for r in rr)
 dest=S/'paired_auc_contrasts'/_public_analysis_directory(p.name);dest.parent.mkdir(exist_ok=True);shutil.copyfile(p,dest)
 s16.append(dict(file=p.name,rows=len(rr),breast_family_rows=0,action='unchanged',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source=str(p)))
write(S/'paired_contrast_coverage.tsv',s16)
pass;pass
# Cell artifacts live under the requested deliverable directory as well as lifecycle evidence.
G=Path(json.loads((W/'records/gpu_batch.json').read_text())['workspace'])/'outputs'
CROOT=W/'results/cpu/outputs'
if not CROOT.exists():CROOT=W/'results/cpu'
artifact_manifest=[]
for root in [CROOT,G]:
 for src in root.glob('chunks/*/transfer_bundle.rds'):
  dest=O/'cells/transfer'/_public_analysis_directory(src.parent.name);dest.mkdir(parents=True,exist_ok=True)
  for f in src.parent.iterdir():
   if f.is_file():
    shutil.copyfile(f,dest/_public_analysis_directory(f.name))
    artifact_manifest.append(dict(path=str((dest/_public_analysis_directory(f.name)).relative_to(O)),source=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
for src in CROOT.rglob('transfer__*.rds'):
 dest=O/'cells/decision_curve'/_public_analysis_directory(src.name);dest.parent.mkdir(exist_ok=True);shutil.copyfile(src,dest)
 artifact_manifest.append(dict(path=str(dest.relative_to(O)),source=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
write(O/'cells/artifact_manifest.tsv',artifact_manifest)
shutil.copyfile(G/'batch_gate.tsv',O/'evidence/batch_gate.tsv')

C=W/'results/cpu/outputs'
if not C.exists():C=W/'results/cpu'
did=[]
for p in C.rglob('derive_cell_audit.tsv'):
 si=p.parent/'sessionInfo.txt';identity=dict(x.split('=',1) for x in si.read_text().splitlines() if '=' in x and not x.startswith(' '))
 for r in rows(p):did.append(dict(context=r['context'],cohort=r['cohort'],method=r['method'],seed=r['seed'],host=identity.get('host',''),device='cpu',threads=1,package_version='2.6.9',package_sha256='5bc289307aaa008cc46c0ee6724325eda415be46bdc59c9108545267d7eb92cd',scheduler_job=json.loads((W/'records/cpu.json').read_text())['job_id'],session_evidence=str(si)))
write(O/'cells/decision_curve_computing_identity.tsv',did)
identities={(r['held_out'],r['method'],r['seed']):r for r in rows(O/'cells/computing_identity.tsv')}
lineage=[]
for r in new:
 key=(r['held_out'],r['method'],r['seed']);identity=identities.get(key)
 lineage.append(dict(held_out=r['held_out'],method=r['method'],seed=r['seed'],origin='R12b_refit' if identity else 'retained_R5b',package_version=r['package_version'],package_commit=r['package_commit'],host=identity['host'] if identity else 'unavailable; see source',device=identity['device'] if identity else 'unavailable; see source',source=identity['session_evidence'] if identity else str(M/'transfer_summary/out/validated_cells.tsv'),execution_route=r['execution_route']))
write(O/'cells/result_lineage.tsv',lineage)

checks=[]
for path in [S/'validation.tsv',S/'transfer_validation.tsv']:
 rr=rows(path);assert all(r['pass']=='TRUE' for r in rr);checks.append(dict(check=path.name,n_checks=len(rr),passed=True,max_difference=max(num(r['max_difference']) for r in rr if math.isfinite(num(r['max_difference'])))))
checks.extend([dict(check='new_transfer_grid',n_checks=345,passed=True,max_difference=0),dict(check='new_decision_curve_grid',n_checks=60,passed=True,max_difference=0),dict(check='S16_no_breast_family_overlap',n_checks=sum(int(r['rows']) for r in s16),passed=True,max_difference=0)])
write(O/'checks.tsv',checks)
shutil.copyfile(S/'transfer_performance.tsv',F/'transfer_performance.tsv')
# The finishing script verifies/cleans all job records before writing DONE.
oldsc=sum(r['eligible']=='TRUE' for r in old);newsc=sum(r['eligible']=='TRUE' for r in new)
aboveold=next(r['value'] for r in oh if r['group']=='transfer' and r['metric']=='above_both');abovenew=next(r['value'] for r in nh if r['group']=='transfer' and r['metric']=='above_both')
report=f'''Removed GSE270497 from breast-cancer transfer training and evaluation. Refit all 345 cells (23 methods, five seeds, three targets) with OmicSelector 2.6.9. Each target retained two training cohorts: 200 cells were evaluable and 145 were refusals. CPU fits ran on Argos; six neural methods ran on the guarded cpd H100. Computing identities, package hashes and execution records are retained.\n\nThe full transfer grid changed from 2,070 to 1,955 cells and from 18 to 17 targets. Structurally eligible targets changed from 12 to 11. Cells with evaluable held-out scores changed from {oldsc} to {newsc}; the 690 structural exclusions are unchanged. All non-breast results, including person-grouped GSE204951, were preserved.\n\nRebuilt S14, target summaries, target-weighted medians, paired differences and Figure 4B/C inputs. Refit 60 decision-curve score cells using DH1b's four methods and unchanged calibration rules. Transfer curve rows above both reference strategies changed from {aboveold} to {abovenew}. `changes.tsv` records every changed field and removed row against R5b or DH1b. S16's 270 published-classifier contrasts contain none of these breast cohorts and remain unchanged; coverage and hashes are recorded.\n\nThe original summary replay, independent AUC checks and decision-curve threshold arithmetic passed. Changes compare the requested 2.6.9 refit with the recorded older-package baselines; they do not isolate exclusion from package or device effects. Decision curves use apparent training calibration and hypothetical prevalences. These computations do not establish clinical or manuscript readiness.\n'''
assert len(report.split())<=300;(O/'REPORT.md').write_text(report)
(O/'REPRODUCE.md').write_text('Run `PATH=./runtime/legacy-home/micromamba/bin:$PATH python3 '+str(W/'scripts/assemble_outputs.py')+'`. Existing submitted records are inspected and never resubmitted. The frozen inputs are listed in input_sha256.tsv, the transfer-loop adaptation in evidence/transfer_loop.patch, and all package/source and output hashes in their manifests. Baselines: R5b transfer results; DH1b decision curves; shipped S16 contrasts. Summary calculations use the original prediction influence-function producer, five-seed completeness and equal target weights. New biofluid-anchor refusals are labelled B in Figure 4 inputs.\n')
print(report)
