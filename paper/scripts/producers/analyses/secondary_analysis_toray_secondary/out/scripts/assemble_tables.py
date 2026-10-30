from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,json,math,hashlib,gzip,difflib,datetime
B=Path('./analyses'); O=B/'secondary_analysis_toray_secondary/out'; ROOT=B.parent.parent
checks=[]
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
 rows=list(rows);fields=fields or (list(rows[0]) if rows else [])
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def number(v):
 try:return float(v)
 except:return math.nan
def same(a,b,tol=1e-10):
 if a==b:return True
 if a in ('','NA','NaN') and b in ('','NA','NaN'):return True
 x,y=number(a),number(b)
 return math.isfinite(x) and math.isfinite(y) and abs(x-y)<=tol

def compare(name,old,new,keys,table='',target=None):
 aa={tuple('' if r.get(k,'') in ('NA','NaN','') else r.get(k,'') for k in keys):r for r in old};bb={tuple('' if r.get(k,'') in ('NA','NaN','') else r.get(k,'') for k in keys):r for r in new}
 assert len(aa)==len(old),(name,table,'old duplicate keys');assert len(bb)==len(new),(name,table,'new duplicate keys')
 out=[]
 for k in sorted(aa.keys()|bb.keys()):
  a=aa.get(k,{});b=bb.get(k,{})
  if target and not target(a or b):continue
  for v in (a.keys()|b.keys())-set(keys):
   av=a.get(v,'NA');bv=b.get(v,'NA');x,y=number(av),number(bv)
   out.append(dict(analysis=name,table=table,**dict(zip(keys,k)),quantity=v,old=av,new=bv,delta=y-x if math.isfinite(x) and math.isfinite(y) else 'NA',changed=not same(av,bv,0)))
 return sorted(out,key=lambda r:tuple(str(r.get(k,'')) for k in ['table']+keys+['quantity']))

def check(name,passed,detail):checks.append(dict(check=name,pass_=passed,evidence=detail))
def validate_against_saved(name,computed,saved,keys,metrics=None,tol=1e-10):
 aa={tuple('' if r.get(k,'') in ('NA','NaN','') else r.get(k,'') for k in keys):r for r in saved};ds=[];bad=[];n=0
 for r in computed:
  key=tuple('' if r.get(k,'') in ('NA','NaN','') else r.get(k,'') for k in keys)
  if key not in aa:bad.append(str(key)+':missing');continue
  a=aa[key]
  for v in metrics or (r.keys()&a.keys())-set(keys):
   if v not in a or v not in r:continue
   n+=1
   if not same(r[v],a[v],tol):bad.append(str(key)+':'+v+':'+r[v]+' vs '+a[v])
   x,y=number(r[v]),number(a[v])
   if math.isfinite(x) and math.isfinite(y):ds.append(abs(x-y))
 check(name,not bad,f'{n} values; tolerance {tol}; max absolute difference {max(ds,default=0)}; mismatches '+ '; '.join(bad[:12]))
 return bad

ex=read(B/'cohort_duplicates_toray_dedup_check/evaluation_exclusions.tsv');excluded={r['row_name'] for r in ex};toray='TORAY_PUBLIC_CANONICAL';prj='PRJEB30542';units={toray,prj}
base=B/'score_identity_identity/out/analysis/related_software'
old=read(base/'results.tsv');new=read(O/'results.tsv')
assert len(new)==6 and all(r['unit']==toray and r['status']=='completed' for r in new)
bytool={r['tool']:r for r in new}
combined=[bytool[r['tool']] if r['unit']==toray else r for r in old]
write(O/'related_software_results.tsv',combined,list(old[0]))
check('GSE94533_results_unchanged',[r for r in combined if r['unit']!=toray]==[r for r in old if r['unit']!=toray],'All non-Toray fields retained from supplied release results; no GSE94533 refit')
check('related_seeds_unchanged',all(r['fit_seed']==bytool[r['tool']]['fit_seed'] and r['audit_seed']==bytool[r['tool']]['audit_seed'] for r in old if r['unit']==toray),'All six fit and audit seeds; original ci=2 retained')
a=read(O/'split_assignments.tsv');check('excluded_absent_from_related_pool',not ({r['specimen'] for r in a}&excluded),f'{len(a)} retained pool profiles; 394 removed')
check('related_train_test_cap',sum(r['role']=='train' for r in a)==1200 and sum(r['role']=='test' for r in a)==800,'2,000 selected; 1,200 train and 800 test')
rows=[]
for f,keys in [('results.tsv',['unit','tool']),('cohorts.tsv',['unit']),('agreement_ktsp.tsv',['unit']),('audit_challenges.tsv',['unit','tool','challenge','instance']),('missing_feature.tsv',['unit','tool','variant'])]:
 rows+=compare('related_software',read(base/_public_analysis_directory(f)),read(O/_public_analysis_directory(f)),keys,f,lambda r:r['unit']==toray)
# A stable common schema retains all emitted scalar summaries and challenge outcomes.
fields=['analysis','table','unit','tool','challenge','instance','variant','quantity','old','new','delta','changed']
write(O/'related_software_toray_reference_comparison.tsv',({k:r.get(k,'') for k in fields} for r in rows),fields)
for f in ['cohorts.tsv','agreement_ktsp.tsv','audit_challenges.tsv','missing_feature.tsv']:
 rr=read(base/_public_analysis_directory(f));merged=[r for r in rr if r['unit']!=toray]+read(O/_public_analysis_directory(f));write(O/(_public_analysis_directory('related_software_'+f)),merged,list(rr[0]))
# The true producer, not rounded TSV scores, supplies all AUCs; precision receipts validate replay.
receipt=read(O/'full_precision_receipt.tsv');check('related_full_precision_replay',len(receipt)==6 and all(r['heldout_decisions_unchanged']=='TRUE' and r['singleton_decisions_unchanged']=='TRUE' and same(r['reported_auc'],r['replayed_full_precision_auc']) for r in receipt),'Six frozen-model replays and exact native decision checks')

for f,keys in [('operating_points.tsv',['context','cohort','method','seed','target_specificity']),('calibration.tsv',['context','cohort','method','seed','calibrator']),('decision_curves.tsv',['context','cohort','method','calibrator','target_prevalence','risk_threshold'])]:
 src=B/'score_identity_identity/out/analysis/decision_curves';orig=read(O/'decision_curves/original'/_public_analysis_directory(f));filt=read(O/'decision_curves/filtered'/_public_analysis_directory(f));saved=read(src/_public_analysis_directory(f))
 selected=lambda rr:[r for r in rr if r['cohort'] in units]
 validate_against_saved('decision_'+f,selected(orig),selected(saved),keys)
 write(O/(_public_analysis_directory(f.stem+'_reference_comparison.tsv')) if isinstance(f,Path) else O/(_public_analysis_directory(f[:-4]+'_reference_comparison.tsv')),compare('decision_curves',selected(orig),selected(filt),keys,f))
 fm={tuple(r[k] for k in keys):r for r in selected(filt)}
 write(O/(_public_analysis_directory('filtered_'+f)),[fm.get(tuple(r[k] for k in keys),r) for r in saved],list(saved[0]))
for f,keys in [('results.tsv',['outcome','tool','scope','stratum']),('auc.tsv',['tool','unit','platform_family','split']),('decision_rule.tsv',['part','scope','stratum','comparator'])]:
 a=read(O/'shared_fold_unusable/original'/_public_analysis_directory(f));b=read(O/'shared_fold_unusable/filtered'/_public_analysis_directory(f))
 validate_against_saved('shared_'+f,a,read(B/'score_identity_identity/out/analysis/shared_fold_unusable'/_public_analysis_directory(f)),keys)
 write(O/(_public_analysis_directory('shared_fold_'+f[:-4]+'_reference_comparison.tsv')),compare('shared_fold_unusable',a,b,keys,f))
for p in sorted(O.glob('*_original_metrics.tsv')):
 name=p.name.removesuffix('_original_metrics.tsv');keys=['unit'] if name.startswith('mirbit') else ['cohort','method'] if name=='hemolysis' else ['scenario','seed','method','cohort'] if name=='panel_restriction' else ['seed','method','cohort']
 write(O/(_public_analysis_directory(name+'_reference_comparison.tsv')),compare(name,read(p),read(O/(_public_analysis_directory(name+'_filtered_metrics.tsv'))),keys,'metrics'))
p=O/'panel_restriction_original_comparisons.tsv'
if p.exists():
 write(O/'panel_restriction_comparisons_reference_comparison.tsv',compare('panel_restriction',read(p),read(O/'panel_restriction_filtered_comparisons.tsv'),['scenario','seed','method','cohort'],'comparisons'))
 validate_against_saved('panel_comparisons',read(p),read(B/'input_robustness_analysis/out/part_b/analysis/panel_restriction/panel_sensitivity_comparisons.tsv'),['scenario','seed','method','cohort'],tol=1e-7)
write(O/'validation_checks.tsv',checks)
# Reconcile the older rounded miRBiT transport output with the complete fresh score replay.
arch=read(ROOT/'analysis/reproduction_secondary/part_a/mirbit/transport.tsv')
fresh={r['unit']:r for r in read(O/'mirbit_original_metrics.tsv')};filtered={r['unit']:r for r in read(O/'mirbit_filtered_metrics.tsv')}
rec=[]
for r in arch:
 if r['dataset'] not in units:continue
 for m in ('AUC','cancer_calls','total_cancer'):
  u=r['dataset'];rec.append(dict(unit=u,quantity=m,archived_reported=r[m],fresh_prefilter=fresh[u][m],fresh_filtered=filtered[u][m],fresh_exclusion_delta=number(filtered[u][m])-number(fresh[u][m]),archived_exclusion_only_new='NA',limitation='Archived extracted inputs lack feature rows; fresh replay has a different baseline and is not an exclusion-only correction of the archived number.'))
write(O/'mirbit_baseline_reconciliation.tsv',rec)
# Hemolysis validation uses the original log-index rows; unavailable estimates stay unavailable.
hs=read(ROOT/'analysis/opportunity/strengthening/platform_concordance/hemolysis_strat.tsv')
ho=read(O/'hemolysis_original_metrics.tsv');hv=[]
for r in ho:
 a=[z for z in hs if z['cohort']==r['cohort'] and (z['record_type']=='index_versus_label' if r['method']=='index' else z['record_type']=='cohort_hemolysis' and z['method']==r['method'])]
 if not a:continue
 for k in ('n_profiles','n_case','n_control','auc_all','auc_below_median','auc_at_or_above_median','spearman_score_vs_index','spearman_p'):
  if a[0].get(k,'') not in ('','NA'):
   hv.append(dict(cohort=r['cohort'],method=r['method'],quantity=k,saved=a[0][k],recomputed=r[k],pass_=same(a[0][k],r[k],5.1e-5 if r['method']=='index' else 1e-10)))
write(O/'hemolysis_validation.tsv',hv);check('hemolysis_index_original_reproduced',all(r['pass_'] for r in hv if r['method']=='index'),'Original S3 index association reproduced; archived AUC has four decimals. Ancillary disease-score differences retained in hemolysis_validation.tsv.')
src=B/'score_identity_identity/out/analysis/unusable_inputs/results.tsv';rr=read(src);target=[r for r in rr if r['unit']==toray]
write(O/'unusable_inputs_reference_comparison.tsv',compare('unusable_inputs',target,target,['unit','tool','policy','condition','batch_or_single'],'results.tsv'))
# reference rank-tree postprocessing is a separate lineage, never a release-fit substitute.
rh=O/'ranktree_reference';sh=ROOT/'analysis/ranktree_external_competitor_sensitivity/final_execution_amend3_6ae8c249_run2/assembled'
keymap={
'method_unit_performance':['unit_id','method_id'],
'method_stratum_performance':['method_id','unit_id','seed','fold','split_id','stratum_id'],
'pairwise_unit':['pair_id','unit_id'],
'pairwise_strata':['pair_id','unit_id','stratum_id'],
'pairwise_meta':['analysis_population','pair_id'],
'coverage_performance':['analysis_population','method_id'],
'complete_support_table':['analysis_population','threshold','item_type','method_id','unit_id'],
'complete_support_rectangle':['analysis_population','threshold','method_id','unit_id']}
if (rh/'filtered').exists():
 for name,keys in keymap.items():
  a=read(rh/'original'/(_public_analysis_directory(name+'.tsv')));b=read(rh/'filtered'/(_public_analysis_directory(name+'.tsv')))
  savedname={'complete_support_table':'complete_support_panels','complete_support_rectangle':'complete_support_performance'}.get(name,name)
  validate_against_saved('historical_ranktree_'+name,a,read(sh/(_public_analysis_directory(savedname+'.tsv'))),keys)
  write(O/(_public_analysis_directory('ranktree_historical_'+name+'_reference_comparison.tsv')),compare('ranktree_historical',a,b,keys,name))
write(O/'validation_checks.tsv',checks)
print('Validation:',sum(r['pass_'] for r in checks),'/',len(checks))
for r in checks:
 if not r['pass_']:print(r)

for name,basepath,keys in [
 ('panel_restriction',B/'input_robustness_analysis/out/part_b/analysis/panel_restriction/panel_sensitivity_cells.tsv',['scenario','seed','method','cohort']),
 ('panel_restriction_comparisons',B/'input_robustness_analysis/out/part_b/analysis/panel_restriction/panel_sensitivity_comparisons.tsv',['scenario','seed','method','cohort'])]:
 p=O/(_public_analysis_directory(name+'_reference_comparison.tsv'));rr=read(p);saved={tuple(z.get(k,'') for k in keys):z for z in read(basepath)}
 for z in rr:
  sr=saved[tuple(z.get(k,'') for k in keys)];m=z['quantity'];before=z['old'];z['recomputed_before_filter']=before;z['original_reported']=sr.get(m,'NA')
  z['baseline_serialization_difference']=number(before)-number(sr.get(m,'NA'))
  if m in sr:z['old']=sr[m]
  x,y=number(z['old']),number(z['new']);z['delta']=y-x if math.isfinite(x) and math.isfinite(y) else 'NA'
  z['score_precision']='saved decimal TSV; native release checkpoints absent'
 write(p,rr)
# Reported hemolysis index is the reassay_design quantity. Ancillary disease-score strata retain the replay discrepancy.
reference_hemolysis=read(ROOT/'analysis/opportunity/strengthening/platform_concordance/hemolysis_strat.tsv')
rows=read(O/'hemolysis_reference_comparison.tsv')
for z in rows:
 a=[r for r in reference_hemolysis if r['cohort']==z['cohort'] and (r['record_type']=='index_versus_label' if z['method']=='index' else r['record_type']=='cohort_hemolysis' and r['method']==z['method'])]
 z['recomputed_before_filter']=z['old'];z['original_reported']=a[0].get(z['quantity'],'NA') if a else 'NA'
 z['baseline_replay_difference']=number(z['old'])-number(z['original_reported'])
 if z['original_reported'] not in ('','NA'):z['old']=z['original_reported']
 x,y=number(z['old']),number(z['new']);z['delta']=y-x if math.isfinite(x) and math.isfinite(y) else 'NA'
write(O/'hemolysis_reference_comparison.tsv',rows)
