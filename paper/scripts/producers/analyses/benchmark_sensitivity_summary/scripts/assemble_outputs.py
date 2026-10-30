"""Benchmark summary tables, analysis comparisons, and interpretation."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,hashlib,json,sys,shutil
import numpy as np
import pandas as pd
from derived import derived, normalization_metrics, rd, save
B=Path(__file__).resolve().parents[1];M=B.parent;O=B/'out';S=O/'synthesis';N=O/'normalization';A=M/'benchmark_reference_summary/out';OLD=A/'synthesis'
GM={'reo-pairratio','inv-scatter','sig-path','proto-net','nc-conformal-anom'};GC={'GSE204951','GSE221088'}
# Read-only verification that promotion preserves held-out participants and reference scores.
import subprocess
with (O/'logs/alignment.log').open('w') as log:
 subprocess.run(['Rscript','--vanilla',str(B/'scripts/alignment.R')],check=True,stdout=log,stderr=subprocess.STDOUT)
checks=[]
def valid(name,ok,detail=''):
 checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
 if not ok:
  save(pd.DataFrame(checks),O/'validation.tsv');raise AssertionError((name,detail))
valid('grouped_outer_alignment',rd(S/'evidence/grouped_outer_alignment.tsv').pass_check.all())
ind=rd(O/'independent_check.tsv');valid('independent_checks',ind.pass_check.all(),f'{len(ind)} checks; maximum error {ind.absolute_error.max():.12g}')
primary_error=ind.loc[ind.scope=='Independent primary S6 REML/KH/BY','absolute_error'].max()
oldmetrics=derived(OLD,'old');newmetrics=derived(S,'primary');fullmetrics=derived(S/'full32_sensitivity','full')
save(oldmetrics,O/'current_headline_metrics.tsv');nm=normalization_metrics(N);save(nm,N/'printed_metrics.tsv');save(pd.concat([newmetrics,nm],ignore_index=True),O/'headline_metrics.tsv')
for tag,root,nroot in [('current',OLD,A/'normalization'),('eligible',S,N)]:
 pa=rd(nroot/'companion_sensitivity/computing_bundle/partA_pooled.tsv');ct=rd(nroot/'companion_sensitivity/computing_bundle/partA_contrasts_dataset.tsv')
 keep=['CN-STD minus FZ-STD','CN-QN minus FZ-QN','CN-COMBAT minus FZ-RAW','SS-CLR minus CN-STD','SS-CLR minus CN-QN','SS-CLR minus CN-COMBAT']
 pp=pa[pa.contrast.isin(keep)].copy();dd=ct[ct.contrast.isin(keep)].copy()
 for d,field in [(pp,'pooled_lift'),(dd,'estimate')]:
  neg=d.type.eq('A1');d.loc[neg,field]=-d.loc[neg,field];lo=d.loc[neg,'ci_lo'].copy();d.loc[neg,'ci_lo']=-d.loc[neg,'ci_hi'];d.loc[neg,'ci_hi']=-lo;d['figure_direction']='single-sample minus cohort-normalized'
 bs=rd(nroot/'companion_sensitivity/computing_bundle/partB_summary.tsv');heat=bs[bs.metric.isin(['B1','B3'])].pivot(index=['dataset','arm'],columns='metric',values='estimate').reset_index();heat['difference']=heat.B3-heat.B1
 heat['printed']=heat.difference.map(lambda x:'' if pd.isna(x) else '0.00' if f'{x:+.2f}' in ['+0.00','-0.00'] else f'{x:+.2f}')
 s7=rd(root/'pooled_effects.tsv');s7['eligible_universe']=30;lasso=pa[(pa.type=='A3')&(pa.contrast=='FZ-STD minus R')].copy();lasso['eligible_universe']=30
 for name,d in [('figure2C_pooled',pp),('figure2C_dataset',dd),('figure2D',heat),('figure3_methods',s7),('figure3_frozen_lasso',lasso)]:save(d,O/'figure_values'/_public_analysis_directory(tag)/(_public_analysis_directory(name+'.tsv')))
# Matched comparison: update the missing partial pair-ratio result and source evidence.
prov=O/'source_provenance';cells=rd(S/'cells_union_complete.tsv')
pcols=['method','cohort','seed','eligible','package_version','package_commit','analysis_code_id','source_bundle','source_layer','source_host','source_device','source_threads']
save(cells[pcols],O/'final_cell_provenance.tsv')
selected=set(rd(O/'matched_baseline_dataset.tsv').method)
src=cells[cells.method.isin(selected)][['method','cohort','seed','eligible','source_bundle','source_device','source_status']].copy()
oldsrc=rd(A/'source_provenance/benchmark_saved_scores_benchmark_sources.tsv');hashes=oldsrc[['source_bundle','source_bundle_sha256']].drop_duplicates().set_index('source_bundle').source_bundle_sha256.to_dict()
source_key=cells.loc[cells.source_layer=='person_grouped_inner_primary','source_bundle'].iloc[0];source=Path(source_key);hashes[source_key]=hashlib.sha256(source.read_bytes()).hexdigest()
src['source_bundle_sha256']=src.source_bundle.map(hashes);valid('matched_source_hashes_present',src.source_bundle_sha256.notna().all());save(src,O/'benchmark_sources.tsv')
folds=rd(A/'source_provenance/benchmark_saved_scores_matched_fold_differences.tsv');folds=folds[~(folds.method.eq('reo-pairratio')&folds.cohort.isin(GC))];save(folds,O/'matched_fold_differences.tsv')
align=rd(A/'source_provenance/benchmark_saved_scores_fold_row_alignment.tsv');mask=align.method.eq('reo-pairratio')&align.cohort.isin(GC)
align.loc[mask,['method_rows','baseline_rows_for_method_seeds','merged_rows']]=0;align.loc[mask,'row_fold_group_label_train_test_match']=np.nan;align.loc[mask,'status']='method_not_evaluable';save(align,O/'fold_row_alignment.tsv')
shutil.copy2(A/'source_provenance/benchmark_saved_scores_intercept_only_folds.tsv',O/'intercept_only_folds.tsv')
shutil.copy2(A/'final_inputs.json',O/'final_inputs_benchmark_reference_summary.json')
finalinputs=json.loads((A/'final_inputs.json').read_text());finalinputs['R12c_grouped_primary']=dict(methods=sorted(GM),cohorts=sorted(GC),cells=50,bundle=str(source),bundle_sha256=hashes[source_key],source_package_version='2.6.7',source_commit='633d7d5a645323c14e8ed874bde28eb59c58f7e6',refits=0,sensitivity='Ungrouped inner resampling sensitivity (R12a primary cells)')
(O/'final_inputs.json').write_text(json.dumps(finalinputs,indent=2)+'\n')
p=rd(S/'evidence/pooled.tsv');ss=rd(S/'evidence/sensitivity_pooled.tsv');w=rd(S/'evidence/within_performance.tsv');ow=rd(OLD/'evidence/within_performance.tsv')
valid('primary_30_units',w.cohort.nunique()==30 and not w.cohort.isin(['GSE270497','GSE186595','GSE145176','GSE31568']).any())
valid('primary_1500_rows_7500_cells',len(w)==1500 and len(cells)==7500 and not cells.duplicated(['method','cohort','seed']).any())
changed=w.method.isin(GM)&w.cohort.isin(GC)
cols=[c for c in w if c not in ['p_lift_gt_zero_BY','p_lift_gt_margin_0_05_BY']]
x=w.loc[~changed,cols].set_index(['method','cohort']).sort_index();y=ow.loc[~(ow.method.isin(GM)&ow.cohort.isin(GC)),cols].set_index(['method','cohort']).sort_index()
valid('unaffected_cohort_statistics_unchanged',x.equals(y))
valid('complete_grouped_counterparts',len(rd(S/'evidence/grouped_cells_source.tsv'))==50)
valid('40_eligible_10_refused_grouped_cells',int(rd(S/'evidence/grouped_cells_source.tsv').eligible.sum())==40)
valid('grouped_predictions_only_complete_cells',not rd(S/'evidence/grouped_predictions.tsv').method.eq('reo-pairratio').any())
valid('S9_12_blocks',ss.analysis.nunique()==12 and ss.groupby('analysis').size().eq(49).all())
fam=ss.groupby('analysis').agg(method_rows=('method','size'),above_zero_finite_tests=('p_lift_gt_zero','count'),above_margin_finite_tests=('p_lift_gt_margin_0_05','count')).reset_index();save(fam,O/'analysis_families.tsv')
valid('S9_both_48_test_families',fam.above_zero_finite_tests.eq(48).all() and fam.above_margin_finite_tests.eq(48).all())
valid('primary_weights_sum_100',np.allclose(rd(S/'evidence/cohort_weights.tsv').groupby('method').weight_percent.sum(),100,atol=1e-10))
valid('all_required_R12a_synthesis_paths_reissued',all((S/_public_analysis_directory(f.relative_to(OLD))).exists() for f in OLD.rglob('*') if f.is_file()))
for branch in ['', 'full32_sensitivity']:
 q=rd(S/_public_analysis_directory(branch)/'ungrouped_inner_sensitivity/pooled.tsv').set_index('method');ref=rd(OLD/_public_analysis_directory(branch)/'evidence/pooled.tsv').set_index('method')
 cols=ref.select_dtypes('number').columns
 valid('ungrouped_replays_R12a_'+(branch or 'primary'),q[cols].isna().equals(ref[cols].isna()) and np.nanmax(np.abs(q[cols]-ref[cols]))<=1e-6)
 newmatched=rd(S/_public_analysis_directory(branch)/'evidence/matched_baseline_pooled.tsv').set_index('method');oldmatched=rd(OLD/_public_analysis_directory(branch)/'evidence/matched_baseline_pooled.tsv').set_index('method')
 valid('matched_pool_invariance_'+(branch or 'primary'),np.nanmax(np.abs(newmatched.select_dtypes('number')-oldmatched.select_dtypes('number')))<=1e-10)
valid('matched_grouped_pairratio_no_partial_folds',not folds.loc[folds.method.eq('reo-pairratio'),'cohort'].isin(GC).any())
naudit=[]
for f in (A/'normalization').rglob('*.tsv'):
 rel=f.relative_to(A/'normalization');a=rd(f);b=rd(N/_public_analysis_directory(rel))
 numeric=a.select_dtypes('number').columns;err=(a[numeric]-b[numeric]).abs().max().max() if len(numeric) else 0.
 ok=a.shape==b.shape and a[numeric].isna().equals(b[numeric].isna()) and (pd.isna(err) or err<=1e-10)
 valid('normalization_invariance_'+str(rel),ok,err);naudit.append(dict(table=str(rel),max_absolute_error=err,pass_check=ok,reason='Normalization arms and the rCLR reference use the same methods in both comparisons'))
save(pd.DataFrame(naudit),O/'normalization_impact.tsv')
# Headline values include an explicit benchmark_reference_summary value for each reissued numeric field.
head=[]
for label,path,oldpath,keys in [('primary pooled versus rCLR',S/'evidence/pooled.tsv',OLD/'evidence/pooled.tsv',['method']),('matched frozen-lasso pooled',O/'matched_baseline_pooled.tsv',A/'matched_baseline_pooled.tsv',['method']),('sensitivity pooled versus rCLR',S/'evidence/sensitivity_pooled.tsv',OLD/'evidence/sensitivity_pooled.tsv',['analysis','method'])]:
 d=rd(path);a=rd(oldpath)
 if 'analysis' in a:a.analysis=a.analysis.replace({'Primary eligible-set synthesis (30 units)':'Primary eligible-set synthesis (30 units; person-grouped inner)','Person-grouped inner resampling':'Ungrouped inner resampling sensitivity (R12a primary cells)'})
 nums=[c for c in d.select_dtypes('number').columns if c in a];a=a[keys+nums].rename(columns={c:'R12a_'+c for c in nums})
 d=d.merge(a,on=keys,how='left');d.insert(0,'metric',label);d['evidence']=str(path.relative_to(O));head.append(d)
save(pd.concat(head,ignore_index=True),O/'headline_results.tsv')
save(newmetrics[newmetrics.metric.isin(['near_zero_count','near_zero_list'])],S/'evidence/near_zero_summary.tsv');save(p[(p.ci_hi>-.01)&(p.ci_hi<0)],S/'evidence/near_zero_methods.tsv')
printed=rd(A/'printed_values.tsv');printed['old']=printed['new']
for i,r in printed.iterrows():
 if str(r.new).startswith('The eligible-set rebuild recorded '):
  n=len(rd(S/'evidence/reml_diagnostics.tsv'));nb=rd(S/'evidence/reml_diagnostics.tsv').reml_boundary.sum();printed.at[i,'new']=f'The grouped-inner eligible-set rebuild recorded {n:,} REML requests, including {nb:,} estimates at the recorded zero-boundary threshold. This is an execution count, not a count of independent analyses.'
 printed.at[i,'new']=str(printed.at[i,'new']).replace('Primary eligible-set synthesis (30 units)','Primary eligible-set synthesis (30 units; person-grouped inner)').replace('Person-grouped inner resampling','Ungrouped inner resampling sensitivity (R12a primary cells)')
extra=dict(source='S1 Grouped inner tuning',line='',old='Person-grouped inner resampling was an exploratory sensitivity.',new='Person-grouped inner resampling is primary for pair-ratio, scattering, path signature, prototype-network episodes and conformal-anomaly calibration in GSE204951 and GSE221088. The former ungrouped cells form the named sensitivity. Pair-ratio remains unevaluable in both datasets.',evidence='synthesis/evidence/grouped_execution_coverage.tsv; synthesis/sensitivity_effects.tsv',scope='design definition')
printed=pd.concat([printed,pd.DataFrame([extra])],ignore_index=True);save(printed,O/'printed_values.tsv')
# Exact field-level ledger over every analysis table, including seed cells, figures and normalization.
LABELS={'Primary eligible-set synthesis (30 units)':'Primary eligible-set synthesis (30 units; person-grouped inner)','Person-grouped inner resampling':'Ungrouped inner resampling sensitivity (R12a primary cells)'}
rows=[];manifest=[]
def equal(a,b):
 if a==b:return True
 if a.lower() in ['true','false'] and b.lower() in ['true','false']:return a.lower()==b.lower()
 if a in ['','NA','NaN','nan'] and b in ['','NA','NaN','nan']:return True
 try:return float(a)==float(b)
 except ValueError:return False
def compare(old,new,label):
 a=pd.read_csv(old,sep='\t',dtype=str,keep_default_na=False) if old.exists() else pd.DataFrame()
 b=pd.read_csv(new,sep='\t',dtype=str,keep_default_na=False)
 for d in [a,b]:
  for c in ['analysis','Sensitivity analysis']:
   if c in d:d[c]=d[c].replace(LABELS)
 candidates=['analysis','Sensitivity analysis','method','Method identifier','cohort','Cohort','Dataset','dataset','unit','part','type','contrast','metric','expectation','arm','model','seed','fold','sample_id','row_id','group_id','repetition','stage','batch','fraction','omitted_cohort','modality','field','scope','context','check','source','line']
 if new.name=='method_cohort_performance.tsv':keys=list(b.columns[:2])
 elif new.name=='pooled_effects.tsv':keys=[b.columns[0]]
 elif new.name=='sensitivity_effects.tsv':keys=list(b.columns[:2])
 else:keys=[c for c in candidates if c in b and (c in a or a.empty)]
 if not keys:
  a['_row']=np.arange(len(a)).astype(str);b['_row']=np.arange(len(b)).astype(str);keys=['_row']
 if a.empty:a=pd.DataFrame(columns=b.columns)
 for d in [a,b]:d['_occurrence']=d.groupby(keys,dropna=False,sort=False).cumcount().astype(str)
 keycols=keys+['_occurrence'];aa={tuple(r[c] for c in keycols):r for r in a.to_dict('records')};bb={tuple(r[c] for c in keycols):r for r in b.to_dict('records')}
 fields=[c for c in dict.fromkeys(list(a.columns)+list(b.columns)) if c not in keycols];start=len(rows)
 for key in dict.fromkeys(list(aa)+list(bb)):
  ar=aa.get(key);br=bb.get(key)
  for f in fields:
   av=ar.get(f,'') if ar is not None else '';bv=br.get(f,'') if br is not None else ''
   if ar is not None and br is not None and f in ar and f in br and equal(av,bv):continue
   try:delta=float(bv)-float(av)
   except ValueError:delta=''
   status='added' if ar is None or f not in ar else 'removed' if br is None or f not in br else 'rounding' if isinstance(delta,float) and abs(delta)<=1e-12 else 'changed'
   rows.append(dict(table=label,row_key=json.dumps(dict(zip(keycols,key)),sort_keys=True),field=f,old=av,new=bv,delta=delta,status=status,old_source=str(old),new_source=str(new)))
 manifest.append(dict(table=label,reference_rows=len(a),analysis_rows=len(b),changed_fields=len(rows)-start,old_source=str(old),new_source=str(new)))
for folder in ['synthesis','normalization']:
 for new in sorted((O/_public_analysis_directory(folder)).rglob('*.tsv')):
  rel=new.relative_to(O)
  # Newly added raw prediction evidence has no benchmark_reference_summary table and is not a changed estimate.
  if not (A/_public_analysis_directory(rel)).exists():continue
  compare(A/_public_analysis_directory(rel),new,str(rel))
for f in ['matched_baseline_dataset.tsv','matched_baseline_pooled.tsv','headline_metrics.tsv','printed_values.tsv']:
 compare(A/_public_analysis_directory(f),O/_public_analysis_directory(f),f)
for f in (O/'figure_values/eligible').glob('*.tsv'):compare(A/'figure_values/eligible'/_public_analysis_directory(f.name),f,str(f.relative_to(O)))
# The renamed sensitivity label is itself a change, separately from its numerical row matching.
for old,new in LABELS.items():rows.append(dict(table='analysis_labels',row_key=old,field='analysis',old=old,new=new,delta='',status='changed',old_source=str(OLD/'sensitivity_effects.tsv'),new_source=str(S/'sensitivity_effects.tsv')))
pass;pass
# Exact numeric primary shifts for clinical readers.
a=rd(OLD/'evidence/pooled.tsv').set_index('method');z=p.set_index('method');delta=z.pooled_lift-a.pooled_lift;changed_methods=delta[delta.abs()>1e-12].index
sents=[]
for m in changed_methods:sents.append(f'{m}: {a.loc[m,"pooled_lift"]:+.6f} → {z.loc[m,"pooled_lift"]:+.6f}')
maxm=delta.abs().idxmax();e=p[p.pooled_lift.notna()]
report=f'''The primary analysis uses person-grouped inner resampling for five methods in GSE204951 and GSE221088. The 50 source cells include 40 eligible cells. Pair-ratio refused all ten grouped cells and remains unevaluable in both datasets. The ungrouped GSE221088 seed is included only in sensitivity analysis.

Pooled AUC differences versus rCLR compare ungrouped sensitivity with grouped primary analysis: {'; '.join(sents)}. The largest absolute difference is {abs(delta.loc[maxm]):.6f}. Other primary point estimates agree across these analyses. Pooled estimates, intervals and weights are in [headline_results.tsv](headline_results.tsv).

The primary set includes 30 units. No method passes either 48-test Benjamini–Yekutieli family. The largest pooled difference is TabPFN {z.loc['ai-tabpfn','pooled_lift']:+.4f} (95% CI {z.loc['ai-tabpfn','ci_lo']:+.4f} to {z.loc['ai-tabpfn','ci_hi']:+.4f}). S9 includes the ungrouped and 32-unit sensitivities.

All six matched frozen-lasso pools and both six-test families agree between analyses. Pair-ratio partial-fold descriptions are outside the primary analysis. Normalization estimates agree because both analyses use the same methods and reference.

Independent Python calculations reproduce primary estimates and both BY families within {primary_error:.2g}; all {len(ind):,} checks pass. No model was refitted. Existing grouped cells used OmicSelector 2.6.7; the required 2.6.9 tarball was hash-verified. This run used CPU, nice 10 and one numerical thread. The source computing host is unavailable. Prepared-input and repeated-cross-validation limitations remain. This calculation does not establish clinical or manuscript readiness.
'''
valid('report_300_words',len(report.split())<=300,len(report.split()))
(O/'REPORT.md').write_text(report);(S/'REPORT.md').write_text('See [R12c results and limitations](../REPORT.md).\n');(S/'results_summary.md').write_text(report.replace('](','](../'))
(S/'table_notes.md').write_text('The primary design uses person-grouped inner resampling for reo-pairratio, inv-scatter, sig-path, proto-net and nc-conformal-anom in GSE204951 and GSE221088. Other cells use the deposited cohort estimates. within_grouped.tsv contains the grouped primary analysis. within_ungrouped.tsv and ungrouped_inner_sensitivity/ contain the ungrouped sensitivity analysis. grouped_changes.tsv defines delta as grouped primary minus ungrouped sensitivity. Both 48-test pooled BY families are recomputed separately per S9 block; S6 has separate within-cohort global BY families. The 30-unit eligibility set and R12a REML/KH estimator, corrected repeated-CV SEs, and finite-positive-SE gates are retained. Deposited source records identify the matched sources and cells. Pair-ratio refusals are observed outcomes, not missing compute.\n')
(N/'README.md').write_text('Normalization uses the deposited estimators and schemas. No affected method enters these arms or the unchanged rCLR reference. Every table was checked against R12a; see ../normalization_impact.tsv. Primary remains 30 units for Part A and 11 planned units for deployment. Full32, without-Toray, nonempty-fit and deposited-input-restriction outputs are retained.\n')
(O/'REBUILD.md').write_text('Run the public summary entry point with the checksum-verified release bundle. Numerical contracts compare the generated cohort and pooled estimates, uncertainty intervals, multiplicity adjustments and sensitivity analyses with deposited values. Independent calculations implement REML, Knapp–Hartung intervals and Benjamini–Yekutieli adjustment. Prepared-input and repeated-cross-validation limitations apply.\n')
for name in ['harness_change.patch','internal_resampling_inventory.tsv','release_source_manifest.tsv','installed_function_evidence.tsv']:
 shutil.copy2(M.parent/'analyses/cohort_benchmark_dataset_computations/out/grouped_tuning'/_public_analysis_directory(name),prov/(_public_analysis_directory('grouped_'+name)))
save(pd.DataFrame(checks),O/'validation.tsv');print(report)
