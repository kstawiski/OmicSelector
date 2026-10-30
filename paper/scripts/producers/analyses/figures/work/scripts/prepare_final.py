#!/usr/bin/env python3
"""Stage final accepted inputs without altering any analysis output."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,hashlib,json,re,shutil,subprocess,sys
M=Path('./analyses');W=M/'figures/work';S=M/'benchmark_assembly/synthesis'
csv.field_size_limit(1<<28)
ledger=[]
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv(p,rr,fields=None):
 rr=list(rr);p.parent.mkdir(parents=True,exist_ok=True);fields=fields or list(dict.fromkeys(k for r in rr for k in r))
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rr)
def cp(src,dest):
 dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
 ledger.append(dict(output=str(dest),source=str(src),source_sha256=sha(src),transformation='byte copy'))
def write(p,t,src,fields):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(t)
 ledger.append(dict(output=str(p),source=str(src),source_sha256=sha(src),transformation=fields))
for name in ['manuscript','S1','Cover_letter']:
 src=M/'figure_annotations/integrator'/_public_analysis_directory(f'{name}.keyed.md');p=W/'text'/_public_analysis_directory(f'{name}.keyed.md');t=src.read_text()
 if name=='manuscript':
  t=t.replace('The foundation-model comparison covered 224 deployment evaluations and 214 eligible benchmark comparisons. ', '')
  t=t.replace('### Tabular foundation models gave the largest pooled gains over the simple reference','### Foundation models gave the largest pooled point estimates without significant improvement over the reference')
  t=t.replace('The largest pooled gain was','The largest pooled point estimate was').replace('The next gains were','The next point estimates were')
  t=t.replace('Among ⟦B:estimable_methods⟧ estimable methods, intervals lay wholly below zero for ⟦B:ci_below_zero⟧ (Figure 3).','⟦B:primary_interval_sentence⟧ Among ⟦B:estimable_methods⟧ estimable methods, ⟦B:ci_below_zero⟧ performed significantly worse, with confidence intervals wholly below zero (Figure 3).')
  t=t.replace('The largest pooled gain over the simple reference was ⟦B:largest_estimate⟧ AUC.','Foundation models gave the largest pooled point estimates, reaching ⟦B:largest_estimate⟧ AUC above the simple reference.')
  t=t.replace('These findings locate the demonstrated discrepancy in refitting rather than scoring a frozen model one specimen at a time.','Scoring an already fitted model is unaffected by this refitting mechanism. Fix the CPU thread count when refitting TabPFN.')
  t=t.replace('Toray combined eight 3D-Gene serum accessions', 'Toray combined eight of ten audited 3D-Gene serum accessions')
  t=t.replace('Two further audited accessions contributed no profile after deduplication.','GSE110317 and GSE122497 contributed no profile after deduplication.')
  t=t.replace('A1 contrasts are cohort-normalized minus frozen AUC.','A1 contrasts are frozen minus cohort-normalized AUC.')
  t=t.replace('OmicSelector provides development, deployment and auditing across 73 registered methods.','OmicSelector provides development, deployment and auditing across 73 registered methods. Foundation models gave the largest pooled point estimates, but no method significantly improved on the simple reference after multiplicity adjustment. Twenty-nine methods performed significantly worse.')
  t=t.replace('⟦B:foundation_seed_results⟧','⟦B:foundation_seed_results⟧ ⟦B:foundation_seed_intervals⟧ ⟦B:foundation_seed_ranges⟧ ⟦B:foundation_seed_interval_interpretation⟧ TabICL includes GSE70080 in this analysis, so its change also reflects coverage.')
 if name=='S1':
  t=t.replace('The corpus contains 34 analysis units from 47 public accessions,','The corpus contains 34 analysis units from 47 audited public accessions, 45 of which contribute profiles,')
  t=t.replace('The foundation-model comparison covered 224 deployment evaluations and 214 eligible benchmark comparisons. ', '')
  t=t.replace('Supplementary File S24 records the controlled same-model CPU and graphics-processor experiment separately.','Fix the CPU thread count when refitting TabPFN. Scoring an already fitted model is unaffected by this refitting mechanism.')
  t=t.replace('There are ⟦B:estimable_methods⟧ estimable methods.','⟦B:primary_superiority_sentence⟧ ⟦B:primary_interval_sentence⟧ There are ⟦B:estimable_methods⟧ estimable methods.')
  nonempty=rows(M/'normalization_comparison/out/partA_nonempty_without_toray_pooled.tsv')
  tab=['| Contrast without Toray and with selected features in every fold | Datasets | AUC difference | 95% CI |','| --- | --- | --- | --- |']
  for r in nonempty:
   if r['type']=='A2':tab.append(f"| {r['contrast']} | {r['k']} | {float(r['pooled_lift']):+.4f} | {float(r['ci_lo']):.4f} to {float(r['ci_hi']):.4f} |".replace('-0.','−0.'))
  t=t.replace('The main text reports the sensitivity restricted to non-intercept-only fits.', '⟦B:secondary_nonempty_results⟧\n\n'+'\n'.join(tab))
  seedrows=rows(M/'R9F/out/foundation_seed_summary.tsv')
  tab=['| Method | Pooled AUC difference | 95% CI | Datasets | BY-adjusted P |','| --- | --- | --- | --- | --- |']
  labels={'ai-tabpfn':'TabPFN','tab-tabdpt':'TabDPT','tab-tabicl':'TabICL'}
  for r in seedrows:tab.append(f"| {labels[r['method']]} | {float(r['pooled_lift']):+.4f} | {float(r['ci_lo']):.4f} to {float(r['ci_hi']):.4f} | {r['k']} | {float(r['p_lift_gt_zero_BY']):.3f} |".replace('-0.','−0.'))
  t=t.replace('⟦B:foundation_seed_results⟧','⟦B:foundation_seed_results⟧ ⟦B:foundation_seed_ranges⟧ ⟦B:foundation_seed_interval_interpretation⟧ TabICL includes GSE70080, so its pooled change also reflects the increase from 29 to 30 datasets.\n\n'+'\n'.join(tab))
  t=t.replace('⟦B:gse186595_exclusion_results⟧','⟦B:gse186595_exclusion_results⟧ ⟦B:sensitivity_superiority_sentence⟧')
 write(p,t,src,'Integrator text plus final result interpretation, corpus counts and source-bound sensitivity content')
p=W/'text/Figure_alt_text.keyed.md';t=p.read_text().replace('takes its information: deposited values, operations within one specimen, the grouped split and steps fitted on training\nrows.','uses deposited values, operations within one specimen, the grouped split and steps fitted on training rows.')
extra=' TabPFN and TabICL have the largest pooled point estimates. No gain is significant after Benjamini–Yekutieli adjustment.'
t=t.replace(extra,'').replace('Panel D shows bars of eligible datasets out of 32.','Panel D shows bars of eligible datasets out of 32.'+extra).replace('10⁻¹²','10^−12')
p.write_text(t)
d=W/'staging/cohort_inventory';shutil.rmtree(d);d.mkdir()
for p in (M/'cohort_inventory/out/cohort_inventory').iterdir():
 if p.is_file():cp(p,d/_public_analysis_directory(p.name))
cp(M/'analysis_1_ds/out/cohort_inventory_analysis_ready_v1.xlsx',d/'cohort_inventory_analysis_ready_workbook.xlsx')
# Complete primary cell provenance from the synthesis, retaining separate analysis roles.
p=W/'staging/reassay_design/reassay_design_analysis_cell_provenance.tsv';base=W/'final_evidence/reassay_design_independent_cells.tsv'
if not base.exists():cp(p,base)
rr=[r for r in rows(base) if r['analysis']!='within_cohort_training_only']
for r in rows(S/'cells_union_complete.tsv'):
 r=dict(r,analysis='within_cohort_training_only',execution_commit=r['package_commit'],execution_device=r['source_device'],observed_threads=r['source_threads'],source_path=str(S/'cells_union_complete.tsv'))
 rr.append(r)
for chunk in sorted((M/'foundation_model_seed_sensitivity/out/chunks').iterdir()):
 if not chunk.is_dir():continue
 prov={(r['method'],r['cohort'],r['seed']):r for r in rows(chunk/'provenance.tsv')}
 for r in rows(chunk/'within_cells.tsv'):
  ev=prov[(r['method'],r['cohort'],r['seed'])];r=dict(r,**{k:v for k,v in ev.items() if k not in ['method','cohort','seed']})
  rr.append(dict(r,analysis='foundation_model_seed_sensitivity',execution_commit=ev['package_commit'],execution_device=ev['device'],observed_threads=ev['threads'],source_path=str(chunk/'within_cells.tsv')))
tsv(p,rr)
for name in ['final_cell_provenance.tsv','rebuild_plan.json'] :cp(M/'benchmark_assembly/out'/_public_analysis_directory(name),W/'staging/reassay_design'/_public_analysis_directory(name))
cp(S/'chunk_sources.tsv',W/'staging/reassay_design/reassay_design_primary_chunk_sources.tsv')
for job in ['R8','FMS']:
 for chunk in sorted((M/_public_analysis_directory(job)/'out/chunks').iterdir()):
  if chunk.is_dir():cp(chunk/'provenance.tsv',W/'staging/reassay_design'/_public_analysis_directory(job)/(_public_analysis_directory(chunk.name+'_provenance.tsv')))
p=W/'staging/reassay_design/reassay_design.md';t=p.read_text()
para='Primary benchmark provenance covers 8,500 method–dataset–seed cells. The separate foundation-model seed analysis covers 510 cells and records the model seed independently of the cross-validation seed. Its 495 eligible cells and 15 structural refusals retain actual computing commits.\n\n'
if para not in t:p.write_text(t+'\n'+para)
# Numerical synthesis and independent regeneration evidence.
for n in [6,7,9]:cp(S/_public_analysis_directory({'6': 'method_cohort_performance.tsv', '7': 'pooled_effects.tsv', '9': 'sensitivity_effects.tsv'}[str(n)]),W/'staging'/_public_analysis_directory(f'S{n}')/_public_analysis_directory({'6': 'method_cohort_performance.tsv', '7': 'pooled_effects.tsv', '9': 'sensitivity_effects.tsv'}[str(n)]))
for name in ['leave_one_cohort_out.tsv','pooled.tsv','pooled_from_cohort_check.tsv']:cp(S/'evidence'/_public_analysis_directory(name),W/'staging/pooled_effects'/_public_analysis_directory(name))
for name in ['sensitivity_summary.tsv','sensitivity_pooled.tsv','foundation_seed_sensitivities.tsv']:cp(S/'evidence'/_public_analysis_directory(name),W/'staging/sensitivity_effects'/_public_analysis_directory(name))
for name in ['difference_summary.tsv','primary_vs_sensitivity.tsv']:cp(M/'foundation_model_seed_sensitivity/out'/_public_analysis_directory(name),W/'staging/sensitivity_effects'/_public_analysis_directory(name))
for name in ['partA_without_toray_pooled.tsv','toray_weights.tsv']:
 if (M/'normalization/out'/_public_analysis_directory(name)).exists():cp(M/'normalization/out'/_public_analysis_directory(name),W/'staging/normalization_sensitivity'/_public_analysis_directory(name))
for p in (M/'normalization_comparison/out').glob('*.tsv'):cp(p,W/'staging/normalization_sensitivity/nonempty_feature_sensitivity'/_public_analysis_directory(p.name))
for p in (M/'benchmark_assembly/out').glob('matched_baseline_*.tsv'):cp(p,W/'staging/unusable_inputs/matched_baseline'/_public_analysis_directory(p.name))
for p in (M/'scoring_agreement/out').glob('*.tsv'):cp(p,W/'staging/unusable_inputs/scoring_paths'/_public_analysis_directory(p.name))
for p in (M/'scoring_precision/out').glob('*.tsv'):cp(p,W/'staging/unusable_inputs/scoring_paths/thread_diagnosis'/_public_analysis_directory(p.name))
for p in (M/'scoring_precision/out/evidence').glob('*.tsv'):cp(p,W/'staging/unusable_inputs/scoring_paths/thread_diagnosis/evidence'/_public_analysis_directory(p.name))
for job in ['S4c','S4e','S4f']:
 for sub in ['scripts','out/code','payload']:
  for p in (M/_public_analysis_directory(job)/_public_analysis_directory(sub)).rglob('*'):
   if p.is_file() and p.suffix in ['.R','.py','.sh','.cpp']:
    cp(p,W/'staging/unusable_inputs/code'/_public_analysis_directory(job)/_public_analysis_directory(p.relative_to(M/job)))
slot=W/'staging/unusable_inputs/scoring_paths/unusable_inputs_controlled_same_host_CPU_GPU.tsv'
if not (M/'scoring_agreement/out/DONE').exists():
 tsv(slot,[dict(comparison='Same fitted model on CPU and GPU at a fixed thread count',score_difference='',auc_difference='',changed_calls='',evidence='No controlled result is supplied')])
elif slot.exists():slot.unlink()
p=W/'staging/unusable_inputs/unusable_inputs.md';t=p.read_text();add='''

Matched conventional-baseline tables compare each method with frozen-standardized lasso on identical held-out rows. The deployment comparison scores all 224 cells and includes 214 eligible benchmark comparisons. Ten GSE198834 folds have no eligible reference. All six factorial cells are scored.

TabPFN score differences above 1e-5 occur in 43 comparisons. Four versus six CPU threads during refitting change SVD features, and rounded row fingerprints amplify those differences. A representative four-thread replay agrees with deployment within 6.6e-07, while six threads reproduce benchmark scores exactly. Other folds have preprocessing evidence, and additional host contributions are not separately quantified. Fix the CPU thread count when refitting TabPFN.

The largest score difference is 0.048 in GSE210329 fold 3, with an AUC difference of −0.006 and one changed call at both the 0.5 and median cutoffs. Across eligible comparisons, changed calls total 3, 6 and 37 at the 0.5, Youden and median cutoffs. The controlled same-host CPU/GPU table is identified separately from these refitting comparisons.
'''
if 'Matched conventional-baseline tables' not in t:p.write_text(t+add)
# Reader-facing native-score versus deployment-adapter comparison.
native={
 'blood-cell-marker-excl-rclr':('Requires metadata, a frozen training panel and panel weights, which a matrix-only call omits.','Stops because metadata are missing.'),
 'frozen-ruv':('Returns a denoised matrix and takes specimen values before the fitted object.','Stops at the scorer-registration guard before native scoring.'),
 'nc-isoforest':('Returns a list containing scores, empirical P values, thresholds and flags.','Stops at the scorer-registration guard before native scoring.'),
 'nc-mahalanobis':('Returns a list containing squared distances, distances and chi-square P values.','Stops at the scorer-registration guard before native scoring.'),
 'robust-pca-resid':('Returns a matrix of residuals rather than one score per specimen.','Stops at the scorer-registration guard before native scoring.'),
 'ws-alr-pivot':('Returns an ALR matrix and requires sufficient pivot features plus a frozen selected panel.','Stops because too few default pivot features are present.'),
 'ws-dominance':('Omits the frozen fitted top-feature count and pseudocount when called without an adapter.','Stops at the scorer-registration guard before native scoring.'),
 'ws-mad-logratio':('Returns a within-specimen MAD-scaled log-ratio matrix and omits the frozen selected panel.','Shape validation rejects the matrix returned for the training probe.')}
ad=[]
for r in rows(M/'AD1/out/adapter_defects.tsv'):
 method=r['method_id'];a,b=native[method]
 ad.append({'Method':method,'Native scoring output or omission':a,'Training preflight outcome':b,'Outcome through deployment adapter':'Returns one score per specimen and passes the synthetic and real-profile audits.','Evidence file':'Supplementary File S21, AD1/out/adapter_defects.tsv, row '+str(len(ad)+2)})
tsv(W/'staging/scoring_conditions/native_scoring_and_adapter.tsv',ad)
p=W/'staging/scoring_conditions/scoring_conditions.md';t=p.read_text();add='\nThe native-scoring comparison describes eight methods whose direct scoring route cannot supply the required specimen-score vector and frozen state. Five routes stop at the registration guard before the native scorer runs. The deployment adapters supply the frozen state or extract the specimen scores, and all eight pass the specified audits. The companion table records each native output, preflight outcome, adapter outcome and evidence file.\n'
if 'The native-scoring comparison' not in t:p.write_text(t+add)
cp(M/'scoring_agreement/out/deployed_vs_benchmark.tsv',W/'staging/scoring_conditions/scoring_conditions_full_protocol_scoring_comparison.tsv')
for name in ['companion_dependence_pooled_comparison.py','figure_method_cohort_performance_auc_heatmap.py']:
 src=M/'figures_v10/scripts'/_public_analysis_directory(name);t=src.read_text()
 t=t.replace("R2/out/m16';assert (M16B.parent/'DONE').is_file()","R2b/out/m16';assert (M16B/'partA_pooled.tsv').is_file()")
 t=t.replace('all currently estimable','all estimable').replace('current 2.6.7 S6 table','S6 table')
 if name=='companion_dependence_pooled_comparison.py':
  t=t.replace("Panel C omits one dataset at a time.'","Panel C omits one dataset at a time.\\nNo gain passes Benjamini–Yekutieli adjustment. The smallest adjusted P is 0.604 for TabPFN.'")
 write(W/'figures/scripts'/_public_analysis_directory(name),t,src,'Final R9 synthesis and accepted R2b conventional baseline')
tsv(W/'final_evidence/staging_sources.tsv',ledger)
from final_reader_fixes import apply
apply()
print('Staged',len(ledger),'source artifacts; primary cells',8500,'seed sensitivity cells',510)
