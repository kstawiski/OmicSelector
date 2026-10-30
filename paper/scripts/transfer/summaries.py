"""Run final deposited producers. Never substitute deposited outputs for regeneration."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import os,sys,subprocess,time,csv,shutil,hashlib,math,json
PAPER=next(p for p in Path(__file__).resolve().parents if (p/'config/project.json').is_file());ROOT=Path(os.environ['OMICSELECTOR_PAPER_ROOT']).resolve();M=ROOT/'analyses';sys.path.insert(0,str(PAPER/'scripts'))
from run import environment,model_environment,write,sha,rows
TIMES=[]
SKIP_R12A='--skip-r12a' in sys.argv
SKIP_R13='--skip-r13' in sys.argv
def run(path,task,cwd=None,env=None):
 e=environment() if task=='R12a' else model_environment();e.update(R1B=str(M/_public_analysis_directory(task)),OUT=str(M/_public_analysis_directory(task)/'out/synthesis'));e.update(env or {})
 cmd=([str(ROOT/'runtime/R/bin/exec/R'),'--vanilla','--slave','--file='+str(path)] if task=='R12a' else ['Rscript','--vanilla',str(path)]) if path.suffix=='.R' else [sys.executable,str(path)]
 start=time.monotonic();log=ROOT/'runtime'/_public_analysis_directory(f'{task}-{path.name}.log')
 with log.open('w') as f:q=subprocess.run(cmd,cwd=cwd or ROOT,env=e,stdout=f,stderr=subprocess.STDOUT)
 TIMES.append(dict(step=f'{task}/{path.name}',seconds=time.monotonic()-start,exit_code=q.returncode));write(ROOT/'verification/final_times.tsv',TIMES)
 if q.returncode:raise RuntimeError(str(path)+' failed; see '+str(log))
for task in ['R12a','R13','R14']:
 for n in ['out/logs','out/evidence','out/synthesis/logs','work']:(M/_public_analysis_directory(task)/_public_analysis_directory(n)).mkdir(parents=True,exist_ok=True)
checks=rows(PAPER/'expected/final_output_sha256.tsv')
if not SKIP_R12A:
 run(PAPER/'data/prepared_hashes.R','R12a')
 # Every principal benchmark_reference_summary export is removed before running its accepted producers.
 checks=rows(PAPER/'expected/final_output_sha256.tsv')
 for r in checks:
  if not SKIP_R12A and r['task']=='R12a' and ('/synthesis/' in r['path'] or '/normalization/' in r['path'] or '/figure_values/' in r['path'] or Path(r['path']).name in ['analysis_families.tsv','matched_baseline_dataset.tsv','matched_baseline_pooled.tsv','headline_metrics.tsv','headline_results.tsv']):(ROOT/_public_analysis_directory(r['path'])).unlink(missing_ok=True)
 # The identity file is input provenance, while timing and host evidence are newly written above.
 for name in ['repool.R','normalization.R','independent.py','assemble_outputs.py']:run(M/'benchmark_reference_summary/scripts'/_public_analysis_directory(name),'R12a')
 from compare_tables import compare
 r12a_checks=[]
 for r in checks:
  if not SKIP_R12A and r['task']=='R12a' and ('/synthesis/' in r['path'] or '/normalization/' in r['path'] or '/figure_values/' in r['path'] or Path(r['path']).name in ['analysis_families.tsv','matched_baseline_dataset.tsv','matched_baseline_pooled.tsv','headline_metrics.tsv','headline_results.tsv']):
   z=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']));z['path']=r['path'];r12a_checks.append(z)
 write(ROOT/'verification/benchmark_reference_summary_comparison.tsv',r12a_checks)
if not SKIP_R13:
 # These four robustness producers reuse scores and person-bootstrap weights. They perform no classifier fitting.
 for name in ['item1_case_mix_operating_points.tsv','item2_case_mix_admission_and_discordance.tsv','item3_reassay.tsv','item4_labels.tsv','item5_haemolysis.tsv','item6_assay_assay_transfer.tsv']:(M/'robustness/out'/_public_analysis_directory(name)).unlink(missing_ok=True)
 for name in ['case_mix.R','reassay.R','label_summaries.R','haemolysis.R']:run(M/'robustness/scripts'/_public_analysis_directory(name),'R13')
 run(M/'robustness/scripts/assemble_labels.py','R13')
 e=model_environment()
 for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:e[key]='1'
 start=time.monotonic()
 with (ROOT/'runtime/robustness-assay_transfer.log').open('w') as f:q=subprocess.run(['Rscript','--vanilla',str(PAPER/'scripts/reassay/assay_transfer.R')],cwd=ROOT,env=e,stdout=f,stderr=subprocess.STDOUT)
 TIMES.append(dict(step='R13/handoff.R',seconds=time.monotonic()-start,exit_code=q.returncode));write(ROOT/'verification/final_times.tsv',TIMES)
 if q.returncode:raise RuntimeError('R13 handoff failed')
# Reproduce the complete published-rule calculation from public training data and deposited target matrices.
run(M/'published_rules_analysis/src/prepare.py','R14')
r14work=M/'published_rules_analysis/work';(r14work/'inputs').mkdir(exist_ok=True)
shutil.copytree(M/'published_rules_analysis/inputs',r14work/'inputs/inputs',dirs_exist_ok=True)
(r14work/'outputs').mkdir(exist_ok=True)
run(M/'published_rules_analysis/inputs/analyze.R','R14',cwd=r14work)
r12b=M/'transfer_analysis';work=r12b/'work';(work/'inputs').mkdir(parents=True,exist_ok=True);(work/'outputs').mkdir(exist_ok=True)
shutil.copytree(r12b/'summary_payload/payload',work/'inputs/payload',dirs_exist_ok=True)
e=model_environment();e['R12_PACKAGE']=str(ROOT/'runtime/packages/2.6.9/OmicSelector')
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:e[key]='1'
start=time.monotonic()
with (ROOT/'runtime/transfer_analysis-summary.log').open('w') as f:q=subprocess.run(['Rscript','--vanilla','inputs/payload/summary.R'],cwd=work,env=e,stdout=f,stderr=subprocess.STDOUT)
TIMES.append(dict(step='R12b/summary.R',seconds=time.monotonic()-start,exit_code=q.returncode));write(ROOT/'verification/final_times.tsv',TIMES)
if q.returncode:raise RuntimeError('R12b summary failed')
(r12b/'results').mkdir(exist_ok=True);shutil.copytree(work/'outputs',r12b/'results/summary',dirs_exist_ok=True)
for r in checks:
 if r['task']=='R12b' and ('/summaries/' in r['path'] or '/transfer_inputs/' in r['path']):(ROOT/_public_analysis_directory(r['path'])).unlink(missing_ok=True)
run(r12b/'scripts/assemble_tables_tables.py','R12b')
# This marker means producer execution only. Comparison is a separate strict gate.
(ROOT/'verification/FINAL_PRODUCERS_EXECUTED').write_text('EXECUTED\n')

from compare_tables import compare
results=[]
for r in checks:
 if not SKIP_R12A and r['task']=='R12a' and ('/synthesis/' in r['path'] or '/normalization/' in r['path'] or '/figure_values/' in r['path'] or Path(r['path']).name in ['analysis_families.tsv','matched_baseline_dataset.tsv','matched_baseline_pooled.tsv','headline_metrics.tsv','headline_results.tsv']):
  z=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']));z['path']=r['path'];results.append(z)
 if r['task']=='R12b' and ('/summaries/' in r['path'] or '/transfer_inputs/' in r['path']):
  z=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']),ignore_execution=True);z['path']=r['path'];results.append(z)
 if not SKIP_R13 and r['task']=='R13' and Path(r['path']).name in ['item1_case_mix_operating_points.tsv','item2_case_mix_admission_and_discordance.tsv','item3_reassay.tsv','item4_labels.tsv','item5_haemolysis.tsv','item6_assay_assay_transfer.tsv']:
  z=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']),ignore_execution=True);z['path']=r['path'];results.append(z)
for expected in (M/'published_rules_analysis/out/managed').glob('*.tsv'):
 z=compare(r14work/'outputs'/_public_analysis_directory(expected.name),expected,ignore_execution=True);z['path']='R14/managed/'+expected.name;results.append(z)
write(ROOT/'verification/final_comparison.tsv',results)
if not all(r['pass_'] for r in results):raise RuntimeError('Final producer table comparison failed; see verification/final_comparison.tsv')
(ROOT/'verification/FINAL_SUBSET_VALID').write_text('VALID\n')
