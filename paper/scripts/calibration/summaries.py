"""Accepted final calibration, practicality and thread-control summary producers."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import os,sys,csv,subprocess,shutil,tarfile,time
from run import ROOT,PAPER,M,model_environment,rows,write
from compare_tables import compare
base=ROOT/_public_analysis_directory(M);timings=[];checks=[]
def run(command,cwd,env=None):
 e=model_environment();e.update(env or {})
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:e[k]='1'
 start=time.monotonic();log=ROOT/'runtime'/(_public_analysis_directory('late-'+str(len(timings))+'.log'))
 with log.open('w') as f:q=subprocess.run(command,cwd=cwd,env=e,stdout=f,stderr=subprocess.STDOUT)
 timings.append(dict(command=' '.join(command),seconds=time.monotonic()-start,exit_code=q.returncode));write(ROOT/'verification/late_times.tsv',timings)
 if q.returncode:raise RuntimeError(str(log))
expected=rows(PAPER/'expected/final_output_sha256.tsv')
tasks=sys.argv[1:] or ['R15b','PB1','TM1']
tasks=[{'calibration':'R15b','practical_cost':'PB1','thread_control':'TM1'}.get(task,task) for task in tasks]
for task in tasks:
 exp=[r for r in expected if r['task']==task and r['status']=='accepted_completed_source']
 if not exp:raise RuntimeError('Completed producer not integrated: '+task)
 work=base/_public_analysis_directory(task)
 if task=='R15b':
  wd=work/'replay';(wd/'inputs').mkdir(parents=True,exist_ok=True);shutil.copytree(work/'payload',wd/'inputs/payload',dirs_exist_ok=True)
  with tarfile.open(work/'payload/candidate.tar.gz') as t:t.extractall(wd/'outputs',filter='data')
  if (wd/'outputs/package').exists():shutil.rmtree(wd/'outputs/package')
  (wd/'outputs/OmicSelector').rename(wd/'outputs/package')
  run(['Rscript','--vanilla','inputs/payload/test_firth.R'],wd,{'R15_HOST':'public-replay'})
  for r in exp:(wd/'outputs'/_public_analysis_directory(Path(r['path']).name)).unlink(missing_ok=True)
  run(['Rscript','--vanilla','inputs/payload/summarize.R'],wd,{'R15_HOST':'public-replay'})
  for r in exp:
   src=wd/'outputs'/_public_analysis_directory(Path(r['path']).name);dst=ROOT/_public_analysis_directory(r['path']);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
 else:
  for r in exp:(ROOT/_public_analysis_directory(r['path'])).unlink(missing_ok=True)
  run(([sys.executable,str(work/'scripts/aggregate.py')] if task=='PB1' else ['Rscript','--vanilla',str(work/'scripts/analyze.R')]),ROOT)
 for r in exp:
  c=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']),ignore_execution=True);c['path']=r['path'];checks.append(c)
 write(ROOT/'verification/late_comparison.tsv',checks)
if not checks or not all(r['pass_'] for r in checks):raise RuntimeError('Final late-producer output mismatch')
(ROOT/'verification/LATE_SUMMARY_VALID').write_text('VALID\n')
