"""calibration_reference summary replay from saved development and held-out states; no classifier fit."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import os,sys,shutil,subprocess,time,platform
from run import ROOT,PAPER,model_environment,write
from compare_tables import compare
B=ROOT/'analyses/calibration_reference';work=B/'replay_summary';(work/'inputs').mkdir(parents=True,exist_ok=True)
shutil.copytree(B/'summary_payload',work/'inputs/summary_payload',dirs_exist_ok=True)
e=model_environment()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:e[k]='1'
e['R15_HOST']=platform.node();e['TMPDIR']=str(work/'outputs/tmp');(work/'outputs/tmp').mkdir(parents=True,exist_ok=True)
shutil.copytree(ROOT/'runtime/packages/2.6.9/OmicSelector',work/'outputs/package',dirs_exist_ok=True)
start=time.monotonic()
with (ROOT/'runtime/calibration.log').open('w') as f:q=subprocess.run(['Rscript','--vanilla','inputs/summary_payload/summarize.R'],cwd=work,env=e,stdout=f,stderr=subprocess.STDOUT)
write(ROOT/'verification/calibration_time.tsv',[dict(seconds=time.monotonic()-start,exit_code=q.returncode)])
if q.returncode:raise RuntimeError('R15 summary failed; see runtime/calibration.log')
checks=[]
for n in ['calibration','net_benefit','failures','audit','cell_status','score_audit','summary']:
 r=compare(work/'outputs'/_public_analysis_directory(f'{n}.tsv'),PAPER/'expected/calibration_reference'/_public_analysis_directory(f'{n}.tsv'),ignore_execution=True);r['path']=n+'.tsv';checks.append(r)
write(ROOT/'verification/calibration.tsv',checks)
if not all(r['pass_'] for r in checks):raise RuntimeError('R15 table-value verification failed')
(ROOT/'verification/CALIBRATION_VALID').write_text('VALID\n')
