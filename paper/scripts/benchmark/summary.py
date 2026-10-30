"""Execute benchmark's completed prepared-input splice and original synthesis producers."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import os,sys,subprocess,shutil,time,json,csv
from run import ROOT,PAPER,M,analysis_environment,rows,write
from compare_tables import compare
contracts={r['task']:r for r in json.loads((PAPER/'config/active_summary_contracts.json').read_text())};assert 'D2' in contracts
base=ROOT/_public_analysis_directory(M)/'benchmark';payload=base/'summary_payload';out=base/'public-summary';shutil.rmtree(out,ignore_errors=True);out.mkdir(parents=True)
e=analysis_environment();e.update(D2_PAYLOAD=str(payload),D2_CELL=str(out));start=time.monotonic()
with (ROOT/'runtime/benchmark-summary.log').open('w') as log:q=subprocess.run(['Rscript','--vanilla',str(payload/'synthesis_entry.R')],cwd=ROOT,env=e,stdout=log,stderr=subprocess.STDOUT)
write(ROOT/'verification/benchmark_times.tsv',[dict(step='corrected_synthesis',seconds=time.monotonic()-start,exit_code=q.returncode)])
if q.returncode:raise RuntimeError('D2 summary producer failed; see runtime/D2-summary.log')
# Assemble rank, transfer and decision summaries from the completed native records.
with (ROOT/'runtime/benchmark-ancillary.log').open('w') as log:
 q=subprocess.run(['Rscript','--vanilla',str(payload/'ancillary_entry.R')],cwd=ROOT,env=e,stdout=log,stderr=subprocess.STDOUT)
if q.returncode:raise RuntimeError('D2 ancillary producer failed; see runtime/D2-ancillary.log')
raw=base/'out/retrieved/aux3/cells/calibration'
assert raw.is_dir(),'D2 native calibration records are missing'
shutil.copytree(raw,out/'calibration',dirs_exist_ok=True)
code=(base/'scripts/deployment_compare.py').read_text().replace("O=B/'out/deployment'","O=B/'public-summary/deployment'")
exec(compile(code,str(base/'scripts/deployment_compare.py'),'exec'),{'__name__':'d2_public_deployment','__file__':str(base/'scripts/deployment_compare.py'), '_public_analysis_directory': _public_analysis_directory, '_public_table': _public_table})
subprocess.run([sys.executable,str(PAPER/'scripts/benchmark/tables.py')],cwd=ROOT,env=e,check=True)
subprocess.run([sys.executable,str(PAPER/'scripts/benchmark_export_evidence.py')],cwd=ROOT,env=e,check=True)
mapping=json.loads((PAPER/'expected/benchmark_provenance_path_mappings.json').read_text());os.environ['OMICSELECTOR_PROVENANCE_REWRITES']=json.dumps(mapping)
checks=[]
for r in rows(PAPER/'expected/final_output_sha256.tsv'):
 if r['task']!='D2':continue
 rel=Path(r['path']).relative_to(M/'benchmark/out');actual=out/_public_analysis_directory(rel)
 z=compare(actual,PAPER/_public_analysis_directory(r['expected_path']),ignore_execution=True);z['path']=str(rel);checks.append(z)
write(ROOT/'verification/benchmark_comparison.tsv',checks);assert checks and all(r['pass_'] for r in checks),'Corrected synthesis differs'
(ROOT/'verification/benchmark_SUMMARY_VALID').write_text('VALID\n')
