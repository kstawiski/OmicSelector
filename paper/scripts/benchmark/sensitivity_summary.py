"""Rebuild final benchmark_sensitivity_summary synthesis with the approved unchanged statistical producers."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import os,subprocess,sys,time,shutil,csv
from run import ROOT,PAPER,M,environment,rows,write
from compare_tables import compare
base=ROOT/_public_analysis_directory(M)/'benchmark_sensitivity_summary';checks=rows(PAPER/'expected/final_output_sha256.tsv')
def selected(r):
 return r['task']=='R12c' and r['path'].endswith('.tsv') and ('/synthesis/' in r['path'] or '/normalization/' in r['path'] or '/figure_values/' in r['path'] or Path(r['path']).name in ['analysis_families.tsv','matched_baseline_dataset.tsv','matched_baseline_pooled.tsv','headline_metrics.tsv','headline_results.tsv'])
checks=[r for r in checks if selected(r)]
assert checks
for r in checks:(ROOT/_public_analysis_directory(r['path'])).unlink(missing_ok=True)
for n in ['out/logs','out/evidence','out/synthesis/logs','out/normalization','out/figure_values']:(base/_public_analysis_directory(n)).mkdir(parents=True,exist_ok=True)
# The approved producer expects the inherited document layout before it rewrites these three files.
for name in ['REPORT.md','results_summary.md','table_notes.md']:
 shutil.copy2(ROOT/_public_analysis_directory(M)/'benchmark_reference_summary/out/synthesis'/_public_analysis_directory(name),base/'out/synthesis'/_public_analysis_directory(name))
e=environment();e.update(PATH=str(ROOT/'runtime/bin')+os.pathsep+e.get('PATH',''),R1B=str(base),OUT=str(base/'out/synthesis'),OMICSELECTOR_PAPER_ROOT=str(ROOT))
timings=[]
for name in ['repool.R','normalization.R','independent.py','assemble_outputs.py']:
 path=base/'scripts'/_public_analysis_directory(name)
 cmd=[str(ROOT/'runtime/R/bin/exec/R'),'--vanilla','--slave','--file='+str(path)] if path.suffix=='.R' else [sys.executable,str(path)]
 start=time.monotonic()
 with (ROOT/'runtime'/(_public_analysis_directory('R12c-'+name+'.log'))).open('w') as log:q=subprocess.run(cmd,cwd=ROOT,env=e,stdout=log,stderr=subprocess.STDOUT)
 timings.append(dict(step=name,seconds=time.monotonic()-start,exit_code=q.returncode));write(ROOT/'verification/benchmark_sensitivity_summary_times.tsv',timings)
 if q.returncode:
  for f in (base/'out/logs').glob('*.log'):print(f.name, f.read_text(errors='replace')[-20000:])
  raise RuntimeError('R12c '+name+' failed; see runtime logs')
# Serialized reference grouped cells retain their original package location.
# Adapt this provenance column only; all numerical values remain untouched.
path_changes=[]
for r in checks:
 f=ROOT/_public_analysis_directory(r['path'])
 with f.open(newline='') as stream:
  reader=csv.DictReader(stream,delimiter='\t');fields=reader.fieldnames;data=list(reader)
 if 'package_path' not in fields:continue
 changed=False
 for i,row in enumerate(data):
  value=row['package_path'];marker='/software_release_267/'
  if value.startswith('/') and marker in value:
   replacement='./software_release_267/'+value.split(marker,1)[1]
   path_changes.append(dict(path=r['path'],row=i+1,column='package_path',original=value,portable=replacement))
   row['package_path']=replacement;changed=True
 if changed:
  with f.open('w',newline='') as stream:
   writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(data)
if path_changes:write(ROOT/'verification/benchmark_sensitivity_summary_provenance_path_changes.tsv',path_changes)
results=[]
for r in checks:
 z=compare(ROOT/_public_analysis_directory(r['path']),PAPER/_public_analysis_directory(r['expected_path']),ignore_execution=True);z['path']=r['path'];results.append(z)
write(ROOT/'verification/benchmark_sensitivity_summary_comparison.tsv',results)
assert all(r['pass_'] for r in results),'Final synthesis table values differ'
(ROOT/'verification/BENCHMARK_SENSITIVITY_SUMMARY_VALID').write_text('VALID\n')
