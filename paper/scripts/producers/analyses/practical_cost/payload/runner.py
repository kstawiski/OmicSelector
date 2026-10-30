from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import os,sys,time,json,csv,subprocess,hashlib,shutil,platform,importlib.metadata,traceback
from pathlib import Path
import psutil
P=Path('inputs/payload');O=Path('outputs');O.mkdir(exist_ok=True)
lib=O/'library';lib.mkdir(exist_ok=True)
with (O/'package-install.log').open('w') as log:
 subprocess.run(['R','CMD','INSTALL','--no-multiarch','--no-test-load','--library='+str(lib),str(O/'package')],stdout=log,stderr=subprocess.STDOUT,check=True)
mode=os.environ.get('PB1_MODE','cpu');smoke=os.environ.get('PB1_SMOKE','0')=='1'
tasks=[x for x in json.loads((P/'tasks.json').read_text()) if x['device']==mode or (mode=='gpu' and x['device']=='cuda')]
if os.environ.get('PB1_METHODS'):
 tasks=[x for x in tasks if x['method'] in os.environ['PB1_METHODS'].split(',')]
if mode=='gpu':tasks.sort(key=lambda x:(x['method']!='tab-tabdpt',x['method'],x['unit']))
if smoke:
 tasks=[x for x in tasks if (x['method'],x['unit']) in [('ws-rclr-trimmed','GSE94533_PUBLIC'),('ai-tabpfn','GSE64591'),('tab-tabdpt','GSE85589'),('tab-tabicl','GSE64591'),('nc-ecod-copod','GSE64591'),('kit-stratified-rclr','breast_cancer__GSE210329')]]
versions={}
for n in ['numpy','scipy','scikit-learn','torch','torchvision','tabpfn','tabicl','tabdpt','pyod','kymatio','cloudpickle','psutil']:
 try:versions[n]=importlib.metadata.version(n)
 except importlib.metadata.PackageNotFoundError:versions[n]='NOT_INSTALLED'
env=dict(hostname=platform.node(),cpu_model=next((l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),'unknown'),kernel=platform.release(),python=sys.version,versions=versions,mode=mode,threads=4,affinity=os.sched_getaffinity(0).__str__(),environment={k:os.environ.get(k,'') for k in ['ATEN_CPU_CAPABILITY','OPENBLAS_CORETYPE','MKL_ENABLE_INSTRUCTIONS','ONEDNN_MAX_CPU_ISA','DNNL_MAX_CPU_ISA','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','CUDA_VISIBLE_DEVICES','HF_HUB_OFFLINE','TRANSFORMERS_OFFLINE']})
try:env['gpu']=subprocess.check_output(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total','--format=csv,noheader'],text=True).strip()
except Exception:env['gpu']='none exposed'
(O/'environment.json').write_text(json.dumps(env,indent=2))
weights=[]
for base in ['outputs/runtime/tabpfn','outputs/runtime/huggingface']:
 for f in Path(base).rglob('*'):
  if f.is_file() and not f.is_symlink():weights.append(dict(path=str(f.relative_to(O/'runtime')),bytes=f.stat().st_size,sha256=hashlib.file_digest(f.open('rb'),'sha256').hexdigest()))
(O/'weights.json').write_text(json.dumps(weights,indent=2))
(O/'process_logs').mkdir(exist_ok=True)
def runphase(task,rep,phase,d):
 log=O/'process_logs'/_public_analysis_directory(f'{d.name}__{phase}.log');mem=d/_public_analysis_directory(f'{phase}_time.txt')
 cmd=['Rscript','--vanilla',str(P/'cell.R'),phase,task['method'],task['unit'],str(rep),task['device'],str(d)]
 if Path('/usr/bin/time').exists():cmd=['/usr/bin/time','-v','-o',str(mem)]+cmd
 start=time.monotonic();peak=0;stage_peaks={}
 with log.open('w') as f:
  process=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT)
  root=psutil.Process(process.pid)
  while process.poll() is None:
   try:
    rss=0
    for child in [root]+root.children(recursive=True):
     try:rss+=child.memory_info().rss
     except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    peak=max(peak,rss)
    stage_file=d/"measurement_phase.txt"
    stage=stage_file.read_text().strip() if stage_file.exists() else "startup"
    stage_peaks[stage]=max(stage_peaks.get(stage,0),rss)
   except (psutil.NoSuchProcess,psutil.AccessDenied):pass
   time.sleep(.05)
 result=dict(exit_code=process.returncode,session_wall_seconds=time.monotonic()-start,sampled_aggregate_peak_rss_bytes=peak,sampling_interval_seconds=.05,stage_peak_rss_bytes=stage_peaks)
 if mem.exists():
  for l in mem.read_text().splitlines():
   if 'Maximum resident set size (kbytes):' in l:result['kernel_maximum_rss_bytes']=int(l.split(':')[1])*1024
 (d/_public_analysis_directory(f'{phase}_process.json')).write_text(json.dumps(result,indent=2));return result
rows=[]
for task in tasks:
 for rep in ([1] if smoke else [1,2,3]):
  ident=f"{task['method']}__{task['unit']}__{task['device']}__r{rep}"
  d=O/'cells'/_public_analysis_directory(ident);d.mkdir(parents=True,exist_ok=True)
  print(time.strftime('%FT%TZ',time.gmtime()),'START',ident,flush=True)
  fit=runphase(task,rep,'fit',d);load=None
  if fit['exit_code']==0:
   model=d/'model';files=list(model.rglob('*'));manifest=[]
   for f in files:
    if f.is_file():manifest.append(dict(file=str(f.relative_to(model)),bytes=f.stat().st_size,sha256=hashlib.file_digest(f.open('rb'),'sha256').hexdigest()))
   (d/'model_manifest.json').write_text(json.dumps(manifest,indent=2))
   load=runphase(task,rep,'load',d)
   # Model files are measurement intermediates; preserve hashes, sizes and
   # fresh-session validation, then release their space before the next cell.
   if load['exit_code']==0:shutil.rmtree(model)
  row=dict(**task,repeat=rep,fit_exit=fit['exit_code'],load_exit=None if load is None else load['exit_code'],directory=str(d));rows.append(row)
  (O/'progress.json').write_text(json.dumps(dict(expected=len(tasks)*(1 if smoke else 3),completed=len(rows),rows=rows),indent=2))
  print(time.strftime('%FT%TZ',time.gmtime()),'END',ident,row['fit_exit'],row['load_exit'],flush=True)
# This marker validates complete, explicitly labelled observations, including
# documented method refusals. Scientific consistency is evaluated separately.
assert len(rows)==len(tasks)*(1 if smoke else 3)
assert len({r['directory'] for r in rows})==len(rows)
for r in rows:
 d=Path(r['directory']);assert (d/'fit_process.json').exists()
 if r['fit_exit']==0:assert (d/'fit.tsv').exists() and (d/'model_manifest.json').exists() and (d/'load_process.json').exists()
 if r['fit_exit']!=0 or r['load_exit']!=0:assert (d/'error.tsv').exists() or (O/'process_logs'/_public_analysis_directory(f"{d.name}__fit.log")).exists()
(O/'measurement_coverage.json').write_text(json.dumps(dict(expected=len(rows),completed=len(rows),fit_success=sum(r['fit_exit']==0 for r in rows),load_success=sum(r['load_exit']==0 for r in rows),smoke=smoke,interpretation='Execution and evidence completeness only; not a claim that all methods deploy successfully.'),indent=2))
(O/'VALID').write_text('VALID\n')
