from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import json,csv,statistics,math,collections,sys,hashlib
P=Path(__file__).resolve().parents[1];O=P/'out';tasks=json.loads((P/'payload/tasks.json').read_text());roster=json.loads((P/'payload/roster.json').read_text())
def tsv(p):
 return list(csv.DictReader(p.open(),delimiter='\t',escapechar='\\')) if p.exists() else []
def write(name,rows,fields=None):
 if fields is None:fields=list(dict.fromkeys(k for r in rows for k in r))
 with (O/_public_analysis_directory(name)).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def num(x):
 try:r=float(x);return r if math.isfinite(r) else None
 except (ValueError,TypeError):return None
def med(v):
 v=[x for x in v if x is not None];return statistics.median(v) if v else ''
def key(r):return (r['method'],r['unit'],r['device'])
scan={r['method']:r for r in tsv(O/'dependency_source_scan.tsv')}
family_names={'A':'Compositional within-sample baselines','B':'Frozen denoisers','D':'Kit-aware transfer','E':'Anchor/reference denoising transfer','F':'Technology-aware transfer','G':'Biofluid-aware transfer','K':'Rank and balance','L':'Likelihood, transport, kernel and density','M':'Image, invariance and signal descriptors','N':'Neural, domain generalisation and self-supervision','O':'Foundation models and nonlinear compositional','NC':'Negative controls'}
roots=[]
if (P/'measurement_sources.json').exists():
 for x in json.loads((P/'measurement_sources.json').read_text()):
  if x.get('kind')=='retention':continue
  root=P/_public_analysis_directory(x['destination'])
  if (root/'measurement_coverage.json').exists():roots.append(root)
events=collections.defaultdict(list);environments=[];weights={};rpkgs=collections.defaultdict(set);all_attempts=[]
for root in roots:
 root_keys=set()
 env=json.loads((root/'environment.json').read_text());env.setdefault('package_version','2.6.9');env.setdefault('package_source_sha256','5bc289307aaa008cc46c0ee6724325eda415be46bdc59c9108545267d7eb92cd');env['E1_source_image']='sha256:2fa387c63144aed56acdd867e81de36aa192346f183dc36f19d6cddd2acfc030';env['container_engine']=env.get('container_engine','apptainer' if env['mode']=='cpu' or env.get('package_version')=='2.6.10' else 'docker');env['SIF_path']='./runtime/reference-home/scratch/study-e1-20261022/lane-2fa387c6.sif' if env['mode']=='cpu' else ''; env['job_verified']=not (root/'SNAPSHOT.json').exists();infos=list((root/'cells').glob('*/sessionInfo.txt'));env['R_version']=infos[0].read_text().splitlines()[0] if infos else '';environments.append(dict(result_directory=str(root.relative_to(P)),**{k:json.dumps(v,sort_keys=True) if isinstance(v,dict) else v for k,v in env.items()}))
 for w in json.loads((root/'weights.json').read_text()):weights[w['sha256']]=w
 for d in sorted((root/'cells').glob('*')):
  if not d.is_dir():continue
  method,unit,device,rep=d.name.rsplit('__',3) if d.name.count('__')==3 else ('','','','')
  # Unit names may include double underscores; parse from each stored row.
  fit=tsv(d/'fit.tsv');load=tsv(d/'load.tsv');err=tsv(d/'error.tsv')
  base=(fit or err)
  if not base:
   if not (d/'fit_process.json').exists():continue
   prefix,device,rep=d.name.rsplit('__',2);method,unit=prefix.split('__',1)
   base=[dict(method=method,unit=unit,device=device,repeat_id=int(rep[1:]),phase='fit',status='execution_failure',error_code='PROCESS_EXIT_WITHOUT_RESULT',error_message='The process ended without a fit or error record; inspect the retained process log.')]
  r=dict(base[0]);r.update(load[0] if load else {});r['repeat_id']=int(r['repeat_id']);r['result_directory']=str(d.relative_to(P));r['host']=env['hostname'];r['cpu_model']=env['cpu_model'];r['gpu']=env.get('gpu','');r['threads']=4;r['package_version']=env['package_version'];r['package_source_sha256']=env['package_source_sha256'];r['package_commit']=env.get('package_commit','');r['fit_package_version']=r.get('fit_package_version',env['package_version']);r['score_package_version']=r.get('score_package_version',env['package_version']);r['persistence_package_version']=r.get('persistence_package_version',env['package_version']);r['job_verified']=not (root/'SNAPSHOT.json').exists()
  r['backend_execution_settings']='ensemble_batch_size=1;estimators=8;float32;frozen_KV_cache=cpu;one_slice_to_cuda;column_embedding_offload=False' if r['method']=='tab-tabicl' and env.get('tabicl_KV_cache_device')=='cpu' else ('ensemble_batch_size=1;frozen_KV_cache=cuda' if r['method']=='tab-tabicl' and env.get('tabicl_ensemble_batch_size')==1 else '')
  for phase in ['fit','score','persistence']:r[phase+'_package_source_sha256']='5bc289307aaa008cc46c0ee6724325eda415be46bdc59c9108545267d7eb92cd' if r[phase+'_package_version']=='2.6.9' else env['package_source_sha256']
  if err:r['error_code']=err[0].get('code','');r['error_message']=err[0].get('message','')
  refusal=tsv(d/'refusal_evidence.tsv')
  if refusal:
   r['n_train']=refusal[0]['n_train'];r['n_features']=refusal[0].get('n_features',int(refusal[0].get('transformed_dimension','0'))+1);r['n_scored']=0
  fitproc=json.loads((d/'fit_process.json').read_text());loadproc=json.loads((d/'load_process.json').read_text()) if (d/'load_process.json').exists() else {}
  r['fit_session_wall_seconds']=fitproc.get('session_wall_seconds','');r['fit_exit']=fitproc['exit_code'];r['load_exit']=loadproc.get('exit_code','')
  r['fit_peak_rss_bytes']=max(fitproc.get('kernel_maximum_rss_bytes',0),fitproc.get('sampled_aggregate_peak_rss_bytes',0));r['load_peak_rss_bytes']=max(loadproc.get('kernel_maximum_rss_bytes',0),loadproc.get('sampled_aggregate_peak_rss_bytes',0)) if loadproc else ''
  r['load_session_wall_seconds']=loadproc.get('session_wall_seconds','')
  stages=fitproc.get('stage_peak_rss_bytes',{})
  r['fitting_sampled_peak_rss_bytes']=stages.get('fit','');r['deployment_sampled_peak_rss_bytes']=stages.get('deployment','');r['serialization_sampled_peak_rss_bytes']=stages.get('serialization','');r['loading_sampled_peak_rss_bytes']=loadproc.get('stage_peak_rss_bytes',{}).get('load','')
  if (d/'derive_process.json').exists():
   repair=json.loads((d/'derive_process.json').read_text());r['canonical_state_session_seconds']=repair['session_wall_seconds'];r['deployment_sampled_peak_rss_bytes']=repair.get('stage_peak_rss_bytes',{}).get('deployment','');r['serialization_sampled_peak_rss_bytes']=repair.get('stage_peak_rss_bytes',{}).get('serialization','');r['canonical_state_peak_rss_bytes']=repair.get('sampled_aggregate_peak_rss_bytes','')
  manifest=json.loads((d/'model_manifest.json').read_text()) if (d/'model_manifest.json').exists() else []
  r['frozen_model_bytes']=sum(x['bytes'] for x in manifest) if manifest else ''
  if key(r) not in root_keys:events[key(r)]=[];root_keys.add(key(r))
  events[key(r)].append(r)
  all_attempts.append(r)
  for q in tsv(d/'r_packages.tsv'):rpkgs[r['method']].add(q['package']+'=='+q['version'])
def cell_measurement_sources(tasks,jobs):
 # Last explicitly declared replacement owns each exact method/unit/device cell.
 owners={}
 for job in jobs:
  exact={tuple(x) for x in job.get('selected_cells',[])}
  for task in tasks:
   if key(task) in exact or (task['device']==job.get('selected_device','cpu') and task['method'] in job.get('selected_methods',[])):
    owners[key(task)]=job
 return owners
for cell,job in cell_measurement_sources(tasks,json.loads((P/'measurement_sources.json').read_text())).items():
 prefixes=[str(Path(job['destination'])/'cells')]
 events[cell]=[e for e in events[cell] if any(e['result_directory'].startswith(prefix+'/') for prefix in prefixes)]
selected={e['result_directory'] for es in events.values() for e in es}
write('all_attempt_measurements.tsv',[dict(e,selected_for_cost_table=e['result_directory'] in selected) for e in all_attempts],fields=list(dict.fromkeys(k for e in all_attempts for k in e))+['selected_for_cost_table'] if all_attempts else ['method','unit','device','repeat_id','selected_for_cost_table'])
# Explicit registry requirements, supplemented by source references and runtime imports.
modules={'reticulate-tabpfn':['numpy','torch','tabpfn'],'reticulate-tabicl':['numpy','torch','tabicl'],'reticulate-tabdpt':['numpy','torch','tabdpt'],'reticulate-kymatio':['numpy','torch','kymatio'],'reticulate-pyod':['numpy','pyod'],'reticulate-torch-torchvision':['numpy','torch','torchvision'],'reticulate-torch-plus-R':['numpy','torch']}
python_metadata={x['distribution'].lower().replace('_','-'):x['version'] for x in tsv(O/'python_dependency_metadata.tsv')}
deps=[]
for r in roster:
 m=r['method_id'];route=r['dep_route'];py=modules.get(route,[])
 req={'glmnet':'glmnet','selbal-github':'selbal','rvinecopulib':'rvinecopulib','isotree-or-base':'isotree (optional; base fallback)'} .get(route,'reticulate' if r['tier']=='R2' else '')
 relevant=[]
 for w in weights.values():
  path=w['path'].lower()
  if (m=='ai-tabpfn' and 'tabpfn' in path and path.endswith('.ckpt'))or(m=='tab-tabicl' and 'tabicl' in path and '/blobs/' in path)or(m=='tab-tabdpt' and 'tabdpt' in path and '/blobs/' in path):relevant.append(w)
 done=[x for t in tasks if t['method']==m for x in events[key(t)]]
 deps.append(dict(method=m,family=r['family'],family_name=family_names[r['family']],estimand=r['estimand'],role=r['role'],registered_fit_function=r['fit_fn'],registered_score_function=r['score_fn'],registry_notes=r['notes'],registry_route=route,registry_R_extras=req,source_referenced_R_extras=scan.get(m,{}).get('source_referenced_R_extras',''),Python_modules=';'.join(py),Python_environment_pins=';'.join(m+'=='+python_metadata.get(m,'unknown') for m in py),Python_required_at_fit=r['tier']=='R2',Python_required_at_scoring='reticulate' in scan.get(m,{}).get('score_referenced_R_extras','').split(';'),fit_source_R_extras=scan.get(m,{}).get('fit_referenced_R_extras',''),score_source_R_extras=scan.get(m,{}).get('score_referenced_R_extras',''),model_weights=';'.join(w['path'] for w in relevant),model_weight_sha256=';'.join(w['sha256'] for w in relevant),model_weight_bytes=sum(w['bytes'] for w in relevant) if relevant else ('' if m in ['ai-tabpfn','tab-tabicl','tab-tabdpt'] else 0),GPU='optional' if m in ['ai-tabpfn','tab-tabicl','tab-tabdpt'] else ('not required; CPU path' if r['tier']=='R2' else 'not used'),GPU_CPU_verification='CPU fit/scoring observed' if any(x.get('device')=='cpu' and x.get('fit_exit')==0 and x.get('consistency_pass')=='TRUE' for x in done) else 'not yet verified',runtime_successful_fits=sum(x.get('fit_exit')==0 for x in done),runtime_successful_reloads=sum(x.get('load_exit')==0 for x in done),observed_noncore_R_packages=';'.join(sorted(rpkgs[m])),serialization_dependency='cloudpickle for Python-held fitted objects; RDS otherwise',GPU_cost_note='TabICL uses 11 GiB RTX 2080 Ti with ensemble batch size 1, eight estimators and full precision; frozen KV tensors reside in host RAM and all eight ensemble slices are copied to CUDA one at a time for specimen scoring. Column-embedding offloading is disabled. Minimum GPU memory was not determined.' if m=='tab-tabicl' else '',dependency_note='Source references may include optional branches; observed packages include benchmark/loading tools. GPU requirement verified by successful CPU runs where available.'))
adapter=P.parent/'gpu_memory_measurements/host_cache/payload/gpu_memory.py'
for d in deps:d['GPU_measurement_adapter']=str(adapter.relative_to(P.parent))+';SHA256='+hashlib.file_digest(adapter.open('rb'),'sha256').hexdigest() if d['method']=='tab-tabicl' else ''
depmap={d['method']:d for d in deps};write('dependencies.tsv',deps)
metrics=['fit_wall_seconds','registered_fit_seconds','fit_session_wall_seconds','score_median_seconds','score_p95_seconds','model_load_seconds','load_session_wall_seconds','fit_peak_rss_bytes','load_peak_rss_bytes','frozen_model_bytes','fitting_sampled_peak_rss_bytes','deployment_sampled_peak_rss_bytes','serialization_sampled_peak_rss_bytes','loading_sampled_peak_rss_bytes']
rows=[]
for t in tasks:
 es=events[key(t)];d=depmap[t['method']];r=dict(t,family_name=family_names[t['family']],seed=101,fold=1 if t['estimand']=='within' else 'held-out cohort',repeats_expected=3,repeats_observed=len(es),status='pending')
 successful=[e for e in es if e.get('fit_exit')==0 and e.get('load_exit')==0 and e.get('consistency_pass')=='TRUE' and e.get('reload_pass')=='TRUE' and e.get('state_match')=='TRUE']
 if len(es)==3:
  if len(successful)==3:r['status']='complete'
  elif all(e.get('error_code') in ['MISSING_METADATA','INSUFFICIENT_ALR_PIVOTS','MISSING_HEMOLYSIS_ANCHORS','UNSEEN_METADATA_LEVEL','EMPTY_FITTED_SUPPORT','NO_USABLE_SUPPORT','MCD_DIMENSION_EXCEEDS_TRAINING_N','SELBAL_GBM_INSUFFICIENT_POSITIVE_VALUES'] for e in es):r['status']='documented_refusal'
  else:r['status']='failed_or_inconsistent'
 elif es:r['status']='incomplete'
 
 for field in ['package_version','package_source_sha256','package_commit','fit_package_version','score_package_version','persistence_package_version','fit_package_source_sha256','score_package_source_sha256','persistence_package_source_sha256','backend_execution_settings']:r[field]=';'.join(sorted({str(e.get(field,'')) for e in es}))
 r['model_load_protocol']='R package startup excluded; Python initialization and restoration included' if t['tier']=='R2' else 'R package startup excluded; frozen R object restoration';r['avoidable_Python_load_setup']=t['tier']=='R2' and not d['Python_required_at_scoring'];r['GPU_measurement_adapter']=d['GPU_measurement_adapter'] if t['device']=='cuda' else '';r['successful_repeats']=len(successful);r['job_verified']=bool(es) and all(e['job_verified'] for e in es)
 # Median measurements retain successful fit/score phases even if reload fails;
 # row status prevents using those as a complete deployment cost.
 for metric in metrics:r[metric]=med([num(e.get(metric)) for e in es])
 r['maximum_difference']=max([num(e.get('maximum_difference')) for e in es if num(e.get('maximum_difference')) is not None],default='')
 r['reload_maximum_difference']=max([num(e.get('reload_maximum_difference')) for e in es if num(e.get('reload_maximum_difference')) is not None],default='')
 r['score_tolerance']=1e-5;r['threads']=4;r['host']=';'.join(sorted({e['host'] for e in es}));r['cpu_model']=';'.join(sorted({e['cpu_model'] for e in es}));r['gpu_device']=';'.join(sorted({e['gpu'] for e in es}));r['n_train']=es[0].get('n_train','') if es else '';r['n_features']=es[0].get('n_features','') if es else '';r['n_scored_per_repeat']=20 if any(str(e.get('n_scored',''))=='20' for e in es) else 0
 for col in ['registry_R_extras','source_referenced_R_extras','Python_modules','model_weights','model_weight_bytes','model_weight_sha256','GPU']:r[col]=d[col]
 r['error_codes']=';'.join(sorted({e.get('error_code','') for e in es if e.get('error_code')}));r['error_messages']=' | '.join(sorted({e.get('error_message','') for e in es if e.get('error_message')}));r['evidence']=';'.join(e['result_directory'] for e in es);rows.append(r)
write('method_costs.tsv',rows);write('repeat_measurements.tsv',[e for es in events.values() for e in es],fields=list(dict.fromkeys(k for es in events.values() for e in es for k in e)) or ['method','unit','device','repeat_id'])
summary=[]
for fam,estimand,dev in sorted({(r['family'],r['estimand'],r['device']) for r in rows}):
 group=[r for r in rows if (r['family'],r['estimand'],r['device'])==(fam,estimand,dev)];valid=[r for r in group if r['status']=='complete'];s=dict(family=fam,family_name=family_names[fam],estimand=estimand,device=dev,expected_cells=len(group),complete_cells=len(valid),all_jobs_verified=bool(group) and all(r['job_verified'] for r in group),refused_cells=sum(r['status']=='documented_refusal' for r in group))
 for field in ['host','cpu_model','gpu_device','package_version','fit_package_version','score_package_version','persistence_package_version']:s[field]=';'.join(sorted({str(r[field]) for r in group if r[field]}))
 s['hardware_comparison_note']='Measurements span different hosts; family GPU timings cannot isolate method cost from hardware.' if len({r['host'] for r in valid})>1 else ''
 for metric in metrics:
  vals=[num(r[metric]) for r in valid if num(r[metric]) is not None]
  s[metric+'_median']=med(vals);s[metric+'_min']=min(vals,default='');s[metric+'_max']=max(vals,default='')
 summary.append(s)
write('family_summary.tsv',summary);write('environment.tsv',environments,fields=list(dict.fromkeys(k for x in environments for k in x)) or ['result_directory','hostname','cpu_model','threads','gpu','versions'])
counts=collections.Counter(r['status'] for r in rows);complete=[r for r in rows if r['status']=='complete'];ref=[r for r in complete if r['method']=='ws-rclr-trimmed']
# A secondary paired comparison changes hardware and ensemble scheduling together.
# It is not the per-device harness comparison and does not change its 1e-5 gate.
secondary_CPU_CUDA=[]
for task in tasks:
 if task['method']!='tab-tabicl' or task['device']!='cuda':continue
 for current in events[key(task)]:
  reference=next((e for e in events[(task['method'],task['unit'],'cpu')] if e['repeat_id']==current['repeat_id']),None)
  if reference is None:continue
  a=tsv(P/_public_analysis_directory(current['result_directory'])/'scores.tsv');b=tsv(P/_public_analysis_directory(reference['result_directory'])/'scores.tsv')
  if len(a)!=20 or len(b)!=20:continue
  assert [x['sample_id'] for x in a]==[x['sample_id'] for x in b]
  secondary_CPU_CUDA.append(max(abs(float(x['score'])-float(y['score'])) for x,y in zip(a,b)))
secondary_maximum=max(secondary_CPU_CUDA,default=None)
report=[]
if ref:
 report.append('Trimmed rCLR: '+ '; '.join(f"{r['unit']}: fit {r['fit_wall_seconds']:.4g} s, specimen median/p95 {r['score_median_seconds']*1000:.4g}/{r['score_p95_seconds']*1000:.4g} ms, load {r['model_load_seconds']:.4g} s, frozen model {r['frozen_model_bytes']/1024:.4g} KiB" for r in ref)+'.')
ss=[s for s in summary if s['device']=='cpu' and s['complete_cells']+s['refused_cells']==s['expected_cells'] and s['score_median_seconds_median']!='']
if ss:
 fastest=min(ss,key=lambda x:x['score_median_seconds_median']);slowest=max(ss,key=lambda x:x['score_median_seconds_median'])
 report.append(f"Among fully observed CPU families, specimen scoring is fastest for {fastest['family_name']} ({fastest['estimand']}; family median {fastest['score_median_seconds_median']:.4g} s) and slowest for {slowest['family_name']} ({slowest['estimand']}; {slowest['score_median_seconds_median']:.4g} s). Family medians and ranges for fit, scoring, load, memory and model size are in family_summary.tsv.")
py=[r['method_id'] for r in roster if r['tier']=='R2']
report.append('Python is required at fit for '+', '.join(py)+'. TabPFN, TabICL, TabDPT, inv-scatter and img-gasfcnn also use Python at scoring. No method requires a GPU; the three foundation classifiers have additional GPU measurements.')
report.append(f"Coverage: {counts['complete']} complete cells and {counts['documented_refusal']} documented refusals of 182, covering 73 methods. Each complete cell has three repetitions on one node and 20 separate specimen calls after warm-up. CPU timings use argos9 with four threads. GPU models and hosts are recorded per row; timings across GPU hosts are not direct hardware comparisons.")
scored=[r for r in rows if r['maximum_difference']!=''];violations=[r for r in scored if float(r['maximum_difference'])>1e-5]
state_warnings=sorted({e['method'] for es in events.values() for e in es if e.get('state_match')=='FALSE' and e.get('reload_pass')=='TRUE'})
if scored:report.append(f"Maximum harness/deployment difference: {max(float(r['maximum_difference']) for r in scored):.6g} (tolerance 1e-5). Transfer methods and CUDA TabPFN/TabICL use 2.6.10. CPU foundation fit measurements use 2.6.9; their scoring and canonical reload checks use 2.6.10. Versions are recorded per measurement. CUDA TabICL uses eight full-precision ensemble members in batches of one; its frozen KV cache resides in host RAM.")
if secondary_maximum is not None:report.append(f'Supplementary TabICL CPU/CUDA comparisons differed by up to {secondary_maximum:.6g}; hardware and ensemble batching changed together.')
report.append('Fit includes deployment construction. Loads use fresh R sessions and include Python setup for R2 methods. Costs depend on these inputs, hardware and serialization protocol. These computation checks establish neither clinical validity nor manuscript readiness.')
if any(r['status'] not in ['complete','documented_refusal'] or not r['job_verified'] for r in rows):report.append('Provisional: computation or managed verification remains pending; see STATUS.md.')
text='\n\n'.join(report)+'\n';assert len(text.split())<=350,len(text.split());(O/'REPORT.md').write_text(text)
ready=all(r['job_verified'] for r in rows) and all(r['status'] in ['complete','documented_refusal'] for r in rows)
(O/'aggregation_status.json').write_text(json.dumps(dict(counts=counts,ready=ready,secondary_TabICL_CPU_CUDA_maximum_difference=secondary_maximum,secondary_TabICL_CPU_CUDA_comparisons=len(secondary_CPU_CUDA),consistency_failure_methods=sorted({r['method'] for r in violations}),reload_fingerprint_warning_methods=state_warnings),indent=2))
print(json.dumps(dict(counts=counts,ready=ready)))
