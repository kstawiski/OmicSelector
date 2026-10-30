"""Independent validation of observed timing, score, persistence and aggregation evidence."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,json,math,statistics
P=Path(__file__).resolve().parents[1];O=P/'out'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
rows=read(O/'method_costs.tsv');reps=read(O/'repeat_measurements.tsv');expected=json.loads((P/'payload/tasks.json').read_text());checks=[]
key=lambda r:(r['method'],r['unit'],r['device'])
assert len(rows)==182 and {key(r) for r in rows}=={key(r) for r in expected}
assert len({r['method'] for r in rows})==73
assert len({(key(r),r['repeat_id']) for r in reps})==len(reps)
checks.append(dict(check='registered_grid_and_unique_repeats',pass_=True,observations=len(rows)))
for r in reps:
 d=P/_public_analysis_directory(r['result_directory']);scores=d/'scores.tsv'
 if not scores.exists():continue
 sc=read(scores);assert len(sc)==20 and len({q['sample_id'] for q in sc})==20
 times=[float(q['seconds']) for q in sc];assert all(math.isfinite(x) and x>=0 for x in times)
 order=sorted(times);p95=order[18]+.05*(order[19]-order[18])
 # R type 7: h=(n-1)*0.95=18.05 for n=20, zero-based.
 assert math.isclose(statistics.median(times),float(r['score_median_seconds']),rel_tol=1e-11,abs_tol=1e-12)
 assert math.isclose(p95,float(r['score_p95_seconds']),rel_tol=1e-11,abs_tol=1e-12)
 delta=max(abs(float(q['score'])-float(q['benchmark'])) for q in sc)
 stored_delta=max(float(q['absolute_difference']) for q in sc)
 assert math.isclose(stored_delta,float(r['maximum_difference']),rel_tol=1e-12,abs_tol=1e-12)
 # R write.table writes score columns at approximately 15 significant digits.
 # Bound subtraction error from those text representations; do not relax the
 # 1e-5 scientific check, which uses the original full-precision differences.
 for q in sc:
  a=float(q['score']);b=float(q['benchmark']);difference=float(q['absolute_difference'])
  assert abs(abs(a-b)-difference)<=1e-14*max(1,abs(a),abs(b))+1e-12
 if r['fit_exit']=='0':
  manifest=json.loads((d/'model_manifest.json').read_text());assert sum(q['bytes'] for q in manifest)==float(r['frozen_model_bytes'])
 checks.append(dict(check='repeat_'+d.name,pass_=True,observations=20))
for r in rows:
 es=[e for e in reps if key(e)==key(r)]
 assert len(es)==int(r['repeats_observed'])
 for field in ['fit_wall_seconds','score_median_seconds','score_p95_seconds','model_load_seconds','frozen_model_bytes']:
  vals=[float(e[field]) for e in es if e.get(field,'') not in ['','NA']]
  if vals:assert math.isclose(float(r[field]),statistics.median(vals),rel_tol=1e-12,abs_tol=1e-12)
 if r['status']=='documented_refusal':
  assert len(es)==3 and r['n_scored_per_repeat']=='0'
  for e in es:
   assert e['fit_exit']!='0' and e['error_code']
   d=P/_public_analysis_directory(e['result_directory'])
   if e['error_code']=='MCD_DIMENSION_EXCEEDS_TRAINING_N':
    q=read(d/'refusal_evidence.tsv')[0]
    assert int(q['n_train'])<=int(q['transformed_dimension']) and q['robustbase_available']=='TRUE' and q['require_robust']=='TRUE'
    assert 'robust MCD requires n_train > transformed dimension' in e['error_message']
   if e['error_code']=='SELBAL_GBM_INSUFFICIENT_POSITIVE_VALUES':
    q=read(d/'refusal_evidence.tsv')[0]
    assert int(q['features_with_fewer_than_two_positive_values'])>0 and q['selbal_available']=='TRUE' and q['zero_replacement']=='bayes'
    assert 'GBM method: not enough information to compute t hyper-parameter' in e['error_message']
 if r['status']=='complete':
  assert len(es)==3 and all(e['consistency_pass']=='TRUE' and e['reload_pass']=='TRUE' and e['state_match']=='TRUE' and e['fit_exit']=='0' and e['load_exit']=='0' for e in es)
  assert len({e['host'] for e in es})==1
  assert len({e['gpu'] for e in es})==1
  assert all(float(e['maximum_difference'])<=float(r['score_tolerance']) and float(e['reload_maximum_difference'])<=float(r['score_tolerance']) for e in es)
cpu_hosts={r['host'] for r in rows if r['device']=='cpu' and r['host']};assert len(cpu_hosts)<=1
checks.append(dict(check='aggregation_and_same_node',pass_=True,observations=len(reps)))
assert len((O/'REPORT.md').read_text().split())<=350
(O/'validation.json').write_text(json.dumps(dict(checks=checks,passed=len(checks),observed_repeats=len(reps),pending_cells=sum(r['status']=='pending' for r in rows),all_jobs_verified=all(r['job_verified']=='True' for r in rows)),indent=2))
print('PASS',len(checks),'checks on',len(reps),'observed repeats; pending cells are not validated as completed.')
