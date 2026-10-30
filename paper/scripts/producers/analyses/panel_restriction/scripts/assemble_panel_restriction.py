from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,math,collections,hashlib
W=Path(__file__).resolve().parents[1];D=W.parent;P=D/'single_specimen/reproduction_objects_v8/panel_restriction/supplementary_file_panel_restriction_panel_restriction.csv';A=W/'out/panel_restriction'
def read(p,sep='\t'):
 with p.open() as f:return list(csv.DictReader(f,delimiter=sep))
rows=read(P,',');fields=list(rows[0]);r=read(A/'panel_sensitivity_comparisons.tsv');cells=read(A/'panel_sensitivity_cells.tsv');ids={'GSE320031','GSE50013','GSE65071','GSE70080','GSE83977','GSE85589'}
scenarios={'full_panel20':'Full assay; 20-miRNA score','full_panel10':'Full assay; 10-miRNA score','assay100_panel20':'100-miRNA assay; 20-miRNA score','assay50_panel20':'50-miRNA assay; 20-miRNA score','assay20_panel10':'20-miRNA assay; 10-miRNA score'}
g=collections.defaultdict(list);cg=collections.defaultdict(list)
for t,grp in [(r,g),(cells,cg)]:
 for x in t:grp[(x['method'],x['cohort'],scenarios[x['scenario']])].append(x)
colspec={8:('spearman_score','mean'),9:('spearman_score','min'),10:('spearman_score','max'),11:('mean_abs_percentile_rank_drift','mean'),12:('auc_method','mean'),13:('auc_method','min'),14:('auc_method','max'),15:('auc_delta_vs_full_panel20','mean'),16:('auc_delta_vs_full_panel20','min'),17:('auc_delta_vs_full_panel20','max')}
def finite(v):
 try:return math.isfinite(float(v))
 except:return False
def summarize(gr,col,op):
 v=[float(r[col]) for r in gr if finite(r[col])]
 if not v:return 'NOT_ESTIMABLE'
 return format({'mean':lambda:math.fsum(v)/len(v),'min':lambda:min(v),'max':lambda:max(v)}[op](),'.15g')
# Gate the exact same aggregation against the stored panel_restriction affected rows.
B=D/'input_robustness_analysis/out/part_b/analysis/panel_restriction'
og=collections.defaultdict(list)
for x in read(B/'panel_sensitivity_comparisons.tsv'):og[(x['method'],x['cohort'],scenarios[x['scenario']])].append(x)
gates=[]
for row in rows:
 k=(row[fields[0]],row[fields[1]],row[fields[2]])
 if k[1] not in ids:continue
 for i,(col,op) in colspec.items():
  actual=summarize(og[k],col,op);stored=row[fields[i]]
  delta=abs(float(actual)-float(stored)) if finite(actual) and finite(stored) else 0 if actual==stored else float('inf')
  tol=1e-9 if 'AUC' in fields[i] else 1e-10
  gates.append(dict(analysis='panel_restriction',unit=k[1],method=k[0],scenario=k[2],metric=fields[i],stored=stored,reproduced=actual,max_abs_difference=delta,tolerance=tol,passed=delta<=tol))
with (A/'summary_reproduction_gate.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(gates[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(gates)
assert all(x['passed'] for x in gates),'S13 aggregation gate failed'
new=[]
for row in rows:
 n=dict(row);k=(row[fields[0]],row[fields[1]],row[fields[2]])
 if k[1] in ids:
  cc=cg[k];gr=g[k];assert len(cc)==5
  ne=sum(x['eligible']=='TRUE' for x in cc)
  n[fields[3]]='5';n[fields[4]]=str(ne);n[fields[5]]='Yes' if ne==5 else 'No'
  n[fields[6]]='; '.join(sorted({x['ineligible_reason'] for x in cc if x['eligible']!='TRUE' and x['ineligible_reason'] not in ['','NA']}))
  n[fields[7]]=str(max(int(x['n']) for x in gr)) if gr else 'NOT_ESTIMABLE'
  for i,(col,op) in colspec.items():n[fields[i]]=summarize(gr,col,op)
 new.append(n)
assert len(new)==680 and sum(x[fields[1]] in ids for x in new)==120
with (A/_public_analysis_directory(P.name)).open('w') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(new)
with (W/'out/input_sha256.tsv').open('a') as f:f.write(f'{P}\t{hashlib.sha256(P.read_bytes()).hexdigest()}\t{P.stat().st_size}\n')
print('panel_restriction assembled: 680 rows, 120 affected rows, 560 unchanged rows')
