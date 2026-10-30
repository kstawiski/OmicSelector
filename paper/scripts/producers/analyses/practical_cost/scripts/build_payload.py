from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import shutil,hashlib,csv,json
P=Path(__file__).resolve().parents[1];M=P.parent;root=M.parent.parent;D=P/'payload'
files={M/'cohort_inventory/out/OmicSelector_2.6.9.tar.gz':'candidate.tar.gz',root/'code/methods/matched_null_benchmark.R':'matched_null_benchmark.R',root/'code/methods/singlesample_roster_benchmark.R':'benchmark.R',M/'transfer/scripts/run_transfer_nsclc_revision_r3.R':'transfer_reference.R'}
for src,dst in files.items():shutil.copy2(src,D/_public_analysis_directory(dst))
roster=list(csv.DictReader((P/'source/OmicSelector/inst/extdata/singlesample_method_manifest.csv').open()))
assert len(roster)==73 and sum(r['estimand']=='within' for r in roster)==50
(D/'roster.json').write_text(json.dumps(roster,indent=2))
tasks=[]
for r in roster:
 for u in (['GSE94533_PUBLIC','GSE64591','GSE85589'] if r['estimand']=='within' else ['breast_cancer__GSE210329']):
  for device in (['cpu','cuda'] if r['method_id'] in ['ai-tabpfn','tab-tabicl','tab-tabdpt'] else ['cpu']):
   tasks.append(dict(method=r['method_id'],unit=u,device=device,family=r['family'],estimand=r['estimand'],tier=r['tier'],dep_route=r['dep_route']))
# Reference then fast R then Python; keep every repeat of a cell adjacent on one node.
tasks.sort(key=lambda x:(x['method']!='ws-rclr-trimmed',x['tier']=='R2',x['estimand'],x['method'],x['unit']))
(D/'tasks.json').write_text(json.dumps(tasks,indent=2))
with (P/'out/expected_cells.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=tasks[0],delimiter='\t');w.writeheader();w.writerows(tasks)
inputs=list(files)+[Path(x) for x in (P/'out/input_paths.txt').read_text().splitlines()]
with (P/'out/input_sha256.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256'])
 for p in inputs:w.writerow([p,p.stat().st_size,hashlib.file_digest(p.open('rb'),'sha256').hexdigest()])
print(len(tasks),'method/unit/device cells;',len(tasks)*3,'repeat runs')
