from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import csv
from compare_tables import compare
checks=[]
for p in Path('expected').glob('*.tsv'):
 r=compare(Path('outputs')/_public_analysis_directory(p.name),p,ignore_execution=True);r['path']=p.name;checks.append(r)
Path('verification').mkdir(exist_ok=True)
with Path('verification/calibration.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(checks)
if not all(r['pass_'] for r in checks):raise SystemExit('Calibration table comparison failed')
Path('verification/CALIBRATION_VALID').write_text('VALID\n')
