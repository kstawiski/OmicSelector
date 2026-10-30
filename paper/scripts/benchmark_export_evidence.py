from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import sys,shutil
import pandas as pd
from run import ROOT,M as RELM
B=ROOT/_public_analysis_directory(RELM)/'benchmark';M=B.parent;T=ROOT/'analyses';O=B/'public-summary';parts_only=True
def read(p):return pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t',low_memory=False)
def save(d,p):
 p.parent.mkdir(parents=True,exist_ok=True);d=d.copy()
 for c in d.select_dtypes(include=['bool']).columns:d[c]=d[c].map({True:'TRUE',False:'FALSE'})
 _public_table(d).to_csv(p,sep=',' if p.suffix=='.csv' else '\t',index=False,na_rep='NA')
def copy(a,b):
 b.parent.mkdir(parents=True,exist_ok=True)
 if a.is_dir():shutil.copytree(a,b,dirs_exist_ok=True)
 else:shutil.copy2(a,b)
for p in O.glob('matched_baseline_*.tsv'):copy(p,O/'synthesis/evidence/GSE198834'/_public_analysis_directory(p.name))
for pp in list(O.glob('independent*.tsv'))+list(O.glob('optimizer*.tsv')):copy(pp,O/'synthesis'/_public_analysis_directory(pp.name))

for name in ['benchmark_sources.tsv','matched_fold_differences.tsv','fold_row_alignment.tsv']:
 copy(O/_public_analysis_directory(name),O/'synthesis/evidence/GSE198834'/_public_analysis_directory(name))
