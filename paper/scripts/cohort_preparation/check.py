"""Require the prepared public identity bridge and every registered preparation contract."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
import subprocess,sys
from run import ROOT,PAPER,analysis_environment,rows
e=analysis_environment()
subprocess.run([sys.executable,PAPER/'data/reconstruct_metadata.py','GSE198834'],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',PAPER/'data/reconstruct_sources.R',PAPER,'GSE31568','GSE85589','GSE94533_PUBLIC','TORAY_PUBLIC_CANONICAL'],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',PAPER/'data/prepare_all.R',PAPER],cwd=ROOT,env=e,check=True)
r=rows(ROOT/'verification/public_preparation.tsv')
assert len(r)==34 and all(x['pass']=='TRUE' for x in r)
(ROOT/'verification/ALL_PREPARATIONS_VALID').write_text('VALID\n')
