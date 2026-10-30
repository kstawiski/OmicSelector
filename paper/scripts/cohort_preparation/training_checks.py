from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from pathlib import Path
import os,subprocess,sys
from run import ROOT,PAPER,initialize,analysis_environment
import argparse
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path);a=p.parse_args();initialize(a.bundle);e=analysis_environment()
ids=['GSE320031','GSE50013','GSE65071','GSE70080','GSE83977','GSE85589']
for id in ids:subprocess.run([sys.executable,PAPER/'data/acquire.py','--accession',id],env=e,cwd=ROOT,check=True)
subprocess.run([sys.executable,PAPER/'data/fetch_mirbase.py'],env=e,cwd=ROOT,check=True)
subprocess.run(['Rscript','--vanilla',PAPER/'data/prepare_all.R',PAPER,*ids],env=e,cwd=ROOT,check=True)
subprocess.run(['Rscript','--vanilla',PAPER/'data/training_contracts.R',PAPER,*ids],env=e,cwd=ROOT,check=True)
