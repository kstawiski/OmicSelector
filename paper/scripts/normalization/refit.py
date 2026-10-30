#!/usr/bin/env python3
"""Public GEO-to-prediction route for the registered GSE117063 normalization cell."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import argparse,os,shutil,subprocess,sys
from run import ROOT,PAPER,M,initialize,analysis_environment
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,default=PAPER.parent/'bundle');p.add_argument('--initialized',action='store_true');a=p.parse_args()
if not a.initialized:initialize(a.bundle.resolve())
e=analysis_environment();out=ROOT/'refit/normalization';unit=out/'units/GSE117063';unit.mkdir(parents=True,exist_ok=True)
subprocess.run([sys.executable,str(PAPER/'data/acquire.py'),'--accession','GSE117063'],env=e,check=True)
subprocess.run([sys.executable,str(PAPER/'data/fetch_mirbase.py')],env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'data/prepare_all.R'),str(PAPER),'GSE117063'],cwd=ROOT,env=e,check=True)
shutil.copytree(PAPER/'scripts/normalization',out/'scripts',dirs_exist_ok=True)
for name in ['datasets.tsv','partA_reproduction_gate.tsv','exclusions.tsv','expected_pred.tsv.gz']:shutil.copy2(PAPER/'expected/normalization'/_public_analysis_directory(name),out/_public_analysis_directory(name))
shutil.copy2(PAPER/'expected/normalization/gate.tsv',unit/'gate.tsv')
# Remove fitted artifacts so deposited models cannot satisfy a refit.
for f in unit.iterdir():
 if f.name!='gate.tsv':f.unlink()
e['OUT']=str(out)
subprocess.run(['Rscript','--vanilla',str(out/'scripts/partA.R'),'GSE117063'],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'scripts/normalization/check_refit.R')],cwd=ROOT,env=e,check=True)
