#!/usr/bin/env python3
"""Public-data refit of the final within-cell calibration producer."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from pathlib import Path
import argparse,subprocess,sys
from run import ROOT,PAPER,initialize,analysis_environment
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,default=PAPER.parent/'bundle');p.add_argument('--initialized',action='store_true');p.add_argument('--accession',default='GSE210329');p.add_argument('--method',default='reo-ktsp');p.add_argument('--seed',default='101');a=p.parse_args()
if not a.initialized:initialize(a.bundle.resolve())
e=analysis_environment()
subprocess.run([sys.executable,str(PAPER/'data/acquire.py'),'--accession',a.accession],env=e,check=True)
if a.accession in ['GSE145176','GSE198834','GSE210329']:subprocess.run([sys.executable,str(PAPER/'data/reconstruct_metadata.py'),a.accession],env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'data/prepare_all.R'),str(PAPER),a.accession],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'scripts/calibration/check_refit.R'),str(PAPER),a.accession,a.method,a.seed],cwd=ROOT,env=e,check=True)
