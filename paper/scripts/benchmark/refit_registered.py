#!/usr/bin/env python3
"""Explicit accession preparation and registered-cell refitting; no readiness claim."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
import argparse,os,subprocess,sys
from run import ROOT,PAPER,initialize,analysis_environment,rows
p=argparse.ArgumentParser();p.add_argument('--accession',action='append',required=True);p.add_argument('--methods',default='');p.add_argument('--seeds',default='101,202,303,404,505');p.add_argument('--bundle',type=__import__('pathlib').Path,default=PAPER.parent/'bundle');p.add_argument('--initialized',action='store_true');a=p.parse_args()
registry={r['analysis_unit']:r for r in rows(PAPER/'data/cohort_registry.tsv')}
unknown=set(a.accession)-set(registry)
if unknown:raise SystemExit('Unregistered analysis units: '+', '.join(sorted(unknown)))
if not a.initialized:initialize(a.bundle.resolve())
e=analysis_environment()
cmd=[sys.executable,str(PAPER/'data/acquire.py')]
for acc in a.accession:cmd+=['--accession',acc]
subprocess.run(cmd,cwd=ROOT,env=e,check=True)
for acc in a.accession:
 if acc in ['GSE145176','GSE198834','GSE210329']:subprocess.run([sys.executable,str(PAPER/'data/reconstruct_metadata.py'),acc],cwd=ROOT,env=e,check=True)
for acc in sorted(set(a.accession)&{'GSE71008','PRJEB30542'}):subprocess.run([sys.executable,str(PAPER/'data/quantify_public.py'),acc],cwd=ROOT,env=e,check=True)
derived=sorted(set(a.accession)&{'GSE31568','GSE85589','GSE94533_PUBLIC','TORAY_PUBLIC_CANONICAL'})
if derived:subprocess.run(['Rscript','--vanilla',str(PAPER/'data/reconstruct_sources.R'),str(PAPER),*derived],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'data/prepare_all.R'),str(PAPER),*a.accession],cwd=ROOT,env=e,check=True)
affected=sorted(set(a.accession)&{'GSE320031','GSE50013','GSE65071','GSE70080','GSE83977','GSE85589'})
if affected:subprocess.run(['Rscript','--vanilla',str(PAPER/'data/training_contracts.R'),str(PAPER),*affected],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',str(PAPER/'scripts/benchmark/refit_registered.R'),str(PAPER),','.join(a.accession),a.methods,a.seeds],cwd=ROOT,env=e,check=True)
