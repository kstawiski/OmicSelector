#!/usr/bin/env python3
"""Run public-source refits through accepted case-mix, re-assay and transfer producers."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import argparse,csv,gzip,json,os,shutil,subprocess,sys,time
from run import ROOT,PAPER,M,initialize,analysis_environment,rows,write
p=argparse.ArgumentParser();p.add_argument('family',nargs='+',choices=['case_mix','reassay','transfer_target','published_rule','all']);p.add_argument('--bundle',type=Path,default=PAPER.parent/'bundle');p.add_argument('--initialized',action='store_true');a=p.parse_args()
if not a.initialized:initialize(a.bundle.resolve())
e=analysis_environment();code=PAPER/'scripts/secondary';timings=[]
def run(cmd):
 start=time.monotonic();subprocess.run(list(map(str,cmd)),cwd=ROOT,env=e,check=True);timings.append(dict(family=family,command=' '.join(map(str,cmd)),seconds=time.monotonic()-start));write(ROOT/'verification/secondary_times.tsv',timings)
def acquire(acc):run([sys.executable,PAPER/'data/acquire.py','--accession',acc])
def prepare(acc):run(['Rscript','--vanilla',PAPER/'data/prepare_all.R',PAPER,acc])
for family in (['case_mix','reassay','transfer_target','published_rule'] if 'all' in a.family else a.family):
 out=ROOT/'refit'/_public_analysis_directory(family)
 if out.exists():shutil.rmtree(out)
 out.mkdir(parents=True);e.update(OUT=str(out),REFIT_CODE=str(code))
 if family=='case_mix':
  acc='GSE188627';acquire(acc);run([sys.executable,PAPER/'data/fetch_mirbase.py']);prepare(acc)
  shutil.copytree(PAPER/'scripts/normalization',out/'scripts');shutil.copy2(code/'case_mix.R',out/'scripts/partB.R');shutil.copy2(code/'case_common.R',out/'scripts/common.R')
  (out/'scripts/preparation.R').write_text("prepare_partition<-function(id,train_ids,test_ids){z<-get_cohort(id);list(train=z$expr_per_sample[train_ids,,drop=FALSE],test=z$expr_per_sample[test_ids,,drop=FALSE],fit=list(preparation='fixed zero fill; no sample-estimated statistic'))}\n")
  for n in ['datasets.tsv','partA_reproduction_gate.tsv','exclusions.tsv']:shutil.copy2(PAPER/'expected/normalization'/_public_analysis_directory(n),out/_public_analysis_directory(n))
  (out/'units'/_public_analysis_directory(acc)).mkdir(parents=True);e['REFIT_ACCESSION']=acc
  run(['Rscript','--vanilla',code/'case_mix_prepare.R']);run(['Rscript','--vanilla',out/'scripts/partB.R',acc,'1'])
 elif family=='reassay':
  acquire('GSE268524');d=ROOT/'data_bucket/GSE268524';raw=d/'raw/matrix/GSE268524_series_matrix.txt.gz'
  tab=list(csv.reader(gzip.open(raw,'rt'),delimiter='\t'))
  get=lambda k:next(r[1:] for r in tab if r and r[0]=='!Sample_'+k)
  titles=get('title');gsms=get('geo_accession');disease=next(r[1:] for r in tab if r and r[0]=='!Sample_characteristics_ch1' and r[1].startswith('disease: '));disease=[v.split(': ',1)[1] for v in disease]
  grouped={}
  for t,g,ds in zip(titles,gsms,disease):grouped.setdefault(t.rsplit('_',1)[0],[]).append((t,g,ds))
  meta=[]
  for specimen,vs in sorted(grouped.items()):
   vs.sort(key=lambda v:int(v[0].rsplit('_',1)[1]));assert len({v[2] for v in vs})==1
   meta.append(dict(sample_id=specimen,array_titles=','.join(v[0] for v in vs),gsm_ids=','.join(v[1] for v in vs),sample_type='control' if vs[0][2].lower() in ['healthy','healthy control','normal','non-cancer'] else 'cancer'))
  assert len(meta)==477;write(d/'metadata.tsv',meta)
  for n in ['prepared','logs','units']:(out/_public_analysis_directory(n)).mkdir(exist_ok=True)
  run(['Rscript','--vanilla',code/'reassay_prepare.R']);run(['Rscript','--vanilla',code/'reassay.R','GSE268524__1__1__2'])
 elif family=='transfer_target':
  for acc in ['GSE210329','GSE320031','GSE41922']:acquire(acc)
  run([sys.executable,PAPER/'data/fetch_mirbase.py']);run([sys.executable,PAPER/'data/reconstruct_metadata.py','GSE210329'])
  for acc in ['GSE210329','GSE320031','GSE41922']:prepare(acc)
  inp=out/'inputs';inp.mkdir();e['REFIT_TRANSFER_INPUT']=str(inp)
  shutil.copytree(ROOT/'code',inp/'code');shutil.copy2(code/'transfer_prep.R',inp/'prep.R');shutil.copy2(code/'transfer_preparation.R',inp/'preparation.R');(inp/'q').mkdir()
  run(['Rscript','--vanilla',code/'transfer_prepare.R']);run(['Rscript','--vanilla',code/'transfer.R','reo-metaktsp','101',out])
 else:
  acquire('GSE85589');run([sys.executable,PAPER/'data/fetch_mirbase.py'])
  # The public series matrix supplies the expression data; no deposited derived matrix can take priority.
  for parent in [ROOT/'external/Ext_data/GSE85589',ROOT/'data_bucket/GSE85589']:
   for n in ['expression_matrix.tsv','expression_matrix_rawcounts.tsv','expression_matrix_rawCq.tsv']:(parent/n).unlink(missing_ok=True)
  prepare('GSE85589');run(['Rscript','--vanilla',code/'published_rule.R'])
 run(['Rscript','--vanilla',code/'compare.R',family])
