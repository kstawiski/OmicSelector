#!/usr/bin/env python3
"""Acquire deposited source files with a required expected SHA-256.

A URL guessed from an accession is never silently accepted: the downloaded
bytes must match the study receipt. Unresolved entries stop the requested run.
"""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import argparse,csv,hashlib,os,urllib.request
PAPER=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
p=argparse.ArgumentParser();p.add_argument('--accession',action='append');p.add_argument('--verify-only',action='store_true');a=p.parse_args()
root=Path(os.environ.get('OMICSELECTOR_PAPER_ROOT',str(PAPER/'work'))).resolve();rows=list(csv.DictReader((PAPER/'data/raw_downloads.tsv').open(),delimiter='\t'))
accessions=set(a.accession or [])
if 'GSE94533_PUBLIC' in accessions:accessions.add('GSE94533')
selected=[r for r in rows if not a.accession or r['accession'] in accessions]
if not selected:raise SystemExit('Accession is absent from the acquisition manifest')
for r in selected:
 if not r['url']:raise SystemExit('Unresolved acquisition dependency: '+r['accession']+' '+r['path'])
 rel=Path(r['path'])
 if rel.is_absolute() or '..' in rel.parts:raise SystemExit('Invalid destination in manifest')
 dest=root/_public_analysis_directory(rel);dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.is_file() and sha(dest)==r['expected_sha256']:continue
 if a.verify_only:raise SystemExit('Missing or changed raw file: '+r['path'])
 tmp=dest.with_name(dest.name+'.partial')
 try:
  with urllib.request.urlopen(r['url'],timeout=120) as response,tmp.open('wb') as f:
   while True:
    b=response.read(1024*1024)
    if not b:break
    f.write(b)
  if sha(tmp)!=r['expected_sha256']:raise RuntimeError('Downloaded bytes differ from the frozen source: '+r['url'])
  tmp.replace(dest)
  evidence=root/'verification/raw_downloads.tsv';evidence.parent.mkdir(parents=True,exist_ok=True)
  new=not evidence.exists()
  with evidence.open('a') as e:
   if new:e.write('accession\tpath\turl\tsha256\tbytes\n')
   e.write('\t'.join([r['accession'],r['path'],r['url'],sha(dest),str(dest.stat().st_size)])+'\n')
 finally:tmp.unlink(missing_ok=True)
 print(r['accession'],r['expected_sha256'],r['path'])
