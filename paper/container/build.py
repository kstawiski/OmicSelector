#!/usr/bin/env python3
"""Build CPU/CUDA images from the verified E1 image and pinned source archive."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import argparse,csv,hashlib,json,os,shutil,subprocess,tarfile,tempfile
PAPER=Path(__file__).resolve().parents[1];ROOT=Path(os.environ.get('OMICSELECTOR_PAPER_ROOT',str(PAPER/'work'))).resolve()
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,default=PAPER.parent/'bundle');p.add_argument('--target',choices=['cpu','cuda'],action='append');a=p.parse_args();bundle=a.bundle.resolve()
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
base='omicselector-study-e1-lane:2fa387c6';expected='sha256:9527e5812bcfde9045b9bbd7f0c70ec53c13d898d38fa5fe7dc91493dcbfbecd'
assert subprocess.check_output(['docker','image','inspect','--format','{{.Id}}',base],text=True).strip()==expected
meta=json.loads((PAPER/'container/package-source.json').read_text());scratch=ROOT/'runtime/container-build';scratch.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='build-',dir=scratch) as temp:
 context=Path(temp);shutil.copytree(PAPER/'container',context/'container');dest=context/'container/OmicSelector-source.tar.gz';source=bundle/'tree/software/OmicSelector-source.tar.gz'
 if source.is_file():shutil.copy2(source,dest)
 else:
  for row in csv.DictReader((bundle/'assets.tsv').open(),delimiter='\t'):
   asset=bundle/_public_analysis_directory(row['file']);assert sha(asset)==row['sha256']
   with tarfile.open(asset,'r:gz') as t:
    try:m=t.getmember('software/OmicSelector-source.tar.gz')
    except KeyError:continue
    assert m.isfile()
    with t.extractfile(m) as src,dest.open('wb') as dst:shutil.copyfileobj(src,dst)
   break
 assert dest.is_file() and sha(dest)==meta['archive_sha256']
 for target in a.target or ['cpu','cuda']:
  subprocess.run(['docker','build','--network=none','--build-arg','E1_BASE='+base,'--target',target,'-t','omicselector-paper:'+target,'-f',str(context/'container/Dockerfile'),str(context)],check=True)
