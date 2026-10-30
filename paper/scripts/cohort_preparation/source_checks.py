from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from pathlib import Path
import os,subprocess,argparse
from run import ROOT,PAPER,initialize,analysis_environment
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path);a=p.parse_args();initialize(a.bundle);e=analysis_environment()
subprocess.run(['Rscript','--vanilla',PAPER/'data/reconstruct_sources.R',PAPER,'GSE31568','GSE85589','GSE94533_PUBLIC','TORAY_PUBLIC_CANONICAL'],cwd=ROOT,env=e,check=True)
subprocess.run(['Rscript','--vanilla',PAPER/'data/prepare_all.R',PAPER,'GSE31568','GSE85589','GSE94533_PUBLIC','TORAY_PUBLIC_CANONICAL'],cwd=ROOT,env=e,check=True)
