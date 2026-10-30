#!/usr/bin/env python3
"""Fetch the pinned public mirror and retain human mature records.
The original miRBase endpoint was unavailable. The representative matrix was
verified exactly against the deposited preparation with this mirrored reference.
"""
from pathlib import Path
import os,json,base64,hashlib,urllib.request,re,sys
PAPER=Path(__file__).resolve().parents[1];ROOT=Path(os.environ.get('OMICSELECTOR_PAPER_ROOT',str(PAPER/'work'))).resolve()
p=ROOT/'data/reference/mirbase_v22/mature.fa'
if p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()=='18cc850341914106f0a4f6d121a2102df94700bd8ad2f408774fe00ecd1b0345':
 print('Pinned human reference cache verified:',p);sys.exit(0)
url='https://api.github.com/repos/nf-core/test-datasets/git/blobs/f600493177fa12431a6d3e8cc6f72d733a2b7544'
with urllib.request.urlopen(url,timeout=120) as f:d=json.load(f)
b=base64.b64decode(d['content']);assert hashlib.sha256(b).hexdigest()=='3c521fc9bea3c7993e71cf188b11caf73e2d40ac83454e32876455317cc6d342'
h=b''.join(x for x in re.split(b'(?=>)',b) if x.startswith(b'>hsa-'));assert hashlib.sha256(h).hexdigest()=='18cc850341914106f0a4f6d121a2102df94700bd8ad2f408774fe00ecd1b0345'
p=ROOT/'data/reference/mirbase_v22/mature.fa';p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(h)
print('Pinned human reference verified:',p)
