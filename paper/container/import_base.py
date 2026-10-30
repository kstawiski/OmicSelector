#!/usr/bin/env python3
"""Load the exact E1 environment from checksum-verified, size-limited release parts."""
from pathlib import Path
import argparse,csv,hashlib,json,subprocess
p=argparse.ArgumentParser();p.add_argument('--assets',type=Path,required=True);a=p.parse_args()
manifest=a.assets/'container_assets.tsv'
if not manifest.exists():raise SystemExit('Missing container_assets.tsv; obtain the environment release assets.')
rows=list(csv.DictReader(manifest.open(),delimiter='\t'))
assert rows and len({r['image_id'] for r in rows})==1
for r in rows:
 f=a.assets/r['file'];assert f.name==r['file'] and f.stat().st_size==int(r['bytes'])
 with f.open('rb') as src:assert hashlib.file_digest(src,'sha256').hexdigest()==r['sha256'],f.name
proc=subprocess.Popen(['docker','load'],stdin=subprocess.PIPE)
try:
 for r in rows:
  with (a.assets/r['file']).open('rb') as src:
   for b in iter(lambda:src.read(1024*1024),b''):proc.stdin.write(b)
finally:proc.stdin.close()
assert proc.wait()==0
subprocess.run(['docker','tag',rows[0]['image_id'],'omicselector-study-e1-lane:2fa387c6'],check=True)
actual=subprocess.check_output(['docker','image','inspect','--format','{{.Id}}','omicselector-m19-e1-lane:2fa387c6'],text=True).strip()
assert actual==rows[0]['image_id'],(actual,rows[0]['image_id'])
print(json.dumps({'image_id':actual,'parts':len(rows),'verified':True}))
