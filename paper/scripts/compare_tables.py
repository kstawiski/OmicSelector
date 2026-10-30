"""Compare table values with explicit transport and canonical checksums.

Pass requires equal headers and row counts, matching strings and missingness,
and every finite numeric difference <= 1e-10 * max(1, abs(expected)).
Canonical hashes use 12 decimal places and normalize missing tokens; named
execution-identity columns are excluded only when requested. Canonical and raw
hash matches and maximum absolute error are reported independently of pass.
"""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
import csv,hashlib,math,json,os,re
from pathlib import Path
MISSING={'','NA','NaN','nan'}
EXECUTION_COLUMNS={'computing_host','computing_identity','summary_host'}
# Filesystem provenance is normalized only in the audited columns and JSON leaf.
# An absolute prefix ending in /OmicSelector_paper/ denotes the project root,
# independent of the machine or account that produced the table.
PATH_COLUMNS={'source_bundle','source_file','package_path','predictions','artifact','path','bundle','Source path','Producer summary source','source_path'}
PROJECT_PREFIX=re.compile(r'^/(?:[^/]+/)*?OmicSelector_paper/')
def normalize_path(value):
 if not value.startswith(('/', './', '../')):return value
 root=os.environ.get('OMICSELECTOR_PAPER_ROOT','')
 if root and value.startswith(root+'/'):value='./'+value[len(root)+1:]
 value=PROJECT_PREFIX.sub('./',value,count=1)
 for old,new in json.loads(os.environ.get('OMICSELECTOR_PROVENANCE_REWRITES','{}')).items():
  if value.startswith(old):value=new+value[len(old):]
 return value
def normalize_provenance(column,value):
 if column in PATH_COLUMNS:return normalize_path(value)
 if column=='Producer summary record':
  try:record=json.loads(value)
  except (ValueError,TypeError):return value
  if isinstance(record,dict) and isinstance(record.get('source_path'),str):
   path=normalize_path(record['source_path'])
   # Retain the exact spelling and spacing of all other serialized fields.
   updated=re.sub(r'("source_path"\s*:\s*)"(?:\\.|[^"\\])*"',
                 lambda match:match.group(1)+json.dumps(path,ensure_ascii=False),value,count=1)
   parsed=json.loads(updated)
   if parsed==dict(record,source_path=path):return updated
 return value
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compare(actual,expected,ignore_execution=False):
 actual=Path(actual);expected=Path(expected)
 out=dict(path=str(actual),expected_sha256=sha(expected),actual_sha256=sha(actual) if actual.exists() else '',raw_hash_match=False,rows=0,maximum_absolute_error=0.,strings_match=True,missingness_match=True,canonical_hash_match=False,pass_=False)
 if not actual.exists():return out
 out['raw_hash_match']=out['expected_sha256']==out['actual_sha256']
 def read(p):
  with p.open(newline='') as f:return list(csv.reader(f,delimiter=',' if p.suffix=='.csv' else '\t'))
 a=read(actual);b=read(expected)
 if not a or not b or a[0]!=b[0] or len(a)!=len(b):return out
 # Every row of both tables must have exactly the header width.
 if any(len(r)!=len(a[0]) for r in a) or any(len(r)!=len(b[0]) for r in b):return out
 keep=[i for i,k in enumerate(a[0]) if not (ignore_execution and k in EXECUTION_COLUMNS)]
 ca=[a[0]];cb=[b[0]];out['rows']=len(a)-1;numeric=True
 for x,y in zip(a[1:],b[1:]):
  if len(x)!=len(y):return out
  nx=[];ny=[]
  for i in keep:
   u,v=x[i],y[i]
   u=normalize_provenance(a[0][i],u)
   if u in MISSING or v in MISSING:
    out['missingness_match']&=(u in MISSING)==(v in MISSING);nx.append('NA' if u in MISSING else u);ny.append('NA' if v in MISSING else v);continue
   try:fu,fv=float(u),float(v)
   except ValueError:
    out['strings_match']&=u==v;nx.append(u);ny.append(v);continue
   if math.isfinite(fu) and math.isfinite(fv):
    delta=abs(fu-fv);out['maximum_absolute_error']=max(out['maximum_absolute_error'],delta);numeric&=delta<=1e-10*max(1,abs(fv));nx.append(format(0. if abs(fu)<5e-13 else fu,'.12f'));ny.append(format(0. if abs(fv)<5e-13 else fv,'.12f'))
   else:out['strings_match']&=u==v;nx.append(u);ny.append(v)
  ca.append(nx);cb.append(ny)
 def valuehash(v):return hashlib.sha256(json.dumps(v,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
 out['actual_value_sha256']=valuehash(ca);out['expected_value_sha256']=valuehash(cb);out['canonical_hash_match']=out['actual_value_sha256']==out['expected_value_sha256']
 out['pass_']=out['strings_match'] and out['missingness_match'] and numeric
 return out
