"""Read named project evidence and write immutable assembly generations."""
from pathlib import Path
import csv, hashlib, io, json, os, re, stat, zipfile
csv.field_size_limit(32*1024*1024)

PAPER = Path('.')
BASE = PAPER / 'analyses'
PKG = BASE / 'analysis_8_md_change_attribution/out/sources/8fe4f684'
OUT = Path(__file__).resolve().parents[1]
SEED = OUT / 'inherited'
V6 = PAPER / 'submission_nargab_v6/upload/supplementary_data'
INPUTS = {}
ISSUES = []
CHECKS = []

READ_ROOTS=[BASE,PAPER/'analysis',PAPER/'manuscript',PAPER/'submission_nargab_v6',PAPER/'code',PAPER/'software_release_267',Path('./package-source'),PKG]

def authorized_read(path):
    p=Path(path).resolve()
    if not any(p==root or root in p.parents for root in READ_ROOTS):
        raise PermissionError('Input path outside the declared analysis sources: '+str(p))
    return p

def sha(data): return hashlib.sha256(data).hexdigest()

def fingerprint(path):
    p=authorized_read(path); h=hashlib.sha256(); size=0
    with p.open('rb') as f:
        while block:=f.read(1024*1024):
            h.update(block);size+=len(block)
    record={'path':str(p),'bytes':size,'sha256':h.hexdigest()}
    INPUTS[str(p)]=record
    return record

def raw(path):
    p = authorized_read(path)
    data = p.read_bytes()
    INPUTS[str(p)] = {'path': str(p), 'bytes': len(data), 'sha256': sha(data)}
    return data

def read(path): return raw(path).decode('utf-8-sig')
def js(path): return json.loads(read(path))
def rows(path, sep=None):
    p = Path(path)
    return list(csv.DictReader(io.StringIO(read(p)), delimiter=sep or (',' if p.suffix=='.csv' else '\t')))

def zread(path, member):
    data = raw(path)
    with zipfile.ZipFile(io.BytesIO(data)) as z: result = z.read(member)
    INPUTS[str(path)+'::'+member] = {'path': str(path)+'::'+member, 'bytes': len(result), 'sha256': sha(result)}
    return result.decode('utf-8-sig')

def zrows(path, member): return list(csv.DictReader(io.StringIO(zread(path, member))))

def write_tsv(dest, data, fields=None):
    data = list(data)
    if fields is None:
        fields = list(dict.fromkeys(k for row in data for k in row))
    with Path(dest).open('x', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n', extrasaction='raise')
        w.writeheader()
        for row in data:
            w.writerow({k: json.dumps(v,ensure_ascii=False,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def write_text(dest, value):
    with Path(dest).open('x') as f: f.write(value)

def annotated(path, sep=None):
    return [dict(r, source_path=str(path), source_row=i) for i,r in enumerate(rows(path,sep),2)]

def issue(area, detail, source): ISSUES.append(dict(area=area,limitation=detail,source_path=str(source)))
def check(name, passed, detail=''):
    CHECKS.append(dict(check=name,passed=bool(passed),detail=detail))
    if not passed: raise AssertionError(name+': '+str(detail))

def numeric(value):
    try:
        n=float(value)
        return n if n==n and abs(n)!=float('inf') else None
    except (TypeError,ValueError): return None

def joined(values): return ';'.join(dict.fromkeys(str(v) for v in values if v is not None and str(v)))

def code_files(root, max_depth=8):
    """Bounded source scan. Never follow directory symlinks or runtime/data trees."""
    root=Path(root)
    skip={'revisions','builds','logs','jobs','units','library','libraries','Rlib','Rlibs','r_library','site-library','runtime','runtime_cache','rootfs','venv','tmp','cache','model_cache','models','checkpoints','inputs','data','raw','downloads','archives','__pycache__','.git','.records','.mplconfig','node_modules','chunks','cells','cell_status','partial_results','full_results','predictions','outputs','staging'}
    if not root.is_dir(): return []
    result=[]
    for parent,dirs,files in os.walk(root,followlinks=False):
        depth=len(Path(parent).relative_to(root).parts)
        dirs[:]=[d for d in dirs if d not in skip and not d.endswith(('.Rcheck','.Rproj.user')) and not Path(parent,d).is_symlink() and depth<max_depth]
        for name in files:
            p=Path(parent,name)
            if p.suffix.lower() in {'.r','.py','.sh'} and not p.is_symlink():result.append(p)
    return sorted(result)

def paragraphs(path):
    text=read(path)
    section=''
    for m in re.finditer(r'\S[^\n]*(?:\n(?!\s*\n)[^\n]*)*',text):
        block=m.group(0).rstrip()
        if block.startswith('#'): section=block.splitlines()[0]
        yield text[:m.start()].count('\n')+1,section,block
