#!/usr/bin/env python3
"""Rebuild missing public metadata and expression/GEO identity bridges.

GSE198834 labels are joined through public count profiles and GEO sample IDs.
The prepared benchmark checksum becomes active only after benchmark completion.
"""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import argparse,csv,gzip,os,re,json
PAPER=Path(__file__).resolve().parents[1];ROOT=Path(os.environ.get('OMICSELECTOR_PAPER_ROOT',str(PAPER/'work'))).resolve()
p=argparse.ArgumentParser();p.add_argument('accession',choices=['GSE145176','GSE198834','GSE210329']);a=p.parse_args();acc=a.accession;out=ROOT/'data_bucket'/_public_analysis_directory(acc)

def write(name,rows):
 out.mkdir(parents=True,exist_ok=True)
 with (out/_public_analysis_directory(name)).open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
if acc=='GSE210329':
 series=list(csv.reader(gzip.open(out/'raw'/(_public_analysis_directory(acc+'_series_matrix.txt.gz')),'rt'),delimiter='\t'))
 def field(k):return next(r[1:] for r in series if r and r[0]=='!Sample_'+k)
 chars=[r[1:] for r in series if r and r[0]=='!Sample_characteristics_ch1']
 relations=next(r[1:] for r in series if r and r[0]=='!Sample_relation' and r[1].startswith('SRA:'))
 columns={k:field(v) for k,v in [('gsm_id','geo_accession'),('title','title'),('source_name','source_name_ch1'),('organism','organism_ch1'),('platform_id','platform_id'),('instrument_model','instrument_model')]}
 metadata=[]
 for i,gsm in enumerate(columns['gsm_id']):
  state=next(c[i] for c in chars if c[i].startswith('disease state:'))
  assert state in ['disease state: case','disease state: control']
  row={k:v[i] for k,v in columns.items()}
  row.update(characteristics_ch1='; '.join(c[i] for c in chars),sample_type='cancer_prediagnostic' if state=='disease state: case' else 'control',mammography_round=re.search(r'_([0-9]+)$',row['title']).group(1),additional_metadata=relations[i]);metadata.append(row)
 assert len(metadata)==131 and sum(r['sample_type']=='cancer_prediagnostic' for r in metadata)==65
 write('metadata.tsv',metadata)
 raw=ROOT/'external/Ext_data'/_public_analysis_directory(acc)/(_public_analysis_directory(acc+'_miRNA_counts.txt.gz'))
 with gzip.open(raw,'rt') as f:ids=next(csv.reader(f,delimiter='\t'))[1:]
 meta=list(csv.DictReader((out/'metadata.tsv').open(),delimiter='\t'))
 def key(s):return re.sub(r' arr$','',re.sub(r'_[0-9]+$','',s))
 lookup={key(r['title']):r for r in meta};assert len(lookup)==len(meta)
 bridge=[]
 for s in ids:
  r=lookup.get(key(s))
  if r is None:assert s=='TO_1T',s
  else:bridge.append(dict(gsm_id=r['gsm_id'],expr_id=s,title=r['title']))
 assert len(bridge)==130 and len({r['gsm_id'] for r in bridge})==130
 write('sample_id_bridge.tsv',bridge)
else:
 rows=list(csv.reader(gzip.open(out/'raw'/(_public_analysis_directory(acc+'_series_matrix.txt.gz')),'rt'),delimiter='\t'))
 def get(k):return next(r[1:] for r in rows if r and r[0]==k)
 gsm=get('!Sample_geo_accession');titles=get('!Sample_title');chars=[r[1:] for r in rows if r and r[0]=='!Sample_characteristics_ch1']
 meta=[];bridge=[]
 if acc=='GSE145176':
  group=next(r for r in chars if r[0].startswith('group: '))
  for g,t,c in zip(gsm,titles,group):
   assert c in ['group: ALL','group: CONTROL'];meta.append(dict(gsm_id=g,title=t,case_control='case' if c=='group: ALL' else 'control'));bridge.append(dict(gsm_id=g,expr_id=t.split(':')[0]))
 else:
  health=next(r for r in chars if r[0].startswith('health: '));urls=get('!Sample_supplementary_file_1')
  d=list(csv.reader(gzip.open(ROOT/'external/Ext_data'/_public_analysis_directory(acc)/(_public_analysis_directory(acc+'_rawCounts.tsv.gz')),'rt'),delimiter='\t'));chips=d[0];features=[r[0].lower() for r in d[1:]]
  profiles={tuple(int(r[i+1]) for r in d[1:]):chip for i,chip in enumerate(chips)};assert len(profiles)==155
  for g,t,h,url in zip(gsm,titles,health,urls):
   sample=list(csv.reader(gzip.open(out/'raw'/_public_analysis_directory(url.rsplit('/',1)[1]),'rt'),delimiter='\t'));values={r[0].lower():int(r[1]) for r in sample[1:]};profile=tuple(values.get(f,0) for f in features)
   bridge.append(dict(gsm_id=g,chip_id=profiles.get(profile,'UNRESOLVED')))
   assert h in ['health: Healthy','health: Cancer'];meta.append(dict(gsm_id=g,title=t,case_control='case' if h=='health: Cancer' else 'control'))
  missing=[r for r in bridge if r['chip_id']=='UNRESOLVED'];remaining=set(chips)-{r['chip_id'] for r in bridge}
  assert len(missing)==1 and missing[0]['gsm_id']=='GSM5957628' and remaining=={'chip.sample122.tsv'}
  missing[0]['chip_id']=remaining.pop()
  evidence=ROOT/'verification';evidence.mkdir(parents=True,exist_ok=True)
  (evidence/'GSE198834_mapping_exception.json').write_text(json.dumps(dict(exact_profiles=154,unique_remaining_pair=['GSM5957628','chip.sample122.tsv'],note='Public per-sample counts differ from the deposited matrix for this pair; the matrix remains unchanged. Publicly mapped labels differ from deposited labels for 52 of 155 samples.'),indent=2)+'\n')
 write('sample_metadata.tsv',meta);write('sample_id_bridge.tsv',bridge)
print(acc,'metadata and identity bridge reconstructed')
