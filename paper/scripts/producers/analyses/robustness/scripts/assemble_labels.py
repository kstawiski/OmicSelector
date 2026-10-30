import csv,json,hashlib,pathlib,subprocess,os,socket,platform,datetime,shutil
B=pathlib.Path('./analyses/robustness');O=B/'out';M=B.parent

def read(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with open(p,'w') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
assert not (O/'DONE').exists(),'Do not modify a sealed completed bundle.'
ref=read(O/'evidence/label_refit_results.tsv');fixed=read(O/'evidence/label_assignment_summary_results.tsv');primary=read(O/'evidence/label_primary_records.tsv')
assert len(ref)==3*6*11*7
for r in ref:r['row_type']='result';r['model_policy']='OmicSelector 2.6.9 refits use the saved train/test allocation. Label-sensitivity assignments apply to training and test specimens. Balanced batches use those assignments.'
for r in fixed:r['row_type']='result'
write(O/'item4_labels.tsv',primary+ref+fixed)
