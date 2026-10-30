from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import csv,gzip,hashlib,json
from pathlib import Path
r=Path('.');b=r/'analysis/singlesample_public_benchmark';o=b/'revision_r4'
new=list(csv.DictReader(gzip.open(o/'tables/GSE204951_prepared_predictions.tsv.gz','rt'),delimiter='\t'))
retained=hashlib.sha256();oldn=newn=0
with gzip.open(b/'within/within_predictions.tsv.gz','rt') as src,gzip.open(o/'tables/within_predictions.tsv.gz','wt',compresslevel=6) as dst:
 header=src.readline();fields=header.rstrip('\n').split('\t');idx=fields.index('cohort');dst.write(header)
 for line in src:
  cols=line.rstrip('\n').split('\t')
  if cols[idx]=='GSE204951':oldn+=1;continue
  dst.write(line);retained.update(line.encode());newn+=1
 writer=csv.DictWriter(dst,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
 for row in new:
  row.update(cache_package_commit='cc483baf93c63159a4cf22039f4269c0aaa33bbb',execution_route='r4_verified_r3_or_pinned_refit',final_package_version='2.6.5',final_package_commit='6455bfe0937a32822ccad003eb8413e47feb2b05',package_version='2.6.5')
  writer.writerow(row)
import pandas as pd
c=pd.read_csv(b/'within/within_cells.tsv',sep='\t');g=pd.read_csv(o/'tables/GSE204951_prepared_cells.tsv',sep='\t');g['execution_route']='r4_verified_r3_or_pinned_refit';g['final_package_version']='2.6.5';g['final_package_commit']='6455bfe0937a32822ccad003eb8413e47feb2b05'
_public_table(pd.concat([c[c.cohort!='GSE204951'],g],ignore_index=True)).to_csv(o/'tables/within_cells.tsv',sep='\t',index=False)
e=pd.read_csv(b/'cohort_eligibility.tsv',sep='\t');e.loc[e.accession=='GSE204951',['n_groups','n_case_groups','n_control_groups']]=[20,10,10];e.loc[e.accession=='GSE204951','reason']='115 profiles from 20 persons; JM04/JM08, JM05/JM09, JM06/JM10 grouped; all aliquots retained'
_public_table(e).to_csv(o/'tables/cohort_eligibility.tsv',sep='\t',index=False)
pass
print('splice complete',oldn,len(new),newn)
