"""Verify the immutable E1 environment before constructing either runtime variant."""
from pathlib import Path
from importlib.metadata import version
import csv,subprocess,json
p=Path(__file__).parent
for line in (p/'E1-python-freeze.txt').read_text().splitlines():
 if '==' in line:
  name,want=line.split('==',1)
  assert version(name)==want,(name,version(name),want)
q=subprocess.check_output(['Rscript','--vanilla','-e','x<-installed.packages();write.table(x[,c("Package","Version")],stdout(),sep="\\t",row.names=FALSE,quote=FALSE)'],text=True)
got={r['Package']:r['Version'] for r in csv.DictReader(q.splitlines(),delimiter='\t')}
for r in csv.DictReader((p/'E1-r-packages.tsv').open(),delimiter='\t'):
 assert got[r['Package']]==r['Version'],(r['Package'],got.get(r['Package']),r['Version'])
print(json.dumps({'python_packages':'matched','R_packages':'matched','scope':'E1 dependency base before installing the recorded paper package'}))
