"""Pinned original nf-core quantification and frozen specimen selection; verify preparation after execution."""
from pathlib import Path as _PublicPath
import sys as _public_sys
_public_root=next(p for p in _PublicPath(__file__).resolve().parents if (p/'config/project.json').is_file())
_public_sys.path.insert(0,str(_public_root/'scripts'))
from run import analysis_directory as _public_analysis_directory
from pathlib import Path
import argparse,csv,os,subprocess,sys,json,shutil
PAPER=Path(__file__).resolve().parents[1];ROOT=Path(os.environ.get('OMICSELECTOR_PAPER_ROOT',str(PAPER/'work'))).resolve();code=PAPER/'data/smrnaseq'
p=argparse.ArgumentParser();p.add_argument('accession',choices=['GSE71008','PRJEB30542']);a=p.parse_args();acc=a.accession
nextflow=shutil.which('nextflow')
if not nextflow:raise SystemExit('Install Nextflow and Java, and provide Docker or Apptainer; the pinned pipeline is nf-core/smrnaseq 2.4.1.')
out=ROOT/'processed'/_public_analysis_directory(acc)/'smrnaseq';out.mkdir(parents=True,exist_ok=True);sheet=out/'public_samplesheet.csv'
with (code/(_public_analysis_directory(acc+'_samplesheet.csv'))).open() as f:rs=list(csv.DictReader(f))
for r in rs:
 for k in ['fastq_1','fastq_2']:
  if r[k]:r[k]=str(ROOT/_public_analysis_directory(r[k]));assert Path(r[k]).is_file(),r[k]
with sheet.open('w') as f:w=csv.DictWriter(f,fieldnames=rs[0]);w.writeheader();w.writerows(rs)
for n in ['hsa_mature.fa','hsa_hairpin.fa','hsa.gff3']:assert (ROOT/'reference/mirbase_v22'/_public_analysis_directory(n)).is_file(),n
adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' if acc=='GSE71008' else 'ATCACCGACTGCCCATAGAGAGG'
cmd=[nextflow,'run','nf-core/smrnaseq','-r','2.4.1','-profile',os.environ.get('OMICSELECTOR_NF_PROFILE','docker'),'-c',str(code/'nextflow.config'),'--input',str(sheet),'--outdir',str(out),'--genome','GRCh38','--protocol','custom','--clip_r1','0','--three_prime_clip_r1','0','--three_prime_adapter',adapter,'-work-dir',str(ROOT/'work/smrnaseq'/_public_analysis_directory(acc)),'-resume']
subprocess.run(cmd,cwd=ROOT,check=True)
raw=out/'mirna_quant/mirtop/mirna.tsv';assert raw.is_file(),raw
converted=out/'public_counts.tsv';subprocess.run(['Rscript','--vanilla',str(code/'convert_mirtop_to_omicselector.R'),acc,str(raw),str(converted),str(ROOT/'reference/mirbase_v22/hsa_mature.fa')],cwd=ROOT,check=True)
if acc=='GSE71008':
 dest=ROOT/'data_bucket/GSE71008/expression_matrix_smrnaseq.tsv';dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(converted,dest)
else:
 ids=(code/'PRJEB30542_selected_rows.tsv').read_text().splitlines()[1:]
 with converted.open() as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;data={v['sample_id']:v for v in r}
 assert all(i in data for i in ids)
 with (out/'PRJEB30542_dedup116_counts.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data[i] for i in ids)
(ROOT/'verification').mkdir(exist_ok=True);(ROOT/'verification'/(_public_analysis_directory(acc+'_quantification_command.json'))).write_text(json.dumps(cmd,indent=2)+'\n')
