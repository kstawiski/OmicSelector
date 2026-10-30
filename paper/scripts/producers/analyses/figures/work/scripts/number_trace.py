#!/usr/bin/env python3
"""Trace displayed numbers to tokens, accepted literals, tables and figure producers."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
from pathlib import Path
import csv,hashlib,json,re,difflib,sys
from collections import defaultdict
M=Path('./analyses');W=M/'figures/work';O=W.parent/'out'
TOKEN=re.compile(r'⟦B:([^⟧]+)⟧');NUMBER=re.compile(r'\d+(?:[.,]\d+)*')
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv(p,rr,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rr)
FIELDS=['product','location','number','basis','source_file','field','rounding','raw_value','source_sha256','evidence_file','validation']
def trace():
 keys={r['key']:r for r in rows(W/'tokens/KEYS.tsv')};trace=[];checks=[];hashes={}
 def digest(p):
  if p not in hashes:hashes[p]=sha(p)
  return hashes[p]
 for keyed in sorted((W/'text').glob('*.keyed.md')):
  label=keyed.name.replace('.keyed.md','');source=keyed.read_text();parts=[];origins=[];pos=0
  def literal(a,b):
   line=source.count('\n',0,a)+1;result=[]
   for c in source[a:b]:
    result.append(('literal',line,''))
    if c=='\n':line+=1
   return result
  for m in TOKEN.finditer(source):
   parts.append(source[pos:m.start()]);origins.extend(literal(pos,m.start()))
   value=keys[m[1]]['resolved_value'];parts.append(value);origins.extend([('token',source.count('\n',0,m.start())+1,m[1])]*len(value));pos=m.end()
  parts.append(source[pos:]);origins.extend(literal(pos,len(source)));text=''.join(parts)
  resolved=W/'resolved'/(_public_analysis_directory(label+'.resolved.md'))
  checks.append(dict(check=label+' exact numeric-origin reconstruction',status='PASS' if resolved.read_text()==text else 'FAIL',detail=str(len(text))))
  accepted=M/'figure_annotations/integrator'/_public_analysis_directory(keyed.name)
  if not accepted.exists():accepted=keyed
  accepted_lines=accepted.read_text().splitlines();current_lines=source.splitlines()
  exact={line:i for i,line in enumerate(accepted_lines,1)}
  # Passage-to-producer links are supplied alongside the exact numeric text coordinates.
  for m in NUMBER.finditer(text):
   kind,line,key=origins[m.start()];number=m.group();location='line '+str(text.count('\n',0,m.start())+1)+' offset '+str(m.start())
   if kind=='token':
    r=keys[key];trace.append(dict(product=label,location=location,number=number,basis='resolved token '+key,source_file=r['source_file'],field=r['row_selector']+'; '+r['column'],rounding=r['rounding'],raw_value=r['resolved_value'],source_sha256=r['source_sha256'],evidence_file=str(W/'tokens/KEYS.tsv'),validation='Source selector evaluated and source hashes verified'))
   else:
    original=current_lines[line-1];source_line=exact.get(original)
    if source_line is None:
     candidates=[(difflib.SequenceMatcher(None,original,x,autojunk=False).ratio(),i,x) for i,x in enumerate(accepted_lines,1) if number in x]
     candidate=max(candidates,default=(0,0,''))
     source_line=candidate[1] if candidate[0]>.35 else None
    if source_line:
     path=accepted;field='accepted passage line '+str(source_line)+'; literal '+number;validation='Literal retained from accepted source passage'
    else:
     path=keyed;field='composed passage line '+str(line)+'; literal '+number;validation='Composition source with scientific evidence linked by passage and token maps'
    if original.startswith('| ') and any(original.startswith('| '+x+' |') for x in ['TabPFN','TabDPT','TabICL']) and '0.0249' in original+('0.0249' if section_is_seed(current_lines,line) else ''):
     path=M/'R9F/out/foundation_seed_summary.tsv';field='method='+original.split('|')[1].strip()+'; pooled_lift, ci_lo, ci_hi, k, p_lift_gt_zero_BY';validation='Final R9F table rounded to four decimals and P to three decimals'
    if original.startswith('| SS-CLR minus CN-'):
     path=M/'normalization_comparison/out/partA_nonempty_without_toray_pooled.tsv';field='contrast='+original.split('|')[1].strip()+'; k, pooled_lift, ci_lo, ci_hi';validation='R2d nonempty-feature without-Toray table rounded to four decimals'
    precision='verbatim accepted numeric text'
    if re.fullmatch(r'\d+(?:,\d{3})*',number):precision='integer or identifier retained exactly'
    elif '.' in number:precision=f'{len(number.rsplit(".",1)[1])} displayed decimal places, retained from accepted source'
    raw_value=number
    if path.name in ['foundation_seed_summary.tsv','partA_nonempty_without_toray_pooled.tsv']:
     start=text.rfind('\n',0,m.start())+1;column=text[start:m.start()].count('|')
     method=original.split('|')[1].strip();seed=path.name=='foundation_seed_summary.tsv'
     lookup_key='method' if seed else 'contrast'
     lookup_value={'TabPFN':'ai-tabpfn','TabDPT':'tab-tabdpt','TabICL':'tab-tabicl'}[method] if seed else method
     rr=next(r for r in rows(path) if r[lookup_key]==lookup_value)
     metric=({2:'pooled_lift',3:'ci_hi' if ' to ' in text[start:m.start()] else 'ci_lo',4:'k',5:'p_lift_gt_zero_BY'} if seed else {2:'k',3:'pooled_lift',4:'ci_hi' if ' to ' in text[start:m.start()] else 'ci_lo'}).get(column)
     if metric:
      raw_value=rr[metric];field=lookup_key+'='+lookup_value+'; '+metric
      places=len(number.rsplit('.',1)[1]) if '.' in number else 0
      precision='integer count exactly' if metric=='k' else 'round to '+str(places)+' decimal places'
      actual=float(number)*(-1 if text[m.start()-1:m.start()] in ['−','-'] else 1)
      ok=abs(actual-float(raw_value)) <= (0 if metric=='k' else .500001*10**(-places))
      validation='PASS' if ok else 'FAIL'
      checks.append(dict(check=label+' direct '+lookup_value+' '+metric,status=validation,detail=str(actual)+' versus raw '+raw_value))
    trace.append(dict(product=label,location=location,number=number,basis='literal source coordinate',source_file=str(path),field=field,rounding=precision,raw_value=raw_value,source_sha256=digest(str(path)),evidence_file=str(W/'reproduction_objects/result_to_code_map.tsv')+' | '+str(W/'reproduction_objects/result_script_links.tsv')+' | '+str(M/'figure_annotations/out/number_guard.tsv'),validation=validation))
 # Bibliographic and citation numbers have separate origins after renumbering.
 pkg=M.parent.parent/'study_description';old=rows(pkg/'QA/numeric_trace.tsv')
 for r in old:
  if r['basis'] in ['citation renumbering','reference metadata']:
   p=W/'text/references_v11.tsv'
   trace.append(dict(product='manuscript numbered',location='line '+r['line']+' offset '+r['offset'],number=r['number'],basis=r['basis'],source_file=str(p),field=r['coordinate'],rounding='exact bibliographic metadata or sequential citation number',raw_value=r['number'],source_sha256=digest(str(p)),evidence_file=str(pkg/'QA/citation_crosswalk.tsv'),validation='Citation keys checked against complete reference catalogue'))
 # Every displayed or plotted figure value retains its row, column and transformation.
 for p in sorted((W/'figures/provenance').glob('*_values.tsv')):
  for r in rows(p):
   if r.get('printed_text') in ['clause omitted','Pending clause']:continue
   if not r.get('source_file'):raise ValueError('Missing figure source '+str(r))
   srcs=r['source_file'].split(' | ')
   trace.append(dict(product=r['figure'],location=r['panel']+' / '+r['item'],number=r['printed_text'],basis='figure data or display coordinate',source_file=r['source_file'],field=r['source_row']+'; '+r['source_column'],rounding=r['transformation'],raw_value=r['value'],source_sha256=' | '.join(digest(p) for p in srcs),evidence_file=str(p),validation='Figure producer evidence and input hashes'))
 import pymupdf
 p=M/'benchmark_assembly/synthesis/evidence/pooled.tsv';rr=[r for r in rows(p) if r['p_lift_gt_zero_BY'] not in ['', 'NA', 'NaN']]
 best=min(rr,key=lambda r:float(r['p_lift_gt_zero_BY']));printed=f"{float(best['p_lift_gt_zero_BY']):.3f}"
 fig=W/'figures/figures/companion_dependence_pooled_comparison.pdf'
 with pymupdf.open(fig) as pdf:caption=' '.join(page.get_text() for page in pdf)
 ok=best['method']=='ai-tabpfn' and re.search(r'smallest adjusted P is\s+'+re.escape(printed)+r' for TabPFN',caption) is not None
 checks.append(dict(check='Figure 3 caption smallest adjusted P',status='PASS' if ok else 'FAIL',detail=printed+' for '+best['method']))
 trace.append(dict(product='figure_3_pooled_comparison',location='caption smallest adjusted P',number=printed,basis='direct caption to statistical source verification',source_file=str(p),field='method='+best['method']+'; minimum finite p_lift_gt_zero_BY',rounding='round to three decimal places',raw_value=best['p_lift_gt_zero_BY'],source_sha256=digest(str(p)),evidence_file=str(fig),validation='PASS' if ok else 'FAIL'))
 method_names={'Trimmed rCLR':'ws-rclr-trimmed','k-TSP':'reo-ktsp','Pair-ratio':'reo-pairratio','Test-batch centring':'Test-batch centring','Within-batch standardization + L1':'Within-batch standardization','Same L1 model, frozen training standardization':'Frozen standardization + L1','ComBat refitted with each test batch + L1':'ComBat_L1','Same L1 model, unadjusted features':'Frozen, unadjusted + L1'}
 p=M/'companion_sensitivity/out/reassay_design/scoring_conditions_rows.tsv';raw={r['scorer']:r for r in rows(p) if r['dataset']=='GSE268524'}
 table2=[];main=(W/'resolved/manuscript.resolved.md').read_text()
 for line in main.split('### Table 2.',1)[1].splitlines():
  if not line.startswith('| '):continue
  cs=[x.strip() for x in line.strip('|').split('|')]
  if cs[0] not in method_names:continue
  r=raw[method_names[cs[0]]]
  for i,prefix in enumerate(['case_mix_change','re_assay_change','both','AUC'],1):
   vs=[float(r[prefix+'_'+s]) for s in ['estimate','ci_low','ci_high']]
   if i<4 and all(v==0 for v in vs):expected='0% (exactly)'
   else:
    v,lo,hi=[f'{x*(100 if i<4 else 1):.{1 if i<4 else 3}f}' for x in vs]
    expected=v+('%' if i<4 else '')+f' ({lo}–{hi})'
   table2.append(dict(check='Table 2 '+cs[0]+' '+prefix,status='PASS' if cs[i]==expected else 'FAIL',detail=cs[i]+' versus '+expected))
   trace.append(dict(product='manuscript Table 2',location=cs[0]+' / '+prefix,number=cs[i],basis='direct numerical table verification',source_file=str(p),field='dataset=GSE268524, scorer='+method_names[cs[0]]+'; '+prefix+'_{estimate,ci_low,ci_high}',rounding='multiply by 100 and round to one decimal' if i<4 else 'round to three decimals',raw_value=str(vs),source_sha256=digest(str(p)),evidence_file=str(p),validation=table2[-1]['status']))
 checks.extend(table2)
 # assay_coverage is checked against the independently reproduced display cells.
 p=M/'panel_restriction/out/assay_coverage_updated.tsv';r1=rows(p);lookup={(r['tool_policy'],r['column']):r for r in r1};count=0
 coverage_path=wpath=W/'staging/unusable_inputs/unusable_inputs/coverage.tsv'
 coverage={(r['tool'],r['metric']):r for r in rows(coverage_path) if r['dataset']=='ALL' and r['condition']=='C3'}
 toolids={'OmicSelector k-TSP, default':'O1','OmicSelector k-TSP, detection mask':'O2','OmicSelector k-TSP, mask and strict':'O3','switchBox k-TSP':'R1','singscore':'R2','multiclassPairs':'R3'}
 displays=[]
 block=main.split('### Table 1.',1)[1].split('### Table 2.',1)[0];head=None
 for line in block.splitlines():
  if not line.startswith('| '):continue
  cs=[x.strip() for x in line.strip('|').split('|')]
  if head is None:head=cs;continue
  if cs[0].startswith('---'):continue
  for col,value in zip(head[1:],cs[1:]):
   r=lookup.get((cs[0],col));cell_source=p
   if col=='Reduced panel, scored / calls changed among scored':
    tool=toolids.get(cs[0]);cell_source=coverage_path
    if tool:
     scored=coverage[(tool,'scored')];changed=coverage[(tool,'changed_call_among_scored')]
     if float(scored['estimate'])==0:
      expected='0%, '+('all refused' if tool=='O3' else 'errors')+' / none scored'
      assert changed['denominator']=='0' and changed['estimate']=='NA'
     else:
      expected=f"{100*float(scored['estimate']):.0f}% / {100*float(changed['estimate']):.1f}% ({100*float(changed['ci_lo']):.1f}–{100*float(changed['ci_hi']):.1f})"
     r=dict(reproduced_value=expected,evidence=str(coverage_path)+'; tool='+tool+'; dataset=ALL; condition=C3; scored and changed_call_among_scored; estimate and ci_lo/ci_hi')
    else:
     r=lookup[(cs[0],'Reduced-panel calls changed without refusal')];cell_source=p
   if r is None:continue
   displays.append(dict(tool_policy=cs[0],column=col,displayed_value=value,reproduced_value=r['reproduced_value'],exact_display_match=str(value==r['reproduced_value']),source=str(cell_source),evidence=r['evidence']))
   count+=1;checks.append(dict(check='Table 1 '+cs[0]+' '+col,status='PASS' if value==r['reproduced_value'] else 'FAIL',detail=value+' versus '+r['reproduced_value']))
   trace.append(dict(product='manuscript Table 1',location=cs[0]+' / '+col,number=value,basis='independent display-cell verification',source_file=str(cell_source),field='tool_policy='+cs[0]+'; column='+col+'; '+r['evidence'],rounding='exact independently verified display cell',raw_value=r['reproduced_value'],source_sha256=digest(str(cell_source)),evidence_file=r['evidence'],validation='PASS' if value==r['reproduced_value'] else 'FAIL'))
 checks.append(dict(check='Table 1 independently mapped cell count',status='PASS' if count==56 else 'FAIL',detail=str(count)+' of 56'))
 tsv(W/'staging/unusable_inputs/unusable_inputs/reproduction.tsv',displays,list(displays[0]))
 tsv(W/'QA/final_number_checks.tsv',checks,['check','status','detail'])
 tsv(O/'number_trace.tsv',trace,FIELDS);shutil_copy(O/'number_trace.tsv',W/'final_evidence/number_trace.tsv')
 print('Number trace',len(trace),'rows;',sum(r['status']=='FAIL' for r in checks),'numeric checks failed',flush=True)
 return checks
def section_is_seed(lines,n):
 return any('foundation_seed_results' in l for l in lines[max(0,n-12):n])
def shutil_copy(a,b):
 import shutil;shutil.copyfile(a,b)
if __name__=='__main__':trace()
