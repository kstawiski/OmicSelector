from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import csv,gzip,math,json,collections,pathlib,hashlib,socket
B=pathlib.Path('./analyses/robustness');O=B/'out';M=B.parent;ROOT=M.parent.parent
checks=[]
def read(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as f: return list(csv.DictReader(f,delimiter='\t'))
def check(name,passed,detail):
    checks.append(dict(check=name,pass_=bool(passed),detail=str(detail)))
def near(a,b,tol=2e-12): return abs(float(a)-float(b))<=tol
op=read(O/'item1_case_mix_operating_points.tsv');ds={x['dataset'] for x in op if x['experiment']=='five_split_deployment'}
check('12 deployment datasets',len(ds)==12,sorted(ds))
check('operating grid',len(op)==12*11*3*2*6+6*5*3*2,len(op))
for name in ['item1_case_mix_operating_points.tsv','item2_case_mix_admission_and_discordance.tsv','item3_reassay.tsv']:
    r=read(O/_public_analysis_directory(name)); bad=[]
    for x in r:
        if x.get('estimate','') not in ('','NA','NaN'):
            v=float(x['estimate']);lo=x.get('ci_lo','');hi=x.get('ci_hi','')
            if not (-1e-12<=v<=1+1e-12) and x.get('metric')!='paired_difference':bad.append(x)
            if lo not in ('','NA','NaN') and hi not in ('','NA','NaN') and float(lo)>float(hi):bad.append(x)
    check(name+' valid bounds',not bad,len(bad))
# Compare every single-split operating point with the authoritative combined case-mix producer.
orig=read(M/'companion_dependence/out/case_mix/results_all_units.tsv')
mis=[]
for x in op:
    if x['experiment']!='original_case_mix':continue
    ar=x['arm'].replace('_',' ')
    y=next(v for v in orig if v['cohort']==x['dataset'] and v['scorer']==ar and v['case_fraction']==x['fraction'])
    if not near(x['estimate'],y[x['metric']+'_mean']):mis.append((x['dataset'],ar,x['fraction'],x['metric']))
check('all 180 single-split operating points match current combined producer',not mis,mis)
# Independent source calculation of original Toray sensitivity and specificity.
A=ROOT/'analyses/case_mix_duplicates_casemix_dedup/out/predictions/TORAY_PUBLIC_CANONICAL__Test-batch_centring.tsv.gz'
z=read(A)
for fr in ['0.1','0.5','0.9']:
    for cl,mt in [('1','sensitivity'),('0','specificity')]:
        a=[x for x in z if x['case_fraction']==fr and x['outcome']==cl];v=sum((x['decision']=='TRUE')==(cl=='1') for x in a)/len(a)
        ref=next(x for x in op if x['dataset']=='TORAY_PUBLIC_CANONICAL' and x['arm']=='Test-batch_centring' and x['fraction']==fr and x['metric']==mt)
        check('original Toray '+fr+' '+mt,near(v,ref['estimate']),(v,ref['estimate']))
# Independent re-assay same-person difference: aligned 50% replicate 2/3 versus all paired 10/90 appearances.
paths=[ROOT/'analyses/score_identity_identity/out/analysis/reassay_casemix/specimen_decisions.tsv.gz',M/'reassay_design/out/frozen/specimen_decisions.tsv.gz']
paired=collections.defaultdict(dict);ext=collections.defaultdict(lambda:[0,0]);badext=set()
for path in paths:
    with gzip.open(path,'rt') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if x['dataset']!='GSE268524' or x['replicate'] not in ('2','3'):continue
            sc=x['scorer'];r=x['repetition'];id=x['specimen'];rep=x['replicate'];call=x['decision']=='TRUE';valid=x['status']=='ok'
            if x['arm']=='third_replicate_eligible' and x['case_fraction']=='0.5':
                paired[(sc,r,id,x['batch'])][rep]=call if valid else None
            if x['arm']=='joint' and ((rep=='2' and x['case_fraction']=='0.1') or (rep=='3' and x['case_fraction']=='0.9')):
                key=(sc,r,id,rep);ext[key][0]+=call;ext[key][1]+=1
                if not valid:badext.add(key)
only=collections.defaultdict(lambda:[0,0])
for (sc,r,id,b),x in paired.items():
    checkpair=set(x)=={'2','3'}
    if not checkpair:raise ValueError('missing matched key')
    if x['2'] is not None and x['3'] is not None:only[(sc,r,id)][0]+=x['2']!=x['3'];only[(sc,r,id)][1]+=1
person=collections.defaultdict(list)
for (sc,r,id),(n,d) in only.items():
    a=ext[(sc,r,id,'2')];b=ext[(sc,r,id,'3')]
    if d and a[1] and b[1] and (sc,r,id,'2') not in badext and (sc,r,id,'3') not in badext:
        p=a[0]/a[1];q=b[0]/b[1];person[(sc,id)].append(p*(1-q)+(1-p)*q-n/d)
rr=read(O/'item3_reassay.tsv')
for sc in sorted({k[0] for k in person}):
    vals=[sum(v)/len(v) for (s,id),v in person.items() if s==sc];v=sum(vals)/len(vals)
    ref=next(x for x in rr if x['scorer']==sc and x['metric']=='paired_difference' and x['scope']=='person_mean_across_splits')
    check('independent paired re-assay '+sc,len(vals)==416 and near(v,ref['estimate']),(len(vals),v,ref['estimate']))
# Independent any-call numerator after retained-intercept reconstruction.
ix=read(O/'item2_case_mix_admission_and_discordance.tsv');fits=read(M/'companion_sensitivity/out/computing_bundle/partB_fits.tsv')
for ds0 in ['GSE204951','GSE188627','TORAY_PUBLIC_CANONICAL']:
    rp='2';calls=collections.defaultdict(set)
    zero={x['model'] for x in fits if x['dataset']==ds0 and x['repetition']==rp and x['features_selected']=='0'}
    zeroarms=set()
    for k in zero:zeroarms.update({'STD':['FZ-STD','CN-STD'],'QN':['FZ-QN','CN-QN'],'RAW':['FZ-RAW','CN-COMBAT'],'CLR':['SS-CLR']}.get(k,[]))
    for x in read(M/_public_analysis_directory(f'companion_dependence/out/runs/computing_bundle/units/{ds0}/B_{rp}/pred.tsv.gz')):
        calls[(x['arm'],x['specimen'])].add('TRUE' if x['arm'] in zeroarms else x['call'])
    for ar in {k[0] for k in calls}:
        v=[c for (a,id),c in calls.items() if a==ar]
        ref=next(x for x in ix if x['dataset']==ds0 and x['arm']==ar and x['scope']=='split_'+rp and x['metric']=='proportion_unique_specimens_any_changed_call')
        if any('NA' in c or '' in c for c in v):ok=ref['estimate']=='NA'
        else:ok=near(sum(len(c)==2 for c in v)/len(v),ref['estimate'])
        check('independent any-call '+ds0+' '+ar,ok,ref['estimate'])
# Preprocessor assay_transfer.
h=read(O/'item6_assay_assay_transfer.tsv');check('96 handoff calls',len(h)==96,len(h));check('four strict scorers',len({x['scorer'] for x in h})==4,sorted({x['scorer'] for x in h}))
check('all aligned masks',all(x['mask_aligned']=='TRUE' for x in h),'dimension names tested at calls')
check('all frozen-preprocessor extra columns unchanged',all(float(x['score_difference'])==0 for x in h if x['condition']=='extra_background' and x['route']=='frozen_preprocessor'),'all four scorers, both policies')
check('all non-detect inputs refused',all(x['code']=='ALL_NONDETECT' for x in h if x['condition'] in ['all_zero_nondetect','positive_floor_all_nondetect']),'zeros and positive floor, both routes and policies')
check('strict required column removal refused',all(x['code']=='MISSING_REQUIRED_FEATURES' for x in h if x['condition']=='missing_required' and x['strict']=='TRUE'),'four scorers')
# Label refits: complete factorial grid and separation from saved-score label sensitivity.
if (O/'item4_labels.tsv').exists():
    lab=read(O/'item4_labels.tsv');fresh=[x for x in lab if x.get('analysis')=='refitted_five_split_deployment']
    check('complete refit summary grid',len(fresh)==3*6*11*7,len(fresh))
    check('no refit scoring failures',all(x.get('n_failed_scores')=='0' for x in fresh),sorted({x.get('n_failed_scores') for x in fresh}))
    old=[x for x in lab if x.get('analysis')=='fixed_scores_five_split_deployment' and x.get('scenario')=='cached_labels_11']
    dif=[]
    for x in fresh:
        if x['scenario']!='cached_labels_11':continue
        y=next(v for v in old if v['scope']==x['scope'] and v['arm']==x['arm'] and v['metric']==x['metric'])
        dif.append(abs(float(x['estimate'])-float(y['estimate'])))
    check('2.6.9 baseline refits reproduce saved summaries',max(dif)<1e-12,max(dif))
with open(O/'evidence/independent_validation.tsv','w') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=['check','pass_','detail']);w.writeheader();w.writerows(checks)
print(json.dumps({'checks':len(checks),'passed':sum(x['pass_'] for x in checks),'host':socket.gethostname()}))
assert all(x['pass_'] for x in checks)
