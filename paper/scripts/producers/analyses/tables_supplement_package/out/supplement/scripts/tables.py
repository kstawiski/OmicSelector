from common import *

from collections import Counter, defaultdict

import statistics, subprocess

def build_transfer(dest):
    p=BASE/'transfer_transfer_final/out/full_results/validated_cells.tsv'; data=annotated(p)
    key=lambda r:(r['method'],r['held_out'],r['seed'])
    check('S14 complete 2160 unique cells',len(data)==len({key(r) for r in data})==2160)
    check('S14 24 methods 18 targets 5 seeds',tuple(len({r[k] for r in data}) for k in ['method','held_out','seed'])==(24,18,5))
    check('S14 outcome class supplied',all(r['outcome_class'] for r in data))
    write_tsv(dest/'transfer_cells.tsv',data)
    groups=defaultdict(list)
    for r in data:groups[(r['method'],r['held_out'])].append(r)
    # Reuse producer summaries if available; never invent uncertainty formulas.
    summaries={}; summary_path=None
    for p2 in [BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance_r3.tsv',BASE/'transfer_transfer_final/out/full_results/merged/transfer_summary.tsv',BASE/'transfer_transfer_final/out/full_results/merged/transfer_performance.tsv']:
        if p2.exists():
            for r in annotated(p2):summaries[(r.get('method'),r.get('held_out'))]=r
            summary_path=p2
    out=[]
    for (method,target),rr in sorted(groups.items()):
        check('five transfer records '+method+'/'+target,len(rr)==5)
        ok=[r for r in rr if r['outcome_class']=='evaluated' and numeric(r['auc_method']) is not None]
        first=rr[0]
        row={'Method identifier':method,'Cancer comparison':first['disease'],'Training analysis units':first['training_cohorts'],
             'Number of training analysis units':first['n_training_cohorts'],'Held-out analysis unit':target,
             'Evaluable repetitions (of five)':len(ok),'All five repetitions evaluable':len(ok)==5,
             'miRNAs in the prespecified disease-wide feature universe':first['n_common_features'],
             'Required mechanism-engagement signal':joined(r['engagement_metric'] for r in rr),
             'Observed mechanism-engagement value':joined(r['engagement_value'] for r in rr),
             'Required mechanism engaged':all(r['engagement_ok']=='TRUE' for r in rr),
             'Mean held-out AUC':statistics.mean(float(r['auc_method']) for r in ok) if len(ok)==5 else '',
             'AUC 95% CI lower limit':'','AUC 95% CI upper limit':'',
             'Mean AUC difference from transfer reference':statistics.mean(float(r['lift']) for r in ok) if len(ok)==5 else '',
             'AUC-difference 95% CI lower limit':'','AUC-difference 95% CI upper limit':'',
             'Mean reference AUC':statistics.mean(float(r['auc_baseline']) for r in ok) if len(ok)==5 else '',
             'Outcome classes':dict(Counter(r['outcome_class'] for r in rr)),'Refusal codes':joined(sorted({r['refusal_code'] for r in rr if r['refusal_code']})),
             'Non-evaluable reasons':joined(r['ineligible_reason'] for r in rr),'Seeds':joined(sorted(r['seed'] for r in rr)),
             'Execution commits':joined(r['package_commit'] for r in rr),'Source path':str(p),'Source rows':joined(r['source_row'] for r in rr),
             'Aggregation':'Arithmetic means only when all five seeds are evaluated; no uncertainty recomputed from seed means.',
             'Producer summary source':'','Producer summary record':''}
        sr=summaries.get((method,target))
        if sr:
            row['Producer summary source']=sr['source_path'];row['Producer summary record']=sr
            row['Aggregation']='Transfer summaries use transfer_performance_r3.tsv. Outcome counts use validated_cells.tsv.'
            # Only exact known field names are eligible for an uncertainty copy.
            for source,destcol in [('auc_method','Mean held-out AUC'),('auc_method_ci_lo','AUC 95% CI lower limit'),('auc_method_ci_hi','AUC 95% CI upper limit'),('lift','Mean AUC difference from transfer reference'),('lift_ci_lo','AUC-difference 95% CI lower limit'),('lift_ci_hi','AUC-difference 95% CI upper limit')]:
                if source in sr:row[destcol]=sr[source]
        out.append(row)
    write_tsv(dest/'transfer_summary.tsv',out)
    check('S14 432 method-target summaries',len(out)==432)
    if not summary_path: issue('S14','The cell evidence supports five-seed arithmetic means. Confidence limits are blank because no named producer summary was found; they were not reconstructed with a new estimator.',p)
    return data,out
