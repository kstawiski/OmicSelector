from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import pandas as pd,numpy as np,hashlib,json,csv
from pathlib import Path
r=Path('.');b=r/'analysis/singlesample_public_benchmark';o=b/'revision_r4';t=o/'tables'
old=pd.read_csv(t/'parsimony_latest_input.tsv',sep='\t');oldf=old[(old.method=='reo-ktsp')&(old.level=='fold')&(old.scope=='main_34')];oldu=old[(old.method=='reo-ktsp')&(old.level=='unit')&(old.scope=='main_34')]
common=oldu.groupby('cohort').eligible.agg(['all','count']);common=common.index[common['all']&(common['count']==10)]
newu=pd.read_csv(t/'parsimony_prepared_units.tsv',sep='\t');newu=newu[(newu.method=='reo-ktsp')&(newu.scope=='main_34')]
vars=['heldout_auc','dense_reference_auc','auc_loss_vs_dense','features_selected_mean'];old34=oldf.groupby(['cohort','pairs_requested'])[vars].mean().reset_index()
rows=[]
for stage,d in [('old_fixed29',oldu[oldu.cohort.isin(common)]),('old_constant_inclusive_fixed34',old34),('corrected_same29',newu[newu.cohort.isin(common)]),('corrected_fixed34',newu)]:
 for k,g in d.groupby('pairs_requested'):
  rows.append(dict(stage=stage,pairs_requested=k,n_units=len(g),**{v:g[v].mean() for v in vars}))
st=pd.DataFrame(rows);_public_table(st).to_csv(t/'parsimony_change_components.tsv',sep='\t',index=False)
# Full per-generation evidence table, including reconstructed reference statistics.
gens=[]
for label,f in [('legacy_variance_reconstruction','current_predictions_cluster_variance_sensitivity.tsv'),('current_all_primary','current_full_method_summary.tsv'),('current_revised','current_method_summary.tsv'),('corrected','prepared_method_summary.tsv')]:
 d=pd.read_csv(t/_public_analysis_directory(f),sep='\t');d['generation']=label
 d['input_id']=hashlib.sha256((b/(_public_analysis_directory('revision_r3/cohort_cache_r3.rds' if label=='corrected' else 'cohort_cache.rds'))).read_bytes()).hexdigest()
 d['prediction_input']=('revision_r4/tables/within_predictions.tsv.gz' if label=='corrected' else 'within/within_predictions.tsv.gz')
 d['eligibility_rule']='all 5 seeds and all K folds; complete finite scores; SE>0; pooling k>=5; exclude whole blood'+('; exclude GSE145176' if label in ['current_revised','corrected'] else '')
 d['reference_impl']='signed-rCLR, training-only top-20 panel; OmicSelector 2.6.3 source-equivalent/2.6.5; baseline rclr_panel'
 d['variance_estimator']=('group-cluster influence, diagnostic reconstruction' if label=='legacy_variance_reconstruction' else 'corrected repeated CV: (1/(5K)+mean test/train ratio)*variance of fold lifts')+'; REML and Knapp-Hartung'
 d['lineage_status']='legacy headline numerically reproduced; original manifest-bound bytes not recovered' if label=='legacy_variance_reconstruction' else 'current input bytes hashed; estimates reproduced from stated inputs'
 gens.append(d)
g=pd.concat(gens);_public_table(g).to_csv(t/'generation_reconciliation.tsv',sep='\t',index=False)
# Hash only explicitly used input files; complete source snapshot is separately retained.
paths=[b/_public_analysis_directory(x) for x in ['cohort_cache.rds','cohort_cache_digest.tsv','cohort_eligibility.tsv','within/within_cells.tsv','within/within_predictions.tsv.gz','summary/within_performance.tsv','summary/method_summary.tsv','revision_r2/revised_method_summary.tsv','revision_r3/cohort_cache_r3.rds','revision_r3/cohort_cache_digest_r3.tsv','revision_r3/code_hashes.txt','revision_r3/runtime_versions.tsv','revision_r3/within_r3/within_cells.tsv','revision_r3/within_r3/within_predictions.tsv.gz','revision_r3/summary_r3/within_performance.tsv','revision_r3/summary_r3/revised_method_summary.tsv','revision_r3/deployment_audit_unseen.tsv']]
paths += [r/_public_analysis_directory(x) for x in ['results/manifests/analysis_a_singlesample_public_benchmark.json','source_data/within_method_meta.tsv','code/analyses/run_within_GSE204951_revision_r3.R','code/analyses/build_revision_r3_cache.R','code/analyses/summarize_singlesample_public_benchmark.R','code/methods/singlesample_roster_benchmark.R','code/methods/matched_null_benchmark.R','code/methods/paired_auc_diff_se.R','code/methods/singlesample_public_design.R','analysis/opportunity/strengthening/parsimony_gradient/parsimony.tsv','analysis/opportunity/strengthening/parsimony_gradient/METHODS_NOTE.md','analysis/opportunity/strengthening/parsimony_gradient/dependency_derive/candidate_predictions.tsv.gz','analysis/opportunity/strengthening/batch_mitigation/agreement.tsv','analysis/opportunity/strengthening/batch_mitigation/contrast.tsv','analysis/opportunity/strengthening/batch_mitigation/METHODS_NOTE.md']]
rows=[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()) for p in paths]
_public_table(pd.DataFrame(rows)).to_csv(t/'input_manifest.tsv',sep='\t',index=False)
for f in ['cohort_cache_digest.tsv','revision_r3/cohort_cache_digest_r3.tsv']:
 d=pd.read_csv(b/_public_analysis_directory(f),sep='\t').iloc[0];actual=next(v['sha256'] for v in rows if v['path']==d.path);assert actual==d.sha256
print('input digests match; generation table and parsimony decomposition written')
print(st[st.pairs_requested.isin([1,10])].to_string(index=False))
