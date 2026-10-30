source(file.path(Sys.getenv('R2'),'scripts/computing_bundle/common.R'))
H<-if(Sys.getenv('M18_PREP_MODE')=='training')file.path(Sys.getenv('R2'),'out/computing_bundle')else file.path(ROOT,'analyses/single_specimen/computing_bundle_headtohead/out')
cl<-fread(file.path(OUT,'decision_rules_classification.tsv'));keep<-cl[class%in%c('a','b'),unit]
alias<-function(x){x[x=='deduplicated 3D-Gene serum unit']<-'TORAY_PUBLIC_CANONICAL';x[x=='GSE94533']<-'GSE94533_PUBLIC';x[x=='GSE259327 (diagnostic arm)']<-'GSE259327_DIAGNOSTIC';x}
import_functions(file.path(H,'scripts/aggregate.R'),'pool_safe')
restr<-list();pooled<-list()
a<-fread(file.path(H,'partA_contrasts_dataset.tsv'))[dataset%in%keep]
restr[[1]]<-a[,.(analysis='Part A',unit=dataset,metric=contrast,estimate,standard_error=corrected_se,ci_lo,ci_hi,status,source_file=file.path(H,'partA_contrasts_dataset.tsv'))]
pa<-a[,pool_safe(data.table(estimate=estimate,se=corrected_se)),by=.(type,contrast,left,right)]
pa[,`:=`(analysis='Part A',metric=contrast,summary='REML/Knapp-Hartung',estimate=pooled_lift)]
pooled[[1]]<-pa
b<-fread(file.path(H,'partB_contrasts_dataset.tsv'))[dataset%in%keep]
restr[[2]]<-b[,.(analysis='Part B deployment AUC',unit=dataset,metric=contrast,estimate,standard_error=se,ci_lo,ci_hi,status,source_file=file.path(H,'partB_contrasts_dataset.tsv'))]
pb<-b[,pool_safe(.SD),by=.(contrast,left,right)];pb[,`:=`(analysis='Part B deployment AUC',metric=contrast,summary='REML/Knapp-Hartung',estimate=pooled_lift)];pooled[[2]]<-pb
cmfile<-if(Sys.getenv('M18_PREP_MODE')=='training')file.path(Sys.getenv('R2'),'out/case_mix/results_all_units.tsv')else file.path(TASK,'case_mix_duplicates_casemix_dedup/out/results_all_units.tsv');cm<-fread(cmfile)[cohort%in%keep & case_fraction==.1]
cmr<-melt(cm,id.vars=c('cohort','scorer','status'),measure.vars=c('all_batches_decision_change_prop','extremes_pairwise_flip_rate','extremes_any_flip_prop'),variable.name='metric',value.name='estimate')
restr[[3]]<-cmr[,.(analysis='Case mix call changes',unit=cohort,metric=paste(scorer,metric,sep=': '),estimate,status,source_file=cmfile)]
pooled[[3]]<-cmr[,.(k=sum(is.finite(estimate)),estimate=median(estimate,na.rm=TRUE),minimum=min(estimate,na.rm=TRUE),maximum=max(estimate,na.rm=TRUE)),by=.(scorer,metric)][,`:=`(analysis='Case mix call changes',summary='Across-unit summaries give the median and range.')]
s6<-fread(file.path(SYN,'method_cohort_performance.tsv'));setnames(s6,names(s6)[1:7],c('method','unit','repetitions','complete','auc','lo','hi'));s6[,unit:=alias(unit)]
r<-s6[method=='ws-rclr-trimmed'&unit%in%keep]
if(Sys.getenv('M18_PREP_MODE')=='training'){
 import_functions(file.path(SYN,'scripts/published/0_summarize_singlesample_public_benchmark.R'),c('placements','cluster_auc_summary'))
 evidence<-list()
 for(uid in intersect(c('GSE320031','GSE83977'),keep)){
  # First verify the original group-cluster AUC interval against method_cohort_performance.
  calc<-function(base){p<-fread(file.path(base,'units',uid,'partA_pred.tsv.gz'))[arm=='R'];p[,`:=`(score_m=score,score_b=score,y=label,group_id=group,sample_id=row_id)];cluster_auc_summary(p,c('seed','fold'),TRUE,5L*as.integer(srb_expected_outer_k(uid)))}
  orig<-calc(file.path(Sys.getenv('R2'),'out/gate_runs/computing_bundle'));saved<-r[unit==uid]
  delta<-max(abs(c(orig$auc_method-saved$auc,orig$auc_method_ci_lo-saved$lo,orig$auc_method_ci_hi-saved$hi)))
  evidence[[uid]]<-data.table(unit=uid,max_difference=delta,tolerance=1e-9,pass=delta<=1e-9)
  stopifnot(delta<=1e-9)
  fresh<-calc(H);r[unit==uid,`:=`(auc=fresh$auc_method,lo=fresh$auc_method_ci_lo,hi=fresh$auc_method_ci_hi)]
 }
 write_tsv(rbindlist(evidence),'reference_interval_reproduction.tsv')
}

restr[[4]]<-r[,.(analysis='Within-cohort reference AUC',unit,metric='ws-rclr-trimmed',estimate=auc,ci_lo=lo,ci_hi=hi,status=ifelse(complete=='Yes','ok','not_estimable'),source_file=file.path(SYN,'method_cohort_performance.tsv'))]
pooled[[4]]<-r[,.(analysis='Within-cohort reference AUC',metric='all class a/b units',summary='Across-unit summaries give the median and range. Absolute AUCs are not meta-analyzed.',k=sum(is.finite(auc)),estimate=median(auc,na.rm=TRUE),mean=mean(auc,na.rm=TRUE),minimum=min(auc,na.rm=TRUE),maximum=max(auc,na.rm=TRUE))]
pooled[[5]]<-r[unit!='GSE145176',.(analysis='Within-cohort reference AUC',metric='active-cancer class a/b units',summary='Across-unit summaries give the median and range. The active-cancer synthesis excludes GSE145176.',k=sum(is.finite(auc)),estimate=median(auc,na.rm=TRUE),mean=mean(auc,na.rm=TRUE),minimum=min(auc,na.rm=TRUE),maximum=max(auc,na.rm=TRUE))]
rr<-rbindlist(restr,fill=TRUE);if(Sys.getenv('M18_PREP_MODE')=='training')rr[analysis=='Within-cohort reference AUC'&unit%in%c('GSE320031','GSE83977'),source_file:=file.path(H,'units',unit,'partA_pred.tsv.gz')];rr<-merge(rr,cl[,.(unit,class)],by='unit',all.x=TRUE);write_tsv(rr,'benchmark_restriction.tsv')
pp<-rbindlist(pooled,fill=TRUE);pp[,status:=ifelse(is.finite(estimate),'ok','not_estimable_fewer_than_five_or_pooling_failure')];write_tsv(pp,'benchmark_pooled.tsv')
write_tsv(cl[,.(unit,class,included_in_restriction=unit%in%keep,partA_present=unit%in%a$dataset,partB_present=unit%in%b$dataset,case_mix_present=unit%in%cm$cohort,reference_present=unit%in%r$unit)],'benchmark_unit_dispositions.tsv')
print(pp[,.(analysis,metric,k,estimate,ci_lo,ci_hi,minimum,maximum)])
