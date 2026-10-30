source('./analyses/robustness/scripts/common.R')
bad<-c('GSM5687109','GSM5687309');res<-list();labelsets<-list(cached_labels_11=c(1,1),source_label_01=c(0,1),alternative_10=c(1,0),alternative_00=c(0,0),exclude_both=c(NA,NA))
calc<-function(p,s,scenario,analysis,repetition){
 if(scenario=='exclude_both'){p<-p[!specimen%in%bad];s<-s[!specimen%in%bad]}else {for(j in 1:2){p[specimen==bad[j],label:=labelsets[[scenario]][j]];s[specimen==bad[j],label:=labelsets[[scenario]][j]]}}
 rows<-list()
 for(ar in unique(p$arm)){
  a<-p[arm==ar];b<-s[arm==ar];st<-a[,.(p10=mean(call[fraction==.1]),p90=mean(call[fraction==.9]),any_changed=if(anyNA(call))NA_real_ else as.numeric(any(call)&any(!call))),by=specimen];st[,paired:=p10*(1-p90)+(1-p10)*p90]
  for(mt in c('full_set_auc','pooled_extremes_auc','mean_batch_auc_10','mean_batch_auc_50','mean_batch_auc_90','mean_paired_call_discordance_10_vs_90','proportion_unique_specimens_any_changed_call')){
   value<-switch(mt,full_set_auc=auc(b$score,b$label),pooled_extremes_auc=auc(a[fraction%in%c(.1,.9),score],a[fraction%in%c(.1,.9),label]),mean_paired_call_discordance_10_vs_90=mean(st$paired,na.rm=TRUE),proportion_unique_specimens_any_changed_call=mean(st$any_changed),{fr<-as.numeric(sub('.*_','',mt))/100;mean(a[fraction==fr,.(a=auc(score,label)),by=batch]$a)})
   rows[[length(rows)+1]]<-data.table(dataset='GSE188627',analysis,scenario,repetition,scope=if(repetition>0)paste0('split_',repetition)else'original_single_split',arm=ar,metric=mt,estimate=value,n_test=uniqueN(a$specimen),label_109=labelsets[[scenario]][1],label_309=labelsets[[scenario]][2],computing_host=HOST,device='CPU',threads=1L,model_policy='Fitted models, thresholds and memberships are fixed. Label-sensitivity assignments and exclusions apply only to evaluation rows.')
  }
 }
 rbindlist(rows)
}
for(r in 1:5){d<-file.path(M,'companion_dependence/out/runs/computing_bundle/units/GSE188627',paste0('B_',r));p<-rd(file.path(d,'pred.tsv.gz'));s<-rd(file.path(d,'single.tsv.gz'))[context=='full_set'];for(scenario in names(labelsets))res[[length(res)+1]]<-calc(copy(p),copy(s),scenario,'fixed_scores_five_split_deployment',r)}
A<-file.path(ROOT,'analysis/case_mix_invariance');pp<-list()
for(ar in c('ws-rclr-trimmed','reo-ktsp','reo-pairratio','Within-batch_standardization','Test-batch_centring')){p<-rd(file.path(A,'predictions',switch(as.character(ar),'specified_scores_five_split_deployment'='GSE188627__specified_scores_five_split_deployment.tsv.gz','specified_scores_original_case_mix'='GSE188627__specified_scores_original_case_mix.tsv.gz',paste0('GSE188627__',ar,'.tsv.gz'))));setnames(p,c('case_fraction','outcome','decision'),c('fraction','label','call'));p[,arm:=ar];pp[[ar]]<-p}
p<-rbindlist(pp);s<-p[0,.(specimen,label,arm,score,call)]
for(scenario in names(labelsets))res[[length(res)+1]]<-calc(copy(p),copy(s),scenario,'fixed_scores_original_case_mix',0L)
x<-rbindlist(res);x<-x[!(analysis=='fixed_scores_original_case_mix'&metric=='full_set_auc')]
avg<-x[analysis=='fixed_scores_five_split_deployment',.(estimate=if(all(is.finite(estimate)))mean(estimate)else NA_real_,n_test=sum(n_test)),by=.(dataset,analysis,scenario,arm,metric,computing_host,device,threads,model_policy)];avg[,`:=`(repetition=0L,scope='mean_five_splits')];x<-rbindlist(list(x,avg),fill=TRUE)
base<-x[scenario=='cached_labels_11',.(analysis,repetition,scope,arm,metric,baseline=estimate)];x<-merge(x,base,by=c('analysis','repetition','scope','arm','metric'));x[,change_from_cached_labels:=estimate-baseline]
# The index AUC is also affected by these labels.
h<-rd(file.path(M,'figures/work/staging/haemolysis/hemolysis_anchor_values.tsv'))[cohort=='GSE188627'&selected==TRUE];hr<-list()
for(scenario in names(labelsets)){z<-copy(h);if(scenario=='exclude_both')z<-z[!sample_id%in%bad]else for(j in 1:2)z[sample_id==bad[j],y:=labelsets[[scenario]][j]]
 hr[[scenario]]<-data.table(dataset='GSE188627',analysis='haemolysis_index',scenario,scope='all_selected_specimens',repetition=0L,arm='haemolysis_index',metric='direction_oriented_auc',estimate=auc(z$index,z$y),n_test=nrow(z),computing_host=HOST,device='CPU',threads=1L)}
h<-rbindlist(hr);h[,baseline:=h[scenario=='cached_labels_11',estimate]];h[,change_from_cached_labels:=estimate-baseline]
wr(rbindlist(list(x,h),fill=TRUE),'evidence/label_assignment_summary_results.tsv')
primary<-data.table(dataset='GSE188627',row_type='primary_record_evidence',accession=bad,cached_label=1L,source_name=c('healthy','CRC'),disease_field=c('CRC','absent'),stage=c('NA','2'),title=c('NGSC059','CRCNT049'),resolution=c('Primary annotation fields conflict. Source-label sensitivity and exclusion analyses address this uncertainty.','The source name and stage 2 support CRC. A missing disease field does not imply a healthy label.'),source_url=paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',bad),retrieved_utc='2026-10-26',browser='persistent-chrome main, Poland',computing_host=HOST,device='CPU',threads=1L)
wr(primary,'evidence/label_primary_records.tsv');finish('label_summaries')
