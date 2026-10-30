source('./analyses/robustness/scripts/common.R')
R2<-file.path(M,'companion_dependence/out/runs/computing_bundle');R2b<-file.path(M,'companion_sensitivity/out')
bo<-rd(file.path(R2b,'computing_bundle/partB_outcomes.tsv'));bs<-rd(file.path(R2b,'computing_bundle/partB_summary.tsv'));env<-rd(file.path(R2b,'actual_task_environment.tsv'))
sel<-grepl('^B4_(sensitivity|specificity)_(10|50|90)$',bo$metric)
ops<-rbindlist(list(bo[sel][,scope:=paste0('split_',repetition)],bs[grepl('^B4_(sensitivity|specificity)_(10|50|90)$',metric)][,scope:='mean_five_splits']),fill=TRUE)
ops[,`:=`(experiment='five_split_deployment',fraction=as.numeric(sub('.*_','',metric))/100,metric=sub('^B4_(sensitivity|specificity)_.*','\\1',metric),source=file.path(R2b,'computing_bundle'),computing_host=env$host[match(paste0('B_',dataset),env$task)],device='CPU',threads=1L,score_package='Saved companion-analysis fits use retained-intercept reconstruction with OmicSelector 2.6.8.',interval='Intervals use 2000 joint person-bootstrap percentile draws and the equal mean of weighted batch rates. Scores and batches are fixed.')]
# Inclusion is based on specimen rows. GSE204951 has repeated profiles per person.
elig<-rd(file.path(R2,'partB_datasets.tsv'));elig[,`:=`(row_type='admission',experiment='five_split_deployment',rule='At least 18 held-out cases and 18 controls in every split; Part A reproduction gate; 32 disease-discrimination datasets screened',computing_host=HOST,device='CPU',threads=1L)]
fits<-rd(file.path(R2b,'computing_bundle/partB_fits.tsv'))
disc<-rbindlist(list(bo[metric=='B5'][,scope:=paste0('split_',repetition)],bs[metric=='B5'][,scope:='mean_five_splits']),fill=TRUE)
disc[,`:=`(row_type='metric',experiment='five_split_deployment',metric='mean_paired_call_discordance_10_vs_90',definition='For each specimen observed at both extremes. P10*(1-p90)+(1-p10)*p90. Equal specimen mean within split, then equal split mean.',computing_host=env$host[match(paste0('B_',dataset),env$task)],device='CPU',threads=1L)]
anyrows<-list();counts<-list();checks<-list()
for(id in unique(bo$dataset)){
 W<-rr(file.path(R2b,'bc_results',paste0('B_',id),'weights.rds'));pm<-rd(file.path(R2b,'bc_results',paste0('B_',id),'person_map.tsv'));stopifnot(ncol(W)==nrow(pm));colnames(W)<-pm$person
 all<-list()
 for(r in 1:5){p<-rd(file.path(R2,'units',id,paste0('B_',r),'pred.tsv.gz'))
  for(kind in c('STD','QN','RAW','CLR'))if(nrow(fits[dataset==id&repetition==r&model==kind&features_selected==0])){aa<-switch(kind,STD=c('FZ-STD','CN-STD'),QN=c('FZ-QN','CN-QN'),RAW=c('FZ-RAW','CN-COMBAT'),CLR='SS-CLR');p[arm%in%aa,call:=TRUE]}
  p[,person:=as.character(group)]
  if(id=='GSE204951'){co<-rr(file.path(ROOT,'analysis/singlesample_public_benchmark/revision_r4/GSE204951_prepared.rds'));p[,person:=as.character(co$group_id)[match(specimen,rownames(co$expr_per_sample))]]}
  stopifnot(!anyNA(p$person),all(p$person%in%colnames(W)))
  q<-p[,.(value=if(anyNA(call))NA_real_ else as.numeric(any(call)&any(!call)),n_appearances=.N,n_refusals=sum(is.na(call))),by=.(arm,specimen,person)]
  q[,repetition:=r];all[[r]]<-q
  checks[[length(checks)+1]]<-data.table(dataset=id,repetition=r,check='batch size and fraction',pass=all(p[,.(n=.N,cases=sum(label)),by=.(arm,fraction,batch)][,n==20&cases==20*fraction]))
 }
 q<-rbindlist(all)
 for(ar in unique(q$arm)){
  d<-q[arm==ar];draws<-list()
  for(r in 1:5){z<-d[repetition==r];v<-if(anyNA(z$value))rep(NA_real_,nrow(W))else wmean(z$value,z$person,W);draws[[r]]<-v
   anyrows[[length(anyrows)+1]]<-cbind(data.table(dataset=id,arm=ar,scope=paste0('split_',r),repetition=r,metric='proportion_unique_specimens_any_changed_call',n_specimens=nrow(z),n_people=uniqueN(z$person),n_changed=sum(z$value),n_refusals=sum(z$n_refusals)),as.data.table(summary_draw(v)))}
  v<-Reduce('+',draws)/5
  anyrows[[length(anyrows)+1]]<-cbind(data.table(dataset=id,arm=ar,scope='mean_five_splits',metric='proportion_unique_specimens_any_changed_call',n_specimens=uniqueN(d$specimen),n_people=uniqueN(d$person),n_refusals=sum(d$n_refusals)),as.data.table(summary_draw(v)))
  u<-d[,.(value=if(anyNA(value))NA_real_ else max(value),person=person[1]),by=specimen]
  v<-if(anyNA(u$value))rep(NA_real_,nrow(W))else wmean(u$value,u$person,W)
  anyrows[[length(anyrows)+1]]<-cbind(data.table(dataset=id,arm=ar,scope='union_of_five_splits',metric='proportion_unique_specimens_changed_within_any_split',n_specimens=nrow(u),n_people=uniqueN(u$person),n_changed=sum(u$value)),as.data.table(summary_draw(v)))
 }
 cat(id,'any-call analysis complete\n');flush.console()
}
a<-rbindlist(anyrows,fill=TRUE);a[,`:=`(experiment='five_split_deployment',row_type='metric',definition='Discordance includes positive and negative calls across batch appearances of the same specimen within a split. Differences between fitted split models are excluded.',computing_host=HOST,device='CPU',threads=1L,status=ifelse(is.finite(estimate),'ok','not_estimable_due_to_missing_calls'))]
# Original six-dataset producer: independently calculated person-bootstrap intervals.
A<-file.path(ROOT,'analysis/case_mix_invariance');record(file.path(A,'scripts/run_analysis.R'));record(file.path(ROOT,'analyses/single_specimen/computing_bundle_headtohead/PROTOCOL.md'));record(file.path(M,'companion_dependence/scripts/combine_case.R'));orig<-list();od<-list()
ids<-c('TORAY_PUBLIC_CANONICAL','GSE64591','GSE94533_PUBLIC','GSE188627','GSE279209','GSE320031')
for(id in ids){folder<-if(id=='TORAY_PUBLIC_CANONICAL')file.path(ROOT,'analyses/case_mix_duplicates_casemix_dedup/out')else if(id=='GSE320031')file.path(M,'companion_dependence/out/runs/case_mix')else A
 for(ar in c('ws-rclr-trimmed','reo-ktsp','reo-pairratio','Within-batch_standardization','Test-batch_centring')){
  p<-rd(file.path(folder,'predictions',paste0(id,'__',ar,'.tsv.gz')));people<-sort(unique(p$specimen));W<-weights(people,202610260L+sum(utf8ToInt(id)))
  for(fr in c(.1,.5,.9))for(cl in 0:1){z<-p[case_fraction==fr&outcome==cl,.(num=sum(if(cl==1)decision else !decision),den=.N),by=specimen];w<-W[,match(z$specimen,colnames(W)),drop=FALSE];v<-as.numeric(w%*%z$num)/as.numeric(w%*%z$den)
   orig[[length(orig)+1]]<-cbind(data.table(dataset=id,arm=ar,scope='original_single_split',experiment='original_case_mix',fraction=fr,metric=if(cl==1)'sensitivity'else'specificity',n_specimens=nrow(z),n_appearances=sum(z$den),computing_host=HOST,device='CPU',threads=1L,source=folder,score_package=if(id=='GSE320031')'Training-only preparation uses OmicSelector 2.6.7.'else'The analysis uses saved case-mix predictions.',interval='Intervals use 2000 person-bootstrap percentile draws and the appearance-weighted rate. Scores and batches are fixed.'),as.data.table(summary_draw(v)))}
  z<-p[,.(any_changed=as.numeric(any(decision)&any(!decision)),p10=mean(decision[case_fraction==.1]),p90=mean(decision[case_fraction==.9])),by=specimen];z[,paired:=p10*(1-p90)+(1-p10)*p90]
  for(mt in c('paired','any_changed')){v<-wmean(z[[mt]],z$specimen,W);od[[length(od)+1]]<-cbind(data.table(dataset=id,arm=ar,scope='original_single_split',experiment='original_case_mix',row_type='metric',metric=if(mt=='paired')'mean_paired_call_discordance_10_vs_90'else'proportion_unique_specimens_any_changed_call',n_specimens=sum(is.finite(z[[mt]])),computing_host=HOST,device='CPU',threads=1L,source=folder),as.data.table(summary_draw(v)))}
 }
}
wr(rbindlist(c(list(ops),orig),fill=TRUE),'item1_case_mix_operating_points.tsv')
wr(rbindlist(c(list(elig,disc,a),od),fill=TRUE),'item2_case_mix_admission_and_discordance.tsv')
wr(rbindlist(checks),'evidence/case_mix_checks.tsv');stopifnot(all(rbindlist(checks)$pass));finish('case_mix')
