source(file.path(Sys.getenv('R2'),'scripts/computing_bundle/common.R'))
R2<-Sys.getenv('R2');old<-file.path(ROOT,'analysis/case_mix_invariance');dedup<-file.path(TASK,'case_mix_duplicates_casemix_dedup/out');new<-file.path(R2,'out/runs/case_mix');baseline<-Sys.getenv('M18_CASE_BASELINE')=='1';OUT<-file.path(R2,if(baseline)'out/stored_case_mix'else'out/case_mix');dir.create(OUT,showWarnings=FALSE)
for(n in c('results.tsv','cohorts.tsv','batch_metrics.tsv','specimen_stability.tsv','model_summary.tsv','split_assignments.tsv','batch_membership.tsv','verification.tsv','glmnet_fold_assignments.tsv','matrix_checks.tsv','seeds.tsv','model_seeds.tsv')){
 if(!file.exists(file.path(old,n)))next
 x<-fread(file.path(old,n));for(d in if(baseline)dedup else c(dedup,new))if(file.exists(file.path(d,n))){v<-fread(file.path(d,n));if('cohort'%in%names(v)&&nrow(v))x<-rbindlist(list(x[!cohort%in%v$cohort],v),use.names=TRUE,fill=TRUE)}
 write_tsv(x,n);if(n%in%c('results.tsv','cohorts.tsv'))write_tsv(x,sub('.tsv','_all_units.tsv',n,fixed=TRUE))
}
res<-fread(file.path(OUT,'results_all_units.tsv'));bm<-fread(file.path(OUT,'batch_metrics.tsv'))
for(n in c('figure_panel_A_data.tsv','figure_panel_A_all_units.tsv','figure_panel_B_data.tsv')){
 x<-fread(file.path(dedup,n))
 if(!baseline&&n=='figure_panel_B_data.tsv'){
  v<-res[cohort=='GSE320031'&case_fraction==.1];for(i in which(x$cohort=='GSE320031')){method<-gsub('\n',' ',x$scorer[i],fixed=TRUE);method<-c('Trimmed-rCLR'='ws-rclr-trimmed','k-TSP'='reo-ktsp','Pair-ratio'='reo-pairratio','Within-batch standardization'='Within-batch standardization','Test-batch centring'='Test-batch centring')[[method]];a<-v[scorer==method];x$extremes_pairwise_flip_rate[i]<-a$extremes_pairwise_flip_rate;x$n_specimens_both_extremes[i]<-a$n_specimens_both_extremes;x$display[i]<-sprintf('%.1f',100*a$extremes_pairwise_flip_rate)}
 }else if(!baseline&&n=='figure_panel_A_all_units.tsv'){
  for(i in which(x$cohort=='GSE320031')){method<-gsub('\n',' ',x$scorer[i],fixed=TRUE);v<-res[cohort=='GSE320031'&scorer_label==method&case_fraction==x$case_fraction[i]];x$value[i]<-v[[paste0(tolower(x$metric[i]),'_mean')]]}
 }
 write_tsv(x,n)
}
ex<-res[case_fraction==.1,.(expectation='singleton_batch_invariance',cohort,scorer,status=ifelse(singleton_applicability!='applicable','not_applicable',ifelse(singleton_max_abs_difference==0&extremes_pairwise_flip_rate==0,'held','not_held')),value=singleton_max_abs_difference,discordance=extremes_pairwise_flip_rate)]
write_tsv(ex,'expectations.tsv')
# Descriptive original figure pooling: equal-cohort means and empirical batch percentiles.
po<-bm[,.(sensitivity_mean=mean(sensitivity),sensitivity_p025=quantile(sensitivity,.025),sensitivity_p975=quantile(sensitivity,.975),specificity_mean=mean(specificity),specificity_p025=quantile(specificity,.025),specificity_p975=quantile(specificity,.975),auc_mean=mean(auc),auc_p025=quantile(auc,.025),auc_p975=quantile(auc,.975)),by=.(scorer,case_fraction)];write_tsv(po,'pooled_descriptive.tsv')
