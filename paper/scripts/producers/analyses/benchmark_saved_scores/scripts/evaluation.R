# Evaluation-only exclusions. This file is never sourced by gate.R.
exclusion_path<-file.path(TASK,'cohort_duplicates_toray_dedup_check/evaluation_exclusions.tsv')
EXCLUSIONS<-fread(exclusion_path)
stopifnot(nrow(EXCLUSIONS)==396L,!anyDuplicated(EXCLUSIONS[,.(unit,row_name)]),
          EXCLUSIONS[unit==TORAY,.N]==394L,EXCLUSIONS[unit=='PRJEB30542',.N]==2L)
# Toray prediction sample_id is accession-qualified; row_name is the exact identifier.
EXCLUSIONS[,`:=`(cohort=unit,prediction_sample_id=row_name)]
stopifnot(all(EXCLUSIONS$sample_id==sub('^.*::','',EXCLUSIONS$prediction_sample_id)))
# Preserve the number of training groups after evaluation-only filtering.
# The published function body remains untouched on disk and in the v6 gate.
eval_body<-deparse(body(uncorrected_cluster_auc_summary),width.cutoff=500L)
eval_line<-grep('n_train_groups <- total_groups - n_test_groups',eval_body,fixed=TRUE)
stopifnot(length(eval_line)==1L)
eval_body[eval_line]<-'        n_train_groups <- if ("n_train_groups_original" %in% names(x)) unique(x$n_train_groups_original) else total_groups - n_test_groups; if (length(n_train_groups) != 1L || anyNA(n_train_groups) || !is.finite(n_train_groups) || n_train_groups <= 0L) stop("Missing, inconsistent or invalid training-group count in stratum ", s)'
body(uncorrected_cluster_auc_summary)<-parse(text=paste(eval_body,collapse='\n'))[[1L]]
filter_evaluation<-function(p,c,method){
 p<-copy(p); audit<-list(); matches<-list()
 for(cohort_value in unique(EXCLUSIONS$cohort))for(seed_value in SEEDS){
  e<-EXCLUSIONS[cohort==cohort_value];cc<-c[cohort==cohort_value & seed==seed_value]
  pp<-p[cohort==cohort_value & seed==seed_value]
  eligible<-nrow(cc)==1L && isTRUE(cc$eligible)
  if(!eligible){
   stopifnot(nrow(pp)==0L)
   audit[[length(audit)+1L]]<-data.table(method=method,cohort=cohort_value,seed=seed_value,status=if(!nrow(cc)||is.na(cc$eligible))'pending_no_predictions'else'ineligible_no_predictions',expected=nrow(e),matched=0L,pass=TRUE)
   next
  }
  stopifnot(!anyDuplicated(pp[,.(sample_id)]),!anyDuplicated(pp[,.(fold,sample_id)]))
  z<-merge(e[,.(cohort,prediction_sample_id,group_id,pair_id)],pp[,.(cohort,prediction_sample_id=sample_id,prediction_group_id=group_id,fold)],by=c('cohort','prediction_sample_id'),all.x=TRUE,allow.cartesian=TRUE)
  ok<-nrow(z)==nrow(e)&&!anyNA(z$fold)&&all(z$group_id==z$prediction_group_id)
  if(!ok)stop('Exclusions do not match exactly one held-out prediction per method/seed: ',method,' ',cohort_value,' ',seed_value)
  z[,`:=`(method=method,seed=seed_value)];matches[[length(matches)+1L]]<-z
  audit[[length(audit)+1L]]<-data.table(method=method,cohort=cohort_value,seed=seed_value,status='matched_once_in_assigned_fold',expected=nrow(e),matched=nrow(z),pass=ok)
 }
 write_tsv(rbindlist(audit),paste0('evidence/methods/',method,'_exclusion_checks.tsv'))
 hit<-rbindlist(matches,fill=TRUE)
 if(!nrow(hit))hit<-data.table(cohort=character(),prediction_sample_id=character(),group_id=character(),pair_id=character(),prediction_group_id=character(),fold=integer(),method=character(),seed=integer())
 write_tsv(hit,paste0('evidence/methods/',method,'_exclusion_matches.tsv'))
 if(!nrow(p))return(p)
 if(!any(p$cohort%in%EXCLUSIONS$cohort))return(p)
 # Fill every stratum before binding: unaffected folds use total minus test
 # groups; evaluation-filtered folds retain this original training count.
 p[,n_train_groups_original:=uniqueN(group_id),by=.(method,cohort,seed)]
 p[,n_train_groups_original:=n_train_groups_original-uniqueN(group_id),by=.(method,cohort,seed,fold)]
 affected<-p[cohort%in%EXCLUSIONS$cohort]
 before<-affected[,.(n_before=.N,n_train=unique(n_train)),by=.(method,cohort,seed,fold)]
 affected<-affected[!EXCLUSIONS,on=.(cohort,sample_id=prediction_sample_id)]
 affected[,n_test:=.N,by=.(method,cohort,seed,fold)]
 after<-affected[,.(n_after=.N,n_train_after=unique(n_train),n_pos=sum(y==1L),n_neg=sum(y==0L)),by=.(method,cohort,seed,fold)]
 a<-merge(before,after,by=c('method','cohort','seed','fold'),all=TRUE)
 stopifnot(!anyNA(a),all(a$n_train==a$n_train_after),all(a$n_pos>=3L & a$n_neg>=3L),sum(a$n_before-a$n_after)==nrow(hit))
 a[,excluded:=n_before-n_after];write_tsv(a,paste0('evidence/methods/',method,'_evaluation_folds.tsv'))
 rbindlist(list(p[!cohort%in%EXCLUSIONS$cohort],affected),fill=TRUE)
}
