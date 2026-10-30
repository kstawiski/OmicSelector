source(file.path(Sys.getenv('OUT'),'scripts/common.R'))
datasets<-fread(file.path(OUT,'datasets.tsv'));gate<-fread(file.path(OUT,'partA_reproduction_gate.tsv'));info<-list();seeds<-list();jobs<-list()
for(di in seq_len(nrow(datasets))) {
 id<-datasets$dataset[di];if(id!=Sys.getenv('REFIT_ACCESSION'))next;z<-co_data(get_cohort(id));keep<-!z$ids%in%exclusions(id);z$X<-z$X[keep,,drop=FALSE];for(nm in c('y','g','ids'))z[[nm]]<-z[[nm]][keep]
 sizes<-list();splits<-list()
 for(r in 1:5) {
  base<-20261017L+100L*r+di
  sd<-data.table(dataset=id,repetition=r,split=base,train_cap=base+10000L,test_cap=base+20000L,model=base+30000L,inner_folds=base+40000L,batch_10=base+50000L,batch_50=base+60000L,batch_90=base+70000L,bootstrap=base+80000L)
  f<-os_make_grouped_stratified_folds(z$y,group_id=z$g,n_folds=5L,seed=base);te<-sort(unlist(f[1:2]));tr<-setdiff(seq_along(z$y),te)
  if(id=='TORAY_PUBLIC_CANONICAL'){tr<-cap_groups(tr,z$y,z$g,3000L,sd$train_cap);te<-cap_groups(te,z$y,z$g,1500L,sd$test_cap)}
  stopifnot(!length(intersect(z$g[tr],z$g[te])))
  sizes[[r]]<-data.table(dataset=id,repetition=r,n_available=length(z$y),n_excluded=sum(!keep),n_train=length(tr),n_test=length(te),train_cases=sum(z$y[tr]==1),train_controls=sum(z$y[tr]==0),test_cases=sum(z$y[te]==1),test_controls=sum(z$y[te]==0),group_overlap=length(intersect(z$g[tr],z$g[te])))
  splits[[r]]<-list(tr=tr,te=te,seeds=sd,ids=z$ids)
  seeds[[length(seeds)+1L]]<-sd
 }
 q<-rbindlist(sizes);eligible<-all(q$test_cases>=18L&q$test_controls>=18L)&&isTRUE(gate[dataset==id,pass]);q[,`:=`(eligible=eligible,reason=if(!isTRUE(gate[dataset==id,pass]))'reproduction_gate_failed'else if(!eligible)'fewer_than_18_cases_or_controls_in_at_least_one_split'else'ok')]
 info[[di]]<-q
 if(eligible)for(r in 1:5){saveRDS(splits[[r]],file.path(OUT,'units',id,paste0('split_B_',r,'.rds')));jobs[[length(jobs)+1L]]<-data.table(dataset=id,repetition=r)}
}
write_tsv(rbindlist(info),'partB_datasets.tsv');write_tsv(rbindlist(seeds),'seeds.tsv');write_tsv(rbindlist(jobs),'partB_jobs.tsv');log_line('PARTB eligibility complete;',length(jobs),'dataset-repetition jobs')
