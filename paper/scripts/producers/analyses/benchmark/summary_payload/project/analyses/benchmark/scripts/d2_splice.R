# Replace only GSE198834; retain the full benchmark_sensitivity_summary scientific summarization and pooling code.
CH<-file.path(Sys.getenv('D2_PAYLOAD'),'changed');cc<-fread(file.path(CH,'primary_cells.tsv'));pp<-fread(file.path(CH,'primary_predictions.tsv.gz'))
stopifnot(nrow(cc)==250L,uniqueN(cc$method)==50L,setequal(cc$seed,SEEDS),all(cc$cohort=='GSE198834'))
znew<-summarize_cells(cc,pp,CJ(method=unique(cc$method),cohort='GSE198834'))
sc<-fread(file.path(CH,'sensitivity_cells.tsv'));sp<-fread(file.path(CH,'sensitivity_predictions.tsv.gz'));zsens<-summarize_cells(sc,sp,CJ(method=unique(sc$method),cohort='GSE198834'))
write_tsv(znew,'evidence/GSE198834_performance.tsv');write_tsv(zsens,'evidence/GSE198834_seed_sensitivity.tsv')
# Restore source paths to the immutable staged replacement bundles for matched-lasso validation.
cc[,source_bundle:=file.path(CH,'bundles',paste0('primary__',method,'__GSE198834__seed_',seed),'within_bundle.rds')]
cc[,`:=`(source_status='completed',source_layer='count_matched_labels',source_device=device,source_threads=as.character(threads))]
local({source(file.path(B,'scripts/d2_matched.R'),local=TRUE)})
matched_new<-fread(file.path(B,'out/matched_baseline_dataset.tsv'))
d2_splice_rows<-function(d,z,keys=c('method','cohort','seed')){
 z<-copy(z);keys<-intersect(keys,names(d));z<-z[d[,..keys],on=keys,nomatch=0L]
 for(n in setdiff(names(d),names(z)))z[,(n):=d[[n]][NA_integer_]]
 rbindlist(list(d[!z[,..keys],on=keys],z[,names(d),with=FALSE]),use.names=TRUE)
}
d2_read_old<-function(f){
 d<-fread(file.path(OLD,f));if(!'cohort'%in%names(d)||!any(d$cohort=='GSE198834'))return(d)
 if(grepl('foundation_seed/',f)){
  if(grepl('_performance.tsv$',f))return(splice_summary(d,zsens))
  if(grepl('_cells.tsv$',f))return(d2_splice_rows(d,sc))
 }
 if(f=='evidence/matched_baseline_dataset.tsv')return(d2_splice_rows(d,matched_new,c('method','cohort')))
 if(grepl('(^cells_union|within_(performance|linear))',f)){
  if('seed'%in%names(d))return(d2_splice_rows(d,cc))
  return(splice_summary(d,znew))
 }
 d
}
