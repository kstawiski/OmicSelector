source('./analyses/robustness/scripts/common.R')
Q<-file.path(ROOT,'analyses/score_identity_identity/out/analysis/reassay_casemix')
p<-rd(file.path(Q,'specimen_decisions.tsv.gz'))[dataset%in%c('GSE268524','GSE71008')]
f<-rd(file.path(M,'reassay_design/out/frozen/specimen_decisions.tsv.gz'))[dataset%in%c('GSE268524','GSE71008')]
p<-rbindlist(list(p,f),fill=TRUE)
fix<-rd(file.path(M,'companion_sensitivity/out/evidence/reassay_intercept_recovery.tsv'))
for(i in seq_len(nrow(fix))){z<-fix[i];p[dataset==z$dataset&repetition==z$repetition&A==z$A&B==z$B&scorer==z$scorer,`:=`(score=z$recovered_intercept,threshold=-Inf,decision=TRUE,status='ok',reason=NA_character_,rejection_code=NA_character_)]}
stopifnot(all(p[status=='ok',decision==(score>=threshold)]))
res<-list();contrib<-list();validation<-list()
# Means give each person equal weight after averaging their specified, evaluable split contributions.
put<-function(z,ds,sc,design,metric,replicate='',pair='',fraction=NA_real_){
 z<-copy(z);stopifnot(all(c('specimen','repetition','value')%in%names(z)))
 if(!nrow(z))return(invisible(NULL))
 z[,`:=`(dataset=ds,scorer=sc,design=design,metric=metric,replicate=replicate,pair=pair,fraction=fraction)];contrib[[length(contrib)+1]]<<-z
 u<-z[,.(value=if(any(is.finite(value)))mean(value[is.finite(value)])else NA_real_,n_defined=sum(is.finite(value)),n_repetitions=.N),by=specimen]
 W<-weights(sort(unique(z$specimen)),202610261L+sum(utf8ToInt(ds)))
 v<-wmean(u$value,u$specimen,W)
 row<-data.table(dataset=ds,scorer=sc,design,metric,replicate,pair,fraction,scope='person_mean_across_splits',n_unique_specimens=nrow(u),n_evaluable_specimens=sum(is.finite(u$value)),n_undefined_specimens=sum(!is.finite(u$value)),n_partially_evaluable=sum(u$n_defined>0&u$n_defined<u$n_repetitions),n_split_specimen_contributions=nrow(z),n_defined_contributions=sum(is.finite(z$value)),raw_numerator=if('num'%in%names(z))sum(z$num,na.rm=TRUE)else NA_real_,raw_denominator=if('den'%in%names(z))sum(z$den,na.rm=TRUE)else NA_real_)
 res[[length(res)+1]]<<-cbind(row,as.data.table(summary_draw(v)))
 for(r in sort(unique(z$repetition))){a<-z[repetition==r];v<-wmean(a$value,a$specimen,W);rrow<-copy(row);rrow[,`:=`(scope=paste0('split_',r),n_unique_specimens=nrow(a),n_evaluable_specimens=sum(is.finite(a$value)),n_undefined_specimens=sum(!is.finite(a$value)),n_partially_evaluable=0L,n_split_specimen_contributions=nrow(a),n_defined_contributions=sum(is.finite(a$value)),raw_numerator=if('num'%in%names(a))sum(a$num,na.rm=TRUE)else NA_real_,raw_denominator=if('den'%in%names(a))sum(a$den,na.rm=TRUE)else NA_real_)];res[[length(res)+1]]<<-cbind(rrow,as.data.table(summary_draw(v)))}
}
for(ds in unique(p$dataset))for(sc in unique(p[dataset==ds,scorer])){
 d<-p[dataset==ds&scorer==sc]
 # The common matched 50% design is available for every scorer.
 for(des in c('main','third_replicate_eligible','singleton')){
  x<-d[arm==des];if(!nrow(x))next
  for(rp in unique(x$replicate))for(fr in unique(x[replicate==rp,case_fraction])){
   a<-x[replicate==rp & (if(is.na(fr))is.na(case_fraction)else case_fraction==fr)]
   for(cl in 0:1){z<-a[outcome==cl,.(num=sum((if(cl==1)decision else !decision)&status=='ok',na.rm=TRUE),den=sum(status=='ok'&!is.na(decision)),n_failed=sum(status!='ok'|is.na(decision))),by=.(specimen,repetition)];z[,value:=ifelse(den>0,num/den,NA_real_)];put(z,ds,sc,des,if(cl==1)'sensitivity_among_scored'else'specificity_among_scored',replicate=rp,fraction=fr)}
   z<-a[,.(num=sum(status!='ok'|is.na(decision)),den=.N),by=.(specimen,repetition)];z[,value:=num/den];put(z,ds,sc,des,'refusal_or_failed_score_rate',replicate=rp,fraction=fr)
  }
 }
 for(pr in if(ds=='GSE268524')list(c('1','2'),c('1','3'),c('2','3'))else list(c('1','2'))){
  des<-if('3'%in%pr)'third_replicate_eligible'else'main'
  a<-d[arm==des&replicate==pr[1]&case_fraction==.5,.(specimen,repetition,batch,outcome,s1=score,d1=decision,ok1=status=='ok')]
  b<-d[arm==des&replicate==pr[2]&case_fraction==.5,.(specimen,repetition,batch,s2=score,d2=decision,ok2=status=='ok')]
  z<-merge(a,b,by=c('specimen','repetition','batch'),all=TRUE);z[,ok:=ok1&ok2&!is.na(d1)&!is.na(d2)]
  validation[[length(validation)+1]]<-data.table(dataset=ds,scorer=sc,pair=paste(pr,collapse='-'),check='matched identical specimen/repetition/batch keys',pass=nrow(z)==nrow(a)&&nrow(z)==nrow(b)&&!anyNA(z$ok1)&&!anyNA(z$ok2),n_specimens=uniqueN(z$specimen))
  for(cl in 0:1)for(conditional in c(FALSE,TRUE)){
   zz<-z[outcome==cl];zz[,denrow:=ok & (!conditional | if(cl==1)d1 else !d1)];zz[,numrow:=ok&(if(cl==1)d1&!d2 else !d1&d2)]
   v<-zz[,.(num=sum(numrow,na.rm=TRUE),den=sum(denrow,na.rm=TRUE)),by=.(specimen,repetition)];v[,value:=ifelse(den>0,num/den,NA_real_)]
   put(v,ds,sc,des,paste0(if(cl==1)'cases_positive_to_negative'else'controls_negative_to_positive',if(conditional)'_conditional_on_initial_call'else'_among_class'),pair=paste(pr,collapse='-'),fraction=.5)
  }
  v<-z[,.(num=sum(ok&d1!=d2,na.rm=TRUE),den=sum(ok,na.rm=TRUE)),by=.(specimen,repetition)];v[,value:=ifelse(den>0,num/den,NA_real_)]
  put(v,ds,sc,des,'matched_reassay_call_discordance',pair=paste(pr,collapse='-'),fraction=.5)
  if(ds=='GSE268524'&&identical(pr,c('2','3'))){
   only<-v[,.(specimen,repetition,reassay=value)]
   lo<-d[arm=='joint'&replicate=='2'&case_fraction==.1,.(p=if(any(status!='ok'))NA_real_ else mean(decision)),by=.(specimen,repetition)]
   hi<-d[arm=='joint'&replicate=='3'&case_fraction==.9,.(p=if(any(status!='ok'))NA_real_ else mean(decision)),by=.(specimen,repetition)]
   both<-merge(lo,hi,by=c('specimen','repetition'),suffixes=c('10','90'),all=TRUE);both[,combined:=p10*(1-p90)+(1-p10)*p90]
   paired<-merge(only,both,by=c('specimen','repetition'),all=TRUE);stopifnot(uniqueN(paired$specimen)==416L)
   paired[!is.finite(reassay)|!is.finite(combined),c('reassay','combined'):=list(NA_real_,NA_real_)];paired[,difference:=combined-reassay]
   for(mt in c('reassay','combined','difference'))put(paired[,.(specimen,repetition,value=get(mt))],ds,sc,'same_416_specimens_replicate2_vs3',paste0('paired_',mt),pair='2-3',fraction=NA_real_)
  }
 }
 cat(ds,sc,'complete\n');flush.console()
}
x<-rbindlist(res,fill=TRUE);x[,`:=`(computing_host=HOST,device='CPU',threads=1L,bootstrap_seed_rule='202610261 + sum(UTF8(dataset))',interval='Intervals use 2000 joint person-bootstrap percentile draws. Repeated appearances and splits stay with each person. Models, thresholds and batches are fixed.',status=ifelse(is.finite(estimate),'ok','not_estimable'),score_source='Re-assay scores combine frozen comparators and retained-intercept reconstruction.')]
wr(x,'item3_reassay.tsv');wr(rbindlist(contrib,fill=TRUE),'evidence/reassay_person_contributions.tsv.gz');wr(rbindlist(validation),'evidence/reassay_checks.tsv');stopifnot(all(rbindlist(validation)$pass))
wr(p[,.(attempts=.N,unique_specimens=uniqueN(specimen)),by=.(dataset,scorer,arm,replicate,status,reason,rejection_code)],'evidence/reassay_refusal_codes.tsv')
finish('reassay')
