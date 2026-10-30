source('./analyses/robustness/scripts/common.R')
suppressPackageStartupMessages(library(pROC))
A<-file.path(M,'figures/work/staging/haemolysis');x<-rd(file.path(A,'hemolysis_anchor_values.tsv'));old<-rd(file.path(A,'hemolysis_index_auc.tsv'));ids<-c('GSE266859','GSE117063','GSE227778','GSE188627');res<-list()
for(id in ids){z<-x[cohort==id&selected==TRUE&is.finite(index)];stopifnot(!anyDuplicated(z$group_id));ro<-roc(z$y,z$index,levels=c(0,1),direction='<',quiet=TRUE);raw<-as.numeric(pROC::auc(ro));ci<-as.numeric(ci.auc(ro,method='delong'));sgn<-if(raw<.5)-1 else 1
 oriented<-roc(z$y,sgn*z$index,levels=c(0,1),direction='<',quiet=TRUE);oc<-as.numeric(ci.auc(oriented,method='delong'));stopifnot(abs(raw-old[cohort==id,auc])<1e-12)
 res[[id]]<-data.table(dataset=id,n_cases=sum(z$y==1),n_controls=sum(z$y==0),index='log(miR-451a+0.5)-log(miR-23a-3p+0.5)',direction=if(sgn==1)'higher index is case-like'else'lower index is case-like',orientation_multiplier=sgn,raw_high_case_auc=raw,raw_ci_lo=ci[1],raw_ci_hi=ci[3],direction_oriented_auc=as.numeric(pROC::auc(oriented)),ci_lo=oc[1],ci_hi=oc[3],interval='95% DeLong.',orientation_note='Direction follows the observed case-control association. This orientation describes the same data and does not estimate externally selected diagnostic performance.',computing_host=HOST,device='CPU',threads=1L)
}
wr(rbindlist(res),'item5_haemolysis.tsv');finish('haemolysis')
