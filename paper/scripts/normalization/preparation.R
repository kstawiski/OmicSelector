# Unaffected GSE117063: select rows from the publicly regenerated abundance matrix.
prepare_partition<-function(id,train_ids,test_ids){stopifnot(id=='GSE117063');z<-get_cohort(id);list(train=z$expr_per_sample[train_ids,,drop=FALSE],test=z$expr_per_sample[test_ids,,drop=FALSE],fit=list(preparation='approved unchanged unaffected cohort'))}
