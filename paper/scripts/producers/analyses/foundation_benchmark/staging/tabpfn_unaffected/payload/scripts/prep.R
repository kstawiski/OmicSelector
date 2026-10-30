# figures P3: descriptive preparation only. No classifier fitting or scoring.
# X_rows is a samples-by-features numeric matrix with unique row identifiers.
# The specification contains specified file membership and static feature identities;
# it contains no medians or retained-feature decisions from held-out values.
.m18_matrix_fit <- function(spec, X_rows) {
 stopifnot(is.matrix(X_rows),is.numeric(X_rows),!is.null(rownames(X_rows)),!anyDuplicated(rownames(X_rows)),!is.null(colnames(X_rows)))
 if(!is.null(spec$unavailable))stop(spec$unavailable)
 mode<-if(is.null(spec$mode))'training'else spec$mode
 if(!mode%in%c('original','original_linear','training'))stop('Unknown preparation mode')
 blocks<-spec$blocks
 rb<-if(length(blocks)==1L)rep(names(blocks),nrow(X_rows))else unname(spec$row_block[rownames(X_rows)])
 if(anyNA(rb))stop('File membership is missing for input rows')
 fits<-lapply(names(blocks),function(k){
  s<-blocks[[k]];r<-which(rb==k)
  if(!length(r))stop('No fitting rows for source block ',k)
  x<-X_rows[r,s$features,drop=FALSE]
  keep<-colSums(is.finite(x))>0L
  if(identical(s$retain,'complete'))keep<-colSums(!is.finite(x))==0L
  if(!is.null(s$minimum_detection_rate))keep<-keep & colMeans(is.finite(x)&x>0)>=s$minimum_detection_rate
  if(isTRUE(s$no_retention)&&mode%in%c('original','original_linear'))keep[]<-TRUE
  features<-colnames(x)[keep];x<-x[,features,drop=FALSE]
  v<-setNames(rep(0,length(features)),features)
  if(s$fill=='median') v<-vapply(seq_along(features),function(j){m<-suppressWarnings(stats::median(x[,j],na.rm=TRUE));if(is.finite(m))m else 0},0.)
  if(s$fill=='median_finite')v<-vapply(seq_along(features),function(j){m<-stats::median(x[is.finite(x[,j]),j]);if(is.finite(m))m else 0},0.)
  if(s$fill=='fixed40')v[]<-40
  if(s$fill=='fixed_ct')v[]<-s$fixed_value
  if(s$fill%in%c('maxct1','globalmaxct1')){fallback<-max(x[is.finite(x)])+1;v<-if(s$fill=='globalmaxct1')rep(fallback,length(features))else vapply(seq_along(features),function(j){z<-x[is.finite(x[,j]),j];if(length(z))max(z)+1 else fallback},0.)}
  if(s$fill=='halfmin_positive_global'){positive<-x[is.finite(x)&x>0];if(!length(positive))stop('No positive training value for half-minimum floor');v[]<-min(positive)/2}
  if(s$fill=='lod')v<-vapply(seq_along(features),function(j)min(x[is.finite(x[,j]),j])-log(sqrt(2),s$base),0.)
  if(s$fill=='none'&&any(!is.finite(x)))stop('Unspecified non-finite input handling')
  names(v)<-features
  list(features=features,fill=v,spec=s,fit_rows=rownames(x))
 });names(fits)<-names(blocks)
 outfeatures<-lapply(fits,function(f){if(is.null(f$spec$post_map))f$features else {mp<-unname(f$spec$post_map[f$features]);if(anyDuplicated(mp))sort(unique(mp))else mp}})
 common<-if(isTRUE(spec$intersection))Reduce(intersect,outfeatures)else unique(unlist(outfeatures,use.names=FALSE))
 list(spec=spec,mode=mode,blocks=fits,features=common,fit_rows=rownames(X_rows))
}
.m18_matrix_apply <- function(fit, X_rows) {
 stopifnot(is.matrix(X_rows),!is.null(rownames(X_rows)))
 rb<-if(length(fit$blocks)==1L)rep(names(fit$blocks),nrow(X_rows))else unname(fit$spec$row_block[rownames(X_rows)])
 if(anyNA(rb))stop('File membership is missing for application rows')
 ans<-matrix(NA_real_,nrow(X_rows),length(fit$features),dimnames=list(rownames(X_rows),fit$features))
 for(k in names(fit$blocks)){
  idx<-which(rb==k);if(!length(idx))next
  f<-fit$blocks[[k]];s<-f$spec
  if(!all(f$features%in%colnames(X_rows)))stop('Application matrix lacks fitted features')
  x<-X_rows[idx,f$features,drop=FALSE]
  if(!s$fill%in%c('none','post0'))for(j in seq_len(ncol(x))){bad<-!is.finite(x[,j]);if(isTRUE(s$nonpositive_missing))bad<-bad|x[,j]<=0;x[bad,j]<-f$fill[j]}
  conv<-s$convert
  if(isTRUE(fit$spec$T4)&&fit$mode=='training')conv<-'exp'
  if(conv=='ct')x<-s$base^((if(is.null(s$ct_reference))40 else s$ct_reference)-x)else if(conv=='negate')x<-s$base^(-x)else if(conv=='exp')x<-s$base^x
  if(!isTRUE(s$no_clip)){x[!is.finite(x)]<-0;x[x<0]<-0}
  if(isTRUE(fit$spec$T4)&&fit$mode=='original_linear')x<-2^x
  if(!is.null(s$post_map)){
   mp<-unname(s$post_map[colnames(x)])
   if(anyDuplicated(mp)){gr<-split(seq_along(mp),mp);x<-do.call(cbind,lapply(gr,function(j)if(length(j)==1L)x[,j]else rowSums(x[,j,drop=FALSE],na.rm=TRUE)));colnames(x)<-names(gr);rownames(x)<-rownames(X_rows)[idx]}else colnames(x)<-mp
  }
  ans[idx,]<-x[,fit$features,drop=FALSE]
 }
 ans
}
# figure_export operates on abundance-scale matrices after each source cohort has been
# prepared with study_tables–figure_annotations on its own training rows. Target masks identify raw missing
# cells carried through the target's specimen-local conversion.
m18_transfer_fit <- function(training_matrices,features=Reduce(intersect,lapply(training_matrices,colnames))) {
 if(!length(training_matrices))stop('No training cohorts')
 X<-do.call(rbind,lapply(training_matrices,function(x)x[,features,drop=FALSE]))
 med<-vapply(seq_len(ncol(X)),function(j){v<-median(X[,j],na.rm=TRUE);if(is.finite(v))v else 0},0.)
 list(features=features,medians=setNames(med,features),scale='abundance')
}
m18_transfer_apply <- function(fit,X_target,missing_mask=!is.finite(X_target)){
 X<-X_target[,fit$features,drop=FALSE];mask<-missing_mask[,fit$features,drop=FALSE]
 for(j in seq_len(ncol(X)))X[mask[,j],j]<-fit$medians[j]
 X
}
# figure_annotations numeric preparation helpers. The caller must specify which analysis arm
# supplies the operation. Target-estimated ComBat parameters have no defined
# training-only estimate for an unseen batch; no such estimate is fabricated.
m18_constants_fit <- function(X,kind=c('standardize','quantile','median_center')){
 kind<-match.arg(kind)
 if(kind=='standardize'){s<-apply(X,2,sd);s[!is.finite(s)|s==0]<-1;return(list(kind=kind,mean=colMeans(X),sd=s))}
 if(kind=='median_center')return(list(kind=kind,median=apply(X,2,median)))
 list(kind=kind,reference=rowMeans(vapply(seq_len(nrow(X)),function(i)sort(X[i,]),numeric(ncol(X)))))
}
m18_constants_apply <- function(fit,X){
 if(fit$kind=='standardize')return(sweep(sweep(X,2,fit$mean,'-'),2,fit$sd,'/'))
 if(fit$kind=='median_center')return(sweep(X,2,fit$median,'-'))
 ans<-X
 for(i in seq_len(nrow(X))){v<-X[i,];o<-order(v);r<-rle(v[o]);ends<-cumsum(r$lengths);starts<-c(1L,head(ends,-1L)+1L);ans[i,o]<-rep(vapply(seq_along(ends),function(j)mean(fit$reference[starts[j]:ends[j]]),0.),r$lengths)}
 ans
}

# Multi-file units retain each file's original scale and full feature schema.
# Supply a named list of matrices for these units; this avoids inventing cells
# for platform features that were never measured in a different file.
m18_prep_fit <- function(spec,X_rows){
 if(is.matrix(X_rows))return(.m18_matrix_fit(spec,X_rows))
 if(!is.list(X_rows)||!setequal(names(X_rows),names(spec$blocks)))stop('Supply a matrix or the named source-block matrices')
 fits<-lapply(names(spec$blocks),function(k){s<-spec;s$blocks<-spec$blocks[k];s$intersection<-FALSE;.m18_matrix_fit(s,X_rows[[k]])});names(fits)<-names(spec$blocks)
 feats<-lapply(fits,`[[`,'features');common<-if(isTRUE(spec$intersection))Reduce(intersect,feats)else unique(unlist(feats,use.names=FALSE))
 list(multifile=TRUE,spec=spec,blocks=fits,features=common)
}
m18_prep_apply <- function(fit,X_rows){
 if(!isTRUE(fit$multifile))return(.m18_matrix_apply(fit,X_rows))
 z<-lapply(names(X_rows),function(k).m18_matrix_apply(fit$blocks[[k]],X_rows[[k]])[,fit$features,drop=FALSE])
 do.call(rbind,z)
}
m18_subset_rows <- function(input,ids){
 if(is.matrix(input$X))return(input$X[ids,,drop=FALSE])
 lapply(input$X,function(x)x[rownames(x)%in%ids,,drop=FALSE])
}
