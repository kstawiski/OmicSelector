# Firth logistic calibration maximizes l(beta) + 0.5 log|X' W X|.
# The covariate is the same clipped affine score used by OmicSelector 2.6.9 Platt.
# Standardization changes the penalty by a beta-independent constant only.
fit_firth_calibration <- function(probs,y) {
 x<-pmax(1e-10,pmin(1-1e-10,as.numeric(probs)));y<-as.numeric(y)
 stopifnot(length(x)==length(y),all(is.finite(x)),all(y%in%c(0,1)),length(unique(y))==2L)
 center<-mean(x);scale<-sd(x)
 if(!is.finite(scale)||scale<1e-14)stop('FIRTH_UNIDENTIFIABLE_CONSTANT_PREDICTOR')
 X<-cbind(1,(x-center)/scale)
 state<-function(beta){
  eta<-as.numeric(X%*%beta);p<-plogis(eta);w<-p*(1-p);I<-crossprod(X,X*w)
  ch<-tryCatch(chol(I),error=function(e)NULL)
  if(is.null(ch))return(list(value=-Inf,score=c(Inf,Inf)))
  inv<-chol2inv(ch);h<-w*rowSums((X%*%inv)*X)
  ll<-sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))+sum(log(diag(ch)))
  U<-as.numeric(crossprod(X,y-p+h*(.5-p)))
  list(value=ll,score=U,information=I)
 }
 beta<-c(qlogis((sum(y)+.5)/(length(y)+1)),0);converged<-FALSE;iterations<-0L
 for(iter in seq_len(1000L)){
  iterations<-iter;s<-state(beta)
  if(max(abs(s$score))<1e-7){converged<-TRUE;break}
  delta<-as.numeric(solve(s$information,s$score));step<-1
  while(step>2^-40){candidate<-beta+step*delta;ns<-state(candidate);if(is.finite(ns$value)&&ns$value>=s$value-1e-12)break;step<-step/2}
  if(step<=2^-40)break
  beta<-candidate
 }
 # Numerical rescue optimizes the identical Firth objective, never another estimator.
 if(!converged){
  opt<-optim(beta,function(b)-state(b)$value,function(b)-state(b)$score,method='BFGS',control=list(maxit=2000L,reltol=1e-14))
  beta<-opt$par
  for(k in seq_len(100L)){
   s<-state(beta);if(max(abs(s$score))<1e-7)break
   delta<-as.numeric(solve(s$information,s$score));step<-1
   while(step>2^-40){candidate<-beta+step*delta;ns<-state(candidate);if(is.finite(ns$value)&&ns$value>=s$value-1e-12)break;step<-step/2}
   if(step<=2^-40)break
   beta<-candidate
  }
  converged<-max(abs(state(beta)$score))<1e-6
 }
 if(!converged||any(!is.finite(beta)))stop('FIRTH_NONCONVERGENCE')
 coefficients<-c(intercept=beta[1]-beta[2]*center/scale,slope=beta[2]/scale)
 f<-function(new_probs){z<-pmax(1e-10,pmin(1-1e-10,as.numeric(new_probs)));plogis(beta[1]+beta[2]*(z-center)/scale)}
 attr(f,'coefficients')<-coefficients;attr(f,'max_adjusted_score')<-max(abs(state(beta)$score));attr(f,'iterations')<-iterations
 attr(f,'penalized_loglik_standardized')<-state(beta)$value
 f
}
fit_selected_calibrator<-function(cal,x,y){
 if(cal=='firth')return(fit_firth_calibration(x,y))
 fn<-if(cal=='platt')'fit_platt_scaling'else'fit_isotonic_calibration'
 f<-getExportedValue('OmicSelector',fn)(x,y)
 if(cal=='platt'){
  m<-attr(f,'model')
  if(!isTRUE(m$converged))stop('PLATT_NONCONVERGENCE')
  if(isTRUE(m$boundary))stop('PLATT_BOUNDARY')
  if(m$rank<2L||any(!is.finite(coef(m))))stop('PLATT_UNIDENTIFIABLE')
 }
 f
}
