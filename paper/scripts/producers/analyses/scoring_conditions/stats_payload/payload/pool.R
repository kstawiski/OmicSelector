# benchmark_assembly derivations; all destinations are below benchmark_assembly. No model-training package is loaded.
MODEL_AUDIT <- list()
MODEL_CACHE <- new.env(parent=emptyenv())
profile_reml <- function(y,v) {
  scale <- max(stats::var(y),median(v),1e-12)
  objective <- function(u) {
    vv <- v + u*scale; w <- 1/vv; mu <- sum(w*y)/sum(w)
    -.5*(sum(log(vv))+log(sum(w))+sum(w*(y-mu)^2))
  }
  hi <- max(1,diff(range(y))^2/scale)
  while(objective(hi*2)>objective(hi) && hi<1e8)hi<-hi*2
  # REML profiles can have several local maxima. Refine every grid peak and
  # compare them with the zero boundary, rather than assuming unimodality.
  grid<-c(0,exp(seq(log(1e-12),log(hi*2),length.out=257L)))
  values<-vapply(grid,objective,0.0)
  peaks<-which(values>=c(-Inf,head(values,-1L)) & values>=c(tail(values,-1L),-Inf))
  candidates<-grid[peaks]
  for(i in peaks[peaks>1L & peaks<length(grid)]) {
   opt<-optimize(objective,c(grid[i-1L],grid[i+1L]),maximum=TRUE,tol=1e-12)
   candidates<-c(candidates,opt$maximum)
  }
  u<-candidates[which.max(vapply(candidates,objective,0.0))]
  list(tau2=u*scale,loglik=objective(u),objective=function(t)objective(t/scale))
}
robust_rma <- function(yi,sei,context='unspecified',control_sequence=NULL) {
  stopifnot(length(yi)>=2L,all(is.finite(yi)),all(is.finite(sei)&sei>0))
  key <- paste(c(format(c(yi,sei),digits=17,scientific=TRUE),deparse(control_sequence)),collapse='|')
  if(exists(key,MODEL_CACHE,inherits=FALSE)) {
    z<-get(key,MODEL_CACHE);a<-copy(z$audit);a[,context:=context]
    MODEL_AUDIT[[length(MODEL_AUDIT)+1L]]<<-a
    return(z)
  }
  prof<-profile_reml(yi,sei^2); all_warnings<-character();attempts<-0L
  controls<-list(list(),list(stepadj=.5,maxiter=10000L,threshold=1e-10),list(stepadj=.25,maxiter=20000L,threshold=1e-12),list(stepadj=.125,maxiter=40000L,threshold=1e-12))
  if(!is.null(control_sequence))controls<-control_sequence
  status<-'default_converged';reason<-character();fit<-NULL
  for(j in seq_along(controls)) {
    attempts<-j;ww<-character()
    fit<-tryCatch(withCallingHandlers(metafor::rma.uni(yi=yi,sei=sei,method='REML',test='knha',control=controls[[j]]),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
    all_warnings<-c(all_warnings,ww)
    if(inherits(fit,'error')){reason<-c(reason,conditionMessage(fit));next}
    boundary<-fit$tau2<=1e-12
    gap<-prof$loglik-prof$objective(fit$tau2)
    tau_ok<-abs(fit$tau2-prof$tau2)<=max(1e-7,1e-3*prof$tau2) && gap<=1e-6
    bad_warning<-any(grepl('converg|boundary|local maximum|set to (zero|0)',ww,ignore.case=TRUE))
    if(j==1L && (boundary||bad_warning||!tau_ok)) {
      reason<-c(reason,if(boundary)'boundary tau2'else if(!tau_ok)'profile discrepancy'else 'optimizer warning');next
    }
    if(bad_warning||!tau_ok){reason<-c(reason,if(!tau_ok)'profile discrepancy'else 'optimizer warning');next}
    status<-if(j==1L)'default_converged'else'retry_converged';break
  }
  if(inherits(fit,'error') || is.null(fit) || !is.finite(fit$tau2) || abs(fit$tau2-prof$tau2)>max(1e-7,1e-3*prof$tau2) || prof$loglik-prof$objective(fit$tau2)>1e-6 || (exists('bad_warning')&&bad_warning)) {
    # The independent REML profile is the fallback optimizer; KH uses that tau2.
    fit<-metafor::rma.uni(yi=yi,sei=sei,method='REML',test='knha',tau2=prof$tau2)
    status<-'profile_optimizer_converged'
  }
  boundary<-fit$tau2<=1e-12
  a<-data.table(context=context,k=length(yi),convergence_status=status,reml_attempts=attempts,reml_boundary=boundary,tau2=fit$tau2,tau2_profile=prof$tau2,tau2_profile_abs_difference=abs(fit$tau2-prof$tau2),profile_loglik_gap=prof$loglik-prof$objective(fit$tau2),profile_check_pass=TRUE,retry_reason=paste(unique(reason),collapse=' | '),optimizer_warnings=paste(unique(all_warnings),collapse=' | '))
  stopifnot(a$tau2_profile_abs_difference<=max(1e-7,1e-3*prof$tau2),a$profile_loglik_gap<=1e-6)
  z<-list(fit=fit,audit=a);assign(key,z,MODEL_CACHE)
  MODEL_AUDIT[[length(MODEL_AUDIT)+1L]]<<-a
  z
}
random_effects_meta <- function(d,context='unspecified',min_k=5L) {
 x<-d[is.finite(lift)&is.finite(lift_se)&lift_se>0]
 if(nrow(x)<min_k)return(data.table(k=nrow(x),pooled_lift=NA_real_,pooled_se=NA_real_,ci_lo=NA_real_,ci_hi=NA_real_,pi_lo=NA_real_,pi_hi=NA_real_,p_lift_gt_zero=NA_real_,p_lift_gt_margin_0_05=NA_real_,tau2=NA_real_,I2=NA_real_,Q=NA_real_,Q_p=NA_real_,df_knha=NA_real_,convergence_status='not_estimable_insufficient_k',reml_attempts=0L,reml_boundary=NA,tau2_profile=NA_real_,tau2_profile_abs_difference=NA_real_,profile_loglik_gap=NA_real_,profile_check_pass=NA))
 z<-robust_rma(x$lift,x$lift_se,context);f<-z$fit;pr<-predict(f);b<-as.numeric(f$b)
 a<-data.table(k=f$k,pooled_lift=b,pooled_se=f$se,ci_lo=f$ci.lb,ci_hi=f$ci.ub,pi_lo=pr$pi.lb,pi_hi=pr$pi.ub,p_lift_gt_zero=pt(b/f$se,f$ddf,lower.tail=FALSE),p_lift_gt_margin_0_05=pt((b-.05)/f$se,f$ddf,lower.tail=FALSE),tau2=f$tau2,I2=f$I2,Q=f$QE,Q_p=f$QEp,df_knha=as.numeric(f$ddf))
 cbind(a,z$audit[,.(convergence_status,reml_attempts,reml_boundary,tau2_profile,tau2_profile_abs_difference,profile_loglik_gap,profile_check_pass)])
}
