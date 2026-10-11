.make_all_routes_data <- function() {
  set.seed(267)
  feat <- unique(c(ws_default_pivot_pool(), 'hsa-miR-451a','hsa-miR-23a-3p', paste0('F',1:48)))
  X <- matrix(rexp(80*length(feat),rate=0.1),80,dimnames=list(NULL,feat))
  y <- rep(0:1,each=40)
  X[y==1,10:15] <- X[y==1,10:15]+20
  meta <- data.frame(kit=rep(c('A','B'),40), biofluid=rep(c('serum-clot','plasma-EDTA'),each=2,length.out=80),
    technology=rep(c('NGS','Toray','qPCR','NanoString'),each=4,length.out=80),rin=seq(5,9,length.out=80),cohort='training')
  list(X=X,y=y,meta=meta)
}

test_that('all legacy routes score each specimen using frozen state and aligned metadata', {
  dat <- .make_all_routes_data()
  probe <- c(1L,6L,11L,16L,40L,80L)
  for (id in .ss_legacy_ids()) {
    dep <- suppressWarnings(deploy_singlesample(dat$X,dat$y,id,meta_train=dat$meta,verify=FALSE))
    before <- singlesample_state_fingerprint(dep)
    expect_true(is_singlesample_deployable(dep,dat$X[probe,],dat$meta[probe,]), info=id)
    expect_identical(singlesample_state_fingerprint(dep), before, info=id)
    # A changed label in prediction metadata must never enter the fitted score.
    m <- dat$meta[probe,]; m$disease <- 1
    expect_identical(score_specimen(dep,dat$X[probe,],m),score_specimen(dep,dat$X[probe,],dat$meta[probe,]),info=id)
    expect_equal(score_specimen(dep,dat$X[probe,ncol(dat$X):1],dat$meta[probe,]),
      score_specimen(dep,dat$X[probe,],dat$meta[probe,]), tolerance=0, info=id)
    for (nm in names(dep$model$metadata_schema)) {
      err <- tryCatch(score_specimen(dep,dat$X[1,,drop=FALSE]),error=identity)
      expect_identical(err$code,'MISSING_METADATA',info=id)
      m <- dat$meta[1,,drop=FALSE]; m[[nm]] <- NA
      err <- tryCatch(score_specimen(dep,dat$X[1,,drop=FALSE],m),error=identity)
      expect_identical(err$code,'MISSING_METADATA',info=id)
      m[[nm]] <- if(nm=='rin') -1 else 'never_seen'
      err <- tryCatch(score_specimen(dep,dat$X[1,,drop=FALSE],m),error=identity)
      expect_identical(err$code,if(nm=='rin')'INVALID_METADATA' else 'UNSEEN_METADATA_LEVEL',info=id)
    }
  }
})

test_that('hemolysis anchors are frozen and absent anchors fail with a code', {
  d <- .make_all_routes_data()
  h <- hemolysis_index_blondal(log(d$X+0.5))
  derived <- deploy_singlesample(d$X,d$y,'hemolysis-resid',verify=FALSE)
  supplied <- deploy_singlesample(d$X,d$y,'hemolysis-resid',meta_train=data.frame(hemo_score=h),verify=FALSE)
  expect_identical(score_specimen(derived,d$X[1:4,]),score_specimen(supplied,d$X[1:4,],data.frame(hemo_score=h[1:4])))
  X <- d$X[,setdiff(colnames(d$X),'hsa-miR-451a'),drop=FALSE]
  err <- tryCatch(score_specimen(derived,X[1,,drop=FALSE]),error=identity)
  expect_identical(err$code,'MISSING_HEMOLYSIS_ANCHORS')
  err <- tryCatch(deploy_singlesample(X,d$y,'hemolysis-resid',verify=FALSE),error=identity)
  expect_identical(err$code,'MISSING_HEMOLYSIS_ANCHORS')
  dep <- deploy_singlesample(X,d$y,'hemolysis-resid',meta_train=data.frame(hemo_score=h),verify=FALSE)
  err <- tryCatch(score_specimen(dep,X[1,,drop=FALSE]),error=identity)
  expect_identical(err$code,'MISSING_HEMOLYSIS_ANCHORS')
})

test_that('legacy scalar adapters freeze the training-median direction', {
  d <- .make_all_routes_data()
  for (id in c('ws-dominance','nc-mahalanobis','nc-conformal-anom','nc-isoforest')) {
    dep <- suppressWarnings(deploy_singlesample(d$X,d$y,id,verify=FALSE))
    X <- d$X[1:4,]; fit <- dep$model$fit
    raw <- function(X) switch(id, `ws-dominance`=ws_dominance_score(X),
      `nc-mahalanobis`=apply_compositional_mahalanobis(fit,X)$distance,
      `nc-conformal-anom`=os_conformal_anomaly(X,fit)$score,
      `nc-isoforest`=apply_isolation_forest_logratio(fit,X)$score)
    train <- raw(d$X)
    direction <- if (median(train[d$y == 1]) < median(train[d$y == 0])) -1 else 1
    expect_identical(dep$score_direction, direction, info=id)
    expect_identical(score_specimen(dep,X),as.numeric(direction * raw(X)),info=id)
  }
})

test_that('Python fingerprint observes parameters, buffers and context mutation', {
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  torch <- reticulate::import('torch',convert=FALSE)
  net <- torch$nn$BatchNorm1d(3L)
  model <- list(net=net,context=matrix(1:6,2))
  .ss_python_freeze(model)
  expect_false(reticulate::py_to_r(net$training))
  context <- .ss_python_state_module()$inference(267L)
  context$`__enter__`()
  expect_false(reticulate::py_to_r(torch$is_grad_enabled()))
  context$`__exit__`(NULL,NULL,NULL)
  old <- singlesample_state_fingerprint(model)
  expect_identical(old,singlesample_state_fingerprint(model))
  net$running_mean$add_(1)
  expect_false(identical(old,singlesample_state_fingerprint(model)))
  expect_error(singlesample_assert_row_equivariant(function(m,X,meta) {
    m$net$running_mean$add_(1); rowSums(X)
  },model,matrix(1,4,3)), 'model state changed')
  before <- singlesample_state_fingerprint(model)
  model$context[1,1] <- 99
  expect_false(identical(before,singlesample_state_fingerprint(model)))
})

test_that('importable Python methods deploy their exported fitted state', {
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  d <- .make_all_routes_data()
  dep <- deploy_singlesample(d$X,d$y,'coda-deepcoda',verify=FALSE,hp=list(epochs=2L))
  expect_identical(dep$model$hp$epochs,2L)
  expect_true(is_singlesample_deployable(dep,d$X[1:4,]))
})

test_that('prediction copies discard Python cache and buffer writes', {
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  torch <- reticulate::import('torch',convert=FALSE)
  model <- list(net=torch$nn$BatchNorm1d(2L))
  .ss_python_freeze(model)
  before <- singlesample_state_fingerprint(model)
  clone <- .ss_python_clone(model,.ss_python_state_module())
  clone$net$running_mean$add_(7)
  expect_identical(singlesample_state_fingerprint(model),before)
  expect_false(identical(singlesample_state_fingerprint(clone),before))
  np <- reticulate::import('numpy',convert=FALSE)
  p <- reticulate::import('pathlib',convert=FALSE)
  configuration <- list(dtype=np$dtype('float64'),path=p$Path('/a/frozen/checkpoint'))
  expect_identical(singlesample_state_fingerprint(configuration),singlesample_state_fingerprint(configuration))
  builtins <- reticulate::import_builtins(convert=FALSE)
  expect_error(singlesample_state_fingerprint(list(opaque=builtins$object())), 'Unsupported state types')
})

test_that('estimator reconstruction restores torch RNG on success and error', {
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  mod <- .ss_python_state_module()
  probe <- reticulate::py_run_string("
import copy
import torch

class ReconstructedEstimator:
    def __init__(self, fail=False):
        self.fail = fail

    def __deepcopy__(self, memo):
        cloned = type(self)(self.fail)
        cloned.cpu_draw = torch.rand(3)
        cloned.cuda_draws = [torch.rand(3, device=i)
                             for i in range(torch.cuda.device_count())]
        if self.fail:
            raise RuntimeError('reconstruction failed after RNG draw')
        return cloned

def check_clone_rng(module, freeze_state, fail):
    cpu = torch.random.get_rng_state()
    cuda = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else []
    try:
        if fail:
            try:
                module.clone_for_inference(ReconstructedEstimator(True), freeze_state)
            except RuntimeError as exc:
                assert str(exc) == 'reconstruction failed after RNG draw'
            else:
                raise AssertionError('copy error was suppressed')
        else:
            expected = copy.deepcopy(ReconstructedEstimator())
            torch.random.set_rng_state(cpu)
            if cuda:
                torch.cuda.set_rng_state_all(cuda)
            actual = module.clone_for_inference(ReconstructedEstimator(), freeze_state)
            assert torch.equal(expected.cpu_draw, actual.cpu_draw)
            assert all(torch.equal(a, b) for a, b in
                       zip(expected.cuda_draws, actual.cuda_draws))
        assert torch.equal(cpu, torch.random.get_rng_state())
        assert all(torch.equal(a, b) for a, b in
                   zip(cuda, torch.cuda.get_rng_state_all() if cuda else []))
        return True
    finally:
        torch.random.set_rng_state(cpu)
        if cuda:
            torch.cuda.set_rng_state_all(cuda)
", local=FALSE, convert=FALSE)
  set.seed(99L)
  before <- .Random.seed
  for (freeze_state in c(FALSE, TRUE)) {
    for (fail in c(FALSE, TRUE)) {
      expect_true(reticulate::py_to_r(probe$check_clone_rng(mod, freeze_state, fail)))
      expect_identical(.Random.seed, before)
    }
  }
})

test_that('R RNG guards restore existing and absent state on success and error', {
  for (seeded in c(FALSE, TRUE)) {
    for (fail in c(FALSE, TRUE)) {
      set.seed(99L)
      before <- .Random.seed
      if (!seeded) rm('.Random.seed', envir=globalenv())
      result <- tryCatch(.ss_with_r_rng({
        runif(1)
        if (fail) stop('RNG guard error')
        7L
      }), error=identity)
      present <- exists('.Random.seed', envir=globalenv(), inherits=FALSE)
      expect_identical(present, seeded)
      if (seeded) expect_identical(.Random.seed, before)
      if (fail) expect_match(conditionMessage(result), 'RNG guard error') else expect_identical(result, 7L)
    }
  }
  early_return <- function() .ss_with_r_rng({ runif(1); return(8L) })
  expect_identical(early_return(), 8L)
  expect_identical(.Random.seed, before)
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  net <- reticulate::import('torch',convert=FALSE)$nn$Linear(2L,2L)
  rm('.Random.seed', envir=globalenv())
  clone <- .ss_python_clone(net,.ss_python_state_module())
  expect_false(exists('.Random.seed', envir=globalenv(), inherits=FALSE))
})

test_that('deployment rejects malformed hyperparameters and unknown training fluid', {
  d <- .make_all_routes_data()
  expect_error(deploy_singlesample(d$X,d$y,hp=list(1)), 'unique non-empty names')
  expect_error(deploy_singlesample(d$X,d$y,'ws-dominance',hp=list(top_k=1)), 'legacy adapters do not accept hp')
  m <- d$meta; m$biofluid <- 'never_seen'
  err <- tryCatch(deploy_singlesample(d$X,d$y,'biofluid-stratified-rclr',meta_train=m,verify=FALSE),error=identity)
  expect_identical(err$code,'UNSEEN_METADATA_LEVEL')
})

test_that('estimator copies can preserve inference flags without sharing state', {
  skip_if_not_installed('reticulate')
  skip_if_not(reticulate::py_module_available('torch'))
  torch <- reticulate::import('torch',convert=FALSE)
  model <- list(net=torch$nn$BatchNorm1d(2L))
  before <- singlesample_state_fingerprint(model)
  mod <- .ss_python_state_module()
  preserved <- .ss_python_clone(model,mod,freeze_state=FALSE)
  frozen <- .ss_python_clone(model,mod)
  expect_true(reticulate::py_to_r(preserved$net$training))
  expect_true(reticulate::py_to_r(preserved$net$weight$requires_grad))
  expect_false(reticulate::py_to_r(frozen$net$training))
  expect_false(reticulate::py_to_r(frozen$net$weight$requires_grad))
  expect_identical(singlesample_state_fingerprint(preserved),before)
  preserved$net$running_mean$add_(7)
  expect_identical(singlesample_state_fingerprint(model),before)
  expect_false(identical(singlesample_state_fingerprint(preserved),before))
})


test_that('complete fingerprints record inference configuration without omitting fitted state', {
  m <- list(hp=list(score_batch=FALSE),context_rclr=matrix(1:6,2),context_y=c(0L,1L))
  changed <- m; changed$hp$score_batch <- TRUE
  for (fn in list(.tabpfn_model_digest,.tabicl_model_digest,.tabdpt_model_digest)) {
    expect_identical(fn(m),fn(m))
    expect_false(identical(fn(m),fn(changed)))
    mutated <- m; mutated$context_rclr[1,1] <- 9
    expect_false(identical(fn(m),fn(mutated)))
  }
})
