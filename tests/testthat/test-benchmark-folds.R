test_that('benchmark folds retain first-row group order and RNG semantics', {
  y <- c(1,0,1,0,1,0,1,0,1,0,1,0)
  groups <- c('z','b','z','d','x','f','w','h','v','j','u','l')
  f <- os_make_benchmark_folds(y,groups,k=3L,seed=101L)
  expect_equal(sort(unlist(f)),seq_along(y))
  expect_true(any(vapply(f,function(x)all(c(1L,3L)%in%x),logical(1))))
  after <- .Random.seed
  expect_identical(f,os_make_benchmark_folds(y,groups,k=3L,seed=101L))
  expect_identical(.Random.seed,after)
  expect_error(os_make_benchmark_folds(y,groups[-1]),'length')
  expect_error(os_make_benchmark_folds(c(0,2),c('a','b')),'binary')
  expect_error(os_make_benchmark_folds(y,groups,k=2.5),'integer')
  # Existing row-count balancing helper remains an independent interface.
  expect_s3_class(os_make_grouped_stratified_folds(y,groups,n_folds=3),'os_grouped_folds')
})
