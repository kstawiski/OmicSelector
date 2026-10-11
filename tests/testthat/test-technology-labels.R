.technology_labels_data <- function() {
  set.seed(267)
  features <- c(.ta_default_pivot_pool(), paste0('F', 1:12))
  X <- matrix(rexp(80 * length(features)) + 1, 80,
              dimnames = list(NULL, features))
  y <- rep(0:1, each = 40)
  X[y == 1, 7:9] <- X[y == 1, 7:9] + 3
  list(X = X, y = y, meta = data.frame(
    technology = rep(c('small-RNA-seq', 'OpenArray qPCR'), 40)))
}

test_that('all transfer cache technology levels have explicit outcomes', {
  # Every distinct sample_meta() technology in cohort_cache_r3.rds (34 cohorts).
  labels <- c('small-RNA-seq', 'qPCR', 'Exiqon qPCR', 'TaqMan TLDA qPCR',
              'OpenArray qPCR', 'TLDA qPCR', 'TaqMan OpenArray qPCR',
              '3D-Gene miRNA microarray', 'HTG EdgeSeq miRNA WTA counts')
  expected <- c('NGS', rep('qPCR', 6), 'Toray', 'NGS')
  expect_identical(.ta_clean_technology(labels), expected)
  expect_identical(.ta_clean_technology(paste0(' ', toupper(labels), ' ')), expected)
  d <- .technology_labels_data()
  for (label in c('five-platform repeated assay', 'microarray',
                  'Affymetrix miRNA microarray', 'unknown assay',
                  'unvalidated qPCR kit', 'novel small-RNA-seq')) {
    expect_error(.ta_clean_technology(label), 'unsupported values', fixed = TRUE)
    m <- d$meta; m$technology[1] <- label
    err <- tryCatch(deploy_singlesample(d$X, d$y, 'tech-stratified-rclr',
                                       meta_train = m, verify = FALSE), error = identity)
    expect_s3_class(err, 'singlesample_input_error')
    expect_identical(err$code, 'UNSEEN_METADATA_LEVEL')
  }
})

test_that('documented assay families and legacy aliases keep their classes', {
  qpcr <- c('TaqMan', 'TaqMan qPCR', 'TLDA', 'OpenArray', 'TaqMan OpenArray',
            'TaqMan TLDA', 'Exiqon', 'LNA qPCR', 'Exiqon LNA qPCR',
            'Exiqon/LNA qPCR', 'SYBR', 'SYBR qPCR', 'SYBR Green qPCR',
            'relative-quantity qPCR', 'relative quantity qPCR', 'ddCt qPCR',
            'delta-delta-Ct qPCR', 'delta delta Ct qPCR', '\u0394\u0394Ct qPCR',
            'q-pcr', 'rt-pcr', 'rt-qPCR')
  ngs <- c('NGS', 'miRNA-seq', 'miRNA_seq', 'small_rna_seq', 'small RNA-seq',
           'small RNA sequencing', 'small-RNA sequencing', 'miRNA sequencing',
           'miRNA-sequencing', 'microRNA sequencing', 'HTG EdgeSeq',
           'HTG EdgeSeq counts', 'HTG EdgeSeq miRNA counts')
  toray <- c('Toray', 'array-toray3d', 'array_toray3d', 'toray-3dgene',
             '3D-Gene', '3D-Gene array', '3D-Gene arrays', 'Toray 3D-Gene',
             'Toray 3D-Gene miRNA microarray')
  nano <- c('NanoString', 'nano-string', 'array-other')
  for (pair in list(list(qpcr, 'qPCR'), list(ngs, 'NGS'), list(toray, 'Toray'),
                   list(nano, 'NanoString'))) {
    expect_identical(.ta_clean_technology(pair[[1]]), rep(pair[[2]], length(pair[[1]])))
  }
  expect_identical(.ta_clean_technology(c(NA, '', 'unknown')), rep('unknown', 3))
})

test_that('default anchors require known classes present only in training', {
  d <- .technology_labels_data(); panel <- colnames(d$X)[7:9]
  a <- select_cross_tech_anchor_features(d$X, d$meta)
  b <- select_cross_tech_anchor_features(d$X, d$meta,
                                         required_technologies = c('NGS', 'qPCR'))
  expect_identical(a, b)
  expect_true(all(a$selected_yes == 'yes'))
  fit <- fit_cross_tech_harmonized_rclr(d$X, d$meta, panel)
  expect_identical(fit$required_technologies, c('NGS', 'qPCR'))
  expect_identical(fit$anchor_table, b)
  probe <- d$meta[1:3, , drop = FALSE]; probe$technology <- 'Toray'
  before <- singlesample_state_fingerprint(fit)
  expect_equal(score_cross_tech_harmonized_rclr(d$X[1:3, ], probe, panel,
                train_expr = d$X, train_meta = d$meta),
               predict_cross_tech_harmonized_rclr(fit, d$X[1:3, ], probe), tolerance = 0)
  expect_identical(before, singlesample_state_fingerprint(fit))
  # Input order and unused factor levels cannot create required classes.
  m <- d$meta; m$technology <- factor(m$technology,
    levels = c('Toray', 'OpenArray qPCR', 'small-RNA-seq'))
  expect_identical(select_cross_tech_anchor_features(d$X[80:1, ], m[80:1, , drop=FALSE]), a)
})

test_that('default refuses fewer than two classes while explicit vectors retain meaning', {
  d <- .technology_labels_data(); panel <- colnames(d$X)[7:9]
  code <- function(expr) tryCatch(force(expr), error = function(e) e$code)
  for (labels in list(rep('qPCR', 80), rep(c('TLDA qPCR', 'OpenArray qPCR'), 40),
                      rep(c('NGS', 'unknown'), 40), rep('unknown', 80))) {
    m <- d$meta; m$technology <- labels
    expect_identical(code(select_cross_tech_anchor_features(d$X, m)),
                     'INSUFFICIENT_TECHNOLOGY_CLASSES')
    expect_identical(code(fit_cross_tech_harmonized_rclr(d$X, m, panel)),
                     'INSUFFICIENT_TECHNOLOGY_CLASSES')
  }
  m <- d$meta; m$technology <- 'NGS'
  expect_identical(code(deploy_singlesample(d$X, d$y, 'cross-tech-harmonized-rclr',
                    meta_train = m, verify = FALSE)), 'INSUFFICIENT_TECHNOLOGY_CLASSES')
  expect_true(all(select_cross_tech_anchor_features(d$X, m,
    required_technologies = 'NGS')$selected_yes == 'yes'))
  explicit <- c('NGS', 'Toray', 'qPCR', 'NanoString')
  fit <- fit_cross_tech_harmonized_rclr(d$X, d$meta, panel, required_technologies=explicit)
  expect_identical(fit$required_technologies, explicit)
  expect_length(fit$anchor_features, 0)
  expect_identical(code(predict_cross_tech_harmonized_rclr(fit, d$X, d$meta)),
                   'INSUFFICIENT_CROSS_TECH_ANCHORS')
})

test_that('deployment accepts raw aliases and freezes the training-derived default', {
  d <- .technology_labels_data()
  d$meta$kit_family <- rep(c('kit-A', 'kit-B'), 40)
  d$meta$biofluid <- 'plasma extracellular vesicles'
  original <- d$meta
  for (id in c('tech-stratified-rclr', 'cross-tech-harmonized-rclr', 'tech-residualized-alr')) {
    dep <- deploy_singlesample(d$X, d$y, id, meta_train = d$meta, verify = FALSE)
    expect_identical(dep$model$metadata_schema$technology, c('NGS', 'qPCR'))
    if (id == 'cross-tech-harmonized-rclr')
      expect_identical(dep$model$fit$required_technologies, c('NGS', 'qPCR'))
    before <- singlesample_state_fingerprint(dep)
    m <- d$meta[1:5, ]; m$technology <- c('miRNA sequencing', 'TLDA qPCR',
                                        'HTG EdgeSeq miRNA WTA counts', 'Exiqon qPCR', 'NGS')
    batch <- score_specimen(dep, d$X[1:5, ], m)
    one <- vapply(1:5, function(i)as.numeric(score_specimen(dep, d$X[i,,drop=FALSE],
                                                          m[i,,drop=FALSE])), numeric(1))
    expect_true(all(is.finite(batch)))
    expect_equal(as.numeric(batch), one, tolerance = 1e-10)
    expect_identical(before, singlesample_state_fingerprint(dep))
    m$technology[1] <- 'Toray'
    err <- tryCatch(score_specimen(dep, d$X[1:5, ], m), error = identity)
    expect_identical(err$code, 'UNSEEN_METADATA_LEVEL')
    expect_identical(d$meta, original)
  }
})
