test_that("requirements cover all built-ins and preserve core routes", {
  x <- method_requirements()
  expect_equal(nrow(x), 73L)
  expect_setequal(x$method_id, singlesample_method_roster()$method_id)
  expect_true(all(x$requirements_known))
  expect_identical(method_requirements("ws-rclr-trimmed")$r_packages, "")
  expect_identical(method_requirements("reo-pairratio")$r_packages, "glmnet")
  expect_true(all(is.na(x$fit_seconds) | x$fit_seconds >= 0))
  expect_error(method_requirements("unknown"), "Unknown methods")
  expect_error(method_requirements(families = "unknown"), "Unknown families")
  expect_error(method_requirements(NA_character_), "without NA")
  expect_equal(nrow(method_requirements(character())), 0L)
  expect_setequal(method_requirements("ai-tabpfn", "A")$method_id,
                  c("ai-tabpfn", x$method_id[x$family == "A"]))
})

test_that("availability distinguishes missing from unchecked dependencies", {
  expect_true(available_methods("ws-rclr-trimmed", python = FALSE)$available)
  local_mocked_bindings(.method_r_available = function(package) FALSE)
  x <- available_methods(c("ai-tabpfn", "reo-pairratio"), python = FALSE)
  expect_false(any(x$available))
  expect_identical(x$missing_r, c("glmnet", "reticulate")[match(x$method_id, c("reo-pairratio", "ai-tabpfn"))])
  expect_match(x$install_command[1], "install_method_dependencies", fixed = TRUE)
  local_mocked_bindings(.method_r_available = function(package) TRUE)
  expect_true(is.na(available_methods("ai-tabpfn", python = FALSE)$available))
})

test_that("Python and weight failures are reported exactly", {
  local_mocked_bindings(.method_r_available = function(package) TRUE,
    .method_python_probe = function(module) if (module == "torch") "torch: import failed" else "",
    .method_weight = function(method, download = FALSE) "weights absent")
  x <- available_methods("ai-tabpfn")
  expect_false(x$available)
  expect_match(x$missing_python, "torch: import failed", fixed = TRUE)
  local_mocked_bindings(.method_python_probe = function(module) "")
  expect_match(available_methods("ai-tabpfn")$missing_weights, "weights absent")
})

test_that("dry-run installation has no Python or download side effects", {
  local_mocked_bindings(.method_python_probe = function(...) stop("must not probe"),
    .method_weight = function(...) stop("must not download"))
  x <- install_method_dependencies(c("ai-tabpfn", "tab-tabicl"), dry_run = TRUE)
  expect_true(all(c("tabpfn==8.0.8", "tabicl==2.1.1") %in% x$python_packages))
  expect_length(x$weights, 2L)
  expect_length(install_method_dependencies("ws-rclr-trimmed", dry_run = TRUE)$python_packages, 0L)
  expect_length(install_method_dependencies("ai-tabpfn", python = FALSE, dry_run = TRUE)$weights, 0L)
  expect_silent(install_method_dependencies("ws-rclr-trimmed"))
})

test_that("error hints preserve signatures, results and condition classes", {
  f <- function(x, y = 2) x + y
  g <- .method_dependency_wrapper(f, "reo-pairratio")
  expect_identical(formals(g), formals(f))
  expect_identical(g(3), f(3))
  h <- .method_dependency_wrapper(function(...) stop("package glmnet is not installed"), "reo-pairratio")
  expect_error(h(), 'install_method_dependencies\\(methods = "reo-pairratio"\\)')
  e <- structure(list(message = "bad input", call = NULL, code = "BAD"), class = c("original_error", "error", "condition"))
  h <- .method_dependency_wrapper(function(...) stop(e), "reo-pairratio")
  expect_error(h(), class = "original_error")
  expect_error(h(), "^bad input$")
})

test_that("wrappers retain caller promises and variadic arguments", {
  f <- function(x, ..., y = 4) c(x, list(...), y)
  g <- .method_dependency_wrapper(f, "reo-pairratio")
  local({
    caller_value <- 17
    expect_identical(g(caller_value, a = caller_value * 2), f(caller_value, a = caller_value * 2))
  })
  expect_identical(.method_dependency_wrapper(function(x = 3) x, "reo-pairratio")(), 3)
})
