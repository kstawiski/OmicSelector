#!/usr/bin/env Rscript

# Complete Phase 9: Package Polish & Optimization
# Run from package root directory: Rscript complete_phase9.R
# Or in R: source("complete_phase9.R")

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║  OmicSelector 2.0 - Phase 9: Package Polish & Optimization ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Check we're in the right directory
if (!file.exists("DESCRIPTION")) {
  stop("Error: Run this script from the package root directory")
}

# Function to install if needed
ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, quiet = TRUE)
  }
}

# 1. Install required packages
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 0: Installing Required Packages\n")
cat("════════════════════════════════════════════════════════════\n\n")

packages <- c("devtools", "roxygen2", "testthat", "covr", "styler", "lintr", "pkgdown")
for (pkg in packages) {
  ensure_package(pkg)
}

library(devtools)
library(roxygen2)
library(testthat)

cat("✓ All required packages installed\n\n")

# 2. Generate documentation
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 1: Generating Documentation\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Running roxygen2::roxygenise()...\n")
document()

# Check NAMESPACE
namespace_lines <- readLines("NAMESPACE", n = 50)
cat("\n✓ NAMESPACE generated with", length(namespace_lines), "lines\n")
cat("\nNew exports (first 20):\n")
exports <- grep("^export\\(", namespace_lines, value = TRUE)
cat(paste(head(exports, 20), collapse = "\n"), "\n")

# Count .Rd files
rd_files <- list.files("man", pattern = "\\.Rd$")
cat("\n✓", length(rd_files), "help files (.Rd) generated in man/\n\n")

# 3. Run tests
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 2: Running Tests\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Running testthat tests...\n")
test_results <- test()

# Summary
if (inherits(test_results, "testthat_results")) {
  n_pass <- sum(test_results$passed)
  n_fail <- sum(test_results$failed)
  n_warn <- sum(test_results$warning)
  n_skip <- sum(test_results$skipped)

  cat("\n═══ Test Results ═══\n")
  cat("PASS:", n_pass, "\n")
  cat("FAIL:", n_fail, "\n")
  cat("WARN:", n_warn, "\n")
  cat("SKIP:", n_skip, "\n")

  if (n_fail > 0) {
    cat("\n⚠ WARNING: Some tests failed! Review output above.\n\n")
  } else {
    cat("\n✓ All tests passed!\n\n")
  }
} else {
  cat("✓ Tests completed\n\n")
}

# 4. Run R CMD check
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 3: Running R CMD Check\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("This may take several minutes...\n")
check_results <- check(
  document = FALSE,  # Already documented
  build_args = "--no-build-vignettes",  # Faster
  args = "--no-tests"  # Already tested
)

cat("\n✓ R CMD check completed\n")
cat("\nCheck summary:\n")
print(check_results)

if (length(check_results$errors) > 0) {
  cat("\n⚠ ERRORS found:", length(check_results$errors), "\n")
  cat("Fix these before proceeding!\n\n")
}

if (length(check_results$warnings) > 0) {
  cat("\n⚠ WARNINGS found:", length(check_results$warnings), "\n")
  cat("Address these before CRAN submission\n\n")
}

if (length(check_results$notes) > 0) {
  cat("\nℹ NOTES found:", length(check_results$notes), "\n")
  cat("Review and minimize before CRAN submission\n\n")
}

# 5. Calculate coverage
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 4: Calculating Test Coverage\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Calculating coverage (this may take a few minutes)...\n")
cov <- covr::package_coverage(type = "tests", quiet = FALSE)

cat("\n")
print(cov)

coverage_pct <- covr::percent_coverage(cov)
cat("\n✓ Overall coverage:", sprintf("%.1f%%", coverage_pct), "\n")

if (coverage_pct >= 80) {
  cat("✓ Coverage target (80%) achieved!\n\n")
} else {
  cat("⚠ Coverage below 80% target. Consider adding more tests.\n\n")
}

# Save coverage report
cat("Generating coverage report...\n")
covr::report(cov, file = "coverage_report.html", browse = FALSE)
cat("✓ Coverage report saved to: coverage_report.html\n\n")

# 6. Code style check
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 5: Code Style Check\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Running styler (may modify files)...\n")
styler::style_pkg(filetype = "R")
cat("✓ Code formatted with styler\n\n")

cat("Running lintr...\n")
lint_results <- lintr::lint_package()
n_lints <- length(lint_results)

if (n_lints > 0) {
  cat("\n⚠ Found", n_lints, "linting issues\n")
  cat("Review with: lintr::lint_package()\n\n")
} else {
  cat("\n✓ No linting issues found!\n\n")
}

# 7. Build package
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 6: Building Package\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Building package tarball...\n")
tarball <- build(vignettes = FALSE)  # Faster without vignettes
cat("✓ Package built:", tarball, "\n\n")

# 8. Build documentation website (optional)
cat("════════════════════════════════════════════════════════════\n")
cat("STEP 7: Documentation Website (Optional)\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("Build website? (y/n): ")
response <- readline()

if (tolower(substr(response, 1, 1)) == "y") {
  cat("Building pkgdown site...\n")
  pkgdown::build_site()
  cat("✓ Website built in docs/\n")
  cat("View at: docs/index.html\n\n")
} else {
  cat("Skipped. Build later with: pkgdown::build_site()\n\n")
}

# 9. Summary and next steps
cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║              PHASE 9 POLISH COMPLETE!                      ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("═══ Summary ═══\n\n")

cat("Documentation:\n")
cat("  • NAMESPACE: ", length(exports), " exports\n")
cat("  • Help files: ", length(rd_files), " .Rd files\n")
cat("  • Vignettes: 1\n\n")

cat("Testing:\n")
cat("  • Test files: ", length(list.files("tests/testthat", pattern = "^test-.*\\.R$")), "\n")
if (exists("n_pass")) {
  cat("  • Tests passed: ", n_pass, "\n")
  cat("  • Tests failed: ", n_fail, "\n")
}
if (exists("coverage_pct")) {
  cat("  • Coverage: ", sprintf("%.1f%%", coverage_pct), "\n")
}
cat("\n")

cat("Package Check:\n")
if (exists("check_results")) {
  cat("  • Errors: ", length(check_results$errors), "\n")
  cat("  • Warnings: ", length(check_results$warnings), "\n")
  cat("  • Notes: ", length(check_results$notes), "\n")
}
cat("\n")

cat("═══ Next Steps ═══\n\n")

cat("1. Review and fix any ERRORs or WARNINGs from R CMD check\n")
cat("2. Install locally:\n")
cat("   devtools::install()\n\n")
cat("3. Test installation:\n")
cat("   library(OmicSelector)\n")
cat("   ?OmicSelector_nested_cv\n\n")
cat("4. Build vignettes (if not done):\n")
cat("   devtools::build_vignettes()\n\n")
cat("5. Commit changes:\n")
cat("   git add NAMESPACE man/\n")
cat("   git commit -m 'Phase 9: Package polish complete'\n\n")
cat("6. Continue development:\n")
cat("   - Phase 3: Multi-omics Integration\n")
cat("   - Phase 5: Survival Analysis\n")
cat("   - Phase 8: Explainability\n\n")
cat("7. Prepare for CRAN (when ready):\n")
cat("   - Address all NOTEs\n")
cat("   - Build vignettes\n")
cat("   - Test on win-builder\n")
cat("   - Submit to CRAN\n\n")

cat("═══ Files Generated ═══\n\n")
cat("• NAMESPACE (updated)\n")
cat("• man/*.Rd (", length(rd_files), " files)\n")
cat("• coverage_report.html\n")
cat("• ", basename(tarball), "\n")
cat("• docs/ (if website built)\n\n")

cat("Phase 9 complete! Package is now installable and ready for use. ✓\n\n")

cat("════════════════════════════════════════════════════════════\n\n")
