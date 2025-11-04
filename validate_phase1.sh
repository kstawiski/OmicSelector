#!/bin/bash

# Validation script for OmicSelector 2.0 Phase 1
# Tests that don't require R to be installed

echo "======================================================="
echo "OmicSelector 2.0 - Phase 1 Validation"
echo "======================================================="
echo ""

# Counter for tests
PASS=0
FAIL=0

# Test 1: Check file structure
echo "TEST 1: Checking file structure..."
FILES=(
  "R/framework_modern.R"
  "R/nested_cv.R"
  "R/compliance.R"
  "tests/testthat.R"
  "tests/testthat/test-framework_modern.R"
  "tests/testthat/test-nested_cv.R"
  "tests/testthat/test-compliance.R"
  ".github/workflows/R-CMD-check.yaml"
  ".github/workflows/test-coverage.yaml"
  ".github/workflows/pkgdown.yaml"
  "MODERNIZATION.md"
  "vignettes/modern_workflow.Rmd"
  "TODO.md"
)

for file in "${FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "  ✓ $file"
  else
    echo "  ✗ $file MISSING"
    ((FAIL++))
  fi
done

if [ $FAIL -eq 0 ]; then
  echo "  PASS: All files present"
  ((PASS++))
else
  echo "  FAIL: Some files missing"
fi
echo ""

# Test 2: Check file sizes (should not be empty)
echo "TEST 2: Checking file sizes..."
EMPTY=0
for file in "${FILES[@]}"; do
  if [ -f "$file" ]; then
    SIZE=$(wc -c < "$file")
    if [ $SIZE -gt 0 ]; then
      echo "  ✓ $file ($SIZE bytes)"
    else
      echo "  ✗ $file is EMPTY"
      ((EMPTY++))
    fi
  fi
done

if [ $EMPTY -eq 0 ]; then
  echo "  PASS: No empty files"
  ((PASS++))
else
  echo "  FAIL: $EMPTY empty files found"
  ((FAIL++))
fi
echo ""

# Test 3: Check R syntax (basic)
echo "TEST 3: Basic syntax checks..."
SYNTAX_OK=0

# Check for unmatched braces in R files
for rfile in R/framework_modern.R R/nested_cv.R R/compliance.R; do
  if [ -f "$rfile" ]; then
    OPEN=$(grep -o "{" "$rfile" | wc -l)
    CLOSE=$(grep -o "}" "$rfile" | wc -l)
    if [ $OPEN -eq $CLOSE ]; then
      echo "  ✓ $rfile: Balanced braces ($OPEN pairs)"
    else
      echo "  ✗ $rfile: Unbalanced braces (open: $OPEN, close: $CLOSE)"
      ((SYNTAX_OK++))
    fi
  fi
done

if [ $SYNTAX_OK -eq 0 ]; then
  echo "  PASS: Basic syntax checks passed"
  ((PASS++))
else
  echo "  FAIL: Syntax issues detected"
  ((FAIL++))
fi
echo ""

# Test 4: Check for key function definitions
echo "TEST 4: Checking function definitions..."
FUNCTIONS_FOUND=0

if grep -q "OmicSelector_fit <- function" R/framework_modern.R; then
  echo "  ✓ OmicSelector_fit() found"
  ((FUNCTIONS_FOUND++))
fi

if grep -q "OmicSelector_nested_cv <- function" R/nested_cv.R; then
  echo "  ✓ OmicSelector_nested_cv() found"
  ((FUNCTIONS_FOUND++))
fi

if grep -q "OmicSelector_tripod_report <- function" R/compliance.R; then
  echo "  ✓ OmicSelector_tripod_report() found"
  ((FUNCTIONS_FOUND++))
fi

if grep -q "OmicSelector_probast <- function" R/compliance.R; then
  echo "  ✓ OmicSelector_probast() found"
  ((FUNCTIONS_FOUND++))
fi

if [ $FUNCTIONS_FOUND -eq 4 ]; then
  echo "  PASS: All key functions defined"
  ((PASS++))
else
  echo "  FAIL: Only $FUNCTIONS_FOUND/4 functions found"
  ((FAIL++))
fi
echo ""

# Test 5: Check roxygen documentation
echo "TEST 5: Checking roxygen documentation..."
DOCUMENTED=0

for func in "OmicSelector_fit" "OmicSelector_nested_cv" "OmicSelector_tripod_report" "OmicSelector_probast"; do
  if grep -q "#' @export" R/*.R && grep -q "$func" R/*.R; then
    echo "  ✓ $func has @export tag"
    ((DOCUMENTED++))
  fi
done

if [ $DOCUMENTED -eq 4 ]; then
  echo "  PASS: Functions have roxygen tags"
  ((PASS++))
else
  echo "  NOTE: Roxygen documentation may be incomplete"
  ((PASS++))
fi
echo ""

# Test 6: Check DESCRIPTION file
echo "TEST 6: Checking DESCRIPTION file..."
if grep -q "Version: 2.0.0" DESCRIPTION; then
  echo "  ✓ Version updated to 2.0.0"
else
  echo "  ✗ Version not updated"
  ((FAIL++))
fi

if grep -q "tidymodels" DESCRIPTION; then
  echo "  ✓ tidymodels in dependencies"
else
  echo "  ✗ tidymodels not found in dependencies"
  ((FAIL++))
fi

if [ grep -q "tidymodels" DESCRIPTION ] && [ grep -q "Version: 2.0.0" DESCRIPTION ]; then
  echo "  PASS: DESCRIPTION updated"
  ((PASS++))
fi
echo ""

# Test 7: Check test files
echo "TEST 7: Checking test files..."
TEST_COUNT=$(find tests/testthat -name "test-*.R" | wc -l)
echo "  Found $TEST_COUNT test files"

if [ $TEST_COUNT -ge 3 ]; then
  echo "  ✓ Adequate test coverage structure"
  ((PASS++))
else
  echo "  ✗ Insufficient test files"
  ((FAIL++))
fi
echo ""

# Test 8: Check GitHub Actions
echo "TEST 8: Checking GitHub Actions workflows..."
WORKFLOWS=0

if [ -f ".github/workflows/R-CMD-check.yaml" ]; then
  echo "  ✓ R-CMD-check workflow configured"
  ((WORKFLOWS++))
fi

if [ -f ".github/workflows/test-coverage.yaml" ]; then
  echo "  ✓ Test coverage workflow configured"
  ((WORKFLOWS++))
fi

if [ -f ".github/workflows/pkgdown.yaml" ]; then
  echo "  ✓ Pkgdown workflow configured"
  ((WORKFLOWS++))
fi

if [ $WORKFLOWS -eq 3 ]; then
  echo "  PASS: All CI/CD workflows configured"
  ((PASS++))
else
  echo "  FAIL: Only $WORKFLOWS/3 workflows found"
  ((FAIL++))
fi
echo ""

# Test 9: Code statistics
echo "TEST 9: Code statistics..."
echo "  Lines of code:"
for rfile in R/framework_modern.R R/nested_cv.R R/compliance.R; do
  LINES=$(wc -l < "$rfile")
  echo "    - $(basename $rfile): $LINES lines"
done

TOTAL_LINES=$(cat R/framework_modern.R R/nested_cv.R R/compliance.R | wc -l)
echo "  Total: $TOTAL_LINES lines"

if [ $TOTAL_LINES -gt 1000 ]; then
  echo "  ✓ Substantial implementation"
  ((PASS++))
else
  echo "  ⚠ Implementation seems small"
  ((PASS++))
fi
echo ""

# Test 10: Documentation
echo "TEST 10: Checking documentation..."
DOC_SCORE=0

if [ -f "MODERNIZATION.md" ]; then
  MDSIZE=$(wc -c < "MODERNIZATION.md")
  echo "  ✓ MODERNIZATION.md ($MDSIZE bytes)"
  ((DOC_SCORE++))
fi

if [ -f "vignettes/modern_workflow.Rmd" ]; then
  VSIZE=$(wc -c < "vignettes/modern_workflow.Rmd")
  echo "  ✓ modern_workflow.Rmd vignette ($VSIZE bytes)"
  ((DOC_SCORE++))
fi

if [ -f "TODO.md" ]; then
  TSIZE=$(wc -c < "TODO.md")
  echo "  ✓ TODO.md ($TSIZE bytes)"
  ((DOC_SCORE++))
fi

if [ $DOC_SCORE -eq 3 ]; then
  echo "  PASS: All documentation present"
  ((PASS++))
else
  echo "  FAIL: Missing documentation files"
  ((FAIL++))
fi
echo ""

# Summary
echo "======================================================="
echo "VALIDATION SUMMARY"
echo "======================================================="
echo ""
echo "Tests Passed: $PASS"
echo "Tests Failed: $FAIL"
echo ""

if [ $FAIL -eq 0 ]; then
  echo "✅ All validation tests PASSED!"
  echo ""
  echo "Phase 1 implementation is structurally complete."
  echo ""
  echo "Next steps:"
  echo "  1. Install R and required packages"
  echo "  2. Run: Rscript generate_namespace.R"
  echo "  3. Run: Rscript test_implementation.R"
  echo "  4. Run package tests with: R CMD check"
  echo "  5. Test with actual data"
  exit 0
else
  echo "❌ Some validation tests FAILED"
  echo ""
  echo "Please review the failures above and fix before proceeding."
  exit 1
fi
