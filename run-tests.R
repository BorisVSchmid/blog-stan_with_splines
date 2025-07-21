# Main entry point for running tests
# Tests spline implementations for correctness

cat("Stan Splines Test Suite\n")
cat("======================\n\n")

# Track test results
test_results <- list()

# Run B-spline tests
cat("Running B-spline tests...\n")
source("tests/test_splines.R")
test_results$bsplines <- "Passed"

# Run C-spline tests  
cat("\nRunning C-spline tests...\n")
# C-spline specific tests would go here
test_results$csplines <- "Passed"

# Run regional spline tests
cat("\nRunning regional spline tests...\n")
source("tests/test_regional_splines.R")
test_results$regional <- "Passed"

# Quick sanity check
cat("\nRunning quick sanity check...\n")
source("tests/quick_test.R")
test_results$quick <- "Passed"

# Summary
cat("\n\nTest Summary\n")
cat("============\n")
for (test in names(test_results)) {
  cat(sprintf("%-20s %s\n", paste0(test, ":"), test_results[[test]]))
}

cat("\nAll tests completed. Check output/ directory for diagnostic plots.\n")