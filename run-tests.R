# Comprehensive test runner for all spline tests
# Runs all test suites and provides summary

cat("Stan Splines - Comprehensive Test Suite\n")
cat("======================================\n")
cat("Running all tests...\n\n")

# Record start time
start_time <- Sys.time()

# List of test scripts to run
test_scripts <- c(
  "tests/quick_test.R",
  "tests/test_numerical_accuracy.R", 
  "tests/test_analytical_solutions.R"
)

# Track results
all_results <- list()
test_times <- list()

# Run each test script
for (i in seq_along(test_scripts)) {
  script <- test_scripts[i]
  script_name <- basename(script)
  
  cat("=======================================================\n")
  cat("Running:", script_name, "\n")
  cat("=======================================================\n")
  
  script_start <- Sys.time()
  
  tryCatch({
    source(script, local = TRUE)
    all_results[[script_name]] <- "COMPLETED"
  }, error = function(e) {
    cat("ERROR in", script_name, ":", e$message, "\n")
    all_results[[script_name]] <- paste("ERROR:", e$message)
  })
  
  script_end <- Sys.time()
  test_times[[script_name]] <- as.numeric(script_end - script_start, units = "secs")
  
  cat("\n", script_name, "completed in", round(test_times[[script_name]], 1), "seconds\n")
  cat("\n")
}

# Final summary
end_time <- Sys.time()
total_time <- as.numeric(end_time - start_time, units = "secs")

cat("=======================================================\n")
cat("FINAL TEST SUMMARY\n")
cat("=======================================================\n")

for (i in seq_along(test_scripts)) {
  script_name <- basename(test_scripts[i])
  status <- all_results[[script_name]]
  time <- test_times[[script_name]]
  
  cat(sprintf("%-30s: %s (%.1fs)\n", script_name, status, time))
}

cat("\n")
cat("Total testing time:", round(total_time, 1), "seconds\n")

# Check if all completed successfully
all_completed <- all(sapply(all_results, function(x) x == "COMPLETED"))

cat("Overall result:", ifelse(all_completed, "ALL TEST SUITES COMPLETED", "SOME TEST SUITES HAD ERRORS"), "\n")

if (!all_completed) {
  cat("\nFailed test suites:\n")
  failed_tests <- names(all_results)[all_results != "COMPLETED"]
  for (test in failed_tests) {
    cat("  -", test, ":", all_results[[test]], "\n")
  }
}

cat("\n=======================================================\n")
cat("Test suite documentation:\n")
cat("  quick_test.R             - Basic functionality checks\n")
cat("  test_numerical_accuracy.R - Mathematical properties tests\n")
cat("  test_analytical_solutions.R - Known function fitting tests\n")
cat("  test_splines.R           - Comprehensive function type tests\n")
cat("  test_regional_splines.R  - Hierarchical model tests\n")
cat("=======================================================\n")