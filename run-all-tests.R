# Comprehensive test runner for ALL spline tests
# Runs every test in the tests/ directory

cat("Stan Splines - Complete Test Suite\n")
cat("==================================\n")
cat("Running ALL tests...\n\n")

# Record start time
start_time <- Sys.time()

# Get all test scripts
all_test_files <- list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE)

# Order tests by complexity/speed (fast tests first)
test_order <- c(
  "test_basic_splines.R",
  "test_3_knots.R",
  "test_smoothing_strength.R",
  "test_numerical_accuracy.R",
  "test_analytical_solutions.R",
  "test_edf_knot_response.R",
  "test_diagnostics_both_splines.R",
  "test_splines.R",
  "test_regional_splines_simple.R",
  "test_regional_splines.R"
)

# Reorder based on preferred order if files exist
ordered_tests <- character()
for (test in test_order) {
  full_path <- file.path("tests", test)
  if (full_path %in% all_test_files) {
    ordered_tests <- c(ordered_tests, full_path)
  }
}

# Add any tests not in the ordered list
remaining_tests <- setdiff(all_test_files, ordered_tests)
all_tests <- c(ordered_tests, remaining_tests)

cat("Found", length(all_tests), "test files\n\n")

# Track results
test_results <- list()
test_times <- list()
test_details <- list()

# Run each test
for (i in seq_along(all_tests)) {
  test_file <- all_tests[i]
  test_name <- basename(test_file)
  
  cat(sprintf("\n[%d/%d] Running %s\n", i, length(all_tests), test_name))
  cat(rep("-", 50), "\n", sep = "")
  
  test_start <- Sys.time()
  
  # Capture output
  output <- capture.output({
    result <- tryCatch({
      source(test_file, local = TRUE)
      "PASSED"
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      paste("FAILED:", e$message)
    }, warning = function(w) {
      cat("WARNING:", w$message, "\n")
      "PASSED WITH WARNINGS"
    })
  })
  
  test_end <- Sys.time()
  test_time <- as.numeric(test_end - test_start, units = "secs")
  
  # Store results
  test_results[[test_name]] <- result
  test_times[[test_name]] <- test_time
  test_details[[test_name]] <- output
  
  # Quick status
  status_symbol <- if (grepl("PASSED", result)) "✓" else "✗"
  cat(sprintf("%s %s completed in %.1f seconds\n", status_symbol, test_name, test_time))
}

# Calculate summary statistics
total_time <- as.numeric(Sys.time() - start_time, units = "secs")
passed_tests <- sum(grepl("PASSED", unlist(test_results)))
failed_tests <- sum(grepl("FAILED", unlist(test_results)))
total_tests <- length(test_results)

# Print detailed summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("TEST SUITE SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

cat(sprintf("\nTotal tests run: %d\n", total_tests))
cat(sprintf("Passed: %d\n", passed_tests))
cat(sprintf("Failed: %d\n", failed_tests))
cat(sprintf("Total time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))

# Detailed results table
cat("\nDetailed Results:\n")
cat(rep("-", 60), "\n", sep = "")
cat(sprintf("%-35s %-15s %10s\n", "Test File", "Status", "Time (s)"))
cat(rep("-", 60), "\n", sep = "")

for (test_name in names(test_results)) {
  status <- test_results[[test_name]]
  time <- test_times[[test_name]]
  
  # Color coding for status (simplified for R output)
  status_display <- if (grepl("PASSED", status)) {
    "PASSED"
  } else if (grepl("WARNING", status)) {
    "PASSED+WARN"
  } else {
    "FAILED"
  }
  
  cat(sprintf("%-35s %-15s %10.1f\n", test_name, status_display, time))
}

cat(rep("-", 60), "\n", sep = "")

# Show failed tests details if any
if (failed_tests > 0) {
  cat("\nFAILED TEST DETAILS:\n")
  cat(rep("=", 60), "\n", sep = "")
  
  for (test_name in names(test_results)) {
    if (grepl("FAILED", test_results[[test_name]])) {
      cat("\n", test_name, ":\n", sep = "")
      cat(test_results[[test_name]], "\n")
      
      # Show last few lines of output
      output <- test_details[[test_name]]
      if (length(output) > 0) {
        relevant_output <- tail(output[output != ""], 10)
        if (length(relevant_output) > 0) {
          cat("Last output lines:\n")
          cat(paste("  ", relevant_output), sep = "\n")
        }
      }
    }
  }
}

# Test categories summary
cat("\n\nTest Categories:\n")
cat(rep("-", 60), "\n", sep = "")
cat("Core functionality:  test_basic_splines.R, test_3_knots.R\n")
cat("Mathematical tests:  test_numerical_accuracy.R, test_analytical_solutions.R\n")
cat("Parameter tests:     test_smoothing_strength.R, test_edf_knot_response.R\n")
cat("Diagnostic tests:    test_diagnostics_both_splines.R\n")
cat("Advanced tests:      test_splines.R, test_regional_splines*.R\n")

# Final status
cat("\n", rep("=", 60), "\n", sep = "")
if (failed_tests == 0) {
  cat("✓ ALL TESTS PASSED!\n")
} else {
  cat("✗ SOME TESTS FAILED - Please review the output above\n")
}
cat(rep("=", 60), "\n\n", sep = "")

# Return summary invisibly
invisible(list(
  total = total_tests,
  passed = passed_tests,
  failed = failed_tests,
  time = total_time,
  results = test_results
))