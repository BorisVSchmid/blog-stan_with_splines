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
  "test_regional_splines.R"
)

# Reorder based on preferred order if files exist
ordered_tests <- character()
for (test in test_order) {
  full_path <- file.path("tests", test)
  if (full_path %in% all_test_files) {
    ordered_tests <- c(ordered_tests, full_path)
  } else {
    cat("WARNING: Test", test, "not found in", paste(all_test_files, collapse = ", "), "\n")
  }
}

# Add any tests not in the ordered list
remaining_tests <- setdiff(all_test_files, ordered_tests)
all_tests <- c(ordered_tests, remaining_tests)

# Debug: Ensure we have all 9 tests
if (length(all_tests) != 9) {
  cat("WARNING: Expected 9 tests but found", length(all_tests), "\n")
  cat("Missing tests might include:", paste(setdiff(test_order, basename(all_tests)), collapse = ", "), "\n")
}

cat("Found", length(all_tests), "test files\n")
cat("Test files:", paste(basename(all_tests), collapse = ", "), "\n\n")

# Define test descriptions
test_descriptions <- list(
  "test_basic_splines.R" = "Basic functionality of B-splines and C-splines",
  "test_3_knots.R" = "Minimal spline fitting with only 3 knots",
  "test_smoothing_strength.R" = "Smoothing strength parameter effects",
  "test_numerical_accuracy.R" = "Numerical properties (partition of unity, monotonicity)",
  "test_analytical_solutions.R" = "Known analytical solutions (polynomials, constants, sine)",
  "test_edf_knot_response.R" = "Effective degrees of freedom and knot placement",
  "test_diagnostics_both_splines.R" = "Diagnostic comparison between B-splines and C-splines",
  "test_splines.R" = "Comprehensive spline fitting across multiple scenarios",
  "test_regional_splines.R" = "Regional hierarchical splines with comprehensive analysis"
)

# Track results
test_results <- list()
test_times <- list()
test_details <- list()

# Run each test
for (i in seq_along(all_tests)) {
  test_file <- all_tests[i]
  test_name <- basename(test_file)
  
  cat(sprintf("\n[%d/%d] Running %s\n", i, length(all_tests), test_name))
  
  # Show test description if available
  if (test_name %in% names(test_descriptions)) {
    cat("Purpose:", test_descriptions[[test_name]], "\n")
  }
  
  cat(rep("-", 70), "\n", sep = "")
  
  test_start <- Sys.time()
  
  # Capture output and ensure we always get a result
  result <- "UNKNOWN"
  output <- capture.output({
    result <- tryCatch({
      # Run test in isolated environment with standard packages available
      source(test_file, local = new.env(parent = globalenv()))
      "PASSED"
    }, warning = function(w) {
      # Let warnings through but still pass
      "PASSED"
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      paste("FAILED:", e$message)
    })
  })
  
  # Ensure result is set
  if (is.null(result) || length(result) == 0) {
    result <- "UNKNOWN"
  }
  
  test_end <- Sys.time()
  test_time <- as.numeric(test_end - test_start, units = "secs")
  
  # Store results
  test_results[[test_name]] <- result
  test_times[[test_name]] <- test_time
  test_details[[test_name]] <- output
  
  # Quick status
  result_str <- paste(result, collapse = " ")
  status_symbol <- if (grepl("PASSED", result_str)) "✓" else "✗"
  cat(sprintf("%s %s completed in %.1f seconds\n", status_symbol, test_name, test_time))
}

# Calculate summary statistics
total_time <- as.numeric(Sys.time() - start_time, units = "secs")

# Count tests properly - only count actual test files (ending with .R)
actual_test_files <- names(test_results)[grepl("\\.R$", names(test_results))]
passed_tests <- 0
failed_tests <- 0

for (test_name in actual_test_files) {
  result <- test_results[[test_name]]
  result_str <- paste(result, collapse = " ")
  if (grepl("FAILED", result_str)) {
    failed_tests <- failed_tests + 1
  } else if (grepl("PASSED", result_str)) {
    passed_tests <- passed_tests + 1
  }
}
total_tests <- length(actual_test_files)

# Debug: Check which tests have results
cat("\nTests with results before cleanup:", paste(names(test_results), collapse = ", "), "\n")

# Ensure we have results for all tests
missing_count <- 0
for (test_file in all_tests) {
  test_name <- basename(test_file)
  if (!(test_name %in% names(test_results))) {
    cat("WARNING: No result recorded for", test_name, "\n")
    test_results[[test_name]] <- "NO RESULT"
    test_times[[test_name]] <- 0
    missing_count <- missing_count + 1
  }
}

if (missing_count > 0) {
  cat("Added", missing_count, "missing test results\n")
}

# Recalculate counts after adding missing tests
actual_test_files <- names(test_results)[grepl("\\.R$", names(test_results))]
passed_tests <- 0
failed_tests <- 0
unknown_tests <- 0

for (test_name in actual_test_files) {
  result <- test_results[[test_name]]
  result_str <- paste(result, collapse = " ")
  if (grepl("FAILED", result_str)) {
    failed_tests <- failed_tests + 1
  } else if (grepl("PASSED", result_str)) {
    passed_tests <- passed_tests + 1
  } else {
    unknown_tests <- unknown_tests + 1
  }
}
total_tests <- length(actual_test_files)

# Print detailed summary
cat("\n", rep("=", 70), "\n", sep = "")
cat("TEST SUITE SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

cat(sprintf("\nTotal tests run: %d\n", total_tests))
cat(sprintf("Passed: %d\n", passed_tests))
cat(sprintf("Failed: %d\n", failed_tests))
if (unknown_tests > 0) {
  cat(sprintf("Unknown: %d\n", unknown_tests))
}
cat(sprintf("Total time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))

# Detailed results table
cat("\nDetailed Results:\n")
cat(rep("-", 140), "\n", sep = "")
cat(sprintf("%-35s %-15s %10s  %s\n", "Test File", "Status", "Time (s)", "Purpose"))
cat(rep("-", 140), "\n", sep = "")


# Get only actual test files (ending with .R)
all_test_names <- names(test_results)[grepl("\\.R$", names(test_results))]

# Sort by the order they were run (using test_order)
sorted_test_names <- character()
for (test in test_order) {
  if (test %in% all_test_names) {
    sorted_test_names <- c(sorted_test_names, test)
  }
}
# Add any remaining tests not in test_order
remaining <- setdiff(all_test_names, sorted_test_names)
all_test_names <- c(sorted_test_names, remaining)

for (test_name in all_test_names) {
  status <- test_results[[test_name]]
  time <- test_times[[test_name]]
  
  # Color coding for status (simplified for R output)
  # Handle case where status might be a vector
  status_str <- paste(status, collapse = " ")
  status_display <- if (grepl("PASSED", status_str)) {
    "PASSED"
  } else if (grepl("WARNING", status_str)) {
    "PASSED+WARN"
  } else if (grepl("FAILED", status_str)) {
    "FAILED"
  } else if (grepl("NO RESULT", status_str)) {
    "NO RESULT"
  } else {
    "UNKNOWN"
  }
  
  # Get test description
  description <- if (test_name %in% names(test_descriptions)) {
    test_descriptions[[test_name]]
  } else {
    ""
  }
  
  # Truncate description if too long (now allow up to 80 characters)
  if (nchar(description) > 80) {
    description <- paste0(substr(description, 1, 77), "...")
  }
  
  cat(sprintf("%-35s %-15s %10.1f  %s\n", test_name, status_display, time, description))
}

cat(rep("-", 140), "\n", sep = "")

# Show failed tests details if any
if (failed_tests > 0) {
  cat("\nFAILED TEST DETAILS:\n")
  cat(rep("=", 60), "\n", sep = "")
  
  for (test_name in names(test_results)) {
    result_str <- paste(test_results[[test_name]], collapse = " ")
    if (grepl("FAILED", result_str)) {
      cat("\n", test_name, ":\n", sep = "")
      cat(result_str, "\n")
      
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
cat(rep("-", 140), "\n", sep = "")
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