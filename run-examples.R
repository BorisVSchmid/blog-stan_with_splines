# Run all spline examples
# Demonstrates various capabilities of the B-spline and C-spline implementations

cat("Stan Splines - Running All Examples\n")
cat("===================================\n\n")

# Track timing
start_time <- Sys.time()

# List of examples to run
examples <- list(
  list(name = "Spline comparison plots", file = "examples/generate_plots_extended.R"),
  list(name = "B-spline smoothing comparison", file = "examples/compare_smoothing.R"),
  list(name = "Smoothing strength guide", file = "examples/smoothing_strength_guide.R"),
  list(name = "Diagnostics demonstration", file = "examples/diagnostics_demo.R"),
  list(name = "Regional splines with hierarchical priors", file = "examples/run_regional_splines.R")
)

# Run each example
for (i in seq_along(examples)) {
  example <- examples[[i]]
  
  cat(sprintf("\n[%d/%d] %s\n", i, length(examples), example$name))
  cat(rep("-", 50), "\n", sep = "")
  
  example_start <- Sys.time()
  
  tryCatch({
    source(example$file)
    cat(sprintf("\n✓ Completed in %.1f seconds\n", 
                as.numeric(Sys.time() - example_start, units = "secs")))
  }, error = function(e) {
    cat(sprintf("\n✗ ERROR: %s\n", e$message))
  })
}

# Summary
total_time <- as.numeric(Sys.time() - start_time, units = "secs")

cat("\n", rep("=", 50), "\n", sep = "")
cat(sprintf("All examples completed in %.1f seconds\n", total_time))
cat("Check output/ directory for generated plots and results\n")
cat(rep("=", 50), "\n", sep = "")