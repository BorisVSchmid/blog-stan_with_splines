# Run hierarchical regional splines example
# Demonstrates how multiple regions can share information through hierarchical priors

cat("Stan Splines - Running Hierarchical Example\n")
cat("==========================================\n\n")

# Track timing
start_time <- Sys.time()

# List of examples to run
examples <- list(
  list(name = "Hierarchical regional splines", file = "examples/hierarchical_regional_splines.R")
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