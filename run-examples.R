# Run extended examples and demonstrations
# More complex analyses showing spline capabilities

cat("Stan Splines - Extended Examples\n")
cat("================================\n\n")

# Available examples
examples <- c(
  "1. Spline comparison plots (examples/generate_plots_extended.R)",
  "2. Regional splines with hierarchical priors (examples/run_regional_splines.R)", 
  "3. B-spline smoothing comparison (examples/compare_smoothing.R)",
  "4. All examples"
)

cat("Available examples:\n")
cat(paste(examples, collapse = "\n"), "\n\n")

# For non-interactive use, run all examples
choice <- 4

if (interactive()) {
  choice <- as.numeric(readline("Select example (1-4): "))
}

# Run selected examples
if (choice == 1 || choice == 4) {
  cat("\n--- Running spline comparison plots ---\n")
  source("examples/generate_plots_extended.R")
}

if (choice == 2 || choice == 4) {
  cat("\n--- Running regional splines demonstration ---\n")
  source("examples/run_regional_splines.R")
}

if (choice == 3 || choice == 4) {
  cat("\n--- Running B-spline smoothing comparison ---\n")
  source("examples/compare_smoothing.R")
}

cat("\nExamples complete. Check output/ directory for results.\n")