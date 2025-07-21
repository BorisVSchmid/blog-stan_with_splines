# Main entry point for spline analyses
# Run key examples and generate outputs

cat("Stan Splines Analysis\n")
cat("====================\n\n")

# Available analyses
analyses <- c(
  "1. Compare B-splines vs C-splines (tests/test_splines.R)",
  "2. Boundary behavior analysis (examples/generate_plots_extended.R)", 
  "3. Regional splines demonstration (examples/run_regional_splines_fixed.R)",
  "4. All of the above"
)

cat("Available analyses:\n")
cat(paste(analyses, collapse = "\n"), "\n\n")

# For non-interactive use, run all analyses
choice <- 4

if (interactive()) {
  choice <- as.numeric(readline("Select analysis (1-4): "))
}

# Run selected analyses
if (choice == 1 || choice == 4) {
  cat("\n--- Running spline comparison ---\n")
  source("tests/test_splines.R")
}

if (choice == 2 || choice == 4) {
  cat("\n--- Running boundary behavior analysis ---\n")
  source("examples/generate_plots_extended.R")
}

if (choice == 3 || choice == 4) {
  cat("\n--- Running regional splines demonstration ---\n")
  source("examples/run_regional_splines_fixed.R")
}

cat("\nAnalysis complete. Check output/ directory for results.\n")