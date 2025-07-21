# Run minimal code examples
# Simple entry point for basic spline demonstrations

cat("Stan Splines - Minimal Examples\n")
cat("==============================\n\n")

# Run B-spline example
cat("--- Running B-spline example with key features ---\n")
source("code/bsplines.R")

# Run C-spline example
cat("\n--- Running C-spline minimal example ---\n")
source("code/csplines.R")

cat("\nMinimal examples complete. Check output/ directory for plots.\n")