# Simple test script without groundhog
# Load required packages
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(123)

# Create output directory if it doesn't exist
dir.create("output", showWarnings = FALSE)

# Source the functions from the main test script
source("test_splines.R", local = TRUE)

# Generate test data
data_sine <- generate_test_data(n = 30, func_type = "sine", noise_sd = 0.1)

# Test B-splines
cat("Testing B-spline model...\n")
fit_b <- fit_bspline(data_sine, num_knots = 7, spline_degree = 3, 
                     chains = 2, iter_warmup = 500, iter_sampling = 500)

# Test C-splines  
cat("\nTesting C-spline model...\n")
fit_c <- fit_cspline(data_sine, num_knots = 7,
                     chains = 2, iter_warmup = 500, iter_sampling = 500)

# Create plots
cat("\nCreating plots...\n")
plots_b <- plot_spline_fit(fit_b, data_sine, "B-spline", show_basis = TRUE)
plots_c <- plot_spline_fit(fit_c, data_sine, "C-spline")

# Save plots
if (!is.null(plots_b$fit) && !is.null(plots_c$fit)) {
  comparison_plot <- plots_b$fit | plots_c$fit
  ggsave("output/spline_comparison_simple.png", comparison_plot, 
         width = 12, height = 6, dpi = 300)
  cat("Saved comparison plot to output/spline_comparison_simple.png\n")
  
  if (!is.null(plots_b$basis)) {
    ggsave("output/bspline_basis_simple.png", plots_b$basis,
           width = 8, height = 6, dpi = 300)
    cat("Saved basis functions to output/bspline_basis_simple.png\n")
  }
}

cat("\nDone!\n")