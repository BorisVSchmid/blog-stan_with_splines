# Test smoothing diagnostics with both B-splines and C-splines

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)
source("code/smoothing_diagnostics.R")

# Generate test data
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x
y <- y_true + rnorm(n, 0, 0.15)

cat("Testing Smoothing Diagnostics for Both Spline Types\n")
cat("===================================================\n\n")

# Test 1: B-splines with various settings
cat("Test 1: B-splines with normal smoothing\n")
cat("---------------------------------------\n")

# Prepare B-spline data
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 8,
  spline_degree = 3,
  smoothing_strength = 0,
  prior_scale = 2 * sd(y)
)

# Compile and fit B-spline model
model_b <- cmdstan_model("code/bsplines.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

# Run diagnostics
diagnosis_b <- diagnose_smoothing(fit_b, x, y, stan_data_b, "bspline")
cat("\nB-spline Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_b)

# Test 2: C-splines
cat("\n\nTest 2: C-splines\n")
cat("-----------------\n")

# Prepare C-spline data (no prior_scale or tau_smooth)
stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 8
)

# Compile and fit C-spline model
model_c <- cmdstan_model("code/csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

# Run diagnostics
diagnosis_c <- diagnose_smoothing(fit_c, x, y, stan_data_c, "cspline")
cat("\nC-spline Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_c)

# Test 3: B-splines with overfitting (many knots)
cat("\n\nTest 3: B-splines with potential overfitting (many knots)\n")
cat("----------------------------------------------------------\n")

stan_data_b_overfit <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 20,  # Many knots for n=40
  spline_degree = 3,
  smoothing_strength = 0,
  prior_scale = 5 * sd(y)  # Large prior scale
)

fit_b_overfit <- model_b$sample(
  data = stan_data_b_overfit,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

diagnosis_b_overfit <- diagnose_smoothing(fit_b_overfit, x, y, stan_data_b_overfit, "bspline")
cat("\nB-spline (overfitting) Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_b_overfit)

# Test 4: C-splines with few knots (over-smoothing)
cat("\n\nTest 4: C-splines with potential over-smoothing (few knots)\n")
cat("------------------------------------------------------------\n")

stan_data_c_smooth <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 3  # Very few knots
)

fit_c_smooth <- model_c$sample(
  data = stan_data_c_smooth,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

diagnosis_c_smooth <- diagnose_smoothing(fit_c_smooth, x, y, stan_data_c_smooth, "cspline")
cat("\nC-spline (over-smoothing) Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_c_smooth)

# Summary
cat("\n\nSummary of Diagnostic Differences\n")
cat("=================================\n")
cat("B-splines have access to:\n")
cat("  - prior_scale adjustments\n")
cat("  - tau_smooth for random walk smoothing\n")
cat("  - num_knots control\n\n")
cat("C-splines only have:\n")
cat("  - num_knots control\n")
cat("  - Natural boundary conditions (built-in)\n\n")
cat("The diagnostics now provide appropriate advice for each spline type.\n")