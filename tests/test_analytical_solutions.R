# Test splines against known analytical solutions
# Tests specific cases where we know the expected behavior

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::var)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

set.seed(123)

cat("Testing splines against analytical solutions...\n")

# Test 1: Polynomial fitting
# Cubic splines should perfectly fit cubic polynomials
test_cubic_polynomial <- function() {
  cat("\nTest 1: Cubic polynomial fitting\n")
  cat("================================\n")
  
  # Generate cubic polynomial data: y = x^3 - 2x^2 + x + 1
  n <- 15
  x <- seq(0, 3, length.out = n)
  y_true <- x^3 - 2*x^2 + x + 1
  y <- y_true + rnorm(n, 0, 0.05)  # Small amount of noise
  
  # Test with C-splines (should handle cubics well)
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 8
  )
  
  model <- cmdstan_model("code/csplines.stan")
  fit <- model$sample(
    data = stan_data,
    chains = 1,
    iter_warmup = 300,
    iter_sampling = 600,
    refresh = 0
  )
  
  draws <- fit$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  
  # Compare with true cubic
  y_true_plot <- x_plot^3 - 2*x_plot^2 + x_plot + 1
  rmse <- sqrt(mean((y_plot - y_true_plot)^2))
  
  cat("  RMSE from true cubic function:", round(rmse, 4), "\n")
  
  # Should fit cubic well
  tolerance <- 0.3
  fits_well <- rmse < tolerance
  
  cat("  Fits cubic well (RMSE <", tolerance, "):", fits_well, "\n")
  
  return(fits_well)
}

# Test 2: Constant function
# Both splines should perfectly fit constant functions
test_constant_function <- function() {
  cat("\nTest 2: Constant function fitting\n")
  cat("=================================\n")
  
  n <- 12
  x <- seq(0, 10, length.out = n)
  constant_value <- 2.5
  y <- rep(constant_value, n) + rnorm(n, 0, 0.02)
  
  results <- list()
  
  # Test B-splines
  cat("  Testing B-splines...\n")
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 6,
    spline_degree = 3,
    smoothing_strength = 0,
    prior_scale = 1
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 400,
    refresh = 0
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  variance_b <- var(y_plot_b)
  mean_diff_b <- abs(mean(y_plot_b) - constant_value)
  
  # Test C-splines
  cat("  Testing C-splines...\n")
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 6
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 400,
    refresh = 0
  )
  
  draws_c <- fit_c$draws(format = "matrix")
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  
  variance_c <- var(y_plot_c)
  mean_diff_c <- abs(mean(y_plot_c) - constant_value)
  
  cat("  B-spline variance:", round(variance_b, 4), "\n")
  cat("  B-spline mean difference:", round(mean_diff_b, 4), "\n")
  cat("  C-spline variance:", round(variance_c, 4), "\n")
  cat("  C-spline mean difference:", round(mean_diff_c, 4), "\n")
  
  # Should have low variance and close mean
  var_tolerance <- 0.01
  mean_tolerance <- 0.1
  
  b_good_variance <- variance_b < var_tolerance
  b_good_mean <- mean_diff_b < mean_tolerance
  c_good_variance <- variance_c < var_tolerance
  c_good_mean <- mean_diff_c < mean_tolerance
  
  results$b_spline <- b_good_variance && b_good_mean
  results$c_spline <- c_good_variance && c_good_mean
  
  cat("  B-spline fits constant well:", results$b_spline, "\n")
  cat("  C-spline fits constant well:", results$c_spline, "\n")
  
  return(results)
}

# Test 3: Step function approximation
# Test how well splines approximate discontinuous functions
test_step_function <- function() {
  cat("\nTest 3: Step function approximation\n")
  cat("===================================\n")
  
  # Create step function data
  n <- 20
  x <- seq(0, 10, length.out = n)
  y_true <- ifelse(x < 5, 1, 3)  # Step at x = 5
  y <- y_true + rnorm(n, 0, 0.1)
  
  # Test with both splines
  results <- list()
  
  # C-splines
  cat("  Testing C-splines...\n")
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 8
  )
  
  model <- cmdstan_model("code/csplines.stan")
  fit <- model$sample(
    data = stan_data,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 400,
    refresh = 0
  )
  
  draws <- fit$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  
  # Check if spline captures the step reasonably
  # Values before step should be closer to 1, after step closer to 3
  before_step <- y_plot[x_plot < 5]
  after_step <- y_plot[x_plot > 5]
  
  mean_before <- mean(before_step)
  mean_after <- mean(after_step)
  step_height <- mean_after - mean_before
  
  cat("  Mean value before step (should be ~1):", round(mean_before, 2), "\n")
  cat("  Mean value after step (should be ~3):", round(mean_after, 2), "\n")
  cat("  Step height (should be ~2):", round(step_height, 2), "\n")
  
  # Check if step is reasonably approximated
  reasonable_before <- abs(mean_before - 1) < 0.5
  reasonable_after <- abs(mean_after - 3) < 0.5
  reasonable_step <- abs(step_height - 2) < 0.8
  
  step_approximated <- reasonable_before && reasonable_after && reasonable_step
  
  cat("  Step function reasonably approximated:", step_approximated, "\n")
  
  return(step_approximated)
}

# Test 4: Sine wave fitting
# Test smooth periodic function fitting
test_sine_wave <- function() {
  cat("\nTest 4: Sine wave fitting\n")
  cat("=========================\n")
  
  # Generate sine wave data
  n <- 25
  x <- seq(0, 4*pi, length.out = n)
  y_true <- sin(x)
  y <- y_true + rnorm(n, 0, 0.1)
  
  # Test with both splines
  results <- list()
  
  # B-splines with sufficient knots
  cat("  Testing B-splines...\n")
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 12,  # Need enough knots for oscillating function
    spline_degree = 3,
    smoothing_strength = 1,  # Mild smoothing for sine wave
    prior_scale = 2
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 1,
    iter_warmup = 300,
    iter_sampling = 500,
    refresh = 0
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  # Calculate correlation with true sine
  y_true_plot_b <- sin(x_plot_b)
  correlation_b <- cor(y_plot_b, y_true_plot_b)
  
  # C-splines
  cat("  Testing C-splines...\n")
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 12
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 1,
    iter_warmup = 300,
    iter_sampling = 500,
    refresh = 0
  )
  
  draws_c <- fit_c$draws(format = "matrix")
  x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  
  y_true_plot_c <- sin(x_plot_c)
  correlation_c <- cor(y_plot_c, y_true_plot_c)
  
  cat("  B-spline correlation with true sine:", round(correlation_b, 3), "\n")
  cat("  C-spline correlation with true sine:", round(correlation_c, 3), "\n")
  
  # Both should have high correlation with sine wave
  min_correlation <- 0.8
  b_fits_sine <- correlation_b > min_correlation
  c_fits_sine <- correlation_c > min_correlation
  
  cat("  B-spline fits sine well (corr >", min_correlation, "):", b_fits_sine, "\n")
  cat("  C-spline fits sine well (corr >", min_correlation, "):", c_fits_sine, "\n")
  
  return(list(b_spline = b_fits_sine, c_spline = c_fits_sine))
}

# Run all analytical tests
cat("Stan Splines - Analytical Solutions Test Suite\n")
cat("==============================================\n")

test_results <- list()

# Run tests with error handling
tryCatch({
  test_results$cubic_polynomial <- test_cubic_polynomial()
}, error = function(e) {
  cat("Error in cubic polynomial test:", e$message, "\n")
  test_results$cubic_polynomial <<- FALSE
})

tryCatch({
  test_results$constant_function <- test_constant_function()
}, error = function(e) {
  cat("Error in constant function test:", e$message, "\n")
  test_results$constant_function <<- list(b_spline = FALSE, c_spline = FALSE)
})

tryCatch({
  test_results$step_function <- test_step_function()
}, error = function(e) {
  cat("Error in step function test:", e$message, "\n")
  test_results$step_function <<- FALSE
})

tryCatch({
  test_results$sine_wave <- test_sine_wave()
}, error = function(e) {
  cat("Error in sine wave test:", e$message, "\n")
  test_results$sine_wave <<- list(b_spline = FALSE, c_spline = FALSE)
})

# Summary
cat("\n\nAnalytical Solutions Test Summary\n")
cat("=================================\n")
cat("Cubic polynomial (C-splines):", ifelse(test_results$cubic_polynomial, "PASS", "FAIL"), "\n")
cat("Constant function (B-splines):", ifelse(test_results$constant_function$b_spline, "PASS", "FAIL"), "\n")
cat("Constant function (C-splines):", ifelse(test_results$constant_function$c_spline, "PASS", "FAIL"), "\n")
cat("Step function approximation:", ifelse(test_results$step_function, "PASS", "FAIL"), "\n")
cat("Sine wave (B-splines):", ifelse(test_results$sine_wave$b_spline, "PASS", "FAIL"), "\n")
cat("Sine wave (C-splines):", ifelse(test_results$sine_wave$c_spline, "PASS", "FAIL"), "\n")

# Overall result
all_passed <- test_results$cubic_polynomial &&
              test_results$constant_function$b_spline &&
              test_results$constant_function$c_spline &&
              test_results$step_function &&
              test_results$sine_wave$b_spline &&
              test_results$sine_wave$c_spline

cat("\nOverall Result:", ifelse(all_passed, "ALL TESTS PASSED", "SOME TESTS FAILED"), "\n")

if (!all_passed) {
  cat("\nNote: Some tests may fail due to the challenging nature of the functions\n")
  cat("or the stochastic nature of MCMC sampling.\n")
}