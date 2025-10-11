# Test splines against known analytical solutions
# Tests specific cases where we know the expected behavior

# Suppress default graphics device to prevent Rplots.pdf
pdf(NULL)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)  # Use base R sd, not posterior::sd

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork")
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
  
  results <- list()
  
  # Test with B-splines
  cat("  Testing B-splines...\n")
  adaptive_knots_b <- max(4, min(round(n/2), 40))
  adaptive_prior <- 2 * sd(y)
  
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_b,
    spline_degree = 3,
    smoothing_strength = 0.1,  # Default mild smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  cat("  [Test 1 - B-spline cubic polynomial fitting] Fitting model...\n")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    adapt_delta = 0.99,
    max_treedepth = 15
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  # Compare with true cubic
  y_true_plot_b <- x_plot_b^3 - 2*x_plot_b^2 + x_plot_b + 1
  rmse_b <- sqrt(mean((y_plot_b - y_true_plot_b)^2))
  
  cat("  B-spline RMSE from true cubic:", round(rmse_b, 4), "\n")
  
  # Test with C-splines (should handle cubics well)
  cat("  Testing C-splines...\n")
  adaptive_knots_c <- max(4, min(round(n/4), 20))
  
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_c
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  cat("  [Test 1 - C-spline cubic polynomial fitting] Fitting model...\n")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    max_treedepth = 15,
    adapt_delta = 0.99,
    refresh = 0
  )
  
  draws_c <- fit_c$draws(format = "matrix")
  x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  
  # Compare with true cubic
  y_true_plot_c <- x_plot_c^3 - 2*x_plot_c^2 + x_plot_c + 1
  rmse_c <- sqrt(mean((y_plot_c - y_true_plot_c)^2))
  
  cat("  C-spline RMSE from true cubic:", round(rmse_c, 4), "\n")
  
  # Should fit cubic well
  tolerance <- 0.3
  results$b_spline <- rmse_b < tolerance
  results$c_spline <- rmse_c < tolerance
  
  cat("  B-spline fits cubic well (RMSE <", tolerance, "):", results$b_spline, "\n")
  cat("  C-spline fits cubic well (RMSE <", tolerance, "):", results$c_spline, "\n")
  
  # Store fits and parameters in parent environment for plotting
  assign("fit_cubic_b", fit_b, envir = parent.frame())
  assign("fit_cubic_c", fit_c, envir = parent.frame())
  assign("cubic_b_knots", stan_data_b$num_knots, envir = parent.frame())
  assign("cubic_b_smoothing", stan_data_b$smoothing_strength, envir = parent.frame())
  assign("cubic_c_knots", stan_data_c$num_knots, envir = parent.frame())
  
  return(results)
}

# Test 2: Constant function
# Both splines should perfectly fit constant functions
test_constant_function <- function() {
  cat("\nTest 2: Constant function fitting\n")
  cat("=================================\n")
  
  n <- 12
  x <- seq(0, 10, length.out = n)
  constant_value <- 2.5
  y <- rep(constant_value, n) + rnorm(n, 0, 0.02)  # Small noise
  
  results <- list()
  
  # Use adaptive knot selection (n/2 for B-splines)
  adaptive_knots_b <- max(4, min(round(n/2), 40))
  adaptive_prior <- 2 * sd(y)
  
  # Test B-splines
  cat("  Testing B-splines...\n")
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_b,
    spline_degree = 3,
    smoothing_strength = 0.1,  # Default mild smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  cat("  [Test 2 - B-spline constant function fitting] Fitting model...\n")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    adapt_delta = 0.99,
    max_treedepth = 12,
    init = function() list(
      alpha_0 = mean(y),
      alpha_raw = rnorm(stan_data_b$num_knots + stan_data_b$spline_degree - 1, 0, 0.1),
      sigma = 0.1  # Initialize sigma away from 0
    )
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  variance_b <- stats::var(y_plot_b)
  mean_diff_b <- abs(mean(y_plot_b) - constant_value)
  
  # Use adaptive knot selection (n/4 for C-splines)
  adaptive_knots_c <- max(4, min(round(n/4), 20))
  
  # Test C-splines
  cat("  Testing C-splines...\n")
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_c
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  cat("  [Test 2 - C-spline constant function fitting] Fitting model...\n")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    adapt_delta = 0.99, 
    max_treedepth = 12,
    init = function() list(
      y_at_knots = rep(mean(y), stan_data_c$num_knots),
      sigma = 0.1  # Initialize sigma away from 0
    )
  )
  
  draws_c <- fit_c$draws(format = "matrix")
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  
  variance_c <- stats::var(y_plot_c)
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
  
  # Store fits and parameters in parent environment for plotting
  assign("fit_constant_b", fit_b, envir = parent.frame())
  assign("fit_constant_c", fit_c, envir = parent.frame())
  assign("constant_b_knots", stan_data_b$num_knots, envir = parent.frame())
  assign("constant_b_smoothing", stan_data_b$smoothing_strength, envir = parent.frame())
  assign("constant_c_knots", stan_data_c$num_knots, envir = parent.frame())
  
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
  
  # Test B-splines
  cat("  Testing B-splines...\n")
  adaptive_knots_b <- max(4, min(round(n/2), 40))
  adaptive_prior <- 2 * sd(y)
  
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_b,
    spline_degree = 3,
    smoothing_strength = 0.1,  # Default mild smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  # Check if B-spline captures the step reasonably
  before_step_b <- y_plot_b[x_plot_b < 5]
  after_step_b <- y_plot_b[x_plot_b > 5]
  
  mean_before_b <- mean(before_step_b)
  mean_after_b <- mean(after_step_b)
  step_height_b <- mean_after_b - mean_before_b
  
  cat("  B-spline - Mean before step (should be ~1):", round(mean_before_b, 2), "\n")
  cat("  B-spline - Mean after step (should be ~3):", round(mean_after_b, 2), "\n")
  cat("  B-spline - Step height (should be ~2):", round(step_height_b, 2), "\n")
  
  # Test C-splines
  cat("  Testing C-splines...\n")
  adaptive_knots_c <- max(4, min(round(n/4), 20))
  
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_c
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0
  )
  
  draws_c <- fit_c$draws(format = "matrix")
  x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  
  # Check if C-spline captures the step reasonably
  before_step_c <- y_plot_c[x_plot_c < 5]
  after_step_c <- y_plot_c[x_plot_c > 5]
  
  mean_before_c <- mean(before_step_c)
  mean_after_c <- mean(after_step_c)
  step_height_c <- mean_after_c - mean_before_c
  
  cat("  C-spline - Mean before step (should be ~1):", round(mean_before_c, 2), "\n")
  cat("  C-spline - Mean after step (should be ~3):", round(mean_after_c, 2), "\n")
  cat("  C-spline - Step height (should be ~2):", round(step_height_c, 2), "\n")
  
  # Check if step is reasonably approximated
  reasonable_before_b <- abs(mean_before_b - 1) < 0.5
  reasonable_after_b <- abs(mean_after_b - 3) < 0.5
  reasonable_step_b <- abs(step_height_b - 2) < 0.8
  
  reasonable_before_c <- abs(mean_before_c - 1) < 0.5
  reasonable_after_c <- abs(mean_after_c - 3) < 0.5
  reasonable_step_c <- abs(step_height_c - 2) < 0.8
  
  results$b_spline <- reasonable_before_b && reasonable_after_b && reasonable_step_b
  results$c_spline <- reasonable_before_c && reasonable_after_c && reasonable_step_c
  
  cat("  B-spline step function approximated well:", results$b_spline, "\n")
  cat("  C-spline step function approximated well:", results$c_spline, "\n")
  
  # Store fits and parameters in parent environment for plotting
  assign("fit_step_b", fit_b, envir = parent.frame())
  assign("fit_step_c", fit_c, envir = parent.frame())
  assign("step_b_knots", stan_data_b$num_knots, envir = parent.frame())
  assign("step_b_smoothing", stan_data_b$smoothing_strength, envir = parent.frame())
  assign("step_c_knots", stan_data_c$num_knots, envir = parent.frame())
  
  return(results)
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
  
  # Use adaptive knot selection
  adaptive_knots_b <- max(4, min(round(n/2), 40))
  adaptive_prior <- 2 * sd(y)
  
  # B-splines with sufficient knots
  cat("  Testing B-splines...\n")
  stan_data_b <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_b,  # Adaptive selection
    spline_degree = 3,
    smoothing_strength = 0.1,  # Default mild smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  draws_b <- fit_b$draws(format = "matrix")
  x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  
  # Calculate correlation with true sine
  y_true_plot_b <- sin(x_plot_b)
  correlation_b <- cor(y_plot_b, y_true_plot_b)
  
  # Use adaptive knot selection (n/4 for C-splines)
  adaptive_knots_c <- max(4, min(round(n/4), 20))
  
  # C-splines
  cat("  Testing C-splines...\n")
  stan_data_c <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots_c
  )
  
  model_c <- cmdstan_model("code/csplines.stan")
  fit_c <- model_c$sample(
    data = stan_data_c,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    adapt_delta = 0.99,
    max_treedepth = 12
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
  
  # Store fits and parameters in parent environment for plotting
  assign("fit_sine_b", fit_b, envir = parent.frame())
  assign("fit_sine_c", fit_c, envir = parent.frame())
  assign("sine_b_knots", stan_data_b$num_knots, envir = parent.frame())
  assign("sine_b_smoothing", stan_data_b$smoothing_strength, envir = parent.frame())
  assign("sine_c_knots", stan_data_c$num_knots, envir = parent.frame())
  
  return(list(b_spline = b_fits_sine, c_spline = c_fits_sine))
}

# Run all analytical tests
cat("Stan Splines - Analytical Solutions Test Suite\n")
cat("==============================================\n")

analytical_test_results <- list()

# Run tests with error handling
tryCatch({
  analytical_test_results$cubic_polynomial <- test_cubic_polynomial()
}, error = function(e) {
  cat("Error in cubic polynomial test:", e$message, "\n")
  analytical_test_results$cubic_polynomial <<- list(b_spline = FALSE, c_spline = FALSE)
})

tryCatch({
  analytical_test_results$constant_function <- test_constant_function()
}, error = function(e) {
  cat("Error in constant function test:", e$message, "\n")
  analytical_test_results$constant_function <<- list(b_spline = FALSE, c_spline = FALSE)
})

tryCatch({
  analytical_test_results$step_function <- test_step_function()
}, error = function(e) {
  cat("Error in step function test:", e$message, "\n")
  analytical_test_results$step_function <<- list(b_spline = FALSE, c_spline = FALSE)
})

tryCatch({
  analytical_test_results$sine_wave <- test_sine_wave()
}, error = function(e) {
  cat("Error in sine wave test:", e$message, "\n")
  analytical_test_results$sine_wave <<- list(b_spline = FALSE, c_spline = FALSE)
})

# Create visualization of all analytical tests
cat("\n\nCreating analytical solutions visualization...\n")

library(ggplot2)
library(patchwork)

# Store all fits for plotting
all_fits <- list()

# 1. Cubic polynomial plots (B-splines and C-splines)
# B-spline cubic
if (exists("fit_cubic_b")) {
  tryCatch({
    draws_b <- fit_cubic_b$draws(format = "matrix")
    x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b_lower <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
    y_plot_b_upper <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)
    
    # True cubic function
    y_true_plot_b <- x_plot_b^3 - 2*x_plot_b^2 + x_plot_b + 1
    
    df_cubic_b <- data.frame(
      x = x_plot_b,
      y_fitted = y_plot_b,
      y_lower = y_plot_b_lower,
      y_upper = y_plot_b_upper,
      y_true = y_true_plot_b
    )
    
    # Get data points from test
    x_data <- seq(0, 3, length.out = 15)
    # Regenerate with noise to match test
    set.seed(123)
    y_data <- x_data^3 - 2*x_data^2 + x_data + 1 + rnorm(15, 0, 0.05)
    
    p_cubic_b <- ggplot(df_cubic_b) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
      geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Cubic Polynomial: y = x³ - 2x² + x + 1 (B-splines)",
           subtitle = paste0("Knots = ", cubic_b_knots, ", smoothing = ", cubic_b_smoothing),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$cubic_b <- p_cubic_b
  }, error = function(e) {
    cat("  Warning: Could not create B-spline cubic polynomial plot:", e$message, "\n")
  })
}

# C-spline cubic
if (exists("fit_cubic_c")) {
  tryCatch({
    draws_c <- fit_cubic_c$draws(format = "matrix")
    x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c_lower <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
    y_plot_c_upper <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)
    
    # True cubic function
    y_true_plot_c <- x_plot_c^3 - 2*x_plot_c^2 + x_plot_c + 1
    
    df_cubic_c <- data.frame(
      x = x_plot_c,
      y_fitted = y_plot_c,
      y_lower = y_plot_c_lower,
      y_upper = y_plot_c_upper,
      y_true = y_true_plot_c
    )
    
    # Use same data points
    p_cubic_c <- ggplot(df_cubic_c) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
      geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Cubic Polynomial: y = x³ - 2x² + x + 1 (C-splines)",
           subtitle = paste0("Knots = ", cubic_c_knots),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$cubic_c <- p_cubic_c
  }, error = function(e) {
    cat("  Warning: Could not create C-spline cubic polynomial plot:", e$message, "\n")
  })
}

# 2. Constant function plots (B-splines and C-splines)
p_constant_b <- NULL
p_constant_c <- NULL

# B-spline constant
if (exists("fit_constant_b")) {
  tryCatch({
    draws_b <- fit_constant_b$draws(format = "matrix")
  y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
  y_plot_b_lower <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
  y_plot_b_upper <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)
  
  df_const_b <- data.frame(
    x = seq(0, 10, length.out = length(y_plot_b)),
    y_fitted = y_plot_b,
    y_lower = y_plot_b_lower,
    y_upper = y_plot_b_upper
  )
  
  # Get data points from test
  x_data <- seq(0, 10, length.out = 12)
  set.seed(123)
  y_data <- rep(2.5, 12) + rnorm(12, 0, 0.02)
  
  p_constant_b <- ggplot(df_const_b) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = 2.5, color = "black", linetype = "dashed", linewidth = 1) +
    geom_point(data = data.frame(x = x_data, y = y_data), 
               aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
    labs(title = "Constant Function: y = 2.5 (B-splines)",
         subtitle = paste0("Knots = ", constant_b_knots, ", smoothing = ", constant_b_smoothing),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8)) +
    ylim(2.2, 2.8)
  
  all_fits$constant_b <- p_constant_b
  }, error = function(e) {
    cat("  Warning: Could not create B-spline constant plot:", e$message, "\n")
  })
}

# C-spline constant
if (exists("fit_constant_c")) {
  tryCatch({
    draws_c <- fit_constant_c$draws(format = "matrix")
  y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
  y_plot_c_lower <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
  y_plot_c_upper <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)
  
  df_const_c <- data.frame(
    x = seq(0, 10, length.out = length(y_plot_c)),
    y_fitted = y_plot_c,
    y_lower = y_plot_c_lower,
    y_upper = y_plot_c_upper
  )
  
  # Use same data points
  p_constant_c <- ggplot(df_const_c) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
    geom_hline(yintercept = 2.5, color = "black", linetype = "dashed", linewidth = 1) +
    geom_point(data = data.frame(x = x_data, y = y_data), 
               aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
    labs(title = "Constant Function: y = 2.5 (C-splines)",
         subtitle = paste0("Knots = ", constant_c_knots),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8)) +
    ylim(2.2, 2.8)
  
  all_fits$constant_c <- p_constant_c
  }, error = function(e) {
    cat("  Warning: Could not create C-spline constant plot:", e$message, "\n")
  })
}

# Don't combine constant plots - let patchwork handle the layout

# 3. Step function plots (B-splines and C-splines)
# B-spline step
if (exists("fit_step_b")) {
  tryCatch({
    draws_b <- fit_step_b$draws(format = "matrix")
    x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b_lower <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
    y_plot_b_upper <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)
    
    # True step function
    y_true_step_b <- ifelse(x_plot_b < 5, 1, 3)
    
    df_step_b <- data.frame(
      x = x_plot_b,
      y_fitted = y_plot_b,
      y_lower = y_plot_b_lower,
      y_upper = y_plot_b_upper,
      y_true = y_true_step_b
    )
    
    # Get data points from test
    x_data <- seq(0, 10, length.out = 20)
    set.seed(123)
    y_data <- ifelse(x_data < 5, 1, 3) + rnorm(20, 0, 0.1)
    
    p_step_b <- ggplot(df_step_b) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
      geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
      geom_step(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Step Function Approximation (B-splines)",
           subtitle = paste0("Knots = ", step_b_knots, ", smoothing = ", step_b_smoothing),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$step_b <- p_step_b
  }, error = function(e) {
    cat("  Warning: Could not create B-spline step function plot:", e$message, "\n")
  })
}

# C-spline step
if (exists("fit_step_c")) {
  tryCatch({
    draws_c <- fit_step_c$draws(format = "matrix")
    x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c_lower <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
    y_plot_c_upper <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)
    
    # True step function
    y_true_step_c <- ifelse(x_plot_c < 5, 1, 3)
    
    df_step_c <- data.frame(
      x = x_plot_c,
      y_fitted = y_plot_c,
      y_lower = y_plot_c_lower,
      y_upper = y_plot_c_upper,
      y_true = y_true_step_c
    )
    
    # Use same data points
    p_step_c <- ggplot(df_step_c) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
      geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
      geom_step(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Step Function Approximation (C-splines)",
           subtitle = paste0("Knots = ", step_c_knots),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$step_c <- p_step_c
  }, error = function(e) {
    cat("  Warning: Could not create C-spline step function plot:", e$message, "\n")
  })
}

# 4. Sine wave plots
# B-spline sine wave
if (exists("fit_sine_b")) {
  tryCatch({
    draws_b <- fit_sine_b$draws(format = "matrix")
    x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
    y_plot_b_lower <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
    y_plot_b_upper <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)
    
    df_sine_b <- data.frame(
      x = x_plot_b,
      y_fitted = y_plot_b,
      y_lower = y_plot_b_lower,
      y_upper = y_plot_b_upper,
      y_true = sin(x_plot_b)
    )
    
    # Get data points from test
    x_data <- seq(0, 4*pi, length.out = 25)
    set.seed(123)
    y_data <- sin(x_data) + rnorm(25, 0, 0.1)
    
    p_sine_b <- ggplot(df_sine_b) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
      geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Sine Wave: y = sin(x) (B-splines)",
           subtitle = paste0("Knots = ", sine_b_knots, ", smoothing = ", sine_b_smoothing),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$sine_b <- p_sine_b
  }, error = function(e) {
    cat("  Warning: Could not create B-spline sine wave plot:", e$message, "\n")
  })
}

# C-spline sine wave
if (exists("fit_sine_c")) {
  tryCatch({
    draws_c <- fit_sine_c$draws(format = "matrix")
    x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
    y_plot_c_lower <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
    y_plot_c_upper <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)
    
    df_sine_c <- data.frame(
      x = x_plot_c,
      y_fitted = y_plot_c,
      y_lower = y_plot_c_lower,
      y_upper = y_plot_c_upper,
      y_true = sin(x_plot_c)
    )
    
    # Use same data points
    p_sine_c <- ggplot(df_sine_c) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
      geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      geom_point(data = data.frame(x = x_data, y = y_data), 
                 aes(x = x, y = y), color = "black", size = 2, alpha = 0.7) +
      labs(title = "Sine Wave: y = sin(x) (C-splines)",
           subtitle = paste0("Knots = ", sine_c_knots),
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    all_fits$sine_c <- p_sine_c
  }, error = function(e) {
    cat("  Warning: Could not create C-spline sine wave plot:", e$message, "\n")
  })
}

# Combine all plots - always try to create output
cat("\nPlots available for combining:\n")
cat("  Number of plots:", length(all_fits), "\n")
cat("  Plot names:", paste(names(all_fits), collapse = ", "), "\n\n")

if (length(all_fits) > 0) {
  combined_plot <- wrap_plots(all_fits, ncol = 2) +
      plot_annotation(
        title = "Analytical Solutions Test: How Well Do Splines Handle Known Functions?",
        subtitle = "Black dashed lines show true functions, colored lines show spline fits (B-splines=blue, C-splines=green)",
        caption = "Note: All tests use adaptive knot selection (B-splines: n/2 knots, C-splines: n/4 knots) with default smoothing=0.1 for B-splines"
      )
  
  # Save the plot
  dir.create("output", showWarnings = FALSE)
  ggsave("output/test-analytical_solutions.png", combined_plot, 
         width = 10, height = 8, dpi = 300)
  cat("Saved analytical solutions plot to output/test-analytical_solutions.png\n")
} else {
  # Create a minimal plot if no tests succeeded
  cat("No successful tests to plot, creating placeholder image...\n")
  empty_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "All analytical solution tests failed\nCheck error messages above", 
             size = 8, hjust = 0.5, vjust = 0.5) +
    theme_void()
  
  dir.create("output", showWarnings = FALSE)
  ggsave("output/test-analytical_solutions.png", empty_plot, 
         width = 10, height = 8, dpi = 300)
  cat("Saved placeholder plot to output/test-analytical_solutions.png\n")
}

# Close the null device at the end
dev.off()

# Summary
cat("\n\nAnalytical Solutions Test Summary\n")
cat("=================================\n")
cat("Cubic polynomial (B-splines):", ifelse(analytical_test_results$cubic_polynomial$b_spline, "PASS", "FAIL"), "\n")
cat("Cubic polynomial (C-splines):", ifelse(analytical_test_results$cubic_polynomial$c_spline, "PASS", "FAIL"), "\n")
cat("Constant function (B-splines):", ifelse(analytical_test_results$constant_function$b_spline, "PASS", "FAIL"), "\n")
cat("Constant function (C-splines):", ifelse(analytical_test_results$constant_function$c_spline, "PASS", "FAIL"), "\n")
cat("Step function (B-splines):", ifelse(analytical_test_results$step_function$b_spline, "PASS", "FAIL"), "\n")
cat("Step function (C-splines):", ifelse(analytical_test_results$step_function$c_spline, "PASS", "FAIL"), "\n")
cat("Sine wave (B-splines):", ifelse(analytical_test_results$sine_wave$b_spline, "PASS", "FAIL"), "\n")
cat("Sine wave (C-splines):", ifelse(analytical_test_results$sine_wave$c_spline, "PASS", "FAIL"), "\n")

# Overall result
all_passed <- analytical_test_results$cubic_polynomial$b_spline &&
              analytical_test_results$cubic_polynomial$c_spline &&
              analytical_test_results$constant_function$b_spline &&
              analytical_test_results$constant_function$c_spline &&
              analytical_test_results$step_function$b_spline &&
              analytical_test_results$step_function$c_spline &&
              analytical_test_results$sine_wave$b_spline &&
              analytical_test_results$sine_wave$c_spline

cat("\nOverall Result:", ifelse(all_passed, "ALL TESTS PASSED", "SOME TESTS FAILED"), "\n")

if (!all_passed) {
  cat("\nNote: Some tests may fail due to the challenging nature of the functions\n")
  cat("or the stochastic nature of MCMC sampling.\n")
}
