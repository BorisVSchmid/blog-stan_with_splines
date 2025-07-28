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
  
  # Use adaptive knot selection (n/4 for C-splines)
  adaptive_knots <- max(4, min(round(n/4), 20))
  
  # Test with C-splines (should handle cubics well)
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots
  )
  
  model <- cmdstan_model("code/csplines.stan")
  cat("  [Test 1 - C-spline cubic polynomial fitting] Fitting model...\n")
  fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
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
  
  # Store fit and parameters in parent environment for plotting
  assign("fit_cubic", fit, envir = parent.frame())
  assign("cubic_knots", stan_data$num_knots, envir = parent.frame())
  
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
  y <- rep(constant_value, n) + rnorm(n, 0, 0.05)  # Increased noise to avoid sigma approaching 0
  
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
    smoothing_strength = 2,  # Default smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  cat("  [Test 2 - B-spline constant function fitting] Fitting model...\n")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 600,
    refresh = 0,
    adapt_delta = 0.98,  # Even higher for constant function with low noise
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
    iter_warmup = 200,
    iter_sampling = 400,
    refresh = 0,
    adapt_delta = 0.95,  # Higher to handle low noise case
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
  
  # Use adaptive knot selection (n/4 for C-splines)
  adaptive_knots <- max(4, min(round(n/4), 20))
  
  # C-splines
  cat("  Testing C-splines...\n")
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = adaptive_knots
  )
  
  model <- cmdstan_model("code/csplines.stan")
  fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
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
  
  # Store fit and parameters in parent environment for plotting
  assign("fit_step", fit, envir = parent.frame())
  assign("step_knots", stan_data$num_knots, envir = parent.frame())
  
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
    smoothing_strength = 2,  # Default smoothing
    prior_scale = adaptive_prior
  )
  
  model_b <- cmdstan_model("code/bsplines.stan")
  fit_b <- model_b$sample(
    data = stan_data_b,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 700,
    refresh = 0,
    adapt_delta = 0.95
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
    iter_sampling = 700,
    refresh = 0,
    adapt_delta = 0.95
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
  analytical_test_results$cubic_polynomial <<- FALSE
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
  analytical_test_results$step_function <<- FALSE
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

# 1. Cubic polynomial plot
if (exists("fit_cubic")) {
  tryCatch({
    draws <- fit_cubic$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  # True cubic function
  y_true_plot <- x_plot^3 - 2*x_plot^2 + x_plot + 1
  
  df_cubic <- data.frame(
    x = x_plot,
    y_fitted = y_plot,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper,
    y_true = y_true_plot
  )
  
  p_cubic <- ggplot(df_cubic) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
    geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
    labs(title = "Cubic Polynomial: y = x³ - 2x² + x + 1",
         subtitle = paste0("C-splines (knots = ", cubic_knots, ") should fit cubic polynomials perfectly"),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
    all_fits$cubic <- p_cubic
  }, error = function(e) {
    cat("  Warning: Could not create cubic polynomial plot:", e$message, "\n")
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
  
  p_constant_b <- ggplot(df_const_b) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = 2.5, color = "black", linetype = "dashed", linewidth = 1) +
    labs(title = "Constant Function: y = 2.5 (B-splines)",
         subtitle = paste0("Knots = ", constant_b_knots, ", smoothing = ", constant_b_smoothing, " (should be flat)"),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
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
  
  p_constant_c <- ggplot(df_const_c) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
    geom_hline(yintercept = 2.5, color = "black", linetype = "dashed", linewidth = 1) +
    labs(title = "Constant Function: y = 2.5 (C-splines)",
         subtitle = paste0("Knots = ", constant_c_knots, " (should be flat)"),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  all_fits$constant_c <- p_constant_c
  }, error = function(e) {
    cat("  Warning: Could not create C-spline constant plot:", e$message, "\n")
  })
}

# Don't combine constant plots - let patchwork handle the layout

# 3. Step function plot
if (exists("fit_step")) {
  tryCatch({
    draws <- fit_step$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  # True step function
  y_true_step <- ifelse(x_plot < 5, 1, 3)
  
  df_step <- data.frame(
    x = x_plot,
    y_fitted = y_plot,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper,
    y_true = y_true_step
  )
  
  p_step <- ggplot(df_step) +
    geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
    geom_step(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
    labs(title = "Step Function Approximation (C-splines)",
         subtitle = paste0("Knots = ", step_knots, " (creates smooth transition)"),
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  all_fits$step <- p_step
  }, error = function(e) {
    cat("  Warning: Could not create step function plot:", e$message, "\n")
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
    
    p_sine_b <- ggplot(df_sine_b) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
      geom_line(aes(x = x, y = y_fitted), color = "blue", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
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
    
    p_sine_c <- ggplot(df_sine_c) +
      geom_ribbon(aes(x = x, ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "darkgreen") +
      geom_line(aes(x = x, y = y_fitted), color = "darkgreen", linewidth = 1.2) +
      geom_line(aes(x = x, y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
      labs(title = "Sine Wave: y = sin(x) (C-splines)",
           subtitle = paste0("Knots = ", sine_c_knots, " (smooth periodic function)"),
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
      subtitle = "Black dashed lines show true functions, colored lines show spline fits"
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
cat("Cubic polynomial (C-splines):", ifelse(analytical_test_results$cubic_polynomial, "PASS", "FAIL"), "\n")
cat("Constant function (B-splines):", ifelse(analytical_test_results$constant_function$b_spline, "PASS", "FAIL"), "\n")
cat("Constant function (C-splines):", ifelse(analytical_test_results$constant_function$c_spline, "PASS", "FAIL"), "\n")
cat("Step function approximation:", ifelse(analytical_test_results$step_function, "PASS", "FAIL"), "\n")
cat("Sine wave (B-splines):", ifelse(analytical_test_results$sine_wave$b_spline, "PASS", "FAIL"), "\n")
cat("Sine wave (C-splines):", ifelse(analytical_test_results$sine_wave$c_spline, "PASS", "FAIL"), "\n")

# Overall result
all_passed <- analytical_test_results$cubic_polynomial &&
              analytical_test_results$constant_function$b_spline &&
              analytical_test_results$constant_function$c_spline &&
              analytical_test_results$step_function &&
              analytical_test_results$sine_wave$b_spline &&
              analytical_test_results$sine_wave$c_spline

cat("\nOverall Result:", ifelse(all_passed, "ALL TESTS PASSED", "SOME TESTS FAILED"), "\n")

if (!all_passed) {
  cat("\nNote: Some tests may fail due to the challenging nature of the functions\n")
  cat("or the stochastic nature of MCMC sampling.\n")
}