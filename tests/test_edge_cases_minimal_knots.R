# Test edge cases with minimal knots for both B-splines and C-splines
# Tests:
# 1. B-splines with 3 knots (minimum that works)
# 2. C-splines with 3 knots (minimum that works)
# 3. Error handling for insufficient knots (2 knots)
# 4. Behavior with very few basis functions

# Suppress default graphics device to prevent Rplots.pdf
pdf(NULL)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)
source("code/smoothing_diagnostics.R")

# Generate test data with complex function
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x 
y <- y_true + rnorm(n, 0, 0.15)

cat("Testing Edge Cases with Minimal Knots\n")
cat("=====================================\n\n")

# Test 1: Testing minimum knot configurations
cat("Test 1: Testing minimum knot configurations\n")
cat("-----------------------------------------------------------------\n")

# B-splines with 2 knots (minimum allowed)
stan_data_2knots <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 2,  # Minimum allowed for B-splines
  spline_degree = 3,
  smoothing_strength = 0.1,
  prior_scale = 2 * sd(y)
)

model_b <- cmdstan_model("code/bsplines.stan")

cat("Testing B-splines with 2 knots (minimum allowed)... ")
fit_b_2knots <- model_b$sample(
  data = stan_data_2knots,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)
cat("SUCCESS: B-splines accepts 2 knots (boundary knots only)\n")
cat("Number of basis functions: 2 + 3 - 1 = 4\n")

# Print diagnostic summary
cat("\nB-spline (2 knots) MCMC Diagnostic Summary:\n")
print(fit_b_2knots$diagnostic_summary())

# Run smoothing diagnostics for 2-knot model
cat("\n")
diagnosis_b_2knots_early <- diagnose_smoothing(fit_b_2knots, x, y, stan_data_2knots, "bspline")
print_smoothing_diagnostics(diagnosis_b_2knots_early)

# C-splines with 2 knots (should fail - requires at least 3)
stan_data_c_2knots <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 2  # Should fail - C-splines require at least 3
)

model_c <- cmdstan_model("code/csplines.stan")

cat("\nTesting C-splines with 2 knots (below minimum)... ")
fit_c_2knots <- NULL
suppressWarnings({
  tryCatch({
    fit_c_2knots <- model_c$sample(
      data = stan_data_c_2knots,
      chains = 1,
      iter_warmup = 100,
      iter_sampling = 100,
      refresh = 0,
      show_messages = FALSE
    )
  }, error = function(e) {
    # Stan parameter constraint errors might not always be caught here
  })
})

if (is.null(fit_c_2knots) || !exists("fit_c_2knots")) {
  cat("Expected: C-splines requires at least 3 knots\n")
} else {
  # Check if the fit actually has valid samples
  tryCatch({
    fit_c_2knots$summary()
    cat("UNEXPECTED: C-splines with 2 knots succeeded!\n")
  }, error = function(e) {
    cat("Expected: C-splines requires at least 3 knots (failed during sampling)\n")
  })
}

# Test 2: B-splines with 3 knots
cat("\n\nTest 2: B-splines with 3 knots (low flexibility configuration)\n")
cat("-----------------------------------------------------------------\n")

stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 3,  # One interior knot plus boundaries
  spline_degree = 3,
  smoothing_strength = 0.1,
  prior_scale = 2 * sd(y)
)

cat("Number of basis functions: 3 + 3 - 1 = 5\n")
cat("Fitting B-spline model...\n")

fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Print diagnostic summary
cat("\nB-spline MCMC Diagnostic Summary:\n")
print(fit_b$diagnostic_summary())

# Run smoothing diagnostics
cat("\n")
diagnosis_b <- diagnose_smoothing(fit_b, x, y, stan_data_b, "bspline")
print_smoothing_diagnostics(diagnosis_b)

# Test 3: C-splines with 3 knots
cat("\n\nTest 3: C-splines with 3 knots (minimal configuration)\n")
cat("-------------------------------------------------------\n")

stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 3  # Minimal for natural cubic splines
)

cat("Fitting C-spline model...\n")

fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Print diagnostic summary
cat("\nC-spline MCMC Diagnostic Summary:\n")
print(fit_c$diagnostic_summary())

# Run smoothing diagnostics
cat("\n")
diagnosis_c <- diagnose_smoothing(fit_c, x, y, stan_data_c, "cspline")
print_smoothing_diagnostics(diagnosis_c)

# Create comparison plots
library(ggplot2)
library(patchwork)

# Extract results for B-splines with 2 knots
draws_b_2knots <- fit_b_2knots$draws(format = "matrix")
x_plot_b_2knots <- colMeans(draws_b_2knots[, grep("x_plot\\[", colnames(draws_b_2knots), value = TRUE)])
y_plot_b_2knots <- colMeans(draws_b_2knots[, grep("y_plot\\[", colnames(draws_b_2knots), value = TRUE)])
y_plot_lower_b_2knots <- apply(draws_b_2knots[, grep("y_plot\\[", colnames(draws_b_2knots), value = TRUE)], 2, quantile, 0.025)
y_plot_upper_b_2knots <- apply(draws_b_2knots[, grep("y_plot\\[", colnames(draws_b_2knots), value = TRUE)], 2, quantile, 0.975)

# Extract results for B-splines with 3 knots
draws_b <- fit_b$draws(format = "matrix")
x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_lower_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
y_plot_upper_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)
sigma_b <- mean(draws_b[, "sigma"])
y_hat_b <- colMeans(draws_b[, grep("y_hat\\[", colnames(draws_b), value = TRUE)])

# Extract results for C-splines
draws_c <- fit_c$draws(format = "matrix")
x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_lower_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
y_plot_upper_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)
sigma_c <- mean(draws_c[, "sigma"])
y_hat_c <- colMeans(draws_c[, grep("y_hat\\[", colnames(draws_c), value = TRUE)])

# High-resolution true function
x_true_hires <- seq(min(x), max(x), length.out = 1000)
y_true_hires <- sin(x_true_hires) + 0.4 * cos(3*x_true_hires) + 0.2*x_true_hires

# B-spline with 2 knots plot
fit_data_b_2knots <- data.frame(
  x = x_plot_b_2knots,
  y_hat = y_plot_b_2knots,
  y_hat_lower = y_plot_lower_b_2knots,
  y_hat_upper = y_plot_upper_b_2knots
)

# Calculate fitted values at data points for diagnostics
y_hat_b_2knots <- colMeans(draws_b_2knots[, grep("y_hat\\[", colnames(draws_b_2knots), value = TRUE)])
sigma_b_2knots <- mean(draws_b_2knots[, "sigma"])

# Run diagnostics for 2-knot model
diagnosis_b_2knots <- diagnose_smoothing(fit_b_2knots, x, y, stan_data_2knots, "bspline")

p_b_2knots <- ggplot() +
  geom_ribbon(data = fit_data_b_2knots, aes(x = x, ymin = y_hat_lower, ymax = y_hat_upper), 
              alpha = 0.3, fill = "skyblue") +
  geom_line(data = fit_data_b_2knots, aes(x = x, y = y_hat), color = "skyblue", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
             color = "black", size = 2, alpha = 0.7) +
  geom_line(data = data.frame(x = x_true_hires, y = y_true_hires), 
            aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "B-spline with 2 Knots (4 Basis Functions) - Minimal Case",
    subtitle = paste0("RMSE = ", round(sqrt(mean((y_hat_b_2knots - y_true)^2)), 3),
                      ", Sigma = ", round(sigma_b_2knots, 3),
                      ", Autocorr = ", round(diagnosis_b_2knots$residual_autocor, 3)),
    x = "x",
    y = "y"
  ) +
  theme_bw()

# B-spline with 3 knots plot
fit_data_b <- data.frame(
  x = x_plot_b,
  y_hat = y_plot_b,
  y_hat_lower = y_plot_lower_b,
  y_hat_upper = y_plot_upper_b
)

p_b <- ggplot() +
  geom_ribbon(data = fit_data_b, aes(x = x, ymin = y_hat_lower, ymax = y_hat_upper), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = fit_data_b, aes(x = x, y = y_hat), color = "blue", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
             color = "black", size = 2, alpha = 0.7) +
  geom_line(data = data.frame(x = x_true_hires, y = y_true_hires), 
            aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "B-spline with 3 Knots (5 Basis Functions)",
    subtitle = paste0("RMSE = ", round(sqrt(mean((y_hat_b - y_true)^2)), 3),
                      ", Sigma = ", round(sigma_b, 3),
                      ", Autocorr = ", round(diagnosis_b$residual_autocor, 3)),
    x = "x",
    y = "y"
  ) +
  theme_bw()

# C-spline plot
fit_data_c <- data.frame(
  x = x_plot_c,
  y_hat = y_plot_c,
  y_hat_lower = y_plot_lower_c,
  y_hat_upper = y_plot_upper_c
)

p_c <- ggplot() +
  geom_ribbon(data = fit_data_c, aes(x = x, ymin = y_hat_lower, ymax = y_hat_upper), 
              alpha = 0.3, fill = "darkgreen") +
  geom_line(data = fit_data_c, aes(x = x, y = y_hat), color = "darkgreen", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
             color = "black", size = 2, alpha = 0.7) +
  geom_line(data = data.frame(x = x_true_hires, y = y_true_hires), 
            aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "C-spline with 3 Knots",
    subtitle = paste0("RMSE = ", round(sqrt(mean((y_hat_c - y_true)^2)), 3),
                      ", Sigma = ", round(sigma_c, 3),
                      ", Autocorr = ", round(diagnosis_c$residual_autocor, 3)),
    x = "x",
    y = "y"
  ) +
  theme_bw()

# Combine plots
combined_plot <- p_b_2knots / p_b / p_c +
  plot_annotation(
    title = "Edge Case Testing: Minimal Knot Configurations",
    subtitle = "All spline configurations struggle with complex function due to limited flexibility",
    caption = "Red dashed = true function, Black points = data, Colored lines = fitted splines with 95% CI"
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-edge_cases_minimal_knots.png", combined_plot, width = 8, height = 12, dpi = 300)
cat("\nSaved comparison plot to output/test-edge_cases_minimal_knots.png\n")

# Summary
cat("\n\nSummary of Edge Case Testing\n")
cat("============================\n")
cat("1. Absolute minimum knots: B-splines = 2 knots, C-splines = 3 knots\n")
cat("2. B-splines with 2 knots: Creates 4 basis functions (boundary knots only)\n")
cat("3. B-splines with 3 knots: Creates 5 basis functions (one interior knot)\n")
cat("4. C-splines with 3 knots: Minimal natural cubic spline configuration\n")
cat("5. Both struggle with complex functions when using minimal knots\n")
cat("\nRecommendation: Use adaptive knot selection (n/2 for B-splines, n/4 for C-splines)\n")

# Close the null device
invisible(dev.off())
