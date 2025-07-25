# Compare B-splines vs C-splines smoothing approaches
# B-splines: Control smoothing via smoothing_strength parameter
# C-splines: Control smoothing via number of knots

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

# Generate complex test data
set.seed(123)
n <- 60
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x
y <- y_true + rnorm(n, 0, 0.15)

cat("Comparing B-splines vs C-splines Smoothing Approaches\n")
cat("=====================================================\n")
cat("Complex function: sin(x) + 0.4*cos(3x) + 0.2*x\n")
cat("Data: n =", n, "points with noise SD = 0.15\n\n")

# Compile models once
cat("Compiling models...\n")
model_b <- cmdstan_model("code/bsplines.stan")
model_c <- cmdstan_model("code/csplines.stan")

# Storage for results
results_b <- list()
results_c <- list()

# B-splines: Fixed knots (default n/2), varying smoothing_strength
cat("\n1. B-splines with fixed knots (default n/2), varying smoothing_strength\n")
cat("---------------------------------------------------------------\n")
# Choose smoothing values to show progression from no smoothing to strong smoothing
fixed_knots_b <- max(4, min(round(n/2), 40))  # Using default n/2 rule
num_basis_b <- fixed_knots_b + 3 - 1  # degree = 3
# Progressive smoothing: 0=none, 2=mild, 4=moderate, 6=moderate-strong, 8=strong
smoothing_values <- c(0, 2, 4, 6, 8)

for (smooth in smoothing_values) {
  cat("  Smoothing strength =", smooth, "... ")
  
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = fixed_knots_b,
    spline_degree = 3,
    smoothing_strength = smooth,
    prior_scale = 2 * sd(y)
  )
  
  fit <- model_b$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Get diagnostics
  diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "bspline")
  
  # Extract fitted values
  draws <- fit$draws(format = "matrix")
  
  # Get fitted values at data points for RMSE calculation
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  
  # Get smooth plotting grid for visualization
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  results_b[[as.character(smooth)]] <- list(
    smoothing = smooth,
    x_plot = x_plot,
    y_plot = y_plot,
    y_plot_lower = y_plot_lower,
    y_plot_upper = y_plot_upper,
    edf = diagnosis$edf,
    sigma = diagnosis$sigma_estimate,
    rmse = sqrt(mean((y_hat - y_true)^2))
  )
  
  cat("RMSE =", round(results_b[[as.character(smooth)]]$rmse, 3), "\n")
}

# C-splines: Varying number of knots (reversed order for display)
cat("\n2. C-splines with varying number of knots\n")
cat("-----------------------------------------\n")
knot_values <- c(13, 11, 9, 7, 5)  # Reversed to show high to low

for (k in knot_values) {
  cat("  Number of knots =", k, "... ")
  
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = k
  )
  
  fit <- model_c$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Get diagnostics
  diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "cspline")
  
  # Extract fitted values
  draws <- fit$draws(format = "matrix")
  
  # Get fitted values at data points for RMSE calculation
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  
  # Get smooth plotting grid for visualization
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  results_c[[as.character(k)]] <- list(
    knots = k,
    x_plot = x_plot,
    y_plot = y_plot,
    y_plot_lower = y_plot_lower,
    y_plot_upper = y_plot_upper,
    edf = diagnosis$edf,
    sigma = diagnosis$sigma_estimate,
    rmse = sqrt(mean((y_hat - y_true)^2))
  )
  
  cat("RMSE =", round(results_c[[as.character(k)]]$rmse, 3), "\n")
}

# Create comparison plots
cat("\n3. Creating comparison plots...\n")

# Prepare data for plotting
library(ggplot2)
library(patchwork)

# True function for reference
x_true_hires <- seq(min(x), max(x), length.out = 200)
y_true_hires <- sin(x_true_hires) + 0.4 * cos(3*x_true_hires) + 0.2*x_true_hires

# B-spline plots
plots_b <- list()
for (i in seq_along(smoothing_values)) {
  smooth <- smoothing_values[i]
  res <- results_b[[as.character(smooth)]]
  
  df_fit <- data.frame(
    x = res$x_plot,
    y = res$y_plot,
    ymin = res$y_plot_lower,
    ymax = res$y_plot_upper
  )
  
  p <- ggplot() +
    geom_ribbon(data = df_fit, aes(x = x, ymin = ymin, ymax = ymax), 
                alpha = 0.2, fill = "blue") +
    geom_line(data = df_fit, aes(x = x, y = y), color = "blue", linewidth = 1) +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
               alpha = 0.3, size = 0.8) +
    geom_line(data = data.frame(x = x_true_hires, y = y_true_hires), 
              aes(x = x, y = y), color = "red", linetype = "dashed", alpha = 0.7) +
    labs(title = paste0("B-spline: smoothing_strength = ", smooth),
         subtitle = paste0("RMSE = ", round(res$rmse, 3))) +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9))
  
  plots_b[[i]] <- p
}

# C-spline plots
plots_c <- list()
for (i in seq_along(knot_values)) {
  k <- knot_values[i]
  res <- results_c[[as.character(k)]]
  
  df_fit <- data.frame(
    x = res$x_plot,
    y = res$y_plot,
    ymin = res$y_plot_lower,
    ymax = res$y_plot_upper
  )
  
  p <- ggplot() +
    geom_ribbon(data = df_fit, aes(x = x, ymin = ymin, ymax = ymax), 
                alpha = 0.2, fill = "darkgreen") +
    geom_line(data = df_fit, aes(x = x, y = y), color = "darkgreen", linewidth = 1) +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
               alpha = 0.3, size = 0.8) +
    geom_line(data = data.frame(x = x_true_hires, y = y_true_hires), 
              aes(x = x, y = y), color = "red", linetype = "dashed", alpha = 0.7) +
    labs(title = paste0("C-spline: knots = ", k),
         subtitle = paste0("RMSE = ", round(res$rmse, 3))) +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9))
  
  plots_c[[i]] <- p
}

# Combine all plots in vertical columns
combined_plot <- (wrap_plots(plots_b, ncol = 1) | wrap_plots(plots_c, ncol = 1)) +
  plot_annotation(
    title = "B-splines vs C-splines: Different Approaches to Smoothing",
    subtitle = paste0("B-splines: Fixed ", fixed_knots_b, " knots (", num_basis_b, " basis), varying smoothing | C-splines: Varying number of knots"),
    caption = "Red dashed line = true function | Blue/Green = fitted spline with 95% CI"
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-flexibility_comparison.png", combined_plot, width = 10, height = 14, dpi = 300)
cat("Saved comparison plot to output/test-flexibility_comparison.png\n")

# Create summary comparison plot
cat("\n4. Creating summary comparison...\n")

# Extract summary statistics
summary_b <- data.frame(
  method = "B-spline",
  parameter = smoothing_values,
  parameter_name = "smoothing_strength",
  edf = sapply(results_b, function(x) x$edf),
  rmse = sapply(results_b, function(x) x$rmse)
)

summary_c <- data.frame(
  method = "C-spline",
  parameter = knot_values,
  parameter_name = "num_knots",
  edf = sapply(results_c, function(x) x$edf),
  rmse = sapply(results_c, function(x) x$rmse)
)
# Sort C-spline summary back to ascending order for plotting
summary_c <- summary_c[order(summary_c$parameter), ]

# RMSE comparison is more meaningful than EDF for regularized splines

# RMSE comparison
# Calculate appropriate y-axis limits based on actual RMSE range
rmse_min <- min(c(summary_b$rmse, summary_c$rmse))
rmse_max <- max(c(summary_b$rmse, summary_c$rmse))
rmse_range <- rmse_max - rmse_min
rmse_ylim <- c(max(0, rmse_min - 0.1 * rmse_range), rmse_max + 0.1 * rmse_range)

p_rmse_b <- ggplot(summary_b, aes(x = parameter, y = rmse)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 3) +
  labs(title = "B-splines: RMSE vs Smoothing Strength",
       x = "Smoothing Strength",
       y = "RMSE") +
  theme_bw() +
  coord_cartesian(ylim = rmse_ylim)

p_rmse_c <- ggplot(summary_c, aes(x = parameter, y = rmse)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 3) +
  labs(title = "C-splines: RMSE vs Number of Knots",
       x = "Number of Knots",
       y = "RMSE") +
  theme_bw() +
  coord_cartesian(ylim = rmse_ylim)

summary_plot <- (p_rmse_b | p_rmse_c) +
  plot_annotation(
    title = "RMSE Comparison: B-splines vs C-splines Smoothing Control",
    subtitle = "Complex function: sin(x) + 0.4*cos(3x) + 0.2*x with noise"
  )

ggsave("output/test-flexibility_comparison_summary.png", summary_plot, width = 10, height = 8, dpi = 300)
cat("Saved summary plot to output/test-flexibility_comparison_summary.png\n")

# Print summary table
cat("\n5. Summary Table\n")
cat("================\n\n")
cat(paste0("B-splines (", fixed_knots_b, " knots, ", num_basis_b, " basis functions, varying smoothing_strength):\n"))
print(summary_b[, c("parameter", "rmse")], row.names = FALSE)

cat("\nC-splines (varying number of knots):\n")
print(summary_c[, c("parameter", "rmse")], row.names = FALSE)

# Find optimal settings
opt_b <- summary_b[which.min(summary_b$rmse), ]
opt_c <- summary_c[which.min(summary_c$rmse), ]

cat("\nOptimal settings for this complex function:\n")
cat("-------------------------------------------\n")
cat("B-splines: smoothing_strength =", opt_b$parameter, 
    "(RMSE =", round(opt_b$rmse, 3), ")\n")
cat("C-splines: num_knots =", opt_c$parameter, 
    "(RMSE =", round(opt_c$rmse, 3), ")\n")

cat("\nNote: B-spline smoothing progression:\n")
cat("  Smoothing 0 = No smoothing (maximum flexibility)\n")
cat("  Smoothing 2 = Mild smoothing (recommended default)\n")
cat("  Smoothing 4 = Moderate smoothing\n")
cat("  Smoothing 6 = Moderate-strong smoothing\n")
cat("  Smoothing 8 = Strong smoothing\n")
cat("\nC-spline knot progression:\n")
cat("  13 knots = High flexibility\n")
cat("  11 knots = Moderate-high flexibility\n")
cat("  9 knots = Moderate flexibility\n")
cat("  7 knots = Limited flexibility\n")
cat("  5 knots = Minimal flexibility\n")

# Close the null device
dev.off()