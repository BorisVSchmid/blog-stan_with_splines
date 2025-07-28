# Minimal B-spline example demonstrating key features

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

# Source smoothing diagnostics
source("code/smoothing_diagnostics.R")

# Note on knot selection:
# We use simple rules based on sample size rather than trying to estimate
# function complexity from noisy data. The key insight is that the smoothing
# prior in B-splines is scale-invariant in x-space (operates knot-to-knot),
# so we can use a generous number of knots and let the smoothing prior
# control flexibility.

# Generate test data with complex function
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.25*x 
y <- y_true + rnorm(n, 0, 0.15)

# Key features demonstration:
# 1. Simple knot selection based on sample size
# Use n/2 knots for B-splines (with smoothing to control flexibility)
adaptive_knots <- max(4, min(round(n/2), 40))

# 2. Adaptive prior scale based on data variance
adaptive_prior <- 2 * sd(y)

# 3. Default smoothing for more stable fits
# smoothing_strength: 0=none, 0.05-0.1=mild, 0.1-0.2=strong

# Prepare data for Stan
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = adaptive_knots,     # Use adaptive knot selection (n/2)
  spline_degree = 3,              # Cubic splines (most common)
  smoothing_strength = 0.1,       # Default mild smoothing with new scaling formula
  prior_scale = adaptive_prior    # Data-driven prior
)

# Compile and fit model
cat("Compiling B-spline model...\n")
model <- cmdstan_model("code/bsplines.stan")

cat("Fitting model with:\n")
cat("  - Number of knots:", stan_data$num_knots, "(n/2 rule)\n")
cat("  - Prior scale:", round(adaptive_prior, 2), "(based on data variance)\n")
cat("  - Smoothing strength:", stan_data$smoothing_strength, " (0.1 is default)\n")

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Print diagnostic summary
cat("\nMCMC Diagnostic Summary:\n")
print(fit$diagnostic_summary())

# Run smoothing diagnostics
cat("\n")
diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "bspline")
print_smoothing_diagnostics(diagnosis)

# Generate and save diagnostic plots
cat("\nGenerating diagnostic plots...\n")
diagnostic_plot <- plot_diagnostic_residuals(fit, x, y)
ggsave("output/code-bspline_diagnostics.png", diagnostic_plot, width = 12, height = 8, dpi = 300)
cat("Diagnostic plots saved to output/code-bspline_diagnostics.png\n")

# Extract and plot results
draws <- fit$draws(format = "matrix")

# Extract smooth plotting grid from generated quantities
x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)

# Extract sigma for model diagnostics
sigma <- mean(draws[, "sigma"])

# Create data frames for plotting
# Smooth fit data
fit_data <- data.frame(
  x = x_plot,
  y_hat = y_plot,
  y_hat_lower = y_plot_lower,
  y_hat_upper = y_plot_upper
)

# Original data with true function
data_points <- data.frame(
  x = x,
  y = y,
  y_true = y_true
)

# High-resolution true function for smooth plotting
x_true_hires <- seq(min(x), max(x), length.out = 1000)
y_true_hires <- sin(x_true_hires) + 0.4 * cos(3*x_true_hires) + 0.25*x_true_hires
true_function_data <- data.frame(
  x = x_true_hires,
  y_true = y_true_hires
)

# Create plot with ggplot2
p <- ggplot() +
  geom_ribbon(data = fit_data, aes(x = x, ymin = y_hat_lower, ymax = y_hat_upper, fill = "95% CI"), alpha = 0.3) +
  geom_line(data = fit_data, aes(x = x, y = y_hat, color = "B-spline fit"), linewidth = 1.2) +
  geom_point(data = data_points, aes(x = x, y = y, color = "Data"), size = 2, alpha = 0.7) +
  geom_line(data = true_function_data, aes(x = x, y = y_true, color = "True function"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("Data" = "black", "B-spline fit" = "blue", "True function" = "red")) +
  scale_fill_manual(values = c("95% CI" = "blue")) +
  labs(
    title = "B-spline Fit Demonstrating Key Features",
    subtitle = "True function: sin(x) + 0.4*cos(3x) + 0.25*x",
    caption = paste0("Parameters: knots = ", stan_data$num_knots, 
                     ", smoothing_strength = ", stan_data$smoothing_strength,
                     ", Estimated sigma = ", round(sigma, 3)),
    x = "x",
    y = "y"
  ) +
  ylim(-1.6, 3.8) +  # Adjusted range for 0.2*x linear term
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(order = 2), color = guide_legend(order = 1))

# Display the plot
print(p)

# Save the plot
dir.create("output", showWarnings = FALSE)
ggsave("output/code-bspline_minimal_example.png", p, width = 8, height = 6, dpi = 300)

# The diagnostic function above already prints all necessary information
