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

# Generate test data with complex function
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x 
y <- y_true + rnorm(n, 0, 0.15)

# Key features demonstration:
# 1. Adaptive knot selection based on data size
# Rule of thumb: n/8 for B-splines to avoid overfitting
# (more conservative than n/4 due to B-spline flexibility)
adaptive_knots <- max(4, min(round(n/8), 20))

# 2. Adaptive prior scale based on data variance
adaptive_prior <- 2 * sd(y)

# 3. Choose between smoothing and flexibility
use_smoothing <- FALSE  # Set to TRUE for smoother curves

# Prepare data for Stan
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = adaptive_knots,     # Using adaptive selection
  spline_degree = 3,              # Cubic splines (most common)
  tau_smooth = if(use_smoothing) 0.2 else 0,  # Smoothing on/off
  prior_scale = adaptive_prior    # Data-driven prior
)

# Compile and fit model
cat("Compiling B-spline model...\n")
model <- cmdstan_model("code/bsplines.stan")

cat("Fitting model with:\n")
cat("  - Adaptive knots:", adaptive_knots, "(based on n =", n, ")\n")
cat("  - Prior scale:", round(adaptive_prior, 2), "(based on data variance)\n")
cat("  - Smoothing:", if(use_smoothing) "enabled (τ = 0.2)" else "disabled", "\n\n")

fit <- model$sample(
  data = stan_data,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Run smoothing diagnostics
cat("\n")
diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "bspline")
print_smoothing_diagnostics(diagnosis)

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
y_true_hires <- sin(x_true_hires) + 0.4 * cos(3*x_true_hires) + 0.2*x_true_hires
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
    subtitle = "True function: sin(x) + 0.4*cos(3x) + 0.2*x",
    caption = paste0("Parameters: Adaptive knots = ", stan_data$num_knots, 
                     ", Adaptive prior = ", round(stan_data$prior_scale, 1),
                     ", Smoothing = ", if(use_smoothing) paste0("tau = ", stan_data$tau_smooth) else "off",
                     ", Estimated sigma = ", round(sigma, 3)),
    x = "x",
    y = "y"
  ) +
  ylim(-1.6, 3.6) +  # Adjusted range for 0.2*x linear term
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(order = 2), color = guide_legend(order = 1))

# Display the plot
print(p)

# Save the plot
dir.create("output", showWarnings = FALSE)
ggsave("output/bspline_minimal_example.png", p, width = 8, height = 6, dpi = 300)

# Print summary
cat("\nModel Summary:\n")
cat("==============\n")
cat("Estimated noise (σ):", round(sigma, 3), "\n")
cat("Number of knots specified:", adaptive_knots, "\n")
cat("Number of basis functions:", adaptive_knots + stan_data$spline_degree - 1, "(knots + degree - 1)\n")
cat("Effective smoothing:", if(use_smoothing) "random walk prior" else "independent priors", "\n")

# Diagnostic check
cat("\nTo try different settings, modify:\n")
cat("  - use_smoothing: TRUE for smooth curves, FALSE for flexible fit\n")
cat("  - n: Sample size (affects adaptive knot selection)\n")
cat("  - Add more noise to y to see adaptive prior scaling\n")
