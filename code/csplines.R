# Minimal viable C-spline (natural cubic spline) implementation in Stan

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
# Use groundhog for dependencies
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

# Load cmdstanr after its dependencies
library(cmdstanr)

# Source smoothing diagnostics
source("code/smoothing_diagnostics.R")

# Generate test data with complex function
set.seed(123)
n <- 30
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.25*x
y <- y_true + rnorm(n, 0, 0.15)

# Prepare data for Stan
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 5       # Number of knots for natural cubic spline
                      # - Use 4-6 knots for simple smooth curves
                      # - Use 7-10 knots for more complex patterns
                      # - C-splines are smoother than B-splines, so often need fewer knots
                      # - Too many knots can lead to numerical instability
)
# Note: C-splines will extrapolate linearly beyond the data range

# Compile and fit model
model <- cmdstan_model("code/csplines.stan")

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
diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "cspline")
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
y_true_hires <- sin(x_true_hires) + 0.4 * cos(3*x_true_hires) + 0.25*x_true_hires
true_function_data <- data.frame(
  x = x_true_hires,
  y_true = y_true_hires
)

# Create plot with ggplot2
p <- ggplot() +
  geom_ribbon(data = fit_data, aes(x = x, ymin = y_hat_lower, ymax = y_hat_upper, fill = "95% CI"), alpha = 0.3) +
  geom_line(data = fit_data, aes(x = x, y = y_hat, color = "C-spline fit"), linewidth = 1.2) +
  geom_point(data = data_points, aes(x = x, y = y, color = "Data"), size = 2, alpha = 0.7) +
  geom_line(data = true_function_data, aes(x = x, y = y_true, color = "True function"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("Data" = "black", "C-spline fit" = "darkgreen", "True function" = "red")) +
  scale_fill_manual(values = c("95% CI" = "darkgreen")) +
  labs(
    title = "C-spline (Natural Cubic Spline) Fit with 95% Credible Interval",
    subtitle = "True function: sin(x) + 0.4*cos(3x) + 0.25*x",
    caption = paste0("Parameters: num_knots = ", stan_data$num_knots, 
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
ggsave("output/code-cspline_minimal_example.png", p, width = 8, height = 6, dpi = 300)
