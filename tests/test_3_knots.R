# Test B-splines with 3 knots

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

# Force 3 knots
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 3,  # Minimal knots
  spline_degree = 3,
  smoothing_strength = 1.0,
  prior_scale = 2 * sd(y)
)

# Compile and fit model
cat("Testing B-splines with 3 knots\n")
cat("==============================\n")
cat("Number of basis functions: 3 + 3 - 1 = 5\n\n")

model <- cmdstan_model("code/bsplines.stan")

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

# Extract sigma and y_hat for diagnostics
sigma <- mean(draws[, "sigma"])
y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])

# Create data frames for plotting
fit_data <- data.frame(
  x = x_plot,
  y_hat = y_plot,
  y_hat_lower = y_plot_lower,
  y_hat_upper = y_plot_upper
)

data_points <- data.frame(
  x = x,
  y = y,
  y_true = y_true,
  y_hat = y_hat
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
    title = "B-spline Fit with 3 Knots (5 Basis Functions)",
    subtitle = "True function: sin(x) + 0.4*cos(3x) + 0.2*x",
    caption = paste0("EDF = ", round(diagnosis$edf, 1), 
                     ", Sigma = ", round(sigma, 3),
                     ", Autocorrelation = ", round(diagnosis$residual_autocor, 3)),
    x = "x",
    y = "y"
  ) +
  ylim(-1.6, 3.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(order = 2), color = guide_legend(order = 1))

# Display the plot
print(p)

# Save the plot
ggsave("claude-tmp/bspline_3_knots.png", p, width = 8, height = 6, dpi = 300)

# Print summary
cat("\nModel Summary:\n")
cat("==============\n")
cat("Estimated noise (σ):", round(sigma, 3), "\n")
cat("Number of knots specified:", 3, "\n")
cat("Number of basis functions:", 5, "(knots + degree - 1)\n")
cat("RMSE from true function:", round(sqrt(mean((data_points$y_hat - data_points$y_true)^2)), 3), "\n")