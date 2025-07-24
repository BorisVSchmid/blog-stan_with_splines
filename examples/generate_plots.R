# Standalone script to generate plots
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)
dir.create("output", showWarnings = FALSE)

# Generate test data
generate_data <- function(n = 30) {
  x <- seq(0, 10, length.out = n)
  y_true <- sin(x)
  y <- y_true + rnorm(n, 0, 0.1)
  list(x = x, y = y, y_true = y_true)
}

# Fit and plot
data <- generate_data(30)

# B-spline
cat("Fitting B-spline...\n")
stan_data_b <- list(
  n_data = length(data$x),
  x = data$x,
  y = data$y,
  num_knots = 7,
  spline_degree = 3
)

model_b <- cmdstan_model("test_bsplines.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

# C-spline
cat("Fitting C-spline...\n")
stan_data_c <- list(
  n_data = length(data$x),
  x = data$x,
  y = data$y,
  num_knots = 7
)

model_c <- cmdstan_model("test_csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

# Extract results and plot
cat("Creating plots...\n")
draws_b <- fit_b$draws(format = "matrix")
draws_c <- fit_c$draws(format = "matrix")

# B-spline plot
x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
y_lower_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
y_upper_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)

plot_b <- ggplot() +
  geom_ribbon(aes(x = x_plot_b, ymin = y_lower_b, ymax = y_upper_b), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = x_plot_b, y = y_plot_b), color = "blue", size = 1) +
  geom_point(data = data.frame(x = data$x, y = data$y), aes(x, y), alpha = 0.6) +
  geom_line(data = data.frame(x = data$x, y = data$y_true), aes(x, y), 
            color = "red", linetype = "dashed") +
  labs(title = "B-spline fit", x = "x", y = "y") +
  theme_minimal()

# C-spline plot
x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
y_lower_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
y_upper_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)

plot_c <- ggplot() +
  geom_ribbon(aes(x = x_plot_c, ymin = y_lower_c, ymax = y_upper_c), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = x_plot_c, y = y_plot_c), color = "blue", size = 1) +
  geom_point(data = data.frame(x = data$x, y = data$y), aes(x, y), alpha = 0.6) +
  geom_line(data = data.frame(x = data$x, y = data$y_true), aes(x, y), 
            color = "red", linetype = "dashed") +
  labs(title = "C-spline (Natural Cubic Spline) fit", x = "x", y = "y") +
  theme_minimal()

# Basis functions for B-spline
basis_cols <- grep("basis_functions\\[", colnames(draws_b), value = TRUE)
if (length(basis_cols) > 0) {
  n_basis <- 7  # We know this from the model
  basis_data <- NULL
  for (i in 1:n_basis) {
    cols <- grep(paste0("basis_functions\\[\\d+,", i, "\\]"), colnames(draws_b), value = TRUE)
    if (length(cols) > 0) {
      basis_vals <- colMeans(draws_b[, cols])
      temp_df <- data.frame(
        x = x_plot_b,
        y = basis_vals,
        basis = factor(i)
      )
      basis_data <- rbind(basis_data, temp_df)
    }
  }
  
  plot_basis <- ggplot(basis_data, aes(x = x, y = y, color = basis)) +
    geom_line(size = 1) +
    labs(title = "B-spline basis functions", x = "x", y = "Basis value") +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave("output/example-bspline_basis_functions.png", plot_basis, width = 8, height = 6, dpi = 300)
  cat("Saved basis functions plot\n")
}

# Combined plot
combined <- plot_b | plot_c
combined_with_title <- combined + 
  plot_annotation(
    title = "B-splines vs C-splines: Sine Function Fitting",
    subtitle = "Blue: fitted spline with 95% CI, Red dashed: true function, Points: noisy data"
  )

ggsave("output/example-spline_comparison.png", combined_with_title, width = 12, height = 6, dpi = 300)
cat("Saved comparison plot\n")

cat("\nPlots saved to output/ directory\n")