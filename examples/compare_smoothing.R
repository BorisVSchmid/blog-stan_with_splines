# Compare B-splines with different levels of smoothing

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

set.seed(123)
dir.create("output", showWarnings = FALSE)

# Generate test data with noise
n <- 50
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.3 * cos(3*x)  # More complex function
y <- y_true + rnorm(n, 0, 0.15)

# Compile model
cat("Compiling model...\n")
model <- cmdstan_model("code/bsplines.stan")  # Single model handles both smoothing and non-smoothing

# Fit standard B-spline (no smoothing)
cat("Fitting standard B-spline (no smoothing)...\n")
stan_data_standard <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 15,  # Many knots to show overfitting
  spline_degree = 3,
  smoothing_strength = 0  # No smoothing for standard fit
)

fit_standard <- model$sample(
  data = stan_data_standard,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Fit B-splines with different smoothing levels
# Note: smoothing_strength scale: 0=none, 1=mild, 10=strong, 100=very strong
smoothing_levels <- c(1, 5, 10, 100)
fits_smooth <- list()

for (i in seq_along(smoothing_levels)) {
  strength <- smoothing_levels[i]
  cat(sprintf("Fitting B-spline with smoothing_strength = %.0f...\n", strength))
  
  stan_data_smooth <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 15,
    spline_degree = 3,
    smoothing_strength = strength
  )
  
  fits_smooth[[i]] <- model$sample(
    data = stan_data_smooth,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0
  )
}

# Extract results and create plots
cat("Creating comparison plots...\n")

# Function to extract fitted values
extract_fit <- function(fit) {
  draws <- fit$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  data.frame(
    x = x_plot,
    y_mean = y_plot,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper
  )
}

# Extract standard fit
df_standard <- extract_fit(fit_standard)
df_standard$model <- "No smoothing"

# Extract smoothed fits
df_smooth_list <- list()
for (i in seq_along(smoothing_levels)) {
  df_smooth_list[[i]] <- extract_fit(fits_smooth[[i]])
  df_smooth_list[[i]]$model <- sprintf("Strength = %.0f", smoothing_levels[i])
}

# Combine all results
df_all <- bind_rows(df_standard, df_smooth_list)
df_all$model <- factor(df_all$model, 
                      levels = c("No smoothing", sprintf("Strength = %.0f", smoothing_levels)))

# Create faceted plot
p_smoothing <- ggplot(df_all, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = y_mean), color = "blue", linewidth = 1) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
             alpha = 0.5, size = 1) +
  geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y_true), 
            color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ model, ncol = 2, scales = "fixed") +
  labs(title = "B-splines with Different Smoothing Strengths",
       subtitle = "Blue: fitted spline with 95% CI, Red dashed: true function, Points: data",
       x = "x", y = "y") +
  theme_bw()

# Calculate effective degrees of freedom for each model
calc_edf <- function(fit) {
  draws <- fit$draws(format = "matrix")
  y_hat_cols <- grep("y_hat\\[", colnames(draws), value = TRUE)
  if (length(y_hat_cols) > 0) {
    y_hat <- colMeans(draws[, y_hat_cols])
    # Approximate EDF by looking at wiggliness
    sum(diff(diff(y_hat))^2)
  } else {
    NA
  }
}

# Create summary plot showing trade-off
edf_standard <- calc_edf(fit_standard)
edf_smooth <- sapply(fits_smooth, calc_edf)

df_summary <- data.frame(
  strength = c(0, smoothing_levels),
  wiggliness = c(edf_standard, edf_smooth),
  model = c("No smoothing", sprintf("Strength = %.0f", smoothing_levels))
)

p_summary <- ggplot(df_summary, aes(x = strength, y = wiggliness)) +
  geom_point(size = 3) +
  geom_line() +
  geom_text(aes(label = model), hjust = -0.1, vjust = -0.5, size = 3) +
  labs(title = "Smoothing Effect on Wiggliness",
       subtitle = "Lower values indicate smoother fits",
       x = "Smoothing strength", 
       y = "Wiggliness (sum of squared second differences)") +
  theme_bw() +
  scale_x_continuous(limits = c(-10, 110))

# Combine plots
combined_plot <- p_smoothing / p_summary + 
  plot_layout(heights = c(2, 1))

# Save plots
ggsave("output/bspline_smoothing_comparison.png", combined_plot, 
       width = 10, height = 12, dpi = 300)

cat("\nPlot saved: output/bspline_smoothing_comparison.png\n")

# Print diagnostics
cat("\nModel comparison:\n")
cat("=================\n")
for (i in seq_along(c("No smoothing", sprintf("Strength = %.0f", smoothing_levels)))) {
  model_name <- c("No smoothing", sprintf("Strength = %.0f", smoothing_levels))[i]
  wiggliness <- df_summary$wiggliness[i]
  cat(sprintf("%-20s: Wiggliness = %.3f\n", model_name, wiggliness))
}

cat("\nInterpretation:\n")
cat("- No smoothing (0): Flexible fit that may overfit the data\n")
cat("- Mild smoothing (1): Slight smoothing, still quite flexible\n")
cat("- Moderate smoothing (5): Good balance between fit and smoothness\n")  
cat("- Strong smoothing (10): Smoother fit, may miss fine details\n")
cat("- Very strong (100): Very smooth, may underfit complex patterns\n")