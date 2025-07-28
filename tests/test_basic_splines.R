# Quick test of spline implementations

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

# Generate simple test data
n <- 20
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.1)

# Calculate default knot counts
# B-splines: n/2 knots (capped at 40)
# C-splines: n/4 knots (capped at 20)
bspline_knots <- max(4, min(round(n/2), 40))
cspline_knots <- max(4, min(round(n/4), 20))

cat(sprintf("Using default adaptive knot selection:\n"))
cat(sprintf("  B-splines: %d knots (n/2 rule)\n", bspline_knots))
cat(sprintf("  C-splines: %d knots (n/4 rule)\n\n", cspline_knots))

# Test B-spline
cat("Testing B-spline model...\n")
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = bspline_knots,
  spline_degree = 3,
  smoothing_strength = 1.0, # Good default for sine waves with adaptive knots
  prior_scale = 2 * sd(y)
)

model_b <- cmdstan_model("code/bsplines.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

cat("\nB-spline MCMC diagnostics:\n")
print(fit_b$diagnostic_summary())

cat("B-spline fit summary:\n")
fit_b$summary(c("sigma", "lp__"))

# Test C-spline
cat("\nTesting C-spline model...\n")
stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = cspline_knots
)

model_c <- cmdstan_model("code/csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

cat("\nC-spline MCMC diagnostics:\n")
print(fit_c$diagnostic_summary())

cat("C-spline fit summary:\n")
fit_c$summary(c("sigma", "lp__"))

cat("\nBoth models ran successfully!\n")

# Create comparison plot
cat("\nCreating comparison plot...\n")

library(ggplot2)
library(patchwork)

# Extract B-spline results
draws_b <- fit_b$draws(format = "matrix")
x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_b_lower <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
y_plot_b_upper <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)

# Extract C-spline results
draws_c <- fit_c$draws(format = "matrix")
x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_c_lower <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
y_plot_c_upper <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)

# Create B-spline plot
df_b <- data.frame(
  x = x_plot_b,
  y = y_plot_b,
  ymin = y_plot_b_lower,
  ymax = y_plot_b_upper
)

# Create C-spline plot
df_c <- data.frame(
  x = x_plot_c,
  y = y_plot_c,
  ymin = y_plot_c_lower,
  ymax = y_plot_c_upper
)

# Calculate common y-axis limits
y_min <- min(c(df_b$ymin, df_c$ymin, y, y_true))
y_max <- max(c(df_b$ymax, df_c$ymax, y, y_true))
y_range <- y_max - y_min
y_limits <- c(y_min - 0.05 * y_range, y_max + 0.05 * y_range)

p_b <- ggplot() +
  geom_ribbon(data = df_b, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.3, fill = "blue") +
  geom_line(data = df_b, aes(x = x, y = y), color = "blue", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), alpha = 0.6) +
  geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y), 
            color = "red", linetype = "dashed", alpha = 0.7) +
  labs(title = "B-spline Fit", 
       subtitle = paste0(stan_data_b$num_knots, " knots, smoothing = ", stan_data_b$smoothing_strength),
       x = "x", y = "y") +
  theme_bw() +
  ylim(y_limits)

p_c <- ggplot() +
  geom_ribbon(data = df_c, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.3, fill = "darkgreen") +
  geom_line(data = df_c, aes(x = x, y = y), color = "darkgreen", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), alpha = 0.6) +
  geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y), 
            color = "red", linetype = "dashed", alpha = 0.7) +
  labs(title = "C-spline Fit", 
       subtitle = paste0(stan_data_c$num_knots, " knots"),
       x = "x", y = "y") +
  theme_bw() +
  ylim(y_limits)

# Combine plots
combined_plot <- p_b + p_c
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Basic Spline Comparison: B-splines vs C-splines",
    subtitle = sprintf("Fitting sin(x) with default adaptive knots (B-splines: %d, C-splines: %d). Red dashed = true function, points = data", 
                      bspline_knots, cspline_knots),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-basic_spline_comparison.png", combined_plot, width = 10, height = 5, dpi = 300)
cat("Saved comparison plot to output/test-basic_spline_comparison.png\n")
