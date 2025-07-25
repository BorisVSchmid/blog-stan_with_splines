# Quick test of spline implementations
library(cmdstanr)
set.seed(123)

# Generate simple test data
n <- 20
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.1)

# Test B-spline
cat("Testing B-spline model...\n")
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 5,
  spline_degree = 3,
  smoothing_strength = 0.0, # if we only do 5 knots, then do little to no smoothing.
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
  num_knots = 5
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

p_b <- ggplot() +
  geom_ribbon(data = df_b, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.3, fill = "blue") +
  geom_line(data = df_b, aes(x = x, y = y), color = "blue", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), alpha = 0.6) +
  geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y), 
            color = "red", linetype = "dashed", alpha = 0.7) +
  labs(title = "B-spline Fit", 
       subtitle = paste0(stan_data_b$num_knots, " knots, smoothing = ", stan_data_b$smoothing_strength),
       x = "x", y = "y") +
  theme_bw()

# Create C-spline plot
df_c <- data.frame(
  x = x_plot_c,
  y = y_plot_c,
  ymin = y_plot_c_lower,
  ymax = y_plot_c_upper
)

p_c <- ggplot() +
  geom_ribbon(data = df_c, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.3, fill = "darkgreen") +
  geom_line(data = df_c, aes(x = x, y = y), color = "darkgreen", linewidth = 1.2) +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), alpha = 0.6) +
  geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y), 
            color = "red", linetype = "dashed", alpha = 0.7) +
  labs(title = "C-spline Fit", 
       subtitle = paste0(stan_data_c$num_knots, " knots"),
       x = "x", y = "y") +
  theme_bw()

# Combine plots
combined_plot <- p_b | p_c
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Basic Spline Comparison: B-splines vs C-splines",
    subtitle = "Both fitting sin(x) with 5 knots. Red dashed = true function, points = data",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-basic_spline_comparison.png", combined_plot, width = 10, height = 5, dpi = 300)
cat("Saved comparison plot to output/test-basic_spline_comparison.png\n")
