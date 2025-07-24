# Full run of regional splines with visualization - fixed version
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# Generate synthetic data for 3 regions
generate_regional_data <- function(n_per_region = 40) {
  # Region 1: Standard sine wave
  x1 <- seq(0, 10, length.out = n_per_region)
  y1_true <- sin(x1)
  y1 <- y1_true + rnorm(n_per_region, 0, 0.15)
  
  # Region 2: Shifted and scaled sine wave
  x2 <- seq(0, 10, length.out = n_per_region)
  y2_true <- 0.7 * sin(x2 + 0.5) + 0.3
  y2 <- y2_true + rnorm(n_per_region, 0, 0.15)
  
  # Region 3: Different frequency
  x3 <- seq(0, 10, length.out = n_per_region)
  y3_true <- sin(1.5 * x3) - 0.2
  y3 <- y3_true + rnorm(n_per_region, 0, 0.15)
  
  # Combine into single dataset
  data.frame(
    x = c(x1, x2, x3),
    y = c(y1, y2, y3),
    y_true = c(y1_true, y2_true, y3_true),
    region = rep(1:3, each = n_per_region),
    region_name = factor(rep(c("Region A", "Region B", "Region C"), each = n_per_region))
  )
}

# Generate data
data <- generate_regional_data(40)

# Prepare Stan data
stan_data <- list(
  n_total = nrow(data),
  n_regions = 3,
  x = data$x,
  y = data$y,
  region = data$region,
  num_knots = 7,
  spline_degree = 3,
  smoothing_strength = 5.0,  # Moderate smoothing
  prior_scale = 2 * sd(data$y)  # Adaptive prior scale based on data variance
)

# Compile and fit model
cat("Compiling regional splines model...\n")
model <- cmdstan_model("tests/test_regional_splines.stan")

cat("Fitting model with hierarchical priors...\n")
fit <- model$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  refresh = 0,
  seed = 123
)

# Extract draws
draws <- fit$draws(format = "matrix")

# Check what variables we have
cat("\nAvailable variables in draws:\n")
var_names <- colnames(draws)
y_plot_vars <- grep("y_plot", var_names, value = TRUE)
cat("y_plot variables found:", length(y_plot_vars), "\n")

# Create simpler plots
dir.create("output", showWarnings = FALSE)

# Plot 1: Raw data
p_data <- ggplot(data, aes(x = x, y = y, color = region_name)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_line(aes(y = y_true), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ region_name, ncol = 1, scales = "free_y") +
  labs(title = "Regional Data and True Functions",
       subtitle = "Points: observed data, Dashed lines: true functions",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

ggsave("output/example-regional_splines_data.png", p_data, width = 8, height = 10, dpi = 300)

# Plot 2: Extract fitted values with credible intervals
y_hat_cols <- grep("y_hat\\[", colnames(draws), value = TRUE)
y_hat_mean <- colMeans(draws[, y_hat_cols])
y_hat_lower <- apply(draws[, y_hat_cols], 2, quantile, 0.025)
y_hat_upper <- apply(draws[, y_hat_cols], 2, quantile, 0.975)

data$y_fitted <- y_hat_mean
data$y_lower <- y_hat_lower
data$y_upper <- y_hat_upper

p_fitted <- ggplot(data, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = region_name), alpha = 0.3) +
  geom_point(aes(y = y, color = region_name), alpha = 0.5) +
  geom_line(aes(y = y_fitted, color = region_name), linewidth = 1.2) +
  geom_line(aes(y = y_true), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ region_name, ncol = 1, scales = "free_y") +
  labs(title = "Regional B-Splines Fit with 95% Credible Intervals",
       subtitle = "Points: data, Solid: fitted splines, Shaded: 95% CI, Dashed: true functions",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

ggsave("output/example-regional_splines_fitted.png", p_fitted, width = 8, height = 10, dpi = 300)

# Plot 3: Global mean coefficients
mu_alpha_cols <- grep("mu_alpha\\[", colnames(draws), value = TRUE)
mu_alpha_mean <- colMeans(draws[, mu_alpha_cols])
mu_alpha_sd <- apply(draws[, mu_alpha_cols], 2, sd)

p_global <- data.frame(
  basis = 1:length(mu_alpha_mean),
  mean = mu_alpha_mean,
  lower = mu_alpha_mean - 2 * mu_alpha_sd,
  upper = mu_alpha_mean + 2 * mu_alpha_sd
) %>%
  ggplot(aes(x = basis)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "darkblue") +
  geom_line(aes(y = mean), color = "darkblue", linewidth = 1.2) +
  geom_point(aes(y = mean), color = "darkblue", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Global Mean Spline Coefficients",
       subtitle = "Hierarchical prior centers regional coefficients around these values",
       x = "Basis Function", y = "Coefficient Value") +
  theme_bw()

ggsave("output/example-regional_splines_global_coefficients.png", p_global, width = 8, height = 6, dpi = 300)

# Plot 4: Regional coefficients
alpha_data <- list()
for (r in 1:3) {
  alpha_cols <- grep(paste0("alpha\\[", r, ","), colnames(draws), value = TRUE, fixed = TRUE)
  if (length(alpha_cols) > 0) {
    alpha_means <- colMeans(draws[, alpha_cols])
    alpha_data[[r]] <- data.frame(
      basis = 1:length(alpha_means),
      coefficient = alpha_means,
      region_name = c("Region A", "Region B", "Region C")[r]
    )
  }
}

if (length(alpha_data) > 0) {
  alpha_df <- bind_rows(alpha_data)
  
  p_coeffs <- ggplot(alpha_df, aes(x = basis, y = coefficient, color = region_name, group = region_name)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = "Regional Spline Coefficients",
         subtitle = "Each region has unique coefficients but shares information through hierarchical prior",
         x = "Basis Function", y = "Coefficient Value",
         color = "Region") +
    theme_bw() +
    scale_color_brewer(palette = "Set2")
  
  ggsave("output/example-regional_splines_coefficients_comparison.png", p_coeffs, width = 10, height = 6, dpi = 300)
}

# Plot 5: Variance components
var_between <- mean(draws[, "var_between"])
var_within <- mean(draws[, "var_within"])
icc <- mean(draws[, "icc"])

var_df <- data.frame(
  component = c("Between Regions", "Within Regions"),
  variance = c(var_between, var_within),
  percentage = c(var_between, var_within) / (var_between + var_within) * 100
)

p_var <- ggplot(var_df, aes(x = component, y = variance, fill = component)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            vjust = -0.5, size = 4) +
  labs(title = "Variance Components",
       subtitle = sprintf("ICC = %.3f (%.1f%% of variance is between regions)", 
                         icc, icc * 100),
       x = "", y = "Variance") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

ggsave("output/example-regional_splines_variance_components.png", p_var, width = 6, height = 5, dpi = 300)

# Print summary
cat("\n======================================\n")
cat("Regional B-Splines Model Summary\n")
cat("======================================\n")

# Regional effects
beta_region <- colMeans(draws[, grep("beta_region\\[", colnames(draws), value = TRUE)])
cat("\nRegional baseline effects:\n")
for (i in 1:3) {
  cat(sprintf("  Region %s: %6.3f\n", c("A", "B", "C")[i], beta_region[i]))
}

# Global mean coefficients
cat("\nGlobal mean spline coefficients:\n")
for (i in 1:length(mu_alpha_mean)) {
  cat(sprintf("  Basis %d: %6.3f (SD: %.3f)\n", i, mu_alpha_mean[i], mu_alpha_sd[i]))
}

# Variance components
cat(sprintf("\nVariance decomposition:\n"))
cat(sprintf("  Between regions: %.3f (%.1f%%)\n", var_between, var_between/(var_between + var_within)*100))
cat(sprintf("  Within regions:  %.3f (%.1f%%)\n", var_within, var_within/(var_between + var_within)*100))
cat(sprintf("  ICC: %.3f\n", icc))

# Model diagnostics
cat("\nModel diagnostics:\n")
cat(sprintf("  Number of divergences: %d\n", sum(fit$sampler_diagnostics()[,,"divergent__"])))
cat(sprintf("  Max R-hat: %.3f\n", max(fit$summary()$rhat, na.rm = TRUE)))
cat(sprintf("  Min ESS Bulk: %.0f\n", min(fit$summary()$ess_bulk, na.rm = TRUE)))

cat("\nOutputs saved:\n")
cat("- output/example-regional_splines_data.png\n")
cat("- output/example-regional_splines_fitted.png\n")
cat("- output/example-regional_splines_global_coefficients.png\n")
cat("- output/example-regional_splines_coefficients_comparison.png\n")
cat("- output/example-regional_splines_variance_components.png\n")

cat("\nInterpretation:\n")
cat("- The hierarchical model allows regions to have unique spline shapes\n")
cat("- Regions borrow strength from each other through shared priors\n")
cat("- The ICC shows what proportion of variance is between vs within regions\n")
cat("- Higher ICC means regions are more different from each other\n")