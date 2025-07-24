# Demo: Regional B-splines with hierarchical priors
# Shows how 3 regions can have different spline shapes while sharing information

library(cmdstanr)
library(ggplot2)
library(dplyr)

set.seed(42)

# Generate synthetic data for 3 regions
n_per_region <- 30
x_vals <- seq(0, 10, length.out = n_per_region)

# Region A: Standard sine wave
y_A <- sin(x_vals) + rnorm(n_per_region, 0, 0.15)

# Region B: Shifted and scaled sine wave  
y_B <- 0.7 * sin(x_vals + 0.5) + 0.3 + rnorm(n_per_region, 0, 0.15)

# Region C: Different frequency
y_C <- sin(1.5 * x_vals) - 0.2 + rnorm(n_per_region, 0, 0.15)

# Combine data
data <- data.frame(
  x = rep(x_vals, 3),
  y = c(y_A, y_B, y_C),
  region = rep(1:3, each = n_per_region),
  region_name = factor(rep(c("Region A", "Region B", "Region C"), each = n_per_region))
)

# Plot raw data
p_data <- ggplot(data, aes(x = x, y = y, color = region_name)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ region_name, ncol = 1, scales = "free_y") +
  labs(title = "Regional Data: 3 Different Patterns",
       subtitle = "Each region follows a different nonlinear pattern",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/example-regional_data.png", p_data, width = 8, height = 8, dpi = 300)

# Prepare Stan data
stan_data <- list(
  n_total = nrow(data),
  n_regions = 3,
  x = data$x,
  y = data$y,
  region = data$region,
  num_knots = 7,
  spline_degree = 3
)

cat("\nRegional Splines Model Demo\n")
cat("===========================\n")
cat("\nData summary:\n")
cat("- 3 regions with different patterns\n")
cat("- 30 observations per region\n")
cat("- Using B-splines with 7 knots\n")
cat("\nKey features:\n")
cat("- Each region gets its own spline coefficients\n")
cat("- Coefficients share hierarchical priors (partial pooling)\n")  
cat("- Regions can have different baseline levels\n")
cat("- Shared noise variance across regions\n")

cat("\nModel structure:\n")
cat("alpha[r] ~ Normal(mu_alpha, tau_alpha)  # Regional spline coefficients\n")
cat("beta[r] ~ Normal(mu_beta, sigma_beta)   # Regional baselines\n")
cat("y ~ Normal(spline(x) + beta[region], sigma)\n")

cat("\nThis allows regions to:\n")
cat("1. Have their own unique nonlinear patterns\n")
cat("2. Borrow strength from each other through the hierarchical prior\n")
cat("3. Maintain region-specific baseline levels\n")

cat("\nExample use cases:\n")
cat("- Disease rates across different regions with similar seasonal patterns\n")
cat("- Economic indicators across countries with shared global trends\n")
cat("- Environmental measurements across sites with common underlying processes\n")

cat("\nTo fit the full model, run:\n")
cat("model <- cmdstan_model('test_regional_splines.stan')\n")
cat("fit <- model$sample(data = stan_data, chains = 4)\n")

cat("\nPlot saved: output/example-regional_data.png\n")