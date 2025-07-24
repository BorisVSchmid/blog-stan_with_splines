# Test regional B-splines with hierarchical priors
# Demonstrates fitting splines to multiple regions that share information

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "tidyr", "patchwork")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

set.seed(42)

# Generate synthetic data for 3 regions
generate_regional_data <- function(n_per_region = 50) {
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
data <- generate_regional_data(50)

# Prepare Stan data
stan_data <- list(
  n_total = nrow(data),
  n_regions = 3,
  x = data$x,
  y = data$y,
  region = data$region,
  num_knots = 8,
  spline_degree = 3,
  smoothing_strength = 1.0,  # Mild smoothing
  prior_scale = 2 * sd(data$y)  # Data-driven prior
)

# Compile and fit model
cat("Compiling regional splines model...\n")
model <- cmdstan_model("tests/test_regional_splines.stan")

cat("Fitting model with hierarchical priors...\n")
fit <- model$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 500,
  seed = 123,
  adapt_delta = 0.98,  # Higher to reduce divergences
  max_treedepth = 15   # Increase to avoid hitting treedepth limit
)

# Print summary of key parameters
cat("\nParameter summary:\n")
fit$summary(c("mu_alpha", "tau_alpha", "beta_region", "sigma", "icc")) %>%
  print(n = 20)

# Extract draws
draws <- fit$draws(format = "matrix")

# Extract predictions for each region
plot_data <- list()


for (r in 1:3) {
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  # Stan outputs array[n_regions] vector[1000] as y_plot[r,i] format
  # Use regex pattern to match y_plot[r,*] where * is any number
  y_plot_cols <- grep(paste0("^y_plot\\[", r, ",\\d+\\]$"), colnames(draws), value = TRUE)
  
  if (length(y_plot_cols) == 0) {
    stop(paste("ERROR: No y_plot columns found for region", r))
  }
  
  y_plot_mean <- colMeans(draws[, y_plot_cols])
  y_plot_lower <- apply(draws[, y_plot_cols], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, y_plot_cols], 2, quantile, 0.975)
  
  plot_data[[r]] <- data.frame(
    x = x_plot,
    y_mean = y_plot_mean,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper,
    region = r,
    region_name = c("Region A", "Region B", "Region C")[r]
  )
}

plot_df <- bind_rows(plot_data)

# Create plots for each region
p1 <- ggplot() +
  # Fitted splines with uncertainty
  geom_ribbon(data = plot_df, 
              aes(x = x, ymin = y_lower, ymax = y_upper, fill = region_name), 
              alpha = 0.3) +
  geom_line(data = plot_df, 
            aes(x = x, y = y_mean, color = region_name), 
            linewidth = 1.2) +
  # Original data
  geom_point(data = data, 
             aes(x = x, y = y, color = region_name), 
             alpha = 0.6, size = 1.5) +
  # True functions
  geom_line(data = data, 
            aes(x = x, y = y_true, group = region_name), 
            linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ region_name, scales = "free_y", ncol = 1) +
  labs(title = "Regional Splines with Hierarchical Priors",
       subtitle = "Each region has its own spline while sharing information through hierarchical structure",
       x = "x", y = "y") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

# Extract and plot the global mean spline coefficients
mu_alpha_cols <- grep("mu_alpha\\[", colnames(draws), value = TRUE)
mu_alpha_mean <- colMeans(draws[, mu_alpha_cols])
mu_alpha_sd <- apply(draws[, mu_alpha_cols], 2, sd)

p2 <- data.frame(
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
  theme_minimal()

# Plot variance components
var_between <- mean(draws[, "var_between"])
var_within <- mean(draws[, "var_within"])
icc <- mean(draws[, "icc"])

var_df <- data.frame(
  component = c("Between Regions", "Within Regions"),
  variance = c(var_between, var_within),
  percentage = c(var_between, var_within) / (var_between + var_within) * 100
)

p3 <- ggplot(var_df, aes(x = component, y = variance, fill = component)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            vjust = -0.5, size = 4) +
  labs(title = "Variance Components",
       subtitle = sprintf("ICC = %.3f (%.1f%% of variance is between regions)", 
                         icc, icc * 100),
       x = "", y = "Variance") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save outputs
dir.create("output", showWarnings = FALSE)

# Combine plots using patchwork
combined_plot <- (p1 | (p2 / p3)) + 
  plot_layout(widths = c(2, 1))

# Save combined visualization
ggsave("output/test-regional_splines_comprehensive.png", combined_plot, 
       width = 14, height = 10, dpi = 300)

# Also save individual plots for detailed examination
ggsave("output/test-regional_splines_fits.png", p1, 
       width = 8, height = 10, dpi = 300)

ggsave("output/test-regional_splines_coefficients.png", p2, 
       width = 6, height = 5, dpi = 300)

ggsave("output/test-regional_splines_variance.png", p3, 
       width = 6, height = 5, dpi = 300)

cat("\nOutputs saved:\n")
cat("- output/test-regional_splines_comprehensive.png (combined visualization)\n")
cat("- output/test-regional_splines_fits.png (individual region fits)\n")
cat("- output/test-regional_splines_coefficients.png (global coefficients)\n")
cat("- output/test-regional_splines_variance.png (variance decomposition)\n")

# Print regional effects
beta_region <- colMeans(draws[, grep("beta_region\\[", colnames(draws), value = TRUE)])
cat("\nRegional baseline effects:\n")
for (i in 1:3) {
  cat(sprintf("  Region %s: %.3f\n", c("A", "B", "C")[i], beta_region[i]))
}

# Print smoothness parameters
lambda_means <- colMeans(draws[, grep("lambda\\[", colnames(draws), value = TRUE)])
cat("\nSmoothness parameters (lambda):\n")
for (i in 1:3) {
  cat(sprintf("  Region %s: %.3f\n", c("A", "B", "C")[i], lambda_means[i]))
}

# Additional diagnostics from the full version
# Extract basis function values for detailed analysis
cat("\n\nDetailed Basis Function Analysis:\n")
cat("==================================\n")

for (r in 1:3) {
  # Extract basis function values for this region
  basis_cols <- grep(paste0("basis_functions\\[\\d+,", r, "\\]"), 
                     colnames(draws), value = TRUE)
  
  if (length(basis_cols) > 0) {
    cat(sprintf("\nRegion %s:\n", c("A", "B", "C")[r]))
    cat(sprintf("  Number of basis functions: %d\n", length(basis_cols)))
    
    # Check sparsity (how many basis functions are effectively zero)
    basis_means <- colMeans(draws[, basis_cols])
    near_zero <- sum(abs(basis_means) < 0.01)
    cat(sprintf("  Basis functions near zero: %d (%.1f%%)\n", 
                near_zero, 100 * near_zero / length(basis_cols)))
  }
}

# Model diagnostics
cat("\n\nModel Diagnostics:\n")
cat("==================\n")
diagnostics <- fit$diagnostic_summary()
cat("Divergences:", sum(diagnostics$num_divergent), "\n")
cat("Max treedepth exceeded:", sum(diagnostics$num_max_treedepth), "\n")
cat("E-BFMI:", round(min(diagnostics$ebfmi), 3), "\n")

# Print effective sample sizes for key parameters
ess <- fit$summary(c("mu_alpha", "tau_alpha", "sigma"))$ess_bulk
cat("\nEffective sample sizes (bulk):\n")
cat("  mu_alpha:", round(min(ess[grep("mu_alpha", names(ess))]), 0), "-",
    round(max(ess[grep("mu_alpha", names(ess))]), 0), "\n")
cat("  tau_alpha:", round(min(ess[grep("tau_alpha", names(ess))]), 0), "-",
    round(max(ess[grep("tau_alpha", names(ess))]), 0), "\n")
cat("  sigma:", round(ess["sigma"], 0), "\n")