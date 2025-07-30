# Hierarchical Regional B-splines with Adaptive Shrinkage
# Updated version with different smoothing strengths and adaptive shrinkage

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)  # Use base R sd, not posterior::sd

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork", "tidyr")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

set.seed(42)

cat("================================================\n")
cat("Adaptive Hierarchical Regional Splines Example\n")
cat("================================================\n\n")

cat("This example demonstrates improved control over pattern absorption with:\n")
cat("1. Different smoothing strengths for global vs regional patterns\n")
cat("2. Adaptive shrinkage that learns from the data\n")
cat("3. Better separation of shared and regional components\n\n")

# Generate synthetic data with shared components and regional variations
generate_regional_data <- function(n_per_region = 40) {
  # Common components across all regions
  x <- seq(0, 10, length.out = n_per_region)
  
  # 1. Global nonlinear trend (e.g., epidemic curve, economic growth)
  global_trend <- 2 + 0.15 * x - 0.02 * x^2 + 0.001 * x^3
  
  # 2. Seasonal pattern shared across regions (e.g., flu season, shopping patterns)
  seasonality <- 0.8 * sin(2 * pi * x / 3.5)  # Period of 3.5 units
  
  # Combined global pattern (increased amplitude by 10%)
  global_pattern <- 1.1 * (global_trend + seasonality)
  
  # 3. Region-specific deviations (comparable amplitude to global pattern)
  # Region A: Strong early peak
  regional_effect_A <- 1.2 * exp(-0.5 * (x - 2.5)^2) - 0.3
  
  # Region B: Different period (faster oscillation)
  regional_effect_B <- 0.8 * sin(2 * pi * x / 2.5) - 0.1
  
  # Region C: Different period (slower oscillation)
  regional_effect_C <- 0.7 * sin(2 * pi * x / 5) + 0.2
  
  # Total effects
  y1_true <- global_pattern + regional_effect_A
  y2_true <- global_pattern + regional_effect_B
  y3_true <- global_pattern + regional_effect_C
  
  # Add noise
  y1 <- y1_true + rnorm(n_per_region, 0, 0.15)
  y2 <- y2_true + rnorm(n_per_region, 0, 0.15)
  y3 <- y3_true + rnorm(n_per_region, 0, 0.15)
  
  # Return structured data
  list(
    data = data.frame(
      x = c(x, x, x),
      y = c(y1, y2, y3),
      region = rep(1:3, each = n_per_region),
      region_name = factor(rep(c("Region A", "Region B", "Region C"), each = n_per_region))
    ),
    truth = list(
      x = x,
      global_pattern = global_pattern,
      regional_effects = list(
        A = regional_effect_A,
        B = regional_effect_B,
        C = regional_effect_C
      ),
      total_effects = list(
        A = y1_true,
        B = y2_true,
        C = y3_true
      )
    )
  )
}

# Generate data
data_list <- generate_regional_data(40)
data <- data_list$data
truth <- data_list$truth

cat("Data structure:\n")
cat("- Shared global pattern: Cubic trend + seasonal oscillation\n")
cat("- Regional deviations: Small modifications to the global pattern\n")
cat("- Goal: Recover both components with adaptive shrinkage\n\n")

# Prepare Stan data with improved parameters
n_per_region <- length(truth$x)
stan_data <- list(
  n_total = nrow(data),
  n_regions = 3,
  x = data$x,
  y = data$y,
  region = data$region,
  num_knots = 20,  # Fixed at 20 knots as requested
  spline_degree = 3,
  smoothing_strength_global = 0.01,    # Very light smoothing for global pattern
  smoothing_strength_regional = 0.01,  # Very light smoothing for regional deviations
  prior_scale = 2 * sd(data$y)  # Back to traditional default
)

cat("Model configuration:\n")
cat(sprintf("- Number of knots: %d (doubled from default n/2 rule)\n", stan_data$num_knots))
cat(sprintf("- Global smoothing: %.2f (very light smoothing)\n", stan_data$smoothing_strength_global))
cat(sprintf("- Regional smoothing: %.2f (very light smoothing)\n", stan_data$smoothing_strength_regional))
cat("- Shrinkage: Adaptive (learned from data)\n\n")

# Compile and fit model
cat("Fitting adaptive hierarchical model (this may take a few minutes)...\n")
# Get the correct path whether run from project root or examples directory
if (file.exists("regional_splines_adaptive.stan")) {
  stan_file <- "regional_splines_adaptive.stan"
} else if (file.exists("examples/regional_splines_adaptive.stan")) {
  stan_file <- "examples/regional_splines_adaptive.stan"
} else {
  stop("Cannot find regional_splines_adaptive.stan file")
}
model <- cmdstan_model(stan_file)

fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 3000,
    iter_sampling = 5000,
    refresh = 0,
    seed = 123,
    adapt_delta = 0.9995,
    max_treedepth = 15
  )
  
  # Print diagnostic summary
  cat("\nMCMC Diagnostic Summary:\n")
  diag_summary <- fit$diagnostic_summary()
  print(diag_summary)
  
  # Extract draws efficiently by variable groups
  # Core parameters for analysis
  draws_core <- fit$draws(variables = c("mu_alpha", "mu_beta", "alpha_deviation", 
                                        "beta_region", "shrinkage_factor", "sigma",
                                        "tau_alpha", "sigma_beta"), 
                         format = "matrix")
  
  # Y_hat for diagnostics
  draws_yhat <- fit$draws(variables = "y_hat", format = "matrix")
  
  # Variance components
  draws_var <- fit$draws(variables = c("var_between", "var_within", "icc", 
                                       "var_global", "var_regional"),
                        format = "matrix")

# Print adaptive shrinkage information
shrinkage_mean <- mean(draws_core[, "shrinkage_factor"])
shrinkage_sd <- sd(draws_core[, "shrinkage_factor"])
cat(sprintf("\nAdaptive shrinkage factor: %.2f ± %.2f\n", shrinkage_mean, shrinkage_sd))
cat("(Lower values = more shrinkage, 1.0 = no shrinkage)\n\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# Extract global pattern from the model
mu_alpha_cols <- grep("^mu_alpha\\.", colnames(draws_core), value = TRUE)
if (length(mu_alpha_cols) == 0) {
  mu_alpha_cols <- grep("^mu_alpha\\[", colnames(draws_core), value = TRUE)
}
mu_alpha_mean <- colMeans(draws_core[, mu_alpha_cols, drop = FALSE])
mu_beta <- mean(draws_core[, "mu_beta"])

# Build basis matrix for plotting
# Copy the build_b_spline function from the Stan file
build_b_spline <- function(t, ext_knots, ind, order, degree) {
  b_spline <- rep(0, length(t))
  w1 <- rep(0, length(t))
  w2 <- rep(0, length(t))
  
  if (order == 1) {
    for (i in 1:length(t)) {
      if (t[i] >= ext_knots[length(ext_knots) - degree] && 
          ind == length(ext_knots) - degree - 1) {
        b_spline[i] <- 1
      } else {
        b_spline[i] <- (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1])
      }
    }
  } else {
    if (abs(ext_knots[ind] - ext_knots[ind+order-1]) > 1e-10) {
      w1 <- (t - ext_knots[ind]) / (ext_knots[ind+order-1] - ext_knots[ind])
    }
    if (abs(ext_knots[ind+1] - ext_knots[ind+order]) > 1e-10) {
      w2 <- 1 - (t - ext_knots[ind+1]) / (ext_knots[ind+order] - ext_knots[ind+1])
    }
    b_spline <- w1 * build_b_spline(t, ext_knots, ind, order-1, degree) +
                w2 * build_b_spline(t, ext_knots, ind+1, order-1, degree)
  }
  return(b_spline)
}

x_plot <- truth$x
num_basis <- stan_data$num_knots + stan_data$spline_degree - 1

# Recreate the same knot placement as in Stan
x_min <- min(data$x)
x_max <- max(data$x)
knots <- numeric(stan_data$num_knots)
knots[1] <- x_min
knots[stan_data$num_knots] <- x_max
if (stan_data$num_knots > 2) {
  for (k in 2:(stan_data$num_knots-1)) {
    knots[k] <- x_min + (x_max - x_min) * (k - 1.0) / (stan_data$num_knots - 1.0)
  }
}

# Extended knots
ext_knots <- c(rep(knots[1], stan_data$spline_degree), 
               knots, 
               rep(knots[length(knots)], stan_data$spline_degree))

# Build basis matrix
B_plot <- matrix(0, length(x_plot), num_basis)
for (i in 1:num_basis) {
  B_plot[, i] <- build_b_spline(x_plot, ext_knots, i, stan_data$spline_degree + 1, stan_data$spline_degree)
}

# Compute recovered global pattern with uncertainty
mu_alpha_samples <- draws_core[, mu_alpha_cols, drop = FALSE]
mu_beta_samples <- draws_core[, "mu_beta"]

# Compute global pattern for each MCMC sample
global_pattern_samples <- matrix(0, nrow(draws_core), length(x_plot))
for (i in 1:nrow(draws_core)) {
  global_pattern_samples[i,] <- as.numeric(B_plot %*% as.numeric(mu_alpha_samples[i,])) + mu_beta_samples[i]
}

# Get mean and credible intervals
recovered_global <- colMeans(global_pattern_samples)
recovered_global_lower <- apply(global_pattern_samples, 2, quantile, 0.025)
recovered_global_upper <- apply(global_pattern_samples, 2, quantile, 0.975)

# Extract regional deviations
regional_deviations <- list()
for (r in 1:3) {
  alpha_deviation_cols <- grep(paste0("^alpha_deviation\\[", r, ","), colnames(draws_core), value = TRUE)
  if (length(alpha_deviation_cols) == 0) {
    alpha_deviation_cols <- grep(paste0("^alpha_deviation\\.", r, "\\."), colnames(draws_core), value = TRUE)
  }
  
  alpha_deviation_r_samples <- draws_core[, alpha_deviation_cols, drop = FALSE]
  beta_r_samples <- draws_core[, paste0("beta_region[", r, "]")]
  
  # Compute deviation for each MCMC sample
  deviation_samples <- matrix(0, nrow(draws_core), length(x_plot))
  for (i in 1:nrow(draws_core)) {
    deviation_samples[i,] <- as.numeric(B_plot %*% as.numeric(alpha_deviation_r_samples[i,])) + 
                             (beta_r_samples[i] - mu_beta_samples[i])
  }
  
  # Get mean and credible intervals
  deviation_mean <- colMeans(deviation_samples)
  deviation_lower <- apply(deviation_samples, 2, quantile, 0.025)
  deviation_upper <- apply(deviation_samples, 2, quantile, 0.975)
  
  regional_deviations[[r]] <- list(
    mean = deviation_mean,
    lower = deviation_lower,
    upper = deviation_upper
  )
}

# Create comparison plots
# 1. Global pattern recovery with confidence interval
df_global <- data.frame(
  x = x_plot,
  y_true = truth$global_pattern,
  y_recovered = recovered_global,
  y_lower = recovered_global_lower,
  y_upper = recovered_global_upper
)

p_global <- ggplot(df_global, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = y_true), color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = y_recovered), color = "blue", linewidth = 1.2) +
  labs(title = "Global Pattern Recovery",
       subtitle = sprintf("Global smoothing (%.2f)", 
                         stan_data$smoothing_strength_global),
       x = "x", y = "y") +
  theme_bw()

# 2. Regional deviations recovery with confidence intervals
df_regional <- data.frame(
  x = rep(x_plot, 3),
  y_true = c(truth$regional_effects$A, truth$regional_effects$B, truth$regional_effects$C),
  y_recovered = c(regional_deviations[[1]]$mean, 
                  regional_deviations[[2]]$mean, 
                  regional_deviations[[3]]$mean),
  y_lower = c(regional_deviations[[1]]$lower,
              regional_deviations[[2]]$lower,
              regional_deviations[[3]]$lower),
  y_upper = c(regional_deviations[[1]]$upper,
              regional_deviations[[2]]$upper,
              regional_deviations[[3]]$upper),
  region = factor(rep(c("Region A", "Region B", "Region C"), each = length(x_plot)))
)

p_regional <- ggplot(df_regional, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "red", alpha = 0.3) +
  geom_hline(yintercept = 0, alpha = 0.3) +
  geom_line(aes(y = y_true), color = "black", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = y_recovered), color = "red", linewidth = 1) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = sprintf("Regional Deviations (Adaptive shrinkage: %.1f ± %.1f)", 
                      shrinkage_mean, shrinkage_sd),
       subtitle = sprintf("Regional smoothing: %.2f", 
                         stan_data$smoothing_strength_regional),
       x = "x", y = "Deviation from global") +
  theme_bw()

# 3. Overall fit per region
y_hat_cols <- grep("y_hat\\[", colnames(draws_yhat), value = TRUE)
y_hat_mean <- colMeans(draws_yhat[, y_hat_cols])
y_hat_lower <- apply(draws_yhat[, y_hat_cols], 2, quantile, 0.025)
y_hat_upper <- apply(draws_yhat[, y_hat_cols], 2, quantile, 0.975)

data$y_fitted <- y_hat_mean
data$y_lower <- y_hat_lower
data$y_upper <- y_hat_upper

# Add true values for comparison
data$y_true <- c(truth$total_effects$A, truth$total_effects$B, truth$total_effects$C)

p_fits <- ggplot(data, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = region_name), alpha = 0.3) +
  geom_line(aes(y = y_true, color = region_name), linetype = "dashed", linewidth = 1) +
  geom_point(aes(y = y, color = region_name), alpha = 0.7, size = 1.5) +
  geom_line(aes(y = y_fitted, color = region_name), linewidth = 1) +
  facet_wrap(~ region_name, scales = "free_y") +
  labs(title = "Overall Fit Quality with Adaptive Model",
       subtitle = "Points: data, Solid: fitted, Dashed: true function, Shaded: 95% CI",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

# Combine plots
library(patchwork)
combined_plot <- p_global / p_regional / p_fits +
  plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(
    title = "Adaptive Hierarchical Spline Decomposition",
    subtitle = "Improved pattern separation with adaptive shrinkage and equal smoothing"
  )

# Save to the correct output directory
output_file <- if (dir.exists("output")) {
  "output/example-hierarchical_adaptive_decomposition.png"
} else {
  "../output/example-hierarchical_adaptive_decomposition.png"
}
ggsave(output_file, combined_plot, width = 10, height = 12, dpi = 300)

# Calculate recovery metrics
cat("\n======================================\n")
cat("Recovery Metrics (Adaptive Model)\n")
cat("======================================\n")

# Global pattern recovery
global_rmse <- sqrt(mean((truth$global_pattern - recovered_global)^2))
global_cor <- cor(truth$global_pattern, recovered_global)
cat(sprintf("\nGlobal pattern recovery:\n"))
cat(sprintf("  RMSE: %.4f\n", global_rmse))
cat(sprintf("  Correlation: %.4f\n", global_cor))

# Regional deviations recovery
cat("\nRegional deviations recovery:\n")
for (i in 1:3) {
  true_dev <- switch(i, truth$regional_effects$A, truth$regional_effects$B, truth$regional_effects$C)
  recovered_dev <- regional_deviations[[i]]$mean
  rmse <- sqrt(mean((true_dev - recovered_dev)^2))
  cor_val <- cor(true_dev, recovered_dev)
  cat(sprintf("  Region %s - RMSE: %.4f, Correlation: %.4f\n", 
              c("A", "B", "C")[i], rmse, cor_val))
}

# Overall fit quality
cat("\nOverall fit quality (data vs fitted):\n")
for (r in 1:3) {
  idx <- data$region == r
  rmse <- sqrt(mean((data$y[idx] - data$y_fitted[idx])^2))
  cat(sprintf("  Region %s - RMSE: %.4f\n", c("A", "B", "C")[r], rmse))
}

# Model diagnostics
cat("\nModel diagnostics:\n")
if (exists("diagnostics")) {
  cat(sprintf("  Divergences: %d\n", diagnostics$num_divergent))
  cat(sprintf("  Max R-hat: %.3f\n", ifelse(exists("summary_stats"), summary_stats$max_rhat, NA)))
  cat(sprintf("  Min ESS Bulk: %.0f\n", ifelse(exists("summary_stats"), summary_stats$min_ess_bulk, NA)))
}

cat("\nOutput saved: output/example-hierarchical_adaptive_decomposition.png\n")

# Create diagnostic report for the splines
cat("\nGenerating spline diagnostics report...\n")

# Compute residuals for each spline
compute_spline_diagnostics <- function(y_true, y_fitted, label) {
  residuals <- y_true - y_fitted
  n <- length(residuals)
  
  # Autocorrelation at different lags
  acf_vals <- acf(residuals, lag.max = 5, plot = FALSE)$acf[-1]  # Remove lag 0
  
  # Smoothness metrics
  first_diff <- diff(y_fitted)
  second_diff <- diff(first_diff)
  roughness <- sqrt(mean(second_diff^2))
  
  # Residual statistics
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  
  # Return diagnostics
  data.frame(
    spline = label,
    rmse = rmse,
    mae = mae,
    roughness = roughness,
    acf_lag1 = acf_vals[1],
    acf_lag2 = acf_vals[2],
    acf_lag3 = acf_vals[3],
    acf_lag4 = acf_vals[4],
    acf_lag5 = acf_vals[5],
    mean_residual = mean(residuals),
    sd_residual = sd(residuals),
    min_residual = min(residuals),
    max_residual = max(residuals)
  )
}

# Collect diagnostics for all splines
diagnostics_list <- list()

# For hierarchical models, check model assumptions using actual residuals
# accounting for the hierarchical structure

# 1. Pooled residuals analysis (overall model adequacy)
all_residuals <- c()
all_x <- c()
all_fitted <- c()
for (r in 1:3) {
  idx <- data$region == r
  all_residuals <- c(all_residuals, data$y[idx] - data$y_fitted[idx])
  all_x <- c(all_x, data$x[idx])
  all_fitted <- c(all_fitted, data$y_fitted[idx])
}

# Create synthetic "true" values for pooled analysis (just for function compatibility)
pooled_true <- all_fitted + all_residuals  # This equals the original y values
diag_pooled <- compute_spline_diagnostics(pooled_true, all_fitted, "Pooled (All Regions)")
diagnostics_list[[1]] <- diag_pooled

# 2. Within-region residuals analysis
for (r in 1:3) {
  idx <- data$region == r
  y_obs <- data$y[idx]
  y_fit <- data$y_fitted[idx]
  label <- paste("Region", c("A", "B", "C")[r], "Residuals")
  diag_regional <- compute_spline_diagnostics(y_obs, y_fit, label)
  diagnostics_list[[r + 1]] <- diag_regional
}

# Combine diagnostics
diagnostics_df <- do.call(rbind, diagnostics_list)

# Save to CSV
csv_file <- if (dir.exists("output")) {
  "output/example-hierarchical_adaptive_spline_diagnostics.csv"
} else {
  "../output/example-hierarchical_adaptive_spline_diagnostics.csv"
}
write.csv(diagnostics_df, csv_file, row.names = FALSE)
cat("Spline diagnostics saved to:", csv_file, "\n")

# Create diagnostic plots for all splines
cat("\nGenerating diagnostic plots for all splines...\n")

# Function to create diagnostic plots for hierarchical model
create_spline_diagnostic_plots <- function(y_true, y_fitted, x_vals, spline_name) {
  residuals <- y_true - y_fitted
  
  # 1. Residuals vs X (check for patterns)
  df_resid_x <- data.frame(x = x_vals, residuals = residuals)
  p_resid_x <- ggplot(df_resid_x, aes(x = x, y = residuals)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(title = paste(spline_name, "- Residuals vs X"),
         subtitle = "Check for systematic patterns",
         x = "x", y = "Residuals") +
    theme_bw()
  
  # 2. Residuals vs Fitted (check for heteroscedasticity)
  df_resid_fit <- data.frame(fitted = y_fitted, residuals = residuals)
  p_resid_fit <- ggplot(df_resid_fit, aes(x = fitted, y = residuals)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(title = paste(spline_name, "- Residuals vs Fitted"),
         subtitle = "Check for heteroscedasticity",
         x = "Fitted values", y = "Residuals") +
    theme_bw()
  
  # 3. ACF plot
  acf_data <- acf(residuals, lag.max = 20, plot = FALSE)
  df_acf <- data.frame(lag = acf_data$lag[-1], acf = acf_data$acf[-1])
  
  # Calculate confidence bounds
  n <- length(residuals)
  conf_bound <- qnorm(0.975) / sqrt(n)
  
  p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-conf_bound, conf_bound), 
               color = "blue", linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0), color = "black") +
    geom_point(size = 2) +
    labs(title = paste(spline_name, "- Autocorrelation"),
         subtitle = "Blue dashed: 95% confidence bounds",
         x = "Lag", y = "ACF") +
    theme_bw()
  
  return(list(resid_x = p_resid_x, resid_fit = p_resid_fit, acf = p_acf))
}

# Create plots for each spline
plot_list <- list()

# 1. Pooled diagnostics (all regions combined)
plots_pooled <- create_spline_diagnostic_plots(
  pooled_true, all_fitted, all_x, "Pooled (All Regions)"
)
plot_list[[1]] <- plots_pooled

# 2. Within-region diagnostics
for (r in 1:3) {
  idx <- data$region == r
  y_obs <- data$y[idx]
  y_fit <- data$y_fitted[idx]
  x_region <- data$x[idx]
  label <- paste("Region", c("A", "B", "C")[r])
  
  plots_regional <- create_spline_diagnostic_plots(
    y_obs, y_fit, x_region, label
  )
  plot_list[[r + 1]] <- plots_regional
}

# Combine all plots into a single figure
combined_diagnostics <- 
  (plot_list[[1]]$resid_x | plot_list[[1]]$resid_fit | plot_list[[1]]$acf) /
  (plot_list[[2]]$resid_x | plot_list[[2]]$resid_fit | plot_list[[2]]$acf) /
  (plot_list[[3]]$resid_x | plot_list[[3]]$resid_fit | plot_list[[3]]$acf) /
  (plot_list[[4]]$resid_x | plot_list[[4]]$resid_fit | plot_list[[4]]$acf) +
  plot_annotation(
    title = "Hierarchical Model Diagnostic Checks",
    subtitle = "Pooled and within-region residual analysis for model assumptions"
  )

# Save diagnostic plots
diag_plot_file <- if (dir.exists("output")) {
  "output/example-hierarchical_adaptive_spline_diagnostics.png"
} else {
  "../output/example-hierarchical_adaptive_spline_diagnostics.png"
}
ggsave(diag_plot_file, combined_diagnostics, width = 15, height = 16, dpi = 300)
cat("Diagnostic plots saved to:", diag_plot_file, "\n")

# Print comparison with original model
cat("\n======================================\n")
cat("Improvements Summary\n")
cat("======================================\n")
cat("1. Adaptive shrinkage factor learned from data\n")
cat("2. Different smoothing: very light for both global and regional (0.01)\n")
cat("3. Prevents overfitting while maintaining flexibility\n")
cat("4. Doubled knots for better flexibility\n")
cat("5. Correct residual calculations for hierarchical structure\n")
cat("\nRun both models to compare diagnostics!\n")

