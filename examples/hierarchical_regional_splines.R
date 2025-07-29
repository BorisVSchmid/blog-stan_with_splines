# Hierarchical Regional B-splines with Stan
# Demonstrates decomposition of shared and regional patterns

library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

cat("======================================\n")
cat("Hierarchical Regional Splines Example\n")
cat("======================================\n\n")

cat("This example demonstrates how hierarchical models can:\n")
cat("1. Separate shared global patterns from regional variations\n")
cat("2. Recover both components accurately\n")
cat("3. Provide good fits while borrowing strength across regions\n\n")

# Generate synthetic data with shared components and regional variations
generate_regional_data <- function(n_per_region = 40) {
  # Common components across all regions
  x <- seq(0, 10, length.out = n_per_region)
  
  # 1. Global nonlinear trend (e.g., epidemic curve, economic growth)
  global_trend <- 2 + 0.15 * x - 0.02 * x^2 + 0.001 * x^3
  
  # 2. Seasonal pattern shared across regions (e.g., flu season, shopping patterns)
  seasonality <- 0.8 * sin(2 * pi * x / 3.5)  # Period of 3.5 units
  
  # Combined global pattern
  global_pattern <- global_trend + seasonality
  
  # 3. Region-specific deviations (comparable amplitude to global pattern)
  # Region A: Strong early peak
  regional_effect_A <- 1.2 * exp(-0.5 * (x - 2.5)^2) - 0.3
  
  # Region B: Strong phase shift
  regional_effect_B <- 0.7 * sin(2 * pi * (x + 0.5) / 3.5) - 0.1
  
  # Region C: Strong amplitude modulation
  regional_effect_C <- 0.6 * (1 + 0.1 * x) * sin(2 * pi * x / 3.5) + 0.2
  
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
cat("- Goal: Recover both components from noisy observations\n\n")

# Prepare Stan data
n_per_region <- length(truth$x)
stan_data <- list(
  n_total = nrow(data),
  n_regions = 3,
  x = data$x,
  y = data$y,
  region = data$region,
  num_knots = max(4, min(round(n_per_region), 80)),  # Doubled: n rule instead of n/2
  spline_degree = 3,
  smoothing_strength_global = 0.01,    # Reduced smoothing for more flexibility
  smoothing_strength_regional = 0.01,  # Reduced smoothing for more flexibility
  prior_scale = 2 * sd(data$y)  # Adaptive prior scale
)

# Check if fitted model exists
# Get the correct path whether run from project root or examples directory
if (dir.exists("output")) {
  model_cache_file <- "output/example-hierarchical_model_draws.rds"
} else if (dir.exists("../output")) {
  model_cache_file <- "../output/example-hierarchical_model_draws.rds"
} else {
  stop("Cannot find output directory")
}

if (file.exists(model_cache_file)) {
  cat("Loading cached model draws...\n")
  cache_data <- readRDS(model_cache_file)
  draws <- cache_data$draws
  diagnostics <- cache_data$diagnostics
  summary_stats <- cache_data$summary
  
  # Print cached diagnostics
  cat("\nCached MCMC Diagnostic Summary:\n")
  cat(sprintf("  Divergences: %d\n", diagnostics$num_divergent))
  cat(sprintf("  Max treedepth hits: %d\n", diagnostics$num_max_treedepth))
  cat(sprintf("  Min EBFMI: %.3f\n", diagnostics$min_ebfmi))
} else {
  # Compile and fit model
  cat("Fitting hierarchical model (this may take a few minutes)...\n")
  # Get the correct path whether run from project root or examples directory
  if (file.exists("regional_splines.stan")) {
    stan_file <- "regional_splines.stan"
  } else if (file.exists("examples/regional_splines.stan")) {
    stan_file <- "examples/regional_splines.stan"
  } else {
    stop("Cannot find regional_splines.stan file")
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
    adapt_delta = 0.999,
    max_treedepth = 15
  )
  
  # Print diagnostic summary
  cat("\nMCMC Diagnostic Summary:\n")
  diag_summary <- fit$diagnostic_summary()
  print(diag_summary)
  
  # Extract draws before fit object is destroyed
  draws <- fit$draws(format = "matrix")
  
  # Save draws and diagnostics for future runs
  cat("Saving model draws for future runs... Takes a while.\n")
  cache_data <- list(
    draws = draws,
    diagnostics = list(
      num_divergent = sum(diag_summary$num_divergent),
      num_max_treedepth = sum(diag_summary$num_max_treedepth),
      min_ebfmi = min(diag_summary$ebfmi)
    ),
    summary = list(
      max_rhat = max(fit$summary()$rhat, na.rm = TRUE),
      min_ess_bulk = min(fit$summary()$ess_bulk, na.rm = TRUE)
    )
  )
  saveRDS(cache_data, model_cache_file)
}

# Create output directory
dir.create("output", showWarnings = FALSE)

# Extract global pattern from the model
# The global pattern is mu_alpha * B + mu_beta
# mu_alpha is stored as a row vector, so we need to extract columns like mu_alpha.1, mu_alpha.2, etc.
mu_alpha_cols <- grep("^mu_alpha\\.", colnames(draws), value = TRUE)
if (length(mu_alpha_cols) == 0) {
  # Try alternative format
  mu_alpha_cols <- grep("^mu_alpha\\[", colnames(draws), value = TRUE)
}
if (length(mu_alpha_cols) == 0) {
  # Debug: print column names to understand the format
  cat("\nDEBUG: First 50 column names in draws:\n")
  cat(paste(head(colnames(draws), 50), collapse = ", "), "\n\n")
  stop("Cannot find mu_alpha columns in draws")
}
mu_alpha_mean <- colMeans(draws[, mu_alpha_cols, drop = FALSE])
mu_beta <- mean(draws[, "mu_beta"])

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
# B_plot is n_plot × num_basis, mu_alpha_mean is num_basis × 1
# Extract all mu_alpha samples to compute credible intervals
mu_alpha_samples <- draws[, mu_alpha_cols, drop = FALSE]
mu_beta_samples <- draws[, "mu_beta"]

# Compute global pattern for each MCMC sample
global_pattern_samples <- matrix(0, nrow(draws), length(x_plot))
for (i in 1:nrow(draws)) {
  global_pattern_samples[i,] <- as.numeric(B_plot %*% as.numeric(mu_alpha_samples[i,])) + mu_beta_samples[i]
}

# Get mean and credible intervals
recovered_global <- colMeans(global_pattern_samples)
recovered_global_lower <- apply(global_pattern_samples, 2, quantile, 0.025)
recovered_global_upper <- apply(global_pattern_samples, 2, quantile, 0.975)

# Extract regional deviations
# For each region: (alpha[r] - mu_alpha) * B + (beta[r] - mu_beta)
regional_deviations <- list()
for (r in 1:3) {
  # Alpha_deviation is stored as alpha_deviation[r,j] where r is region and j is basis function
  # The grep needs to escape the brackets and comma properly
  alpha_deviation_cols <- grep(paste0("^alpha_deviation\\[", r, ","), colnames(draws), value = TRUE)
  if (length(alpha_deviation_cols) == 0) {
    # Try alternative format with dots
    alpha_deviation_cols <- grep(paste0("^alpha_deviation\\.", r, "\\."), colnames(draws), value = TRUE)
  }
  if (length(alpha_deviation_cols) == 0) {
    # Debug: print column names to understand the format
    cat(paste0("\nDEBUG: Cannot find alpha_deviation columns for region ", r, "\n"))
    cat("Looking for pattern: alpha_deviation[", r, ",j] or alpha_deviation.", r, ".j\n")
    alpha_dev_like <- grep("^alpha_deviation", colnames(draws), value = TRUE)
    cat("First 20 alpha_deviation columns: ", paste(head(alpha_dev_like, 20), collapse = ", "), "\n\n")
    stop("Cannot find alpha_deviation columns for region ", r)
  }
  # Get all samples for this region's deviations
  alpha_deviation_r_samples <- draws[, alpha_deviation_cols, drop = FALSE]
  beta_r_samples <- draws[, paste0("beta_region[", r, "]")]
  
  # Compute deviation for each MCMC sample
  deviation_samples <- matrix(0, nrow(draws), length(x_plot))
  for (i in 1:nrow(draws)) {
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
       subtitle = "Shared trend and seasonality across all regions (shaded: 95% CI)",
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
  labs(title = "Regional Deviations Recovery",
       subtitle = "Region-specific patterns with 95% CI (dashed: true, solid: recovered)",
       x = "x", y = "Deviation from global") +
  theme_bw()

# 3. Overall fit per region
y_hat_cols <- grep("y_hat\\[", colnames(draws), value = TRUE)
y_hat_mean <- colMeans(draws[, y_hat_cols])
y_hat_lower <- apply(draws[, y_hat_cols], 2, quantile, 0.025)
y_hat_upper <- apply(draws[, y_hat_cols], 2, quantile, 0.975)

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
  labs(title = "Overall Fit Quality",
       subtitle = "Points: data (matching color), Solid: fitted, Dashed: true function, Shaded: 95% CI",
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
    title = "Hierarchical Spline Decomposition",
    subtitle = "How well does the model separate global patterns from regional variations?"
  )

# Save to the correct output directory
output_file <- if (dir.exists("output")) {
  "output/example-hierarchical_decomposition.png"
} else {
  "../output/example-hierarchical_decomposition.png"
}
ggsave(output_file, combined_plot, width = 10, height = 12, dpi = 300)

# Calculate recovery metrics
cat("\n======================================\n")
cat("Recovery Metrics\n")
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
  # Use cached diagnostics
  cat(sprintf("  Divergences: %d\n", diagnostics$num_divergent))
  cat(sprintf("  Max R-hat: %.3f\n", ifelse(exists("summary_stats"), summary_stats$max_rhat, NA)))
  cat(sprintf("  Min ESS Bulk: %.0f\n", ifelse(exists("summary_stats"), summary_stats$min_ess_bulk, NA)))
} else {
  # Calculate from fit object if available
  cat(sprintf("  Divergences: %d\n", sum(fit$sampler_diagnostics()[,,"divergent__"])))
  cat(sprintf("  Max R-hat: %.3f\n", max(fit$summary()$rhat, na.rm = TRUE)))
  cat(sprintf("  Min ESS Bulk: %.0f\n", min(fit$summary()$ess_bulk, na.rm = TRUE)))
}

cat("\nOutput saved: output/example-hierarchical_decomposition.png\n")

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

# Global pattern diagnostics
diag_global <- compute_spline_diagnostics(truth$global_pattern, recovered_global, "Global")
diagnostics_list[[1]] <- diag_global

# Regional deviation diagnostics
for (i in 1:3) {
  true_dev <- switch(i, truth$regional_effects$A, truth$regional_effects$B, truth$regional_effects$C)
  recovered_dev <- regional_deviations[[i]]$mean
  label <- paste("Region", c("A", "B", "C")[i], "Deviation")
  diag_regional <- compute_spline_diagnostics(true_dev, recovered_dev, label)
  diagnostics_list[[i + 1]] <- diag_regional
}

# Combine diagnostics
diagnostics_df <- do.call(rbind, diagnostics_list)

# Save to CSV
csv_file <- if (dir.exists("output")) {
  "output/example-hierarchical_spline_diagnostics.csv"
} else {
  "../output/example-hierarchical_spline_diagnostics.csv"
}
write.csv(diagnostics_df, csv_file, row.names = FALSE)
cat("Spline diagnostics saved to:", csv_file, "\n")

# Create diagnostic plots for all splines
cat("\nGenerating diagnostic plots for all splines...\n")

# Function to create diagnostic plots for a single spline
create_spline_diagnostic_plots <- function(y_true, y_fitted, x_vals, spline_name) {
  residuals <- y_true - y_fitted
  
  # 1. Fitted vs True
  df_fit <- data.frame(x = x_vals, true = y_true, fitted = y_fitted)
  p_fit <- ggplot(df_fit, aes(x = x)) +
    geom_line(aes(y = true), color = "black", linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = fitted), color = "blue", linewidth = 1) +
    labs(title = paste(spline_name, "- Fit"),
         subtitle = "Dashed: True, Solid: Fitted",
         x = "x", y = "y") +
    theme_bw()
  
  # 2. Residuals vs x
  df_resid <- data.frame(x = x_vals, residuals = residuals)
  p_resid <- ggplot(df_resid, aes(x = x, y = residuals)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(title = paste(spline_name, "- Residuals"),
         subtitle = "Loess smooth with 95% CI",
         x = "x", y = "Residuals") +
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
  
  return(list(fit = p_fit, resid = p_resid, acf = p_acf))
}

# Create plots for each spline
plot_list <- list()

# 1. Global pattern
plots_global <- create_spline_diagnostic_plots(
  truth$global_pattern, recovered_global, x_plot, "Global Pattern"
)
plot_list[[1]] <- plots_global

# 2. Regional deviations
for (i in 1:3) {
  true_dev <- switch(i, truth$regional_effects$A, 
                     truth$regional_effects$B, 
                     truth$regional_effects$C)
  recovered_dev <- regional_deviations[[i]]$mean
  label <- paste("Region", c("A", "B", "C")[i], "Deviation")
  
  plots_regional <- create_spline_diagnostic_plots(
    true_dev, recovered_dev, x_plot, label
  )
  plot_list[[i + 1]] <- plots_regional
}

# Combine all plots into a single figure
library(patchwork)
combined_diagnostics <- 
  (plot_list[[1]]$fit | plot_list[[1]]$resid | plot_list[[1]]$acf) /
  (plot_list[[2]]$fit | plot_list[[2]]$resid | plot_list[[2]]$acf) /
  (plot_list[[3]]$fit | plot_list[[3]]$resid | plot_list[[3]]$acf) /
  (plot_list[[4]]$fit | plot_list[[4]]$resid | plot_list[[4]]$acf) +
  plot_annotation(
    title = "Hierarchical Spline Diagnostics",
    subtitle = "Fit quality, residual patterns, and autocorrelation for all components"
  )

# Save diagnostic plots
diag_plot_file <- if (dir.exists("output")) {
  "output/example-hierarchical_spline_diagnostics.png"
} else {
  "../output/example-hierarchical_spline_diagnostics.png"
}
ggsave(diag_plot_file, combined_diagnostics, width = 15, height = 16, dpi = 300)
cat("Diagnostic plots saved to:", diag_plot_file, "\n")
