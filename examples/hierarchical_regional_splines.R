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
  
  # 3. Region-specific deviations (smaller to emphasize shared pattern)
  # Region A: Small early peak
  regional_effect_A <- 0.3 * exp(-0.5 * (x - 2.5)^2) - 0.1
  
  # Region B: Small phase shift
  regional_effect_B <- 0.2 * sin(2 * pi * (x + 0.5) / 3.5) - 0.05
  
  # Region C: Small amplitude modulation
  regional_effect_C <- 0.15 * (1 + 0.1 * x) * sin(2 * pi * x / 3.5) + 0.1
  
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
  num_knots = max(4, min(round(n_per_region/2), 40)),  # n/2 rule
  spline_degree = 3,
  smoothing_strength = 1,        # Default smoothing
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
    adapt_delta = 0.995,
    max_treedepth = 12
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
mu_alpha_cols <- grep("mu_alpha\\[", colnames(draws), value = TRUE)
mu_alpha_mean <- colMeans(draws[, mu_alpha_cols])
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

# Compute recovered global pattern
# B_plot is n_plot × num_basis, mu_alpha_mean is num_basis × 1
recovered_global <- as.numeric(B_plot %*% mu_alpha_mean) + mu_beta

# Extract regional deviations
# For each region: (alpha[r] - mu_alpha) * B + (beta[r] - mu_beta)
regional_deviations <- list()
for (r in 1:3) {
  alpha_cols <- grep(paste0("alpha\\[", r, ","), colnames(draws), value = TRUE, fixed = TRUE)
  alpha_r_mean <- colMeans(draws[, alpha_cols])
  beta_r <- mean(draws[, paste0("beta_region[", r, "]")])
  
  deviation <- as.numeric(B_plot %*% (alpha_r_mean - mu_alpha_mean)) + (beta_r - mu_beta)
  regional_deviations[[r]] <- deviation
}

# Create comparison plots
# 1. Global pattern recovery
df_global <- data.frame(
  x = rep(x_plot, 2),
  y = c(truth$global_pattern, recovered_global),
  type = factor(rep(c("True", "Recovered"), each = length(x_plot)))
)

p_global <- ggplot(df_global, aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("True" = "black", "Recovered" = "blue")) +
  scale_linetype_manual(values = c("True" = "dashed", "Recovered" = "solid")) +
  labs(title = "Global Pattern Recovery",
       subtitle = "Shared trend and seasonality across all regions",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.title = element_blank())

# 2. Regional deviations recovery
df_regional <- data.frame(
  x = rep(x_plot, 6),
  y = c(truth$regional_effects$A, regional_deviations[[1]],
        truth$regional_effects$B, regional_deviations[[2]], 
        truth$regional_effects$C, regional_deviations[[3]]),
  type = factor(rep(c("True", "Recovered"), 3, each = length(x_plot))),
  region = factor(rep(c("Region A", "Region B", "Region C"), each = 2*length(x_plot)))
)

p_regional <- ggplot(df_regional, aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, alpha = 0.3) +
  facet_wrap(~ region, scales = "free_y") +
  scale_color_manual(values = c("True" = "black", "Recovered" = "red")) +
  scale_linetype_manual(values = c("True" = "dashed", "Recovered" = "solid")) +
  labs(title = "Regional Deviations Recovery",
       subtitle = "Region-specific patterns after removing global component",
       x = "x", y = "Deviation from global") +
  theme_bw() +
  theme(legend.title = element_blank())

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
  geom_point(aes(y = y, color = region_name), alpha = 0.5, size = 1) +
  geom_line(aes(y = y_fitted, color = region_name), linewidth = 1) +
  geom_line(aes(y = y_true), linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ region_name, scales = "free_y") +
  labs(title = "Overall Fit Quality",
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
  recovered_dev <- regional_deviations[[i]]
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