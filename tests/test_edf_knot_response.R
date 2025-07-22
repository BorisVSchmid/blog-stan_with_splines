# Test if EDF (Effective Degrees of Freedom) responds to knot number changes

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)
source("code/smoothing_diagnostics.R")

# Generate test data
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x
y <- y_true + rnorm(n, 0, 0.15)

cat("Testing EDF Response to Different Knot Numbers\n")
cat("==============================================\n")
cat("Data: n =", n, "points\n\n")

# Compile model once
model <- cmdstan_model("code/bsplines.stan")

# Test different knot numbers
knot_numbers <- c(3, 5, 7, 10, 15, 20)
results <- data.frame(
  num_knots = integer(),
  num_basis = integer(),
  edf = numeric(),
  sigma = numeric(),
  autocor = numeric()
)

for (k in knot_numbers) {
  cat("\nTesting with", k, "knots...\n")
  
  # Prepare data
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = k,
    spline_degree = 3,
    smoothing_strength = 0,
    prior_scale = 2 * sd(y)
  )
  
  # Fit model
  fit <- model$sample(
    data = stan_data,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 500,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Get diagnostics
  diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "bspline")
  
  # Store results
  results <- rbind(results, data.frame(
    num_knots = k,
    num_basis = k + 3 - 1,
    edf = diagnosis$edf,
    sigma = diagnosis$sigma_estimate,
    autocor = diagnosis$residual_autocor
  ))
  
  cat("  Basis functions:", k + 3 - 1, "\n")
  cat("  EDF:", round(diagnosis$edf, 1), "\n")
  cat("  Sigma:", round(diagnosis$sigma_estimate, 3), "\n")
}

# Display summary
cat("\n\nSummary of Results\n")
cat("==================\n")
print(results)

# Calculate EDF calculation details for one fit
cat("\n\nDetailed EDF Calculation Check\n")
cat("==============================\n")

# Refit with moderate knots
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 7,
  spline_degree = 3,
  smoothing_strength = 0,
  prior_scale = 2 * sd(y)
)

fit <- model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

# Manual EDF calculation check
draws <- fit$draws(format = "matrix")
y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])

# Method 1: From our diagnostic
influence_matrix_trace <- sum((y_hat - mean(y))^2) / sum((y - mean(y))^2) * n
cat("\nMethod 1 (scaled variance):", round(influence_matrix_trace, 1), "\n")

# Method 2: Alternative calculation
residuals <- y - y_hat
rss <- sum(residuals^2)
tss <- sum((y - mean(y))^2)
r_squared <- 1 - rss/tss
alt_edf <- r_squared * (n - 1) + 1
cat("Method 2 (R-squared based):", round(alt_edf, 1), "\n")

# Check if y_hat values are reasonable
cat("\nFitted values check:\n")
cat("  Range of y:", round(range(y), 2), "\n")
cat("  Range of y_hat:", round(range(y_hat), 2), "\n")
cat("  Variance of y:", round(stats::var(y), 3), "\n")
cat("  Variance of y_hat:", round(stats::var(y_hat), 3), "\n")

# Plot EDF vs knots
library(ggplot2)
p <- ggplot(results, aes(x = num_knots, y = edf)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = n, linetype = "dashed", color = "red") +
  geom_text(aes(label = paste0("n=", n)), x = min(knot_numbers), y = n + 1, color = "red") +
  labs(
    title = "EDF Response to Number of Knots",
    x = "Number of Knots",
    y = "Effective Degrees of Freedom",
    subtitle = "Testing if EDF is stuck near 35 or responds to knot changes"
  ) +
  theme_bw()

ggsave("claude-tmp/edf_vs_knots.png", p, width = 8, height = 6)