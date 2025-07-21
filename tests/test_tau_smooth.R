# Test different tau_smooth values

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

# Generate test data
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x 
y <- y_true + rnorm(n, 0, 0.15)

# Test different tau values
tau_values <- c(0, 0.1, 0.5, 1.0, 2.0)
model <- cmdstan_model("code/bsplines.stan")

# Store results
results <- list()

for (tau in tau_values) {
  cat("\nTesting tau_smooth =", tau, "\n")
  
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 7,
    spline_degree = 3,
    tau_smooth = tau,
    prior_scale = 2 * sd(y)
  )
  
  fit <- model$sample(
    data = stan_data,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 500,
    refresh = 0
  )
  
  draws <- fit$draws(format = "matrix")
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  
  # Check variation in fitted values
  cat("  Range of y_hat:", round(range(y_hat), 2), "\n")
  cat("  SD of y_hat:", round(sd(y_hat), 3), "\n")
  
  # Store for plotting
  results[[as.character(tau)]] <- data.frame(
    x = x,
    y = y,
    y_hat = y_hat,
    tau = tau
  )
}

# Create comparison plot
library(ggplot2)
library(patchwork)

# Combine all results
all_results <- do.call(rbind, results)
all_results$tau_label <- paste0("tau = ", all_results$tau)

p <- ggplot(all_results, aes(x = x)) +
  geom_point(aes(y = y), alpha = 0.3) +
  geom_line(aes(y = y_hat, color = as.factor(tau)), size = 1) +
  facet_wrap(~ tau_label, ncol = 3) +
  labs(
    title = "Effect of tau_smooth on B-spline fits",
    subtitle = "tau = 0.1 appears to over-smooth dramatically",
    y = "y",
    color = "tau_smooth"
  ) +
  theme_bw()

print(p)
ggsave("claude-tmp/tau_smooth_comparison.png", p, width = 10, height = 6)

# Check the actual coefficient values for tau = 0.1
cat("\n\nDetailed check for tau = 0.1:\n")
stan_data_check <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 7,
  spline_degree = 3,
  tau_smooth = 0.1,
  prior_scale = 2 * sd(y)
)

fit_check <- model$sample(
  data = stan_data_check,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

draws_check <- fit_check$draws(format = "matrix")
alpha_cols <- grep("alpha\\[", colnames(draws_check), value = TRUE)
if (length(alpha_cols) > 0) {
  alpha_values <- colMeans(draws_check[, alpha_cols])
  cat("Alpha coefficients:", round(alpha_values, 3), "\n")
  cat("Range of alphas:", round(range(alpha_values), 3), "\n")
}

# The issue: with tau=0.1, the random walk is too constrained
cat("\nThe problem: tau_smooth scales the innovation standard deviation\n")
cat("With tau=0.1, each coefficient can only differ by ~0.1 from the previous\n")
cat("This is far too restrictive for typical spline fitting\n")