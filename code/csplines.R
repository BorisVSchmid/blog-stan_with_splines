# Minimal viable C-spline (natural cubic spline) implementation in Stan

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
# Use groundhog for dependencies
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

# Load cmdstanr after its dependencies
library(cmdstanr)

# Generate simple test data
n <- 30
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.1)

# Prepare data for Stan
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 7
)

# Compile and fit model
model <- cmdstan_model("code/csplines.stan")

fit <- model$sample(
  data = stan_data,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Extract and plot results
draws <- fit$draws(format = "matrix")
y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])

# Simple plot
plot(x, y, main = "C-spline (Natural Cubic Spline) Fit", xlab = "x", ylab = "y")
lines(x, y_hat, col = "darkgreen", lwd = 2)
lines(x, y_true, col = "red", lty = 2)
legend("topright", c("Data", "C-spline fit", "True function"), 
       col = c("black", "darkgreen", "red"), 
       pch = c(1, NA, NA), lty = c(NA, 1, 2))