# Minimal viable B-spline implementation in Stan

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
  num_knots = 7,      # Number of knots (including boundaries)
                      # - Use 4-6 knots for simple smooth curves
                      # - Use 7-10 knots for more complex patterns
                      # - More knots = more flexibility but risk of overfitting
  spline_degree = 3   # Degree of B-spline (3 = cubic, most common)
                      # - Degree 1: piecewise linear (connect-the-dots)
                      # - Degree 2: piecewise quadratic (smooth but can be wobbly)
                      # - Degree 3: piecewise cubic (smooth, most common choice)
                      # - Higher degrees rarely needed in practice
)

# Compile and fit model
model <- cmdstan_model("code/bsplines.stan")

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
plot(x, y, main = "B-spline Fit", xlab = "x", ylab = "y")
lines(x, y_hat, col = "blue", lwd = 2)
lines(x, y_true, col = "red", lty = 2)
legend("topright", c("Data", "B-spline fit", "True function"), 
       col = c("black", "blue", "red"), 
       pch = c(1, NA, NA), lty = c(NA, 1, 2))