# Quick test of spline implementations
library(cmdstanr)
set.seed(123)

# Generate simple test data
n <- 20
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.1)

# Test B-spline
cat("Testing B-spline model...\n")
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 5,
  spline_degree = 3
)

model_b <- cmdstan_model("test_bsplines.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

cat("B-spline fit summary:\n")
fit_b$summary(c("sigma", "lp__"))

# Test C-spline
cat("\nTesting C-spline model...\n")
stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 5
)

model_c <- cmdstan_model("test_csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

cat("C-spline fit summary:\n")
fit_c$summary(c("sigma", "lp__"))

cat("\nBoth models ran successfully!\n")