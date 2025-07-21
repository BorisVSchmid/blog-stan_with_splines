# Extended plotting script with better knot placement and extended C-spline range
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)
dir.create("output", showWarnings = FALSE)

# Generate test data
generate_data <- function(n = 30) {
  x <- seq(0, 10, length.out = n)
  y_true <- sin(x)
  y <- y_true + rnorm(n, 0, 0.1)
  list(x = x, y = y, y_true = y_true)
}

# Create modified B-spline model with boundary knots at data extremes
create_bspline_boundary_knots <- function() {
  stan_code <- '
// B-spline with knots including boundaries
functions {
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    
    if (order==1)
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=1> n_data;
  array[n_data] real x;
  vector[n_data] y;
  int<lower=1> num_knots;
  int<lower=1> spline_degree;
}

transformed data {
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  vector[num_knots] knots;
  array[2*spline_degree + num_knots] real ext_knots_arr;
  matrix[num_basis, n_data] B;
  
  {
    array[n_data] real x_sorted = sort_asc(x);
    real x_min = x_sorted[1];
    real x_max = x_sorted[n_data];
    
    // Include boundaries in knots
    knots[1] = x_min;
    knots[num_knots] = x_max;
    
    // Interior knots at equal spacing
    if (num_knots > 2) {
      for (k in 2:(num_knots-1)) {
        knots[k] = x_min + (x_max - x_min) * (k - 1.0) / (num_knots - 1.0);
      }
    }
  }
  
  // Extended knots
  for (i in 1:spline_degree) {
    ext_knots_arr[i] = knots[1];
  }
  for (i in 1:num_knots) {
    ext_knots_arr[spline_degree + i] = knots[i];
  }
  for (i in 1:spline_degree) {
    ext_knots_arr[spline_degree + num_knots + i] = knots[num_knots];
  }
  
  // Build basis
  for (i in 1:num_basis) {
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order));
  }
  
  print("B-spline with boundary knots:");
  print("  Knot locations: ", knots);
}

parameters {
  row_vector[num_basis] alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n_data] y_hat = to_vector(alpha * B);
}

model {
  alpha ~ normal(0, 10);
  sigma ~ normal(0, 1);
  y ~ normal(y_hat, sigma);
}

generated quantities {
  array[100] real x_plot;
  matrix[num_basis, 100] B_plot;
  vector[100] y_plot;
  
  real x_min = min(x) - 0.1 * (max(x) - min(x));
  real x_max = max(x) + 0.1 * (max(x) - min(x));
  
  for (i in 1:100) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 99.0;
  }
  
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order));
  }
  
  y_plot = to_vector(alpha * B_plot);
  
  matrix[100, num_basis] basis_functions;
  for (i in 1:num_basis) {
    basis_functions[:, i] = to_vector(B_plot[i, :]);
  }
}
'
  write(stan_code, "test_bsplines_boundary.stan")
}

# Fit and plot
data <- generate_data(30)

# Create modified B-spline model
create_bspline_boundary_knots()

# B-spline with boundary knots
cat("Fitting B-spline with boundary knots...\n")
stan_data_b <- list(
  n_data = length(data$x),
  x = data$x,
  y = data$y,
  num_knots = 7,
  spline_degree = 3
)

model_b <- cmdstan_model("test_bsplines_boundary.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

# C-spline with extended range
cat("Fitting C-spline with extended evaluation range...\n")
stan_data_c <- list(
  n_data = length(data$x),
  x = data$x,
  y = data$y,
  num_knots = 7
)

model_c <- cmdstan_model("test_csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0
)

# Extract and plot
cat("Creating plots...\n")
draws_b <- fit_b$draws(format = "matrix")
draws_c <- fit_c$draws(format = "matrix")

# B-spline plot
x_plot_b <- colMeans(draws_b[, grep("x_plot\\[", colnames(draws_b), value = TRUE)])
y_plot_b <- colMeans(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)])
y_lower_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.025)
y_upper_b <- apply(draws_b[, grep("y_plot\\[", colnames(draws_b), value = TRUE)], 2, quantile, 0.975)

# For C-spline, we need to evaluate on extended range (-2 to 12)
# Extract the fitted spline parameters
y_at_knots <- colMeans(draws_c[, grep("y_at_knots\\[", colnames(draws_c), value = TRUE)])
knot_locs <- stan_data_c$n_data  # We know knot locations from the model

# Create extended x range for plotting
x_extended <- seq(-2, 12, length.out = 200)

# For visualization, we'll show the natural behavior at boundaries
# C-splines are linear outside the knot range

plot_b <- ggplot() +
  geom_ribbon(aes(x = x_plot_b, ymin = y_lower_b, ymax = y_upper_b), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = x_plot_b, y = y_plot_b), color = "blue", linewidth = 1) +
  geom_point(data = data.frame(x = data$x, y = data$y), aes(x, y), alpha = 0.6, size = 2) +
  geom_line(data = data.frame(x = data$x, y = data$y_true), aes(x, y), 
            color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "B-spline fit (boundary knots at data extremes)", x = "x", y = "y") +
  theme_bw() +
  coord_cartesian(xlim = c(-0.5, 10.5))

# C-spline plot with extended range
x_plot_c <- colMeans(draws_c[, grep("x_plot\\[", colnames(draws_c), value = TRUE)])
y_plot_c <- colMeans(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)])
y_lower_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.025)
y_upper_c <- apply(draws_c[, grep("y_plot\\[", colnames(draws_c), value = TRUE)], 2, quantile, 0.975)

# Create true function for extended range
x_true_extended <- seq(-2, 12, length.out = 200)
y_true_extended <- sin(x_true_extended)

plot_c <- ggplot() +
  geom_ribbon(aes(x = x_plot_c, ymin = y_lower_c, ymax = y_upper_c), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = x_plot_c, y = y_plot_c), color = "blue", linewidth = 1) +
  geom_point(data = data.frame(x = data$x, y = data$y), aes(x, y), alpha = 0.6, size = 2) +
  geom_line(data = data.frame(x = x_true_extended, y = y_true_extended), aes(x, y), 
            color = "red", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(0, 10), linetype = "dotted", alpha = 0.5) +
  annotate("text", x = 5, y = -1.5, label = "Data range", size = 3) +
  annotate("text", x = -1, y = -1.5, label = "Linear\nextrapolation", size = 3) +
  annotate("text", x = 11, y = -1.5, label = "Linear\nextrapolation", size = 3) +
  labs(title = "C-spline (Natural Cubic Spline) - Extended view", 
       subtitle = "Note: C-splines extrapolate linearly beyond knot boundaries",
       x = "x", y = "y") +
  theme_bw() +
  coord_cartesian(xlim = c(-2, 12), ylim = c(-2, 2))

# Save individual plots
ggsave("output/bspline_boundary_knots.png", plot_b, width = 8, height = 6, dpi = 300)
ggsave("output/cspline_extended_range.png", plot_c, width = 10, height = 6, dpi = 300)

# Combined comparison
combined <- plot_b / plot_c
combined_with_title <- combined + 
  plot_annotation(
    title = "B-splines vs C-splines: Boundary Behavior",
    subtitle = "B-spline now includes boundary knots; C-spline shows linear extrapolation"
  )

ggsave("output/spline_comparison_extended.png", combined_with_title, width = 10, height = 12, dpi = 300)

cat("\nPlots saved:\n")
cat("- output/bspline_boundary_knots.png\n")
cat("- output/cspline_extended_range.png\n")
cat("- output/spline_comparison_extended.png\n")