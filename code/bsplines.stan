// B-splines with scale-invariant random walk smoothing prior
// 
// Key feature: The smoothing strength parameter works consistently across different data scales.
// A smoothing_strength of 4.0 provides the same level of regularization whether your data
// ranges from -1 to 1 or -1000 to 1000. This is achieved by scaling the random walk step
// size by the data's standard deviation.
functions {
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order, int degree) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    
    if (order==1)
      for (i in 1:size(t)) {
        // For the rightmost data point, assign it to the last non-zero interval
        if (t[i] >= ext_knots[size(ext_knots) - degree] && 
            ind == size(ext_knots) - degree - 1) {
          b_spline[i] = 1;
        } else {
          b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
        }
      }
    else {
      if (abs(ext_knots[ind] - ext_knots[ind+order-1]) > 1e-10)
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (abs(ext_knots[ind+1] - ext_knots[ind+order]) > 1e-10)
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1, degree) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1, degree);
    }
    return b_spline;
  }
}

data {
  int<lower=0> n_data;
  array[n_data] real x;
  vector[n_data] y;
  int<lower=1> num_knots;
  int<lower=1> spline_degree;
  real<lower=0> smoothing_strength;  // Smoothing strength (0=none, 1-2=mild, 5-10=strong)
  real<lower=0> prior_scale;  // Scale for coefficient priors
}

transformed data {
  // Input validation
  if (num_knots < 2) reject("num_knots must be at least 2");
  if (spline_degree < 1) reject("spline_degree must be at least 1");
  if (n_data < num_knots) reject("insufficient data points for given number of knots");
  
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  vector[num_knots] knots;
  array[2*spline_degree + num_knots] real ext_knots_arr;
  matrix[num_basis, n_data] B;
  
  // Transform smoothing_strength to tau_smooth (random walk SD)
  // smoothing_strength: 0=none, 1-2=mild, 5-10=strong (scales with num_basis)
  // tau_smooth: 0=special case for independent, small=smooth, large=flexible
  //
  // CRITICAL: Scale-invariant smoothing implementation
  // The random walk prior is: alpha[i] = alpha[i-1] + alpha_raw[i] * tau_smooth
  // To make smoothing strength consistent across different data scales AND different
  // numbers of knots, we scale tau_smooth by both prior_scale (data variance) and
  // number of basis functions.
  // Without this scaling, the same smoothing_strength would:
  // 1. Underfit large-scale data while overfitting small-scale data
  // 2. Have different effects with different numbers of knots
  real tau_smooth;
  if (smoothing_strength == 0) {
    tau_smooth = 0;  // Special case: triggers independent coefficients
  } else {
    tau_smooth = prior_scale / sqrt(smoothing_strength * num_basis);  // Scale by data variance and number of basis functions
  }
  
  {
    array[n_data] real x_sorted = sort_asc(x);
    real x_min = x_sorted[1];
    real x_max = x_sorted[n_data];
    
    knots[1] = x_min;
    knots[num_knots] = x_max;
    
    if (num_knots > 2) {
      for (k in 2:(num_knots-1)) {
        knots[k] = x_min + (x_max - x_min) * (k - 1.0) / (num_knots - 1.0);
      }
    }
  }
  
  for (i in 1:spline_degree) {
    ext_knots_arr[i] = knots[1];
  }
  for (i in 1:num_knots) {
    ext_knots_arr[spline_degree + i] = knots[i];
  }
  for (i in 1:spline_degree) {
    ext_knots_arr[spline_degree + num_knots + i] = knots[num_knots];
  }
  
  for (i in 1:num_basis) {
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order, spline_degree));
  }
}

parameters {
  real alpha_0;  // Intercept
  vector[num_basis] alpha_raw;  // Raw coefficients for non-centered parameterization
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_basis] alpha;  // Actual spline coefficients
  vector[n_data] y_hat;
  
  // Apply smoothing based on tau_smooth value
  if (tau_smooth == 0) {
    // No smoothing - independent coefficients scaled by prior_scale
    alpha = alpha_raw * prior_scale;
  } else {
    // Random walk smoothing with scale-invariant steps
    // First coefficient starts at the data scale
    alpha[1] = alpha_raw[1] * prior_scale;
    // Subsequent coefficients follow a random walk with steps proportional to data scale
    // This ensures consistent smoothing behavior regardless of whether y ranges from
    // -1 to 1 or -1000 to 1000
    for (i in 2:num_basis) {
      alpha[i] = alpha[i-1] + alpha_raw[i] * tau_smooth;
    }
  }
  
  y_hat = alpha_0 + to_vector(alpha' * B);
}

model {
  // Priors
  alpha_0 ~ normal(mean(y), fmax(sd(y), 0.1 * fmax(abs(mean(y)), 1.0)));  // Scale adapts to data magnitude, min 0.1
  alpha_raw ~ std_normal();  // Standard normal for non-centered parameterization
  sigma ~ exponential(2);  // Better prior for positive-constrained parameter
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  array[1000] real x_plot;
  matrix[num_basis, 1000] B_plot;
  vector[1000] y_plot;
  
  real x_min = min(x);
  real x_max = max(x);
  
  for (i in 1:1000) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 999.0;
  }
  
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order, spline_degree));
  }
  
  y_plot = alpha_0 + to_vector(alpha' * B_plot);
}
