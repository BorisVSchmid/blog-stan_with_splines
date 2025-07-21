// B-splines with random walk smoothing prior
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
  int<lower=0> n_data;
  array[n_data] real x;
  vector[n_data] y;
  int<lower=1> num_knots;
  int<lower=1> spline_degree;
  real<lower=0> tau_fixed;  // Fixed smoothing parameter
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
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order));
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
  
  // Random walk prior with fixed smoothing
  alpha[1] = alpha_raw[1];
  for (i in 2:num_basis) {
    alpha[i] = alpha[i-1] + alpha_raw[i] * tau_fixed;
  }
  
  y_hat = alpha_0 + to_vector(alpha' * B);
}

model {
  // Priors
  alpha_0 ~ normal(mean(y), sd(y));
  alpha_raw ~ std_normal();  // Standard normal for non-centered parameterization
  sigma ~ normal(0, 1);
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  array[100] real x_plot;
  matrix[num_basis, 100] B_plot;
  vector[100] y_plot;
  
  real x_min = min(x);
  real x_max = max(x);
  
  for (i in 1:100) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 99.0;
  }
  
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order));
  }
  
  y_plot = alpha_0 + to_vector(alpha' * B_plot);
}