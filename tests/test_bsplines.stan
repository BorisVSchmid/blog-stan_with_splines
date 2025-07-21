// Minimal B-spline implementation for testing
// Based on: https://mc-stan.org/learn-stan/case-studies/splines_in_stan.html
// Original implementation by: Milad Kharratzadeh (2017)
// This is a derivative work - please check the original source for licensing terms

functions {
  // B-spline basis function - Direct implementation from Stan documentation
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    
    if (order==1)
      // B-splines of order 1 are piece-wise constant
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of 
      // two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=1> n_data;             // Number of data points
  array[n_data] real x;            // x values
  vector[n_data] y;                // y values to fit
  
  // Spline parameters
  int<lower=1> num_knots;          // Number of interior knots
  int<lower=1> spline_degree;      // Degree of spline (3 for cubic)
}

transformed data {
  // Calculate dimensions
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  
  // Set up knots
  vector[num_knots] knots;
  array[2*spline_degree + num_knots] real ext_knots_arr;
  
  // B-spline basis matrix
  matrix[num_basis, n_data] B;  // Basis functions as rows, data points as columns
  
  // Place knots at quantiles of the x data
  {
    array[n_data] real x_sorted = sort_asc(x);
    real x_min = x_sorted[1];
    real x_max = x_sorted[n_data];
    
    // Place interior knots at equally spaced quantiles
    for (k in 1:num_knots) {
      real p = k * 1.0 / (num_knots + 1);
      int idx = to_int(floor(p * n_data + 0.5));
      if (idx < 1) idx = 1;
      if (idx > n_data) idx = n_data;
      knots[k] = x_sorted[idx];
    }
  }
  
  // Extended knots: repeat boundary knots spline_degree times
  for (i in 1:spline_degree) {
    ext_knots_arr[i] = knots[1];
  }
  for (i in 1:num_knots) {
    ext_knots_arr[spline_degree + i] = knots[i];
  }
  for (i in 1:spline_degree) {
    ext_knots_arr[spline_degree + num_knots + i] = knots[num_knots];
  }
  
  // Build B-spline basis matrix
  for (i in 1:num_basis) {
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order));
  }
  
  // Print diagnostics
  print("B-spline setup:");
  print("  Number of data points: ", n_data);
  print("  Number of interior knots: ", num_knots);
  print("  Spline degree: ", spline_degree);
  print("  Number of basis functions: ", num_basis);
  print("  Knot locations: ", knots);
  print("  Data range: [", min(x), ", ", max(x), "]");
}

parameters {
  // Spline coefficients
  row_vector[num_basis] alpha;
  
  // Noise
  real<lower=0> sigma;
}

transformed parameters {
  // Compute fitted values
  vector[n_data] y_hat = to_vector(alpha * B);
}

model {
  // Priors
  alpha ~ normal(0, 10);
  sigma ~ normal(0, 1);
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // For plotting: evaluate spline on a fine grid
  int n_plot = 100;
  array[100] real x_plot;  // Must use literal for array size
  matrix[num_basis, 100] B_plot;
  vector[100] y_plot;
  
  // Create plotting grid
  real x_min = min(x) - 0.1 * (max(x) - min(x));
  real x_max = max(x) + 0.1 * (max(x) - min(x));
  
  for (i in 1:n_plot) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / (n_plot - 1.0);
  }
  
  // Build basis for plotting
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order));
  }
  
  // Evaluate spline
  y_plot = to_vector(alpha * B_plot);
  
  // Also output individual basis functions for visualization
  matrix[100, num_basis] basis_functions;  // Must use literal for array size
  for (i in 1:num_basis) {
    basis_functions[:, i] = to_vector(B_plot[i, :]);
  }
}
