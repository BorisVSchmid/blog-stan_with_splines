// Regional B-splines: Multiple splines for different regions with shared priors
// Simplified version with better priors

functions {
  // B-spline basis function
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
  int<lower=1> n_total;                    // Total number of observations
  int<lower=1> n_regions;                  // Number of regions
  array[n_total] real x;                   // x values
  vector[n_total] y;                       // y values to fit
  array[n_total] int<lower=1,upper=n_regions> region;  // Region indicator
  
  // Spline parameters
  int<lower=1> num_knots;                  // Number of interior knots
  int<lower=1> spline_degree;              // Degree of spline (3 for cubic)
}

transformed data {
  // Calculate dimensions
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  
  // Set up knots
  vector[num_knots] knots;
  array[2*spline_degree + num_knots] real ext_knots_arr;
  
  // B-spline basis matrix
  matrix[num_basis, n_total] B;
  
  // Place knots based on overall data range
  {
    array[n_total] real x_sorted = sort_asc(x);
    real x_min = x_sorted[1];
    real x_max = x_sorted[n_total];
    
    // Include boundaries
    knots[1] = x_min;
    knots[num_knots] = x_max;
    
    // Interior knots
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
  
  // Build B-spline basis matrix
  for (i in 1:num_basis) {
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order));
  }
  
  print("Regional B-spline setup:");
  print("  Number of regions: ", n_regions);
  print("  Number of basis functions: ", num_basis);
  print("  Knot locations: ", knots);
}

parameters {
  // Regional spline coefficients
  array[n_regions] row_vector[num_basis] alpha;
  
  // Hierarchical priors for spline coefficients
  row_vector[num_basis] mu_alpha;              // Global mean
  vector<lower=0>[num_basis] tau_alpha;        // Global SD
  
  // Shared noise parameter
  real<lower=0> sigma;
  
  // Regional baseline shifts
  vector[n_regions] beta_region;
  real mu_beta;
  real<lower=0> sigma_beta;
}

transformed parameters {
  // Compute fitted values for each observation
  vector[n_total] y_hat;
  
  for (i in 1:n_total) {
    y_hat[i] = alpha[region[i]] * B[:, i] + beta_region[region[i]];
  }
}

model {
  // Hierarchical priors for spline coefficients
  mu_alpha ~ normal(0, 2);
  tau_alpha ~ normal(0, 1);
  
  for (r in 1:n_regions) {
    alpha[r] ~ normal(mu_alpha, tau_alpha);
  }
  
  // Priors for regional effects
  mu_beta ~ normal(0, 1);
  sigma_beta ~ normal(0, 0.5);
  beta_region ~ normal(mu_beta, sigma_beta);
  
  // Shared noise prior
  sigma ~ normal(0, 0.5);
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // Predictions for each region on a fine grid
  array[100] real x_plot;
  matrix[num_basis, 100] B_plot;
  array[n_regions] vector[100] y_plot;
  
  // Variance components
  real var_between = variance(beta_region);
  real var_within = square(sigma);
  real icc = var_between / (var_between + var_within);
  
  // Create plotting grid
  real x_min = min(x);
  real x_max = max(x);
  
  for (i in 1:100) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 99.0;
  }
  
  // Build basis for plotting
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order));
  }
  
  // Evaluate spline for each region
  for (r in 1:n_regions) {
    y_plot[r] = to_vector(alpha[r] * B_plot) + beta_region[r];
  }
}