// Regional B-splines: Multiple splines for different regions with shared priors
// Simplified version with better priors

functions {
  // B-spline basis function
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
  int<lower=1> n_total;                    // Total number of observations
  int<lower=1> n_regions;                  // Number of regions
  array[n_total] real x;                   // x values
  vector[n_total] y;                       // y values to fit
  array[n_total] int<lower=1,upper=n_regions> region;  // Region indicator
  
  // Spline parameters
  int<lower=1> num_knots;                  // Number of interior knots
  int<lower=1> spline_degree;              // Degree of spline (3 for cubic)
  real<lower=0> smoothing_strength;        // Smoothing strength (0=none, 1-2=mild, 5-10=strong)
  real<lower=0> prior_scale;               // Scale for coefficient priors
}

transformed data {
  // Calculate dimensions
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  
  // Transform smoothing_strength to tau_smooth
  real tau_smooth;
  if (smoothing_strength == 0) {
    tau_smooth = 0;  // Special case: independent coefficients
  } else {
    tau_smooth = 1 / sqrt(smoothing_strength);  // Convert to SD scale
  }
  
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
    B[i, :] = to_row_vector(build_b_spline(x, ext_knots_arr, i, spline_order, spline_degree));
  }
  
  print("Regional B-spline setup:");
  print("  Number of regions: ", n_regions);
  print("  Number of basis functions: ", num_basis);
  print("  Knot locations: ", knots);
}

parameters {
  // Raw coefficients for non-centered parameterization
  array[n_regions] row_vector[num_basis] alpha_raw;
  
  // Hierarchical priors for spline coefficients
  row_vector[num_basis] mu_alpha_raw;          // Global mean (raw)
  vector<lower=0>[num_basis] tau_alpha;        // Global SD
  
  // Shared noise parameter
  real<lower=0> sigma;
  
  // Regional baseline shifts
  vector[n_regions] beta_region;
  real mu_beta;
  real<lower=0> sigma_beta;
}

transformed parameters {
  // Apply random walk smoothing to get actual coefficients
  array[n_regions] row_vector[num_basis] alpha;
  row_vector[num_basis] mu_alpha;
  
  // Apply smoothing based on tau_smooth value
  if (tau_smooth == 0) {
    // No smoothing - independent coefficients
    mu_alpha = mu_alpha_raw * prior_scale;
    for (r in 1:n_regions) {
      alpha[r] = alpha_raw[r] .* tau_alpha';
    }
  } else {
    // Random walk smoothing
    mu_alpha[1] = mu_alpha_raw[1] * prior_scale;
    for (j in 2:num_basis) {
      mu_alpha[j] = mu_alpha[j-1] + mu_alpha_raw[j] * tau_smooth;
    }
    
    for (r in 1:n_regions) {
      alpha[r][1] = alpha_raw[r][1] * tau_alpha[1];
      for (j in 2:num_basis) {
        alpha[r][j] = alpha[r][j-1] + alpha_raw[r][j] * tau_smooth * tau_alpha[j];
      }
    }
  }
  
  // Compute fitted values for each observation
  vector[n_total] y_hat;
  
  for (i in 1:n_total) {
    y_hat[i] = alpha[region[i]] * B[:, i] + beta_region[region[i]];
  }
}

model {
  // Priors for raw coefficients (non-centered parameterization)
  mu_alpha_raw ~ std_normal();
  tau_alpha ~ normal(0, prior_scale);
  
  for (r in 1:n_regions) {
    alpha_raw[r] ~ std_normal();
  }
  
  // Priors for regional effects
  mu_beta ~ normal(mean(y), fmax(sd(y), 0.1 * fmax(abs(mean(y)), 1.0)));
  sigma_beta ~ exponential(2);
  beta_region ~ normal(mu_beta, sigma_beta);
  
  // Shared noise prior
  sigma ~ exponential(2);
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // Predictions for each region on a fine grid
  array[1000] real x_plot;
  matrix[num_basis, 1000] B_plot;
  array[n_regions] vector[1000] y_plot;
  
  // Variance components
  real var_between = variance(beta_region);
  real var_within = square(sigma);
  real icc = var_between / (var_between + var_within);
  
  // Create plotting grid
  real x_min = min(x);
  real x_max = max(x);
  
  for (i in 1:1000) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 999.0;
  }
  
  // Build basis for plotting
  for (i in 1:num_basis) {
    B_plot[i, :] = to_row_vector(build_b_spline(x_plot, ext_knots_arr, i, spline_order, spline_degree));
  }
  
  // Evaluate spline for each region
  for (r in 1:n_regions) {
    y_plot[r] = to_vector(alpha[r] * B_plot) + beta_region[r];
  }
}

