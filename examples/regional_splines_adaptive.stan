// Regional B-splines with adaptive shrinkage
// Improved version with data-driven shrinkage for regional deviations

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
  real<lower=0> smoothing_strength_global;  // Smoothing strength for global pattern
  real<lower=0> smoothing_strength_regional; // Smoothing strength for regional deviations
  real<lower=0> prior_scale;               // Scale for coefficient priors
}

transformed data {
  // Calculate dimensions
  int num_basis = num_knots + spline_degree - 1;
  int spline_order = spline_degree + 1;
  
  // Transform smoothing_strength to tau_smooth for global and regional
  real tau_smooth_global;
  real tau_smooth_regional;
  
  if (smoothing_strength_global == 0) {
    tau_smooth_global = 0;  // Special case: independent coefficients
  } else {
    tau_smooth_global = prior_scale / (smoothing_strength_global * pow(num_basis, sqrt(2)));
  }
  
  if (smoothing_strength_regional == 0) {
    tau_smooth_regional = 0;  // Special case: independent coefficients
  } else {
    tau_smooth_regional = prior_scale / (smoothing_strength_regional * pow(num_basis, sqrt(2)));
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
  
  print("Adaptive Regional B-spline setup:");
  print("  Number of regions: ", n_regions);
  print("  Number of basis functions: ", num_basis);
  print("  Knot locations: ", knots);
}

parameters {
  // Global spline coefficients (raw for smoothing)
  row_vector[num_basis] mu_alpha_raw;          // Global pattern (raw)
  
  // Regional deviations from global pattern (raw for non-centered)
  array[n_regions] row_vector[num_basis] alpha_deviation_raw;
  
  // Hierarchical SD for regional deviations (now with adaptive shrinkage)
  vector<lower=0>[num_basis] tau_alpha_base;   // Base SD of deviations
  real<lower=0> shrinkage_factor;              // Adaptive shrinkage parameter
  
  // Shared noise parameter
  real<lower=0> sigma;
  
  // Regional baseline shifts
  vector[n_regions] beta_region;
  real mu_beta;
  real<lower=0> sigma_beta;
}

transformed parameters {
  // Apply random walk smoothing to get actual coefficients
  row_vector[num_basis] mu_alpha;
  array[n_regions] row_vector[num_basis] alpha_deviation;
  array[n_regions] row_vector[num_basis] alpha;  // Total regional coefficients
  vector[num_basis] tau_alpha;                   // Adaptively shrunk tau values
  
  // Apply smoothing to global pattern
  if (tau_smooth_global == 0) {
    // No smoothing - independent coefficients
    mu_alpha = mu_alpha_raw * prior_scale;
  } else {
    // Random walk smoothing
    mu_alpha[1] = mu_alpha_raw[1] * prior_scale;
    for (j in 2:num_basis) {
      mu_alpha[j] = mu_alpha[j-1] + mu_alpha_raw[j] * tau_smooth_global;
    }
  }
  
  // Apply adaptive shrinkage to tau_alpha
  tau_alpha = tau_alpha_base * (prior_scale / shrinkage_factor);
  
  // Transform regional deviations
  for (r in 1:n_regions) {
    if (tau_smooth_regional == 0) {
      // No smoothing - independent deviations
      alpha_deviation[r] = alpha_deviation_raw[r] .* tau_alpha';
    } else {
      // Random walk smoothing for deviations
      alpha_deviation[r][1] = alpha_deviation_raw[r][1] * tau_alpha[1];
      for (j in 2:num_basis) {
        alpha_deviation[r][j] = alpha_deviation[r][j-1] + 
                                alpha_deviation_raw[r][j] * tau_smooth_regional * tau_alpha[j];
      }
    }
    
    // Total regional coefficients = global + deviation
    alpha[r] = mu_alpha + alpha_deviation[r];
  }
  
  // Compute fitted values for each observation
  vector[n_total] y_hat;
  
  for (i in 1:n_total) {
    y_hat[i] = alpha[region[i]] * B[:, i] + beta_region[region[i]];
  }
}

model {
  // Prior for global pattern coefficients
  mu_alpha_raw ~ std_normal();
  
  // Prior for base deviation scales
  tau_alpha_base ~ std_normal();  // Unit scale, will be scaled by shrinkage
  
  // Adaptive shrinkage prior - weakly informative, centered at 5 (current fixed value)
  shrinkage_factor ~ normal(5, 2);
  
  // Priors for regional deviations (non-centered)
  for (r in 1:n_regions) {
    alpha_deviation_raw[r] ~ std_normal();
  }
  
  // Priors for regional baseline effects
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
  
  // Variance explained by global vs regional patterns
  real var_global = 0;
  real var_regional = 0;
  
  // Diagnostic: effective degrees of freedom for regional deviations
  real edf_regional = sum(tau_alpha) / mean(tau_alpha);
  
  {
    // Compute variance of global pattern across plotting grid
    vector[1000] global_pattern;
    array[n_regions] vector[1000] regional_patterns;
    
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
    
    // Compute patterns
    global_pattern = to_vector(mu_alpha * B_plot) + mu_beta;
    for (r in 1:n_regions) {
      regional_patterns[r] = to_vector(alpha_deviation[r] * B_plot);
      y_plot[r] = to_vector(alpha[r] * B_plot) + beta_region[r];
    }
    
    var_global = variance(global_pattern);
    for (r in 1:n_regions) {
      var_regional += variance(regional_patterns[r]) / n_regions;
    }
  }
}