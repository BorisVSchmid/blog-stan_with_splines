// Minimal natural cubic spline implementation for testing
// Using stan-splines library by Sergey Koposov (University of Edinburgh)
// https://github.com/segasai/stan-splines
// Please cite: https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract
// DOI: https://doi.org/10.5281/zenodo.7193910
// This is a derivative work - please check the original source for licensing terms

functions {
  // Include the spline functions from the library
  #include "cspline_library/spline.stan"
}

data {
  int<lower=1> n_data;             // Number of data points
  vector[n_data] x;                // x values (as vector for spline functions)
  vector[n_data] y;                // y values to fit
  
  // Spline parameters
  int<lower=3> num_knots;          // Number of knots for natural cubic spline
}

transformed data {
  // Input validation
  if (num_knots < 2) reject("num_knots must be at least 2");
  if (n_data < num_knots) reject("insufficient data points for given number of knots");
  
  // Knot locations
  vector[num_knots] knot_locations;
  
  // Find position of each data point relative to knots
  array[n_data] int x_knot_positions;
  
  // Place knots at quantiles of the x data
  {
    vector[n_data] x_sorted = sort_asc(x);
    real x_min = x_sorted[1];
    real x_max = x_sorted[n_data];
    real buffer = 0.001 * (x_max - x_min);  // Small buffer to ensure all points are inside
    
    // First and last knots must encompass all data
    knot_locations[1] = x_min - buffer;
    knot_locations[num_knots] = x_max + buffer;
    
    // Interior knots at quantiles
    if (num_knots > 2) {
      for (k in 2:(num_knots-1)) {
        real position = (n_data - 1.0) * (k - 1.0) / (num_knots - 1.0) + 1.0;
        int idx_low = to_int(floor(position));
        int idx_high = to_int(ceil(position));
        if (idx_high > n_data) idx_high = n_data;
        
        if (idx_low == idx_high) {
          knot_locations[k] = x_sorted[idx_low];
        } else {
          real weight = position - idx_low;
          knot_locations[k] = x_sorted[idx_low] * (1 - weight) + x_sorted[idx_high] * weight;
        }
      }
    }
  }
  
  // Find which knot interval each data point belongs to
  x_knot_positions = spline_findpos(knot_locations, x);
  
  // Print diagnostics
  print("Natural cubic spline setup:");
  print("  Number of data points: ", n_data);
  print("  Number of knots: ", num_knots);
  print("  Knot locations: ", knot_locations);
  print("  Data range: [", min(x), ", ", max(x), "]");
}

parameters {
  // Values at knots
  vector[num_knots] y_at_knots;
  
  // Noise
  real<lower=0> sigma;
}

transformed parameters {
  // Get spline coefficients
  vector[num_knots] spline_coeffs = spline_getcoeffs(knot_locations, y_at_knots);
  
  // Evaluate spline at data points
  vector[n_data] y_hat = spline_eval(knot_locations, y_at_knots, spline_coeffs, x, x_knot_positions);
}

model {
  // Priors
  y_at_knots ~ normal(0, 10);
  sigma ~ normal(0, 1);
  
  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // For plotting: evaluate spline on a fine grid
  int n_plot = 1000;
  vector[1000] x_plot;  // Must use literal for array size
  array[1000] int x_plot_positions;
  vector[1000] y_plot;
  
  // Create plotting grid
  real x_min = min(x) - 0.1 * (max(x) - min(x));
  real x_max = max(x) + 0.1 * (max(x) - min(x));
  
  // Stay within knot range to avoid extrapolation issues
  if (x_min < knot_locations[1]) {
    x_min = knot_locations[1];
  }
  if (x_max > knot_locations[num_knots]) {
    x_max = knot_locations[num_knots];
  }
  
  for (i in 1:1000) {
    x_plot[i] = x_min + (x_max - x_min) * (i - 1.0) / 999.0;
  }
  
  // Find positions and evaluate spline
  x_plot_positions = spline_findpos(knot_locations, x_plot);
  y_plot = spline_eval(knot_locations, y_at_knots, spline_coeffs, x_plot, x_plot_positions);
  
  // Output the spline coefficients for diagnostics
  vector[num_knots] coefficients = spline_coeffs;
}
