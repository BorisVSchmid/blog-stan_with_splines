
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

