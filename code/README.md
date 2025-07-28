# Stan Splines Core Implementations

This directory contains the core Stan models and R functions for B-splines and C-splines.

For detailed differences between our implementations and the original sources, see [IMPLEMENTATION_DIFFERENCES.md](IMPLEMENTATION_DIFFERENCES.md).

## Stan Models

### bsplines.stan
Self-contained B-spline implementation with:
- Cox-de Boor recursive algorithm for basis function construction
- Adaptive knot placement at quantiles
- Scale-invariant smoothing prior (operates knot-to-knot)
- Boundary handling via polynomial continuation
- Number of basis functions = `num_knots + spline_degree - 1`

Parameters:
- `num_knots`: Number of interior knots
- `spline_degree`: Polynomial degree (default 3 for cubic)
- `smoothing_strength`: Controls smoothness (0=none, 0.05-0.1=mild, 0.1-0.2=strong)
- `prior_scale`: Scale for coefficient priors

### csplines.stan
Natural cubic spline implementation that:
- Depends on stan-math spline library via `#include`
- Uses matrix algebra for efficient computation
- Enforces natural boundary conditions (linear extrapolation)
- Requires all evaluation points within knot range
- Global support (each basis affects entire curve)

Parameters:
- `num_knots`: Total number of knots including boundaries

## R Scripts

### bsplines.R
Demonstrates B-spline usage with:
- Adaptive knot selection (n/2 rule)
- Data-driven prior scaling
- MCMC diagnostics output
- Smoothing diagnostics
- Visualization with uncertainty intervals

Key features:
- Runs 4 chains in parallel
- Prints MCMC diagnostic summary
- Generates publication-quality plots

### csplines.R
Demonstrates C-spline usage with:
- Simple knot configuration
- Natural boundary conditions
- MCMC diagnostics output
- Smoothing diagnostics
- Visualization with uncertainty intervals

Key features:
- Runs 4 chains in parallel
- Prints MCMC diagnostic summary
- Generates publication-quality plots

### smoothing_diagnostics.R
Advanced diagnostic functions for spline models:
- Effective degrees of freedom (EDF) calculation
- Residual analysis and autocorrelation
- Cross-validation metrics
- Model complexity assessment
- Automated recommendations for parameter tuning

Functions:
- `diagnose_smoothing()`: Main diagnostic function
- `print_smoothing_diagnostics()`: Pretty-print results
- `fit_with_diagnostics()`: Wrapper for fitting with diagnostics


## Usage Examples

### Basic B-spline fitting
```r
source("code/bsplines.R")
# Generates test data, fits model, shows diagnostics and plots
```

### Basic C-spline fitting
```r
source("code/csplines.R")
# Generates test data, fits model, shows diagnostics and plots
```

### Custom B-spline analysis
```r
library(cmdstanr)
source("code/smoothing_diagnostics.R")

# Your data
x <- your_x_values
y <- your_y_values

# Compile model
model <- cmdstan_model("code/bsplines.stan")

# Set up data
stan_data <- list(
  n_data = length(x),
  x = x,
  y = y,
  num_knots = 10,
  spline_degree = 3,
  smoothing_strength = 0.1,
  prior_scale = 2 * sd(y)
)

# Fit with diagnostics
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000
)

# Check diagnostics
print(fit$diagnostic_summary())
diagnosis <- diagnose_smoothing(fit, x, y, stan_data, "bspline")
print_smoothing_diagnostics(diagnosis)
```

## Key Differences: B-splines vs C-splines

### B-splines
- Local support (each basis affects only nearby region)
- More flexible with many knots
- Smoothing prior helps control wiggliness
- Better for complex, multi-scale features
- Polynomial extrapolation at boundaries

### C-splines
- Global support (each basis affects entire curve)
- Naturally smooth (C² continuous)
- Fewer knots needed
- Better for smooth, global trends
- Linear extrapolation at boundaries

## Parameter Guidelines

### B-splines
- `num_knots`: Start with n/2 for flexibility, adjust based on diagnostics
- `smoothing_strength`: 0=none, 0.05-0.1=mild (0.1 default), 0.1-0.2=strong (scales with num_basis^sqrt(2))
- `spline_degree`: 3 (cubic) is standard, 2 (quadratic) for simpler curves

### C-splines
- `num_knots`: Start with n/4, rarely need more than 10-15
- No smoothing parameter (inherently smooth)

## Output Files

Running the example scripts creates:
- `output/code-bspline_minimal_example.png`: B-spline demonstration
- `output/code-cspline_minimal_example.png`: C-spline demonstration

## Dependencies

All scripts use:
- `cmdstanr` for Stan interface
- `groundhog` for reproducible package management
- `conflicted` for namespace management
- Standard tidyverse packages for data manipulation and plotting