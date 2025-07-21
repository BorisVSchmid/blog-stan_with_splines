# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Structure

- `code/` - Core Stan models and R functions
- `tests/` - Test scripts for verifying implementations
- `examples/` - Example analyses and demonstrations
- `output/` - Generated plots and results
- `claude-tmp/` - Temporary/experimental scripts created by Claude (not part of main project)

## Commands

### Running Tests
```r
# Run comprehensive test suite (tests multiple functions, creates plots)
source("test_splines.R")

# Run quick test (simple sine function test)
source("quick_test.R")
```

### Building Stan Models
Models are compiled automatically when running tests. To manually compile:
```r
library(cmdstanr)
model_b <- cmdstan_model("test_bsplines.stan")
model_c <- cmdstan_model("test_csplines.stan")
```

### Testing Individual Functions
```r
# Generate test data
data <- generate_test_data(n = 30, func_type = "sine", noise_sd = 0.1)
# func_type options: "sine", "polynomial", "step", "exponential", "complex", "linear"

# Fit B-spline
fit_b <- fit_bspline(data, num_knots = 7, spline_degree = 3)

# Fit C-spline
fit_c <- fit_cspline(data, num_knots = 7)
```

## Architecture

### Spline Implementations

The project provides two independent spline implementations in Stan:

1. **B-splines** (`test_bsplines.stan`)
   - Self-contained implementation using Cox-de Boor recursive algorithm
   - Function `build_b_spline()` constructs basis functions recursively
   - Knots placed at quantiles with extended boundary knots
   - Number of basis functions = `num_knots + spline_degree - 1`

2. **C-splines** (`test_csplines.stan`)
   - Depends on `spline.stan` library via `#include`
   - Uses matrix algebra for natural cubic splines
   - Requires all evaluation points within knot range
   - Functions: `spline_getcoeffs()`, `spline_eval()`, `spline_findpos()`

### Key Architectural Differences

- **B-splines**: Local support means coefficients (`alpha`) directly multiply basis functions
- **C-splines**: Global support with values specified at knots (`y_at_knots`), then transformed to coefficients

### Testing Framework

`test_splines.R` orchestrates comprehensive testing:
1. Generates synthetic data with known functions
2. Fits both spline types with various configurations
3. Extracts posterior draws and constructs uncertainty intervals
4. Creates multi-panel plots using patchwork
5. Saves all outputs to `output/` directory

## Critical Stan Requirements

### Array Syntax
Always use modern array syntax:
```stan
array[N] real x;           // Correct
real x[N];                 // Wrong - old syntax
```

### Array Size Constraints
In generated quantities, array sizes must be literals:
```stan
generated quantities {
  array[100] real x_plot;  // Correct - literal size
  array[n_plot] real x_plot; // Wrong - variable size
}
```

### File Formatting
All Stan files must end with a newline character.

## Package Dependencies

Managed via groundhog with date "2025-06-01":
- cmdstanr (v0.9.0+) - installed from Stan repo
- Core packages: posterior, dplyr, tidyr, ggplot2, patchwork
- Namespace conflicts resolved via conflicted library

## Important Constraints

1. **No Extrapolation**: Both spline types should not extrapolate beyond data range
2. **Fixed Knots**: Knot positions determined by data quantiles, not optimized
3. **No Shrinkage**: Neither implementation includes automatic regularization
4. **Boundary Behavior**: 
   - B-splines: Polynomial continuation at boundaries
   - C-splines: Linear extrapolation (natural boundary conditions)
5. **C-spline Strictness**: The stan-splines library rejects points outside knot range