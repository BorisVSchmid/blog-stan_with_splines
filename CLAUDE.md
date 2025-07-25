# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About This Project

This project was built collaboratively with Claude Code to:
- Extract minimal, working spline implementations from complex epidemiological models
- Create comprehensive test frameworks for B-splines and C-splines in Stan
- Implement advanced features like smoothing priors and regional hierarchical models
- Document best practices and usage patterns for Stan splines

## Project Structure

- `code/` - Core Stan models and R functions (bsplines.stan, csplines.stan, smoothing_diagnostics.R)
- `tests/` - Test scripts for verifying implementations (9 test files)
- `examples/` - Example analyses and demonstrations (hierarchical regional splines)
- `output/` - Generated plots and results (all test outputs saved here)

## Commands

### Running the Project
```r
# Run minimal code examples (simple demonstrations)
source("run-code.R")

# Run extended examples (complex analyses, comparisons)
source("run-examples.R")

# Run complete test suite (all 9 tests)
source("run-tests.R")
```

### Running Individual Components
```r
# Core implementations
source("code/bsplines.R")   # B-spline minimal example
source("code/csplines.R")   # C-spline minimal example

# Run specific tests
source("tests/test_basic_splines.R")         # Quick functionality check
source("tests/test_splines.R")               # Comprehensive testing
source("tests/test_regional_splines.R")      # Hierarchical models

# Run examples
source("examples/hierarchical_regional_splines.R")
```

### Building Stan Models
Models are compiled automatically when running scripts. To manually compile:
```r
library(cmdstanr)
model_b <- cmdstan_model("code/bsplines.stan")
model_c <- cmdstan_model("code/csplines.stan")
model_regional <- cmdstan_model("examples/regional_splines.stan")
```

## Architecture

### Spline Implementations

The project provides two independent spline implementations in Stan:

1. **B-splines** (`code/bsplines.stan`)
   - Self-contained implementation using Cox-de Boor recursive algorithm
   - Function `build_b_spline()` constructs basis functions recursively
   - Knots placed at quantiles with extended boundary knots
   - Number of basis functions = `num_knots + spline_degree - 1`
   - Includes `smoothing_strength` parameter for scale-invariant smoothing

2. **C-splines** (`code/csplines.stan`)
   - Depends on `spline.stan` library via `#include`
   - Uses matrix algebra for natural cubic splines
   - Requires all evaluation points within knot range
   - Functions: `spline_getcoeffs()`, `spline_eval()`, `spline_findpos()`

### Key Architectural Differences

- **B-splines**: Local support means coefficients (`alpha`) directly multiply basis functions
- **C-splines**: Global support with values specified at knots (`y_at_knots`), then transformed to coefficients

### MCMC Configuration

All models now run with:
- 4 chains in parallel (`chains = 4, parallel_chains = 4`)
- Diagnostic summaries printed after each fit (`fit$diagnostic_summary()`)
- Adaptive parameters for complex models (`adapt_delta`, `max_treedepth`)

### Testing Framework

The test suite (`run-tests.R`) includes 9 tests in order of complexity:
1. `test_basic_splines.R` - Basic functionality check
2. `test_3_knots.R` - Minimal knot configuration
3. `test_smoothing_strength.R` - Smoothing parameter effects
4. `test_numerical_accuracy.R` - Mathematical properties
5. `test_analytical_solutions.R` - Known function fitting
6. `test_edf_knot_response.R` - Effective degrees of freedom
7. `test_diagnostics_both_splines.R` - Diagnostic comparison
8. `test_splines.R` - Comprehensive scenarios
9. `test_regional_splines.R` - Hierarchical models

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
3. **No Automatic Shrinkage**: Neither implementation includes automatic regularization (B-splines have optional smoothing_strength)
4. **Boundary Behavior**: 
   - B-splines: Polynomial continuation at boundaries
   - C-splines: Linear extrapolation (natural boundary conditions)
5. **C-spline Strictness**: The stan-splines library rejects points outside knot range

## Diagnostics and Output

### Smoothing Diagnostics
The `code/smoothing_diagnostics.R` file provides:
- `diagnose_smoothing()`: Calculates EDF, residual patterns, cross-validation metrics
- `print_smoothing_diagnostics()`: Pretty-prints diagnostic results with recommendations
- Automated parameter tuning suggestions based on diagnostic metrics

### Output Organization
- All plots saved to `output/` directory with consistent naming:
  - `code-*.png`: Output from code examples
  - `test-*.png`: Output from test scripts
  - `example-*.png`: Output from example scripts
- No automatic PDF generation (removed `print()` calls that created Rplots.pdf)

## Regional/Hierarchical Models

The `examples/` directory contains advanced hierarchical implementations:
- `regional_splines.stan`: Hierarchical B-splines with region-specific deviations
- Variance component estimation (between vs within-region)
- Model caching for repeated runs
- Comprehensive visualization of regional patterns