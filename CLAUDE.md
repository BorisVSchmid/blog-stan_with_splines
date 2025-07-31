# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About This Project

This project was built collaboratively with Claude Code to:
- Extract minimal, working spline implementations from complex epidemiological models
- Create comprehensive test frameworks for B-splines and C-splines in Stan
- Implement advanced features like smoothing priors and regional hierarchical models
- Document best practices and usage patterns for Stan splines

The development process involved iterative refinement, with Claude Code helping to:
- Debug Stan compilation issues and syntax requirements
- Implement proper uncertainty quantification
- Create visualizations using ggplot2 with appropriate themes
- Ensure proper attribution and licensing compliance
- Organize code following R package conventions

## Project Structure

- `code/` - Core Stan models and R functions
- `tests/` - Test scripts for verifying implementations
- `examples/` - Example analyses and demonstrations
- `output/` - Generated plots and results
- `claude-tmp/` - Temporary/experimental scripts created by Claude (not part of main project)

## Commands

### Running the Project
```r
# Run minimal code examples (simple demonstrations)
source("run-code.R")

# Run extended examples (complex analyses, comparisons)
source("run-examples.R")

# Run complete test suite (all 5 tests)
source("run-tests.R")
```

### Running Individual Components
```r
# Core implementations
source("code/bsplines.R")   # B-spline minimal example
source("code/csplines.R")   # C-spline minimal example

# Run specific tests
source("tests/test_basic_splines.R")              # Quick functionality check
source("tests/test_analytical_solutions.R")       # Known function fitting
source("tests/test_diagnostic_recommendations.R") # Smoothing diagnostics

# Run examples
source("examples/hierarchical_regional_splines_adaptive.R")
```

### Building Stan Models
Models are compiled automatically when running scripts. To manually compile:
```r
library(cmdstanr)
model_b <- cmdstan_model("code/bsplines.stan")
model_c <- cmdstan_model("code/csplines.stan")
model_regional_adaptive <- cmdstan_model("examples/regional_splines_adaptive.stan")
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
   - Smoothing formula: `tau = prior_scale / (smoothing_strength * num_basis^sqrt(2))`
   - Default smoothing_strength = 0.1 (mild smoothing)

2. **C-splines** (`code/csplines.stan`)
   - Depends on `spline.stan` library via `#include`
   - Uses matrix algebra for natural cubic splines
   - Requires all evaluation points within knot range
   - Functions: `spline_getcoeffs()`, `spline_eval()`, `spline_findpos()`

### Key Architectural Differences

- **B-splines**: Local support means coefficients (`alpha`) directly multiply basis functions
- **C-splines**: Global support with values specified at knots (`y_at_knots`), then transformed to coefficients

### Default Knot Selection

The project uses adaptive knot selection as defaults:
- **B-splines**: `max(4, min(round(n/2), 40))` - half the data points, capped at 40
- **C-splines**: `max(4, min(round(n/4), 20))` - quarter of data points, capped at 20

These defaults balance flexibility with stability across different sample sizes.

### MCMC Configuration

All models now run with:
- 4 chains in parallel (`chains = 4, parallel_chains = 4`)
- Diagnostic summaries printed after each fit (`fit$diagnostic_summary()`)
- Adaptive parameters for complex models (`adapt_delta`, `max_treedepth`)

### Testing Framework

The test suite (`run-tests.R`) includes 6 tests in order of complexity:
1. `test_basic_splines.R` - Basic functionality check
2. `test_edge_cases_minimal_knots.R` - Minimal knot configuration (2-3 knots)
3. `test_numerical_accuracy.R` - Mathematical properties
4. `test_analytical_solutions.R` - Known function fitting
5. `test_sine_wave_smoothing.R` - B-spline smoothing effects on sine waves
6. `test_diagnostic_recommendations.R` - Tests smoothing diagnostics with multiple scenarios

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

## Smoothing Parameter Guide

### B-spline Smoothing Strength
- **0**: No smoothing (independent coefficients)
- **0.05-0.1**: Mild smoothing (0.1 is default)
- **0.1-0.2**: Strong smoothing

The smoothing uses the formula `tau = prior_scale / (smoothing_strength * num_basis^sqrt(2))`, which provides:
- Scale invariance: consistent behavior regardless of y-axis scale
- Knot-count invariance: consistent smoothing across different numbers of basis functions

### Recent Updates (2025)
- Changed smoothing formula from `sqrt(num_basis)` to `num_basis^sqrt(2)` for better scaling
- Updated default smoothing_strength from 1.0 to 0.1
- Adjusted recommended ranges to reflect the new scaling

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
- `fit_with_diagnostics()`: Wrapper function for fitting with automatic diagnostics
- `plot_diagnostic_residuals()`: Creates diagnostic residual plots

Additional diagnostic documentation:
- `code/DIAGNOSTIC_INTERPRETATION.md`: Guide for interpreting diagnostic plots
- `code/IMPLEMENTATION_DIFFERENCES.md`: Differences from original spline implementations

### Output Organization
- All plots saved to `output/` directory with consistent naming:
  - `code-*.png`: Output from code examples
  - `test-*.png`: Output from test scripts
  - `example-*.png`: Output from example scripts
- No automatic PDF generation (removed `print()` calls that created Rplots.pdf)

## Regional/Hierarchical Models

The `examples/` directory contains advanced hierarchical implementations:
- `regional_splines_adaptive.stan`: Hierarchical B-splines with adaptive shrinkage
- `hierarchical_regional_splines_adaptive.R`: Example demonstrating the adaptive model
- Features:
  - Separate smoothing strengths for global vs regional patterns
  - Adaptive shrinkage factor that learns from data
  - Variance component estimation (between vs within-region)
  - Better separation of shared and regional components
  - Model caching for repeated runs
  - Comprehensive visualization of regional patterns