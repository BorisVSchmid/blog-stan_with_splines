# Stan Splines Test Suite

Comprehensive test suite for B-spline and C-spline implementations in Stan.

## Running Tests

### Run all tests
From the project root directory:
```bash
Rscript run-tests.R
```

### Run individual test
```bash
cd tests
Rscript test_basic_splines.R
```

## Test Files (5 total)

### Core Functionality Tests (2)
1. **test_basic_splines.R** - Essential functionality for both B-splines and C-splines
   - Fits simple sine function with both spline types
   - Creates comparison plot showing fitted curves with uncertainty intervals
   - Outputs: `test-basic_spline_comparison.png`

2. **test_edge_cases_minimal_knots.R** - Edge case testing with minimal knot configurations
   - Tests both B-splines and C-splines with 2 and 3 knots
   - Verifies error handling for insufficient knots
   - Shows limitations of minimal flexibility
   - Outputs: `test-edge_cases_2knots_bspline.png`, `test-edge_cases_3knots_comparison.png`

### Mathematical Property Tests (2)
3. **test_numerical_accuracy.R** - Verifies mathematical properties
   - Partition of unity test for B-splines
   - Linear function interpolation for both spline types
   - Interpolation accuracy at data points
   - Monotonicity preservation
   - Outputs: `test-numerical_accuracy_properties.png`

4. **test_analytical_solutions.R** - Tests against known analytical solutions
   - Polynomials (degree 1-3)
   - Constant functions
   - Sine functions
   - Tests both spline types against each function

### Parameter Test (1)
5. **test_sine_wave_smoothing.R** - B-spline smoothing effects on sine waves
   - Tests different smoothing values (0, 0.05, 0.10, 0.15, 0.20)
   - Shows how smoothing affects oscillatory function fitting
   - Compares with C-splines at different knot counts
   - Creates 5x4 grid visualization
   - Outputs: `test-sine_wave_smoothing.png`

## Test Categories

1. **Quick Tests** (<5 seconds each)
   - Basic functionality
   - Minimal configurations
   - Parameter validation

2. **Mathematical Tests** (5-20 seconds each)
   - Numerical accuracy
   - Analytical solutions
   - Mathematical properties

3. **Comprehensive Tests** (>20 seconds each)
   - Multiple scenarios
   - Regional/hierarchical models
   - Extensive parameter sweeps

## Dependencies

Tests use the main implementations:
- `code/bsplines.stan` - B-spline implementation with smoothing
- `code/csplines.stan` - Natural cubic spline implementation
- `code/smoothing_diagnostics.R` - Diagnostic functions

The only Stan file in the tests directory is:
- `test_regional_splines.stan` - Hierarchical spline model for testing

## Expected Output

Successful test runs will show:
- Number of tests per file
- Purpose/description of each test
- Pass/fail status for each test  
- Timing information
- MCMC diagnostic summaries after each model fit
- Detailed results table with all 5 tests
- Summary statistics

### Generated Plots

All tests now generate visualization plots saved to the `output/` directory:

1. **Basic functionality**: `test-basic_spline_comparison.png`
2. **Edge cases**: `test-edge_cases_2knots_bspline.png`, `test-edge_cases_3knots_comparison.png`
3. **Numerical accuracy**: `test-numerical_accuracy_properties.png`
4. **Flexibility comparison**: `test-flexibility_comparison.png`, `test-flexibility_comparison_summary.png`
5. **Diagnostic recommendations**: `test-diagnostic_recommendations_scenarios.png`, `test-diagnostic_metrics_comparison.png`
6. **Various target functions**: `test-spline_comparison.png`, `test-diagnostics_summary.csv`
7. **Regional splines**: Multiple plots showing variance components and regional fits

### MCMC Diagnostics
All test files now output MCMC diagnostic summaries after model fitting, showing:
- Number of divergent transitions per chain
- Max treedepth hits per chain
- EBFMI (Energy Bayesian Fraction of Missing Information) values
- Any diagnostic warnings

Example diagnostic output:
```
MCMC Diagnostic Summary:
$num_divergent
[1] 0 0 0 0

$num_max_treedepth
[1] 0 0 0 0

$ebfmi
[1] 1.234 1.156 1.089 1.201
```

### Parallel Chains
All models are now configured to run 4 chains in parallel for better sampling and faster execution.

Example output:
```
[1/8] Running test_basic_splines.R
Purpose: Basic functionality of B-splines and C-splines
----------------------------------------------------------------------
Testing B-spline model...

B-spline MCMC diagnostics:
$num_divergent
[1] 0 0 0 0

...

✓ test_basic_splines.R completed in 12.3 seconds

...

======================================================================
TEST SUITE SUMMARY
======================================================================

Total tests run: 8
Passed: 8
Failed: 0
Total time: 352.5 seconds (5.9 minutes)
```

## Troubleshooting

If tests fail:
1. Check that cmdstanr is properly installed
2. Verify Stan models compile correctly
3. Check for missing R packages
4. Review specific error messages in test output

Common issues:
- **Compilation errors**: Update cmdstanr or check Stan syntax
- **Numeric differences**: May need to adjust tolerance in tests
- **Missing functions**: Source required files or install packages