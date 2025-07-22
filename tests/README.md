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

## Test Files (10 total)

### Core Functionality Tests (2)
1. **test_basic_splines.R** - Essential functionality for both B-splines and C-splines
2. **test_3_knots.R** - Edge case testing with minimal knot configuration

### Mathematical Property Tests (2)
3. **test_numerical_accuracy.R** - Verifies mathematical properties (partition of unity, interpolation)
4. **test_analytical_solutions.R** - Tests against known analytical solutions (constant, cubic, sine)

### Parameter and Feature Tests (3)
5. **test_smoothing_strength.R** - Tests the smoothing_strength parameter behavior
6. **test_edf_knot_response.R** - Tests effective degrees of freedom calculations  
7. **test_diagnostics_both_splines.R** - Tests diagnostic functionality

### Advanced Tests (3)
8. **test_splines.R** - Extensive tests with multiple function types and performance comparison
9. **test_regional_splines_simple.R** - Simplified regional model for testing
10. **test_regional_splines.R** - Full hierarchical model implementation

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
- Pass/fail status for each test
- Timing information
- Summary statistics

Example output:
```
[1/10] Running test_basic_splines.R
--------------------------------------------------
✓ test_basic_splines.R completed in 12.3 seconds
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