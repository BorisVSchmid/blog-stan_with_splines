# Stan Splines Test Suite

This directory contains comprehensive tests for the B-spline and C-spline implementations.

## Test Organization

### Core Tests (Run by default)
- **test_basic_splines.R** - Basic functionality checks for both B-splines and C-splines
- **test_numerical_accuracy.R** - Mathematical properties tests (partition of unity, linear interpolation)
- **test_analytical_solutions.R** - Known function fitting tests

### Diagnostic Tests
- **test_diagnostics_both_splines.R** - Tests the smoothing diagnostics for both spline types
- **test_edf_knot_response.R** - Tests effective degrees of freedom calculation
- **test_smoothing_strength.R** - Demonstrates smoothing_strength parameter behavior

### Special Case Tests
- **test_3_knots.R** - Tests minimal knot configuration

### Complex Tests (Run separately)
- **test_splines.R** - Comprehensive tests with multiple function types
- **test_regional_splines.R** - Hierarchical model tests
- **test_regional_splines_simple.R** - Simplified regional model tests

## Running Tests

### Quick validation
```bash
Rscript run-tests.R
```

### Simple test without groundhog dependencies
```bash
cd tests
Rscript run_tests_simple.R
```

### Individual test
```bash
cd tests
Rscript test_basic_splines.R
```

## Stan Model Files

The tests use Stan models from the main `code/` directory:
- `../code/bsplines.stan` - B-spline implementation
- `../code/csplines.stan` - C-spline implementation

The only Stan file kept in this directory is test_regional_splines.stan which is used by the hierarchical model tests.