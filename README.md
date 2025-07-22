# Stan Spline Implementations: B-splines and C-splines

This repository contains minimal implementations for testing B-splines and natural cubic splines (C-splines) in Stan, with comprehensive testing and comparison tools.

**Built with [Claude Code](https://claude.ai/code)** - An AI-powered coding assistant that helped create, test, and document this project.

## Overview

Two spline implementations are provided:
- **B-splines** (`test_bsplines.stan`): Based on Stan documentation and case studies
- **C-splines** (`test_csplines.stan`): Natural cubic splines using the stan-splines library

## Files

### Core Implementations
- `test_bsplines.stan` - B-spline implementation
- `test_csplines.stan` - Natural cubic spline implementation  
- `spline.stan` - Natural cubic spline functions from stan-splines library

### Testing Scripts
- `test_splines.R` - Comprehensive testing script
- `quick_test.R` - Simple test for quick verification
- `generate_plots_extended.R` - Extended boundary behavior visualization

### Regional/Hierarchical Models
- `test_regional_splines.stan` - Full hierarchical model with smoothness penalties
- `test_regional_splines_v2.stan` - Simplified hierarchical model
- `run_regional_splines.R` - Complete regional analysis with visualizations
- `demo_regional_splines.R` - Simple demonstration of regional concepts

### Output
- `output/` - Directory for all generated plots and results

## Usage

Three main entry points are provided:

```r
# Run minimal code examples (simple demonstrations)
source("run-code.R")

# Run extended examples (complex analyses, comparisons)
source("run-examples.R")

# Run test suite (verify implementations)
source("run-tests.R")       # Core tests only
source("run-all-tests.R")   # Complete test suite
```

For direct access to specific functionality:
```r
# Minimal examples
source("code/bsplines.R")   # B-spline minimal example
source("code/csplines.R")   # C-spline minimal example

# Extended examples
source("examples/compare_smoothing.R")      # B-spline smoothing comparison
source("examples/run_regional_splines.R")   # Regional hierarchical models
```

## Capabilities Comparison

### B-splines

**Advantages:**
- ✅ Local support - changes to data affect only nearby basis functions
- ✅ Flexible degree control (linear, quadratic, cubic, etc.)
- ✅ Can represent polynomial functions exactly (up to degree)
- ✅ Stable numerical properties
- ✅ Flexible boundary behavior

**Limitations:**
- ❌ No automatic shrinkage/smoothing
- ❌ Knot positions are fixed (quantiles of data)
- ❌ No built-in warnings for poor coverage
- ❌ More parameters than C-splines for same flexibility

**Implementation Details:**
- Uses recursive Cox-de Boor algorithm
- Knots placed at quantiles of input data
- Extended knots repeat boundary knots `degree` times
- Number of basis functions = `num_knots + degree - 1`

### C-splines (Natural Cubic Splines)

**Advantages:**
- ✅ Natural boundary conditions (linear beyond endpoints)
- ✅ Fewer parameters than B-splines
- ✅ Optimal interpolating spline (minimizes curvature)
- ✅ Smooth second derivatives
- ✅ More efficient for smooth functions

**Limitations:**
- ❌ Global support - changes affect entire curve
- ❌ Fixed cubic degree
- ❌ No automatic shrinkage
- ❌ Knot positions are fixed
- ❌ Implementation has strict boundary checking

**Implementation Details:**
- Uses matrix algebra for coefficient computation
- Natural boundary conditions (second derivative = 0 at endpoints)
- Requires all evaluation points within knot range
- Based on stan-splines library by Sergey Koposov

## Key Differences

| Feature | B-splines | C-splines |
|---------|-----------|-----------|
| **Degree flexibility** | Any degree (2, 3, 4, ...) | Fixed cubic (degree 3) |
| **Support** | Local | Global |
| **Boundary behavior** | Polynomial continuation | Linear at boundaries |
| **Parameters** | More (num_basis coefficients) | Fewer (num_knots values) |
| **Smoothness** | C^(degree-1) continuous | C^2 continuous |
| **Computation** | Recursive evaluation | Matrix solution |
| **Best for** | Flexible local fitting | Smooth global curves |
| **Extrapolation** | Not recommended | Not recommended |

## Shrinkage and Regularization

Neither implementation includes automatic shrinkage or regularization:

- **No built-in penalties** on roughness or complexity
- **No adaptive smoothing** parameters
- Smoothness controlled only by:
  - Number of knots (fewer = smoother)
  - Degree (for B-splines)
  - Prior distributions on coefficients

To add shrinkage, you could:
1. Add a smoothing penalty in the model block
2. Use hierarchical priors on spline coefficients
3. Implement penalized splines (P-splines)

## Knot Placement

Both implementations use **fixed knot positions**:

- Knots placed at quantiles of input data
- No optimization of knot locations
- No adaptive knot placement

For B-splines:
- Interior knots at data quantiles
- Boundary knots extended beyond data range

For C-splines:
- Knots must encompass all data points
- Small buffer added to prevent boundary errors

## Coverage Warnings

Neither implementation provides automatic warnings for:
- ❌ Too few knots for data complexity
- ❌ Poor knot placement
- ❌ Overfitting with too many knots

Best practices:
- Start with 4-6 knots for simple curves
- Use 8-12 knots for complex patterns
- Check residual plots for systematic patterns
- Compare different knot numbers
- Keep predictions within the data range
- Avoid extrapolation with splines

## Performance

From testing on a 100-point dataset with 8 knots:
- B-splines and C-splines have similar computation time
- Both scale well with data size
- C-splines slightly faster for simple curves
- B-splines more stable for complex patterns

## Recommendations

**Use B-splines when:**
- You need local control
- Degree flexibility is important
- Working with non-smooth data
- Need to minimize parameters affected by local changes

**Use C-splines when:**
- Smoothness is paramount
- Working within data range
- Fewer parameters preferred
- Natural boundary behavior desired

## Example Results

See `output/` directory for:
- `spline_comparison.pdf` - Side-by-side comparisons
- `bspline_basis_functions.pdf` - Visualization of basis functions
- `bspline_degree_comparison.pdf` - Effect of spline degree
- `diagnostics_summary.rds` - Full MCMC diagnostics
- `timing_results.csv` - Performance benchmarks

## Regional/Hierarchical Splines

The repository now includes implementations for fitting splines to multiple regions with shared hierarchical priors:

### Model Structure
```stan
// Each region has its own spline coefficients
array[n_regions] row_vector[num_basis] alpha;

// Coefficients share hierarchical priors
row_vector[num_basis] mu_alpha;        // Global mean
vector<lower=0>[num_basis] tau_alpha;  // Global SD

// Regional coefficients centered around global mean
alpha[r] ~ normal(mu_alpha, tau_alpha);
```

### Key Features
- **Partial pooling**: Regions borrow strength from each other
- **Regional baselines**: Each region can have different intercepts
- **Variance decomposition**: ICC shows between vs within-region variance
- **Optional smoothness penalties**: Control wiggliness per region

### Example Use Cases
- Disease rates across regions with similar seasonal patterns
- Economic indicators across countries with shared trends
- Environmental measurements across sites with common processes

### Running Regional Models
```r
# Simple demonstration
source("demo_regional_splines.R")

# Full analysis with visualizations
source("run_regional_splines.R")
```

## Future Enhancements

Potential improvements not currently implemented:
1. **Adaptive knot placement** - Optimize knot locations
2. **Penalized splines** - Add roughness penalties
3. **Automatic smoothing** - Cross-validation for smoothing parameters
4. **Coverage diagnostics** - Warn about poor knot coverage
5. **Monotonic constraints** - Ensure monotonic fits
6. **Periodic splines** - For cyclic data
7. **Tensor product splines** - For 2D/3D smoothing

## Attribution and Licensing

See [ATTRIBUTION.md](ATTRIBUTION.md) for important information about third-party code and required citations.

## References

- B-splines: [Stan Case Study](https://mc-stan.org/learn-stan/case-studies/splines_in_stan.html) by Milad Kharratzadeh
- C-splines: [stan-splines library](https://github.com/segasai/stan-splines) by Sergey Koposov
  - **Required citation**: Koposov et al. (2019), MNRAS, 485, 4726
- General reference: "The Elements of Statistical Learning" by Hastie, Tibshirani, and Friedman

## Development

This project was developed with [Claude Code](https://claude.ai/code), an AI coding assistant that helped:
- Extract and adapt minimal spline implementations from complex outbreak models
- Create comprehensive test suites for both B-splines and C-splines
- Implement smoothing/regularization for B-splines using random walk priors
- Develop regional spline models with hierarchical priors
- Generate visualization code with proper uncertainty quantification
- Document capabilities, limitations, and usage patterns
- Ensure proper attribution and licensing compliance

## Testing

This project includes a comprehensive test suite to verify the correctness and properties of both spline implementations.

### Test Files

| Test Script | Purpose | Coverage |
|-------------|---------|----------|
| `tests/quick_test.R` | Basic functionality check | Quick verification of both implementations |
| `tests/test_numerical_accuracy.R` | Mathematical properties | Partition of unity, interpolation accuracy, monotonicity |
| `tests/test_analytical_solutions.R` | Known function fitting | Polynomials, constants, step functions, sine waves |
| `tests/test_splines.R` | Comprehensive function tests | Multiple function types and edge cases |
| `tests/test_regional_splines.R` | Hierarchical models | Multi-region splines with shared priors |

### Running Tests

```r
# Run all tests with summary
source("run-tests.R")

# Run individual test suites
source("tests/quick_test.R")                    # Fast basic check
source("tests/test_numerical_accuracy.R")       # Mathematical properties
source("tests/test_analytical_solutions.R")     # Known solutions
source("tests/test_splines.R")                  # Comprehensive tests
source("tests/test_regional_splines.R")         # Hierarchical models
```

### Test Coverage

Our test suite covers:

#### Mathematical Properties
- **Partition of Unity**: B-splines sum to 1 at any point
- **Interpolation Accuracy**: Closeness to data points
- **Linear Function Fitting**: Perfect fit for linear data
- **Monotonicity Preservation**: Maintaining trends in data

#### Analytical Solutions
- **Polynomial Fitting**: Cubic splines fitting cubic polynomials
- **Constant Functions**: Both splines fitting flat data
- **Step Function Approximation**: Handling discontinuous data
- **Sine Wave Fitting**: Smooth periodic function approximation

#### Edge Cases and Robustness
- **Minimal Data**: Testing with few data points (n=5)
- **Different Spline Degrees**: B-splines with degree 2, 3, 4
- **Various Function Types**: Linear, polynomial, exponential, complex
- **Boundary Behavior**: Knot placement and extrapolation limits

#### Performance and Comparison
- **Timing Comparisons**: B-splines vs C-splines execution time
- **Regional Models**: Hierarchical splines for multiple regions
- **Parameter Sensitivity**: Different knot numbers and smoothing levels

### Test Philosophy

The tests follow the validation approach of the original repositories:

- **C-splines**: Based on numerical accuracy tests from [segasai/stan-splines](https://github.com/segasai/stan-splines)
- **B-splines**: Extended beyond original case study to include comprehensive validation
- **Pure R implementation**: All tests run in R without external dependencies

Tests are designed to be:
- **Robust** to MCMC stochasticity
- **Fast** enough for regular validation
- **Comprehensive** in coverage
- **Clear** in output and failure reporting

## TODO

### Comparison with R's splines package
- Validate B-spline implementation against `splines::bs()`
- Compare C-spline results with `splines::ns()` 
- Benchmark performance differences
- Document any discrepancies in knot placement or boundary behavior

### Performance Optimizations
- Implement iterative B-spline algorithm to avoid recursion overhead
- Add memoization for repeated basis function evaluations
- Profile and optimize for large datasets (10,000+ points)

### Enhanced Functionality
- Add derivative computation for both spline types
- Implement alternative boundary conditions for B-splines (not-a-knot, clamped)
- Add monotonic and shape-constrained spline options
- Support for periodic/cyclic splines

### Model Selection and Diagnostics
- Implement automatic knot selection via cross-validation
- Add WAIC and LOO-CV computation for model comparison
- Create diagnostic plots for knot placement effectiveness
- Add residual analysis tools

### Documentation and Examples
- Create vignette comparing spline types for different use cases
- Add examples with real-world datasets
- Document computational complexity and memory usage
- Create interactive Shiny app for spline exploration