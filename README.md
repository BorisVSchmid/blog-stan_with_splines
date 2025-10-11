# Stan Spline Implementations: B-splines and C-splines

This repository contains minimal implementations for testing B-splines and natural cubic splines (C-splines) in Stan, with comprehensive testing and comparison tools.

**Built with [Claude Code](https://claude.ai/code)** - An AI-powered coding assistant that helped create, test, and document this project. So be careful there.
I spend quite some time with this code, so I trust it to a decent degree. There is an accompanying blog post on boris.earth.

## Overview

Two spline implementations are provided:
- **B-splines** (`code/bsplines.stan`): Based on Stan documentation with scale-invariant smoothing
- **C-splines** (`code/csplines.stan`): Natural cubic splines using the stan-splines library

## Files

### Core Implementations
- `code/bsplines.stan` - B-spline implementation with scale-invariant smoothing
- `code/csplines.stan` - Natural cubic spline implementation  
- `code/spline.stan` - Natural cubic spline functions from stan-splines library (dependency)
- `code/smoothing_diagnostics.R` - Diagnostic tools for detecting over/under-smoothing

### Testing Suite (8 tests)
- `tests/test_basic_splines.R` - Basic functionality of both spline types
- `tests/test_edge_cases_minimal_knots.R` - Edge cases with 2-3 knots
- `tests/test_flexibility_bspline_vs_cspline.R` - Smoothing comparison
- `tests/test_numerical_accuracy.R` - Mathematical properties verification
- `tests/test_analytical_solutions.R` - Known function fitting
- `tests/test_diagnostic_recommendations.R` - Diagnostic system testing
- `tests/test_various_target_functions.R` - Six different target functions
- `tests/test_regional_splines.R` - Hierarchical regional models

### Example Scripts
- `examples/hierarchical_regional_splines.R` - Advanced hierarchical decomposition
- `examples/regional_splines.stan` - Stan model for regional hierarchical splines

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
source("run-tests.R")
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

## When to Use B-splines vs C-splines

### Choose B-splines when you need:

**1. Local Control**
- Changes to data in one region shouldn't affect the entire curve
- Modeling data with local features, sharp peaks, or discontinuities
- Robustness to outliers (they only affect nearby regions)
- Example: Epidemic curves with distinct outbreak peaks

**2. Flexible Degree Selection**
- Linear (degree 1) for piecewise linear fits
- Quadratic (degree 2) for simple smooth curves
- Cubic (degree 3) for standard smooth modeling
- Higher degrees for special requirements

**3. Explicit Smoothing Control**
- The `smoothing_strength` parameter provides fine control
- 0 = no smoothing (independent coefficients)
- 0.05-0.1 = mild smoothing (0.1 is default)
- 0.1-0.2 = strong smoothing
- Scale-invariant smoothing that adapts to the y-axis standard deviation and the number of basis functions (num_basis^sqrt(2))
- Example: Adjusting smoothness based on data quality or sample size

**4. Hierarchical Modeling**
- Local basis functions work well with region-specific effects
- Each region can have different coefficient values
- Natural for multi-level models

### Choose C-splines when you need:

**1. Natural Smoothness**
- Automatic C² continuity (smooth second derivatives)
- Optimal in the sense of minimizing total curvature
- Natural boundary conditions with linear extrapolation
- Example: Long-term mortality trends or baseline hazards

**2. Parsimonious Models**
- Fewer parameters to estimate (one per knot)
- More efficient for simple smooth curves
- Better when you have limited data
- Example: Age-specific patterns that change gradually

**3. Global Smoothness Constraints**
- Ensures smooth transitions across entire domain
- No risk of artificial wiggles between knots
- Appropriate for inherently smooth processes
- Example: Temperature trends or growth curves

### Practical Examples for Epidemiology

**B-splines are ideal for:**
```r
# Modeling epidemic curves with potential sharp features
epidemic_fit <- fit_bspline(outbreak_data, 
                           num_knots = 15,        # More knots for flexibility
                           smoothing_strength = 0.1)  # Default mild smoothing

# Regional models with local variation
regional_model <- stan("regional_splines.stan", 
                       data = list(spline_degree = 3,
                                  num_knots = 10))
```

**C-splines are ideal for:**
```r
# Smooth baseline mortality trends
baseline_fit <- fit_cspline(mortality_data,
                           num_knots = 7)  # Fewer knots needed

# Age-specific infection fatality rates
age_ifr_fit <- fit_cspline(age_data,
                          num_knots = 5)  # Smooth age progression
```

### Mixed Approaches

In complex models, you might use both:
- B-splines for the main time-varying effects (with local features)
- C-splines for smooth confounders or baseline trends
- Example: COVID-19 model with B-splines for daily cases and C-splines for age effects

## Smoothing and Regularization

### B-splines: Scale-Invariant Smoothing

B-splines now include a sophisticated scale-invariant smoothing system:

```r
# The smoothing parameter adapts to your data automatically
tau_smooth = prior_scale / (smoothing_strength * num_basis^sqrt(2))
```

**Key features:**
- **Scale-invariant**: Adapts to the standard deviation of your y-values
- **Basis-aware**: Accounts for the number of basis functions
- **Intuitive control**: 
  - `smoothing_strength = 0`: No smoothing (independent coefficients)
  - `smoothing_strength = 0.05-0.1`: Mild smoothing (0.1 is recommended default)
  - `smoothing_strength = 0.1-0.2`: Strong smoothing

The implementation uses a random walk prior on coefficients:
```stan
alpha[i] ~ normal(alpha[i-1], tau_smooth)
```

### C-splines: Knot-Based Smoothing

C-splines control smoothness through the number of knots:
- Fewer knots = smoother curves
- More knots = more flexible fitting
- Natural boundary conditions ensure smooth behavior

### Diagnostic System

The package includes `code/smoothing_diagnostics.R` which provides:
- Automatic detection of over/under-smoothing
- Context-aware recommendations
- Key metrics: residual autocorrelation, runs test, EDF approximation
- Actionable parameter adjustments based on diagnostics

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

See `output/` directory for generated plots:

### Code Examples
- `code-bspline_minimal_example.png` - Basic B-spline fit
- `code-cspline_minimal_example.png` - Basic C-spline fit

### Test Outputs  
- `test-basic_spline_comparison.png` - Side-by-side comparison
- `test-edge_cases_*.png` - Minimal knot configurations
- `test-flexibility_comparison*.png` - Smoothing parameter effects
- `test-numerical_accuracy_properties.png` - Mathematical properties
- `test-diagnostic_recommendations_*.png` - Diagnostic scenarios
- `test-spline_comparison.png` - Six function types comparison
- `test-diagnostics_summary.csv` - Diagnostic metrics summary

### Example Outputs
- `example-hierarchical_decomposition.png` - Regional pattern decomposition
- `example-hierarchical_model_draws.rds` - Cached model draws

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
# Hierarchical regional spline example
source("examples/hierarchical_regional_splines.R")

# Test hierarchical models
source("tests/test_regional_splines.R")
```

## Using Splines in Your Stan Models

This repository provides two ways to incorporate splines into your Stan models:

### Option 1: Copy-paste the functions (Recommended for B-splines)

For B-splines, you can copy the `build_b_spline` function directly into your Stan model:

```stan
functions {
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order, int degree) {
    // Copy the function from code/bsplines.stan
  }
}

data {
  // Your data
  int<lower=0> n_data;
  array[n_data] real x;
  vector[n_data] y;
  
  // Spline parameters
  int<lower=1> num_knots;
  int<lower=1> spline_degree;
  real<lower=0> smoothing_strength;
  real<lower=0> prior_scale;
}

transformed data {
  // Set up B-spline basis (see code/bsplines.stan for full implementation)
  int num_basis = num_knots + spline_degree - 1;
  // ... knot placement and basis construction ...
}

parameters {
  real alpha_0;  // Intercept
  vector[num_basis] alpha_raw;  // Spline coefficients
  // ... other parameters ...
}

model {
  // Use the spline in your model
  vector[n_data] y_hat = alpha_0 + to_vector(alpha' * B);
  y ~ normal(y_hat, sigma);
}
```

### Option 2: Include the spline library (Required for C-splines)

For C-splines, you must include the stan-splines library:

```stan
functions {
  #include spline.stan
}

data {
  int<lower=0> n_data;
  array[n_data] real x;
  vector[n_data] y;
  int<lower=1> num_knots;
}

transformed data {
  // Set up knots and get spline coefficients
  vector[num_knots] knots = // ... set knot positions ...
  matrix[4, num_knots - 1] spline_coeffs = spline_coeff(knots, y_at_knots);
}

parameters {
  vector[num_knots] y_at_knots;  // Values at knot positions
  // ... other parameters ...
}

model {
  // Evaluate spline at data points
  vector[n_data] y_hat;
  for (i in 1:n_data) {
    y_hat[i] = spline_eval(x[i], knots, spline_coeffs);
  }
  y ~ normal(y_hat, sigma);
}
```

### Key Implementation Details

#### B-splines
1. **Knot placement**: Place knots at quantiles of your x data
2. **Extended knots**: Repeat boundary knots `spline_degree` times
3. **Basis construction**: Use `build_b_spline` to create basis matrix
4. **Smoothing**: Apply random walk prior with `tau_smooth = prior_scale / (smoothing_strength * num_basis^sqrt(2))`

#### C-splines
1. **Include directive**: Must have `spline.stan` file in same directory or in include path
2. **Knot coverage**: Ensure knots cover full data range
3. **Functions**: Use `spline_coeff()` once, then `spline_eval()` for each point
4. **Natural boundaries**: Automatic linear extrapolation at endpoints

### Minimal Working Examples

See these files for complete, working implementations:
- **B-splines**: `code/bsplines.stan` - Full implementation with smoothing
- **C-splines**: `code/csplines.stan` - Natural cubic splines with stan-splines library
- **Regional models**: `examples/regional_splines.stan` - Hierarchical splines for multiple regions

### Common Pitfalls

1. **Array syntax**: Use modern Stan array syntax (`array[N] real x` not `real x[N]`)
2. **Knot coverage**: Ensure knots span your data range (especially for C-splines)
3. **Include paths**: Place `spline.stan` where Stan can find it
4. **Boundary behavior**: Both spline types can extrapolate poorly - keep predictions within data range

### Choosing Spline Parameters

**Number of knots**:
- B-splines: Start with `n/2` (where n is number of data points)
- C-splines: Start with `n/4` (they need fewer knots)
- Adjust based on complexity of your function

**Smoothing (B-splines only)**:
- `smoothing_strength = 0`: No smoothing (independent coefficients)
- `smoothing_strength = 0.05-0.1`: Mild smoothing (0.1 is recommended default)
- `smoothing_strength = 0.1-0.2`: Strong smoothing
- The parameter automatically scales with data variance and number of basis functions (num_basis^sqrt(2))

**Prior scale**:
- Set to `2 * sd(y)` for automatic scaling
- Adjusts coefficient priors to match your data scale

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

This project includes a comprehensive test suite with 8 tests to verify the correctness and properties of both spline implementations.

### Test Files

| Test Script | Purpose | Coverage |
|-------------|---------|----------|
| `tests/test_basic_splines.R` | Basic functionality check | Quick verification with visualization |
| `tests/test_edge_cases_minimal_knots.R` | Minimal knot configurations | Tests with 2-3 knots, error handling |
| `tests/test_flexibility_bspline_vs_cspline.R` | Smoothing comparison | B-spline smoothing vs C-spline knots |
| `tests/test_numerical_accuracy.R` | Mathematical properties | Partition of unity, interpolation, monotonicity |
| `tests/test_analytical_solutions.R` | Known function fitting | Polynomials, constants, sine waves |
| `tests/test_diagnostic_recommendations.R` | Diagnostic system | Over/under-smoothing detection |
| `tests/test_various_target_functions.R` | Comprehensive tests | Six different function types |
| `tests/test_regional_splines.R` | Hierarchical models | Multi-region splines with shared priors |

### Running Tests

```r
# Run all tests with summary (recommended)
source("run-tests.R")

# Run individual test suites
source("tests/test_basic_splines.R")            # Fast basic check
source("tests/test_numerical_accuracy.R")       # Mathematical properties
source("tests/test_diagnostic_recommendations.R") # Diagnostic system
source("tests/test_various_target_functions.R")  # Comprehensive tests
source("tests/test_regional_splines.R")         # Hierarchical models
```

All tests now generate visualization plots saved to the `output/` directory.

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
