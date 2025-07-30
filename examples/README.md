# Stan Splines Examples

This directory contains advanced examples demonstrating hierarchical and regional spline models.

## Files

### hierarchical_regional_splines_adaptive.R
Hierarchical B-splines example with adaptive shrinkage:
- 3 regions with shared global trend and regional deviations
- Adaptive shrinkage that learns from data
- Doubled knots: `2 * max(4, min(round(n/2), 40))`
- Equal smoothing (0.05) for both global and regional patterns
- Proper residual diagnostics for hierarchical models

### regional_splines_adaptive.stan
Stan model implementing adaptive hierarchical B-splines:
- Global spline coefficients with region-specific deviations
- Adaptive shrinkage: `shrinkage_factor ~ normal(5, 2)`
- Transforms to: `tau_alpha = tau_alpha_base * (prior_scale / shrinkage_factor)`
- Variance decomposition in generated quantities

## Running the Examples

```r
# Run the original model (fixed shrinkage)
source("examples/hierarchical_regional_splines.R")

# Run the adaptive model (recommended)
source("examples/hierarchical_regional_splines_adaptive.R")
```

Both examples will:
1. Generate synthetic data for 3 regions with known patterns
2. Fit the hierarchical model (or load from cache)
3. Print MCMC diagnostics and shrinkage information
4. Create comprehensive visualizations
5. Save output plots to `output/` directory
6. Generate diagnostic plots with correct residuals

## Output Files

### Original Model
- `output/example-hierarchical_model_draws.rds`: Cached model draws
- `output/example-hierarchical_decomposition.png`: Pattern decomposition
- `output/example-hierarchical_spline_diagnostics.png`: Diagnostic plots
- `output/example-hierarchical_spline_diagnostics.csv`: Diagnostic metrics

### Adaptive Model
- `output/example-hierarchical_adaptive_model_draws.rds`: Cached model draws
- `output/example-hierarchical_adaptive_decomposition.png`: Improved decomposition
- `output/example-hierarchical_adaptive_spline_diagnostics.png`: Better diagnostics
- `output/example-hierarchical_adaptive_spline_diagnostics.csv`: Diagnostic metrics

## Model Details

### Parameters (Original Model)
- `mu_alpha`: Global spline coefficient means
- `tau_alpha`: Between-region standard deviation (fixed at `prior_scale/5`)
- `alpha_deviation_raw`: Standardized region deviations
- `sigma`: Observation noise
- `beta_region`: Region-specific intercepts

### Additional Parameters (Adaptive Model)
- `shrinkage_factor`: Learned shrinkage parameter (~5 typical)
- `tau_alpha_base`: Base standard deviations before shrinkage
- Transforms to: `tau_alpha = tau_alpha_base * (prior_scale / shrinkage_factor)`

### Key Differences
- **Fixed model**: Regional variance is `prior_scale/5` (hardcoded)
- **Adaptive model**: Regional variance adapts to data via `shrinkage_factor`

### Sampling Configuration
- 4 chains run in parallel
- 3000 warmup iterations
- 5000 sampling iterations
- `adapt_delta = 0.999` for challenging posterior
- `max_treedepth = 15` for complex hierarchies

## Interpretation

The hierarchical model allows:
1. **Borrowing Strength**: Regions with less data borrow information from others
2. **Shrinkage**: Regional deviations shrink toward global mean
3. **Uncertainty Quantification**: Proper accounting for all sources of variation
4. **Variance Decomposition**: Separate global vs regional variation

## Choosing Between Models

### Use the Adaptive Model When:
- You don't know how much regional variation to expect
- Regional effects may vary across datasets
- You want the model to learn the appropriate balance
- You need better diagnostics (correct residuals)


## Key Parameters to Tune

### Smoothing Strength
- Controls wiggliness (shape) of patterns
- **0**: No smoothing (maximum flexibility)
- **0.05**: Mild smoothing (recommended default)
- **0.1+**: Strong smoothing (very smooth patterns)

### Number of Knots
- More knots = more flexibility
- Default: `max(4, min(round(n/2), 40))`
- Adaptive model doubles this for better flexibility

### Adaptive Shrinkage
- Controls amplitude of regional deviations
- Learned from data, typically 2-10
- Higher values = more shrinkage toward global

## Troubleshooting

Common issues:

1. **Divergent transitions**: Increase `adapt_delta` (up to 0.999)
2. **Max treedepth**: Increase `max_treedepth` (up to 15)
3. **Slow sampling**: Reduce `num_knots` or use stronger priors
4. **Poor mixing**: Check parameterization, consider longer chains

## References

This example demonstrates concepts from:
- Hierarchical modeling (Gelman & Hill, 2007)
- Functional data analysis (Ramsay & Silverman, 2005)
- Regional statistics (Cressie & Wikle, 2011)