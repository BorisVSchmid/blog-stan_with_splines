# Stan Splines Examples

This directory contains advanced examples demonstrating hierarchical and regional spline models.

## Files

### hierarchical_regional_splines_adaptive.R
Hierarchical B-splines example with adaptive shrinkage:
- 3 regions with shared global trend and regional deviations
- Adaptive shrinkage that learns from data
- Uses adaptive knot selection: `max(4, min(round(n/2), 40))`
- Different smoothing for global (0.05) vs regional (0.10) patterns
- Comprehensive diagnostics and visualization

### regional_splines_adaptive.stan
Stan model implementing adaptive hierarchical B-splines:
- Global spline coefficients with region-specific deviations
- Adaptive shrinkage: `shrinkage_factor ~ normal(1.0, 0.5)`
- Transforms to: `tau_alpha = tau_alpha_base * prior_scale * shrinkage_factor`
- Variance decomposition in generated quantities
- Self-contained implementation (no dependencies on code/ directory)

## Running the Example

```r
# Run the adaptive hierarchical model
source("examples/hierarchical_regional_splines_adaptive.R")
```

This example will:
1. Generate synthetic data for 3 regions with known patterns
2. Fit the hierarchical model using Stan
3. Print MCMC diagnostics and shrinkage information
4. Create comprehensive visualizations
5. Save output plots to `output/` directory
6. Generate diagnostic plots with residuals and metrics

## Output Files
- `output/example-hierarchical_adaptive_decomposition.png`: Pattern decomposition visualization
- `output/example-hierarchical_adaptive_spline_diagnostics.png`: Diagnostic plots
- `output/example-hierarchical_adaptive_spline_diagnostics.csv`: Diagnostic metrics

## Model Details

### Parameters
- `mu_alpha_raw`: Raw global spline coefficients (transformed via random walk)
- `alpha_deviation_raw`: Raw regional deviations (non-centered parameterization)
- `tau_alpha_base`: Base standard deviations before shrinkage
- `shrinkage_factor`: Learned shrinkage parameter (typically ~1, prior centered at 1.0)
- `sigma`: Observation noise
- `beta_region`: Region-specific intercepts
- `mu_beta`, `sigma_beta`: Hierarchical parameters for regional intercepts

### Key Transform
- `tau_alpha = tau_alpha_base * prior_scale * shrinkage_factor`
- This allows the model to learn appropriate regional variation from the data

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

## Model Features

The adaptive model provides:
- Data-driven learning of regional variation strength
- Separate smoothing control for global vs regional patterns
- Proper uncertainty quantification for all parameters
- Variance decomposition between global and regional components


## Key Parameters to Tune

### Smoothing Strength
- Controls wiggliness (shape) of patterns
- **0**: No smoothing (maximum flexibility)
- **0.05**: Mild smoothing (recommended default)
- **0.1+**: Strong smoothing (very smooth patterns)

### Number of Knots
- More knots = more flexibility
- Default: `max(4, min(round(n/2), 40))`

### Adaptive Shrinkage
- Controls amplitude of regional deviations
- Learned from data with prior centered at 1.0
- Higher values = larger regional deviations allowed
- Lower values = more shrinkage toward global pattern

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