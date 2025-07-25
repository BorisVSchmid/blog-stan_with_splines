# Stan Splines Examples

This directory contains advanced examples demonstrating hierarchical and regional spline models.

## Files

### hierarchical_regional_splines.R
Comprehensive example of hierarchical B-splines with:
- Multiple regions with shared global trend
- Region-specific deviations
- Variance component estimation
- Model caching for repeated runs
- Advanced visualization including:
  - Individual region fits
  - Global trend extraction
  - Variance decomposition
  - Coefficient shrinkage patterns

Key features:
- 4 regions with different patterns
- Hierarchical structure for borrowing strength
- MCMC diagnostics with 4 parallel chains
- Publication-quality multi-panel plots

### regional_splines.stan
Stan model implementing hierarchical B-splines:
- Global spline coefficients with region-specific deviations
- Proper hierarchical priors
- Variance component parameters
- Generated quantities for prediction and model assessment

Model structure:
```
y[i,j] ~ normal(global_spline[i] + regional_deviation[i,j], sigma)
regional_deviation[i,j] ~ normal(0, tau_alpha)
```

## Running the Example

```r
# Run the hierarchical regional splines example
source("examples/hierarchical_regional_splines.R")
```

This will:
1. Generate synthetic data for 4 regions
2. Fit the hierarchical model (or load from cache)
3. Print MCMC diagnostics
4. Create comprehensive visualizations
5. Save output plots to `output/` directory

## Output Files

The example generates:
- `output/example-hierarchical_model_fit.rds`: Cached model fit
- `output/example-regional_splines_data.png`: Raw data visualization
- `output/example-regional_splines_fitted.png`: Fitted curves by region
- `output/example-regional_splines_global_coefficients.png`: Coefficient patterns
- `output/example-regional_splines_variance_components.png`: Variance decomposition

## Model Details

### Parameters
- `mu_alpha`: Global spline coefficient means
- `tau_alpha`: Between-region standard deviation
- `alpha_raw`: Standardized region deviations
- `sigma`: Observation noise
- `beta_region`: Region-specific intercepts

### Derived Quantities
- `alpha[r,k]`: Actual coefficients for region r, basis k
- `y_global`: Global trend (average across regions)
- `icc`: Intraclass correlation coefficient

### Sampling Configuration
- 4 chains run in parallel
- 3000 warmup iterations
- 5000 sampling iterations
- `adapt_delta = 0.995` for challenging posterior
- `max_treedepth = 12` for complex hierarchies

## Interpretation

The hierarchical model allows:
1. **Borrowing Strength**: Regions with less data borrow information from others
2. **Shrinkage**: Regional deviations shrink toward global mean
3. **Uncertainty Quantification**: Proper accounting for all sources of variation
4. **Variance Decomposition**: Separate global vs regional variation

## Extending the Example

To adapt for your data:

1. **Change number of regions**:
   ```r
   n_regions <- 6  # Your number of regions
   ```

2. **Modify data generation**:
   ```r
   # Replace synthetic data with your real data
   data_list <- list(
     region1 = data.frame(x = x1, y = y1),
     region2 = data.frame(x = x2, y = y2),
     ...
   )
   ```

3. **Adjust spline parameters**:
   ```r
   num_knots <- 12  # More knots for complex patterns
   spline_degree <- 3  # Keep cubic for smoothness
   ```

4. **Customize priors**:
   ```stan
   // In regional_splines.stan
   tau_alpha ~ normal(0, 1);  // Tighter prior for less variation
   ```

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