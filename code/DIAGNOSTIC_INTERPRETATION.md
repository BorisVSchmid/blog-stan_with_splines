# Spline Diagnostic Interpretation Guide

This guide helps interpret the diagnostic plots and metrics from `smoothing_diagnostics.R`.

## 1. Heteroscedasticity (Non-constant Variance)

**What to look for:** Residuals that fan out, funnel in, or show other patterns of changing spread.

**What it means:**
- The variance of residuals changes across the range of x values
- Your model assumes constant variance but the data violates this assumption
- Standard errors and confidence intervals may be unreliable

**Common patterns:**
- **Increasing variance** (fan out): Often seen with exponential growth, count data, or percentages near boundaries
- **Decreasing variance** (funnel in): Sometimes seen with saturation effects or measurement precision improving at higher values
- **U-shaped variance**: Higher variance at extremes, common with extrapolation uncertainty

**Implications for splines:**
- The spline is fitting the mean correctly, but uncertainty quantification is wrong
- Some regions will be over-smoothed while others under-smoothed
- Credible intervals will be too narrow where variance is high, too wide where it's low

**Solutions:**
1. **Transform the response**: log(y), sqrt(y), or Box-Cox transformation
2. **Use a different error distribution**: e.g., `y ~ lognormal()` or `y ~ gamma()`
3. **Model the variance**: Add a variance function `sigma = sigma_0 * f(x)`
4. **Weighted regression**: If you know the variance pattern, use weights
5. **Robust methods**: Use Student-t errors instead of normal

## 2. Residuals vs X Plot

**What to look for:** Any systematic patterns in how residuals relate to x values.

**Ideal pattern:**
- Random scatter around zero
- No trends, curves, or clusters
- Constant spread (homoscedasticity)

**Problem patterns and interpretations:**

### Curved patterns
- **U-shaped or inverted U**: Spline is too stiff, increase knots or decrease smoothing
- **S-shaped**: May need more knots at inflection points
- **Multiple waves**: Periodic pattern not captured, add more knots

### Linear trends
- **Upward/downward slope**: Boundary effects or missing linear trend
- **Segmented trends**: Knots may be poorly placed

### Clusters or gaps
- **Vertical bands**: Discrete x values or measurement artifacts
- **Horizontal bands**: Discrete y values or rounding
- **Empty regions**: Extrapolation or sparse data regions

### Variance patterns
- **See heteroscedasticity section above**

**Solutions:**
- Increase number of knots for complex patterns
- Decrease smoothing strength for more flexibility
- Check for data issues (outliers, gaps)
- Consider different knot placement strategies

## 3. Residual Autocorrelation

**What to look for:** Correlation between consecutive residuals (assuming data is ordered by x).

**What it means:**
- **Positive autocorrelation**: Residuals tend to have the same sign as their neighbors
  - Indicates under-fitting - the spline is too smooth
  - Missing important features in the data
  
- **Negative autocorrelation**: Residuals alternate signs rapidly
  - Indicates over-fitting - the spline is too wiggly
  - Fitting noise rather than signal

**Interpretation of ACF plot:**
- Blue dashed lines: 95% confidence bounds for white noise
- Bars outside bounds: Significant autocorrelation
- Lag 1 most important for spline fitting

**Common patterns:**
- **Decaying positive ACF**: Smooth underlying process not fully captured
- **Alternating +/- ACF**: Over-fitting with too many knots
- **Periodic ACF**: Missing seasonal/cyclical component

**Solutions for positive autocorrelation:**
1. Add more knots to capture finer features
2. Reduce smoothing strength (smaller smoothing_strength parameter)
3. Check if there's a missing covariate or time trend
4. Consider modeling the correlation structure directly

**Solutions for negative autocorrelation:**
1. Reduce number of knots
2. Increase smoothing strength
3. Check for measurement error or artificial discretization

## Quick Reference Table

| Diagnostic | Good | Problem | Action |
|------------|------|---------|--------|
| Heteroscedasticity | Constant spread | Fan/funnel pattern | Transform y or model variance |
| Residuals vs X | Random scatter | Curves/trends | Adjust knots/smoothing |
| Autocorrelation | Near zero | Significant lag 1 | More knots if positive, fewer if negative |

## Using with smoothing_diagnostics.R

The diagnostic function automatically detects these patterns and provides recommendations:

```r
diagnosis <- diagnose_smoothing(x, y, fitted_values)
print_smoothing_diagnostics(diagnosis)
```

The function will flag:
- Heteroscedasticity via Breusch-Pagan test
- Residual patterns via runs test
- Autocorrelation via Durbin-Watson test

And suggest specific parameter adjustments based on the combination of issues detected.