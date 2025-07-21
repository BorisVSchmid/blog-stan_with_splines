# Attribution and Third-Party Code

This project includes code derived from the following sources:

## B-spline Implementation

- **Source**: [Stan Case Study: Splines in Stan](https://mc-stan.org/learn-stan/case-studies/splines_in_stan.html)
- **GitHub Repository**: [https://github.com/milkha/Splines_in_Stan](https://github.com/milkha/Splines_in_Stan)
- **Author**: Milad Kharratzadeh
- **Date**: October 24, 2017
- **License**: Not explicitly stated
- **Files affected**: `bsplines.stan`, `test_bsplines.stan` (build_b_spline function)

The B-spline implementation follows the algorithm described in the Stan case study. Users should consult the original source for any licensing terms.

### Penalization/Smoothing Approach

Both Milad's original implementation and our implementation use a random walk penalization approach:
- **Original (b_spline_penalized.stan)**: Uses `a[i] = a[i-1] + a_raw[i]*tau`
- **Our implementation**: Uses the same approach with `alpha[i] = alpha[i-1] + alpha_raw[i] * tau_smooth`

This creates a smoothing effect by constraining successive B-spline coefficients to be similar, with the smoothing strength controlled by the tau parameter.

## Natural Cubic Spline Implementation

- **Source**: [stan-splines library](https://github.com/segasai/stan-splines)
- **Author**: Sergey Koposov (University of Edinburgh, Institute for Astronomy)
- **License**: Listed as "Other (Open)" on Zenodo with no further description
- **License Status**: ⚠️ **Unclear** - No LICENSE file in GitHub repository as of January 2025
- **Files affected**: `code/cspline_library/spline.stan` (entire file), `csplines.stan` and `test_csplines.stan` (include spline.stan)
- **Modifications**: Added boundary case handling in spline_findpos function

### Required Citation

When using the natural cubic spline implementation, please cite:

```
Koposov et al. (2019), MNRAS, 485, 4726
https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract
DOI: https://doi.org/10.5281/zenodo.7193910
```

## License Compatibility Note

### B-splines
The Stan case study does not specify a license. As it's part of Stan's educational materials, it's likely intended for public use, but no explicit license is stated.

### C-splines (stan-splines)
**⚠️ LICENSE UNCLEAR**: 
- Zenodo lists as "Other (Open)" without specification
- No LICENSE file in GitHub repository (checked January 2025)
- Author contact: Sergey Koposov (skoposov AT ed DOT ac DOT uk)

**Recommendations for C-spline usage:**
1. For academic/research use: Cite as required (see above)
2. For commercial use: Contact the author for explicit permission
3. For distribution: Consider implementing your own natural cubic spline functions

This project's BSD 3-Clause License applies only to original code, not to third-party implementations.

## Original Work

All other code in this repository (test scripts, documentation, testing framework) is original work licensed under the BSD 3-Clause License as stated in the LICENSE file.