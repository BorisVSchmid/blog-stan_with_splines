# Attribution and Third-Party Code

This project includes code derived from the following sources:

## B-spline Implementation

- **Source**: [Stan Case Study: Splines in Stan](https://mc-stan.org/learn-stan/case-studies/splines_in_stan.html)
- **Author**: Milad Kharratzadeh
- **Date**: October 24, 2017
- **License**: Not explicitly stated
- **Files affected**: `test_bsplines.stan` (build_b_spline function)

The B-spline implementation follows the algorithm described in the Stan case study. Users should consult the original source for any licensing terms.

## Natural Cubic Spline Implementation

- **Source**: [stan-splines library](https://github.com/segasai/stan-splines)
- **Author**: Sergey Koposov (University of Edinburgh, Institute for Astronomy)
- **License**: Listed as "Other (Open)" on Zenodo, specific license not found
- **Files affected**: `spline.stan` (entire file), `test_csplines.stan` (includes spline.stan)
- **Modifications**: Added boundary case handling in spline_findpos function

### Required Citation

When using the natural cubic spline implementation, please cite:

```
Koposov et al. (2019), MNRAS, 485, 4726
https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract
DOI: https://doi.org/10.5281/zenodo.7193910
```

## License Compatibility Note

Neither the B-spline case study nor the stan-splines library have clearly stated standard open-source licenses. While this project is distributed under the BSD 3-Clause License, this applies only to the original code written for this project. 

Users should:
1. Check the original sources for their licensing terms
2. Ensure proper attribution when using code derived from these sources
3. Contact the original authors if clarification is needed for commercial use

## Original Work

All other code in this repository (test scripts, documentation, testing framework) is original work licensed under the BSD 3-Clause License as stated in the LICENSE file.