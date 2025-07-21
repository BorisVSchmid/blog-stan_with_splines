# Natural Cubic Spline Library Source Information

## Origin
The `spline.stan` file is from the stan-splines library:
- **Repository**: https://github.com/segasai/stan-splines/tree/v2.0.0/stan
- **Author**: Sergey Koposov (University of Edinburgh, Institute for Astronomy)
- **Version**: v2.0.0
- **License**: Listed as "Other (Open)" on Zenodo (see below)

## Description
This library provides efficient implementations of natural cubic splines (also known as C-splines) in Stan. Natural cubic splines are piecewise cubic polynomials that are:
- Continuous up to the second derivative at knot points
- Linear beyond the boundary knots (natural boundary conditions)
- Globally smooth across the entire domain

## Functions Provided
- `spline_geths()`: Calculates spacings between knot points
- `spline_getcoeffs()`: Computes spline coefficients using natural boundary conditions
- `spline_findpos()`: Finds position of a point relative to knots
  - **⚠️ Extrapolation enabled**: Points outside knot range are assigned to boundary intervals
  - Points < first knot use first interval
  - Points ≥ last knot use last interval
- `spline_eval()`: Evaluates the spline at given points (including extrapolated values)

## Modifications
The version in this repository includes minor modifications:
- Enhanced boundary case handling in `spline_findpos()` to enable extrapolation
- Removed unreachable error condition after boundary checks
- Added buffer zones for knot placement to improve numerical stability
- Documented extrapolation behavior in comments

## Citation Requirements
When using this natural cubic spline implementation, please cite:

```bibtex
@article{Koposov2019,
  author = {Koposov, S. E. and others},
  title = {RVSpecFit: Radial velocity and stellar atmospheric parameter fitting},
  journal = {Monthly Notices of the Royal Astronomical Society},
  volume = {485},
  number = {4},
  pages = {4726-4742},
  year = {2019},
  doi = {10.1093/mnras/stz613}
}
```

**Zenodo DOI**: https://doi.org/10.5281/zenodo.7193910  
**ADS Link**: https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract

## Usage in This Project
The library is included in:
- `code/csplines.stan`: Minimal C-spline implementation
- `tests/test_csplines.stan`: Comprehensive C-spline tests

To use in your Stan programs:
```stan
functions {
  #include "cspline_library/spline.stan"
}
```

## Key Differences from B-splines
Unlike B-splines which have local support (each basis function is non-zero only over a limited range), natural cubic splines use global basis functions. This means:
- C-splines generally require fewer parameters for smooth fits
- They naturally extrapolate linearly beyond the data range
- They can be more numerically stable for certain applications
- They may be less flexible for highly local features

## License Note
While the repository lists the license as "Other (Open)" on Zenodo, the specific terms are not clearly stated. Users should:
1. Always provide proper attribution as shown above
2. Check the original repository for any updates to licensing terms
3. Contact the author for clarification if needed for commercial use