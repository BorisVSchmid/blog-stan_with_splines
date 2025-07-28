# Implementation Differences from Original Sources

This document details the differences between our spline implementations and the original sources they were derived from.

## B-Splines (compared to Milad Kharratzadeh's Stan Case Study)

Original source: https://github.com/milkha/Splines_in_Stan

### Major Additions and Enhancements

1. **Scale-Invariant Smoothing with Inverted Parameter (Major Addition)**
   - **Original**: Used fixed `tau` parameter for random walk: `a[i] = a[i-1] + a_raw[i]*tau`
     - Larger `tau` = more flexibility (less smooth)
   - **Our Implementation**: Inverted the relationship and made it scale-invariant:
     ```stan
     tau_smooth = prior_scale / sqrt(smoothing_strength * num_basis);
     ```
   - **Key innovation**: Higher `smoothing_strength` = MORE smoothing (smoother curves)
     - `smoothing_strength = 0`: No smoothing (independent coefficients)
     - `smoothing_strength = 0.05-0.1`: Mild smoothing (0.1 is default)
     - `smoothing_strength = 0.1-0.2`: Strong smoothing
   - This ensures consistent smoothing behavior regardless of data scale or number of knots
   - The inverted relationship is more intuitive: larger values mean more smoothing

2. **Boundary Handling for Rightmost Data Point**
   - **Original**: Standard B-spline definition could miss the rightmost data point
   - **Our Implementation**: Added explicit handling in `build_b_spline`:
     ```stan
     if (t[i] >= ext_knots[size(ext_knots) - degree] && 
         ind == size(ext_knots) - degree - 1) {
       b_spline[i] = 1;
     }
     ```

3. **Zero Division Protection**
   - **Original**: Direct division without checking for zero denominators
   - **Our Implementation**: Added epsilon checks:
     ```stan
     if (abs(ext_knots[ind] - ext_knots[ind+order-1]) > 1e-10)
     ```

### Data and Parameter Structure Changes

4. **Data Block**
   - **Original**: Required pre-computed knots and used old array syntax
   - **Our Implementation**: 
     - Automatically computes knots from data
     - Added `smoothing_strength` parameter for user control
     - Added input validation with informative error messages
     - Uses modern array syntax: `array[n_data] real x` instead of `real X[num_data]`

5. **Transformed Data Block**
   - **Original**: Minimal setup
   - **Our Implementation**: 
     - Comprehensive knot placement at quantiles
     - Automatic calculation of extended knots
     - Transform smoothing_strength to tau_smooth with scaling

### Prior and Model Improvements

6. **Priors**
   - **Original**: Generic normal priors
   - **Our Implementation**: 
     - Data-adaptive intercept prior: `normal(mean(y), fmax(sd(y), 0.1 * fmax(abs(mean(y)), 1.0)))`
     - Exponential prior for sigma: `exponential(2)` (better for positive parameters)
     - Standard normal for raw coefficients in non-centered parameterization

7. **Parameterization**
   - **Original**: Centered parameterization
   - **Our Implementation**: Non-centered parameterization for better sampling efficiency

### Generated Quantities

8. **Plotting Support**
   - **Original**: No built-in plotting support
   - **Our Implementation**: Generates 1000-point smooth grid for visualization

### Documentation

9. **Comments and Documentation**
   - **Original**: Minimal comments
   - **Our Implementation**: Extensive documentation explaining scale-invariant smoothing

## C-Splines (compared to Sergey Koposov's stan-splines)

Original source: https://github.com/segasai/stan-splines

### Modifications to spline.stan

1. **Boundary Case Handling in spline_findpos**
   - **Original**: Would fail with error if point outside knot range
   - **Our Implementation**: Added graceful handling with extrapolation to nearest interval:
     ```stan
     if (success==0) {
       // Handle boundary cases - assign to nearest knot interval
       if (x[i] < nodes[1]) {
         ret[i] = 1;  // Extrapolate using first interval
       } else {
         ret[i] = n_nodes - 1;  // Extrapolate using last interval
       }
     }
     ```

### Wrapper Implementation (csplines.stan)

2. **Automatic Knot Placement**
   - **Original**: User must provide knot locations
   - **Our Implementation**: Automatically places knots at quantiles with buffer:
     ```stan
     real buffer = 0.001 * (x_max - x_min);
     knot_locations[1] = x_min - buffer;
     knot_locations[num_knots] = x_max + buffer;
     ```

3. **Input Validation**
   - **Original**: No validation
   - **Our Implementation**: Comprehensive validation with informative error messages

4. **Data Types**
   - **Original**: Uses both array and vector types
   - **Our Implementation**: Consistent use of vectors for spline functions

5. **Prior Specification**
   - **Original**: Not specified (library only)
   - **Our Implementation**: 
     - Normal priors for knot values: `y_at_knots ~ normal(0, 10)`
     - Exponential prior for sigma: `sigma ~ exponential(2)`

6. **Generated Quantities**
   - **Original**: Not provided
   - **Our Implementation**: 
     - 1000-point plotting grid
     - Boundary checking to prevent extrapolation beyond knot range
     - Diagnostic output of spline coefficients

7. **Modern Stan Syntax**
   - Uses modern array declaration syntax
   - Fixed-size arrays in generated quantities (Stan requirement)

## Summary of Key Differences

### Philosophy Changes
1. **User-Friendliness**: Both implementations prioritize ease of use with automatic knot placement
2. **Scale-Invariance**: B-spline implementation ensures consistent behavior across data scales
3. **Robustness**: Added input validation and boundary handling to prevent runtime errors

### Technical Improvements
1. **Numerical Stability**: Added epsilon checks for division by zero
2. **Sampling Efficiency**: Non-centered parameterization for B-splines
3. **Flexibility**: Smoothing strength parameter that adapts to data characteristics

### Practical Enhancements
1. **Visualization**: Built-in plotting grid generation
2. **Diagnostics**: Comprehensive diagnostic output
3. **Documentation**: Extensive inline documentation explaining design choices