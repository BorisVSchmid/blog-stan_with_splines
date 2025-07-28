# Smoothing diagnostics for B-splines and C-splines
# Provides automatic detection of over/under-smoothing with actionable advice

# NOTE: These diagnostics and recommendations are an EXPERIMENTAL FEATURE
# The suggested parameter adjustments are based on heuristics and may need
# fine-tuning for your specific application. Always validate the results.

# GUIDANCE SECTION FOR CLAUDE CODE:
#
# IMPORTANT: Diagnostic advice should be data-driven, not hard-coded
#
# WRONG - Hard-coded multipliers:
#   "Increase prior_scale (try multiplying by 2-5)"
#   "Set smoothing_strength to 40-100"
#
# RIGHT - Data-driven suggestions:
#   "Increase prior_scale (autocorrelation 0.45 suggests multiplying by 2.4)"
#   "Set smoothing_strength > 0 (noise level suggests strength=20)"
#
# The advice should be calculated from the actual diagnostic metrics,
# not predetermined constants
#.
# END OF GUIDANCE SECTION.

library(dplyr)

# Diagnostic function to analyze smoothing level
diagnose_smoothing <- function(fit, x, y, stan_data, model_type = "bspline") {
  
  # Extract draws and compute diagnostics
  draws <- fit$draws(format = "matrix")
  
  # Get fitted values at data points
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  sigma <- mean(draws[, "sigma"])
  
  # Calculate residuals
  residuals <- y - y_hat
  n <- length(y)
  
  # 1. Effective degrees of freedom (approximate)
  # WARNING: This is a crude approximation that doesn't account for smoothing priors!
  # For regularized splines, true EDF is much lower than the number of parameters.
  # This calculation assumes all parameters contribute equally, which is false
  # when smoothing_strength > 0. Use with caution!
  
  # Simple approximation based on model structure
  # For B-splines: num_basis = num_knots + spline_degree - 1
  # Total parameters = num_basis + 1 (for intercept alpha_0)
  if (model_type == "bspline" && "spline_degree" %in% names(stan_data)) {
    max_edf <- stan_data$num_knots + stan_data$spline_degree - 1 + 1  # +1 for intercept
  } else {
    # For C-splines: parameters = num_knots (y_at_knots)
    max_edf <- stan_data$num_knots
  }
  
  # Approximate EDF using ratio of residual variance to total variance
  # This gives a measure of model complexity relative to saturated model
  residual_var <- sum(residuals^2) / (n - max_edf)
  total_var <- sum((y - mean(y))^2) / (n - 1)
  
  # EDF approximation: interpolates between 1 (just intercept) and max_edf
  # Based on how much variance the model explains
  r_squared <- 1 - sum(residuals^2) / sum((y - mean(y))^2)
  edf <- 1 + (max_edf - 1) * r_squared
  
  # Ensure EDF is within reasonable bounds
  edf <- min(max(edf, 1), max_edf)
  
  # 2. Residual analysis
  residual_sd <- sd(residuals)
  residual_autocor <- cor(residuals[-n], residuals[-1])
  
  # 3. Cross-validation approximation (leave-one-out residuals)
  # Approximate using influence diagonal
  h_approx <- pmax(0.01, pmin(0.99, abs(y_hat - mean(y)) / (abs(y - mean(y)) + 1e-6)))
  loo_residuals <- residuals / (1 - h_approx)
  loo_rmse <- sqrt(mean(loo_residuals^2))
  
  # 4. Smoothness measure - second differences
  if (length(unique(x)) > 3) {
    x_order <- order(x)
    y_hat_ordered <- y_hat[x_order]
    second_diff <- diff(diff(y_hat_ordered))
    smoothness <- sd(second_diff)
  } else {
    smoothness <- NA
  }
  
  # 5. Bias-variance indicators
  # High bias: systematic patterns in residuals
  # High variance: fitted values too wiggly
  runs_test <- sum(diff(sign(residuals)) != 0) / (n - 1)  # Proportion of sign changes
  
  # Calculate number of basis functions for B-splines
  if (model_type == "bspline" && "spline_degree" %in% names(stan_data)) {
    num_basis <- stan_data$num_knots + stan_data$spline_degree - 1
  } else {
    # C-splines use num_knots directly
    num_basis <- stan_data$num_knots
  }
  
  # Diagnose smoothing level
  diagnosis <- list(
    edf = edf,
    sigma_estimate = sigma,
    residual_sd = residual_sd,
    residual_autocor = residual_autocor,
    loo_rmse = loo_rmse,
    smoothness = smoothness,
    runs_proportion = runs_test,
    num_knots = stan_data$num_knots,
    num_basis = num_basis,
    smoothing_strength = ifelse("smoothing_strength" %in% names(stan_data), stan_data$smoothing_strength, NA),
    prior_scale = ifelse("prior_scale" %in% names(stan_data), stan_data$prior_scale, NA),
    has_bspline_params = "smoothing_strength" %in% names(stan_data) && "prior_scale" %in% names(stan_data)
  )
  
  # Determine if over-smoothed or overfitted
  over_smoothed <- FALSE
  overfitted <- FALSE
  warnings <- character()
  suggestions <- character()
  
  # Analyze primary indicators with thresholds
  # Good ranges:
  # - Autocorrelation: < 0.2
  # - Runs test: 0.4-0.6 (roughly half sign changes)
  # - Smoothness: < 1.5 * residual_sd
  # Note: EDF thresholds removed - not meaningful for regularized splines
  
  # EDF thresholds based on proportion of available parameters used
  edf_ratio <- edf / max_edf
  high_edf <- edf_ratio > 0.9  # Using >90% of available parameters
  low_edf <- edf_ratio < 0.3   # Using <30% of available parameters
  high_autocor <- residual_autocor > 0.3
  moderate_autocor <- residual_autocor > 0.2
  few_runs <- runs_test < 0.3
  high_smoothness <- !is.na(smoothness) && smoothness > 2 * residual_sd
  
  # Check if basis functions are too close to n
  basis_ratio <- diagnosis$num_basis / n
  too_many_basis <- basis_ratio > 0.25  # More than 25% of data points
  
  # Calculate severity scores for prioritization
  autocor_severity <- max(0, (residual_autocor - 0.2) / 0.3)
  edf_severity <- if (high_edf) (edf_ratio - 0.8) / 0.2 else if (low_edf) (0.4 - edf_ratio) / 0.4 else 0
  smoothness_severity <- if (!is.na(smoothness)) max(0, (smoothness - residual_sd) / residual_sd) else 0
  
  # Determine smoothing state
  # High autocorrelation with many parameters might indicate model misspecification
  if (high_edf && high_autocor) {
    warnings <- c(warnings, 
      "High autocorrelation despite many parameters suggests possible model misspecification",
      "The function may have features that the current spline configuration cannot capture")
    # Don't give contradictory advice in this case
  } else if (high_autocor || few_runs || low_edf) {
    # Classic over-smoothing (too smooth, not flexible enough)
    over_smoothed <- TRUE
    if (high_autocor) warnings <- c(warnings, "High residual autocorrelation suggests over-smoothing")
    if (few_runs) warnings <- c(warnings, "Few sign changes in residuals suggest systematic bias")
    if (low_edf) warnings <- c(warnings, "Very low effective degrees of freedom suggests over-smoothing")
  } else if (high_edf || high_smoothness) {
    # Classic overfitting (too flexible, not smooth enough)
    overfitted <- TRUE
    if (high_edf) warnings <- c(warnings, "Model may be using too many parameters - consider more smoothing")
    if (high_smoothness) warnings <- c(warnings, "High variability in second differences suggests overfitting")
  }
  
  # Generate suggestions based on diagnostic metrics
  # Note: All parameter changes are made in the R script (stan_data list), not in the Stan file
  
  # Special case: model misspecification
  if (high_edf && high_autocor) {
    suggestions <- c(suggestions,
      "Possible model misspecification:"
    )
    
    if (model_type == "bspline" && diagnosis$has_bspline_params) {
      # Check if we already have many knots relative to data
      if (diagnosis$num_basis >= n * 0.4) {
        # Many basis functions but still poor fit - smoothing is likely the issue
        if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 4) {
          suggestions <- c(suggestions,
            sprintf("  - Try reducing smoothing_strength (current: %.1f) to allow more flexibility", 
                    diagnosis$smoothing_strength),
            "  - Consider smoothing_strength = 2-4 for mild smoothing",
            "  - Check for outliers or structural breaks in the data"
          )
        } else if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 0) {
          suggestions <- c(suggestions,
            sprintf("  - Smoothing is already mild (%.1f), try setting to 0 for no smoothing", 
                    diagnosis$smoothing_strength),
            "  - Check for outliers or structural breaks in the data",
            "  - The data may have features that smooth splines cannot capture"
          )
        } else {
          suggestions <- c(suggestions,
            sprintf("  - With %d basis functions and no smoothing, check data quality", diagnosis$num_basis),
            "  - Look for outliers, measurement errors, or model misspecification",
            "  - The data may have features that splines cannot capture well"
          )
        }
      } else {
        suggestions <- c(suggestions,
          sprintf("  - Try reducing smoothing_strength (current: %.1f) for more flexibility", 
                  ifelse(!is.na(diagnosis$smoothing_strength), diagnosis$smoothing_strength, 0)),
          sprintf("  - OR increase num_knots (current: %d) if function is very complex", 
                  diagnosis$num_knots)
        )
      }
    } else {
      # C-spline only has knot control
      # Check if already have very few knots
      if (diagnosis$num_knots <= 5) {
        suggestions <- c(suggestions,
          sprintf("  - With only %d knots, the model may be too inflexible", diagnosis$num_knots),
          "  - Consider increasing num_knots to capture the data complexity",
          "  - OR consider B-splines which offer smoothing control"
        )
      } else {
        suggestions <- c(suggestions,
          sprintf("  - Try adjusting num_knots (current: %d)", diagnosis$num_knots),
          "  - Fewer knots for smoother fit, more knots for flexibility"
        )
      }
    }
  } else if (over_smoothed) {
    suggestions <- c(suggestions,
      "To reduce over-smoothing:"
    )
    
    if (model_type == "bspline" && diagnosis$has_bspline_params) {
      # B-spline specific recommendations
      
      # Check if we already have many knots
      if (diagnosis$num_basis >= n * 0.4) {
        # Already have lots of basis functions
        if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 0) {
          # Calculate suggested smoothing based on current value
          suggested_smooth <- max(0, diagnosis$smoothing_strength / 2)
          suggestions <- c(suggestions,
            sprintf("  - Reduce smoothing_strength to %.1f (current: %.1f)",
                    suggested_smooth, diagnosis$smoothing_strength),
            "  - OR set smoothing_strength = 0 for no smoothing"
          )
        } else {
          suggestions <- c(suggestions,
            "  - Already using no smoothing with many basis functions",
            "  - Over-smoothing diagnosis may be incorrect",
            "  - Check if residual patterns are due to data issues"
          )
        }
      } else {
        # Can potentially add more knots
        min_knots <- max(4, min(round(n/2), 40))
        
        if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 5) {
          suggested_strength <- max(2, diagnosis$smoothing_strength / 2)
          suggestions <- c(suggestions,
            sprintf("  - Reduce smoothing_strength to %.0f (current: %.1f)",
                    suggested_strength, diagnosis$smoothing_strength)
          )
        } else if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 0) {
          suggestions <- c(suggestions,
            sprintf("  - Reduce smoothing_strength (current: %.1f) or set to 0 for no smoothing",
                    diagnosis$smoothing_strength)
          )
        }
        
        if (diagnosis$num_knots < min_knots) {
          suggestions <- c(suggestions,
            sprintf("  - Increase num_knots to %d (current: %d) for more flexibility", 
                    min_knots, diagnosis$num_knots)
          )
        }
      }
    } else {
      # C-spline or other spline types - only knot-based control
      if (diagnosis$num_knots <= 5) {
        suggestions <- c(suggestions,
          sprintf("  - With only %d knots, increasing knots is the main option", diagnosis$num_knots),
          sprintf("  - Try num_knots = %d for more flexibility", min(10, round(n/4)))
        )
      } else {
        suggestions <- c(suggestions,
          sprintf("  - Add more knots - increase num_knots (current: %d)", diagnosis$num_knots)
        )
      }
      
      if (model_type == "cspline") {
        suggestions <- c(suggestions,
          "  - C-splines have limited smoothing control - consider B-splines for finer control"
        )
      }
      
      suggestions <- c(suggestions,
        "  - Check if your data has sufficient variation to support the fitted curve"
      )
    }
  }
  
  if (overfitted) {
    suggestions <- c(suggestions,
      "To reduce overfitting (in order of priority):"
    )
    
    if (model_type == "bspline" && diagnosis$has_bspline_params) {
      # B-spline specific recommendations
      # Prioritize based on severity
      if (edf_severity > 0.5 || too_many_basis) {
        # Severe overfitting - strong measures needed
        # Calculate reasonable knot range
        min_knots <- max(4, min(round(n/2), 40))
        suggested_knots <- max(min_knots, round(diagnosis$num_knots * 0.75))
        
        if (too_many_basis) {
          # Suggest specific smoothing values based on current setting
          current_smooth <- ifelse(!is.na(diagnosis$smoothing_strength), diagnosis$smoothing_strength, 0)
          if (current_smooth < 2) {
            suggested_smooth <- 4
          } else if (current_smooth < 5) {
            suggested_smooth = current_smooth * 2
          } else {
            suggested_smooth = current_smooth + 5
          }
          
          suggestions <- c(suggestions,
            sprintf("  - Increase smoothing_strength to %.0f (current: %.1f) for smoother fit", 
                    suggested_smooth, current_smooth),
            sprintf("  - OR reduce num_knots to %d (currently %d basis functions for %d data points)", 
                    suggested_knots, diagnosis$num_basis, n)
          )
        } else {
          current_smooth <- ifelse(!is.na(diagnosis$smoothing_strength), diagnosis$smoothing_strength, 0)
          suggested_smooth <- ifelse(current_smooth < 2, 2, current_smooth * 1.5)
          
          suggestions <- c(suggestions,
            sprintf("  - Try smoothing_strength = %.0f (current: %.1f) for smoother fit", 
                    suggested_smooth, current_smooth),
            sprintf("  - OR reduce num_knots to %d (current: %d)", 
                    suggested_knots, diagnosis$num_knots)
          )
        }
        
        # Remove prior_scale recommendations entirely
      } else {
        # Mild overfitting - gentler approach
        min_knots <- max(4, min(round(n/2), 40))
        suggested_knots <- max(min_knots, round(diagnosis$num_knots * 0.85))
        
        suggestions <- c(suggestions,
          sprintf("  - Try increasing smoothing_strength slightly (current: %.1f)", 
                  ifelse(!is.na(diagnosis$smoothing_strength), diagnosis$smoothing_strength, 0)),
          sprintf("  - OR reduce num_knots to %d", suggested_knots)
        )
      }
    } else {
      # C-spline or other spline types - only knot-based control
      if (diagnosis$num_knots <= 5) {
        # Very few knots - can't reduce further
        suggestions <- c(suggestions,
          sprintf("  - With only %d knots, overfitting is unlikely", diagnosis$num_knots),
          "  - The wiggliness may be inherent to natural cubic splines",
          "  - Consider B-splines for better smoothing control"
        )
      } else if (edf_severity > 0.5 || too_many_basis) {
        # Severe overfitting with many knots
        min_knots <- max(4, min(round(n/2), 40))  # Use consistent formula
        suggested_knots <- max(min_knots, round(diagnosis$num_knots * 0.7))
        
        suggestions <- c(suggestions,
          sprintf("  - Try reducing num_knots to %d (current: %d)", 
                  suggested_knots, diagnosis$num_knots),
          sprintf("  - %d basis functions for %d data points may be excessive", 
                  diagnosis$num_basis, n)
        )
      } else {
        # Mild overfitting
        min_knots <- max(4, min(round(n/2), 40))
        suggested_knots <- max(min_knots, round(diagnosis$num_knots * 0.85))
        
        suggestions <- c(suggestions,
          sprintf("  - Consider reducing num_knots to %d (current: %d)", 
                  suggested_knots, diagnosis$num_knots)
        )
      }
      
      if (model_type == "cspline") {
        suggestions <- c(suggestions,
          "  Note: C-splines only support knot-based smoothing control"
        )
      }
    }
  }
  
  if (!over_smoothed && !overfitted && !high_edf && !high_autocor) {
    # Check if in good range
    if (residual_autocor < 0.2 && runs_test > 0.3 && runs_test < 0.7) {
      suggestions <- c(suggestions, 
        "Smoothing level appears appropriate:",
        sprintf("  - Residual autocorrelation %.3f is low (< 0.2)", residual_autocor),
        sprintf("  - Runs test %.3f indicates good randomness", runs_test),
        "  - Visual inspection confirms reasonable fit"
      )
    } else {
      # Borderline case - give gentle suggestions
      suggestions <- c(suggestions, "Smoothing level is acceptable but could be fine-tuned:")
      
      if (moderate_autocor && model_type == "bspline" && !is.na(diagnosis$smoothing_strength)) {
        suggestions <- c(suggestions,
          sprintf("  - Slight autocorrelation (%.3f) - consider minor increase in smoothing_strength", residual_autocor)
        )
      }
      
      if (edf_ratio > 0.85) {
        suggestions <- c(suggestions,
          sprintf("  - EDF %.1f of max %.0f (%.0f%%) is moderately high - consider slightly fewer knots", 
                  edf, max_edf, edf_ratio * 100)
        )
      }
    }
  }
  
  diagnosis$over_smoothed <- over_smoothed
  diagnosis$overfitted <- overfitted
  diagnosis$warnings <- warnings
  diagnosis$suggestions <- suggestions
  
  return(diagnosis)
}

# Print diagnostic results
print_smoothing_diagnostics <- function(diagnosis) {
  # Model Summary first
  cat("\nModel Summary:\n")
  cat("==============\n")
  cat(sprintf("Number of knots: %d\n", diagnosis$num_knots))
  cat(sprintf("Basis functions: %d", diagnosis$num_basis))
  if (!is.na(diagnosis$num_basis) && !is.na(diagnosis$num_knots)) {
    cat(" (knots + degree - 1)\n")
  } else {
    cat("\n")
  }
  
  # Add total parameters line for clarity
  if (!is.na(diagnosis$num_basis) && diagnosis$has_bspline_params) {
    cat(sprintf("Total parameters: %d (basis + intercept)\n", diagnosis$num_basis + 1))
  } else if (!is.na(diagnosis$num_knots)) {
    cat(sprintf("Total parameters: %d (knot values)\n", diagnosis$num_knots))
  }
  
  if (!is.na(diagnosis$smoothing_strength)) {
    cat(sprintf("Smoothing strength: %.1f (2.0 is default)\n", diagnosis$smoothing_strength))
  }
  cat(sprintf("Estimated noise σ: %.3f\n", diagnosis$sigma_estimate))
  
  # Smoothing Diagnostics
  cat("\nSmoothing Diagnostics (Experimental):\n")
  cat("====================================\n")
  
  # Add help option
  if (!is.null(diagnosis$show_help) && diagnosis$show_help) {
    cat("\nDiagnostic Metrics Explained:\n")
    cat("- Residual autocorrelation: Correlation between adjacent residuals. High values suggest over-smoothing\n")
    cat("- Residual SD: Standard deviation of residuals. Should be close to true noise level\n")
    cat("- LOO RMSE: Leave-one-out cross-validation error. Lower is better\n")
    cat("- Smoothness (2nd diff SD): Variability in curvature. High values suggest wiggliness\n")
    cat("- Runs test: Proportion of sign changes in residuals. ~0.5 is ideal (random pattern)\n\n")
  }
  
  # Primary diagnostics first
  # Autocorrelation with interpretation
  autocor_level <- ifelse(diagnosis$residual_autocor > 0.3, "high", 
                          ifelse(diagnosis$residual_autocor > 0.2, "moderate", "low"))
  cat(sprintf("Residual autocorrelation:      %5.3f (%s)      → Pattern in residuals (high = over-smoothed)\n", 
              diagnosis$residual_autocor, autocor_level))
  
  runs_status <- ifelse(abs(diagnosis$runs_proportion - 0.5) < 0.1, "good",
                       ifelse(diagnosis$runs_proportion < 0.3, "poor",
                              ifelse(diagnosis$runs_proportion > 0.7, "poor", "fair")))
  cat(sprintf("Runs test:                     %5.3f (%s)      → Randomness of residuals (0.5 = ideal)\n", 
              diagnosis$runs_proportion, runs_status))
  
  cat(sprintf("Residual SD:                   %5.3f             → Spread of residuals (compare to σ=%.3f)\n", 
              diagnosis$residual_sd, diagnosis$sigma_estimate))
  
  if (!is.na(diagnosis$smoothness)) {
    cat(sprintf("Smoothness (2nd diff SD):      %5.3f             → Curvature variation (high = wiggly)\n", 
                diagnosis$smoothness))
  }
  
  cat(sprintf("LOO RMSE:                      %5.3f            → Cross-validation error (lower = better)\n", 
              diagnosis$loo_rmse))
  
  # EDF at the end with clear caveat
  if (!is.na(diagnosis$num_basis) && diagnosis$has_bspline_params) {
    max_edf <- diagnosis$num_basis + 1  # B-spline with intercept
  } else {
    max_edf <- diagnosis$num_knots  # C-spline (y_at_knots includes intercept)
  }
  
  cat("\nModel complexity (use with caution):\n")
  cat(sprintf("  Parameters: %d basis functions + 1 intercept = %d total\n", 
              diagnosis$num_basis, diagnosis$num_basis + 1))
  if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 0) {
    cat(sprintf("  Approximate EDF: %.1f (unreliable for regularized splines)\n", diagnosis$edf))
  } else {
    cat(sprintf("  Approximate EDF: %.1f (no smoothing applied)\n", diagnosis$edf))
  }
  
  # Assessment section
  cat("\nAssessment: ")
  if (diagnosis$over_smoothed) {
    cat("Model appears over-smoothed\n")
  } else if (diagnosis$overfitted) {
    cat("Model appears overfitted\n")
  } else if (length(diagnosis$warnings) > 0 && grepl("misspecification", diagnosis$warnings[1])) {
    cat("Model shows signs of misspecification\n")
  } else {
    cat("Smoothing level appears appropriate\n")
  }
  
  if (length(diagnosis$warnings) > 0) {
    for (warning in diagnosis$warnings) {
      cat("-", warning, "\n")
    }
  }
  
  cat("\nRecommendations:\n")
  if (length(diagnosis$suggestions) > 0) {
    for (suggestion in diagnosis$suggestions) {
      cat(" ", suggestion, "\n")
    }
    # Add general advice about parameter adjustment
    cat("\nParameter adjustment guide:\n")
    cat("- smoothing_strength: 0=none, 1-2=mild, 5-10=strong smoothing (scales with number of basis functions)\n")
    cat("- num_knots: Keep high enough to capture function complexity, adjust smoothing for regularization\n")
    
    cat("\nInterpreting the metrics:\n")
    cat("- Good fit: Autocorr < 0.2, runs test ~0.5, reasonable visual appearance\n")
    cat("- Overfitting: High smoothness measure, very low residual SD, wiggly appearance\n")
    cat("- Over-smoothing: High autocorr (>0.3), runs test far from 0.5, systematic residual patterns\n")
  } else {
    cat("No adjustments needed - model fit appears appropriate.\n")
  }
}

# Wrapper function to add diagnostics to model fitting
fit_with_diagnostics <- function(model, stan_data, x, y, model_type = "bspline", ...) {
  
  cat("Fitting model...\n")
  
  # Fit the model
  fit <- model$sample(
    data = stan_data,
    ...
  )
  
  # Run diagnostics
  cat("\nRunning smoothing diagnostics...\n")
  diagnosis <- diagnose_smoothing(fit, x, y, stan_data, model_type)
  
  # Print diagnostics
  print_smoothing_diagnostics(diagnosis)
  
  # Return both fit and diagnostics
  return(list(
    fit = fit,
    diagnostics = diagnosis
  ))
}

# Example usage function
example_diagnostics <- function() {
  cat("Example: Detecting over-smoothing\n")
  cat("================================\n")
  
  # Generate data
  set.seed(123)
  n <- 50
  x <- seq(0, 10, length.out = n)
  y_true <- sin(x) + 0.4 * cos(3*x)
  y <- y_true + rnorm(n, 0, 0.15)
  
  # Fit with too much smoothing
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 4,      # Few knots
    spline_degree = 3,
    smoothing_strength = 10  # Strong smoothing
  )
  
  cat("\nFitting with restrictive settings (likely to over-smooth)...\n")
  
  # This would run if model is available:
  # model <- cmdstan_model("code/bsplines.stan")
  # result <- fit_with_diagnostics(model, stan_data, x, y, "bspline",
  #                               chains = 2, iter_warmup = 500, 
  #                               iter_sampling = 1000, refresh = 0)
  
  cat("\nThe diagnostics would detect over-smoothing and suggest:\n")
  cat("  - Reducing smoothing_strength\n")
  cat("  - Adding more knots\n")
}

# Function to plot diagnostic information
plot_diagnostic_residuals <- function(fit, x, y) {
  require(ggplot2)
  require(patchwork)
  
  draws <- fit$draws(format = "matrix")
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  residuals <- y - y_hat
  
  df <- data.frame(
    x = x,
    y = y,
    y_hat = y_hat,
    residuals = residuals
  )
  
  # Residuals vs x
  p1 <- ggplot(df, aes(x = x, y = residuals)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    labs(title = "Residuals vs X", 
         subtitle = "Systematic patterns indicate over-smoothing",
         x = "X", y = "Residuals") +
    theme_bw()
  
  # Residuals vs fitted
  p2 <- ggplot(df, aes(x = y_hat, y = residuals)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Residuals vs Fitted", 
         subtitle = "Check for heteroscedasticity",
         x = "Fitted values", y = "Residuals") +
    theme_bw()
  
  # ACF plot approximation
  acf_data <- data.frame(
    lag = 1:(min(20, length(residuals)-1)),
    acf = sapply(1:(min(20, length(residuals)-1)), 
                 function(k) cor(residuals[1:(length(residuals)-k)], 
                               residuals[(k+1):length(residuals)]))
  )
  
  p3 <- ggplot(acf_data, aes(x = lag, y = acf)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = c(-2/sqrt(length(residuals)), 2/sqrt(length(residuals))), 
               linetype = "dashed", color = "red") +
    labs(title = "Residual Autocorrelation",
         subtitle = "High autocorrelation suggests over-smoothing",
         x = "Lag", y = "ACF") +
    theme_bw()
  
  # Combine plots
  combined <- p1 + p2 + p3 + plot_layout(ncol = 2)
  
  return(combined)
}
