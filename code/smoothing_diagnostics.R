# Smoothing diagnostics for B-splines and C-splines
# Provides automatic detection of over/under-smoothing with actionable advice

# NOTE: These diagnostics and recommendations are an EXPERIMENTAL FEATURE
# The suggested parameter adjustments are based on heuristics and may need
# fine-tuning for your specific application. Always validate the results.

# IMPORTANT: Diagnostic advice should be data-driven, not hard-coded
#
# WRONG - Hard-coded multipliers:
#   "Increase prior_scale (try multiplying by 2-5)"
#   "Set smoothing_strength to 10-100"
#
# RIGHT - Data-driven suggestions:
#   "Increase prior_scale (autocorrelation 0.45 suggests multiplying by 2.4)"
#   "Set smoothing_strength > 0 (noise level suggests strength=20)"
#
# The advice should be calculated from the actual diagnostic metrics,
# not predetermined constants.

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
  # For splines, EDF approximates the number of parameters used
  # Maximum EDF = number of basis functions + intercept
  
  # Simple approximation based on model structure
  if (model_type == "bspline" && "spline_degree" %in% names(stan_data)) {
    max_edf <- stan_data$num_knots + stan_data$spline_degree - 1 + 1  # +1 for intercept
  } else {
    max_edf <- stan_data$num_knots + 1  # +1 for intercept  
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
  # - EDF: between 40-80% of max possible EDF
  # - Autocorrelation: < 0.2
  # - Runs test: 0.4-0.6 (roughly half sign changes)
  # - Smoothness: < 1.5 * residual_sd
  
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
  # High EDF with high autocorrelation might indicate model misspecification
  if (high_edf && high_autocor) {
    warnings <- c(warnings, 
      "High EDF with high autocorrelation suggests possible model misspecification",
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
    if (high_edf) warnings <- c(warnings, sprintf("Very high effective degrees of freedom (%.1f of %d possible) suggests overfitting", edf, max_edf))
    if (high_smoothness) warnings <- c(warnings, "High variability in second differences suggests overfitting")
  }
  
  # Generate suggestions based on diagnostic metrics
  # Note: All parameter changes are made in the R script (stan_data list), not in the Stan file
  
  # Special case: model misspecification
  if (high_edf && high_autocor) {
    suggestions <- c(suggestions,
      "Possible solutions for model misspecification:"
    )
    
    if (model_type == "bspline" && diagnosis$has_bspline_params) {
      suggestions <- c(suggestions,
        sprintf("  - Try increasing smoothing_strength (current: %.1f) for smoother fit", 
                ifelse(!is.na(diagnosis$smoothing_strength), diagnosis$smoothing_strength, 0)),
        sprintf("  - OR increase num_knots (current: %d) to capture more complexity", 
                diagnosis$num_knots)
      )
    } else {
      # C-spline only has knot control
      suggestions <- c(suggestions,
        sprintf("  - Try reducing num_knots (current: %d) for smoother fit", 
                diagnosis$num_knots),
        "  - OR increase num_knots to capture more complexity"
      )
    }
  } else if (over_smoothed) {
    suggestions <- c(suggestions,
      "To reduce smoothing (in order of priority):"
    )
    
    if (model_type == "bspline" && diagnosis$has_bspline_params) {
      # B-spline specific recommendations - use smoothing_strength control
      if (autocor_severity > 0.5 || low_edf) {
        # Strong over-smoothing - reduce smoothing strength
        if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 1) {
          suggested_strength <- max(0.1, diagnosis$smoothing_strength / 2)
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce smoothing_strength to %.1f (current: %.1f)",
                    suggested_strength, diagnosis$smoothing_strength)
          )
        }
        if (edf < n/4) {
          suggestions <- c(suggestions,
            sprintf("  2. In R script: Add more knots - increase num_knots (current EDF %.1f is low)", edf)
          )
        }
      } else {
        # Mild over-smoothing - gentle adjustment
        if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength > 0.5) {
          suggested_strength <- max(0.5, diagnosis$smoothing_strength * 0.7)
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce smoothing_strength to %.1f", suggested_strength)
          )
        }
      }
      
      if (!is.na(diagnosis$smoothing_strength) && diagnosis$smoothing_strength < 10) {
        suggestions <- c(suggestions,
          sprintf("  - In R script: Increase smoothing_strength (current %.1f) for more smoothing", 
                  diagnosis$smoothing_strength))
      }
    } else {
      # C-spline or other spline types - only knot-based control
      suggestions <- c(suggestions,
        sprintf("  1. In R script: Add more knots - increase num_knots (current %d with EDF %.1f needs more flexibility)", 
                diagnosis$num_knots, edf)
      )
      
      if (model_type == "cspline") {
        suggestions <- c(suggestions,
          "  2. C-splines have limited smoothing control - consider B-splines for more options"
        )
      }
      
      suggestions <- c(suggestions,
        "  3. Check if your data has sufficient variation to support the fitted curve"
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
        if (too_many_basis) {
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce num_knots significantly (%d basis functions for %d data points is excessive)", 
                    diagnosis$num_basis, n),
            sprintf("     Try %d knots → %d basis functions", 
                    round(n/10), round(n/10) + stan_data$spline_degree - 1)
          )
        } else {
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce num_knots (current %d with EDF %.1f is too flexible)", 
                    diagnosis$num_knots, edf)
          )
        }
        
        if (is.na(diagnosis$smoothing_strength) || diagnosis$smoothing_strength < 1) {
          noise_ratio = residual_sd / sd(y)
          suggested_strength = max(1, min(100, 10 / noise_ratio))
          suggestions <- c(suggestions,
            sprintf("  2. In R script: Set smoothing_strength = %.0f for smoother transitions", suggested_strength)
          )
        }
        
        if (!is.na(diagnosis$prior_scale) && !is.na(smoothness) && smoothness > residual_sd) {
          scale_divisor = min(smoothness / residual_sd, 3)
          suggestions <- c(suggestions,
            sprintf("  3. In R script: Decrease prior_scale (divide by %.1f)", scale_divisor)
          )
        }
      } else {
        # Mild overfitting - gentler approach
        suggestions <- c(suggestions,
          sprintf("  1. In R script: Slightly reduce num_knots to %d", 
                  round(diagnosis$num_knots * 0.8))
        )
        
        if (!is.na(diagnosis$prior_scale) && diagnosis$prior_scale > sd(y)) {
          suggestions <- c(suggestions,
            "  2. In R script: Decrease prior_scale to match data scale"
          )
        }
      }
    } else {
      # C-spline or other spline types - only knot-based control
      if (edf_severity > 0.5 || too_many_basis) {
        # Severe overfitting
        if (too_many_basis && model_type == "bspline") {
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce num_knots (%d basis functions for %d data points is excessive)", 
                    diagnosis$num_basis, n),
            sprintf("     Recommended: %d-%d knots", round(n/12), round(n/8))
          )
        } else {
          suggestions <- c(suggestions,
            sprintf("  1. In R script: Reduce num_knots significantly (current %d with EDF %.1f is overfitting)", 
                    diagnosis$num_knots, edf),
            sprintf("  2. Try reducing to %d knots", max(3, round(diagnosis$num_knots * 0.6)))
          )
        }
      } else {
        # Mild overfitting
        suggestions <- c(suggestions,
          sprintf("  1. In R script: Reduce num_knots to %d (current EDF %.1f is moderately high)", 
                  round(diagnosis$num_knots * 0.8), edf)
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
    if (edf_ratio >= 0.4 && edf_ratio <= 0.8 && residual_autocor < 0.2 && runs_test > 0.3 && runs_test < 0.7) {
      suggestions <- c(suggestions, 
        "Smoothing level appears appropriate:",
        sprintf("  - EDF %.1f of max %.0f (%.0f%% utilization)", edf, max_edf, edf_ratio * 100),
        sprintf("  - Residual autocorrelation %.3f is low (< 0.2)", residual_autocor),
        sprintf("  - Runs test %.3f indicates good randomness", runs_test)
      )
    } else {
      # Borderline case - give gentle suggestions
      suggestions <- c(suggestions, "Smoothing level is acceptable but could be fine-tuned:")
      
      if (moderate_autocor) {
        suggestions <- c(suggestions,
          sprintf("  - Slight autocorrelation (%.3f) - consider minor increase in prior_scale", residual_autocor)
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
  if (!is.na(diagnosis$smoothing_strength)) {
    cat(sprintf("Smoothing strength: %.1f\n", diagnosis$smoothing_strength))
  }
  cat(sprintf("Estimated noise σ: %.3f\n", diagnosis$sigma_estimate))
  
  # Smoothing Diagnostics
  cat("\nSmoothing Diagnostics (Experimental):\n")
  cat("====================================\n")
  
  # Calculate max EDF for context
  max_edf <- ifelse(!is.na(diagnosis$num_basis), diagnosis$num_basis + 1, diagnosis$num_knots + 1)
  edf_percent <- (diagnosis$edf / max_edf) * 100
  cat(sprintf("Effective degrees of freedom: %.1f of %d possible (%.0f%%)\n", 
              diagnosis$edf, max_edf, edf_percent))
  # Autocorrelation with interpretation
  autocor_level <- ifelse(diagnosis$residual_autocor > 0.3, " (high)", 
                          ifelse(diagnosis$residual_autocor > 0.2, " (moderate)", " (low)"))
  cat(sprintf("Residual autocorrelation: %.3f%s\n", diagnosis$residual_autocor, autocor_level))
  cat(sprintf("Residual SD: %.3f\n", diagnosis$residual_sd))
  cat(sprintf("LOO RMSE: %.3f\n", diagnosis$loo_rmse))
  
  if (!is.na(diagnosis$smoothness)) {
    cat(sprintf("Smoothness (2nd diff SD): %.3f\n", diagnosis$smoothness))
  }
  
  cat(sprintf("Runs test: %.3f\n", diagnosis$runs_proportion))
  
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
    cat("- smoothing_strength: 0=none, 1=mild, 10=strong smoothing\n")
    cat("- num_knots: Keep high enough to capture function complexity, adjust smoothing for regularization\n")
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
    smoothing_strength = 0,
    prior_scale = 0.5   # Restrictive prior
  )
  
  cat("\nFitting with restrictive settings (likely to over-smooth)...\n")
  
  # This would run if model is available:
  # model <- cmdstan_model("code/bsplines.stan")
  # result <- fit_with_diagnostics(model, stan_data, x, y, "bspline",
  #                               chains = 2, iter_warmup = 500, 
  #                               iter_sampling = 1000, refresh = 0)
  
  cat("\nThe diagnostics would detect over-smoothing and suggest:\n")
  cat("  - Increasing prior_scale\n")
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