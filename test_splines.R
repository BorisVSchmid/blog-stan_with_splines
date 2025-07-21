# Comprehensive test script for B-splines and C-splines in Stan
# Tests both implementations with various functions and edge cases

# Package management
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
# Note: cmdstanr must be installed from Stan repo
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Use groundhog for dependencies
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

# Load cmdstanr after its dependencies
library(cmdstanr)

set.seed(123)

# Create output directory if it doesn't exist
dir.create("output", showWarnings = FALSE)

# Function to generate test data
generate_test_data <- function(n = 50, func_type = "sine", noise_sd = 0.1) {
  x <- seq(0, 10, length.out = n)
  
  y_true <- switch(func_type,
    "sine" = sin(x),
    "polynomial" = 0.1 * x^3 - x^2 + 2*x + 1,
    "step" = ifelse(x < 5, 0, 1),
    "exponential" = exp(-x/3),
    "complex" = sin(x) + 0.5 * cos(3*x),
    "linear" = 2 * x + 1
  )
  
  y <- y_true + rnorm(n, 0, noise_sd)
  
  list(x = x, y = y, y_true = y_true, func_type = func_type)
}

# Function to fit B-spline model
fit_bspline <- function(data, num_knots = 5, spline_degree = 3, 
                        chains = 2, iter_warmup = 500, iter_sampling = 1000) {
  
  stan_data <- list(
    n_data = length(data$x),
    x = data$x,
    y = data$y,
    num_knots = num_knots,
    spline_degree = spline_degree
  )
  
  # Compile model if needed
  model <- cmdstan_model("test_bsplines.stan")
  
  # Fit model
  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  return(fit)
}

# Function to fit C-spline model
fit_cspline <- function(data, num_knots = 5, 
                        chains = 2, iter_warmup = 500, iter_sampling = 1000) {
  
  stan_data <- list(
    n_data = length(data$x),
    x = data$x,
    y = data$y,
    num_knots = num_knots
  )
  
  # Compile model if needed
  model <- cmdstan_model("test_csplines.stan")
  
  # Fit model
  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  return(fit)
}

# Function to extract and plot results
plot_spline_fit <- function(fit, data, spline_type, show_basis = FALSE) {
  draws <- fit$draws(format = "matrix")
  
  # Extract fitted values at plot points
  x_plot_cols <- grep("x_plot\\[", colnames(draws), value = TRUE)
  y_plot_cols <- grep("y_plot\\[", colnames(draws), value = TRUE)
  
  if (length(x_plot_cols) > 0 && length(y_plot_cols) > 0) {
    x_plot <- colMeans(draws[, x_plot_cols])
    y_plot_mean <- colMeans(draws[, y_plot_cols])
    y_plot_lower <- apply(draws[, y_plot_cols], 2, quantile, 0.025)
    y_plot_upper <- apply(draws[, y_plot_cols], 2, quantile, 0.975)
    
    plot_df <- data.frame(
      x = x_plot,
      y_mean = y_plot_mean,
      y_lower = y_plot_lower,
      y_upper = y_plot_upper
    )
    
    p <- ggplot() +
      geom_ribbon(data = plot_df, aes(x = x, ymin = y_lower, ymax = y_upper), 
                  alpha = 0.2, fill = "blue") +
      geom_line(data = plot_df, aes(x = x, y = y_mean), color = "blue", size = 1) +
      geom_point(data = data.frame(x = data$x, y = data$y), 
                 aes(x = x, y = y), alpha = 0.6) +
      geom_line(data = data.frame(x = data$x, y = data$y_true), 
                aes(x = x, y = y), color = "red", linetype = "dashed", alpha = 0.5) +
      labs(title = paste(spline_type, "fit"),
           subtitle = data$func_type,
           x = "x", y = "y") +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 9))
    
    # Add basis functions for B-splines if requested
    if (show_basis && spline_type == "B-spline") {
      basis_cols <- grep("basis_functions\\[", colnames(draws), value = TRUE)
      if (length(basis_cols) > 0) {
        n_basis <- length(grep("basis_functions\\[1,", colnames(draws), value = TRUE))
        basis_df <- NULL
        
        for (i in 1:n_basis) {
          cols <- grep(paste0("basis_functions\\[\\d+,", i, "\\]"), 
                      colnames(draws), value = TRUE)
          if (length(cols) > 0) {
            basis_vals <- colMeans(draws[, cols])
            temp_df <- data.frame(
              x = x_plot,
              y = basis_vals,
              basis = factor(i)
            )
            basis_df <- rbind(basis_df, temp_df)
          }
        }
        
        if (!is.null(basis_df)) {
          p_basis <- ggplot(basis_df, aes(x = x, y = y, color = basis)) +
            geom_line() +
            labs(title = "B-spline basis functions",
                 x = "x", y = "Basis value") +
            theme_bw() +
            theme(legend.position = "none",
                  plot.title = element_text(size = 10))
          
          return(list(fit = p, basis = p_basis))
        }
      }
    }
    
    return(list(fit = p, basis = NULL))
  } else {
    warning("Could not extract plot data from fit")
    return(list(fit = NULL, basis = NULL))
  }
}

# Function to run diagnostics
check_diagnostics <- function(fit, model_name) {
  cat("\n=== Diagnostics for", model_name, "===\n")
  
  # Check for divergences
  sampler_diags <- fit$sampler_diagnostics()
  divergences <- sum(sampler_diags[,,"divergent__"])
  cat("Divergences:", divergences, "\n")
  
  # Check Rhat
  summary_df <- fit$summary()
  high_rhat <- summary_df %>% filter(rhat > 1.01)
  if (nrow(high_rhat) > 0) {
    cat("Parameters with Rhat > 1.01:", nrow(high_rhat), "\n")
  } else {
    cat("All Rhat values < 1.01\n")
  }
  
  # Check ESS
  low_ess <- summary_df %>% filter(ess_bulk < 400 | ess_tail < 400)
  if (nrow(low_ess) > 0) {
    cat("Parameters with low ESS:", nrow(low_ess), "\n")
  } else {
    cat("All ESS values adequate\n")
  }
  
  # Return diagnostics summary
  list(
    divergences = divergences,
    high_rhat = nrow(high_rhat),
    low_ess = nrow(low_ess),
    summary = summary_df
  )
}

# Main testing routine
cat("=== Testing B-splines and C-splines in Stan ===\n\n")

# Store all plots
all_plots <- list()

# Test 1: Simple sine function
cat("Test 1: Sine function\n")
data_sine <- generate_test_data(n = 30, func_type = "sine", noise_sd = 0.1)

cat("  Fitting B-spline model...\n")
fit_b_sine <- fit_bspline(data_sine, num_knots = 7, spline_degree = 3)
diag_b_sine <- check_diagnostics(fit_b_sine, "B-spline (sine)")
plots_b_sine <- plot_spline_fit(fit_b_sine, data_sine, "B-spline", show_basis = TRUE)

cat("\n  Fitting C-spline model...\n")
fit_c_sine <- fit_cspline(data_sine, num_knots = 7)
diag_c_sine <- check_diagnostics(fit_c_sine, "C-spline (sine)")
plots_c_sine <- plot_spline_fit(fit_c_sine, data_sine, "C-spline")

# Test 2: Polynomial function
cat("\n\nTest 2: Polynomial function\n")
data_poly <- generate_test_data(n = 40, func_type = "polynomial", noise_sd = 0.5)

cat("  Fitting B-spline model...\n")
fit_b_poly <- fit_bspline(data_poly, num_knots = 5, spline_degree = 3)
diag_b_poly <- check_diagnostics(fit_b_poly, "B-spline (polynomial)")
plots_b_poly <- plot_spline_fit(fit_b_poly, data_poly, "B-spline")

cat("\n  Fitting C-spline model...\n")
fit_c_poly <- fit_cspline(data_poly, num_knots = 5)
diag_c_poly <- check_diagnostics(fit_c_poly, "C-spline (polynomial)")
plots_c_poly <- plot_spline_fit(fit_c_poly, data_poly, "C-spline")

# Test 3: Complex function
cat("\n\nTest 3: Complex function (sin + cos)\n")
data_complex <- generate_test_data(n = 50, func_type = "complex", noise_sd = 0.1)

cat("  Fitting B-spline model...\n")
fit_b_complex <- fit_bspline(data_complex, num_knots = 10, spline_degree = 3)
diag_b_complex <- check_diagnostics(fit_b_complex, "B-spline (complex)")
plots_b_complex <- plot_spline_fit(fit_b_complex, data_complex, "B-spline")

cat("\n  Fitting C-spline model...\n")
fit_c_complex <- fit_cspline(data_complex, num_knots = 10)
diag_c_complex <- check_diagnostics(fit_c_complex, "C-spline (complex)")
plots_c_complex <- plot_spline_fit(fit_c_complex, data_complex, "C-spline")

# Create combined plots using patchwork
if (!is.null(plots_b_sine$fit) && !is.null(plots_c_sine$fit)) {
  # Main comparison plot
  main_plot <- (plots_b_sine$fit | plots_c_sine$fit) / 
               (plots_b_poly$fit | plots_c_poly$fit) / 
               (plots_b_complex$fit | plots_c_complex$fit) +
    plot_annotation(
      title = "B-splines vs C-splines: Function Fitting Comparison",
      subtitle = "Blue: fitted spline with 95% CI, Red dashed: true function, Points: data"
    )
  
  ggsave("output/spline_comparison.png", main_plot, width = 10, height = 12, dpi = 300)
  cat("\n\nMain comparison plot saved to output/spline_comparison.png\n")
  
  # Basis functions plot
  if (!is.null(plots_b_sine$basis)) {
    basis_plot <- plots_b_sine$basis +
      plot_annotation(title = "B-spline Basis Functions (Sine Example)")
    ggsave("output/bspline_basis_functions.png", basis_plot, width = 8, height = 6, dpi = 300)
    cat("Basis functions plot saved to output/bspline_basis_functions.png\n")
  }
}

# Test edge cases
cat("\n\n=== Testing edge cases ===\n")

# Test 4: Minimal data (n = 5)
cat("\nTest 4: Minimal data points (n=5)\n")
data_minimal <- generate_test_data(n = 5, func_type = "linear", noise_sd = 0.05)

cat("  Fitting B-spline model with 3 knots...\n")
tryCatch({
  fit_b_minimal <- fit_bspline(data_minimal, num_knots = 3, spline_degree = 2)
  diag_b_minimal <- check_diagnostics(fit_b_minimal, "B-spline (minimal)")
  plots_b_minimal <- plot_spline_fit(fit_b_minimal, data_minimal, "B-spline")
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

cat("\n  Fitting C-spline model with 3 knots...\n")
tryCatch({
  fit_c_minimal <- fit_cspline(data_minimal, num_knots = 3)
  diag_c_minimal <- check_diagnostics(fit_c_minimal, "C-spline (minimal)")
  plots_c_minimal <- plot_spline_fit(fit_c_minimal, data_minimal, "C-spline")
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# Test 5: Different spline degrees for B-splines
cat("\n\nTest 5: Different B-spline degrees\n")
data_test <- generate_test_data(n = 30, func_type = "sine", noise_sd = 0.1)

degree_plots <- list()
for (degree in c(2, 3, 4)) {
  cat(paste0("\n  Testing B-spline with degree ", degree, "...\n"))
  tryCatch({
    fit_b_degree <- fit_bspline(data_test, num_knots = 6, spline_degree = degree, 
                                chains = 2, iter_warmup = 300, iter_sampling = 500)
    diag <- check_diagnostics(fit_b_degree, paste0("B-spline (degree ", degree, ")"))
    plots <- plot_spline_fit(fit_b_degree, data_test, paste0("B-spline (degree ", degree, ")"))
    if (!is.null(plots$fit)) {
      degree_plots[[as.character(degree)]] <- plots$fit
    }
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
  })
}

# Create degree comparison plot
if (length(degree_plots) > 0) {
  degree_comparison <- wrap_plots(degree_plots, ncol = 3) +
    plot_annotation(title = "B-spline Degree Comparison")
  ggsave("output/bspline_degree_comparison.png", degree_comparison, width = 12, height = 5, dpi = 300)
  cat("\nDegree comparison plot saved to output/bspline_degree_comparison.png\n")
}

# Save diagnostics summary
diagnostics_summary <- list(
  sine = list(
    bspline = diag_b_sine,
    cspline = diag_c_sine
  ),
  polynomial = list(
    bspline = diag_b_poly,
    cspline = diag_c_poly
  ),
  complex = list(
    bspline = diag_b_complex,
    cspline = diag_c_complex
  ),
  timestamp = Sys.time()
)

saveRDS(diagnostics_summary, "output/diagnostics_summary.rds")

# Performance comparison
cat("\n\n=== Performance comparison ===\n")
cat("Comparing computation time for B-splines vs C-splines\n")

data_perf <- generate_test_data(n = 100, func_type = "complex", noise_sd = 0.1)

time_b <- system.time({
  fit_b_perf <- fit_bspline(data_perf, num_knots = 8, chains = 2, 
                            iter_warmup = 500, iter_sampling = 1000)
})

time_c <- system.time({
  fit_c_perf <- fit_cspline(data_perf, num_knots = 8, chains = 2, 
                            iter_warmup = 500, iter_sampling = 1000)
})

cat("\nB-spline time:", time_b["elapsed"], "seconds\n")
cat("C-spline time:", time_c["elapsed"], "seconds\n")

# Save timing results
timing_results <- data.frame(
  method = c("B-spline", "C-spline"),
  elapsed_time = c(time_b["elapsed"], time_c["elapsed"]),
  n_data = 100,
  n_knots = 8,
  chains = 2,
  total_iterations = 1500
)

write.csv(timing_results, "output/timing_results.csv", row.names = FALSE)

cat("\n\n=== Testing complete ===\n")
cat("Results saved to output/ directory:\n")
cat("  - spline_comparison.png: Main comparison plots\n")
cat("  - bspline_basis_functions.png: B-spline basis visualization\n")
cat("  - bspline_degree_comparison.png: Effect of spline degree\n")
cat("  - diagnostics_summary.rds: Full diagnostics results\n")
cat("  - timing_results.csv: Performance comparison\n")