# Test smoothing diagnostics with both B-splines and C-splines

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)
source("code/smoothing_diagnostics.R")

# Generate test data
set.seed(123)
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.4 * cos(3*x) + 0.2*x
y <- y_true + rnorm(n, 0, 0.15)

cat("Testing Smoothing Diagnostics for Both Spline Types\n")
cat("===================================================\n\n")

# Test 1: B-splines with various settings
cat("Test 1: B-splines with normal smoothing\n")
cat("---------------------------------------\n")

# Prepare B-spline data
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 8,
  spline_degree = 3,
  smoothing_strength = 0,
  prior_scale = 2 * sd(y)
)

# Compile and fit B-spline model
model_b <- cmdstan_model("code/bsplines.stan")
fit_b <- model_b$sample(
  data = stan_data_b,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

# Run diagnostics
diagnosis_b <- diagnose_smoothing(fit_b, x, y, stan_data_b, "bspline")
cat("\nB-spline Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_b)

# Test 2: C-splines
cat("\n\nTest 2: C-splines\n")
cat("-----------------\n")

# Prepare C-spline data (no prior_scale or tau_smooth)
stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 8
)

# Compile and fit C-spline model
model_c <- cmdstan_model("code/csplines.stan")
fit_c <- model_c$sample(
  data = stan_data_c,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

# Run diagnostics
diagnosis_c <- diagnose_smoothing(fit_c, x, y, stan_data_c, "cspline")
cat("\nC-spline Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_c)

# Test 3: B-splines with overfitting (many knots)
cat("\n\nTest 3: B-splines with potential overfitting (many knots)\n")
cat("----------------------------------------------------------\n")

stan_data_b_overfit <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 20,  # Many knots for n=40
  spline_degree = 3,
  smoothing_strength = 0,
  prior_scale = 5 * sd(y)  # Large prior scale
)

fit_b_overfit <- model_b$sample(
  data = stan_data_b_overfit,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

diagnosis_b_overfit <- diagnose_smoothing(fit_b_overfit, x, y, stan_data_b_overfit, "bspline")
cat("\nB-spline (overfitting) Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_b_overfit)

# Test 4: C-splines with few knots (over-smoothing)
cat("\n\nTest 4: C-splines with potential over-smoothing (few knots)\n")
cat("------------------------------------------------------------\n")

stan_data_c_smooth <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 3  # Very few knots
)

fit_c_smooth <- model_c$sample(
  data = stan_data_c_smooth,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)

diagnosis_c_smooth <- diagnose_smoothing(fit_c_smooth, x, y, stan_data_c_smooth, "cspline")
cat("\nC-spline (over-smoothing) Diagnostics:\n")
print_smoothing_diagnostics(diagnosis_c_smooth)

# Create visualization plots
cat("\n\nCreating diagnostic visualization plots...\n")

library(patchwork)

# Function to create plot for each scenario
create_diagnostic_plot <- function(fit, x, y, y_true, title, subtitle, color = "blue") {
  draws <- fit$draws(format = "matrix")
  
  # Get fitted values for plot
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  # Get fitted values at data points for residuals
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  residuals <- y - y_hat
  
  # Main fit plot
  df_fit <- data.frame(
    x = x_plot,
    y = y_plot,
    ymin = y_plot_lower,
    ymax = y_plot_upper
  )
  
  p_fit <- ggplot() +
    geom_ribbon(data = df_fit, aes(x = x, ymin = ymin, ymax = ymax), 
                alpha = 0.3, fill = color) +
    geom_line(data = df_fit, aes(x = x, y = y), color = color, linewidth = 1.2) +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
               alpha = 0.6, size = 1.5) +
    geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y), 
              color = "red", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    labs(title = title, subtitle = subtitle, x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  # Residual plot
  df_resid <- data.frame(x = x, residuals = residuals)
  
  p_resid <- ggplot(df_resid, aes(x = x, y = residuals)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.6, color = color) +
    geom_smooth(method = "loess", se = TRUE, color = "black", 
                fill = "gray80", alpha = 0.3, size = 0.5) +
    labs(title = "Residuals", x = "x", y = "Residuals") +
    theme_bw() +
    theme(plot.title = element_text(size = 9))
  
  return(list(fit = p_fit, resid = p_resid))
}

# Create plots for all 4 scenarios
plots_b_normal <- create_diagnostic_plot(fit_b, x, y, y_true, 
                                         "B-spline: Normal Settings",
                                         sprintf("8 knots, smoothing=0, autocorr=%.3f", 
                                                 diagnosis_b$residual_autocor),
                                         "blue")

plots_c_normal <- create_diagnostic_plot(fit_c, x, y, y_true,
                                         "C-spline: Normal Settings", 
                                         sprintf("8 knots, autocorr=%.3f",
                                                 diagnosis_c$residual_autocor),
                                         "darkgreen")

plots_b_overfit <- create_diagnostic_plot(fit_b_overfit, x, y, y_true,
                                          "B-spline: Overfitting",
                                          sprintf("20 knots, smoothing=0, autocorr=%.3f",
                                                  diagnosis_b_overfit$residual_autocor),
                                          "darkred")

plots_c_smooth <- create_diagnostic_plot(fit_c_smooth, x, y, y_true,
                                         "C-spline: Over-smoothing",
                                         sprintf("3 knots, autocorr=%.3f",
                                                 diagnosis_c_smooth$residual_autocor),
                                         "darkorange")

# Combine all plots
combined_plot <- (plots_b_normal$fit | plots_c_normal$fit | 
                  plots_b_overfit$fit | plots_c_smooth$fit) /
                 (plots_b_normal$resid | plots_c_normal$resid | 
                  plots_b_overfit$resid | plots_c_smooth$resid)

combined_plot <- combined_plot + 
  plot_annotation(
    title = "Diagnostic Recommendations Test: Four Scenarios",
    subtitle = "Top: Fitted curves (red dashed = true function) | Bottom: Residual patterns",
    caption = "Diagnostics detect over-smoothing (high autocorr) and overfitting (wiggly fit)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-diagnostic_recommendations_scenarios.png", combined_plot, 
       width = 14, height = 8, dpi = 300)
cat("Saved diagnostic scenarios plot to output/test-diagnostic_recommendations_scenarios.png\n")

# Create a second plot showing diagnostic metrics
create_metrics_plot <- function(diagnosis_list, names) {
  # Extract key metrics
  metrics_df <- data.frame(
    Scenario = factor(names, levels = names),
    Autocorrelation = sapply(diagnosis_list, function(d) d$residual_autocor),
    Runs_Test = sapply(diagnosis_list, function(d) d$runs_proportion),
    Residual_SD = sapply(diagnosis_list, function(d) d$residual_sd),
    Smoothness = sapply(diagnosis_list, function(d) ifelse(is.na(d$smoothness), 0, d$smoothness))
  )
  
  # Convert to long format for plotting
  metrics_long <- tidyr::pivot_longer(metrics_df, 
                                      cols = -Scenario,
                                      names_to = "Metric",
                                      values_to = "Value")
  
  # Create the plot
  p <- ggplot(metrics_long, aes(x = Scenario, y = Value, fill = Scenario)) +
    geom_col() +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("blue", "darkgreen", "darkred", "darkorange")) +
    labs(title = "Diagnostic Metrics Comparison",
         subtitle = "Key metrics that drive the diagnostic recommendations",
         x = "", y = "Metric Value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  return(p)
}

# Create metrics comparison plot
diagnosis_list <- list(diagnosis_b, diagnosis_c, diagnosis_b_overfit, diagnosis_c_smooth)
scenario_names <- c("B-spline Normal", "C-spline Normal", 
                    "B-spline Overfit", "C-spline Over-smooth")

metrics_plot <- create_metrics_plot(diagnosis_list, scenario_names)

ggsave("output/test-diagnostic_metrics_comparison.png", metrics_plot,
       width = 10, height = 8, dpi = 300)
cat("Saved diagnostic metrics plot to output/test-diagnostic_metrics_comparison.png\n")

# Summary
cat("\n\nSummary of Diagnostic Differences\n")
cat("=================================\n")
cat("B-splines have access to:\n")
cat("  - prior_scale adjustments\n")
cat("  - tau_smooth for random walk smoothing\n")
cat("  - num_knots control\n\n")
cat("C-splines only have:\n")
cat("  - num_knots control\n")
cat("  - Natural boundary conditions (built-in)\n\n")
cat("The diagnostics now provide appropriate advice for each spline type.\n")