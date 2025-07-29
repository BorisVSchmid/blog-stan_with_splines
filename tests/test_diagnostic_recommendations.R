# Test smoothing diagnostics with both B-splines and C-splines

# Suppress default graphics device to prevent Rplots.pdf
pdf(NULL)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork")
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

# Prepare B-spline data with default settings
default_bspline_knots <- max(4, min(round(n/2), 40))  # n/2 rule
stan_data_b <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = default_bspline_knots,
  spline_degree = 3,
  smoothing_strength = 0.07,  # Lighter smoothing for normal case
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

# Prepare C-spline data with default settings
stan_data_c <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 11
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
  num_knots = 2 * default_bspline_knots, 
  spline_degree = 3,
  smoothing_strength = 0, # no smoothing.
  prior_scale = 2 * sd(y) # default
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
  num_knots = 4  # Few knots for over-smoothing
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
               alpha = 0.6, linewidth = 1.5) +
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
                fill = "gray80", alpha = 0.3, linewidth = 0.5) +
    labs(title = "Residuals", x = "x", y = "Residuals") +
    theme_bw() +
    theme(plot.title = element_text(size = 9))
  
  return(list(fit = p_fit, resid = p_resid))
}

# Create plots for all 4 scenarios
plots_b_normal <- create_diagnostic_plot(fit_b, x, y, y_true, 
                                         "B-spline: Fitting Settings",
                                         sprintf("%d knots, smoothing=%.2f, autocorr=%.3f", 
                                                 stan_data_b$num_knots,
                                                 stan_data_b$smoothing_strength,
                                                 diagnosis_b$residual_autocor),
                                         "blue")

plots_c_normal <- create_diagnostic_plot(fit_c, x, y, y_true,
                                         "C-spline: Fitting Settings", 
                                         sprintf("%d knots, autocorr=%.3f",
                                                 stan_data_c$num_knots,
                                                 diagnosis_c$residual_autocor),
                                         "darkgreen")

plots_b_overfit <- create_diagnostic_plot(fit_b_overfit, x, y, y_true,
                                          "B-spline: Overfitting",
                                          sprintf("%d knots, smoothing=%.2f, autocorr=%.3f",
                                                  stan_data_b_overfit$num_knots,
                                                  stan_data_b_overfit$smoothing_strength,
                                                  diagnosis_b_overfit$residual_autocor),
                                          "darkred")

plots_c_smooth <- create_diagnostic_plot(fit_c_smooth, x, y, y_true,
                                         "C-spline: Over-smoothing",
                                         sprintf("%d knots, autocorr=%.3f",
                                                 stan_data_c_smooth$num_knots,
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
  # Extract key metrics for each diagnosis
  plot_list <- list()
  
  for (i in seq_along(diagnosis_list)) {
    d <- diagnosis_list[[i]]
    
    # Create data frame for this scenario (removed Residual SD)
    metrics_df <- data.frame(
      Metric = c("Autocorrelation", "Runs Test", "Smoothness"),
      Value = c(d$residual_autocor, d$runs_proportion, 
                ifelse(is.na(d$smoothness), 0, d$smoothness))
    )
    
    # Set colors based on metric type and values (using lighter colors)
    colors <- c(
      ifelse(abs(d$residual_autocor) > 0.3, "#ff6b6b",  # Light red for high positive or negative
             ifelse(abs(d$residual_autocor) > 0.2, "#ffa726", "#66bb6a")),  # Light orange/green
      ifelse(abs(d$runs_proportion - 0.5) > 0.2, "#ff6b6b", 
             ifelse(abs(d$runs_proportion - 0.5) > 0.1, "#ffa726", "#66bb6a")),  # Runs
      ifelse(!is.na(d$smoothness) && d$smoothness > d$residual_sd * 1.5, 
             "#ff6b6b", "#66bb6a")  # Smoothness
    )
    
    # Create individual plot
    p <- ggplot(metrics_df, aes(x = Metric, y = Value, fill = Metric)) +
      geom_col() +
      scale_fill_manual(values = colors) +
      ylim(-0.6, 1.0) +  # Consistent y-axis range
      labs(title = names[i], y = "Value") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 10, hjust = 0.5))
    
    plot_list[[i]] <- p
  }
  
  # Combine all plots
  combined <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] +
    plot_layout(nrow = 1) +
    plot_annotation(
      title = "Diagnostic Metrics by Scenario",
      subtitle = "Green = good, Orange = warning, Red = problem. Each subplot shows one model configuration.",
      caption = "Target values: |Autocorrelation| < 0.2, Runs Test within 0.1 of 0.5, Smoothness < 1.5 × Residual SD",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  return(combined)
}

# Create metrics comparison plot
diagnosis_list <- list(diagnosis_b, diagnosis_c, diagnosis_b_overfit, diagnosis_c_smooth)
scenario_names <- c(
  sprintf("B-spline Fitting\n(%d knots, smoothing=%.2f)", diagnosis_b$num_knots, diagnosis_b$smoothing_strength),
  sprintf("C-spline Fitting\n(%d knots)", diagnosis_c$num_knots),
  sprintf("B-spline Overfit\n(%d knots, smoothing=%.2f)", diagnosis_b_overfit$num_knots, diagnosis_b_overfit$smoothing_strength),
  sprintf("C-spline Over-smooth\n(%d knots)", diagnosis_c_smooth$num_knots)
)

metrics_plot <- create_metrics_plot(diagnosis_list, scenario_names)

ggsave("output/test-diagnostic_metrics_comparison.png", metrics_plot,
       width = 10, height = 8, dpi = 300)
cat("Saved diagnostic metrics plot to output/test-diagnostic_metrics_comparison.png\n")

# Summary
cat("\n\nSummary of Spline Control Parameters\n")
cat("===================================\n")
cat("B-splines:\n")
cat("  - smoothing_strength: Primary control (default 0.1) - adjust this first!\n")
cat("  - num_knots: Auto-selected using n/2 rule (rarely needs adjustment)\n")
cat("  - prior_scale: Auto-scaled to 2*sd(y) (rarely needs adjustment)\n\n")
cat("C-splines:\n")
cat("  - num_knots: Only control parameter (default n/4 rule)\n")
cat("  - Natural boundary conditions built-in\n\n")
cat("The diagnostics guide you to appropriate parameter adjustments.\n")

# Close the null device
dev.off()
