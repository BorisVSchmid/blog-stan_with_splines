# Test B-spline smoothing behavior on sine waves
# Explores how smoothing strength affects oscillatory function fitting

# Suppress default graphics device to prevent Rplots.pdf
pdf(NULL)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::sd)

library(groundhog)
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2", "patchwork", "tidyr")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

library(cmdstanr)

set.seed(123)

cat("Testing B-spline smoothing on sine waves\n")
cat("=========================================\n\n")

# Function to fit B-spline with specific configuration
fit_bspline_config <- function(x, y, num_knots, smoothing_strength) {
  stan_data <- list(
    n_data = length(x),
    x = x,
    y = y,
    num_knots = num_knots,
    spline_degree = 3,
    smoothing_strength = smoothing_strength,
    prior_scale = 2 * sd(y)
  )
  
  model <- cmdstan_model("code/bsplines.stan")
  fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    max_treedepth = 12,
    refresh = 0,
    adapt_delta = 0.99
  )
  
  # Extract results
  draws <- fit$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  list(
    x_plot = x_plot,
    y_plot = y_plot,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper,
    stan_data = stan_data
  )
}

# Function to fit C-spline
fit_cspline_config <- function(x, y, num_knots) {
  stan_data <- list(
    n_data = length(x),
    x = x,
    y = y,
    num_knots = num_knots
  )
  
  model <- cmdstan_model("code/csplines.stan")
  fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    max_treedepth = 12,
    adapt_delta = 0.99,
    refresh = 0
  )
  
  # Extract results
  draws <- fit$draws(format = "matrix")
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  y_plot_lower <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
  y_plot_upper <- apply(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)
  
  list(
    x_plot = x_plot,
    y_plot = y_plot,
    y_lower = y_plot_lower,
    y_upper = y_plot_upper,
    stan_data = stan_data
  )
}

# Test parameters
smoothing_values <- c(0, 0.1, 0.2, 0.3, 0.4)
cspline_knots <- c(11, 9, 7, 5, 3)
sample_sizes <- c(20, 40, 60)  # For B-spline columns

# Store all results
all_results <- list()

# Generate data for each sample size
for (n in sample_sizes) {
  x <- seq(0, 10, length.out = n)
  y_true <- sin(x)
  y <- y_true + rnorm(n, 0, 0.1)
  all_results[[paste0("data_n", n)]] <- list(x = x, y = y, y_true = y_true)
}

# Column 0: C-splines with n=20 and varying knots
cat("\nProcessing column 0: C-splines with n=20 and varying knots...\n")
data_n20 <- all_results[["data_n20"]]
for (i in seq_along(cspline_knots)) {
  knots <- cspline_knots[i]
  cat(sprintf("  C-spline with %d knots\n", knots))
  result <- fit_cspline_config(data_n20$x, data_n20$y, num_knots = knots)
  all_results[[paste0("cspline_knots", i)]] <- result
  all_results[[paste0("cspline_knots", i, "_value")]] <- knots
}

# Columns 1-3: B-splines with different sample sizes and varying smoothing
for (j in seq_along(sample_sizes)) {
  n <- sample_sizes[j]
  cat(sprintf("\nProcessing column %d: B-splines with n=%d...\n", j, n))
  data <- all_results[[paste0("data_n", n)]]
  
  # Use adaptive knots based on sample size (n/2, capped at 40)
  adaptive_knots <- max(4, min(round(n/2), 40))
  cat(sprintf("  Using %d knots (n/2 rule)\n", adaptive_knots))
  
  for (i in seq_along(smoothing_values)) {
    smooth <- smoothing_values[i]
    cat(sprintf("  B-spline with smoothing = %.1f\n", smooth))
    result <- fit_bspline_config(data$x, data$y, num_knots = adaptive_knots, smoothing_strength = smooth)
    all_results[[paste0("bspline_n", n, "_smooth", i)]] <- result
    all_results[[paste0("bspline_n", n, "_smooth", i, "_value")]] <- smooth
    all_results[[paste0("bspline_n", n, "_knots")]] <- adaptive_knots
  }
}

# Create the 5x4 plot (4 columns, 5 rows)
cat("\n\nCreating comparison plot...\n")

library(ggplot2)
library(patchwork)

# Create plots list
plots <- list()
plot_idx <- 1

# Create plots row by row for better organization
for (row in 1:5) {
  # Column 1: C-splines
  data_n20 <- all_results[["data_n20"]]
  cspline_result <- all_results[[paste0("cspline_knots", row)]]
  knots_value <- all_results[[paste0("cspline_knots", row, "_value")]]
  
  if (!is.null(cspline_result)) {
    df_cspline <- data.frame(
      x = cspline_result$x_plot,
      y = cspline_result$y_plot,
      y_lower = cspline_result$y_lower,
      y_upper = cspline_result$y_upper
    )
    
    p_cspline <- ggplot() +
      geom_ribbon(data = df_cspline, aes(x = x, ymin = y_lower, ymax = y_upper), 
                  alpha = 0.3, fill = "darkgreen") +
      geom_line(data = df_cspline, aes(x = x, y = y), color = "darkgreen", linewidth = 1) +
      geom_point(data = data.frame(x = data_n20$x, y = data_n20$y), 
                 aes(x = x, y = y), size = 1.5, alpha = 0.6) +
      geom_line(data = data.frame(x = data_n20$x, y = data_n20$y_true), 
                aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.8) +
      labs(title = if(row == 1) "C-spline (n=20)" else "",
           subtitle = paste0("knots=", knots_value)) +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 9)) +
      ylim(-1.5, 1.5)
    
    plots[[plot_idx]] <- p_cspline
  }
  plot_idx <- plot_idx + 1
  
  # Columns 2-4: B-splines for each sample size
  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    data_info <- all_results[[paste0("data_n", n)]]
    bspline_result <- all_results[[paste0("bspline_n", n, "_smooth", row)]]
    smooth_value <- smoothing_values[row]
    
    if (!is.null(bspline_result)) {
      df_bspline <- data.frame(
        x = bspline_result$x_plot,
        y = bspline_result$y_plot,
        y_lower = bspline_result$y_lower,
        y_upper = bspline_result$y_upper
      )
      
      p_bspline <- ggplot() +
        geom_ribbon(data = df_bspline, aes(x = x, ymin = y_lower, ymax = y_upper), 
                    alpha = 0.3, fill = "blue") +
        geom_line(data = df_bspline, aes(x = x, y = y), color = "blue", linewidth = 1) +
        geom_point(data = data.frame(x = data_info$x, y = data_info$y), 
                   aes(x = x, y = y), size = 1.5, alpha = 0.6) +
        geom_line(data = data.frame(x = data_info$x, y = data_info$y_true), 
                  aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.8) +
        labs(title = if(row == 1) paste0("B-spline (n=", n, ", knots=", all_results[[paste0("bspline_n", n, "_knots")]], ")") else "",
             subtitle = paste0("smoothing=", smooth_value)) +
        theme_bw() +
        theme(plot.title = element_text(size = 10, face = "bold"),
              plot.subtitle = element_text(size = 9)) +
        ylim(-1.5, 1.5)
      
      plots[[plot_idx]] <- p_bspline
    }
    plot_idx <- plot_idx + 1
  }
}

# Ensure we have exactly 20 plots (fill with empty if needed)
if (length(plots) < 20) {
  empty_plot <- ggplot() + theme_void()
  for (i in (length(plots)+1):20) {
    plots[[i]] <- empty_plot
  }
}

# Combine all plots
combined_plot <- wrap_plots(plots, ncol = 4, nrow = 5) +
  plot_annotation(
    title = "Spline Smoothing Effects on Sine Wave Fitting",
    subtitle = "Black dashed = true sin(x), points = data. Col 0: C-splines (n=20) varying knots. Cols 1-3: B-splines (n/2 knots) varying smoothing",
    caption = "Key finding: B-splines need minimal smoothing (0-0.4) for oscillatory functions. C-splines robust across knot counts."
  )

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/test-sine_wave_smoothing.png", combined_plot, width = 14, height = 12.5, dpi = 300)
cat("Saved plot to output/test-sine_wave_smoothing.png\n")

# Print summary
cat("\n\nSummary of Results\n")
cat("==================\n")

# C-spline results
cat("\nC-splines (n=20): Knot effects on fit quality\n")
for (i in seq_along(cspline_knots)) {
  knots <- cspline_knots[i]
  result <- all_results[[paste0("cspline_knots", i)]]
  y_true_plot <- sin(result$x_plot)
  rmse <- sqrt(mean((result$y_plot - y_true_plot)^2))
  cat(sprintf("  Knots = %d: RMSE = %.3f\n", knots, rmse))
}

# B-spline results for each sample size
for (n in sample_sizes) {
  knots <- all_results[[paste0("bspline_n", n, "_knots")]]
  cat(sprintf("\nB-splines (n=%d, %d knots): Smoothing effects on fit quality\n", n, knots))
  for (i in seq_along(smoothing_values)) {
    smooth <- smoothing_values[i]
    result <- all_results[[paste0("bspline_n", n, "_smooth", i)]]
    y_true_plot <- sin(result$x_plot)
    rmse <- sqrt(mean((result$y_plot - y_true_plot)^2))
    cat(sprintf("  Smoothing = %.1f: RMSE = %.3f\n", smooth, rmse))
  }
}

# Close the null device
invisible(dev.off())

cat("\nTest complete!\n")
