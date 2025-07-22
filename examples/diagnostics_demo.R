# Demonstration of B-spline and C-spline diagnostic capabilities
# Shows how to use the built-in diagnostics to detect over/underfitting

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::var)

library(cmdstanr)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source diagnostic functions
source("code/smoothing_diagnostics.R")

set.seed(42)
dir.create("output", showWarnings = FALSE)

# Generate complex test data
n <- 50
x <- seq(0, 10, length.out = n)
y_true <- sin(x) + 0.3 * cos(3*x)  # Complex function
y <- y_true + rnorm(n, 0, 0.15)

cat("==========================================\n")
cat("Spline Diagnostics Demonstration\n")
cat("==========================================\n\n")

# 1. B-spline with too little smoothing (overfitting)
cat("1. B-spline with insufficient smoothing (overfitting)\n")
cat("------------------------------------------------\n")
bspline_overfit <- fit_bspline(
  x = x, 
  y = y,
  num_knots = 20,  # Many knots
  smoothing_strength = 0.5  # Very little smoothing
)

# 2. B-spline with appropriate smoothing
cat("\n2. B-spline with appropriate smoothing\n")
cat("------------------------------------------------\n")
bspline_good <- fit_bspline(
  x = x, 
  y = y,
  num_knots = 10,
  smoothing_strength = 5.0  # Moderate smoothing
)

# 3. B-spline with too much smoothing (underfitting)
cat("\n3. B-spline with excessive smoothing (underfitting)\n")
cat("------------------------------------------------\n")
bspline_underfit <- fit_bspline(
  x = x, 
  y = y,
  num_knots = 10,
  smoothing_strength = 100  # Very strong smoothing
)

# 4. C-spline with too many knots (overfitting)
cat("\n4. C-spline with too many knots (overfitting)\n")
cat("------------------------------------------------\n")
cspline_overfit <- fit_cspline(
  x = x, 
  y = y,
  num_knots = 15  # Too many knots for C-spline
)

# 5. C-spline with appropriate knots
cat("\n5. C-spline with appropriate knots\n")
cat("------------------------------------------------\n")
cspline_good <- fit_cspline(
  x = x, 
  y = y,
  num_knots = 5  # Good for C-spline
)

# Create comparison plots
plots <- list()

# Extract fitted values for plotting
extract_fitted <- function(result) {
  data.frame(
    x = result$x_plot,
    y_mean = result$y_plot,
    y_lower = result$y_plot_lower,
    y_upper = result$y_plot_upper
  )
}

# B-spline plots
scenarios <- list(
  list(result = bspline_overfit, title = "B-spline: Overfitting\n(20 knots, strength=0.5)"),
  list(result = bspline_good, title = "B-spline: Good fit\n(10 knots, strength=5)"),
  list(result = bspline_underfit, title = "B-spline: Underfitting\n(10 knots, strength=100)"),
  list(result = cspline_overfit, title = "C-spline: Overfitting\n(15 knots)"),
  list(result = cspline_good, title = "C-spline: Good fit\n(5 knots)")
)

for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  df_fit <- extract_fitted(scenario$result)
  
  p <- ggplot() +
    geom_ribbon(data = df_fit, aes(x = x, ymin = y_lower, ymax = y_upper), 
                alpha = 0.3, fill = "blue") +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), 
               alpha = 0.5, size = 1.5) +
    geom_line(data = df_fit, aes(x = x, y = y_mean), 
              color = "blue", linewidth = 1) +
    geom_line(data = data.frame(x = x, y = y_true), aes(x = x, y = y_true), 
              color = "red", linetype = "dashed", alpha = 0.7) +
    labs(title = scenario$title,
         x = "x", y = "y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
  
  plots[[i]] <- p
}

# Combine all plots
combined_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) / 
                 (plots[[4]] | plots[[5]] | plot_spacer()) +
  plot_annotation(
    title = "Spline Fitting Diagnostics Demonstration",
    subtitle = "Blue: fitted spline with 95% CI, Red dashed: true function, Points: data"
  )

ggsave("output/diagnostics_demo_comparison.png", combined_plot, 
       width = 12, height = 8, dpi = 300)

# Print summary of diagnostics
cat("\n\n==========================================\n")
cat("Summary of Diagnostic Outputs\n")
cat("==========================================\n")

cat("\nKey observations:\n")
cat("1. Overfitting indicators:\n")
cat("   - Effective degrees of freedom (EDF) close to number of basis functions\n")
cat("   - Large wiggliness measure\n")
cat("   - Recommendations to increase smoothing or reduce knots\n\n")

cat("2. Underfitting indicators:\n")
cat("   - Low EDF relative to data complexity\n")
cat("   - Systematic patterns in residuals\n")
cat("   - Recommendations to decrease smoothing or increase knots\n\n")

cat("3. Good fit characteristics:\n")
cat("   - Balanced EDF (not too high or low)\n")
cat("   - Residuals appear random\n")
cat("   - No diagnostic warnings\n\n")

cat("Output saved: output/diagnostics_demo_comparison.png\n")