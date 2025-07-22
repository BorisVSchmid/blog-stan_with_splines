# Guide to understanding the smoothing_strength parameter
# Shows the intuitive scale: 0 = no smoothing, 10 = strong, 100 = very strong

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)
conflicts_prefer(stats::var)

library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)
dir.create("output", showWarnings = FALSE)

# Generate test data - simple sine wave with noise
n <- 40
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.15)

# Compile model once
cat("Compiling B-spline model...\n")
model <- cmdstan_model("code/bsplines.stan")

# Test different smoothing strengths
smoothing_values <- c(0, 1, 2, 5, 10, 20, 50, 100)
results <- list()

cat("\nFitting models with different smoothing strengths:\n")
cat("==================================================\n")

for (i in seq_along(smoothing_values)) {
  strength <- smoothing_values[i]
  cat(sprintf("Smoothing strength = %3.0f", strength))
  
  stan_data <- list(
    n_data = n,
    x = x,
    y = y,
    num_knots = 12,  # Fixed number of knots
    spline_degree = 3,
    smoothing_strength = strength
  )
  
  fit <- model$sample(
    data = stan_data,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Extract fitted values
  draws <- fit$draws(format = "matrix")
  y_plot <- colMeans(draws[, grep("y_plot\\[", colnames(draws), value = TRUE)])
  x_plot <- colMeans(draws[, grep("x_plot\\[", colnames(draws), value = TRUE)])
  
  # Calculate roughness metric (second derivative proxy)
  roughness <- sum(diff(diff(y_plot))^2)
  
  # Extract effective degrees of freedom
  y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
  
  # Simple EDF approximation based on fitted values complexity
  edf_approx <- sum((y_hat - mean(y_hat))^2) / var(y)
  
  cat(sprintf(" -> Roughness: %6.3f, Approx EDF: %4.1f\n", roughness, edf_approx))
  
  results[[i]] <- list(
    strength = strength,
    x_plot = x_plot,
    y_plot = y_plot,
    roughness = roughness,
    edf = edf_approx,
    description = case_when(
      strength == 0 ~ "No smoothing",
      strength == 1 ~ "Very mild",
      strength == 2 ~ "Mild",
      strength == 5 ~ "Moderate",
      strength == 10 ~ "Strong",
      strength == 20 ~ "Stronger",
      strength == 50 ~ "Very strong",
      strength == 100 ~ "Extremely strong"
    )
  )
}

# Create visualization
plot_data <- map_df(results, function(res) {
  data.frame(
    x = res$x_plot,
    y = res$y_plot,
    strength = res$strength,
    label = sprintf("%s (λ=%g)", res$description, res$strength)
  )
})

# Order labels
plot_data$label <- factor(plot_data$label, 
                         levels = unique(plot_data$label[order(plot_data$strength)]))

# Main comparison plot
p_main <- ggplot(plot_data, aes(x = x)) +
  geom_point(data = data.frame(x = x, y = y), aes(y = y), 
             alpha = 0.3, size = 0.8) +
  geom_line(aes(y = y, color = as.factor(strength)), linewidth = 1) +
  geom_line(data = data.frame(x = x, y = y_true), aes(y = y_true),
            color = "black", linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ label, ncol = 4, scales = "fixed") +
  labs(title = "B-spline Smoothing Strength Guide",
       subtitle = "Points: data, Colored lines: fitted splines, Black dashed: true function",
       x = "x", y = "y") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 9))

# Summary metrics plot
summary_df <- map_df(results, function(res) {
  data.frame(
    strength = res$strength,
    roughness = res$roughness,
    edf = res$edf,
    description = res$description
  )
})

p_roughness <- ggplot(summary_df, aes(x = strength, y = roughness)) +
  geom_line(linewidth = 1.2, color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_text(aes(label = description), hjust = -0.1, vjust = -0.5, 
            size = 3, angle = 45) +
  scale_x_continuous(trans = "log1p", breaks = smoothing_values) +
  scale_y_log10() +
  labs(title = "Roughness vs Smoothing Strength",
       subtitle = "Lower roughness = smoother fit",
       x = "Smoothing strength (log scale)", 
       y = "Roughness (log scale)") +
  theme_bw() +
  expand_limits(x = c(0, 150))

# Combine plots
library(patchwork)
combined_plot <- p_main / p_roughness + 
  plot_layout(heights = c(2, 1))

ggsave("output/smoothing_strength_guide.png", combined_plot, 
       width = 12, height = 12, dpi = 300)

# Print recommendations
cat("\n\nSmoothing Strength Recommendations:\n")
cat("====================================\n")
cat("  0: No smoothing - Use when you want maximum flexibility\n")
cat("  1-2: Very mild - Slight smoothing, preserves most details\n")
cat("  5: Moderate (default) - Good balance for most applications\n")
cat("  10: Strong - Smooth curves, filters out noise\n")
cat("  20-50: Very strong - Very smooth, only major trends\n")
cat("  100+: Extremely strong - Nearly linear behavior\n")

cat("\nChoosing the right value:\n")
cat("- Start with the default (5) and adjust based on results\n")
cat("- If fit is too wiggly: increase smoothing_strength\n")
cat("- If fit is too smooth: decrease smoothing_strength\n")
cat("- Use diagnostics to guide your choice\n")

cat("\nNote: The actual smoothing effect also depends on:\n")
cat("- Number of knots (more knots need more smoothing)\n")
cat("- Data noise level\n")
cat("- Sample size\n")

cat("\nOutput saved: output/smoothing_strength_guide.png\n")