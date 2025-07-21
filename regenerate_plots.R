# Quick script to regenerate key plots with white backgrounds
library(ggplot2)
library(dplyr)

# First, let's regenerate the regional data plot from demo_regional_splines.R
set.seed(42)

# Generate synthetic data for 3 regions
n_per_region <- 30
x_vals <- seq(0, 10, length.out = n_per_region)

# Region A: Standard sine wave
y_A <- sin(x_vals) + rnorm(n_per_region, 0, 0.15)

# Region B: Shifted and scaled sine wave  
y_B <- 0.7 * sin(x_vals + 0.5) + 0.3 + rnorm(n_per_region, 0, 0.15)

# Region C: Different frequency
y_C <- sin(1.5 * x_vals) - 0.2 + rnorm(n_per_region, 0, 0.15)

# Combine data
data <- data.frame(
  x = rep(x_vals, 3),
  y = c(y_A, y_B, y_C),
  region = rep(1:3, each = n_per_region),
  region_name = factor(rep(c("Region A", "Region B", "Region C"), each = n_per_region))
)

# Plot raw data with white background
p_data <- ggplot(data, aes(x = x, y = y, color = region_name)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ region_name, ncol = 1, scales = "free_y") +
  labs(title = "Regional Data: 3 Different Patterns",
       subtitle = "Each region follows a different nonlinear pattern",
       x = "x", y = "y") +
  theme_bw() +  # White background with grid lines
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

# Save plot
dir.create("output", showWarnings = FALSE)
ggsave("output/regional_data.png", p_data, width = 8, height = 8, dpi = 300)

cat("Regenerated regional_data.png with white background\n")

# If you want to regenerate other specific plots, run:
# source("run_regional_splines_fixed.R")  # This will refit models and regenerate all regional plots
# source("generate_plots_extended.R")      # This will regenerate extended boundary plots