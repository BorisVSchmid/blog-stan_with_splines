# Script to update all plots to have white backgrounds
# This will re-run the main plotting scripts with theme_bw() for pure white backgrounds

library(ggplot2)

# Set global theme for all plots in this session
theme_set(theme_bw())

# Update the default theme to have pure white panel background
theme_update(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  legend.background = element_rect(fill = "white", color = NA),
  legend.box.background = element_rect(fill = "white", color = NA)
)

cat("Updated global ggplot2 theme to use white backgrounds\n")
cat("Now re-running main plotting scripts...\n\n")

# Re-run the main plotting scripts
cat("1. Running test_splines.R...\n")
source("test_splines.R")

cat("\n2. Running generate_plots_extended.R...\n")
source("generate_plots_extended.R")

cat("\n3. Running run_regional_splines_fixed.R...\n")
source("run_regional_splines_fixed.R")

cat("\nAll plots have been regenerated with white backgrounds!\n")