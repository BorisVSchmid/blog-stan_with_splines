# Helper function to create a consistent white background theme
white_theme <- function() {
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "grey80"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "grey95", color = "grey80")
  )
}

# Export the function
cat("White theme function created. Use it in your plots like this:\n")
cat("  ggplot(data, aes(...)) + \n")
cat("    geom_point() + \n")
cat("    white_theme()\n\n")
cat("Or combine with other theme modifications:\n")
cat("  ggplot(data, aes(...)) + \n")
cat("    geom_point() + \n")
cat("    white_theme() + \n")
cat("    theme(legend.position = 'none')\n")