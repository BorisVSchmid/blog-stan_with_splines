# Example showing B-spline with smoothing enabled

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

library(groundhog)
# Use groundhog for dependencies
stan_pkgs <- c("posterior", "checkmate", "R6", "jsonlite", "processx")
pkgs <- c("dplyr", "ggplot2")
groundhog.library(c(stan_pkgs, pkgs), "2025-06-01")

# Load cmdstanr after its dependencies
library(cmdstanr)

# Generate simple test data with more noise
n <- 30
x <- seq(0, 10, length.out = n)
y_true <- sin(x)
y <- y_true + rnorm(n, 0, 0.2)  # More noise

# Prepare data for Stan with smoothing
stan_data <- list(
  n_data = n,
  x = x,
  y = y,
  num_knots = 10,     # More knots but with smoothing
  spline_degree = 3,  
  smoothing_strength = 5.0    # Moderate smoothing (was tau_smooth = 0.2)
)

# Compile and fit model
model <- cmdstan_model("code/bsplines.stan")

fit <- model$sample(
  data = stan_data,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# Extract and plot results
draws <- fit$draws(format = "matrix")
y_hat <- colMeans(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)])
y_hat_lower <- apply(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)], 2, quantile, 0.025)
y_hat_upper <- apply(draws[, grep("y_hat\\[", colnames(draws), value = TRUE)], 2, quantile, 0.975)

# Create data frame for plotting
plot_data <- data.frame(
  x = x,
  y = y,
  y_hat = y_hat,
  y_hat_lower = y_hat_lower,
  y_hat_upper = y_hat_upper,
  y_true = y_true
)

# Create plot with ggplot2
p <- ggplot(plot_data, aes(x = x)) +
  geom_ribbon(aes(ymin = y_hat_lower, ymax = y_hat_upper), alpha = 0.3, fill = "blue") +
  geom_point(aes(y = y), color = "black", size = 2, alpha = 0.7) +
  geom_line(aes(y = y_hat), color = "blue", linewidth = 1.2) +
  geom_line(aes(y = y_true), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "B-spline Fit with Smoothing",
    subtitle = paste0("num_knots = ", stan_data$num_knots, 
                      ", degree = ", stan_data$spline_degree,
                      ", smoothing_strength = ", stan_data$smoothing_strength),
    x = "x",
    y = "y"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(p)

# Save the plot
dir.create("output", showWarnings = FALSE)
ggsave("output/example-bspline_smooth_example.png", p, width = 8, height = 6, dpi = 300)