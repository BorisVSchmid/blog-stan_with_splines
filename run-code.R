if (!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted")
if (!requireNamespace("groundhog", quietly = TRUE)) install.packages("groundhog")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("cmdstanr", quietly = TRUE)) install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Might have to restart a few times if this is your first time using groundhog.

# Might have to run the install script for cmdstanr if you don't have it installed already.
#library(cmdstanr)
#install_cmdstan(overwrite=TRUE)

# Pay attention to the messages.
# For example, I get a please add /c/Users/mail/.cmdstan/cmdstan-2.36.0/stan/lib/stan_math/lib/tbb to your PATH.

# Run minimal code examples
# Simple entry point for basic spline demonstrations

cat("Stan Splines - Minimal Examples\n")
cat("==============================\n\n")

# Run B-spline example
cat("--- Running B-spline example with key features ---\n")
source("code/bsplines.R")

# Run C-spline example
cat("\n--- Running C-spline minimal example ---\n")
source("code/csplines.R")

cat("\nMinimal examples complete. Check output/ directory for plots.\n")