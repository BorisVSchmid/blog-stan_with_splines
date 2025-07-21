# Install required packages
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"), 
                 repos = "https://cloud.r-project.org")

# Install cmdstanr
install.packages("cmdstanr", 
                 repos = c("https://stan-dev.r-universe.dev", 
                          "https://cloud.r-project.org"))

cat("Packages installed. Now you can run the test scripts.\n")