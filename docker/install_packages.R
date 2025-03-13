#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"))
# install requirements from file
packages_to_install = read.csv("my_packages.csv")
# Check which packages are not installed yet
not_installed <- packages_to_install$Package[!packages_to_install$Package %in% installed.packages()[,"Package"]]

# Install only the packages that are not yet installed
if(length(not_installed) > 0) {
  install.packages(not_installed)
  cat("Installed packages:", paste(not_installed, collapse=", "), "\n")
} else {
  cat("All packages are already installed.\n")
}
library("devtools")
install_github("carbocation/aberrant")


