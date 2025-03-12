#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"))

install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repos = NULL, type="source")
# install requirements from file
packages_to_install = read.csv("my_packages.csv")
install.packages(packages_to_install$Package)
library("devtools")
install_github("carbocation/aberrant")