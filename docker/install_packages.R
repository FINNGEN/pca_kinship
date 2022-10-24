#!/usr/bin/env Rscript

install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repos = NULL, type="source")


req_packages <- c("plotly","optparse","dplyr","ggplot2","data.table","tidyr","htmlwidgets","MCMCpack","mvtnorm", "ellipse","igraph")
for (pack in req_packages) {
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack, repos = "http://cran.us.r-project.org")
    }
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack)
    }
}