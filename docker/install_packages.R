#!/usr/bin/env Rscript

req_packages <- c("optparse","dplyr","ggplot2","data.table","tidyr","plotly","htmlwidgets","MCMCpack","mvtnorm", "ellipse","igraph")
for (pack in req_packages) {
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack, repos = "http://cran.us.r-project.org")
    }
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack)
    }
}