#!/usr/bin/env Rscript

#req_packages <- c("optparse","dplyr","ggplot2","data.table","tidyr","plotly","htmlwidgets","MCMCpack","mvtnorm", "ellipse")
#for (pack in req_packages) {
#  if(!require(pack,character.only = TRUE)) {
#    install.packages(pack, repos = c(CRAN = "http://cran.r-project.org"))
#  }
#  require( pack , character.only = TRUE )
#}

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}


script_path <- getScriptPath()
source( paste0(script_path,"/classify_outliers_functions.R"))

if(!require("aberrant")) {
  system( paste0("R CMD INSTALL ",script_path,"/../aberrant/"))
}

library(optparse)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(aberrant)
library(plotly)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="pca file as output from plink without header. ", metavar="character"),
  make_option(c("-s", "--sample_annot"), type="character", default=NULL, 
              help="Two column file where first is individual id and second is a label for the sample. Use an arbitraty header . Used for plotting", metavar="character"),
  make_option(c("-e", "--eigen"), type="character", default=NULL, 
              help="eigenvalues as output from plink", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="output file prefix", metavar="character"),
  make_option(c("-n", "--n_pcs"), type="integer", default=3, 
              help="Number of PCs to use for outlier detection",
              metavar="number"),
  make_option(c("--pc_list"), type="character", default=NULL, 
              help="Comma separated list of PCs to use. Overrides --n_pcs setting"),
  make_option(c("-p", "--skip_plotting"), action="store_true", default=FALSE,
           help="Skip plotting [default]"),
  
  make_option(c("--n_iterations"), type="integer", default=10000,
              help="number of iterations to run MCMC [default]")
  
make_option(c("--lambda"), type="integer", default=20,
              help="number of iterations to run MCMC [default]")
  

  )
  

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


pca <- fread( opt$f)
setkey(pca,IID )
colnames(pca) <- c("FID","IID", paste0("PC", seq( ncol(pca)-2)) )

batch_col <- NULL
pca$batch <- ""
if( ! is.null( opt$sample_annot  )) {
  samples <- fread(opt$sample_annot, header=T)
  setkeyv(samples, colnames(samples)[1])
  pca <- merge(pca, samples, by.x="IID", by.y= colnames(samples)[1], all.x=T)
  
  print(colnames(samples)[2])
  batch_col <- colnames(samples)[2]
  pca$batch <- subset(pca, select=batch_col)
  
}

eigen <- NULL
if( ! is.null( opt$eigen  )) {
  eigen <- read.table(opt$eigen)
}

pca$outlier <- FALSE

outcols <- c("IID","outlier")

pc_vec <- c()

if(! is.null(opt$pc_list) ) {
  pc_vec <- unlist(strsplit( opt$pc_list, "," ))
} else {
  pc_vec <- paste0("PC", 1:3)
}



for( i in 1:length(pc_vec)) {
  if( i==length(pc_vec) &  length(pc_vec)>1) {
    break
  }

  for (j in (i+1):length(pc_vec)) {
    cols <- c(pc_vec[i],pc_vec[j])
    if ( j> length(pc_vec) ) {
      ## only one pc given
      cols <- c(pc_vec[i],pc_vec[i])
    }

    print(cols)
    dat <-  subset(pca, select=cols)
    print(summary(dat))
    outliers_lambda20 <- aberrant( dat,  lambda=opt$lambda, niter = opt$n_iterations)
    pca[outliers_lambda20$group==1, ]$outlier <- T
    outcols[[length(outcols)+1]] <- paste0("OUTLIERS_", cols,  collapse="_")
    pca[, paste0("OUTLIERS_", cols, collapse = "_")] <- outliers_lambda20$group==1
  }
}

if(! opt$skip_plotting) {
  pdf( paste0( opt$o,"_outlier_pcas.pdf" ))
  basic_plot_pcs(pca, eigen_values=eigen$V1, colour_by ="outlier", n_first=4)
  if( ! is.null(batch_col)) {
    basic_plot_pcs(pca, eigen_values=eigen$V1, colour_by =batch_col, n_first=4)  
  }
  dev.off()
  
  
  if(system('pandoc -v', show.output.on.console = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)==0) {
    results <- tryCatch({
      pl <- plot_ly(pca, x=~PC1, y=~PC2, color=~outlier,hoverinfo='text', text=~paste0(IID,"<BR>",batch))
      saveWidget(pl, file=paste0( opt$o,"_outlier_pcas.html"))
    }, 
    error = function(x) print( paste0(x,". Install pandoc to enable interactive html plotting." ) )  
    )
  } else  {
    print("For interactive plots install pandoc in your system. Skipping interactive plot generation")
  }

}


if( ! is.null( batch_col)) {
  
  sums <- full_join( pca %>% filter(outlier) %>% group_by_(.dots = batch_col) %>% summarize( n_outlier_lambda20=n() ) %>% replace_na( list(n_outlier_lambda20=0) ), 
             pca %>% filter(!outlier) %>% group_by_(.dots = batch_col) %>% summarize( n_inlier_lambda20=n() ) %>% replace_na( list(n_inlier_lambda20=0) ), by=c(batch_col) ) 
  
  write.table(sums,file=paste0( opt$o,"_outlier_summaries.tsv") , row.names = F, quote=F, sep="\t")
}

if(! is.null( batch_col) ) {
  outcols <- c(outcols, batch_col)
}

write.table(  subset(pca, select=outcols), file=paste0( opt$o,"_outlier_samples.tsv"), row.names = F, quote=F, sep="\t" )


