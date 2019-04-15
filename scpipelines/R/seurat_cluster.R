#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")
library("Seurat")

option_list = list(
  make_option(c("-w", "--workingdir"), type="character", default=".", 
              help="working directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default="sce.rds", 
              help="input file containing the seurat rds object"),
  make_option(c("-o", "--out"), type="character", default="sce.rds", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--reductionmethod"), default="pca",
                help="Type of dimensionality reduction to perform. Options include 'pca', 'ica'"),
  make_option(c("--pccomponents"), type="integer", default=10,
                help="The number of principle components to use"),
  make_option(c("--resolution"), type="double", default=1,
                help="cluster resolution"),
  make_option(c("--algorithm"), type="integer", default=3,
                help="1=original Louvain, 2=Louvain multilevel, 3=SLM")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

so <- readRDS(opt$input)

comps <- 1:as.numeric(opt$pccomponents)


so <-  FindClusters(so,
                  reduction.type=opt$pccomponents,
                  dims.use = comps,
                  resolution = opt$resolution,
                  algorithm = opt$algorithm,
                  print.output = FALSE,
                  save.SNN = F)

nclusters <- length(unique(so@ident))


so <- buildClusterTree(so,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:11)

# TODO: impliment functions that allow more indepth interrogation
# of clustering relationship. For example, cluster correlations heatmap


saveRDS(object = so, file = opt$out)