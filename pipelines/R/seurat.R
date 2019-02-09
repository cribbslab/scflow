#!/usr/bin/env Rscript

library("tidyverse")
library("SingleCellExperiment")
library("optparse")
library("Seurat")

option_list = list(
  make_option(c("-w", "--workingdir"), type="character", default=NULL, 
              help="working directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default="sce.rds", 
              help="input file (quants_mat.gz path) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="sce.rds", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wd = opt$workingdir

setwd(wd)
out <- opt$out
input <- opt$input

sce_object = readRDS(input)

seurat_object = Convert(from = sce_object, to = "seurat")


saveRDS(object = seurat_object, file = opt$out)