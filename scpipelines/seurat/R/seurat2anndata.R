#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library("optparse")
library(stringr)


option_list = list(
  make_option(c("-i", "--input"), type="character", default="", 
              help="input seurat object [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="", 
              help="output file name for anndata [default= %default]", metavar="character")
	      ); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

so <- readRDS(opt$input)

name <- str_replace(opt$input, "RDS_objects.dir/", "")
name <- str_replace(name, "filtered_SeuratObject.rds", "")

so_name <- paste0("Anndata.dir/", name, ".h5Seurat")
h5_name <- paste0("Anndata.dir/", name, ".h5ad")

SaveH5Seurat(so, filename = so_name)
Convert(so_name, dest = h5_name)