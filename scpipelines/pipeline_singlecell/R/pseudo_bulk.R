#!/usr/bin/env Rscript

library("tidyverse")
library("SingleCellExperiment")
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character",  
              help="input files (a list of sce.rds files) ", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


out <- opt$out
input_list_string <- opt$input

input <- gsub(" ", "", input_list_string)
input <- unlist(strsplit(input, ","))
sample_name_list <- sapply(input, function(x) unlist(strsplit(x, "/"))[2], USE.NAMES = FALSE)

sce <- readRDS(input[1])

pseudo_bulk <- tibble("genes" = rownames(sce))
for(i in 1:length(input)){
  sce.pass <- readRDS(input[i])
  sample_name <- sample_name_list[i]
  pseudo_bulk[[sample_name]] <- as.integer(as.vector(rowSums(counts(sce.pass))))
  
}

pseudo_bulk_df <- column_to_rownames(as.data.frame(pseudo_bulk), var = "genes")
saveRDS(pseudo_bulk_df, file= out)



