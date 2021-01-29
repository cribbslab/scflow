library(Seurat)
library(tidyverse)
library(optparse)
library(scRNAseq)
library(celldex)
library(clustifyr)


option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd integrated Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
sample_name <- opt$sample
reference_sce_loc <- opt$reference

# Load seurat object
seurat_object <- readRDS(input_file)

# Load in reference sce object
ref_sce <- readRDS("reference_sce_loc")
ref_mat <- as.matrix(assay(ref_sce))

res <- clustify(input=seurat_object, ref_mat= ref_mat, cluster_col="seurat_clusters", obj_out=TRUE)