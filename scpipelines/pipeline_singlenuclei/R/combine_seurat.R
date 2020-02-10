#!/usr/bin/env Rscript

library("Seurat")
library("Matrix")
library("tidyverse")
library("SingleCellExperiment")
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character",  
              help="input files (a list of sce.rds files) ", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="output file name", metavar="character"),
  make_option(c("-p", "--pseudo"), type="character", default="alevin", 
              help="pseudoaligner [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


out <- opt$out
input_list_string <- opt$input
pseudo <- opt$pseudo


input <- gsub(" ", "", input_list_string)
input <- unlist(strsplit(input, ","))
print(input)
sample_name_list <- sapply(input, function(x) unlist(strsplit(x, "/"))[2], USE.NAMES = FALSE)

seurat_list <- vector("list", length(input))
gene_list <- c()
for(j in 1:length(input)){
  seurat_object_path <- input[j]
  seurat_object <- readRDS(seurat_object_path)
  #seurat_object <- UpdateSeuratObject(seurat_object)
  seurat_object@meta.data$stim <- sample_name_list[j]
  
  # Normalize data
  seurat_object <- NormalizeData(object = seurat_object)
  
  seurat_object <- ScaleData(seurat_object)

  # detect highly variable genes
  seurat_object <- FindVariableGenes(object = seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
  g <- head(rownames(seurat_object@hvg.info), 1000)

  gene_list <- unique(c(gene_list, g))
  gene_list <- intersect(gene_list, rownames(seurat_object@scale.data))

  seurat_list[[j]] <- seurat_object
                             
}



cca_combined <- RunMultiCCA(object.list = seurat_list, genes.use = gene_list, add.cell.ids = sample_name_list)

saveRDS(cca_combined, file= out)
