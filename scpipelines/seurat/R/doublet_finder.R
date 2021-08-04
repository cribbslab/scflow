library(tidyverse)
library(Seurat)
library(optparse)
library(scDblFinder)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds file and path, filtered Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name")
)


opt <- parse_args(OptionParser(option_list=option_list))

sample_name <- opt$sample
sample_name <- gsub(" ", "", sample_name, fixed = TRUE)

# Read in RDS files (may take some time)
seurat_filtered_rds <- opt$input
filtered_seurat_object <- readRDS(seurat_filtered_rds)

sce <- as.SingleCellExperiment(filtered_seurat_object)
sce <- scDblFinder(sce, dbr=0.1)


outfile <- paste0("Doublet_Figures.dir/", opt$sample, "_doublets.csv")

write.csv(table(sce$scDblFinder.class), file=outfile, row.names = F)

sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

metadata <- sce.seurat@meta.data

out_meta <- paste0("Doublet_Figures.dir/", opt$sample, "_metadata.csv")

write.csv(metadata, file=out_meta, row.names = T)
