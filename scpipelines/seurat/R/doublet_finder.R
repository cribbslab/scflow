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

# save output on number of doublets identified
outfile <- paste0("Doublet_Figures.dir/", opt$sample, "_doublets.csv")

write.csv(table(sce$scDblFinder.class), file=outfile, row.names = F)

# save the RDS file
saveRDS(sce, gsub("SAMPLE_FILE",sample_name ,"RDS_objects.dir/SAMPLE_FILE_filtered_clustered_doublet_SeuratObject.rds"))
