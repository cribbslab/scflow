library(Seurat)
library(tidyverse)
library(optparse)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds file and path, filtered Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option(c("-m", "--minPct"), default=0.1, type="double",
			help="Genes only tested if found in minimum percentage of cells in either population [default %default]"),
		make_option(c("-l", "--logfc"), default=0.25, type ="double",
			help="Limit testing to genes which have (on average) a log fold change greater than this threshold [default %default]"),
        make_option(c("-t", "--testuse"), default="wilcox", type="character",
			help="Test to use. Options: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2. [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

sample_name <- read.table(opt$sample,header = TRUE, fill = T)
min_pct <- opt$minPct
logfc_threshold <- opt$logfc
test_use <- opt$testuse
seurat_object_path <- read.table(opt$input,header = TRUE, fill = T)

# Read in RDS files (may take some time)
seurat_object <- readRDS(seurat_object_path

# Might need to add more stat stuff

# Something like this... need to look at object to see meta data.
for(i in 1:length(clusters)){
    so_markers <- cluster1.markers <- FindMarkers(seurat_object, ident.1 = i, min.pct = min_pct,
                                                  logfc.threshold = logfc_threshold, test.use = test_use)
    assign(paste("so_markers", i, sep = "."), so_markers)

	# cbind all data ? rownames to column for genes. Add column for cluster number
	# Then save as tsv file
}


