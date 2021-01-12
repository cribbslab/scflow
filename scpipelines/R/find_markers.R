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

sample_name <- opt$sample
min_pct <- opt$minPct
logfc_threshold <- opt$logfc
test_use <- opt$testuse
seurat_object_path <- opt$input

# Read in RDS files (may take some time)
seurat_object <- readRDS(seurat_object_path)

# Get total number of clusters
num_clusters <- nlevels(seurat_object$seurat_clusters)
final_cluster <- num_clusters - 1

# Start clustering from 1 not 0
#if(!(as.numeric(levels(seurat_object$seurat_clusters)[1]))){
#	 seurat_object$seurat_clusters <- as.factor(as.numeric(as.character(seurat_object$seurat_clusters)) + 1)
#
#}

# Might need to add more stat stuff

# Empty vector to combine results into
combined <- c()

for(i in 0:final_cluster){
    so_markers <-  FindMarkers(seurat_object, ident.1 = i, min.pct = min_pct,
                                                  logfc.threshold = logfc_threshold, test.use = test_use)
    assign(paste("so_markers", i, sep = "."), so_markers)

	# Get back, probably don't need
	so <- get(gsub("CLUSTER",i , "so_markers.CLUSTER"))
	so$ensembl <- rownames(so)
	# Get rid of .1 .12 etc. at end of ensembl name
	so$ensembl_short <-   gsub("\\.\\d+", "", so$ensembl)
	so$cluster <- i

	so <- as_tibble(so)
	combined <- rbind(combined, so)


}

combined <- combined %>% dplyr::select(ensembl, ensembl_short, cluster, everything())

combined_logfc_order <- combined %>% dplyr::arrange(cluster, desc(abs(avg_logFC)), p_val_adj)
tsv_filename <- paste0("clustering_markers.dir/", sample_name, "_markers.tsv")
tsv_filename_log <- paste0("clustering_markers.dir/", sample_name, "logfc_ordered_markers.tsv")


write_csv(combined, tsv_filename)
write_csv(combined_logfc_order, tsv_filename_log)

