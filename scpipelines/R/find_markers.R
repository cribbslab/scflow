library(Seurat)
library(tidyverse)
library(optparse)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option(c("-m", "--minPct"), default=0.1, type="double",
			help="Genes only tested if found in minimum percentage of cells in either population [default %default]"),
		make_option(c("-l", "--logfc"), default=0.25, type ="double",
			help="Limit testing to genes which have (on average) a log fold change greater than this threshold [default %default]"),
        make_option(c("-t", "--testuse"), default="wilcox", type="character",
			help="Test to use. Options: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2. [default %default]"),
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

# Percentage differences
genes <- rownames(seurat_object@assays$RNA@data)[Matrix::rowSums(seurat_object@assays$RNA@data)>0]

# Function
expmeanDCG <- function(mat,iter) {
	df<-NULL
	for(n in seq(0,nrow(mat),iter)){
		init <- n+1
		f<-n+iter
		if(f>nrow(mat)){ f<- nrow(mat)}
		message(paste0("range [", init, ":",f,"]") )
		df[init:f] <- apply(mat[init:f, ],1,FUN=ExpMean )
		}
	names(df) <- rownames(mat)
	return(df)
}

# Empty vector to combine results into
combined <- c()
filter_stats_combined <- c()
markers_filter_stats_combined <- c()

for(i in 0:final_cluster){

	# Which cells are in each cluster
	cluster_cells <- WhichCells(seurat_object,idents=i)
	other_clusters <- setdiff(0:final_cluster,i)
	other_cells <- WhichCells(seurat_object,idents=other_clusters)

	# Percentages and differences
	cluster_pct <- round(Matrix::rowSums(seurat_object@assays$RNA@data[genes,cluster_cells, drop=F]>0)/length(cluster_cells),digits=3)
	other_pct <- round(Matrix::rowSums(seurat_object@assays$RNA@data[genes,other_cells, drop=F]>0)/length(other_cells),digits=3)

	pcts <- cbind(cluster_pct,other_pct)
	max_pct <- apply(pcts,1,max)
	min_pct <- apply(pcts,1,min)
	diff_pct <- max_pct - min_pct

	cluster_mean <- expmeanDCG(seurat_object@assays$RNA@data[genes,cluster_cells],2000)
	other_mean <- expmeanDCG(seurat_object@assays$RNA@data[genes,other_cells],2000)
	diff_mean <- abs(cluster_mean - other_mean)


	so_filter_stats <- data.frame(rownames = genes, ensembl = genes, cluster = i, cluster_mean=signif(cluster_mean,4), other_mean=signif(other_mean,4),cluster_pct = cluster_pct,other_pct = other_pct)
	filter_stats_combined <- rbind(filter_stats_combined, so_filter_stats)

	take_min_pct <- max_pct > min_pct

    so_markers <-  FindMarkers(seurat_object, ident.1 = i, min.pct = min_pct, features = genes,
                                                  logfc.threshold = logfc_threshold, test.use = test_use)

	if(nrow(so_markers) > 0){
		so_markers$ensembl <- rownames(so_markers)
		# Get rid of .1 .12 etc. at end of ensembl name
		so_markers$ensembl_short <-   gsub("\\.\\d+", "", so_markers$ensembl)
		so_markers$cluster <- i
	}

	so <- as_tibble(so_markers)
	combined <- rbind(combined, so)

	full_marker_stats <- as_tibble(cbind(so_markers,so_filter_stats[rownames(so_markers),]))
	markers_filter_stats_combined <- rbind(markers_filter_stats_combined, full_marker_stats )


}

# BH corrected p-value, change their p_val_adj to bonferroni
combined$p.adj <- p.adjust(combined$p_val, method="BH")
combined$p.adj.bonferroni <- combined$p_val_adj
combined <- combined %>% dplyr::select(ensembl, ensembl_short, cluster, p.adj, p_val, p.adj.bonferroni, avg_logFC, pct.1, pct.2)

# Rearrange and select columns
combined_padj_order <- combined %>% dplyr::arrange(cluster, p.adj, desc(abs(avg_logFC)))
combined_logfc_order <- combined %>% dplyr::arrange(cluster, desc(abs(avg_logFC)), p_val_adj)

# All together, test and see which output is best. Need to look out for bugs in this
markers_filter_stats_combined$p.adj <- p.adjust(markers_filter_stats_combined$p_val, method="BH")
markers_filter_stats_combined$p.adj.bonferroni <- markers_filter_stats_combined$p_val_adj
markers_filter_stats_combined <- markers_filter_stats_combined %>% dplyr::select(ensembl, ensembl_short, cluster, p.adj, p_val, p.adj.bonferroni, avg_logFC, pct.1, pct.2,
	cluster_mean, other_mean, cluster_pct ,other_pct)
markers_filter_stats_combined <- markers_filter_stats_combined %>% dplyr::arrange(cluster, p.adj, desc(abs(avg_logFC)))


# Make file names
csv_filename <- paste0("clustering_markers.dir/", sample_name, "_markers.tsv")
csv_filename_log <- paste0("clustering_markers.dir/", sample_name, "logfc_ordered_markers.tsv")
filter_stats_filename <- paste0("clustering_markers.dir/", sample_name, "_filter_stats.tsv")
csv_filename_all <- paste0("clustering_markers.dir/", sample_name, "_markers_filter_stats.tsv")

# Save CSVs
write_csv(combined_padj_order, csv_filename)
write_csv(combined_logfc_order, csv_filename_log)
write_csv(filter_stats_combined, filter_stats_filename)
write_csv(markers_filter_stats_combined, csv_filename_all)