library(Seurat)
library(tidyverse)
library(optparse)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option(c("-c", "--maxClusters"), default=0, type="integer",
                        help="If you want to set a maximum number of clusters to find markers for. 0 = Do it for all clusters [default %default]"),
        make_option(c("-g", "--group"), default="stim", type = "character",
			help="Name of grouping variable. [default %default]"),
        make_option(c("-m", "--meta"), default="metadata.csv", type = "character",
			help="Name of meta data csv file [default %default]"),
        make_option(c("--de"), type = "character",
			help="Differential expression versus conditions.")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

seurat_object_path <- opt$input
sample_name <- opt$sample
max_clusters <- opt$maxClusters
group_var <- opt$group
meta_data_path <- opt$meta

de_conditions_long <- opt$de
de_conditions <- str_split(de_conditions_long, "-")[[1]] # In for loop split by underscore to get conditions

###############################
# Read in and format meta data
###############################



# Read in RDS files (may take some time)
seurat_object <- readRDS(seurat_object_path)

# Get total number of clusters
num_clusters <- nlevels(seurat_object$seurat_clusters)

# Get minimum between number of clusters selected and total number of clusters
if(max_clusters){
	num_clusters_to_use <- min(c(max_clusters, final_cluster))
	final_cluster <- num_clusters_to_use - 1
} else{
	final_cluster <- num_clusters - 1
}

DefaultAssay(seurat_object) <- "RNA"

############################################
# Identify conserved cell type markers
############################################

combined <- c()
combined_topmarkers <- c()

for(i in 0:final_cluster){

  conserved_markers <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = group_var, verbose = FALSE)
  top9 <- rownames(conserved_markers)[1:9]
  top2 <- rownames(conserved_markers)[1:2]
  combined_topmarkers <- c(combined_topmarkers, top2)

  name<- paste0("Annotation_Figures.dir/FeaturePlot_Top9ConservedMarkers_", sample_name, "_cluster_", i, ".eps")
  postscript(name)
  print(FeaturePlot(seurat_object, features = top9, min.cutoff = "q9"))
  dev.off()

  conserved_markers$ensembl <- rownames(conserved_markers)
  conserved_markers$cluster <- i

  conserved_markers_tib <- as_tibble(so_markers)
  combined <- rbind(combined, conserved_markers_tib)


}

combined_topmarkers <- unique(combined_topmarkers)
# Dotplot showing 2 top strong marker genes for each cluster
dotplot <- DotPlot(seurat_object, features = rev(combined_topmarkers), cols = c("blue", "red"), dot.scale = 8,
    split.by = group_var) + RotatedAxis()
name <- paste0("Annotation_Figures.dir/DotPlot_TopConservedMarkers_", sample_name, ".eps")
postscript(name)
print(dotplot)
dev.off()


name_file <- paste0("Annotation_stats.dir/ConservedMarkers_", sample_name, ".csv")

combined <- combined %>% dplyr::select(ensembl, cluster, everything())
write_csv(combined, name_file)

############################################
# Differential expression across conditions
############################################

