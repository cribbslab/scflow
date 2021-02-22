library(Seurat)
library(tidyverse)
library(optparse)
library(RColorBrewer)

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
			help="Differential expression versus conditions."),
        make_option(c("--predefined"), defualt = NULL,
			help="User pre-defined list of marker genes")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

seurat_object_path <- opt$input
sample_name <- opt$sample
max_clusters <- opt$maxClusters
group_var <- opt$group
meta_data_path <- opt$meta

de_conditions_long <- opt$de
de_conditions <- str_split(de_conditions_long, "---")[[1]] # In for loop split by underscore to get conditions

#################################################
# Read in and format meta data and seurat object
#################################################

# Read in SO RDS files (may take some time)
seurat_object <- readRDS(seurat_object_path)

meta <- read_csv(meta_data_path)
meta <- meta[rowSums(is.na(meta))<ncol(meta),colSums(is.na(meta))<nrow(meta)] # Remove rows and columns that are all NAs

sample_names <- unique(seurat_object@meta.data$sample_name)

if((!setequal(meta$sample_name,sample_names)) | (length(sample_names) != length(meta$sample_name))){
  stop("Check metadata.csv sample_name matches file names of objects")
}

meta <- meta[match(sample_names, meta$sample_name),] # Puts in the same order as SO
cell_value <- as.vector(table(seurat_object@meta.data$sample_name)) # Number of cells of each sample.

for(column in 2:length(colnames(meta))){
  column_name <- colnames(meta)[column]
  col <- meta[[column]]
  big_vec <- rep(col,cell_value)
  seurat_object@meta.data[[column_name]] <- big_vec
}

# Get total number of clusters
num_clusters <- nlevels(seurat_object$seurat_clusters)

# Get minimum between number of clusters selected and total number of clusters
if(max_clusters){
	num_clusters_to_use <- min(c(max_clusters, final_cluster))
	final_cluster <- num_clusters_to_use - 1
} else{
	final_cluster <- num_clusters - 1
}

# Colour vector
n <- length(unique(meta[[group_var]]))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col, n)

############################################
# Identify conserved cell type markers
############################################

DefaultAssay(seurat_object) <- "RNA"

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

  conserved_markers_tib <- as_tibble(conserved_markers)
  combined <- plyr::rbind.fill(combined, conserved_markers_tib)


}

combined_topmarkers <- unique(combined_topmarkers)
# Dotplot showing 2 top strong marker genes for each cluster
dotplot <- DotPlot(seurat_object, features = rev(combined_topmarkers), cols = col_vector, dot.scale = 8,
    split.by = group_var) + RotatedAxis()  # Edit so number of colours reflect number of conditions in group being split by
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


if(!group_var %in% colnames(seurat_object@meta.data)){
  stop("Grouping variable not listed, check spelling against colnames of meta data file")
}

for(comparison in de_conditions){
  conditions <- str_split(comparison, "_v_")[[1]]
  condition1 <- conditions[1]
  condition2 <- conditions[2]

  de_markers <- FindMarkers(seurat_object, ident.1 = condition1, ident.2 = condition2,
                            group.by = group_var) # Overall differences

  de_markers$ensembl <- rownames(de_markers)
  name_file <- paste0("Annotation_stats.dir/DifferentialMarkersAll_", sample_name, "_", condition1, "_vs_",condition2, ".csv")
  de_markers <- de_markers %>% dplyr::select(ensembl, everything())
  write_csv(de_markers, name_file)

  combined <- c()
  combined_topmarkers <- c()

  # Not every cluster is in every group. Subset for clusters that are in group
  subset_so1 <- (seurat_object@meta.data %>% dplyr::filter(!as.symbol(group_var) == condition1))$seurat_clusters %>%
    unique() %>% as.character() %>% as.numeric() %>% sort()
  subset_so2 <- (seurat_object@meta.data %>% dplyr::filter(!as.symbol(group_var) == condition2))$seurat_clusters %>%
    unique() %>% as.character() %>% as.numeric() %>% sort()
  clusters_keep <- intersect(subset_so1, subset_so2)

  for(i in clusters_keep){

    de_markers_cluster <- FindMarkers(seurat_object, ident.1 = condition1, ident.2 = condition2,
                            group.by = group_var, subset.ident = i)

    de_markers_cluster$ensembl <- rownames(de_markers_cluster)
    de_markers_cluster$cluster <- i

    de_markers_tib <- as_tibble(de_markers_cluster)
    combined <- plyr::rbind.fill(combined, de_markers_tib)
  }

  name_file <- paste0("Annotation_stats.dir/DifferentialMarkersPerCluster_", sample_name, "_", condition1, "_vs_",condition2, ".csv")

  combined <- combined %>% dplyr::select(ensembl, cluster, everything())
  write_csv(combined, name_file)

}

DefaultAssay(seurat_object) <- "RNA"
# Pre-defined list of markers
predefined <- opt$predefined
if(!is.null(predefined)){
  predefined_list <- str_split(predefined, "_")[[1]]
  number_markers <- length(predefined_list)

  for(marker in predefined_list){

    vln_plt <- VlnPlot(seurat_object, features = marker,  group.by = "seurat_clusters",
    pt.size = 0)

    ft_plt <- FeaturePlot(seurat_object, features = marker, label=TRUE)

    name<- paste0("Annotation_Figures.dir/Clusters_ViolinPlot_", sample_name, "_", marker ,".eps")
    postscript(name)
    print(vln_plt)
    dev.off()

    name<- paste0("Annotation_Figures.dir/Clusters_FeaturePlot_", sample_name, "_", marker ,".eps")
    postscript(name)
    print(ft_plt)
    dev.off()

  }

}