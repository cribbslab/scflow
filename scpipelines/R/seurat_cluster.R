library(tidyverse)
library(Seurat)
library(optparse)


option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds file and path, filtered Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option(c("-v", "--variableFeatures"), default=2000, type="integer",
			help="Number of variable features  [default %default]"),
		make_option(c("--reddim"), default="pca",
			help="What dimensionality reduction to use. PCA, harmony, zinbwave, CCA etc. [default %default]"),
		  make_option(c("-d","--numdim"), default=NULL,
              help="how many dimensions to use. Default to all calculated ones"),
		  make_option(c("--resolution"), default=0.5,
              help="Resolution for finding clusters. Sets the 'granularity' of the downstream clustering,
              with increased values leading to a greater number of clusters. Between 0.4-1.2 good for ~3K cells.
              [default %default]")
)



opt <- parse_args(OptionParser(option_list=option_list))


sample_name <- opt$sample
num_variable_features <- opt$variableFeatures
num_dimensions <- opt$numdim
reduction_technique <- opt$reddim
resolution <- opt$resolution


seurat_filtered_rds <- opt$input
seurat_unfiltered_rds <- stringr::str_replace(opt$input, "filtered_SeuratObject.rds$", "unfiltered_SeuratObject.rds")

# Read in RDS files (may take some time)
filtered_seurat_object <- readRDS(seurat_filtered_rds)
unfiltered_seurat_object <- readRDS(seurat_unfiltered_rds)

# Normalise the data
filtered_seurat_object <- NormalizeData(filtered_seurat_object)

# Variable features (default = top 2000 genes)
filtered_seurat_object <- FindVariableFeatures(filtered_seurat_object, selection.method = "vst", nfeatures = num_variable_features)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(filtered_seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(filtered_seurat_object)
labelled_variable_feature_plot <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scale the data- linear transformation prior to dimension reduction
# (Might take some time)
all.genes <- rownames(filtered_seurat_object)
filtered_seurat_object <- ScaleData(filtered_seurat_object, features = all.genes)

# Principal component analysis (PCA)
filtered_seurat_object <- RunPCA(filtered_seurat_object, features = VariableFeatures(object = filtered_seurat_object))


# How much to reduce
filtered_seurat_object <- JackStraw(filtered_seurat_object, num.replicate = 100)
filtered_seurat_object <- ScoreJackStraw(filtered_seurat_object, dims = 1:20)


name<- paste0("Clustering_Figures.dir/JackStrawPlot_", sample_name, ".eps")
postscript(name)
print(JackStrawPlot(filtered_seurat_object, dims = 1:15))
dev.off()


name<- paste0("Clustering_Figures.dir/ElbowPlot_", sample_name, ".eps")
postscript(name)
print(ElbowPlot(filtered_seurat_object))
dev.off()

# Number of dimensions to reduce to. Look at Jackdraw / elbow plots, if not user-defined runs Embeddings
if(is.null(num_dimensions)){
  reddim<-Embeddings(filtered_seurat_object, reduction.type = reduction_technique)
  ncol(reddim)->maxN
  dims.use <- seq(1,maxN,1)
} else{
  dims.use <- seq(1,num_dimensions,1) #change this if you want to use less components, or use the significant ones only (jackstraw procedure on PCA)
}

filtered_seurat_object <- FindNeighbors(filtered_seurat_object, dims = dims.use, reduction = reduction_technique)
filtered_seurat_object <- FindClusters(filtered_seurat_object, resolution = resolution)

# Non-linear cluster techniques, e.g. UMAP and tSNE.
# Uses PCA to reduce dimensions.

filtered_seurat_object <- RunUMAP(filtered_seurat_object, dims = dims.use)

name<- paste0("Clustering_Figures.dir/UMAP_", sample_name, ".eps")
postscript(name)
print(DimPlot(filtered_seurat_object, reduction="umap", pt.size = 0.5))
dev.off()

filtered_seurat_object <- RunTSNE(filtered_seurat_object, dims = dims.use)

name<- paste0("Clustering_Figures.dir/tSNE_", sample_name, ".eps")
postscript(name)
print(DimPlot(filtered_seurat_object, reduction="tsne", pt.size = 0.5))
dev.off()

saveRDS(filtered_seurat_object, gsub("SAMPLE_FILE",sample_name ,"RDS_objects.dir/SAMPLE_FILE_clustered_filtered_SeuratObject.rds"))

