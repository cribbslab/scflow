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
              [default %default]"),
		make_option(c("-m", "--metadata"), default="none",
			help="The input metadata file from doublet pipeline")
)


opt <- parse_args(OptionParser(option_list=option_list))


sample_name <- opt$sample
num_variable_features <- opt$variableFeatures
num_dimensions <- opt$numdim
reduction_technique <- opt$reddim
resolution <- opt$resolution

if(opt$metadata != "none"){
   metadata <- read.csv(opt$metadata, row.names=1)
}
seurat_filtered_rds <- opt$input

# Read in RDS files (may take some time)
filtered_seurat_object <- readRDS(seurat_filtered_rds)

if(opt$metadata != "none"){
# add new metadata to so
print(head(opt$metadata))
filtered_seurat_object@meta.data <- metadata

# Filter according to metadata doublet class
filtered_seurat_object <- subset(filtered_seurat_object, subset = scDblFinder.class == "singlet")
}

# Normalise the data
filtered_seurat_object <- NormalizeData(filtered_seurat_object)

# Variable features (default = top 2000 genes)
filtered_seurat_object <- FindVariableFeatures(filtered_seurat_object, selection.method = "vst", nfeatures = num_variable_features)

# Scale the data- linear transformation prior to dimension reduction
# (Might take some time)
all.genes <- rownames(filtered_seurat_object)
filtered_seurat_object <- ScaleData(filtered_seurat_object, features = all.genes)

# Principal component analysis (PCA)
filtered_seurat_object <- RunPCA(filtered_seurat_object, features = VariableFeatures(object = filtered_seurat_object))

# Number of dimensions to reduce to. Look at Jackdraw / elbow plots in Cluster.Rmd, if not user-defined runs Embeddings
if(is.null(num_dimensions)){
  red_dim<-Embeddings(filtered_seurat_object, reduction = reduction_technique)
  ncol(red_dim)->maxN
  dims.use <- seq(1,maxN,1)
} else{
  dims.use <- seq(1,num_dimensions,1) #change this if you want to use less components, or use the significant ones only (jackstraw procedure on PCA)
}

filtered_seurat_object <- FindNeighbors(filtered_seurat_object, dims = dims.use, reduction = reduction_technique)
filtered_seurat_object <- FindClusters(filtered_seurat_object, resolution = resolution)

# Non-linear cluster techniques, e.g. UMAP and tSNE.
# Uses PCA to reduce dimensions.

filtered_seurat_object <- RunUMAP(filtered_seurat_object, dims = dims.use)
filtered_seurat_object <- RunTSNE(filtered_seurat_object, dims = dims.use)

# save RDS object with clustering
saveRDS(filtered_seurat_object, gsub("SAMPLE_FILE",sample_name ,"RDS_objects.dir/SAMPLE_FILE_filtered_clustered_SeuratObject.rds"))

