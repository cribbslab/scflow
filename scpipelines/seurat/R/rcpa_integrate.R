library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(sctransform)
library(patchwork)
library(yaml)

sample_files <- str_replace(Sys.glob("RDS_objects.dir/*_filtered_clustered_SeuratObject.rds"), "RDS_objects.dir/", "")

ini <- read_yaml("pipeline.yml")

pattern <- ini$pattern
colname <- ini$colname
vars_to_regress <- ini$vars_to_regress
num_variable_features <- ini$num_variable_features
resolution <- ini$resolution
num_dimensions <- ini$num_dimensions

for (i in sample_files){
  name <- paste0("RDS_objects.dir/", i)
  so <- readRDS(name)
  so@meta.data$sample_name <- str_replace(i, "_filtered_clustered_SeuratObject.rds", "")
  assign(paste("so", i, sep = "."), so)
}



# Normalise the data by SCTransform

data.list <- mget(ls(pattern = "filtered_clustered_SeuratObject"))

print(data.list)

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Integration by Seurat
options(future.globals.maxSize = 20000 *1024^2)

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = num_variable_features)

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- ScaleData(x, features = data.features, verbose = FALSE)
    x <- RunPCA(x, features = data.features, verbose = FALSE)
})

data.anchors <- FindIntegrationAnchors(object.list = data.list, reduction = "rpca", 
                                       anchor.features = data.features, verbose = FALSE)
data.integrated <- IntegrateData(anchorset = data.anchors, verbose = FALSE, k.weight = ini$kweight)


data.integrated <- ScaleData(data.integrated, verbose=FALSE)

# Principal component analysis (PCA)
data.integrated <- RunPCA(object = data.integrated, npcs = 30, verbose = FALSE)

# Non-linear cluster techniques, e.g. UMAP and tSNE.
# Uses PCA to reduce dimensions.
# UMAP
dims.use <- seq(1,num_dimensions,1)

data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = dims.use)


# Clustering and data visualization in UMAP and tSNE - Seurat

data.integrated <- FindNeighbors(data.integrated, dims = dims.use, reduction = "pca")
data.integrated <- FindClusters(data.integrated, resolution = resolution)


pdf("Integration_Figures.dir/UMAP_rcpa_Integration_allsamples.eps")
print(DimPlot(data.integrated, reduction="umap", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/UMAP_rcpa_Integration_allsamples_colour.eps")
print(DimPlot(data.integrated, reduction="umap", group.by="sample", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/UMAP_rcpa_Integration_persample.eps")
print(DimPlot(data.integrated, reduction="umap", split.by = "sample", ncol = 4, label = TRUE, pt.size = 0.5), width=15, height=15)
dev.off()

# tSNE

data.integrated <- RunTSNE(object = data.integrated, reduction = "pca", dims = dims.use)

pdf("Integration_Figures.dir/tSNE_rcpa_Integration_allsamples.eps")
print(DimPlot(data.integrated, reduction="tsne", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/tSNE_rcpa_Integration_persamples_colour.eps")
print(DimPlot(data.integrated, reduction="tsne", group.by="sample", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/tSNE_rcpa_Integration_persample.eps")
print(DimPlot(data.integrated, reduction="tsne", split.by = "sample", ncol = 4, label = TRUE, pt.size = 0.5), width=15, height=15)
dev.off()

saveRDS(data.integrated, file= "RDS_objects.dir/rcpa_integrated_SeuratObject.rds")

