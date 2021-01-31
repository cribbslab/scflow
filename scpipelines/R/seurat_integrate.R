
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(sctransform)
library(patchwork)
library(yaml)

sample_files <- str_replace(Sys.glob("RDS_objects.dir/*_filtered_SeuratObject.rds"), "RDS_objects.dir/", "")

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
  assign(paste("so", i, sep = "."), so)
}



# Normalise the data by SCTransform

data.list <- mget(ls(pattern = "filtered_SeuratObject"))

for (i in 1:length(data.list)) {
  data.list[[i]] <- PercentageFeatureSet(data.list[[i]], pattern = pattern, col.name = colname)
  data.list[[i]] <- SCTransform(data.list[[i]], vars.to.regress = vars_to_regress, verbose = FALSE)
}

# Integration by Seurat
options(future.globals.maxSize = 20000 *1024^2)

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = num_variable_features)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = FALSE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)


# Principal component analysis (PCA)
data.integrated <- RunPCA(object = data.integrated, npcs = 30, verbose = FALSE)


# Number of dimensions to reduce

dims.use <- seq(1,num_dimensions,1)

# Clustering and data visualization in UMAP and tSNE - Seurat

data.integrated <- FindNeighbors(data.integrated, dims = dims.use, reduction = "pca")
data.integrated <- FindClusters(data.integrated, resolution = resolution)

# Non-linear cluster techniques, e.g. UMAP and tSNE.
# Uses PCA to reduce dimensions.
# UMAP

data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = dims.use)

pdf("Integration_Figures.dir/UMAP_Seurat_Integration_allsamples.eps")
print(DimPlot(data.integrated, reduction="umap", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/UMAP_Seurat_Integration_persample.eps")
print(DimPlot(data.integrated, reduction="umap", split.by = "sample", ncol = 4, label = TRUE, pt.size = 0.5), width=15, height=15)
dev.off()

# tSNE

data.integrated <- RunTSNE(object = data.integrated, reduction = "pca", dims = dims.use)

pdf("Integration_Figures.dir/tSNE_Seurat_Integration_allsamples.eps")
print(DimPlot(data.integrated, reduction="tsne", pt.size = 0.5, label = TRUE))
dev.off()

pdf("Integration_Figures.dir/tSNE_Seurat_Integration_persample.eps")
print(DimPlot(data.integrated, reduction="tsne", split.by = "sample", ncol = 4, label = TRUE, pt.size = 0.5), width=15, height=15)
dev.off()

saveRDS(data.integrated, file= "RDS_objects.dir/SCT_integrated_SeuratObject.rds")

