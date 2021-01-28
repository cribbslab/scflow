
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(sctransform)
library(patchwork)

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), default=NULL,
              help="The input rds file and path, filtered Seurat object"),
  make_option(c("-v", "--variableFeatures"), default=3000, type="integer",
              help="Number of variable features  [default %default]"),
  make_option(c("--resolution"), default=0.5,
              help="Resolution for finding clusters. Sets the 'granularity' of the downstream clustering,
              with increased values leading to a greater number of clusters. Between 0.4-1.2 good for ~3K cells.
              [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))


num_variable_features <- opt$variableFeatures
resolution <- opt$resolution
seurat_filtered_rds <- opt$input

# Read in RDS files (may take some time)
filtered_seurat_object <- readRDS(seurat_filtered_rds)

# Normalise the data by SCTransform

data.list <- mget(ls(pattern = "filtered_SeuratObject"))

for (i in 1:length(data.list)) {
  data.list[[i]] <- PercentageFeatureSet(data.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  data.list[[i]] <- SCTransform(data.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

# Integration by Seurat

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = num_variable_features)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = FALSE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)


# Principal component analysis (PCA)
data.integrated <- RunPCA(object = data.integrated, npcs = 30, verbose = FALSE)

# Clustering and data visualization in UMAP and tSNE - Seurat


data.integrated <- FindNeighbors(data.integrated, dims = 1:30, reduction = "pca")
data.integrated <- FindClusters(data.integrated, resolution = resolution)

# Non-linear cluster techniques, e.g. UMAP and tSNE.
# Uses PCA to reduce dimensions.
# UMAP

data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = 1:30)

name<- paste0("Clustering_Figures.dir/UMAP_Seurat_Integration_allsamples.eps")
postscript(name)
print(DimPlot(data.integrated, reduction="umap", pt.size = 0.5, label = TRUE))
dev.off()

name<- paste0("Clustering_Figures.dir/UMAP_Seurat_Integration_persample.eps")
postscript(name)
print(DimPlot(data.integrated, reduction="umap", split.by = "sample", ncol = 2, label = TRUE, pt.size = 0.5))
dev.off()

# tSNE

data.integrated <- RunTSNE(object = data.integrated, reduction = "pca", dims =  1:30)

name<- paste0("Clustering_Figures.dir/tSNE_Seurat_Integration_allsamples.eps")
postscript(name)
print(DimPlot(data.integrated, reduction="tsne", pt.size = 0.5, label = TRUE))
dev.off()

name<- paste0("Clustering_Figures.dir/tSNE_Seurat_Integration_persample.eps")
postscript(name)
print(DimPlot(data.integrated, reduction="tsne", split.by = "sample", ncol = 2, label = TRUE, pt.size = 0.5))
dev.off()

saveRDS(data.integrated, file= "RDS_objects.dir/SCT_integrated_SeuratObject.rds")

