
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(harmony)
library(cowplot)

# Data processing and integration by Harmony
# Uses harmony to reduce dimensions.

data.integrated <- RunHarmony(data.integrated, group.by.vars = "orig.ident", plot_convergence = TRUE, assay.use = "SCT")

data_harmony <- data.integrated %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()


# Data visualization of harmony reduction (DimPlot) and VlnPlot - Harmony

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data_harmony, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data_harmony, features = "harmony_1", group.by = "sample", pt.size = .1)
plt <- plot_grid(p1,p2)
print(plt)
name<- paste0("Clustering_Figures.dir/DimPlot_Integration_Harmony_Dim_vs_Vlnplots.eps")
postscript(name)
print(plt)
dev.off()

# Clustering and data visualization by Harmony- UMAP and tSNE

options(repr.plot.height = 4, repr.plot.width = 6)

# UMAP

name<- paste0("Clustering_Figures.dir/UMAP_Integration_harmony_all_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "umap", label = TRUE, pt.size = .1))
dev.off()

name<- paste0("Clustering_Figures.dir/UMAP_Integration_harmony_per_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "umap", split.by = "sample", ncol = 2, label = TRUE))
dev.off()


# tSNE

name<- paste0("Clustering_Figures.dir/tSNE_Integration_harmony_all_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "tsne", label = TRUE, pt.size = .1))
dev.off()

name<- paste0("Clustering_Figures.dir/tSNE_Integration_harmony_per_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "tsne", split.by = "sample", ncol = 2, label = TRUE))
dev.off()

saveRDS(data_harmony, file="RDS_objects.dir/Harmony_integrated_SeuratObject.rds")

