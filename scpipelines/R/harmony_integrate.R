
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(sctransform)
library(harmony)
library(patchwork)
library(cowplot)
library(yaml)


ini <- read_yaml("pipeline.yml")

num_variable_features <- ini$num_variable_features
resolution <- ini$resolution
num_dimensions <- ini$num_dimensions
group_by_vars <- ini$group_vars
plot_convergence <- ini$plot_convergence
assay_use <- ini$assay_use

# Number of dimensions to reduce

dims.use <- seq(1,num_dimensions,1)


# Data processing and integration by Harmony
# Uses harmony to reduce dimensions.
data.integrated <- readRDS("RDS_objects.dir/SCT_integrated_SeuratObject.rds")


data.integrated <- RunHarmony(data.integrated, group.by.vars = group_by_vars, plot_convergence = plot_convergence, assay.use = assay_use)

data_harmony <- data.integrated %>% 
  RunUMAP(reduction = "harmony", dims = dims.use) %>% 
  RunTSNE(reduction = "harmony", dims = dims.use) %>% 
  FindNeighbors(reduction = "harmony", dims = dims.use) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()


# Data visualization of harmony reduction (DimPlot) and VlnPlot - Harmony

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data_harmony, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data_harmony, features = "harmony_1", group.by = "sample", pt.size = .1)
plt <- plot_grid(p1,p2)
print(plt)
name<- paste0("Integration_Figures.dir/DimPlot_Harmony_Integration_Dim_vs_Vlnplots.eps")
postscript(name)
print(plt)
dev.off()

# Clustering and data visualization by Harmony- UMAP and tSNE

options(repr.plot.height = 4, repr.plot.width = 6)

# UMAP

name<- paste0("Integration_Figures.dir/UMAP_harmony_Integration_all_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "umap", label = TRUE, pt.size = .1))
dev.off()

name<- paste0("Integration_Figures.dir/UMAP_harmony_Integration_per_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "umap", split.by = "sample", ncol = 2, label = TRUE))
dev.off()


# tSNE

name<- paste0("Integration_Figures.dir/tSNE_harmony_Integration_all_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "tsne", label = TRUE, pt.size = .1))
dev.off()

name<- paste0("Integration_Figures.dir/tSNE_harmony_Integration_per_samples.eps")
postscript(name)
print(DimPlot(data_harmony, reduction = "tsne", split.by = "sample", ncol = 2, label = TRUE))
dev.off()

saveRDS(data_harmony, file="RDS_objects.dir/Harmony_integrated_SeuratObject.rds")

