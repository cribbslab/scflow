
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(sctransform)
library(harmony)
library(patchwork)
library(cowplot)
library(yaml)

sample_files <- str_replace(Sys.glob("RDS_objects.dir/*_filtered_clustered_SeuratObject.rds"), "RDS_objects.dir/", "")

ini <- read_yaml("pipeline.yml")

num_variable_features <- ini$num_variable_features
resolution <- ini$resolution
num_dimensions <- ini$num_dimensions
group_by_vars <- ini$group_vars
plot_convergence <- ini$plot_convergence
assay_use <- ini$assay_use

# Number of dimensions to reduce

dims.use <- seq(1,num_dimensions,1)

for (i in sample_files){
  name <- paste0("RDS_objects.dir/", i)
  so <- readRDS(name)
  so@meta.data$sample_name <- str_replace(i, "_filtered_clustered_SeuratObject.rds", "")
  assign(paste("so", i, sep = "."), so)
}

# Data processing and integration by Harmony
# Uses harmony to reduce dimensions.

data.list <- mget(ls(pattern = "filtered_clustered_SeuratObject"))

data.integrated <- Reduce(function(x, y){merge(x,y)}, data.list)

data.integrated <- data.integrated %>% 
                FindVariableFeatures(selection.method = "vst", nfeatures = num_variable_features) %>% 
                ScaleData(verbose=FALSE) %>% 
                RunPCA(npcs = num_dimensions, verbose=FALSE)



# Data visualization of harmony reduction (DimPlot) and VlnPlot - Harmony
options(repr.plot.height = 5, repr.plot.width = 20)
p1 <- DimPlot(object = data.integrated, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data.integrated, features = "PC_1", group.by = "sample", pt.size = .1)
plt <- plot_grid(p1,p2)
print(plt)
pdf("Integration_Figures.dir/DimPlot_Harmony_Integration_Dim_vs_Vlnplots_before_harmony.eps")
print(plt)
dev.off()


data.integrated <- RunHarmony(data.integrated, group.by.vars = "sample_name", plot_convergence = plot_convergence)

data.integrated <- data.integrated %>% 
  RunUMAP(reduction = "harmony", dims = dims.use) %>% 
  RunTSNE(reduction = "harmony", dims = dims.use) %>% 
  FindNeighbors(reduction = "harmony", dims = dims.use) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()


# Data visualization of harmony reduction (DimPlot) and VlnPlot - Harmony

options(repr.plot.height = 5, repr.plot.width = 20)
p1 <- DimPlot(object = data.integrated, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data.integrated, features = "harmony_1", group.by = "sample", pt.size = .1)
plt <- plot_grid(p1,p2)
print(plt)
pdf("Integration_Figures.dir/DimPlot_Harmony_Integration_Dim_vs_Vlnplots.eps")
print(plt)
dev.off()

# Clustering and data visualization by Harmony- UMAP and tSNE

options(repr.plot.height = 4, repr.plot.width = 6)

# UMAP

pdf("Integration_Figures.dir/UMAP_harmony_Integration_all_samples.eps")
print(DimPlot(data.integrated, reduction = "umap", label = TRUE, pt.size = .1))
dev.off()


pdf("Integration_Figures.dir/UMAP_harmony_Integration_all_samples_colour.eps")
print(DimPlot(data.integrated, reduction = "umap", group.by="sample", label = TRUE, pt.size = .1))
dev.off()


pdf("Integration_Figures.dir/UMAP_harmony_Integration_per_samples.eps")
print(DimPlot(data.integrated, reduction = "umap", split.by = "sample", ncol = 4, label = TRUE), width=15, height=15)
dev.off()


# tSNE

pdf("Integration_Figures.dir/tSNE_harmony_Integration_all_samples.eps")
print(DimPlot(data.integrated, reduction = "tsne", label = TRUE, pt.size = .1))
dev.off()

pdf("Integration_Figures.dir/tSNE_harmony_Integration_allsamples_colour.eps")
print(DimPlot(data.integrated, reduction = "tsne", group.by="sample" ,label = TRUE, pt.size = .1))
dev.off()

pdf("Integration_Figures.dir/tSNE_harmony_Integration_per_samples.eps")
print(DimPlot(data.integrated, reduction = "tsne", split.by = "sample", ncol = 4, label = TRUE), width=15, height=15)
dev.off()

saveRDS(data.integrated, file="RDS_objects.dir/Harmony_integrated_SeuratObject.rds")
