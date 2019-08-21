library(Seurat)
library(velocyto.R)
library(ggplot2)
library(tidyverse)
library("optparse")


option_list = list(
  make_option(c("-l", "--loom"), type="character",
              help="Name of loom file", metavar="character"),
  make_option(c("-d", "--dimensions"), type="integer", default = 6,
              help="Number of dimensions for clustering, [default= %default]. Look ar JackStrawPlot and rerun.", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

loom_file <- opt$loom
max_dimensions <- opt$dimensions

sample_name <- loom_file %>% str_remove(".loom")

ldat <- read.loom.matrices(loom_file)

# Percentage of spliced vs unspliced vs ambiguous
type_vec <- c("spliced", "ambiguous", "unspliced")
splicing <- tibble(type = type_vec, total = c(sum(ldat$spliced),sum(ldat$ambiguous), sum(ldat$unspliced)))
splicing$fraction <- splicing$total / sum(splicing$total)
splicing$type <- factor(type_vec, levels = c("spliced", "ambiguous", "unspliced"))

png(paste(c(sample_name, "_plot_fractions.png"), collapse = ""))
ggplot(data= splicing, aes(x= type, y = fraction)) + geom_bar(stat="identity")
dev.off()

emat <- ldat$spliced
emat <- emat[,colSums(emat)>=1e3]

seurat <- CreateSeuratObject(counts = emat, names.delim = ":")

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

png(paste(c(sample_name, "_violin_plots.png"), collapse = ""))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

counts_features_plot <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(counts_features_plot)

# Subset based on
seurat <- subset(seurat, subset = nFeature_RNA > 100 )

# Normalise
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
print(plot1)

top20 <- head(VariableFeatures(seurat), 20)
png(paste(c(sample_name, "_variable_genes.png"), collapse = ""))
LabelPoints(plot = plot1, points = top20, repel = TRUE)
dev.off()

all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

# Run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

DimPlot(seurat, reduction = "pca")

DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 2, cells = 500, balanced = TRUE)

# Determine dimensionality

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:20)

png(paste(c(sample_name, "_jack_straw_plot.png"), collapse = ""))
JackStrawPlot(seurat, dims = 1:15)
dev.off()

png(paste(c(sample_name, "_elbow_plot.png"), collapse = ""))
ElbowPlot(seurat)
dev.off()

# Clustering 
seurat <- FindNeighbors(seurat, dims = 1:max_dimensions) # Using PCs 1-10
seurat <- FindClusters(seurat, resolution = 0.5)

# UMAP
seurat <- RunUMAP(seurat, dims = 1:max_dimensions)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
png(paste(c(sample_name, "_UMAP.png"), collapse = ""))
DimPlot(seurat, reduction = "umap")
dev.off()

# TSNE
seurat <- RunTSNE(seurat, dims =1:max_dimensions)

png(paste(c(sample_name, "_TSNE.png"), collapse = ""))
DimPlot(seurat, reduction = "tsne")
dev.off()

# Back to velocyto
emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,colnames(seurat)]
nmat <- nmat[,colnames(seurat)]

# Embeddings (tsne or UMAP?)
emb <- seurat@reductions$tsne@cell.embeddings

# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(emb)))

# Not sure of this variable
fit.quantile <- 0.02

# Main velocity estimation
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=24)

png(paste(c(sample_name, "_cell_velocity.png"), collapse = ""))
gg <- DimPlot(seurat, reduction = "tsne")

colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)


show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                      cell.colors=ac(colors,alpha=0.5),
                                      cex=0.8,arrow.scale=2,show.grid.flow=T,
                                      min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                      do.par=F,cell.border.alpha = 0.1,
                                      n.cores=24,main="Cell Velocity")

dev.off()
