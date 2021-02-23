library(Seurat)
library(tidyverse)
library(optparse)
library(scRNAseq)
library(celldex)
library(SingleR)
library(pheatmap)
#library(scuttle)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd integrated Seurat object"),
		make_option(c("-o", "--output"), default=NULL,
			help="The output rds path."),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file"),
        make_option(c("-d", "--DEmethod"), default=0, type="integer",
            help="Whether to apply Wilcoxin pairwise comparison tests between labels. [default %default, i.e. no DE test]"),
        make_option(c("-m", "--method"), default="single", type="character",
            help="SingleR method, annotating by single cells or by cluster. [default %default]"),
        make_option(c("--predefined"), default=NULL,
			help="User pre-defined list of marker genes")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
output_file <- opt$output
sample_name <- opt$sample
reference_sce_loc <- opt$reference

de_wilcoxin <- opt$DEmethod
method <- opt$method

# Load seurat object
seurat_object <- readRDS(input_file)
# Convert to single cell object
sce <- as.SingleCellExperiment(seurat_object)


# Load in reference sce object
ref_sce <- readRDS(reference_sce_loc)

if(method == "cluster"){
  # Character vector of each cluster for each cell in object
  clusters <-  as.vector(seurat_object$seurat_clusters)
}else{
  clusters <- NULL
}

if(de_wilcoxin){
  pred <- SingleR(test=sce, assay.type.ref="logcounts", assay.type.test = "logcounts", ref=ref_sce, labels=ref_sce$label, de.method="wilcox", clusters=clusters)
}else{
  pred <- SingleR(test=sce, assay.type.ref="logcounts", assay.type.test="logcounts", ref=ref_sce, labels=ref_sce$label, clusters=clusters)
}

seurat_object@meta.data[['singleR_labels']] <- pred$labels
seurat_object@meta.data[['singleR_pruned.labels']] <- pred$pruned.labels

saveRDS(seurat_object, output_file) # Save seurat object with labels


# Plots
name<- paste0("Annotation_Figures.dir/singleR_scoreHeatmap_", sample_name, ".eps")
postscript(name)
print(plotScoreHeatmap(pred))
dev.off()

name<- paste0("Annotation_Figures.dir/singleR_scoreDistribution_", sample_name, ".eps")
postscript(name)
print(plotScoreDistribution(pred))
dev.off()

name<- paste0("Annotation_Figures.dir/singleR_deltaDistribution_", sample_name, ".eps")
postscript(name)
print(plotDeltaDistribution(pred, ncol = 3))
dev.off()

DefaultAssay(seurat_object) <- "RNA"

# Labelled UMAP and tSNE dimension plots with singleR labels
umap_plt <- DimPlot(seurat_object, group.by = "singleR_labels", reduction = "umap")
umap_labelled <- DimPlot(seurat_object, group.by = "singleR_labels", label = TRUE, reduction = "umap")

tsne_plt <- DimPlot(seurat_object, group.by = "singleR_labels", reduction = "tsne")
tsne_labelled <- DimPlot(seurat_object, group.by = "singleR_labels", label = TRUE, reduction = "tsne")

name<- paste0("Annotation_Figures.dir/singleR_UMAP_", sample_name, ".eps")
postscript(name)
print(umap_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/singleR_UMAP_labelled_", sample_name, ".eps")
postscript(name)
print(umap_labelled)
dev.off()

name<- paste0("Annotation_Figures.dir/singleR_tSNE_", sample_name, ".eps")
postscript(name)
print(tsne_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/singleR_tSNE_labelled_", sample_name, ".eps")
postscript(name)
print(tsne_labelled)
dev.off()

# Pre-defined list of markers
predefined <- opt$predefined
if(!is.null(predefined)){
  predefined_list <- str_split(predefined, "_")[[1]]
  number_markers <- length(predefined_list)

  for(marker in predefined_list){

    vln_plt <- VlnPlot(seurat_object, features = marker,  group.by = "singleR_labels",
    pt.size = 0)

    name<- paste0("Annotation_Figures.dir/singleR_ViolinPlot_", sample_name, "_", marker ,".eps")
    postscript(name)
    print(vln_plt)
    dev.off()
  }
}