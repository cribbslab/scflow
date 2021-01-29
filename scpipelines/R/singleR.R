library(Seurat)
library(tidyverse)
library(optparse)
library(scRNAseq)
library(celldex)
library(SingleR)
library(scuttle)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd integrated Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file"),
        make_option(c("-d", "--DEmethod"), default=0, type="integer",
            help="Whether to apply Wilcoxin pairwise comparison tests between labels. [default %default, i.e. no DE test]"),
        make_option(c("-m", "--method"), default=0, type="single",
            help="SingleR method, annotating by single cells or by cluster. [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
sample_name <- opt$sample
reference_sce_loc <- opt$reference

de_wilcoxin <- opt$DEmethod
method <- opt$method

# Load seurat object
seurat_object <- readRDS(input_file)
# Might need to convert to sce object instead
# Check if there are logcounts already

# Load in reference sce object
ref_sce <- readRDS("reference_sce_loc")

if(method == "cluster"){
  # Character vector of each cluster for each cell in object
  clusters <-  as.vector(seurat_object$seurat_clusters)
}else{
  clusters <- NULL
}

if(de_wilcoxin){
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=ref_sce$label, de.method="wilcox", clusters=clusters)
}else{
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=ref_sce$label, clusters=clusters)
}

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

# To combine annotations with tSNE and UMAP.
# ...