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
        make_option(c("-r", "--reference"), default="celldex", type = "character"
			help="Which reference package to use, celldex or scRNAseq package"),
		make_option(c("-c", "--referenceNameCellDex"), default=NULL, type="character",
			help="Reference dataset from celldex. E.g. HumanPrimaryCellAtlasData [default %default]."),
        make_option(c("-s", "--referenceNameSC"), default=NULL, type="character",
            help="Reference dataset from scRNAseq package. E.g. KotliarovPBMCData [default %default]. Find datasets here:
			 https://bioconductor.org/packages/devel/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html#available-data-sets"),
        make_option(c("-o", "--scRNAseqOption"), default=NULL, type="character",
            help="Any options for scRNAseq reference, e.g. 'human', [default %default]"),
        make_option(c("-d", "--DEmethod"), default=0, type="integer",
            help="Whether to apply Wilcoxin pairwise comparison tests between labels. [default %default, i.e. no DE test]"),
        make_option(c("-m", "--method"), default=0, type="single",
            help="SingleR method, annotating by single cells or by cluster. [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
sample_name <- opt$sample
reference_package <- opt$reference
celldex_reference_name <- opt$referenceNameCellDex
scRNAseq_reference_name <- opt$referenceNameSC
scRNAseqOption <- opt$scRNAseqOption
de_wilcoxin <- opt$DEmethod
method <- opt$method

# Read in seurat object
seurat_object <- readRDS(input_file)
# Might need to convert to sce object instead
# Check if there are logcounts already

# Cell-dex reference call
if((reference_package == "celldex") & (!is.null(celldex_reference_name))){
  ref <- do.call(celldex_reference_name, list())
  labels <- ref$label.main
}

# scRNAseq package reference call
if((reference_package == "scRNAseq") & (!is.null(scRNAseq_reference_name))){
  # if there are options in the function, e.g. 'human' or 'mouse'
  if(!is.null(scRNAseqOption)){
    ref_sce <- do.call(scRNAseq_reference_name, list(scRNAseqOption))
  }else{
    ref_sce <- do.call(scRNAseq_reference_name, list())
  }

  # Get rid of unlabbelled cells
  ref_sce <- ref_sce[,!is.na(ref_sce$label)]
  ref_sce <- logNormCounts(ref_sce)
  labels <- ref_sce$label
}

if(method == "cluster"){
  # Character vector of each cluster for each cell in object
  clusters <-  as.vector(seurat_object$seurat_clusters)
}else{
  clusters <- NULL
}

if(de_wilcoxin){
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=labels, de.method="wilcox", clusters=clusters)
}else{
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=labels, clusters=clusters)
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