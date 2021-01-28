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

if(de_wilcoxin){
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=labels, de.method="wilcox", method = method)
}else{
  pred <- SingleR(test=seurat_object, assay.type.test=1, ref=ref_sce, labels=labels, method = method)
}
