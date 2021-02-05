library(Seurat)
library(optparse)
library(scRNAseq)
library(celldex)

option_list <- list(
        make_option(c("-r", "--reference"), default="celldex", type = "character",
			help="Which reference package to use, celldex or scRNAseq package"),
		make_option(c("-c", "--referenceNameCellDex"), default=NULL, type="character",
			help="Reference dataset from celldex. E.g. HumanPrimaryCellAtlasData [default %default]."),
        make_option(c("-s", "--referenceNameSC"), default=NULL, type="character",
            help="Reference dataset from scRNAseq package. E.g. KotliarovPBMCData [default %default]. Find datasets here:
			 https://bioconductor.org/packages/devel/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html#available-data-sets"),
        make_option(c("-o", "--scRNAseqOption"), default=NULL, type="character",
            help="Any options for scRNAseq reference, e.g. 'human', [default %default]"),
        make_option(c("--outfile"), default="reference_sce.rds", type="character",
            help="Output rds file, [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

reference_package <- opt$reference
celldex_reference_name <- opt$referenceNameCellDex
scRNAseq_reference_name <- opt$referenceNameSC
scRNAseqOption <- opt$scRNAseqOption
outfile <- opt$outfile


# Cell-dex reference
if((reference_package == "celldex") & (!is.null(celldex_reference_name))){
  ref_sce <- do.call(celldex_reference_name, list(ensembl=TRUE))
  labels <- ref_sce$label.main
  ref_sce$label <- ref_sce$label.main
}

# scRNAseq package reference call
if((reference_package == "scRNAseq") & (!is.null(scRNAseq_reference_name))){
  # if there are options in the function, e.g. 'human' or 'mouse'
  if(!is.null(scRNAseqOption)){
    ref_sce <- do.call(scRNAseq_reference_name, list(scRNAseqOption, ensembl=TRUE))
  }else{
    ref_sce <- do.call(scRNAseq_reference_name, list(ensembl=TRUE))
  }

  # Get rid of unlabbelled cells
  ref_sce <- ref_sce[,!is.na(ref_sce$label)]
  ref_sce <- logNormCounts(ref_sce) # Do I need this??
}

# Save reference matrix to RDS file in main directory
saveRDS(ref_sce, outfile)