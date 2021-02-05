library(Seurat)
library(tidyverse)
library(optparse)
library(scRNAseq)
library(celldex)
library("scClassify")


option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd integrated Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file"),
        make_option(c("-p", "--pretrained"), default=FALSE,
            help="Whether to use a pretrained model, default=0, or give path"),
        make_option(c("-m", "--method"), default="predict", type = "character",
            help="scClassify method to use. predict, ensemble, nonensemble. [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
sample_name <- opt$sample
reference_sce_loc <- opt$reference
pretrain <- opt$pretrained
method <- opt$method


# Load seurat object
seurat_object <- readRDS(input_file)
# Sparse matrix
dgc_mat <- as(object_data(seurat_object, "data"), "dgCMatrix")

if(pretrain==FALSE){
  # Train new model - if statement needed and parameter
  ref_sce <- readRDS(reference_sce_loc) # Load in reference sce object
  ref_mat <- as(assay(ref_sce), "dgCMatrix")
  cellTypes_train <- ref_sce$label

  if(method == "predict"){
  # Train model using reference
    trainClass <- train_scClassify(exprsMat_train = ref_mat,
                                 cellTypes_train = cellTypes_train,
                                 selectFeatures = c("limma", "BI"),
                                 returnList = FALSE)
  }
}else{
  trainClass <- readRDS(pretrain)

}

# Pre trained model
if(method == "predict"){
  pred_res <- predict_scClassify(exprsMat_test = dgc_mat,
                                 trainRes = trainClass, # or pre-made reference using
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = c("pearson", "spearman"),
                                 prob_threshold = 0.7,
                                 verbose = TRUE)
}


# Ensemble classify
if(method == "ensemble"){
  scClassify_res_ensemble <- scClassify(exprsMat_train = ref_mat,
                                        cellTypes_train = cellTypes_train,
                                        exprsMat_test = dgc_mat,
                                        tree = "HOPACH",
                                        algorithm = "WKNN",
                                        selectFeatures = c("limma"),
                                        similarity = c("pearson", "spearman"), # Could do as list parameter with split??
                                        weighted_ensemble = TRUE, # TRUE or false parameterised ???
                                        returnList = FALSE,
                                        verbose = FALSE)
}

# Non-Ensemble classify
if(method == "ensemble"){
  scClassify_res <- scClassify(exprsMat_train = ref_mat,
                                        cellTypes_train = cellTypes_train,
                                        exprsMat_test = dgc_mat,
                                        tree = "HOPACH",
                                        algorithm = "WKNN",
                                        selectFeatures = c("limma"),
                                        similarity = c("pearson"),
                                        returnList = FALSE,
                                        verbose = FALSE)
}