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
 		make_option(c("-o", "--output"), default=NULL,
			help="Output file name"),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file"),
        make_option(c("--pretrained"), default=0,
            help="Whether to use a pretrained model, default=0, or give path"),
        make_option(c("-m", "--method"), default="predict", type = "character",
            help="scClassify method to use. predict, ensemble, nonensemble. [default %default]"),
        make_option(c("--similarity"), default="pearson_spearman", type = "character",
            help="Similarity test to use. If multiple, separated by space. [default %default]")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
output <- opt$output
sample_name <- opt$sample
reference_sce_loc <- opt$reference
method <- opt$method

pretrain <- opt$pretrained
if(pretrain == "0" || method != "predict" ){
  pretrain <- 0
}
similarity_string <- opt$similarity
similarity <- str_split(similarity_string, "_")[[1]]


# Load seurat object
seurat_object <- readRDS(input_file)
# Sparse matrix
dgc_mat <- seurat_object[["RNA"]]@counts

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
  if ("scClassifyTrainModelList" %in% is(trainClass)) { # If train class is a list of multiple models
    pred_res <- predict_scClassifyJoint(exprsMat_test = dgc_mat,
                                 trainRes = trainClass, # or pre-made reference using
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = similarity,
                                 prob_threshold = 0.7,
                                 verbose = FALSE)
  }else{
    pred_res <- predict_scClassify(exprsMat_test = dgc_mat,
                                 trainRes = trainClass, # or pre-made reference using
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = similarity,
                                 prob_threshold = 0.7,
                                 verbose = FALSE)
  }

  pearson_cells_assigned <- as.vector(pred_res$pearson_WKNN_limma$predRes)
  spearmen_cells_assigned <- as.vector(pred_res$spearman_WKNN_limma$predRes)
  pred_results  <- pred_res$ensembleRes
  pred_results$pearson <- pearson_cells_assigned
  pred_results$spearman <- spearmen_cells_assigned
  pred_results$cell_barcode <- rownames(pred_results)

  #name_file <- paste0("scclassify_predict_", sample_name, ".rds")
  #saveRDS(pred_results, name_file)

  name_file <- paste0("Annotation_stats.dir/scclassify_predict_", sample_name, ".csv.gz")
  write_csv(pred_results,name_file)

  seurat_object@meta.data[['scclassify_labels']] <- pred_results$cellTypes

}


# Ensemble/non-ensembl classify. Train and test model in one
if(method != "predict"){

  if (any(table(cellTypes_train) == 1)) { # Sort out stop error thrown by scClassify which doesn't like cell types with 1 cell only
    tab <- table(cellTypes_train)
    good <- cellTypes_train %in% (tab[tab !=1] %>% rownames) # True and falses
    ref_mat <- ref_mat[,good]
    cellTypes_train <- cellTypes_train[good]
  }

  scClassify_res_ensemble <- scClassify(exprsMat_train = ref_mat,
                                        cellTypes_train = cellTypes_train,
                                        exprsMat_test = dgc_mat,
                                        tree = "HOPACH",
                                        algorithm = "WKNN",
                                        selectFeatures = c("limma"),
                                        similarity = similarity, # Could do as list parameter with split??
                                        weighted_ensemble = TRUE, # TRUE or false parameterised ???
                                        returnList = FALSE,
                                        verbose = FALSE)
  #name_file <- paste0("scclassify_ensembl_", sample_name, ".rds")
  #saveRDS(scClassify_res_ensemble,name_file)

  results <- scClassify_res_ensemble$testRes$test$ensembleRes
  for(metric in similarity){
    column_name <- paste0(metric, "_WKNN_limma")
    metric_cells_assigned <- as.vector(scClassify_res_ensemble$testRes$test[[column_name]]$predRes)
    results[[metric]] <- metric_cells_assigned
  }
  results$cell_barcode <- rownames(results)
  name_file <- paste0("Annotation_stats.dir/scclassify_train_", sample_name, ".csv.gz")
  write_csv(results,name_file)


  seurat_object@meta.data[['scclassify_labels']] <- results$cellTypes

}

saveRDS(seurat_object, output) # Save seurat object with labels