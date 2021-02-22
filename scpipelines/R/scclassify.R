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
output_file <- opt$output
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
                                 selectFeatures = c("limma"),
                                 returnList = FALSE)
  }
}else{
  trainClass <- readRDS(pretrain)

}

# Pre trained model
if(method == "predict"){
  if ("scClassifyTrainModelList" %in% is(trainClass)) { # If train class is a list of multiple models
    similarity  <- similarity[1] # Error with scClassify where multiple metrics seem to make predict_scClassifyJoint fail
    pred_res <- predict_scClassifyJoint(exprsMat_test = dgc_mat,
                                 trainRes = trainClass, # or pre-made reference using
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = similarity,
                                 prob_threshold = 0.7,
                                 verbose = FALSE)

    pred_results  <- pred_res$jointRes
    cells_assigned <- pred_results$cellTypes

  }else{
    pred_res <- predict_scClassify(exprsMat_test = dgc_mat,
                                 trainRes = trainClass, # or pre-made reference using
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = similarity,
                                 prob_threshold = 0.7,
                                 verbose = FALSE)

    if(length(similarity) > 1){
      pred_results  <- pred_res$ensembleRes
      cells_assigned <- pred_results$cellTypes

    }else{
      column_name <- paste0(metric, "_WKNN_limma")
      cells_assigned <- as.vector(pred_res[[column_name]]$predRes)

    }

  }

  seurat_object@meta.data[['scclassify_labels']] <- cells_assigned
  saveRDS(seurat_object, output_file) # Save seurat object with labels

}


# Ensemble/non-ensembl classify. Train and test model in one
if(method == "train"){

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

  if(length(similarity) > 1){
    results  <- scClassify_res_ensemble$testRes$test$ensembleRes
    cells_assigned <- results$cellTypes

  }else{
    column_name <- paste0(metric, "_WKNN_limma")
    cells_assigned <- as.vector(scClassify_res_ensemble$testRes$test[[column_name]]$predRes)

  }


  seurat_object@meta.data[['scclassify_labels']] <- cells_assigned
  saveRDS(seurat_object, output_file) # Save seurat object with labels

}

DefaultAssay(seurat_object) <- "RNA"

# Labelled UMAP and tSNE dimension plots with scClassify labels
umap_plt <- DimPlot(seurat_object, group.by = "scclassify_labels", reduction = "umap")
umap_labelled <- DimPlot(seurat_object, group.by = "scclassify_labels", label = TRUE, reduction = "umap")

tsne_plt <- DimPlot(seurat_object, group.by = "scclassify_labels", reduction = "tsne")
tsne_labelled <- DimPlot(seurat_object, group.by = "scclassify_labels", label = TRUE, reduction = "tsne")

name<- paste0("Annotation_Figures.dir/scClassify_UMAP_", sample_name, ".eps")
postscript(name)
print(umap_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/scClassify_UMAP_labelled_", sample_name, ".eps")
postscript(name)
print(umap_labelled)
dev.off()

name<- paste0("Annotation_Figures.dir/scClassify_tSNE_", sample_name, ".eps")
postscript(name)
print(tsne_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/scClassify_tSNE_labelled_", sample_name, ".eps")
postscript(name)
print(tsne_labelled)
dev.off()