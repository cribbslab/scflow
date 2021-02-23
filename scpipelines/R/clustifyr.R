library(Seurat)
library(tidyverse)
library(optparse)
library(scRNAseq)
library(celldex)
library(clustifyr)


option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, filtered clusterd integrated Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option(c("-o", "--output"), default=NULL,
			help="The output rds path."),
        make_option(c("-r", "--reference"), default="reference_sce.rds", type = "character",
			help="Location of reference sce rds file"),
		make_option(c("-d", "--dimReduction"), default="umap", type = "character",
			help="Dimension reduction technique, e.g. umap, pca, tsne."),
		make_option(c("-v", "--varFeatures"), default=1, type = "integer",
			help="Whether to use previously generated list of variable features for clustify function."),
        make_option(c("--predefined"), default=NULL,
			help="User pre-defined list of marker genes")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
output_file <- opt$output
sample_name <- opt$sample
reference_sce_loc <- opt$reference
dim_red <- opt$dimReduction
var_features <- opt$varFeatures

# Load seurat object
seurat_object <- readRDS(input_file)

#seurat_object_mat <- as.matrix(object_data(seurat_object, "data"))
#metadata <- seurat_meta(seurat_object, dr = dr)

# Load in reference sce object
ref_sce <- readRDS(reference_sce_loc)
ref_mat <- as.matrix(assay(ref_sce))
colnames(ref_mat) <- ref_sce$label

# Variable features
if(var_features){
	vargenes <- VariableFeatures(seurat_object)
}else{
	vargenes <- NULL
}

res <- clustify(input=seurat_object, ref_mat= ref_mat, cluster_col="seurat_clusters", dr = dim_red, query_genes = vargenes, obj_out = FALSE)
res_so <- clustify(input=seurat_object, ref_mat= ref_mat, cluster_col="seurat_clusters", dr = dim_red, query_genes = vargenes) # Outputs seurat object

name_file <- paste0(c("Annotation_stats.dir/clustifyr_correlation_matrix_", sample_name,".csv"), collapse="")
write.csv(x = res, file = name_file, row.names = TRUE)

name_file2 <- paste0(c("Annotation_stats.dir/clustifyr_cluster_annotations_", sample_name,".csv"), collapse="")
res2 <- cor_to_call(cor_mat = res, cluster_col = "seurat_clusters")
res2$seurat_clusters <- as.integer(res2$seurat_clusters)
res2 <- res2 %>% dplyr::arrange(seurat_clusters)
write_csv(res2, name_file2)


cell_types <- res_so@meta.data$type
r <- res_so@meta.data$r

seurat_object@meta.data[['clustifyr_labels']] <- cell_types
seurat_object@meta.data[['clustifyr_rvalues']] <- r

saveRDS(seurat_object, output_file) # Save seurat object with labels

DefaultAssay(seurat_object) <- "RNA"

# Labelled UMAP and tSNE dimension plots with clustifyr labels
umap_plt <- DimPlot(seurat_object, group.by = "clustifyr_labels", reduction = "umap")
umap_labelled <- DimPlot(seurat_object, group.by = "clustifyr_labels", label = TRUE, reduction = "umap")

tsne_plt <- DimPlot(seurat_object, group.by = "clustifyr_labels", reduction = "tsne")
tsne_labelled <- DimPlot(seurat_object, group.by = "clustifyr_labels", label = TRUE, reduction = "tsne")

name<- paste0("Annotation_Figures.dir/clustifyr_UMAP_", sample_name, ".eps")
postscript(name)
print(umap_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/clustifyr_UMAP_labelled_", sample_name, ".eps")
postscript(name)
print(umap_labelled)
dev.off()

name<- paste0("Annotation_Figures.dir/clustifyr_tSNE_", sample_name, ".eps")
postscript(name)
print(tsne_plt)
dev.off()

name<- paste0("Annotation_Figures.dir/clustifyr_tSNE_labelled_", sample_name, ".eps")
postscript(name)
print(tsne_labelled)
dev.off()

# Pre-defined list of markers
predefined <- opt$predefined
if(!is.null(predefined)){
  predefined_list <- str_split(predefined, "_")[[1]]
  number_markers <- length(predefined_list)

  for(marker in predefined_list){

    vln_plt <- VlnPlot(seurat_object, features = marker,  group.by = "clustifyr_labels",
    pt.size = 0)

    name<- paste0("Annotation_Figures.dir/clustifyr_ViolinPlot_", sample_name, "_", marker ,".eps")
    postscript(name)
    print(vln_plt)
    dev.off()
  }
}