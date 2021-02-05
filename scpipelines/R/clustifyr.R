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
			help="Dimension reduction technique, e.g. umap, pca, tsne.")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
output_file <- opt$output
sample_name <- opt$sample
reference_sce_loc <- opt$reference
dim_red <- opt$dimReduction

# Load seurat object
seurat_object <- readRDS(input_file)

#seurat_object_mat <- as.matrix(object_data(seurat_object, "data"))
#metadata <- seurat_meta(seurat_object, dr = dr)

# Load in reference sce object
ref_sce <- readRDS(reference_sce_loc)
ref_mat <- as.matrix(assay(ref_sce))
colnames(ref_mat) <- ref_sce$label

# With their function
new_ref_matrix_sce <- object_ref(
  input = ref_sce,               # SCE object
  cluster_col = "label"       # name of column in colData containing cell identities
)


res <- clustify(input=seurat_object, ref_mat= new_ref_matrix_sce, cluster_col="seurat_clusters", dr = dim_red )
res1 <- res
res1$cluster <- rownames(res)

name_file <- paste0(c("clustifyr_correlation_matrix_", sample_name,".csv"), collapse="")
write_csv(res1, name_file)

name_file2 <- paste0(c("clustifyr_cluster_annotations_", sample_name,".csv"), collapse="")
res2 <- cor_to_call(cor_mat = res, cluster_col = seurat_clusters)
write_csv(res2, name_file2)

res3 <- dplyr::select(res2, seurat_clusters, clustifyr_labels = type)
tib <- tibble(seurat_clusters = seurat_object@meta.data$seurat_clusters )

joined <- plyr::join(tib,res3, by="seurat_clusters", type ="left")
seurat_object@meta.data$clustifyr_labels <- joined$clustifyr_labels

saveRDS(seurat_object, output_file) # Save seurat object with labels
