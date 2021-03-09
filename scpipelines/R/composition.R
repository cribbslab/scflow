# Script to do composition analysis

library(tidyverse)
library(Seurat)
library(optparse)
library(RColorBrewer)

option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds path, annotated Seurat object"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
        make_option(c("--seed"), default=NULL,
			help="Set seed for sampling. Default = null, assign integer to set seed.")
)

# Read in options
opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
sample_name <- opt$sample
seed_no <- opt$seed

if(!is.null(seed_no)){
  set.seed(seed_no)
}

seurat_object <- readRDS(input_file)

so <- seurat_object@meta.data
clusters <- unique(so$seurat_cluster)%>% as.character() %>% as.numeric() %>% sort()
samples <- unique(so$sample_name)

composite <- c()

for(sample_test in samples){
    so_col <- dplyr::filter(so, sample_name == sample_test)
    counts <- table(so_col$seurat_cluster) %>% as.vector()
    composite <- cbind(composite,counts)



}
composite_tib <- as_tibble(composite)
colnames(composite_tib) <- samples
composite_tib$cluster <- clusters
composite_tib <- dplyr::select(composite_tib, cluster, everything())

proportion <- prop.table(composite, margin=2) %>% as_tibble()
colnames(proportion) <- samples
proportion$cluster <- clusters
proportion <- dplyr::select(proportion, cluster, everything())

gather_composition <- gather(composite_tib, "Sample", "Number_cells", -cluster)
gather_composition_prop <- gather(proportion, "Sample", "Proportion_cells", -cluster)

gathered <- gather_composition
gathered$Proportion_cells <- gather_composition_prop$Proportion_cells
gathered$cluster <- as.factor(gathered$cluster)

# Barcharts

n <- len(samples)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols=sample(col_vector, n)

bar_prop <- ggplot(data=gathered, aes(x=cluster, y=Proportion_cells, fill=Sample)) +
  geom_bar(stat="identity") + theme_minimal() + scale_fill_manual(values=cols) +
  labs(x= "Cluster", y = "Proportion of cells (for each sample")

bar_number <- ggplot(data=gathered, aes(x=cluster, y=Number_cells, fill=Sample)) +
  geom_bar(stat="identity") + theme_minimal() + scale_fill_manual(values=cols) +
  labs(x= "Cluster", y = "Number of cells")

name_bar_prop <- paste0("Annotation_Figures.dir/Composition_Clusters_Barchart_Proportion_", sample_name ,".eps")
name_bar_numb <- paste0("Annotation_Figures.dir/Composition_Clusters_Barchart_NumberCells_", sample_name ,".eps")

ggsave(filename = name_bar_prop, plot = bar_prop, device = "eps")
ggsave(filename = name_bar_numb, plot = bar_number, device = "eps")

# Dot plots

n <- len(clusters)
set.seed(seed_no + 10) # Doesn't share colours with first pallete
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)

dot_prop <- ggplot(data = gathered, aes(x= Sample, y = cluster, color=cluster)) +
  geom_point(aes(size = Proportion_cells)) + scale_color_manual(values = col) +
  scale_size_continuous(range = c(0,8)) +theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    labs(x= "Sample", y = "Cluster")

dot_numb <- ggplot(data = gathered, aes(x= Sample, y = cluster, color =cluster)) +
  geom_point(aes(size = Number_cells)) + scale_color_manual(values = col) +
  scale_size_continuous(range = c(0,8)) +theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    labs(x= "Sample", y = "Cluster")

name_dot_prop <- paste0("Annotation_Figures.dir/Composition_Clusters_DotPlot_Proportion_", sample_name ,".eps")
name_dot_numb <- paste0("Annotation_Figures.dir/Composition_Clusters_Dotplot_NumberCells_", sample_name ,".eps")

ggsave(filename = name_dot_prop, plot = dot_prop, device = "eps")
ggsave(filename = name_dot_numb, plot = dot_numb, device = "eps")

# To do with annotations