library(BUSpaRse)
library(DropletUtils)
library(Seurat)
library(ggplot2)


library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default="", 
              help="input file GTF file [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
	      ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



res_mat <- read_count_output(opt$input,
                             name = "genes", tcc = FALSE)

tot_counts <- Matrix::colSums(res_mat)

bc_rank <- barcodeRanks(res_mat)

out_pdf = paste0(opt$out, "/kneeplot.pdf")
pdf(out_pdf)
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
dev.off()


res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
dim(res_mat)

library(tidyverse)
gene_species <- ifelse(str_detect(rownames(res_mat), "^ENSMUSG"), "mouse", "human")
mouse_inds <- gene_species == "mouse"
human_inds <- gene_species == "human"
# mark cells as mouse or human
cell_species <- tibble(n_mouse_umi = Matrix::colSums(res_mat[mouse_inds,]),
                       n_human_umi = Matrix::colSums(res_mat[human_inds,]),
                       tot_umi = Matrix::colSums(res_mat),
                       prop_mouse = n_mouse_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)

cell_species <- cell_species %>% 
  mutate(species = case_when(
    prop_mouse > 0.9 ~ "mouse",
    prop_human > 0.9 ~ "human",
    TRUE ~ "mixed"
  ))

out_pdf = paste0(opt$out, "/barnyard.pdf")
pdf(out_pdf)
ggplot(cell_species, aes(n_human_umi, n_mouse_umi, color = species)) +
  geom_point(size = 0.5)
dev.off()
