library(BUSpaRse)
library(DropletUtils)
library(Seurat)
library(ggplot2)

tr2g <- transcript2gene(fasta_file = c("./data/Homo_sapiens.GRCh38.cdna.all.fa.gz", "./data/Mus_musculus.GRCm38.cdna.all.fa.gz"),
                        kallisto_out_path = "./output/out_hgmmDS/")
head(tr2g)

save_tr2g_bustools(tr2g, "./output/tr2g_hgmm.tsv")


res_mat <- read_count_output("./output/out_hgmmDS/genecount",
                             name = "genes", tcc = FALSE)
dim(res_mat)
tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

bc_rank <- barcodeRanks(res_mat)

qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")


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

ggplot(cell_species, aes(n_human_umi, n_mouse_umi, color = species)) +
  geom_point(size = 0.5)
