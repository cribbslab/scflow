library("BUSpaRse")
library("DropletUtils")
library("ggplot2")
library("tidyverse")
library("optparse")


option_list = list(
  make_option(c("-i", "--infile"), type="character", 
              help="Sorted bus text file", metavar="character"),
  make_option(c("--estcells"), type="integer", default = 1000,
              help="Estimated number of cells, [default= %default]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", 
              help="Output file. SCE object", metavar="character"),
  make_option(c("-t", "--t2g"), type="character", default = "transcript2geneMap.tsv",
              help="Transcript to gene file.", metavar="character"),
  make_option(c("-g", "--geneset"), type="character", default = "geneset_all.gtf.gz",
              help="Name of geneset file, [default= %default]", metavar="character")
); 

# make_option(c("-g", "--geneset"), type="character", default = "geneset_all.gtf.gz",
#              help="Name of geneset file, [default= %default]", metavar="character"),

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

geneset <- opt$geneset
bus_text <- opt$infile
t2gmap <- opt$t2g 
outfile <- opt$outfile
estcells <- opt$estcells
bus_folder <- dirname(bus_text)

if (!file.exists("t2g_rename.tsv")){
tr2g_tsv <- read_tsv(t2gmap, col_names = TRUE) %>% as.data.frame()
colnames(tr2g_tsv) <- c("transcript", "gene")
write_tsv(tr2g_tsv, "t2g_rename.tsv")
}


#tr2g <- tr2g_gtf(file= geneset, transcript_version = NULL, gene_version = NULL, gene_name = NULL)
tr2g <- sort_tr2g(file = "t2g_rename.tsv", kallisto_out_path = bus_folder, verbose = FALSE)


TCC <- make_sparse_matrix(bus_text, tr2g = tr2g, est_ncells = estcells, est_ngenes = nrow(tr2g), whitelist = NULL, TCC = TRUE, gene_count = FALSE)
TCC_file <- paste(c(bus_folder, "TCC.mat"), collapse = "/")
saveRDS(TCC, file=TCC_file)

GC_mat <- make_sparse_matrix(bus_text, tr2g = tr2g, est_ncells = estcells, est_ngenes = nrow(tr2g), whitelist = NULL, TCC = FALSE, gene_count = TRUE)
GC_mat_file <- paste(c(bus_folder, "GC.mat"), collapse = "/")
saveRDS(GC_mat, file=GC_mat_file)

tot_counts2 <- Matrix::colSums(GC_mat)
summary(tot_counts2)

bc_rank <- barcodeRanks(GC_mat)

# Rotated knee plot
knee_plot_file <- paste(c(bus_folder, "knee_plot.png"), collapse = "/")
png(knee_plot_file)
qplot(bc_rank$total, bc_rank$rank, geom = "line") + 
  geom_vline(xintercept = bc_rank$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = bc_rank$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(bc_rank$knee, bc_rank$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  annotate("text", y = 500, x = 1.5 * c(bc_rank$knee, bc_rank$inflection),
           label = c(bc_rank$knee, bc_rank$inflection), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
dev.off()
# This is assumed to represent the difference between empty droplets with little RNA and cell-containing dropleta with more RNA. More sophisticated method below:

# Threshold using point of inflection
GC_mat_filt <- GC_mat[, tot_counts2 > bc_rank$inflection]
#dim(GC_mat_filt)
v_1 <- log2(GC_mat_filt + 1)
sce_inflection <- SingleCellExperiment(assays = list(counts = GC_mat_filt, logcounts = v_1))
out_inflection <- paste0(dirname(outfile), "/inflection_sce.rds")
saveRDS(sce_inflection, file= out_inflection)

set.seed(100)
e.out <- emptyDrops(GC_mat)

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm = TRUE)
is.cell_logical <- is.cell %>% replace_na(FALSE)
GC_mat_filt2 <- GC_mat[,is.cell_logical]

logprob_plot_file <- paste(c(bus_folder, "log_probability_plot.png"), collapse = "/")
png(logprob_plot_file)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"), xlab = "Total UMI Count", ylab = "-Log Probability")
dev.off()
#table(Limited = e.out$Limited, Significant = is.cell)

mat <- GC_mat_filt2
v <- log2(mat + 1)
sce <- SingleCellExperiment(assays = list(counts = mat, logcounts = v))
saveRDS(sce, file= outfile)




