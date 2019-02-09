#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")
library("Seurat")

option_list = list(
  make_option(c("-w", "--workingdir"), type="character", default=NULL, 
              help="working directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default="sce.rds", 
              help="input file (quants_mat.gz path) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="sce.rds", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--outdir"), type="character", default="NULL", 
              help="output dir name"),
  make_option(c("--mingenes"), type="character", default="NULL", 
              help="Threshold for the minimum number of genes to filter"),
  make_option(c("--maxmitopercent"), type="character", default="NULL", 
              help="Threshold for the maximum percentage of mitochrondrial to filter")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wd = opt$workingdir

so <- readRDS(opt$input)

# see if any meta data is required to be added:
# so <- AddMetaData(so, metadata)

# Identify MT genes, if mouse then ignoring case will
# handle both situations. 
mito.genes <- grep("^MT-", rownames(so@data), value=TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(so@data[mito.genes, ]) / Matrix::colSums(so@raw.data)

# Add the meta data for percent mito to the object
so <- AddMetaData(so, percent.mito, "percent.mito")

svg(paste(opt$outdir,"/","vln.svg", sep=""),width=14,height=7)
VlnPlot(
    so, c("nGene", "nUMI", "percent.mito"), size.title.use=14, size.x.use=12,
    group.by=opt$groupby, x.lab.rot=TRUE, point.size.use=0.1, nCol=3
    )
dev.off()

plot_fn <- function() {
    par(mfrow=c(1, 2))
    GenePlot(s, "nUMI", "percent.mito", cex.use=1)
    GenePlot(s, "nUMI", "nGene", cex.use=1)
    par(mfrow=c(1, 1))


################### Filter cells #########################

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
so <- FilterCells(object = so,
   subset.names = c("nGene", "percent.mito"), 
   low.thresholds = c(opt$mingenes, -Inf),
   high.thresholds = c(Inf, opt$maxmitopercent))


#################### Normalise data ######################

so <- NormalizeData(
    object=so,
    normalization.method="LogNormalize",
    scale.factor=10000
    )

#################### Scale data ######################


so <- ScaleData(object = so,
   vars.to.regress = c("nUMI", "percent.mito"))

# TODO: add function to handle cell cycle gene analysis and regression

#################### Find HVG ######################


so <- FindVariableGenes(object = so,
   mean.function = ExpMean,
   dispersion.function = LogVMR, 
   x.low.cutoff = 0.0125,
   x.high.cutoff = 3,
   y.cutoff = 0.5)

############# Perform Linear dim reduction ################


so <- RunPCA(object = so,
   pc.genes = pbmc@var.genes,
   do.print = FALSE,
   pcs.compute-50)


svg(paste(opt$outdir,"/","pca.svg", sep=""),width=14,height=7)
PCHeatmap(so,
	pc.use=1:12,
	cells.use=min(1000, length(s@cell.names)),
	do.balanced=TRUE,
        label.columns=TRUE,
	cexRow=0.8,
	use.full=FALSE)
dev.off()

svg(paste(opt$outdir,"/","pca_elbow.svg", sep=""),width=14,height=7)
PCElbowPlot(s, num.pc=nPCs)
dev.off()

############# Jack straw ################

PC <- min(dim(so@dr$pca@cell.embeddings)[2],50)

so <- JackStraw(so,
   num.replicate=200,
   num.pc = nPCs,
   do.par=TRUE,
   do.print=FALSE)

so <- JackStrawPlot(s0, PCs=1:PC)

svg(paste(opt$outdir,"/","jackstraw.svg", sep=""),width=14,height=7)
so@dr$pca@misc$jackstraw.plot
dev.off()


saveRDS(object = so, file = opt$out)