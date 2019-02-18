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
              help="Threshold for the maximum percentage of mitochrondrial to filter"),
  make_option(c("--reductionmethod"), default="pca",
                help="Type of dimensionality reduction to perform. Options include 'pca', 'ica'"),
  make_option(c("--npcs"), type="character", default="NULL", 
              help="Number of PC components for elbow plot"),
  make_option(c("--pccomponents"), type="integer", default=10,
                help="The number of principle components to use"),
  make_option(c("--resolution"), type="double", default=1,
                help="cluster resolution"),
  make_option(c("--algorithm"), type="integer", default=3,
                help="1=original Louvain, 2=Louvain multilevel, 3=SLM")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wd = opt$workingdir

so <- readRDS(opt$input)

mito.genes <- grep("^MT-", rownames(so@data), value=TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(so@data[mito.genes, ]) / Matrix::colSums(so@raw.data)

# Add the meta data for percent mito to the object
so <- AddMetaData(so, percent.mito, "percent.mito")


svg(paste(wd,"vln.svg", sep=""),width=14,height=7)
VlnPlot(
    so, c("nGene", "nUMI", "percent.mito"), size.title.use=14, size.x.use=12, x.lab.rot=TRUE, point.size.use=0.1, nCol=3
    )
dev.off()



################### Filter cells #########################

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
so <- FilterCells(object = so,
   subset.names = c("nGene", "percent.mito"), 
   low.thresholds = c(as.numeric(opt$mingenes), -Inf),
   high.thresholds = c(Inf, as.numeric(opt$maxmitopercent)))


#################### Normalise data ######################

so <-  NormalizeData(object = so, normalization.method = "LogNormalize", 
                     scale.factor = 10000)


#################### Find HVG ######################

svg(paste(wd,"variable_genes.svg", sep=""),width=14,height=7)
so <- FindVariableGenes(object = so,
   mean.function = ExpMean,
   dispersion.function = LogVMR, 
   x.low.cutoff = 0.0125,
   x.high.cutoff = 3,
   y.cutoff = 0.5)
dev.off()

#################### Scale data ######################


so <- ScaleData(so, genes.use = so@var.genes)

# TODO: add function to handle cell cycle gene analysis and regression



############# Perform Linear dim reduction ################


so <- RunPCA(object = so,
   pc.genes = so@var.genes,
   do.print = FALSE,
   pcs.compute=as.numeric(opt$npcs))


svg(paste(wd,"pca_results_viz.svg", sep=""),width=14,height=7)
VizPCA(so, pcs.use = 1:2)
dev.off()


svg(paste(wd,"pca_results.svg", sep=""),width=14,height=7)
PCAPlot(so, dim.1 = 1, dim.2 = 2)
dev.off()

svg(paste(wd,"pca.svg", sep=""),width=14,height=7)
PCHeatmap(so,
	pc.use=1:3,
	cells.use=min(1000, length(so@cell.names)),
	do.balanced=TRUE,
        label.columns=TRUE,
	cexRow=0.8,
	use.full=FALSE)
dev.off()

svg(paste(wd,"pca_elbow.svg", sep=""),width=14,height=7)
PCElbowPlot(so, num.pc=opt$npcs)
dev.off()


############# Find Clusters #############

so <- ProjectPCA(so, do.print = FALSE)

so <-  FindClusters(so,
                  reduction.type=opt$pccomponents,
                  dims.use = 1:as.numeric(opt$npcs),
                  resolution = opt$resolution,
                  algorithm = opt$algorithm,
                  print.output = FALSE,
                  save.SNN = F)


nclusters <- length(unique(so@ident))


so <- buildClusterTree(so,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:11)


############# Run tsne ##################

so <- RunTSNE(so,
	     dims.use=comps,
	     do.fast=T)


svg(paste(wd,"tsne.svg", sep=""),width=14,height=7)
TSNEPlot(so)
dev.off()


############# Jack straw ################

#PC <- min(dim(so@dr$pca@cell.embeddings)[2],50)

#so <- JackStraw(so)

#so <- JackStrawPlot(s0)

#svg(paste(wd,"jackstraw.svg", sep=""), width=14, height=7)
#so@dr$pca@misc$jackstraw.plot
#dev.off()


saveRDS(object = so, file = opt$out)