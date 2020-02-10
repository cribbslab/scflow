#!/usr/bin/env Rscript
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(SingleCellExperiment)
library(optparse)
library(rtracklayer)
library(edgeR)
library(tidyverse)

option_list = list(
  make_option(c("--sample1"), type="character",  
              help="Name of control sample ", metavar="character"),
  make_option(c("--sample2"), type="character", 
              help="Name of sample to compare to sample 1", metavar="character"),
  make_option(c("-g", "--geneset"), type="character", default = "geneset_all.gtf.gz",
              help="Name of geneset file [default= %default]", metavar="character"),
  make_option(c("-s", "--species"), type="character", default = "human",
              help="Species, human or mouse, [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

geneset <- opt$geneset
sample1_name <- opt$sample1
sample2_name <- opt$sample2
species <- opt$species

### Functions ###

fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)

ensembl_to_symbol <- function(result.f, col_name = "gene", clust_results = FALSE){
  result <- result.f
  test <- as.data.frame(result)

  data <- as.vector(test[[col_name]])
  if(species == "human"){
    annots <-  AnnotationDbi::select(org.Hs.eg.db, keys=data,
                                   columns="SYMBOL", keytype = "ENSEMBL")
  } else if (species == "mouse"){
    annots <-  AnnotationDbi::select(org.Mm.eg.db, keys=data,
                                   columns="SYMBOL", keytype = "ENSEMBL")
  } else { stop('Species must be human or mouse. Assign in yaml file under sce subheader.')}
  
    
  result <- merge(test, annots, by.x=col_name, by.y= "ENSEMBL")
  if(clust_results){
    result <- result[rev(order(abs(result$avg_logFC))),]
  }
  else{result <- result[order(result$wilcox),]}
  
  return(result)
}

wilcoxin_test <- function(scepass1, scepass2, GC_length_ensembl_id_output, symbol = FALSE){
  groupA.m <- logcounts(scepass1)
  groupA.m <- na.omit(groupA.m)
  groupB.m <- logcounts(scepass2)
  groupB.m <- na.omit(groupB.m)
  output <- GC_length_ensembl_id_output
  
  
  m.A <- match(rownames(groupA.m), rownames(output))
  gene.lengths <- output$output[m.A]
  
  RPKM.A <- rpkm(groupA.m, gene.length = gene.lengths,normalized.lib.sizes=TRUE)
  
  m.B <- match(rownames(groupB.m), rownames(output))
  gene.lengths <- output$output[m.B]
  
  RPKM.B <- rpkm(groupB.m, gene.length = gene.lengths,normalized.lib.sizes=TRUE)
  
  
  fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
  
  # set up empty data frame
  result.f<-data.frame()
  
  # loop over nrow groupA.m i.e. for each cell
  for(i in 1:nrow(RPKM.A)){
    
    n<-RPKM.A[i,] # take the i'th sample from groupA
    m<-RPKM.B[i,] # take the i'th sample from groupB
    n.f<-n[n>0] # subset all data for genes with count >0 in groupA
    m.f<-m[m>0] # subset all data for genes with count >0 in groupB
    
    x<-length(n) # calulate length of groupA
    y<-length(m) # calculate length of groupB
    
    x.1<-length(n.f) # calculate length of n.f groupA
    y.1<-length(m.f) # calulate length of m.f of groupB
    x.2<-x-x.1 # calulate difference between expression in groupA
    y.2<-y-y.1 # calulate difference between expression in groupB
    
    m.m <- matrix(c(x.1, x.2, y.1, y.2), ncol = 2) # now set up a matrix (contingency table) of sample expression A vs B
    fisher.test.p <-fisher.test(m.m)$p.value # Now use fishers to test if the comparrison between cell frequency is different in group A vs groupB
    
    wilcox.p<-wilcox.test(n,m,paired = FALSE)$p.value # wilcoxon test for difference in expression. ?? However, shouldnt it now be done on genes not rows? Am I missing something? Clarify with supat?
    
    all.p<-c(fisher.test.p,wilcox.p) # combined test
    fisher.p<-fishersMethod (all.p) # Then perform a fishers method function
    
    my.meanA<-mean(RPKM.A[i,]) # calculate the mean of groupA
    my.meanB<-mean(RPKM.B[i,]) # calulcate the mean of groupB
    my.log2.fc<-log2(my.meanB)-log2(my.meanA) # calulate the difference in mean
    my.test.f<-data.frame(gene=rownames(RPKM.A)[i],meanA=my.meanA,meanB=my.meanB,
                          log2fc=my.log2.fc,nCellA=x,nCellB=y,expCellA=x.1,expCellB=y.1,
                          expFractionCellA=x.1/x,expFractionCellB=y.1/y,fisher.test.p=fisher.test.p,wilcox=wilcox.p,fisher=fisher.p) # make a huge dataframe will all of the analysis in
    result.f<-rbind(result.f,my.test.f) # Then bind the result to my tests dataframe.
  }
  
  result.f<-unique(result.f[order(result.f$wilcox),]) # order the result by fishers test
  result.f$p.adjust<-p.adjust(result.f$fisher, method ="BH") # Then perform benjamini-hochberg adjustment
  if(symbol){
    result.f <- ensembl_to_symbol(result.f, col_name = "gene")
  }
  
  return(result.f)
  
}


### Main script ###

## Geneset part ##
GTFfile <- geneset

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="GRCh38", feature.type="exon")
grl <- GenomicRanges::reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), lengths(grl))

elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC
calc_GC_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length)

output <- as.data.frame(output)


## Samples of interest ##
path <- list.dirs(paste0('./SCE.dir/', sample1_name, sep = ""), recursive = FALSE)[1]
pseudo <- unlist(strsplit(path, "/"))[4]
path_sample1 <- paste("SCE.dir", sample1_name, pseudo, "pass.rds", sep = "/")
path_sample2 <- paste("SCE.dir", sample2_name, pseudo, "pass.rds", sep = "/")

sce.pass1 <- readRDS(path_sample1)
sce.pass2 <- readRDS(path_sample2)


results <- wilcoxin_test(scepass1 = sce.pass1, scepass2 = sce.pass2,  GC_length_ensembl_id_output = output, symbol = TRUE)

# Filter mean A and B columns
result_final <- results %>% 
  filter(meanA > 1 & meanB > 1) %>% filter(!is.na(p.adjust) && p.adjust < 0.05) 

file_name = paste0("DE.dir/", sample2_name, "vs", sample1_name, ".txt")
file_name_DE = paste0("DE.dir/", sample2_name, "vs", sample1_name, "_DE.txt")

write.table(results, file=file_name, append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(result_final,file=file_name_DE, append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)


