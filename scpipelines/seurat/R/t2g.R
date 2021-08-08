#!/usr/bin/env Rscript

library("GenomicFeatures")
library("optparse")


option_list = list(
  make_option(c("-i", "--input"), type="character", default="", 
              help="input file GTF file [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
	      ); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


tsg <- makeTxDbFromGFF(opt$input)

k <- keys(tsg, keytype="GENEID")

df <- select(tsg, keys = k, columns="TXNAME", keytype = "GENEID")

write.table(df, file=opt$out, row.names=FALSE, sep="\t")