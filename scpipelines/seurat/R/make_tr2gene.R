#!/usr/bin/env Rscript

library(BUSpaRse)
library(optparse)

option_list = list(
  make_option(c("-i", "--input1"), type="character", default="", 
              help="input fastafile if there are mixed species then specify -j too", metavar="character"),
  make_option(c("-j", "--input2"), type="character", default=FALSE, 
              help="Second input fastafile", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="", 
              help="output dir name [default= %default]", metavar="character"),
  make_option(c("-f", "--outfile"), type="character", default="", 
              help="output filename for tsv", metavar="character")
	      );
 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


tr2g <- BUSpaRse::tr2g_fasta(file=opt$input1)

write.table(tr2g, file=opt$outfile, quote=FALSE, sep="\t", col.names=F, row.names=F)