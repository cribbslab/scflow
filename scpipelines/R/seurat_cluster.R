library(tidyverse)
library(Seurat)
library(optparse)


option_list <- list(
		make_option(c("-i", "--input"), default=NULL,
			help="The input rds file and path"),
		make_option(c("-s", "--sample"), default=NULL,
			help="Sample name"),
		make_option
		
