"""
===================
Pipeline clustering
===================

The pipeline follows pipeline_kb and pipeline_kb-filter-2
and should be ran in the same directory as these pipelines


Authors
=======


Code
====
""" 

# Load modules
from ruffus import *

import sys
import os 
import re

import cgatcore.pipeline as P 
import cgatcore.experiment as E


# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Root of Rmarkdown folder in pipeline folder
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-cluster-3","Rmarkdown")
# R folder in main directory
R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

# filtered files in format:
# RDS_objects.dir/<SAMPLE>_filtered_SeuratObject.rds"
rds_suffixes = "*_filtered_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_SeuratObject.rds"),
	r"")
def cluster(infile, outfile)
	'''
	Rscript task to run seurat clustering
	'''
	
	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_filtered_SeuratObject.rds', file_name).group(1)

	job_memory = "50G"

	statement = Rscript %(R_PATH)s/seurat_cluster.R -i %(infile)s -s %(sample)s

	P.run(statement)

