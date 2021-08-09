"""
===================================
Pipeline Integration Reciprocal PCA
===================================

This pipeline is for large datasets that cannot be integrated by Seurat CCA.

The pipeline follows pipeline_kb, pipeline_kb-sampleqc-1 and pipeline_kb-filter-2.

It takes filtered Seurat objects and performs integrations based on reciprocal PCA.

the workflow consists of the following steps:

    Create a list of Seurat objects to integrate
    Perform normalization, feature selection, and scaling separately for each dataset
    Run PCA on each object in the list
    Integrate datasets, and proceed with joint analysis

The integrations performed by both the methods could be visualized as tSNE and UMAP dimensional reduction.

=======
Code
=======
""" 

# Load modules
from ruffus import *

import sys
import os 
import re
import glob
import sqlite3

import cgatcore.pipeline as P 
import cgatcore.experiment as E


# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Root of Rmarkdown folder in pipeline folder
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_integration-rcpa","Rmarkdown")
# R folder in main directory

R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

@follows(mkdir("Integration_Figures-rcpa.dir"))
@originate("RDS_objects.dir/rcpa_integrated_SeuratObject.rds")
def integrate(outfile):
	'''
	R script task to run seurat integration
	'''
	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

	job_memory = "70G"

	statement = '''
	Rscript %(R_PATH)s/rcpa_integrate.R '''

	P.run(statement)


@follows(integrate)
@originate("RDS_objects.dir/Harmony_integrated_SeuratObject.rds")
def harmony(outfile):
	'''
	R script task to run harmony integration
	'''
	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

	job_memory = "70G"

	statement = '''
	Rscript %(R_PATH)s/harmony_integrate.R'''

	P.run(statement)

@follows(harmony)
@merge([integrate, harmony],
	"Integration.html")
def integrate_rmarkdown(infile, outfile):
	'''
	R markdown to visualise clustering and dimensional reduction after integration
	'''

	RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_integration-5","Rmarkdown")
	job_memory = "50G"

	statement = ''' 
	cp %(RMD_ROOT)s/Integration.Rmd . &&
	R -e "rmarkdown::render('Integration.Rmd', output_file='Integration.html')" '''

	P.run(statement)

@follows(integrate, harmony, integrate_rmarkdown)
def full():
	pass

def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))
