"""
===================
Pipeline Integration
===================
The pipeline follows pipeline_kb, pipeline_kb-sampleqc-1 and pipeline_kb-filter-2.

It takes filtered Seurat objects and performs integrations by two methods, Seurat and Harmony.
Each dataset is first processed individually by SCTransform normalization and then perfrom integration and clustering by Seurat.
The integrated data could be furthered accepted by Harmony with integration/clustering.  
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
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-integration-4","Rmarkdown")
# R folder in main directory

R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

@follows(mkdir("Integration_Figures.dir"))
@originate("RDS_objects.dir/SCT_integrated_SeuratObject.rds")
def integrate(outfile):
	'''
	R script task to run seurat integration
	'''
	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

	job_memory = "70G"

	statement = '''
	Rscript %(R_PATH)s/seurat_integrate.R '''

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

	RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-integration-4","Rmarkdown")
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
