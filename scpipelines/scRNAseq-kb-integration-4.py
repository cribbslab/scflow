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
====
""" 

# Load modules
from ruffus import *

import sys
import os 
import re
import glob

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

filtered_suffixes = "*_filtered_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

@follows(mkdir("Clustering_Figures.dir"))
@merge(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_SeuratObject.rds"),
	r"RDS_objects.dir/\1_SCT_integrated_SeuratObject.rds")

def integrate(infiles, outfile):
	'''
	R script task to run seurat integration
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))
	infiles = glob("filtered_SeuratObject.rds")
        
        #file_name = os.path.basename(infile)
        #sample = re.match(r'(\S+)_filtered_SeuratObject.rds', file_name).group(1)

	res = PARAMS['resolution']
	nvf = PARAMS['num_variable_features']
	
	job_memory = "70G"

	statement = '''
	Rscript %(R_PATH)s/seurat_integrate.R -i %(infiles)s -v %(nvf)s --resolution %(res)s'''

	P.run(statement)

@follows(integrate)
@transform(integrate, suffix("RDS_objects.dir/SCT_integrated_SeuratObject.rds"),
		   r"RDS_objects.dir/Harmony_integrated_SeuratObject.rds")
def harmony(infile, outfile):
	'''
	R script task to run harmony integration
	'''
	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))
	#file_name = os.path.basename(infile)


	job_memory = "70G"

	statement = '''
	Rscript %(R_PATH)s/harmony_integrate.R -i %(infile)s'''

	P.run(statement)

#@follows(find_markers)
@follows(mkdir("Clustering_Figures.dir"))
@merge(integrate,
	"Integrate.html")
def integrate_rmarkdown(infile, outfile):
	'''
	R markdown to visualise clustering and dimensional reduction
	'''

	RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-integration-4","Rmarkdown")
	job_memory = "50G"

	statement = ''' 
	cp %(RMD_ROOT)s/Integrate.Rmd . &&
	R -e "rmarkdown::render('Integrate.Rmd', output_file='Integrate.html')" '''

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
