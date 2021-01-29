"""
===================
Pipeline annotation
===================

The pipeline follows:
pipeline_kb, pipeline_kb-sampleqc-1 and pipeline_kb-filter-2, pipeline_kb-cluster-3, pipeline_kb-integration-4
Annotation of the clusters are performed using either SingleR, clustifyR or scClassify.
Output figures are visualised on a Rmarkdown notebook and html file.

Authors
=======
Anna James-Bott- anna.james-bott@st-hildas.ox.ac.uk
Dr Adam P Cribbs- adam.cribbs@imm.ox.ac.uk

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
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-annotation-5","Rmarkdown")
# R folder in main directory
R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))


##### Guessing what Chen-yi's output will be #####
filtered_suffixes = "*_filtered_clustered_integrated_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

######################
# Generate reference
######################

@active_if(PARAMS['reference_generate'])
@originate("reference_sce.rds")
def reference_generate(outfile):
	'''
	Generate a SCE rds file from scratch using either celldex or scRNAseq
	'''

	ref = PARAMS['reference_package']
	celldex_ref = PARAMS['reference_celldex_reference_name']
	scRNAseq_ref = PARAMS['reference_scRNAseq_reference_name']
	scRNAseq_option = PARAMS['reference_scRNAseq_option']

	statement = '''
	Rscript %(R_PATH)s/reference_sce.R -r %(ref)s  -c %(celldex_reference_name)s 
	-s %(scRNAseq_ref)s -o %(scRNAseq_option)s --outfile %(outfile)s'''

	P.run(statement)

@active_if(not PARAMS['reference_generate'])
@originate("reference_sce.rds")
def reference_copy(outfile):
	'''
	Copy existing SCE rds file to reference_sce.rds location
	'''

	orig_loc = PARAMS['reference_path']

	statement = '''cp %(orig_loc)s %(outfile)s'''

	P.run(statement)


##########
# singleR
##########

# Unsure about output files... or regex
@follows(mkdir("Annotation_Figures.dir"))
@follows(reference_generate, reference_copy)
@active_if(PARAMS['singler_run'])
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_clustered_integrated_SeuratObject.rds"),
	r"RDS_objects.dir/\1_filtered_clustered_integrated_annotated_SeuratObject.rds")
def singleR(infile, outfile):
	'''
    R script task to run SingleR package for annotation
	'''


	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_filtered_clustered_integrated_SeuratObject.rds', file_name).group(1)

	ref = "reference_sce.rds"
	DE = PARAMS['singler_DEmethod']
	method = PARAMS['singler_method']

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/singleR.R -i %(infile)s -s %(sample)s -r %(ref)s  -d %(DEmethod)s -m %(method)s '''

	P.run(statement)

############
# clustifyr
############

@active_if(PARAMS['clustifyr_run'])
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_clustered_integrated_SeuratObject.rds"),
	r"RDS_objects.dir/\1_filtered_clustered_integrated_annotated_SeuratObject.rds")
def clustifyr(infile, outfile):
	'''
    R script task to run clustifyr package for annotation
	'''


	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_filtered_clustered_integrated_SeuratObject.rds', file_name).group(1)

	ref = "reference_sce.rds"

	DE = PARAMS['singler_DEmethod']
	method = PARAMS['singler_method']

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/clustifyr.R -i %(infile)s -s %(sample)s -r %(ref)s   '''

	P.run(statement)

@follows(singleR, clustifyr, scclassify, )
def full():
	pass

def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))