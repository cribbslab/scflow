"""
===================
Pipeline annotation
===================

The pipeline follows:
pipeline_kb, pipeline_kb-sampleqc-1 and pipeline_kb-filter-2, pipeline_kb-cluster-3, pipeline_kb-integration-4
Annotation of the clusters are performed using either SingleR, clustifyR or scClassify.
Output figures are visualised on a Rmarkdown notebook / html file.

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
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_annotation-6","Rmarkdown")
# R folder in main directory
R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

# Integrated output:
# Harmony_integrated_SeuratObject.rds
# SCT_integrated_SeuratObject.rds
filtered_suffixes = "*_integrated_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

######################
# Seurat markers
######################

@follows(mkdir("Annotation_Figures.dir"))
@follows(mkdir("Annotation_stats.dir"))
@active_if(PARAMS['markerdiff'])
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
	r"Annotation_stats.dir/ConservedMarkers_\1.csv")
def integrated_markers(infile, outfile):
	'''
	Find conserved and differentially expressed markers across conditions
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))
	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_integrated_SeuratObject.rds', file_name).group(1)

	meta = PARAMS['markers_meta_data']
	group = PARAMS['markers_group']
	DE = PARAMS['markers_DE_versus'].replace(" ", "---")

	predefined_list = PARAMS['predefined_list']
	if predefined_list:
		predef_options = "--predefined " + predefined_list
	else:
		predef_options = ""

	statement = '''
	Rscript %(R_PATH)s/conserved_differential_markers.R -i %(infile)s -m %(meta)s -g %(group)s 
	--de %(DE)s %(predef_options)s -s %(sample)s'''

	P.run(statement)

######################
# Generate reference
######################

@active_if(PARAMS['reference_generate'])
@originate("reference_sce.rds")
def reference_generate(outfile):
	'''
	Generate a SCE rds file from scratch using either celldex or scRNAseq
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

	ref = PARAMS['reference_package']
	celldex_ref = PARAMS['reference_celldex_reference_name']
	scRNAseq_ref = PARAMS['reference_scRNAseq_reference_name']

	scRNAseq_option = PARAMS['reference_scRNAseq_option']
	if scRNAseq_option:
		scRNAseq_options = "-o " + scRNAseq_option
	else:
		scRNAseq_options = ""

	statement = '''
	Rscript %(R_PATH)s/reference_sce.R -r %(ref)s  -c %(celldex_ref)s 
	-s %(scRNAseq_ref)s %(scRNAseq_options)s --outfile %(outfile)s'''

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

@active_if(PARAMS['singler_run'])
@follows(reference_generate, reference_copy)
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
	r"RDS_objects.dir/\1_singleR_annotated_SeuratObject.rds")
def singleR(infile, outfile):
	'''
    R script task to run SingleR package for annotation
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_integrated_SeuratObject.rds', file_name).group(1)

	ref = "reference_sce.rds"
	DE = PARAMS['singler_DEmethod']
	method = PARAMS['singler_method']

	predefined_list = PARAMS['predefined_list']
	if predefined_list:
		predef_options = "--predefined " + predefined_list
	else:
		predef_options = ""

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/singleR.R -i %(infile)s -s %(sample)s -r %(ref)s  -d %(DE)s 
	-m %(method)s -o %(outfile)s %(predef_options)s '''

	P.run(statement)

############
# clustifyr
############

@active_if(PARAMS['clustifyr_run'])
@follows(reference_generate, reference_copy)
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
	r"RDS_objects.dir/\1_clustifyr_annotated_SeuratObject.rds")
def clustifyr(infile, outfile):
	'''
    R script task to run clustifyr package for annotation
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_integrated_SeuratObject.rds', file_name).group(1)

	ref = "reference_sce.rds"
	dim_red = PARAMS['clustifyr_dimRed']
	var_features = PARAMS['clustifyr_var_features']

	predefined_list = PARAMS['predefined_list']
	if predefined_list:
		predef_options = "--predefined " + predefined_list
	else:
		predef_options = ""

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/clustifyr.R -i %(infile)s -s %(sample)s -r %(ref)s  -d %(dim_red)s 
	-v %(var_features)s -o %(outfile)s %(predef_options)s'''

	P.run(statement)

#############
# scClassify
#############

@active_if(PARAMS['scclassify_run'])
@follows(reference_generate, reference_copy)
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
	r"RDS_objects.dir/\1_scclassify_annotated_SeuratObject.rds")
def scclassify(infile, outfile):
	'''
    R script task to run scClassify package for annotation
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_integrated_SeuratObject.rds', file_name).group(1)

	ref = "reference_sce.rds"
	pretrained = PARAMS['scclassify_pretrained']
	method = PARAMS['scclassify_method']
	sim = PARAMS['scclassify_similarity']
	sim = sim.replace(" ", "_")

	predefined_list = PARAMS['predefined_list']
	if predefined_list:
		predef_options = "--predefined " + predefined_list
	else:
		predef_options = ""

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/scclassify.R -i %(infile)s -s %(sample)s -r %(ref)s  -m %(method)s
	--pretrained %(pretrained)s --similarity %(sim)s -o %(outfile)s %(predef_options)s'''

	P.run(statement)

@follows(integrated_markers,singleR, clustifyr, scclassify)
@originate("Annotation.html")
def annotation_rmarkdown(outfile):
	'''
	Rmarkdown and html to visualise annotation
	'''

	RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_annotation-6","Rmarkdown")
	job_memory = "50G"

	statement = ''' 
	cp %(RMD_ROOT)s/Annotation.Rmd . &&
	R -e "rmarkdown::render('Annotation.Rmd', output_file='Annotation.html')" '''

	P.run(statement)

@follows(integrated_markers, reference_generate, reference_copy,
		 singleR, clustifyr, scclassify, annotation_rmarkdown)
def full():
	pass

def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))
