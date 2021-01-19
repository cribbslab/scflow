"""
===================
Pipeline clustering
===================

The pipeline follows pipeline_kb, pipeline_kb-sampleqc-1 and pipeline_kb-filter-2
It performs PCA, tSNE and UMAP dimensional reduction. Markers for each cluster are found
Clustering is visualised in an Rmarkdown notebook.

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

filtered_suffixes = "*_filtered_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

@follows(mkdir("Clustering_Figures.dir"))
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_SeuratObject.rds"),
	r"RDS_objects.dir/\1_filtered_clustered_SeuratObject.rds")
def cluster(infile, outfile):
	'''
	R script task to run seurat clustering
	'''

	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))
	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_filtered_SeuratObject.rds', file_name).group(1)

	res = PARAMS['resolution']
	nvf = PARAMS['num_variable_features']
	red = PARAMS['reduction_technique']

	num_dimensions = PARAMS["num_dimensions"]
	if num_dimensions:
		# User defined number of dimensions to reduce to (look at elbow and jackstraw plot post hoc)
		num_dimensions_option = "-d " + str(num_dimensions)
	else:
		# Calculated with embeddings
		num_dimensions_option = " "

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/seurat_cluster.R -i %(infile)s -s %(sample)s -v %(nvf)s --reddim %(red)s
	--resolution %(res)s %(num_dimensions_option)s'''

	P.run(statement)

@follows(mkdir("clustering_markers.dir"))
@transform(cluster,
		   regex("RDS_objects.dir/(\S+)_filtered_clustered_SeuratObject.rds"),
		   r"clustering_markers.dir/\1_markers.csv")
def find_markers(infile, outfile):
	'''
	R script to find markers for each cluster
	'''
	R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))
	file_name = os.path.basename(infile)
	sample = re.match(r'(\S+)_clustered_filtered_SeuratObject.rds', file_name).group(1)

	min_percent = PARAMS['min_percent']
	logfc_thresh = PARAMS['logfc_thresh']
	test_use = PARAMS['test_use']

	job_memory = "50G"

	statement = '''
	Rscript %(R_PATH)s/find_markers.R -i %(infile)s -s %(sample)s  -l %(logfc_thresh)s -t %(test_use)s 
	--minPct %(min_percent)s '''

	P.run(statement)

#@follows(find_markers)
@follows(mkdir("Clustering_Figures.dir"))
@merge(cluster,
	"Cluster.html")
def cluster_rmarkdown(infile, outfile):
	'''
	R markdown to visualise clustering and dimensional reduction
	'''

	RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-cluster-3","Rmarkdown")
	job_memory = "50G"

	statement = ''' 
	cp %(RMD_ROOT)s/Cluster.Rmd . &&
	R -e "rmarkdown::render('Cluster.Rmd', output_file='Cluster.html')" '''

	P.run(statement)

@follows(cluster, find_markers, cluster_rmarkdown)
def full():
	pass

def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))
