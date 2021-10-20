"""
===========
Pipeline visium spaceranger count
===========
Overview
==================
This pipeline peforms spaceranger count on multiple samples. Spaceranger takes a microscope slide image and FASTQ files from spaceranger mkfastq and performs alignment, 
tissue detection, fiducial detection, and barcode/UMI counting. The pipeline uses the Visium spatial barcodes to generate feature-spot matrices, determine clusters, 
and perform gene expression analysis. In addition to fresh-frozen tissues, spaceranger count is also used for FFPE and Targeted Gene Expression libraries. The pipeline 
uses the RNA-seq aligner STAR for fresh-frozen samples or a probe aligner algorithm for FFPE tissues. 
Outputs are delivered in BAM, MEX, CSV, HDF5, TIFF, PNG, JPEG and HTML formats.
See https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger for more details

=====
Configuration
-------------
The pipeline uses CGAT-core and CGAT-apps throughout the pipeline. 

Input files
-----------
The pipeline is ran using fastq files that follow the naming convention Read1: Name_S#_R1_001.fastq.gz
and read2: Name_S#_R2_001.fastq.gz
And brightfield images (jpeg) using the naming convention: SlideName.jpeg

==================
The output of running this pipeline is the generation of a counts matrix and initial analysis. 
The end user can use this as a base to then further develop project specific analyses.

Code
==================
"""
#IMPORT MODULES
#from ruffus import *
#import sys
#import os
#import re
#import sqlite3
#import glob

#import cgatcore.pipeline as P
#import cgatcore.experiment as E



# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


JUPYTER_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_quantcells","Jupyter")

# Determine the location of the input fastq files

try:
    PARAMS['data']
except NameError:
    DATADIR = "."
else:
    if PARAMS['data'] == 0:
        DATADIR = "."
    elif PARAMS['data'] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS['data']

SEQUENCESUFFIXES = ("*.fastq.gz",
		    "*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

############################################
# Run spaceranger
############################################

def spaceranger_count(infile, outfile):
    '''
    spaceranger count is run to make spatial counts matrix
    '''
    if "" in infile
    
    Module load spaceranger/1.2.2
    for x in number_of_samples:
        spaceranger count --id=sample_id \
            --transcriptome= genome_file \
            --fastqs=fastq_file \
            --sample=sample_id \
            --image=image_file \
            --slide=slide_name \
            --area=area_name \
            --localcores=6 \
            --localmem=25


