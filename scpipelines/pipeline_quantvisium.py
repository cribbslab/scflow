"""
===========
Pipeline visium spaceranger count
===========
Overview
==================
This pipeline peforms spaceranger count on multiple
samples.

Spaceranger takes a microscope slide image and FASTQ
files from spaceranger mkfastq and performs alignment, 
tissue detection, fiducial detection, and barcode/UMI
counting. 

The pipeline uses the Visium spatial barcodes
to generate feature-spot matrices, determine clusters, 
and perform gene expression analysis. 

In addition to fresh-frozen tissues, spaceranger count
is also used for FFPE and Targeted Gene Expression libraries.

The pipeline uses the RNA-seq aligner STAR for fresh-frozen
samples or a probe aligner algorithm for FFPE tissues. 
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
from ruffus import *
import sys
import os
import re
import pandas as pd
import cgatcore.pipeline as P
import cgatcore.experiment as E

# Load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# read the desin file and save to global dicts
DESIGN_DICT = pd.read_table("design.tsv", header=0, index_col=0).to_dict()


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

SEQUENCESUFFIXES = ("*_R1_001.fastq.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

############################################
# Run spaceranger
############################################

@transform(SEQUENCEFILES,
           regex("(\S+)_R1_001.fastq.gz"),
           r"/\1/outs/filtered_feature_bc_matrix/")
def spaceranger_count(infile, outfile):
    '''
    spaceranger count is run to make spatial counts matrix
    '''

    name = infile.replace("_R1_001.fastq.gz", "")
    name = name.replace(PARAMS['data']+ "/", "")


    image_loc = DESIGN_DICT["image_loc"][name]
    slide_name = DESIGN_DICT["slide_name"][name]
    slide_area = DESIGN_DICT["slide_area"][name]

    statement = """module load %(spaceranger_module)s &&
                   spaceranger count --id=%(name)s 
                    --transcriptome=%(trans_file)s
                    --fastqs=data.dir
                    --sample=%(name)s
                    --slide=%(slide_name)s
                    --image=%(image_loc)s
                    --area=%(slide_area)s
                    --localcores=5
                    --localmem=20"""


    P.run(statement, job_memory="20G", job_threads=5)



@follows(spaceranger_count)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
