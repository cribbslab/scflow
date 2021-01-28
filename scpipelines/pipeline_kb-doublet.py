"""
================
Pipeline doublet
================

The pipeline should be  exectued within the directory that pipeline_kb.py
was previously ran within and the filtering pipeline should have already been exectured prior to
running this pipeline.
"""

from ruffus import *

import sys
import os
import re
import sqlite3
import glob

import cgatcore.pipeline as P
import cgatcore.experiment as E



# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

filtered_suffixes = "*_filtered_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])


@follows(mkdir("Doublet_Figures.dir"))
@transform(SEURAT_OBJECTS,
	regex("RDS_objects.dir/(\S+)_filtered_SeuratObject.rds"),
	r"RDS_objects.dir/\1_filtered_SeuratObject.rds")
def doublet_finder(infile, outfile):
    '''Runs doublet analysis and outputs a seurat object'''

    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

    file_name = os.path.basename(infile)
    sample = re.match(r'(\S+)_filtered_SeuratObject.rds', file_name).group(1)

    job_memory = "50G"

    statement = '''Rscript %(R_PATH)s/doublet_finder.R -i %(infile)s -s %(sample)s'''

    P.run(statement)


@originate("Doublets.html")
def rmarkdown_doublet(outfile):
    '''

    '''

    RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-filter-2","Rmarkdown")

    job_memory = "50G"

    statement = '''
    cp %(RMD_ROOT)s/Doublets.Rmd . &&
    R -e "rmarkdown::render('Doublets.Rmd', output_file='Doublets.html')"
    '''

    P.run(statement)


@follows(rmarkdown_doublet)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
