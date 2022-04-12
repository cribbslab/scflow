2"""
===============
Pipeline filter
===============

The pipeline should be  exectued within the directory that pipeline_kb.py
was previously ran within.
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


@originate("Filter.html")
def rmarkdown_stats(outfile):
    '''
    Runs a quality report in Rmarkdown to assess the statistics of single-cell
    experiments ran using pipeline_kb.py
    '''

    RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_filter-2","Rmarkdown")

    job_memory = PARAMS['mem']

    statement = '''
    cp %(RMD_ROOT)s/Filter.Rmd . &&
    R -e "rmarkdown::render('Filter.Rmd', output_file='Filter.html')"
    '''

    P.run(statement)


RDSFILES = os.path.join("RDS_objects.dir", "*")

@follows(mkdir("Anndata.dir"))
@follows(rmarkdown_stats)
@transform(RDSFILES,
           regex("RDS_objects.dir/(\S+)_filtered_SeuratObject"),
           r"Anndata.dir/\1.h5ad")
def seurat2anndata(infile, outfile):
    '''
    convert from seurat to anndata 
    '''

    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    job_memory = "20G"

    statement = '''
    Rscript %(R_ROOT)s/seurat2anndata.R --input=%(infile)s
    '''

    P.run(statement)


@follows(seurat2anndata)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
