"""
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


@originate("Filter.html")
def rmarkdown_stats(outfile):
    '''
    Runs a quality report in Rmarkdown to assess the statistics of single-cell
    experiments ran using pipeline_kb.py
    '''

    RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-filter-2","Rmarkdown")

    job_memory = "50G"

    statement = '''
    cp %(RMD_ROOT)s/Filter.Rmd . &&
    R -e "rmarkdown::render('Filter.Rmd', output_file='Filter.html')"
    '''

    P.run(statement)

@follows(rmarkdown_stats)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
