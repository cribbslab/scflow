"""
=================
Pipeline sampleqc
=================

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


@originate("QC.html")
def rmarkdown_stats(outfile):
    '''
    Runs a quality report in Rmarkdown to assess the statistics of single-cell
    experiments ran using pipeline_kb.py
    '''

    RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kb-sampleqc-1","Rmarkdown")

    job_memory = "50G"

    statement = '''
    cp %(RMD_ROOT)s/QC.Rmd . &&
    R -e "rmarkdown::render('QC.Rmd', output_file='QC.html')"
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
