"""
=================
Pipeline cluster
=================

The pipeline should be  exectued within the directory that pipeline_quantcells.py or pipeline_quantnuclei.py was previously ran within.

This pipeline runs the scanpy workflow according to the scanpy tutorial processing and clustering
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


@mkdir("Scanpy_QC.dir")
@transform("kallisto.dir/*/bus/genecount/genes.mtx",
           regex("kallisto.dir/(\S+)/bus/genecount/genes.mtx"),
           r"Scanpy_QC.dir/\1.hd5ad")
def preprocessing(infile, outfile):
    '''
    Runs a quality report in Rmarkdown to assess the statistics of single-cell
    experiments ran using pipeline_kb.py
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_cluster","python")

    job_memory = PARAMS['memory']

    barcode = infile.replace("mtx", "barcodes.txt")
    ec_path = infile.replace("genecount/genes.mtx", "output.bus/matrix.ec")
    gene_path = infile.replace("mtx", "genes.txt")


    statement = '''
    python %(PYTHON_ROOT)s/preprocessing.py --matrix %(infile)s --barcodes %(barcode)s
               --ec %(ec_path)s --genename %(gene_path)s --outfile %(outfile)s
    '''

    P.run(statement)


@mkdir("Scanpy_cluster.dir")
@transform(preprocessing,
           regex("Scanpy_QC.dir/(\S+).hd5ad"),
           r"Scanpy_cluster.dir/\1.html")
def run_jupyter(infile, outfile):
    """Run a jupyter notebook containing a scanpy workflow that implements clustering"""

    JUPYTER_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_cluster","python")

    name = infile.replace(".hd5ad", "")
    name = name.replace("Scanpy_QC.dir/", "")

    results_file = infile.replace("Scanpy_QC.dir/", "Scanpy_cluster.dir")
    results_file = results_file.replace(".hd5ad", "_cluster.hd5ad")

    statement = """cp %(JUPYTER_ROOT)s/run_scanpy.ipynb %(name)s_run_scanpy.ipynb
                   papermill   Scanpy_cluster.dir/%(name)s_run_scanpy.ipynb  \
                   Scanpy_cluster.dir/%(name)s_run_scanpy.ipynb -p h5_input %(infile)s results_file %(results_file)s &&
                   cd Scanpy_cluster.dir && jupyter nbconvert %(name)s_run_scanpy.ipynb --to html
                   """


@follows(run_jupyter)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
