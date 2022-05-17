"""
=======================
Pipeline integ_3_filter
=======================

"""

from ruffus import *

import sys
import os

import cgatcore.pipeline as P
import cgatcore.experiment as E


# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


@mkdir('scanpy_integ_3_filter.dir')
def topy():
    '''
    convert a file from ipynb to py
    '''
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'pipeline_integ_3_filter')

    statement = '''
    jupyter nbconvert --to python %(ROOT_DIR)s/h5.ipynb --output-dir=scanpy_integ_3_filter.dir
    '''

    P.run(statement)


@mkdir('scanpy_integ_3_filter_fig.dir')
@follows(topy)
def run():
    ROOT_DIR = 'scanpy_integ_3_filter.dir'

    job_memory = PARAMS['memory']

    statement = '''
    python %(ROOT_DIR)s/h5.py --work_dir=./
    '''
    P.run(statement)


@follows(run)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))