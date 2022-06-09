"""
===================
Pipeline integ_1_qc
===================


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
print(PARAMS)

@mkdir('scanpy_integ_1_qc.dir')
def topy():
    '''
    convert a file from ipynb to py
    '''
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'pipeline_integ_1_qc')

    statement = '''
    jupyter nbconvert --to python %(ROOT_DIR)s/h5.ipynb --output-dir=scanpy_integ_1_qc.dir
    '''

    P.run(statement)


@mkdir('scanpy_integ_1_qc_fig.dir')
@follows(topy)
def run():
    ROOT_DIR = 'scanpy_integ_1_qc.dir'

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