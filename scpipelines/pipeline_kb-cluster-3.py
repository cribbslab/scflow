"""
===================
Pipeline clustering
===================

The pipeline follows pipeline_kb and pipeline_kb-filter-2
and should be ran in the same directory as these pipelines


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


