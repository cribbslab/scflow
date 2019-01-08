#!/usr/bin/env python

"""
ModuleSC.py - Tasks for running single cell pipleine

"""
import os
import re
import pysam
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import cgatcore.pipeline as P
import cgatcore.database as Database
import pandas as pd
