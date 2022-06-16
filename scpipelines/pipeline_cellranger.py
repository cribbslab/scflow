"""
===================
Pipeline cellranger
===================


Overview
==================

This pipeline is a wrapper for cellranger pipeline

Usage
=====

Configuration
-------------

The pipeline uses CGAT-core and CGAT-apps throughout the pipeline. Please see installation
and setup and installation instructions at `cgat-core documentation <>`_


Input files
-----------

The pipeline is ran using fastq files that follow the naming convention Read1: Name.fastq.1.gz
and read2: Name.fastq.2.gz.

 * a fastq file (single /paired end (always paired end for drop seq methods)
 * a GTF geneset

For drop-seq, the default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains UMI/cellular barcode data and fastq.2.gz contains sequencing reads.
Chromium outputis of the format: samplename_R1.fastq.gz and samplename_R2.fastq.gz so will require conversion to the default file format above.

Pipeline output
==================


Code
==================

"""
from ruffus import *

import sys
import os
import re
import sqlite3
import glob

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as iotools



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

SEQUENCESUFFIXES = ("*.fastq.gz",
		    "*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

############################################
# Build indexes
############################################
# Need to build one index for cDNA and one for intron


@follows(mkdir("counts.dir"))
@transform(SEQUENCEFILES,
           regex("data.dir/(\S+).fastq.1.gz"),
           r"counts.dir/\1/cellranger.count.dummy")
def cellranger_count(infile, outfile):
    '''
    Execute the cell ranger pipleline for each sample.
    '''

    read2 = infile.replace(".1.gz",".2.gz")

    # set key parameters
    transcriptome = PARAMS["cellranger_transcriptome"]

    # set the maximum number of jobs for cellranger
    max_jobs = PARAMS["cellranger_maxjobs"]

    log_file = infile + ".log"

    id_tag = infile.replace(".fastq.1.gz", "")

    job_threads = 3
    job_memory = "24000M"
    statement ='''cellranger count
                   --id %(id_tag)s
                   --fastqs %(infile)s %(read2)s
                   --transcriptome %(transcriptome)s
                   --chemistry %(cellranger_chemistry)s
            &> %(log_file)s
        '''

    P.run(statement)

    iotools.touch_file(outfile)


############################################
# Perform read quality steps
############################################


@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           regex("(\S+).fastq.(\d).gz"),
           r"fastqc_pre.dir/\1.fastq.\2_fastqc.html")
def run_fastqc(infile, outfile):
    '''
    Fastqc is ran to determine the quality of the reads from the sequencer
    '''
    # paired end mode
    if "fastq.1.gz" in infile:
        second_read = infile.replace(".fastq.1.gz", ".fastq.2.gz")
        statement = "fastqc -q -o fastqc_pre.dir/ %(infile)s %(second_read)s"

    else:
        statement = "fastqc -q -o fastqc_pre.dir/ %(infile)s"

    P.run(statement)


@follows(cellranger_count)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
