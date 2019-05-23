"""
=================
Pipeline velocyto
=================


Overview
==================

This pipeline performs velocyto for dropseq methods. First the samples are mapped
using dropseq tools from mccarrol lab and then the velocyto command line tools are
used to process the data into a loom file then the Rmarkdown generated visualisations
of the data.


Usage
=====

Configuration
-------------

The pipeline uses CGAT-core and CGAT-apps throught the pipeline. Please see installation
and setup and installation instructions at `cgat-core documentation <>`_


Input files
-----------

The pipeline is ran using fastq files that follow the naming convention Read1: Name.fastq.1.gz
and read2: Name.fastq.2.gz. 

 
Pipeline output
==================


Code
==================

"""
from ruffus import *

import sys
import os
import sqlite3
import pandas as pd

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgat.GTF as GTF
import cgatcore.iotools as iotools

# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Determine the location of the input fastq files
SEQUENCESUFFIXES = ("*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


############################################
# Build indexes
############################################


@active_if(PARAMS['star_index'] == 1)
@mkdir('star_index.dir')
@transform(PARAMS['geneset'],
           suffix(".gtf.gz"),
           add_inputs(PARAMS['genome']),
           "star_index.dir/Genome")
def index_genome_star(infiles, outfile):
    """
    Generate a star index if no index is supplied
    """

    gtffile, genome = infiles

    statement = """STAR
                   --runThreadN 12 
                   --runMode genomeGenerate 
                   --genomeDir star_index.dir/
                   --genomeFastaFiles %(genome)s
                   --sjdbGTFfile %(gtffile)s"""

    P.run(statement)

@follows(PARAMS['star_index'] == 0)
@mkdir('star_index.dir')
@originate('star_index.dir/Genome')
def index_genome_copy(outfile):

    """
    If a star index has already been generated then it will be copied
    to the star_index.dir directory for tracability
    """

    star_location = PARAMS['star_index_dir']
    statement = """cp -R %(star_location)s star_index.dir/ """

    P.run(statement)

@follows(mkdir('data.dir'))
@transform(SEQUENCEFILES,
           formatter(r"(?P<track>[^/]+).(?P<suffix>fastq.1.gz)"),
           r"data.dir/{track[0]}_unmapped.bam")
def fastq_to_bam(infile, outfile):
    """
    generate an unmapped bam file for further parsing into dropseq tools
    """

    second_read = infile.replace(".fastq.1.gz", ".fastq.2.gz")
    name = infile.replace(".fastq.1.gz", "")

    statement = """picard FastqToSam 
                   F1=%(infile)s 
                   F2=%(second_read)s 
                   O=%(outfile)s 
                   SM=%(name)s"""

    P.run(statement)


@follows()
def quant():
    pass


@follows()
def loom_generation():
    pass


@follows()
def velocyto():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
