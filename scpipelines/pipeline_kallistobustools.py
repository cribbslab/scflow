##############################################################################
#
#   Botnar Resaerch Centre
#
#   $Id$
#
#   Copyright (C) 2018 Adam Cribbs
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""
====================
Pipeline single cell
====================


Overview
==================

This pipeline performs alignment free based quantification of single-nuclei-seq.
The pipeline is based off a biostars response by lior pachter:
https://www.biostars.org/p/397671/.

The nuclei data will be processed using kalliso and a custom built DNA and intron
index for the species of interest. Two matrices will be generated: one for
spliced transcripts and one for unspliced transcripts, which will be summed
to obtain the total nuclear transcripts.

To learn how to generate a cDNA and intron index see the following tutorial:
https://www.kallistobus.tools/velocity_index_tutorial.html.

Important: The mouse cDNA and intron index is about 26GB. Because of this,
building it and processing data with it requires significantly more RAM than
typical kallisto workflows, and we recomend using a machine with at least
64GB RAM for this workflow.

The downstream analysis of the data is then almost identical to the singlecell
pipeline.

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

The output of running this pipeline is the generation of an
Rmarkdown report with summary statistics. The end user can use this as a base
to then further develop project specific analyses.

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

import pandas as pd

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


@mkdir('geneset.dir')
@originate("geneset.dir/index.idx")
def build_kallisto_index(outfile):
    '''
    Builds a kallisto index for the reference transcriptome
    Parameters
    ----------
    '''


    job_memory = "65G"

    statement = '''
    kb ref -i geneset.dir/index.idx -g geneset.dir/t2g.txt -f1 geneset.dir/cdna.fa
    -f2 geneset.dir/intron.fa -c1 geneset.dir/cdna_t2c.txt -c2 geneset.dir/intron_t2c.txt
    --workflow %(kallisto_workflow)s  %(genome_file)s %(geneset)s 2> ref.log
    '''

    P.run(statement)


############################################
# Perform read quality steps
############################################


@follows(build_kallisto_index)
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


############################################
# Pseudoalignment
############################################

if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s.(fastq.gz|fastq.1.gz)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))

    SEQUENCEFILES_KALLISTO_OUTPUT = (
        r"kallisto.dir/%s/bus/output.bus" % (
            PARAMS["merge_pattern_output"].strip()))

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.gz|fastq.1.gz)")

    SEQUENCEFILES_KALLISTO_OUTPUT = (
        r"kallisto.dir/\1/bus/output.bus")


#############################
# Kallisto- Bus
#############################


@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(build_kallisto_index),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def run_kallisto_bus(infiles, outfile):
    '''
    Generates BUS files for single-cell sequencing

    infiles: raw sequencing fastq files, kallisto index

    '''

    fastqfile, index = infiles[0]
    glob_search = "geneset.dir/index*"
    index_files = glob.glob(glob_search)
    index_files = ",".join(index_files)

    read2 = fastqfile.replace(".fastq.1.gz",".fastq.2.gz")
    fastqfiles = " ".join([fastqfile, read2])

    outfolder = outfile.rsplit('/',1)[0]

    statement = '''
    kb count -i %(index_files)s -g geneset.dir/t2g.txt
    -c1 geneset.dir/cdna_t2c.txt -c2 geneset.dir/intron_t2c.txt -x %(kallisto_sctechnology)s
    -o %(outfolder)s --workflow %(kallisto_workflow)s --%(kallisto_output_format)s  %(fastqfiles)s
    2> %(outfolder)s_kblog.log
    '''

    job_memory = '100G'

    P.run(statement)


@transform(run_kallisto_bus,
           regex("kallisto.dir/(\S+)/bus/output.bus"),
           r"kallisto.dir/\1/bus/genecount/genes.mtx")
def merge_mtx(infile, outfile):
    '''
    merge the spliced and unspliced matrix
    '''

    outpath = outfile.replace("genes.mtx", "")
    inpath = outfile.replace("/genecount/genes.mtx", "/counts_unfiltered/")

    if not os.path.exists(outpath):
        os.mkdir(outpath)


    PYTHON_PATH =  os.path.join(os.path.dirname(__file__), "python/")

    statement = """python %(PYTHON_PATH)sMergeSplicedMatrix.py -o %(outpath)s -s %(inpath)sspliced.mtx 
                   -c %(inpath)sspliced.barcodes.txt -a %(inpath)sspliced.genes.txt -u %(inpath)sunspliced.mtx
                   -b %(inpath)sunspliced.barcodes.txt -g %(inpath)sunspliced.genes.txt"""

    P.run(statement)


@follows(merge_mtx)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
