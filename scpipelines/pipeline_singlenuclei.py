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

import cgatcore.pipeline as P
import cgatcore.experiment as E
import scpipelines.ModuleSC as ModuleSC

import pandas as pd

import cgat.GTF as GTF
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
@originate("geneset.dir/index.idx.0")
def build_kallisto_index(infiles, outfile):
    '''
    Builds a kallisto index for the reference transcriptome
    Parameters
    ----------
    '''


    job_memory = "65G"

    statement = '''
    kb ref -i geneset.dir/index.idx -g geneset.dir/t2g.txt -f1 geneset.dir/cdna.fa
    -f2 geneset.dir/intron.fa -c1 geneset.dir/cdna_t2c.txt -c2 geneset.dir/intron_t2c.txt
    --workflow nucleus -n 8 %(genome_file)s %(geneset)s
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
    read2 = fastqfile.replace(".fastq.1.gz",".fastq.2.gz")
    fastqfiles = " ".join([fastqfile, read2])

    outfolder = outfile.rsplit('/',1)[0]

    statement = '''
    kallisto bus -i %(index)s -o %(outfolder)s -x %(kallisto_sctechnology)s
    -t %(kallisto_threads)s %(fastqfiles)s
    '''

    job_memory = '30G'

    P.run(statement)


@transform(run_kallisto_bus,
           regex("kallisto.dir/./(\S+)/bus/output.bus"),
           r"kallisto.dir/\1/bus/output.sort.whitelist.bus")
def bustools_sort_for_whitelist(infile, outfile):
    '''use bustools sort to sort the corrected bus record '''

    statement = '''
    bustools sort -o %(outfile)s  -t 4 %(infile)s
    '''

    P.run(statement)


@collate(bustools_sort_for_whitelist,
         regex("kallisto.dir/(\S+)/bus/output.sort.whitelist.bus"),
         r"\1_whitelist.txt")
def whitelist(infile, outfile):
    '''use umitools to generate whitelist of barcodes'''

    infile = infile[0]

    statement = '''
                bustools whitelist -o %(outfile)s %(infile)s
    '''

    P.run(statement)


@transform(run_kallisto_bus,
           regex("kallisto.dir/(\S+)/bus/output.bus"),
           add_inputs(whitelist),
           r"kallisto.dir/\1/bus/output.correct.bus")
def bustools_correct(infiles, outfile):
    ''' then
    bustools correct function'''

    bus_file = infiles[0]

    # pick up correct whitelist
    match = bus_file.replace("kallisto.dir/", "")
    match = match.replace("/bus/output.bus","")
    whitelist = match + "_whitelist.txt"


    statement = '''

    bustools correct -w %(whitelist)s -o %(outfile)s  %(bus_file)s
    '''

    P.run(statement)


@transform(bustools_correct,
           regex("kallisto.dir/(\S+)/bus/output.correct.bus"),
           r"kallisto.dir/\1/bus/output.sort.bus")
def bustools_sort(infile, outfile):
    '''use bustools sort to sort the corrected bus record '''

    statement = '''
    bustools sort -o %(outfile)s  -t 4 %(infile)s
    '''

    P.run(statement)


@transform(bustools_sort,
           regex("kallisto.dir/(\S+)/bus/output.sort.bus"),
           add_inputs(introns_transcripts_no_version),
           r"kallisto.dir/\1/bus/introns_capture.bus")
def bustools_capture_intron(infiles, outfile):
    '''use bustools capture for cDNA '''

    infile, capture_list = infiles

    matrix = infile.replace("output.sort.bus","matrix.ec")
    trans = infile.replace("output.sort.bus","transcripts.txt")

    statement = '''
    bustools capture -s -o %(outfile)s -c %(capture_list)s  -e %(matrix)s -t %(trans)s  %(infile)s
    '''

    P.run(statement)

# Bustools capture cDNA and then introns
@transform(bustools_sort,
           regex("kallisto.dir/(\S+)/bus/output.sort.bus"),
           add_inputs(capture_list),
           r"kallisto.dir/\1/bus/cDNA_capture.bus")
def bustools_capture_cdna(infiles, outfile):
    '''use bustools capture for cDNA '''

    infile, capture_list = infiles

    matrix = infile.replace("output.sort.bus","matrix.ec")
    trans = infile.replace("output.sort.bus","transcripts.txt")

    statement = '''
    bustools capture -s -o %(outfile)s -c %(capture_list)s  -e %(matrix)s -t %(trans)s  %(infile)s
    '''

    P.run(statement)

@transform(bustools_capture_intron,
           regex("kallisto.dir/(\S+)/bus/introns_capture.bus"),
           r"kallisto.dir/\1/bus/unspliced/unspliced.mtx")
def bustools_count_intron(infile, outfile):

    outfile = outfile.replace(".mtx","")

    matrix = infile.replace("introns_capture.bus","matrix.ec")
    trans = infile.replace("introns_capture.bus","transcripts.txt")

    statement = '''
    bustools count -o %(outfile)s -g kallisto.dir/cDNA_introns_t2g.txt -e %(matrix)s -t %(trans)s --genecounts %(infile)s
    '''

    P.run(statement)


@transform(bustools_capture_cdna,
           regex("kallisto.dir/(\S+)/bus/cDNA_capture.bus"),
           r"kallisto.dir/\1/bus/spliced/spliced.mtx")
def bustools_count_cdna(infile, outfile):

    outfile = outfile.replace(".mtx","")

    matrix = infile.replace("cDNA_capture.bus","matrix.ec")
    trans = infile.replace("cDNA_capture.bus","transcripts.txt")

    statement = '''
    bustools count -o %(outfile)s -g kallisto.dir/cDNA_introns_t2g.txt -e %(matrix)s -t %(trans)s --genecounts %(infile)s
    '''

    P.run(statement)


@merge([bustools_count_intron, bustools_count_cdna],
    "geneset.dir/introns_t2g.txt")
def merge_matrices(infiles, outfile):
    '''use python script to merge the spliced and unspliced matrix'''

    unspliced, spliced = infiles
    unspliced_barcode = unspliced.replace(".mtx",".barcodes.txt")
    unspliced_genes = unspliced.replace(".mtx",".genes.txt")

    spliced_barcode = spliced.replace(".mtx",".barcodes.txt")
    spliced_genes = spliced.replace(".mtx",".genes.txt")

    bus_dir = spliced.replace("spliced/spliced.mtx","")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python")

    statement = '''python %(PYTHON_ROOT)s/MergeSplicedMatrix.py -o %(outfile)s -d %(bus_dir)s -s %(spliced)s
                   -c %(spliced_barcode)s -t %(spliced_genes)s -u %(unspliced)s -b %(spliced_barcode)s -g %(spliced_genes)s'''

    P.run(statement)


#########################
# Multiqc
#########################

@follows(mkdir("MultiQC_report.dir"))
@follows(run_fastqc, bustools_count_cdna)
@originate("MultiQC_report.dir/multiqc_report.html")
def build_multiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)


@follows(merge_matrices, build_multiqc)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
