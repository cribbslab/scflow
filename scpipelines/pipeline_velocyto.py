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
                   --genomeFastaFiles %(genome)s"""

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


@transform(fastq_to_bam,
           regex("(\S+)_unmapped.bam"),
           r"data.dir/\1_tagged_Cell.bam")
def cell_barcode_bam(infile, outfile):
    """
    Extracts bases from the cell barcode encoding read
    (BARCODED_READ), and creates a new BAM tag with those bases on the genome read.
    """

    name = infile.replace("_unmapped.bam","")

    statement = """TagBamWithReadSequenceExtended 
                   INPUT=%(infile)s
                   OUTPUT=%(outfile)s 
                   SUMMARY=%(name)s_tagged_Cellular.bam_summary.txt 
                   BASE_RANGE=1-12 
                   BASE_QUALITY=10 
                   BARCODED_READ=1 
                   DISCARD_READ=False 
                   TAG_NAME=XC 
                   NUM_BASES_BELOW_QUALITY=1 """

    P.run(statement)

@transform(cell_barcode_bam,
           regex("(\S+)_tagged_Cell.bam"),
           r"data.dir/\1_tagged_CellMolecular.bam")
def molecular_barcode_bam(infile, outfile):
    """

    """

    name = infile.replace("_tagged_Cell.bam", "")

    statement = """TagBamWithReadSequenceExtended 
                   INPUT=%(input)s 
                   OUTPUT=%(outfile)s 
                   SUMMARY=%(name)s_tagged_Molecular.bam_summary.txt 
                   BASE_RANGE=13-20 
                   BASE_QUALITY=10 
                   BARCODED_READ=1 
                   DISCARD_READ=True 
                   TAG_NAME=XM 
                   NUM_BASES_BELOW_QUALITY=1 """

    P.run(statement)


@transform(molecular_barcode_bam,
           regex("(\S+)_tagged_CellMolecular.bam"),
           r"data.dir/\1_tagged_filtered.bam")
def filter_bam(infile, outfile):
    """
    filter the bam file to regect XQ tag
    """

    statement = """FilterBam 
                   TAG_REJECT=XQ 
                   INPUT=%(infile)s 
                   OUTPUT=%(outfile)s"""

    P.run(statement)

@transform(filter_bam,
           regex("(\S+)_tagged_filtered.bam"),
           r"data.dir/\1_tagged_trimmed_smart.bam")
def trim_starting_sequence(infile, outfile):
    """
    Trim the starting sequence of each read in the bamfile
    """

    name = infile.replace("_tagged_filtered.bam", "")

    statement = """TrimStartingSequence 
                   INPUT=%(infile)s 
                   OUTPUT=%(outfile)s 
                   OUTPUT_SUMMARY=%(name)s_adapter_trimming_report.txt 
                   SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG 
                   MISMATCHES=0 
                   NUM_BASES=5"""

    P.run(statement)


@transform(trim_starting_sequence,
           regex("(\S+)_tagged_trimmed_smart.bam"),
           r"data.dir/\1_polyA_filtered.bam")
def polyA_trimmer(infile, outfile):
    """
    Remove the poly A tail from each read
    """

    statement = """PolyATrimmer 
                   INPUT=unaligned_tagged_trimmed_smart.bam 
                   OUTPUT=unaligned_mc_tagged_polyA_filtered.bam 
                   OUTPUT_SUMMARY=polyA_trimming_report.txt 
                   MISMATCHES=0 
                   NUM_BASES=6 
                   USE_NEW_TRIMMER=true"""

    P.run(statement)

@follows(mkdir("dropest.dir"))
@transform(trim_starting_sequence,
           regex("(\S+)_ployA_filtered.bam"),
           add_inputs(PARAMS['geneset']),
           r"dropest.dir/\1_dropEst.rds")
def run_dropest(infiles, outfile):
    """
    runs dropEst on the bam file and generated an rds file as output
    """
    bamfile, gtffile = infiles
    config = PARAMS['dropest_config']

    statement = """dropest -m -V -b -f 
                   -g %(gtffile)s
                   -o dropEst_out 
                   -L eiEIBA 
                   -c %(config)s 
                   %(bamfile)s"""

    P.run(statement)


@follows(mkdir("star.dir"))
@transform(trim_starting_sequence,
           regex("(\S+)_ployA_filtered.bam"),
           add_inputs(PARAMS['geneset'],
                      PARAMS['genome']),
           r"star.dir/\1_mapped.bam")
def star_mapping(infile, outfile):
    """
    Perform star mapping
    """

    bamfile, gtffile, genome = infiles
    root, dirs, files = os.walk(outfile)

    statement = """STAR 
                   --runThreadN 12 
                   --runMode genomeGenerate 
                   --genomeDir %(root)s
                   --genomeFastaFiles %(genome)s
                   --sjdbGTFfile %(gtffile)s"""

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
