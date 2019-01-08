#!/usr/bin/env python

"""===========================
Pipeline single cell
===========================

Overview
========

This pipeline was developed to perform mapping of sequencing data obtained from single cell techniques
including DropSeq, 10X and (...). Pseudoalignment is performed on the RNA reads, using kallisto or Alevin
and the resulting data is quantitatvely and qualitatively analysed.



Requires:
 * a fastq file (single/paired end??)
 * an indexed genome
 *

Kallisto is faster than Alevin, however Kallisto and BUStools is not yet supported for 10X V3 data.

Pipeline output
===============

The output of running this software is the generation of a html report.

Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import pandas as pd
import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleSC
import cgat.IndexedFasta as IndexedFasta

import cgat.GTF as GTF
import cgatcore.iotools as iotools

import cgatpipelines.tasks.geneset as geneset
import cgatpipelines.tasks.rnaseq as rnaseq
import cgatpipelines.tasks.tracks as tracks
from cgatpipelines.report import run_report

import cgatpipelines.tasks.expression as Expression

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

# some files are fastq.2.gz etc, do we need * or just 1 like in rnaseqdiffexpression???
SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.*.gz|fastq.gz|sra)")

############################################
# Quality control of the fastq files
############################################

@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           suffix(".fastq.gz"),
           r"fastqc_pre.dir/\1.fastq")
def fastqc_pre(infile, outfile):
    """
    Runs fastQC on each input file
    """

    statement = "fastqc -q -o fastqc_pre.dir/ %(infile)s"

    P.run(statement)

@follows(fastqc_pre)
@follows(mkdir("processed.dir"))
@transform(SEQUENCEFILES,
           suffix(".fastq.gz"),
           r"processed.dir/\1_processed.fastq.gz")
def process_reads(infile, outfile):
    """
    Runs trimmomatic quality related trimming
    """

    if PARAMS["trimmomatic_run"]:

        trimmomatic_options = PARAMS["trimmomatic_options"]

        trimmomatic_options = "ILLUMINACLIP:%s:%s:%s:%s" % (
            PARAMS["trimmomatic_adapter"],
            PARAMS["trimmomatic_mismatches"],
            PARAMS["trimmomatic_p_thresh"],
            PARAMS["trimmomatic_c_thresh"]) + "\t" + trimmomatic_options

        phred = PARAMS["trimmomatic_phred"]

        ModuleTrna.process_trimmomatic(infile, outfile, phred,
                                   trimmomatic_options)
    else:

        statement = "cp %(infile)s %(outfile)s"

        P.run(statement)

@follows(mkdir("fastqc_post.dir"))
@transform(process_reads,
           regex("processed.dir/(\S+)_processed.fastq.gz"),
           r"fastqc_post.dir/\1.fastq")
def fastqc_post(infile, outfile):
    """
    Runs fastQC on each of the processed files
    """

    statement = """fastqc -q -o fastqc_post.dir/ %(infile)s
                """

    P.run(statement)

############################################
# Build indexes
############################################

# Index for kallisto taken from rnaseqdiffexpression pipeline
@mkdir('geneset.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           r"geneset.dir/\1.fa")
def buildReferenceTranscriptome(infile, outfile):
    '''
    Builds a reference transcriptome from the provided GTF geneset - generates
    a fasta file containing the sequence of each feature labelled as
    "exon" in the GTF.
    --fold-at specifies the line length in the output fasta file
    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat %(infile)s |
    awk '$3=="exon"'|
    cgat gff2fasta
    --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
    --log=%(outfile)s.log > %(outfile)s;
    samtools faidx %(outfile)s
    '''

    P.run(statement)

@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".kallisto.index")
def buildKallistoIndex(infile, outfile):
    '''
    Builds a kallisto index for the reference transcriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    kallisto_kmer: int
       :term: `PARAMS` kmer size for Kallisto.  Default is 31.
       Kallisto will ignores transcripts shorter than this.
    outfile: str
       path to output file
    '''

    job_memory = "12G"

    statement = '''
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s %(infile)s
    '''

    P.run(statement)

# Input fastqc

# Pseudoalignment
# BUStools approach
@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex, getTranscript2GeneMap),
         r"kallisto.dir/output.bus")
def runKallistoBus(infiles, outfile):
    '''
    Generates BUS files for single-cell sequencing

    infiles: raw sequencing fastq files, kallisto index, transcript2genemap

    # Probably need to separate sequencing files
    '''
    sequence_files, kallisto_index, t2gmap = infiles
    # Unsure what to call output
    statement = '''
    kallisto bus -i %(kallisto_index)s -o kallisto.dir -x %(kallisto_sctechnology)s
    -t %(kallisto_threads)s %(sequence_files)
    '''

    P.run(statement)

######################
# Process bus file
######################

# Must have bustools installed, see https://github.com/BUStools/bustools

@transform(runKallistoBus,
           suffix(".bus"),
           r"kallisto.dir/\1_sorted.txt")
def busText(infile, outfile):
    '''
    Sort the bus file produced by kallisto and then convert it to a text file.
    '''

    tmp_bus  = P.get_temp_filename(".")

    statement = '''
    bustools sort -o %(tmp_bus)s %(infile)s |
    bustools text -o %(outfile)s tmp_bus
    '''

## Alevin

# Count

# Quality control
# Scater (levin swing???)

## Multi QC
## Generate our own multi QC report using R

# Pseudotime
# clustering
# Velocyte
# Cell cycle (cyclone), blocking

# Create R data object
