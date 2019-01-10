#!/usr/bin/env python

"""===========================
Pipeline single cell
===========================

Overview
========

This pipeline was developed to perform mapping of sequencing data obtained from single cell techniques
including DropSeq, 10X, celseq and gemcode. Pseudoalignment is performed on the RNA reads,
using kallisto or Alevin and the resulting data is quantitatvely and qualitatively analysed.



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

# Determine the location of the input fastq files
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
		    "*.fastq.1.gz",
		    "*.fastq.2.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.[1-2].gz|fastq.gz)")

############################################
# Build indexes
############################################

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
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    '''
    Builds a salmon index for the reference transriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    salmon_kmer: int
       :term: `PARAMS` kmer size for sailfish.  Default is 31.
       Salmon will ignores transcripts shorter than this.
    salmon_index_options: str
       :term: `PARAMS` string to append to the salmon index command to
       provide specific options e.g. --force --threads N
    outfile: str
       path to output file
    '''

    job_memory = "unlimited"
    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s;
    salmon index -k %(salmon_kmer)i %(salmon_index_options)s -t %(infile)s -i %(outfile)s
    -k %(salmon_kmer)s
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
# Does this apply??? Alevin and bustools give generic output, change to
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s.(fastq.gz|fastq.1.gz|fastq.2.gz)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))

    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"kallisto.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"salmon.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.gz|fastq.1.gz|fastq.2.gz)")

    # Need to run bustools to find exact output files
    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/\1/output.bus",
        r"kallisto.dir/\1/matrix.mtx"]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/\1/quants_mat.gz",
        r"salmon.dir/\1/quants_mat_cols.txt",
        r"salmon.dir/\1/quants_mat_rows.txt",
        r"salmon.dir/\1/quants_tier_mat.gz",
        r"salmon.dir/\1/quants_tier_mat.gz"]

# Alevin
# Count matrix, multiple samples? run seperately??? Gene by cell, so sample separate matrix?
@active_if(salmon_alevin)
@follows(mkdir("salmon.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSalmonIndex, getTranscript2GeneMap),
         SEQUENCEFILES_SALMON_OUTPUT)
def runSalmonAlevin(infiles, outfile):
    '''
    Alevin is integrated with salmon to quantify and analyse 3' tagged-end
    single-cell sequencing data. Alevin supports 10Xv1, 10Xv2 and Drop-Seq
    sc technology
    '''

    # Probably need to separate sequencing files
    sequence_files, salmon_index, t2gmap = infiles
    statement = '''
    salmon alevin -l %(salmon_librarytype)s -1 CB_UMI_sequences?? -2  %(sequence_files)s
    --%(salmon_sctechnology)s -i %(salmon_index)s -p %(salmon_threads)s -o salmon.dir
    --tgMap %(t2gmap)s --dumpCsvCounts
    '''

# BUStools approach
@active_if(kallisto_bustools)
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
@active_if(kallisto_bustools)
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


# Count

# Quality control

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

# Scater (levin swing???)

## Multi QC
## Generate our own multi QC report using R

# Pseudotime
# clustering
# Velocyte
# Cell cycle (cyclone), blocking

# Create R data object
