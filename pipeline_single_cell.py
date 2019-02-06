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


# To do:
# merger reads if over multiple lanes

"""===========================
Pipeline single cell
===========================

Overview
========


This pipeline performs alignment free based quantification of drop-seq, 10X and smart-seq2
single-cell seqeucning analysis. Pseudoalignment is performed on the RNA reads,
using kallisto or Alevin and the resulting data is quantitatvely and qualitatively analysed.

The pipeline performs the following analyses:
* Alignment using kallisto or alevin (pert of salmon)
* QC of reads using the scater package


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

 * a fastq file (single /paired end?? - AC:require both (always paired end for drop seq methods and
potentially single end or paired end for smartseq2)
 * a GTF geneset

For chromium v2/ v3 CB and UMI barcodes are included in R1 files and R2 contains the raw sequencing reads.
The default file format assumes the following convention:
<samplename>.fastq.gz (fastq.1.gz (and fastq.2.gz for second read of paired data) are also accepted for raw reads)

Pipeline output
===============

The output of running this software is the generation of a html report, count matrices of gene expression.

Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import pandas as pd
from functools import reduce #you dont seem to use this anywhere in the pipeline?
import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleSC
import cgat.IndexedFasta as IndexedFasta # Again, doesnt seem to be used

import cgat.GTF as GTF
import cgatcore.iotools as iotools

import cgatpipelines.tasks.geneset as geneset # Doesnt seem to be used
import cgatpipelines.tasks.rnaseq as rnaseq # not used
import cgatpipelines.tasks.tracks as tracks # Not used
from cgatpipelines.report import run_report # Again not used

import cgatpipelines.tasks.expression as Expression # Not used

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
		    "*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


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

    job_memory = "12G" # I wouldnt set job limit unlimited because this is bad practice - it could break cluster
    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s;
    salmon index -k %(salmon_kmer)i %(salmon_index_options)s -t %(infile)s -i %(outfile)s
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


@originate("transcript2geneMap.tsv")
def getTranscript2GeneMap(outfile):
    ''' Extract a 1:1 map of transcript_id to gene_id from the geneset '''

    iterator = GTF.iterator(iotools.open_file(PARAMS['geneset']))
    transcript2gene_dict = {}

    for entry in iterator:

        # Check the same transcript_id is not mapped to multiple gene_ids!
        if entry.transcript_id in transcript2gene_dict:
            if not entry.gene_id == transcript2gene_dict[entry.transcript_id]:
                raise ValueError('''multipe gene_ids associated with
                the same transcript_id %s %s''' % (
                    entry.gene_id,
                    transcript2gene_dict[entry.transcript_id]))
        else:
            transcript2gene_dict[entry.transcript_id] = entry.gene_id

    with iotools.open_file(outfile, "w") as outf:
        outf.write("transcript_id\tgene_id\n")
        for key, value in sorted(transcript2gene_dict.items()):
            outf.write("%s\t%s\n" % (key, value))

# Input fastqc

# Pseudoalignment


if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s.(fastq.gz|fastq.1.gz)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))

    SEQUENCEFILES_KALLISTO_OUTPUT = (
        r"kallisto.dir/%s/output.bus" % (
            PARAMS["merge_pattern_output"].strip()))

    SEQUENCEFILES_SALMON_OUTPUT = (
        r"salmon.dir/%s/alevin/quants_mat.gz" % (
            PARAMS["merge_pattern_output"].strip()))

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.gz|fastq.1.gz)")

    SEQUENCEFILES_KALLISTO_OUTPUT = (
        r"kallisto.dir/\1/output.bus")

    SEQUENCEFILES_SALMON_OUTPUT = (
        r"salmon.dir/\1/alevin/quants_mat.gz")

############################################
# Perform read quality steps
############################################


@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           formatter(r"(?P<track>[^/]+).(?P<suffix>fastq.1.gz|fastq.gz)"),
           r"fastqc_pre.dir/{track[0]}.fastqc")
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

#############################
# Salmon- Alevin
#############################

@active_if(PARAMS['salmon_alevin'])
@follows(mkdir("salmon.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSalmonIndex, getTranscript2GeneMap),
         SEQUENCEFILES_SALMON_OUTPUT)
def runSalmonAlevin(infiles, outfile):
    '''
    Alevin is integrated with salmon to quantify and analyse 3' tagged-end
    single-cell sequencing data. Alevin supports 10Xv1, 10Xv2 and Drop-Seq
    sc technology.
    '''

    aligner = 'salmon_alevin'
    infiles = ModuleSC.check_multiple_read_files(infiles)
    fastqfile, index, t2gmap = infiles
    fastqfiles = ModuleSC.check_paired_end(fastqfile, aligner)
    if isinstance(fastqfiles, list):
        CB_UMI_fastq = " ".join(fastqfiles[0])
        reads_fastq = " ".join(fastqfiles[1]) 

    outfolder = outfile.rsplit('/',2)[0]
  
    statement = '''
    salmon alevin -l %(salmon_librarytype)s -1 %(CB_UMI_fastq)s -2  %(reads_fastq)s
    --%(salmon_sctechnology)s -i %(index)s -p %(salmon_threads)s -o %(outfolder)s
    --tgMap %(t2gmap)s --dumpCsvCounts
    '''

    job_memory = "30G"

    P.run(statement)

#############################
# Kallisto- Bus
#############################

@active_if(PARAMS['kallisto_bustools'])
@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex, getTranscript2GeneMap),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def runKallistoBus(infiles, outfile):
    '''
    Generates BUS files for single-cell sequencing

    infiles: raw sequencing fastq files, kallisto index

    '''
    aligner = 'kallisto_bus'
    infiles = ModuleSC.check_multiple_read_files(infiles)
    fastqfile, index, t2gmap = infiles
    fastqfiles = ModuleSC.check_paired_end(fastqfile, aligner)
    fastqfiles = " ".join(fastqfiles)

    outfolder = outfile.rsplit('/',1)[0]

    statement = '''
    kallisto bus -i %(index)s -o %(outfolder)s -x %(kallisto_sctechnology)s
    -t %(kallisto_threads)s %(fastqfiles)s
    '''

    P.run(statement)

######################
# Process bus file
######################

# Must have bustools installed
# https://github.com/BUStools/bustools

@active_if(PARAMS['kallisto_bustools'])
@transform(runKallistoBus,
           regex(r"kallisto.dir/(.*)/output.bus"),
           r"kallisto.dir/\1/\1_sorted.txt")
def busText(infile, outfile):
    '''
    Sort the bus file produced by kallisto and then convert it to a text file.
    '''

    tmp_bus  = P.get_temp_filename(".")
    E.warn("====================================================")
    E.warn(outfile)
    statement = '''
    bustools sort -o %(tmp_bus)s %(infile)s ;
    bustools text -o %(outfile)s %(tmp_bus)s
    '''

    P.run(statement)

#########################
# SCE object  
#########################

@follows(mkdir("R.dir"))
@active_if(PARAMS['salmon_alevin'])
@transform(runSalmonAlevin,
           regex(r"salmon.dir/(.*)/alevin/quants_mat.gz"),
           r"R.dir/\1.rds")
def readAlevinSCE(infile,outfile):
    '''
    Collates alevin count matrices for each sample
    Creates a single cell experiment class in R and saves as and r object
    '''

    working_dir = os.getcwd()
    sc_directory = PARAMS['sc_dir']
    script_loc = sc_directory + "/sce.r"
    
    job_memory = "10G"

    statement = '''
    Rscript %(script_loc)s -w %(working_dir)s -i %(infile)s -o %(outfile)s
    '''
    
    P.run(statement)


@follows(buildReferenceTranscriptome, buildSalmonIndex, buildKallistoIndex,
         getTranscript2GeneMap, runSalmonAlevin, runKallistoBus)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
