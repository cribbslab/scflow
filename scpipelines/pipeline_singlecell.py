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

This pipeline performs alignment free based quantification of drop-seq, 10X and smart-seq2
single-cell sequencing analysis. Pseudoalignment is performed on the RNA reads,
using kallisto or Alevin and the resulting data is quantitatvely and qualitatively analysed.

The pipeline performs the following analyses:
* Alignment using kallisto or alevin (part of salmon)
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

 * a fastq file (single /paired end (always paired end for drop seq methods and
potentially single end or paired end for smartseq2)
 * a GTF geneset

The default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains UMI/cellular barcode data and fastq.2.gz contains sequencing reads.
Chromium outputis of the format: samplename_R1.fastq.gz and samplename_R2.fastq.gz so will require conversion to the default file format above.

Pipeline output
===============

The output of running this software is the generation of a SingleCellExperiment object and further downstream analysis including: clustering, pseudotime analysis, velocity time graphs, quality control analysis.

Code
====

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

R_ROOT = os.path.join(os.path.dirname(__file__), "scpipelines", "pipeline_singlecell","R")

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


@follows(run_fastqc)
@follows(mkdir("fastq_post.dir"))
@transform(SEQUENCEFILES,
           regex("(\S+).fastq.(\d).gz"),
           [r"fastq_post.dir/\1.fastq.1.gz",
            r"fastq_post.dir/\1.fastq.2.gz"])
def run_fastp(infile, outfiles):
    '''
    Fastp will trim the quality of the reads to improve mappability
    '''
    out_first = outfiles[0]
    out_second = outfiles[1]

    report_out = out_first.replace(".fastq.1.gz", ".html")

    # paired end mode
    if "fastq.1.gz" in infile:
        second_read = infile.replace(".fastq.1.gz", ".fastq.2.gz")

        statement = "fastp -i %(infile)s -I %(second_read)s -o %(out_first)s -O %(out_second)s -h %(report_out)s  %(fastp_options)s"

        P.run(statement)


############################################
# Build indexes
############################################


@mkdir('geneset.dir')
@merge([PARAMS['prim_trans1'], PARAMS['prim_trans2']],
           r"geneset.dir/geneset_all.fa")
def buildReferenceSalmon(infiles, outfile):
    '''
    Builds a reference transcriptome and decoy sequneces for alevin and kallisto
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
    prim_trans1, prim_trans2 = infiles
    genome_file1 = PARAMS['genome1']

    if PARAMS['mixed_species']:
        genome_file2 = PARAMS['genome2']
        tmp1_g = P.get_temp_filename('.')
        tmp2_g = P.get_temp_filename('.')
        tmp_trans1 = P.get_temp_filename('.')
        tmp_trans2 = P.get_temp_filename('.')
        statement = '''
                    guzip %(genome_file1)s > %(tmp1_g)s &&
                    gunzip %(genome_file2)s > %(tmp2_g)s &&
                    grep "^>" %(tmp1_g)s %(tmp2_g)s | cut -d " " -f 1 > decoys.txt &&
                    sed -i.bak -e 's/>//g' decoys.txt &&
                    zcat %(prim_trans1)s  > %(tmp_trans1)s &&
                    zcat %(prim_trans2)s  > %(tmp_trans2)s &&
                    cat %(tmp_trans1)s %(tmp_trans2)s %(tmp1_g)s %(tmp2_g)s > %(outfile)s
                    '''
    else:
        statement = '''
                       grep "^>" <(gunzip -c %(genome_file1)s) | cut -d " " -f 1 > decoys.txt &&
                       sed -i.bak -e 's/>//g' decoys.txt &&
                       cat %(prim_trans1)s %(genome_file1)s > %(outfile)s
                       '''

    P.run(statement)

    if PARAMS['mixed_species']:
        os.unlink(tmp1_g)
        os.unlink(tmp2_g)
        os.unlink(tmp_trans1)
        os.unlink(tmp_trans)

@mkdir('geneset.dir')
@merge([PARAMS['prim_trans1'], PARAMS['prim_trans2']],
           r"geneset.dir/geneset_all.fa")
def buildReferenceKallisto(infiles, outfile):
    '''
    Builds a reference transcriptome and decoy sequneces for alevin and kallisto
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
    prim_trans1, prim_trans2 = infiles
    genome_file1 = PARAMS['genome1']

    if PARAMS['mixed_species']:
        genome_file2 = os.path.abspath(
        os.path.join(PARAMS["genome_dir2"], PARAMS["genome2"] + ".fa"))
        tmp1 = P.get_temp_filename('.')
        tmp2 = P.get_temp_filename('.')
        statement = '''
                       zcat %(geneset1)s |
                       awk '$3=="exon"'|
                       cgat gff2fasta
                       --is-gtf
                       --genome-file=%(genome_file1)s
                       --fold-at=60 -v 0
                       --log=%(outfile)s.log > %(tmp1)s &&
                       zcat %(geneset2)s |
                       awk '$3=="exon"'|
                       cgat gff2fasta
                       --is-gtf
                       --genome-file=%(genome_file2)s
                       --fold-at=60 -v 0
                       --log=%(outfile)s.log > %(tmp2)s &&
                       cat %(tmp1)s %(tmp2)s > %(outfile)s &&
                       samtools faidx %(outfile)s
                       '''
    else:
        statement = '''
                       grep "^>" <(gunzip -c %(genome_file1)s) | cut -d " " -f 1 > decoys.txt &&
                       sed -i.bak -e 's/>//g' decoys.txt &&
                       cat %(prim_trans1)s %(genome_file1)s > %(outfile)s
                       '''

    P.run(statement)

    if PARAMS['mixed_species']:
        os.unlink(tmp1)
        os.unlink(tmp2)


@active_if(PARAMS['salmon_alevin'])
@transform(buildReferenceSalmon,
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

    job_memory = "12G"

    statement = '''
    rm -rf %(outfile)s;
    salmon index -p 12 -d decoys.txt %(salmon_index_options)s -t %(infile)s -i %(outfile)s --gencode
    '''

    P.run(statement)

@active_if(PARAMS['kallisto_bustools'])
@transform(buildReferenceKallisto,
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

    geneset1 = PARAMS['geneset1']

    R_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = """Rscript %(R_ROOT)s/t2g.R -i %(geneset1)s -o %(outfile)s"""
    
    tmp = P.get_temp_filename('.')

    if PARAMS['mixed_species']:
        geneset2 = PARAMS['geneset2']

        statement = """cat  %(geneset1)s %(geneset2)s > %(tmp)s"""

        P.run(statement)

        statement = """Rscript %(R_ROOT)s/t2g.R -i %(tmp)s -o %(outfile)s"""

    P.run(statement, to_cluster=False)
    os.unlink(tmp)




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

    SEQUENCEFILES_SALMON_OUTPUT = (
        r"salmon.dir/%s/alevin/quants_mat.gz" % (
            PARAMS["merge_pattern_output"].strip()))

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.gz|fastq.1.gz)")

    SEQUENCEFILES_KALLISTO_OUTPUT = (
        r"kallisto.dir/\1/bus/output.bus")

    SEQUENCEFILES_SALMON_OUTPUT = (
        r"salmon.dir/\1/alevin/quants_mat.gz")

#############################
# Salmon- Alevin
#############################

@active_if(PARAMS['salmon_alevin'])
@follows(mkdir("salmon.dir"))
@transform(run_fastp,
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
    fastqfiles, index, t2gmap = infiles
    if isinstance(fastqfiles, list):
        CB_UMI_fastq = fastqfiles[0]
        reads_fastq = fastqfiles[1]

    outfolder = outfile.rsplit('/',2)[0].replace("salmon.dir/fastq_post.dir/", "salmon.dir/")

    salmon_options = PARAMS['salmon_run_options']

    statement = '''
    salmon alevin -l %(salmon_librarytype)s -1 %(CB_UMI_fastq)s -2  %(reads_fastq)s
    --%(salmon_sctechnology)s -i %(index)s -p %(salmon_threads)s -o %(outfolder)s
    --tgMap %(t2gmap)s --dumpFeatures --dumpUmiGraph %(salmon_options)s
    '''

    job_memory = PARAMS["salmon_job_memory"]
    job_threads = PARAMS['salmon_threads']
    P.run(statement)

#############################
# Kallisto- Bus
#############################

@active_if(PARAMS['kallisto_bustools'])
@follows(mkdir("kallisto.dir"))
@collate(run_fastp,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex, getTranscript2GeneMap),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def runKallistoBus(infiles, outfile):
    '''
    Generates BUS files for single-cell sequencing

    infiles: raw sequencing fastq files, kallisto index

    '''
    aligner = 'kallisto_bus'

    fastqfiles, index, t2gmap = infiles[0]
    if isinstance(fastqfiles, list):
        CM_UMI_fastq = fastqfiles[0]
        reads_fastq = fastqfiles[1]

    outfolder = outfile.rsplit('/',1)[0].replace("fastq_post.dir/", "")

    statement = '''
    kallisto bus -i %(index)s -o %(outfolder)s -x %(kallisto_sctechnology)s
    -t %(kallisto_threads)s %(CM_UMI_fastq)s %(reads_fastq)s
    '''

    job_memory = '20G'

    P.run(statement)

#########################
# Process bus file
#########################

# Must have bustools installed
#https://github.com/BUStools/bustools

@active_if(PARAMS['kallisto_bustools'])
@transform(runKallistoBus,
           regex("kallisto.dir/(\S+)/bus/output.bus"),
           r"kallisto.dir/\1/bus/output.bus.sorted.txt")
def busText(infile, outfile):
    '''
    Sort the bus file produced by kallisto and then convert it to a text file.
    '''

    tmp_bus  = P.get_temp_filename(".")
    job_memory = '10G'

    job_memory = '10G'

    statement = '''
    bustools sort -o %(tmp_bus)s %(infile)s ;
    bustools text -o %(outfile)s %(tmp_bus)s
    '''

    P.run(statement)



@transform(busText,
           suffix(".sorted.txt"),
           add_inputs(getTranscript2GeneMap),
           r"\1_GCcoordmatrix.mtx")
def busCount(infiles, outfile):
    '''
    Takes the sorted BUS file, corresponding ec matrix and transcript text file and generates a count matrix and tag count comparison??
    '''

    sorted_bus, t2gmap = infiles
    folder = sorted_bus.rsplit('/', 1)[0]
    ROOT = os.path.dirname(__file__)
    bus2count = ROOT + "python" + "/bus2count.py"
    exp_cells = PARAMS['kallisto_expectedcells']
    threads = PARAMS['kallisto_threads']

    statement = '''
    rm -rf %(folder)s/bus_count.log;
    python %(bus2count)s --dir %(folder)s --t2gmap %(t2gmap)s --expectedcells %(exp_cells)s --threads %(threads)s -o %(outfile)s
    '''

    job_memory = "30G"

    P.run(statement)

## Kallisto SCE object
@follows(mkdir("SCE.dir"))
@active_if(PARAMS['kallisto_bustools'])
@transform(busCount,
           regex("kallisto.dir/(\S+)/bus/output.bus_GCcoordmatrix.mtx"),
           r"SCE.dir/\1/bus/sce.rds")
def readBusSCE(infile, outfile):
    '''
    Takes in gene count matrices for each sample
    Creates a single cell experiment class in R and saves as an r object
    '''

    working_dir = os.getcwd()
    species = PARAMS['sce_species']
    gene_name = PARAMS['sce_genesymbol']
    pseudo = 'kallisto'

    job_memory = "10G"

    statement = '''
    Rscript %(R_ROOT)s/sce.R -w %(working_dir)s -i %(infile)s -o %(outfile)s --species %(species)s --genesymbol %(gene_name)s --pseudoaligner %(pseudo)s
    '''

    P.run(statement)

## Kallisto SCE object using BUSpaRse R package and emptydrops (DropletUtils function)
@follows(mkdir("SCE.dir"))
@follows(busText)
@active_if(PARAMS['kallisto_bustools'])
@transform(busText,
           regex("kallisto.dir/(\S+)/bus/output.bus.sorted.txt"),
           add_inputs(PARAMS['geneset'], getTranscript2GeneMap),
           r"SCE.dir/\1/bus/sce.rds")
def BUSpaRse(infiles, outfile):
    '''
    Create kallisto SCE object. Use BUSpaRse package to read in bus file and convert to TCC and gene counts matrix.
    Create knee plot and use point of inflection to estimate number of empty droplets and cells.
    Or use emptyDrops function from DropletUtils package to compare to the ambient profile.
    '''

    bus_text, gtf, t2gmap = infiles
    est_cells = PARAMS['kallisto_expectedcells']

    job_memory = '50G'

    statement = '''
    Rscript %(R_ROOT)s/BUSPaRse.R -i %(bus_text)s -o %(outfile)s --estcells %(est_cells)s --t2g %(t2gmap)s -g %(gtf)s
    '''

    P.run(statement)


#########################
# Multiqc
#########################

@follows(mkdir("MultiQC_report.dir"))
@follows(run_fastqc, runSalmonAlevin)
@originate("MultiQC_report.dir/multiqc_report.html")
def build_multiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)


@follows(mkdir("Report.dir"))
@follows(build_multiqc)
@originate("Report.dir/_site.yml")
def copy_report(infile):
    '''Copy the Rmarkdown report to current directory'''

    RMARKDOWN_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_singlecell","Rmarkdown")

    statement = '''cp %(RMARKDOWN_ROOT)s/* Report.dir/'''

    P.run(statement)


@follows(runSalmonAlevin, BUSpaRse, copy_report)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
