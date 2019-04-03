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
import re

import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleVelo

import cgat.GTF as GTF

#import cgatpipelines.tasks.mapping as mapping

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
		    "*.fastq.2.gz")
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


## Align reads (bam files needed)
# STAR

@active_if(PARAMS["star_run"])
@transform(buildReferenceTranscriptome,
           regex("(.*)/(.*).fa"),
           r"\1/star.index")
def buildStarIndex(infile, outfile):
    '''
    Builds a STAR index from the reference transcriptome.

    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    star_threads: str
        :term: `PARAMS` the number of threads for genome generation
    outfile: str
        path to output file
    '''

    os.mkdir("geneset.dir/star.index")

    threads = PARAMS["star_threads"]
    job_memory = PARAMS["star_memory"]
    RAM_limit = str(PARAMS["star_RAM_limit"])

    statement = '''
    STAR --runThreadN %(threads)s --runMode genomeGenerate --genomeDir %(outfile)s --genomeFastaFiles %(infile)s --limitGenomeGenerateRAM=%(RAM_limit)s
    '''

    P.run(statement)

# Merge lanes

if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]: 
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s.(fastq.gz|fastq.2.gz)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))

    SEQUENCEFILES_STAR_OUTPUT = (
        r"STAR.dir/%s_Aligned.sortedByCoord.out.bam" % (
            PARAMS["merge_pattern_output"].strip()))

    SEQUENCEFILES_WHITELIST_OUTPUT = (
        r"UMItools.dir/%s_whitelist.txt"  % (
            PARAMS["merge_pattern_output"].strip()))

    SEQUENCEFILES_ALEVIN_OUTPUT = (
        r"processed.dir/%s/%s_processed.fq" % (
            PARAMS["merge_pattern_output"].strip(),PARAMS["merge_pattern_output"].strip()))

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.gz|fastq.2.gz)")

    SEQUENCEFILES_STAR_OUTPUT = (
        r"STAR.dir/\1_Aligned.sortedByCoord.out.bam")

    SEQUENCEFILES_WHITELIST_OUTPUT = (
        r"UMItools.dir/\1_whitelist.txt")

    SEQUENCEFILES_ALEVIN_OUTPUT = (
        r"processed.dir/\1/\1_processed.fq")

## Use UMI tools to extract CB and UMIs

@collate(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           SEQUENCEFILES_WHITELIST_OUTPUT)
def generateCellularBarcodeWhitelist(infiles, outfile):
    '''
    Generate whitelist of CBs using UMI tools
    '''
    
    infiles = ModuleVelo.check_multiple_read_files(infiles)
    reads =  infiles
    temp_file = P.get_temp_filename(".")
    BC_pattern = PARAMS['UMItools_BCpattern']
    cell_num = PARAMS['UMItools_cell_number']

    if isinstance(reads, list):
        reads = " ".join(reads)

    statement = '''
    cat reads > %(temp_file)s ; 
    umi_tools whitelist --stdin %(temp_file)s  --bc-pattern=%(BC_pattern)s --set-cell-number=%(cell_num)s --log2stderr > %(outfile)s
    '''

## Use alevin to error correct CBs and attach to read file

@active_if(PARAMS['salmon_alevin'])
@follows(mkdir("processed.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         SEQUENCEFILES_ALEVIN_OUTPUT)
def alevinDumpFastq(infiles, outfile):
    '''
    Salmon alevin used with options dumpfq and noQuant to  
    error correct CBs and attach to header file. No pseudoalignment
    performed and no counts matrix generated
    '''

    fastqfile = ModuleVelo.check_multiple_read_files_no_index(infiles)
    fastqfiles = ModuleVelo.check_paired_end(fastqfile)
    tmp_fastq = P.get_temp_filename(".")

    if isinstance(fastqfiles, list):
        CB_UMI_fastq = " ".join(fastqfiles[0])
        reads_fastq = " ".join(fastqfiles[1]) 

    if PARAMS['salmon_default_barcodes']:
        barcode_options = ''
    else:
        barcode_options = '--barcodeLength %s --umiLength %s --end %s'%(PARAMS['salmon_barcode_length'], PARAMS['salmon_umi_length'], PARAMS['salmon_end'])

    outfolder = re.search('(processed.dir/.*)/.*\.fq', outfile)[1]

    # Need to give files for index and transcript2genemap despite not needing them as no quant, have yml files as placeholders
    statement = '''
    salmon alevin -l %(salmon_librarytype)s -1 %(CB_UMI_fastq)s -2  %(reads_fastq)s
    --%(salmon_sctechnology)s --noQuant --dumpfq -i pipeline.yml  --tgMap pipeline.yml -p %(salmon_threads)s %(barcode_options)s -o %(outfolder)s > %(tmp_fastq)s &&
    sed '/@/,$!d' tmp_fastq > %(outfile)s
    '''

    P.run(statement)

#################
# Mapping
#################

@active_if(PARAMS["star_run"])
@follows(mkdir("STAR.dir"))
@collate(SEQUENCEFILES,
           SEQUENCEFILES_REGEX, 
           add_inputs(buildStarIndex),
           SEQUENCEFILES_STAR_OUTPUT)
def mapReadsWithSTAR(infiles, outfile):
    '''
    Mapping reads to the genome/transcriptome.
    '''
    
    infiles = ModuleVelo.check_multiple_read_files(infiles)
    reads, star_index = infiles

    if isinstance(reads, list):
        reads = " ".join(reads)

    sample = re.search('STAR\.dir/(.*)Aligned\.out\.bam', outfile)[1]

    threads = PARAMS["star_threads"]
    job_memory = PARAMS["star_memory"]

    statement = '''
    STAR --runThreadN %(threads)s --genomeDir %(star_index)s --readFilesIn %(reads)s --outFileNamePrefix STAR.dir/%(sample)s --readFilesCommand zcat 
    '''
    # --outSAMtype BAM Unsorted SortedByCoordinate
    P.run(statement)


@active_if(PARAMS["star_run"])
@follows(mkdir("STAR.dir"))
@transform(alevinDumpFastq,
           regex("processed.dir/(.*)/.*_processed.fq"), 
           add_inputs(buildStarIndex),
           r"STAR.dir/\1_Aligned.out.bam")
def mapReadsWithSTAR2(infiles, outfile):
    '''
    Mapping reads to the genome/transcriptome.
    '''
    
    reads, star_index = infiles

    sample = re.search('STAR\.dir/(.*)Aligned\.out\.bam', outfile)[1]
    temp_file = P.get_temp_filename(".")

    threads = PARAMS["star_threads"]
    job_memory = PARAMS["star_memory"]

    statement = '''
    STAR --runThreadN %(threads)s --genomeDir %(star_index)s --readFilesIn %(reads)s --outFileNamePrefix STAR.dir/%(sample)s_star 
    --outSAMtype BAM Unsorted SortedByCoordinate && samtools sort -T %(temp_file)s -o %(outfile)s STAR.dir/%(sample)s_starAligned.out.bam  && samtools index %(outfile)s
    '''

    # --readFilesCommand zcat
    P.run(statement)

@transform(mapReadsWithSTAR2,
           regex("STAR.dir/(.*)\.out\.bam"),
           r"STAR.dir/\1_tag.out.bam")
def tagBAM(infile,outfile):
    ''' 
    Extract CBs and UMIs from read names and assign tags CB and UB in BAM file
    Using pysam module
    '''

    sc_directory = PARAMS['sc_dir']
    tagbam = sc_directory + "/pipelines/TagBAM.py"

    statement = '''
    python %(tagbam)s -i %(infile)s -o %(outfile)s
    '''

    P.run(statement)


## velocyto run


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
