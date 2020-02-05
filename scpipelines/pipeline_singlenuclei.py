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
import scpipelines.ModuleSC

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

# Determine read length - need documentation  for how to download the data according to lior
'''
$ zcat R1.fastq.gz | head -2 ## note on a mac you would do zcat < R1.fastq.gz | head
@SRR8742283.1 NS500422:552:HJ5Y3BGX3:1:11101:21875:1038 length=61
CAGTCNTTTTTTTTAATTTAAAAAAAAAAAAAAGATTTATTAACAGTTTTAGAAGGCAGTT

$ echo -n CAGTCNTTTTTTTTAATTTAAAAAAAAAAAAAAGATTTATTAACAGTTTTAGAAGGCAGTT | wc -c
61
L= 61
'''

@mkdir('geneset.dir')
@originate("geneset.dir/tr2g.txt")
def t2g(outfile):
    '''This function will convert transcript to genes using t2g
       script'''

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__)))

    statement = '''
       zcat < %(geneset)s | %(PY_SRC_PATH)s/t2g.py  > %(outfile)s'''

    P.run(statement)


@mkdir('geneset.dir')
@originate("geneset.dir/introns.fa")
def intron_bed2fa(outfile):
    '''This converts introns bed to introns fa'''

    tmp_bed = P.get_temp_filename(".")

    statement = '''zcat < %(intron_bed)s > %(tmp_bed)s &&
    # index the fa file before next line &&
    bedtools getfasta -name -fo introns.fa -fi %(genome_file)s -bed %(tmp_bed)s'''

    P.run(statement)
    os.unlink(tmp_bed)



@transform(intron_bed2fa,
           regex("(\S+).fa"),
           r"geneset.dir/\1_transcripts.txt")
def introns_transcripts(infile, outfile):
    '''get a list of all of the intronic transcript IDs represented in our FASTA file, with version numbers.'''

    statement = '''cat %(infile)s | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3"."$4}' > %(outfile)s'''

    P.run(statement)


@transform(introns_transcripts,
           regex("geneset.dir/(\S+)_transcripts.txt"),
           r"geneset.dir/\1_transcripts_no_version.txt")
def introns_transcripts_no_version(infile, outfile):
    '''list of all of the intronic transcript IDs represented in our FASTA file, without version numbers.'''

    staement = '''cat %(infile)s | tr "." " " | awk '{print $1}' > %(outfile)s '''

    P.run(statement)


@transform(introns_transcripts,
           regex("geneset.dir/(\S+)_transcripts.txt"),
           r"geneset.dir/\1_transcripts.to_capture.txt")
def add_identifier(infile outfile):
    '''add an identifier to the transcript IDs '''

    statement = '''cat %(infile)s | awk '{print $0"."NR"-I"}' > %(outfile)s'''

    P.run(statement)


@merge([introns_transcripts,t2g],
       r"geneset.dir/introns_t2g.txt")
def map_trans2gene(infiles, outfile):
    '''map the transcripts to their respective genes.'''

    introns_transcripts, tr2g  = infiles

    statement = '''awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' %(tr2g)s %(introns_transcripts)s > %(outfile)s'''

    P.run(statement)


@merge([map_trans2gene,intron_bed2fa],
           r"geneset.dir/\introns.correct_header.fa")
def fix_intron_fasta(infiles, outfile):
    '''fix all of the headers for the introns FASTA file so that they
    contain the transcript ID, an identifier specifying that the transcript
    is an “intronic” transcript, and a unique number to avoid duplicates.'''

    introns_t2g, introns = infiles

    tmp_fasta = P.get_temp_filename(".")

    statement = '''awk '{print ">"$1"."NR"-I"" gene_id:"$2" gene_name:"$3}' %(introns_t2g)s > %(tmp_fasta)s  &&
                awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' %(tmp_fasta)s  %(introns)s > %(outfile)s '''

    P.run(statement)
    os.unlink()


@mkdir('geneset.dir')
@originate("geneset.dir/cDNA_transcripts_no_version.txt")
def capture_list(outfile):
    '''Get the transcripts to capture list and transcripts to genes for cDNA'''

    tmp_cdna = P.get_temp_filename(".")

    statement = '''zcat < %(cdna_fasta)s cDNA.fa | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3}' > %(tmp_cdna)s  &&
                   cat %(tmp_cdna)s  | tr "." " " | awk '{print $1}' > %(outfile)s '''


    P.run(statement)
    os.unlink()


@mkdir('geneset.dir')
@originate("geneset.dir/cDNA_transcripts.to_capture.txt")
def map_tran_gene(outfile):
    ''''Add an identifier to the transcript IDs'''

    tmp_cdna = P.get_temp_filename(".")

    statement = '''zcat < %(cdna_fasta)s | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3}' > %(tmp_cdna)s &&
                   cat %(tmp_cdna)s | awk '{print $0"."NR}' > %(outfile)s'''

    P.run(statement)


@transform(t2g,
           regex("(\S+).txt"),
           r"geneset.dir/cDNA_\1.txt")
def map_tr2gene(infile, outfile):
    '''Map the transcripts to genes.'''

    tmp_cdna = P.get_temp_filename(".")

    statement = '''zcat < %(cdna_fasta)s | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3}' > %(tmp_cdna)s &&
                   awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' %(infile)s %(tmp_cdna)s  > %(outfile)s'''


    P.run(statement)


@transform(map_trans2gene,
           regex("(\S+)_t2g.txt"),
           r"geneset.dir/\1.correct_fasta.fa")
def find_intron_fa_header(infile, outfile):
    '''Fix the INTRONS FASTA header'''

    statement = '''awk '{print ">"$1"."NR" gene_id:"$2" gene_name:"$3}' %(infile)s > geneset.dir/cDNA_fasta_header.txt &&
                   awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' geneset.dir/cDNA_fasta_header.txt $cDNA_fa >
                  %(outfile)s'''

    P.run(statement)


@mkdir('kallisto.dir')
@merge([find_intron_fa_header,fix_intron_fasta, map_tr2gene, map_trans2gene],
           "kallisto.dir/kallisto.index")
def buildKallistoIndex(infiles, outfile):
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
    cDNA_correct_header, introns_correct_header, map_tr2gene, map_trans2gene = infiles

    job_memory = "12G"

    statement = '''
    cat %(cDNA_correct_header)s %(introns_correct_header)s > kallisto.dir/cDNA_introns.fa
    cat %(map_tr2gene)s %(map_trans2gene)s > kallisto.dir/cDNA_introns_t2g.txt
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s kallisto.dir/cDNA_introns.fa
    '''

    P.run(statement)


############################################
# Perform read quality steps
############################################


@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           regex("(\S+).fastq.(\d).gz"),
           r"fastqc_pre.dir/\1.fastq.\2_fastqc.html")
def runFastQC(infile, outfile):
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
         add_inputs(buildKallistoIndex),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def runKallistoBus(infiles, outfile):
    '''
    Generates BUS files for single-cell sequencing

    infiles: raw sequencing fastq files, kallisto index

    '''
    aligner = 'kallisto_bus'
    infiles = ModuleSC.check_multiple_read_files(infiles)
    fastqfile, index = infiles
    fastqfiles = ModuleSC.check_paired_end(fastqfile, aligner)
    fastqfiles = " ".join(fastqfiles)

    outfolder = outfile.rsplit('/',1)[0]

    statement = '''
    kallisto bus -i %(index)s -o %(outfolder)s -x %(kallisto_sctechnology)s
    -t %(kallisto_threads)s %(fastqfiles)s
    '''

    job_memory = '20G'

    P.run(statement)


@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         r"\1_whitelist.txt")
def whitelist(infile, outfile):
    '''use umitools to generate whitelist of barcodes'''

    statement = '''
                umi_tools whitelist --stdin=%(infile)s
                --bc-pattern=%(umitools_barcode_pattern)s
                --extract-method=regex
                --log2stderr
                > %(outfile)s
    '''

    P.run(statement)


@transform(runKallistoBus,
           regex("kallisto.dir/(\S+)/bus/output.bus"),
           add_inputs(whitelist),
           r"kallisto.dir/\1/bus/output.correct.bus")
def bustools_correct(infiles, outfile):
    ''' then
    bustools correct function'''

    # may need a function that collates the whitelists
    # znd selects the correct whitelist for the correct fastq
    bus_file, whitelist = infiles

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


@transform(bustools_correct,
           regex("kallisto.dir/(\S+)/bus/output.sort.bus"),
           add_inputs(introns_transcripts_no_version),
           r"kallisto.dir/\1/bus/introns_capture.bus")
def bustools_capture_intron(infiles, outfile):
    '''use bustools capture for cDNA '''

    infile, capture_list = infiles

    matrix = infile.replace("introns_capture.bus","matrix.ec")
    trans = infile.replace("introns_capture.bus","transcripts.txt")

    statement = '''
    bustools capture -s -o %(outfile)s -c %(capture_list)s  -e %(matrix)s -t %(trans)s  %(infile)s
    '''

    P.run(statement)

# Bustools capture cDNA and then introns
@transform(bustools_correct,
           regex("kallisto.dir/(\S+)/bus/output.sort.bus"),
           add_inputs(capture_list),
           r"kallisto.dir/\1/bus/cDNA_capture.bus")
def bustools_capture_cdna(infiles, outfile):
    '''use bustools capture for cDNA '''

    infile, capture_list = infiles

    matrix = infile.replace("cDNA_capture.bus","matrix.ec")
    trans = infile.replace("cDNA_capture.bus","transcripts.txt")

    statement = '''
    bustools capture -s -o %(outfile)s -c %(capture_list)s  -e %(matrix)s -t %(trans)s  %(infile)s
    '''

    P.run(statement)

@transform(bustools_capture_intron,
           regex("kallisto.dir/(\S+)/bus/introns_capture.bus"),
           r"kallisto.dir/\1/bus/unspliced/unspliced")
def bustools_count_intron(infile, outfile):
# Merge spliced and unspliced

    matrix = infile.replace("introns_capture.bus","matrix.ec")
    trans = infile.replace("introns_capture.bus","transcripts.txt")

    statement = '''
    bustools count -o %(outfile)s -g kallisto.dir/cDNA_introns_t2g.txt -e %(matrix)s -t %(trans)s --genecounts %(infile)s
    '''

    P.run(statement)


@transform(bustools_capture_cdna,
           regex("kallisto.dir/(\S+)/bus/cDNA_capture.bus"),
           r"kallisto.dir/\1/bus/spliced/spliced")
def bustools_count_cdna(infile, outfile):
# Merge spliced and unspliced

    matrix = infile.replace("cDNA_capture.bus","matrix.ec")
    trans = infile.replace("cDNA_capture.bus","transcripts.txt")

    statement = '''
    bustools count -o %(outfile)s -g kallisto.dir/cDNA_introns_t2g.txt -e %(matrix)s -t %(trans)s --genecounts %(infile)s
    '''

    P.run(statement)


#########################
# SCE object
#########################

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
    R_ROOT = os.path.join(os.path.dirname(__file__),"pipeline_singlecell","R")
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
@follows(runFastQC, bustools_count_cdna)
@originate("MultiQC_report.dir/multiqc_report.html")
def build_multiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)


@follows(BUSpaRse, build_multiqc)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
