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

The default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains UMI/cellular barcode data and fastq.2.gz contains sequencing reads. 
Chromium outputis of the format: samplename_R1.fastq.gz and samplename_R2.fastq.gz so will require conversion to the default file format above.

Pipeline output
==================

The output of running this software is the generation of a SingleCellExperiment object and further downstream analysis including: clustering, pseudotime analysis, velocity time graphs, quality control analysis. 

Code
==================

"""
from ruffus import *

import sys
import os
import sqlite3

import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleSC

import pandas as pd
import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleSC

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

    job_memory = "12G"

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

############################################
# Perform read quality steps
############################################


@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           formatter(r"(?P<track>[^/]+).(?P<suffix>fastq.1.gz|fastq.gz)"),
           r"fastqc_pre.dir/{track[0]}.fastqc")
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

#########################
# Process bus file
#########################

# Must have bustools installed
# https://github.com/BUStools/bustools

@active_if(PARAMS['kallisto_bustools'])
@transform(runKallistoBus,
           suffix(".bus"),
           r"\1.bus.sorted.txt")
def busText(infile, outfile):
    '''
    Sort the bus file produced by kallisto and then convert it to a text file.
    '''

    tmp_bus  = P.get_temp_filename(".")

    statement = '''
    bustools sort -o %(tmp_bus)s %(infile)s ;
    bustools text -o %(outfile)s %(tmp_bus)s
    '''

    P.run(statement)



@transform(busText,
           suffix(".sorted.txt"),
           r"\1.count")
def busCount(infile, outfile):
    '''
    Takes the sorted BUS file, corresponding ec matrix and transcript text file and generates a count matrix and tag count comparison??
    ''' 

    folder = infile.rsplit('/', 1)[0]
    sc_directory = PARAMS['sc_dir']
    bus2count = sc_directory + "/pipelines/bus2count.py"

    statement = '''
    python %(bus2count)s --busdir %(folder)s 
    '''

    P.run(statement)

#########################
# SCE object  
#########################

@follows(mkdir("SCE.dir"))
@active_if(PARAMS['salmon_alevin'])
@transform(runSalmonAlevin,
           regex(r"salmon.dir/(.*)/alevin/quants_mat.gz"),
           r"SCE.dir/\1.rds")
def readAlevinSCE(infile,outfile):
    '''
    Collates alevin count matrices for each sample
    Creates a single cell experiment class in R and saves as an r object
    '''

    working_dir = os.getcwd()
    sc_directory = PARAMS['sc_dir']
    script_loc = sc_directory + "/pipelines/R/sce.R"
    
    job_memory = "10G"

    statement = '''
    Rscript %(script_loc)s -w %(working_dir)s -i %(infile)s -o %(outfile)s
    '''
    
    P.run(statement)

#########################
# Multiqc
#########################

@follows(mkdir("MultiQC_report.dir"))
@originate("Mapping_qc.dir/multiqc_report.html")
def build_multiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)

#########################
# QC step  
#########################

@follows(build_multiqc)
@follows(mkdir("QC_report.dir"))
@transform(readAlevinSCE,
           suffix(".rds"),
           "SCE.dir/pass.rds")
def run_qc(infile, outfile):
    """
    Runs an Rmarkdown report that allows users to visualise and set their
    quality parameters according to the data. The aim is for the pipeline
    to generate default thresholds then the user can open the Rmarkdown in
    rstudio and re-run the report, modifying parameters changesto suit the
    data
    """

    NOTEBOOK_ROOT = os.path.join(os.path.dirname(__file__), "Rmarkdown")

    #probably just need to knit one document not render_site
    statement = '''cp %(NOTEBOOK_ROOT)s/Sample_QC.Rmd QC_report.dir &&
                   cd QC_report.dir && R -e "rmarkdown::render('Sample_QC.Rmd',encoding = 'UTF-8')"'''

    P.run(statement)


#########################
# Seurat analysis  
#########################

# tSNE plotting on saved surat object - a range of perplexity choices
# plot tSNE perplexity hyper parameters on tSNE layout
# UMAP analysis
# Diffusion map
# find clusters (findMarkers i think the function is called)
# differential expression of markers in clusters


#########################
# Velocity analysis  
#########################

# Rmarkdown maybe then users can play around with parameters?


#########################
# Visulalisation of selected data  
#########################

# make violin plots from user selected genes
# heatmap of lists of genes

@follows(readAlevinSCE, busText)
def quant():
    pass

# what about adding a knee plot - how does alevin or kallisto handle the
# expected number of cells? - check documentation
@follows(run_qc)
def qc():
    pass


@follows(mkdir("Seurat.dir"))
@transform(readAlevinSCE,
           regex(r"SCE.dir/(.*).rds"),
           r"Seurat.dir/\1.rds")
def seurat_generate(infile,outfile):
    ''' 
    Takes sce object and converts it to a seurat object for further analysis
    '''
    working_dir = os.getcwd()
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    statement = '''
    Rscript %(R_ROOT)s/seurat.R -w %(working_dir)s -i %(infile)s -o %(outfile)s
    '''
    
    P.run(statement)


@transform(seurat_generate,
           regex("Seurat.dir/(\S+).rds"),
           r"Seurat.dir/\1.dim_reduction.rds")
def seurat_dimreduction(infile, outfile):
    '''
    Takes a seurate object and computes a PCA-based dimension reduction
    '''
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    statement = '''Rscript %(R_ROOT)s/seurat_dimreduction.R
    				-w %(working_dir)s
    				-i %(infile)s
    				-o %(outfile)s
    				--mingenes=%(seurat_mingenes)s
    				--maxmitopercent=%(seurat_maxmitopercent)s'''

    P.run(statement)


@transform(seurat_generate,
           regex("Seurat.dir/(\S+).rds"),
           r"Seurat.dir/\1.seurat_cluster.rds")
def seurat_clustering(infile, outfile):
    '''
    Takes sce seurat object and creates a series of cluster visualisations
    over a number of perplexities.
    '''
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")
	
    statement = '''Rscript %(R_ROOT)s/seurat_cluster.R'''

    P.run(statement)


@transform(seurat_clustering,
           regex("Seurat.dir/(\S+).rds"),
           r"Seurat.dir/\1.seurat_marker.rds")
def run_seurat_markdown(infile, outfile):
	'''
    Takes sce seurat object from clustering and generates
    an Rmarkdown report for running tsne, visualising clusters, finding marker
    genes and creating feature plots
    '''

	R_ROOT = os.path.join(os.path.dirname(__file__), "Rmarkdown")
	statement = ''' '''

	P.run(statement)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
