"""
===========
Pipeline kb
===========


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



# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


JUPYTER_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_quantcells","Jupyter")

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


    job_memory = "10G"

    statement = '''
    kallisto index -i %(outfile)s  %(geneset)s 2> index.log
    '''

    P.run(statement)


############################################
# Perform read quality steps
############################################


@follows(build_kallisto_index)
@follows(mkdir("fastqc_pre.dir"))
@transform(SEQUENCEFILES,
           regex("{}/(\S+).fastq.(\d).gz".format(DATADIR)),
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
SEQUENCEFILES_REGEX = regex(
        r"%s/(\S+).(fastq.gz|fastq.1.gz)" % (
            DATADIR))

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


    statement = '''
     kallisto bus -i %(index_files)s -t %(kallisto_threads)s -x %(kallisto_sctechnology)s
    -o %(outfile)s   %(fastqfiles)s
    2> %(outfile)s_kblog.log
    '''

    job_memory = PARAMS['bustools_memory']

    P.run(statement, job_options='-t 24:00:00')


@jobs_limit(1)
@transform(run_kallisto_bus,
           regex("(\S+)/output.bus"),
           r"\1/tr2gene.tsv")
def build_tr2g(infile, outfile):
    """Build a transcript to gene relationship """

    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))
    PYTHON_PATH =  os.path.abspath(os.path.join(os.path.dirname(__file__),"python"))

    if PARAMS['mixed_species']:


        if PARAMS['custom_fasta']:

            statement = """zcat %(input1)s %(input2)s > cat_input.fasta && 
                           python %(PYTHON_PATH)s/tr2gene.py --fasta cat_input.fasta --output %(outfile)s 2> %(outfile)s.log &&
                           rm -rf cat_input.fasta"""

        else:
            input1, input2 = PARAMS['geneset'].split(" ")

            statement = """zcat %(input1)s %(input2)s > cat_input.fasta && 
                       Rscript %(R_PATH)s/make_tr2gene.R -i cat_input.fasta  -o %(infile)s/
                       -f %(outfile)s 2> %(outfile)s.log && rm -rf tx_filtered.fa tr2g.tsv """
    else:
        if PARAMS['custom_fasta']:

            statement = "python %(PYTHON_PATH)s/tr2gene.py --fasta %(geneset)s --output %(outfile)s 2> %(outfile)s.log"

        else:

            statement = """Rscript  %(R_PATH)s/make_tr2gene.R -i %(geneset)s -o %(infile)s -f %(outfile)s"""

    P.run(statement)


@transform(run_kallisto_bus,
           suffix(".bus"),
           "_sorted.bus")
def bustools_sort(infile, outfile):
    """
    Generate a sorted bus file
    """

    tmp = P.get_temp_filename(".")

    statement = """bustools sort -T %(tmp)s -t %(kallisto_threads)s -o %(outfile)s %(infile)s/output.bus"""

    P.run(statement)

@follows(build_tr2g)
@transform(bustools_sort,
           regex("(\S+)/output_sorted.bus"),
           r"\1/genecount/genes.barcodes.txt")
def bustools_count(infile, outfile):
    """Generate a counts file from bus record"""

    out_dir = outfile.replace("genes.barcodes.txt", "genes")
    tr2g = infile.replace("output_sorted.bus","tr2gene.tsv")
    mat = infile.replace("output_sorted.bus","output.bus/matrix.ec")
    trans = infile.replace("output_sorted.bus","output.bus/transcripts.txt")



    statement = """bustools count -o %(out_dir)s -g %(tr2g)s -e %(mat)s -t %(trans)s --genecounts 2> %(out_dir)s.count.log %(infile)s"""
    P.run(statement, job_options='-t 24:00:00')


@transform(bustools_count,
           regex(r"(\S+)/genecount/genes.barcodes.txt"),
           r"\1/adata.h5ad")
def generate_h5ad(infile, outfile):
    """Convert mtx files to h5ad format"""
    
    PYTHON_PATH =  os.path.abspath(os.path.join(os.path.dirname(__file__),"python"))

    genecount_dir = os.path.dirname(infile)
    sample_dir = os.path.dirname(genecount_dir)
    
    matrix_file = os.path.join(genecount_dir, "genes.mtx")
    genes_file = os.path.join(genecount_dir, "genes.genes.txt")
    barcodes_file = os.path.join(genecount_dir, "genes.barcodes.txt")
    
    output_file = outfile

    statement = '''
    python %(PYTHON_PATH)s/parse_mtx_to_h5ad.py --matrix {matrix_file} \
                                --genes {genes_file} \
                                --barcodes {barcodes_file} \
                                --output {output_file}
    '''.format(matrix_file=matrix_file,
               genes_file=genes_file,
               barcodes_file=barcodes_file,
               output_file=output_file)

    P.run(statement)


@active_if(PARAMS['mixed_species'])
@transform(bustools_count,
           regex("(\S+)/genes.barcodes.txt"),
           r"\1")
def barnyard_plot(infile, outfile):
    """Construct a barnyard plot from bus file """

    infile = infile.replace("/genes.barcodes.txt", "")

    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

    statement = """Rscript %(R_PATH)s/plot_barnyard.R -i %(infile)s -o %(outfile)s"""

    P.run(statement)


@follows(generate_h5ad, barnyard_plot)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
