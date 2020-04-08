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

ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kallistobus","python")
R_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_kallistobus", "R")


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

# Determine the location of the input fastq files
SEQUENCESUFFIXES = ("*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

############################################
# Build indexes
############################################

genome_fasta = PARAMS['genome_dir'] + PARAMS['genome'] + ".fa"

@active_if(PARAMS['star_index'] == 1)
@mkdir('star_index.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           add_inputs(genome_fasta),
           r"star_index.dir/Genome")
def index_genome_star(infiles, outfile):
    """
    Generate a star index if no index is supplied
    """

    gtffile, genome = infiles
    threads = PARAMS['star_threads']

    statement = """STAR
                   --runThreadN %(threads)s
                   --runMode genomeGenerate
                   --genomeDir star_index.dir/
                   --genomeFastaFiles %(genome)s"""

    job_memory = 'unlimited'

    P.run(statement)

@active_if(PARAMS['star_index'] == 0)
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

    job_memory = '10G'

    statement = """picard FastqToSam
                   F1=%(infile)s
                   F2=%(second_read)s
                   O=%(outfile)s
                   SM=%(name)s"""

    job_memory = '20G'
    P.run(statement)


@transform(fastq_to_bam,
           regex("(\S+)_unmapped.bam"),
           r"\1_tagged_Cell.bam")
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


    job_memory = '30G'

    P.run(statement)

@transform(cell_barcode_bam,
           regex("(\S+)_tagged_Cell.bam"),
           r"\1_tagged_CellMolecular.bam")
def molecular_barcode_bam(infile, outfile):
    """
    Extracts bases from UMI barcode encoding read (fastq.1.gz)
    Creates a new BAM tag with those bases on the genome read.
    """

    name = infile.replace("_tagged_Cell.bam", "")

    statement = """TagBamWithReadSequenceExtended
                   INPUT=%(infile)s
                   OUTPUT=%(outfile)s
                   SUMMARY=%(name)s_tagged_Molecular.bam_summary.txt
                   BASE_RANGE=13-20
                   BASE_QUALITY=10
                   BARCODED_READ=1
                   DISCARD_READ=True
                   TAG_NAME=XM
                   NUM_BASES_BELOW_QUALITY=1 """

    job_memory = '30G'

    P.run(statement)


@transform(molecular_barcode_bam,
           regex("(\S+)_tagged_CellMolecular.bam"),
           r"\1_tagged_filtered.bam")
def filter_bam(infile, outfile):
    """
    filter the bam file to reject XQ tag
    """

    statement = """FilterBam
                   TAG_REJECT=XQ
                   INPUT=%(infile)s
                   OUTPUT=%(outfile)s"""

    job_memory = '10G'

    P.run(statement)

@transform(filter_bam,
           regex("(\S+)_tagged_filtered.bam"),
           r"\1_tagged_trimmed_smart.bam")
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

    job_memory  ='10G'

    P.run(statement)


@transform(trim_starting_sequence,
           regex("(\S+)_tagged_trimmed_smart.bam"),
           r"\1_polyA_filtered.bam")
def polyA_trimmer(infile, outfile):
    """
    Remove the poly A tail from each read
    """

    name = infile.replace("_tagged_trimmed_smart.bam", "")


    statement = """PolyATrimmer
                   INPUT=%(infile)s
                   OUTPUT=%(outfile)s
                   OUTPUT_SUMMARY=%(name)s_polyA_trimming_report.txt
                   MISMATCHES=0
                   NUM_BASES=6
                   USE_NEW_TRIMMER=true"""

    job_memory = '20G'

    P.run(statement)


@follows(mkdir("fastq_file.dir"))
@transform(polyA_trimmer,
           regex("data.dir/(\S+)_polyA_filtered.bam"),
           r"fastq_file.dir/\1.fastq")
def bam_to_fastq(infile, outfile):
    """
    Convert to fastq file so it can be mapped using star
    """

    statement = """picard SamToFastq
                   INPUT=%(infile)s
                   FASTQ=%(outfile)s"""

    job_memory = '20G'

    P.run(statement)


@follows(mkdir("star.dir"))
@transform(bam_to_fastq,
           regex("fastq_file.dir/(\S+).fastq"),
           add_inputs(PARAMS['geneset'],
                      genome_fasta),
           r"star.dir/\1_Aligned.out.bam")
def star_mapping(infiles, outfile):
    """
    Perform star mapping. Ouputs unsorted BAM file.
    """

    fastq, gtffile, genome = infiles
    dirs, files = os.path.split(outfile)
    name = fastq.replace("fastq_file.dir/","")
    outfile_name = name.replace(".fastq","")
    threads = PARAMS['star_threads']

    statement = """STAR
                   --readFilesIn %(fastq)s
                   --runThreadN %(threads)s
                   --genomeDir star_index.dir
                   --outSAMmultNmax 1
                   --outSAMunmapped Within
                   --outSAMtype BAM Unsorted
                   --outFileNamePrefix star.dir/%(outfile_name)s_"""

    job_memory = 'unlimited'
    P.run(statement)


@transform(star_mapping,
           regex("star.dir/(\S+)_Aligned.out.bam"),
           r"star.dir/\1_Aligned.sorted_qn.out.bam")
def sort_query_name(infile,outfile):
    """
    Sort bam file using picard program SortSam
    """

    statement = """ picard SortSam
                    INPUT=%(infile)s
                    OUTPUT=%(outfile)s
                    SORT_ORDER=queryname """

    job_memory = '20G'
    P.run(statement)


@transform(sort_query_name,
           regex("star.dir/(\S+)_Aligned.sorted_qn.out.bam"),
           add_inputs(genome_fasta),
           r"star.dir/\1_Aligned.tagged.sorted_qn.out.bam")
def recover_tags(infiles, outfile):
    """
    Script using pysam to merge the sorted aligned output from STAR with the unaligned BAM that has been tagged with cellular and molecular barcodes (XC/XM).
    """

    mapped_bam, genome_fasta = infiles
    unmapped_bam = "data.dir/" + os.path.split(mapped_bam)[1].replace("Aligned.sorted_qn.out.bam", "polyA_filtered.bam")
    mergebam = ROOT + "/MergeBam.py"

    statement = ''' python %(mergebam)s --unmapped %(unmapped_bam)s --aligned %(mapped_bam)s -o %(outfile)s '''

    job_memory ='10G'

    P.run(statement)

## Merge lanes

@active_if(PARAMS["merge_pattern_input"])
@collate(recover_tags,
     regex("%s_(\S+)\.bam" % PARAMS["merge_pattern_input"].strip()),
     # the last expression counts number of groups in pattern_input
     r"%s.merged.aligned.coord.bam" % (PARAMS["merge_pattern_output"]))
def mergeBAMFiles(infiles, outfile):
    '''merge BAM files from the same experiment using user-defined regex
    For the mapping stages it is beneficial to perform mapping
    seperately for each sequence read infile(s) per sample so that
    the consistency can be checked. However, for downstream tasks,
    the merged :term:`bam` alignment files are required.
    Parameters
    ----------
    infiles : list
       list of :term:`bam` format alignment files
    outfile : str
       Output filename in :term:`bam` format
    '''

    if "merge_pattern_output" not in PARAMS or \
       not PARAMS["merge_pattern_output"]:
        raise ValueError(
            "no output pattern 'merge_pattern_output' specified")

    if len(infiles) == 1:
        if not os.path.isfile(os.path.join(infiles[0], outfile)):
            E.info(
                "%(outfile)s: only one file for merging - creating "
                "softlink" % locals())
            os.symlink(os.path.basename(infiles[0]), outfile)
            os.symlink(os.path.basename(infiles[0]) + ".bai", outfile + ".bai")
            return
        else:
            E.info(
                "%(outfile)s: only one file for merging - softlink "
                "already exists" % locals())
            return

    infiles = " ".join(infiles)
    tmp_bam  = P.get_temp_filename(".")

    statement = '''
    samtools merge %(tmp_bam)s %(infiles)s >& %(outfile)s_merge.log &&
    samtools sort %(tmp_bam)s -o %(outfile)s &&
    samtools index %(outfile)s
    '''

    job_memory = '20G'

    P.run(statement)


@active_if(not PARAMS["merge_pattern_input"] or "merge_pattern_input" not in PARAMS)
@transform(recover_tags,
           suffix("_Aligned.tagged.sorted_qn.out.bam"),
           ".merged.aligned.coord.bam")
def cpBAMFile(infile, outfile):
    ''' Non-merged files have the same file name for next task'''

    statement = '''cp %(infile)s %(outfile)s'''

    P.run(statement)

@follows(cpBAMFile)
@transform(PARAMS['geneset'],
           suffix(".gtf.gz"),
           ".gtf")
def extract_geneset(infile, outfile):
    ''' Gunzip geneset gtf file '''

    statement = '''gunzip -c %(infile)s > %(outfile)s'''

    P.run(statement)

@follows(mkdir('dropest.dir'))
@follows(cpBAMFile)
@follows(mergeBAMFiles)
@transform(mergeBAMFiles,
           regex("star.dir/(\S+).bam"),
           add_inputs(extract_geneset),
           r"dropest.dir/\1.filtered.bam")
def run_dropest(infiles, outfile):
    """
    Make new config file with location of whitelisted barcodes.
    Runs dropEst on the bam file and generate an rds file as output.
    """

    bamfile, gtffile = infiles
    generic_config_path = PARAMS['velocyto_dropest_config']
    barcode_suffix = PARAMS['velocyto_whitelist_suffix']

    sample_name = os.path.basename(bamfile).replace(".merged.aligned.coord.bam", "")
    new_config = "dropest.dir/" + sample_name + "_config_desc.xml"

    dropEst_out = "dropest.dir/" + sample_name + "_dropEst"
    # dropEst dumps filtered bam file in current directory, move to dropest.dir
    bam_out = sample_name + ".merged.aligned.coord.filtered.bam"

    config_file = ROOT + "/BarcodeFileDropest.py"

    # Not sure if -f (lowercase) is giving the problems. When removed struggled with UMIs. Try without

    statement = """python %(config_file)s --input %(bamfile)s --barcode_suffix %(barcode_suffix)s --config %(generic_config_path)s &&
                   dropest -m -f -F -L eiEIBA -g %(gtffile)s -o %(dropEst_out)s -c %(new_config)s %(bamfile)s &&
                   mv -T %(bam_out)s %(outfile)s"""

    job_memory = '50G'

    P.run(statement)


# Is this needed??? dropEst seems to barcode correct- jupyter notebook. barcodes change between CR and XC
@follows(mkdir('velocyto.dir'))
@transform(run_dropest,
           regex("dropest.dir/(\S+)\.merged\.aligned\.coord\.filtered\.bam"),
           r"velocyto.dir/\1/correct_\1.bam")
def dropest_bc_correct(infile, outfile):
    """
    Velocyto tools to make a new error corrected bam file
    """

    sample_name = infile.replace(".merged.aligned.coord.filtered.bam", "")
    #bam_path = "./velocyto.dir/" + sample_name + "/" + infile
    rds_path = sample_name + "_dropEst.rds"
    # not sure if right rds file, could be _dropEst.matrices.rds

    statement = """velocyto tools dropest-bc-correct -o  %(outfile)s %(infile)s %(rds_path)s"""

    job_memory = '20G'

    P.run(statement)

@follows(mkdir("velocyto.dir"))
@transform(dropest_bc_correct,
           regex("velocyto.dir/(\S+)/correct_(\S+).bam"),
           r"velocyto.dir/\1/cellsorted_correct_\1.bam")
def CB_sort(infile, outfile):
    '''
    First step of velocyto is to sort using samtools by CB.
    Doing this step first helps with parellelisation when running bulk jobs.
    Avoids runtimes errors using velocyto run.
    '''

    statement = """samtools sort -t CB -O BAM -o %(outfile)s %(infile)s"""

    P.run(statement)


@follows(CB_sort)
@transform(dropest_bc_correct,
           regex("velocyto.dir/(\S+)/correct_(\S+).bam"),
           add_inputs(extract_geneset),
           r"velocyto.dir/\1/\2.loom")
def velocyto_run_dropest(infiles, outfile):
    """
    Generate a loom file using BAM file preprocessed with dropEst tools.
    """

    bamfile, gtf_file = infiles
    output_folder = os.path.split(bamfile)[0]

    sample_name = os.path.basename(bamfile).replace(".bam","").replace("correct_", "")

    if PARAMS['velocyto_whitelist_active']:
        barcode_file = sample_name + PARAMS['velocyto_whitelist_suffix']
        whitelist = "-b " + barcode_file
    else:
        whitelist = ""

    if PARAMS['velocyto_mask_active']:
        mask = "-m " + PARAMS['velocyto_mask_file']
    else:
        mask = ""

    statement = """ velocyto run-dropest %(whitelist)s -o %(output_folder)s  -e %(sample_name)s %(mask)s %(bamfile)s %(gtf_file)s >
    %(sample_name)s_loom.log"""

    job_memory = '50G'

    P.run(statement)

@transform(velocyto_run_dropest,
         suffix(".loom"),
         "_cell_velocity.png")
def loom_analysis(infile, outfile):
    '''
    Takes all seurat.rds objects and combines them using RunCCA/MergeSeurat/RunMultiCCA
    into one large seurat object with annotations for each sample, ready for the monocle library.
    '''

    statement = """ Rscript %(R_ROOT)s/velocyto_analysis.R -l %(infile)s """

    job_memory = '20G'

    P.run(statement)

@follows()
def velocyto(run_dropest):
    pass

@follows()
def clear_temps():
    """
    Clears temporary files
    """

    statement = """ rm ctmp* """

    P.run(statement)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
