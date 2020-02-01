
===================
singlecell pipeline
===================

Overview
==================

This pipeline performs alignment free based quantification of drop-seq, 10X
single-cell sequencing analysis using wither kallisto or salmon.
Pseudoalignment is performed on the RNA reads,
using kallisto or Alevin and the resulting data is quantitatively
and qualitatively analysed.

The pipeline performs the following analyses:
* Alignment using kallisto or alevin (part of salmon)
* QC of reads using the scater package

Input files
-----------

The pipeline is ran using fastq files that follow the naming convention Read1: Name.fastq.1.gz
and read2: Name.fastq.2.gz.

 * a fastq file (paired end following the naming convention below)
 * a GTF geneset

The default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains
UMI/cellular barcode data and fastq.2.gz contains sequencing reads.
Chromium output is of the format: samplename_R1.fastq.gz and
samplename_R2.fastq.gz so will require conversion to the default file
format above.

Running the pipeline
--------------------

To run the pipeline you will need to set up the cluster configuration according
to the cluster documentation.

However the pipeline can also be run locally without the cluster using the
commandline flag `--no-cluster`.

The following command will run the pipeline::

   scflow singlecell make full -v5


Report generation
-----------------

The pipeline also generates Rmarkdown reports by running the following command::

   scflow singlecell make build_report -v5


output
------
