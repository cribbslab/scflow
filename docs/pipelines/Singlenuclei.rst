
===========
kallistobus
===========

Overview
========

This pipeline performs alignment free based quantification using the kallisto bus tools
pipeline. For further information on the code that is wrapped up in this pipeline please
refer to to the following `documentation <https://github.com/pachterlab/kb_python>`_

This pipeline runs kallisto bustools and runs alignment, counting and outputs either a 
loom or h5ad file that can be further imported into Seurat or Scanpy.

You have the option of running the following quantification:
* standard - will generate a standard pseudoaligned analysis output
* lamanno -  for RNA- velocity analysis based on la Manno et al 2018
* kite - for feature barcoding
* kite:10xFB - feature barcoding for 10X genomics 

A number of technologies are supported by kallisto bustools:
* 10XV1
* 10XV2
* 10XV3
* CELSEQ
* CELSEQ2
* DROPSEQ
* INDROPSV1
* INDROPSV2
* INDROPSV3
* SCRUBSEQ
* SURECELL


All of these options are set within the pipeline.yml configuration file. Please see
more info below.

Input files
-----------

The pipeline is ran using fastq files that follow the naming convention
Read1: Name.fastq.1.gz and read2: Name.fastq.2.gz.

 * a fastq file (paired end following the naming convention below)
 * a GTF geneset

The default file format assumes the following convention:
fastq.1.gz and fastq.2.gz for paired data, where fastq.1.gz contains
UMI/cellular barcode data and fastq.2.gz contains sequencing reads.
Chromium output is of the format: samplename_R1.fastq.gz and
samplename_R2.fastq.gz so will require conversion to the default file
format above.

Configuring the pipeline
------------------------

To set the input values for the pipeline you need to modify a configuration
file. To generate this yml file run the following::

   scflow kallistobus config -v5

Then open up the pipeline.yml file and modify the default values before running the pipeline.

Running the pipeline
--------------------

To run the pipeline you will need to set up the cluster configuration according
to the cluster documentation.

However the pipeline can also be run locally without the cluster using the
commandline flag `--no-cluster`.

The following command will run the pipeline::

   scflow kallistobus make full -v5



output
------

The output of the piepline is a bus record for each sample than can be combined using downstream tools
such as Seurat or Scanpy.
