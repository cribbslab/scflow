#!/usr/bin/env python

"""
ModuleVelo.py - Tasks for running single cell velocyto pipeline.

"""
import os
import cgatcore.pipeline as P


def check_paired_end(fastqfiles):
    "checks if data is paired end or single end"
    fastqfile1 = []
    fastqfile2 = []
    for fastqfile in fastqfiles:
        if fastqfile.endswith(".fastq.2.gz"):
            bn = P.snip(fastqfile, ".fastq.2.gz")
            infile1 = "%s.fastq.1.gz" % bn
            infile2 = "%s.fastq.2.gz" % bn
            if not os.path.exists(infile1):
                raise ValueError("cant find paired end file "
	                     "'%s' for '%s'" % (infile1, infile2))
            fastqfile1.append(infile1)
            fastqfile2.append(infile2)
        else:
            raise ValueError("Alevin requires UMI/CB file and reads file")

    fastqfiles = [fastqfile1, fastqfile2]
    return(fastqfiles)

def check_multiple_read_files(infiles):
    if isinstance(infiles[0], tuple):
        index = infiles[0][1]

        fastqs = [x[0] for x in infiles]

    else: 
        fastqs = [infiles[0]]
        index = infiles[1]

    output = [fastqs, index]
    return(output)

def check_multiple_read_files_no_index(infiles):
    if isinstance(infiles, tuple):
        fastqs = [x for x in infiles]

    else: 
        fastqs = [infile]


    output = fastqs
    return(output)
