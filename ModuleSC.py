#!/usr/bin/env python

"""
ModuleSC.py - Tasks for running single cell pipleine

"""
import os
import cgatcore.pipeline as P


def check_paired_end(fastqfiles):
    "checks if data is paired end or single end"
    fastqfile1 = []
    fastqfile2 = []
    for fastqfile in fastqfiles:
        if fastqfile.endswith(".fastq.1.gz"):
            bn = P.snip(fastqfile, ".fastq.1.gz")
            infile1 = "%s.fastq.1.gz" % bn
            infile2 = "%s.fastq.2.gz" % bn
            if not os.path.exists(infile2):
                raise ValueError("cant find paired end file "
                             "'%s' for '%s'" % (infile, infile2))
            fastqfile1.append(infile1)
            fastqfile2.append(infile2)
        else:
            raise ValueError("Alevin requires UMI/CB file and reads file")

    fastqfiles = [fastqfile1, fastqfile2]
    return(fastqfiles)


def check_multiple_read_files(infiles):
    if isinstance(infiles[0], tuple):
        index = infiles[0][1]
        t2gmap = infiles[0][2]
        fastqs = [x[0] for x in infiles]
        #fastqs = " ".join(fastqs)

    else: 
        fastqs = [infiles[0]]
        index = infiles[1]
        t2gmap = infiles[2]
    
    output = [fastqs, index, t2gmap]
    return(output)
