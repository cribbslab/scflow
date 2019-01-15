#!/usr/bin/env python

"""
ModuleSC.py - Tasks for running single cell pipleine

"""
import os
import cgatcore.pipeline as P


def check_paired_end(fastqfile):
    "checks if data is paired end or single end"
    if fastqfile.endswith(".fastq.1.gz"):
        bn = P.snip(fastqfile, ".fastq.1.gz")
        infile1 = "%s.fastq.1.gz" % bn
        infile2 = "%s.fastq.2.gz" % bn
        if not os.path.exists(infile2):
            raise ValueError("cant find paired end file "
                             "'%s' for '%s'" % (infile, infile2))
        fastqfile = [infile1, infile2]
        return fastqfile
    else:
        return fastqfile
