'''
bus2count.py - generate count matrices from kallisto BUS output
=================================================================
:Tags: Python

Purpose
-------
Takes input of a sorted bus text files, an ec matrix and transcript text file and generates
count matrices which can be converted to SCE further downstream. 

Usage
-----
Example::
   python bus2count.py --path

Command line options
--------------------
--help (-h) help for the script
--busdir (-b) the path to the folder of the kallisto bus output for the given sample

'''

import os 
import sys
import re
import pandas as pd
import collections
import math
import numpy as np

import cgatcore.experiment as E



def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--busdir", dest="bus_dir", type ="string" ,
                      help="directory bus output is located within")

    parser.add_option("-t", "--t2gmap", dest="t2gmap", type ="string" ,
                      help="Transcript to gene map file")
   
    parser.set_defaults(bus_dir=None, t2gmap = 'transcript2geneMap.tsv')

    (options, args) = E.start(parser)

    bus_dir = options.bus_dir
    t2gmap = options.t2gmap

    matrix_ec = bus_dir + "/matrix.ec"
    transcripts = bus_dir + "/transcripts.txt"
    sorted_text = bus_dir + "/output.bus.sorted.txt"

    tr2g = {}
    trlist = []
    with open(t2gmap) as f:
        next(f)
        for line in f:
            l = line.split()
            tr2g[l[0]] = l[1]
            trlist.append(l[0])

    genes = list(set(tr2g[t] for t in tr2g))

    # Equivalence classes
    ecs = {}
    with open(matrix_ec) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs
        
    def ec2g(ec):
        if ec in ecs:
            return list(set(tr2g[trlist[t]] for t in ecs[ec]))        
        else:
            return []
    # if parameter['whitelist']!=0:
    #whitelist = set(x.strip() for x in open(parameter['whitelist']))

    barcodes=np.array(pd.read_csv(sorted_text ,delimiter='\t',usecols=[0],header=None, dtype=str)).reshape(-1,)

    counts = collections.Counter(barcodes)
    labels, values = zip(*counts.items())
    # sort the values in descending order
    indSort = np.argsort(values)[::-1]
    # rearrange the data
    labels = np.array(labels)[indSort]
    values = np.array(values)[indSort]

    indices = np.arange(len(labels))
    print("NUM_OF_DISTINCT_BARCODES =",len(indices))

    cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
    pbar=None
    pumi=None
    with open(sorted_text) as f:
        gs = set()
        for line in f:
            l = line.split()
            barcode,umi,ec,count = line.split()
            ec = int(ec)
        
            if barcode == pbar:
                # same barcode
                if umi == pumi:
                    # same UMI, let's update with intersection of genelist
                    gl = ec2g(ec)
                    gs.intersection_update(gl)
                else:
                    # new UMI, process the previous gene set
                    for g in gs:
                        cell_gene[barcode][g] += 1.0/len(gs)
                        # record new umi, reset gene set
                        pumi = umi
                        gs = set(ec2g(ec))
            else:
            # work with previous gene list
                for g in gs:
                    cell_gene[pbar][g] += 1.0/len(gs)
            
                if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                    del cell_gene[pbar]
            
                pbar = barcode
                pumi = umi
            
                gs = set(ec2g(ec))
        #remember the last gene
        for g in gs:
            cell_gene[pbar][g] += 1.0/len(gs)
        
        #if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
        #    del cell_gene[pbar]

    print("MATRIX??", cell_gene)
    print(pbar)
    print(pumi)

if __name__ == "__main__":
    sys.exit(main())
