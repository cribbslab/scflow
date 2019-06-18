import os 
import sys
import re
import pandas as pd
import collections
import math
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from itertools import chain, combinations, product,compress

import cgatcore.experiment as E

def main(argv=sys.argv):

##################
# Option parsing 
##################

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--dir", dest="bus_dir", type ="string" ,
                      help="directory bus output is located within")

    parser.add_option("-t", "--t2gmap", dest="t2gmap", type ="string" ,
                      help="Transcript to gene map file")

    parser.add_option("-b", "--barcodethresh", dest="thresh", type ="int" ,
                      help="Minimum threshold for cellular barcodes")   

    parser.add_option("-o", "--out", dest="outfile", type ="string" ,
                      help="Gene count matrix outfile")

    parser.add_option("-e", "--expectedcells", dest="exp_cells", type ="int" ,
                      help="Expected number of cells")

    parser.add_option( "--threads", dest="threads", type ="int" ,
                      help="Number of threads")
    

    parser.set_defaults(bus_dir=None, t2gmap = 'transcript2geneMap.tsv', thresh = 100, outfile = 'kallisto.dir/output.bus.mtx', exp_cells = 1000, threads = 1)

    (options, args) = E.start(parser)

    bus_dir = options.bus_dir
    t2gmap = options.t2gmap
    barcode_thresh = options.thresh
    outfile = options.outfile
    exp_cells = options.exp_cells
    threads = options.threads

    # Bus files

    matrix_ec = bus_dir + "/matrix.ec"
    transcripts = bus_dir + "/transcripts.txt"
    sorted_text = bus_dir + "/output.bus.sorted.txt"


    def t2g_dict(infile):
        d={}
        with open(infile) as f:
            next(f)
            for line in f:
                (key, value)=line.split()
                d[key]=value
        return(d)
 
    # load transcripts        
    trlist = []
    with open(bus_dir+'/transcripts.txt') as f:
        for line in f:
            trlist.append(line.rstrip('\n'))

    # Dictionaries for transcript to gene and for gene to gene symbol
    tr2g = t2g_dict(t2gmap)
 
    # load equivalence classes
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

    # load kallisto bus output dataset

    cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
    pbar=None
    pumi=None
    with open(bus_dir+'/output.bus.sorted.txt') as f:
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
        
        if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
            del cell_gene[pbar]

    barcode_hist = collections.defaultdict(int)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        s = len([cg[g] for g in cg])
        barcode_hist[barcode] += s

    threshold = 0 # this filters the data by gene count
    bcv = [x for b,x in barcode_hist.items() if x > threshold] 
    fig, ax = plt.subplots()
    ax.hist(bcv,bins=40, log=True)
    
    ax.set_ylabel('Number barcodes', color='k')
    ax.set_xlabel('Number of gene counts', color='k')

    fig.savefig(bus_dir+'/Number_barcodes_vs_gene_counts.png')


    barcodes_to_use = [b for b,x in barcode_hist.items() if x > 0 ]

    num_entries = 0
    for barcode in barcodes_to_use:
        num_entries += len([x for x in cell_gene[barcode].values() if round(x)>0])

    genes = list(set(tr2g[t] for t in tr2g))
    gene_to_id = dict((g,i+1) for i,g in enumerate(genes))

    with open(outfile, 'w') as of:
        of.write('%%MatrixMarket matrix coordinate real general\n%\n')
        #number of genes
        of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), num_entries))
        bcid = 0
        for barcode in barcodes_to_use:
            bcid += 1
            cg = cell_gene[barcode]
            gl = [(gene_to_id[g],round(cg[g])) for g in cg if round(cg[g]) > 0]
            gl.sort()
            for x in gl:
                of.write("%d %d %d\n"%(x[0],bcid,x[1]))



if __name__ == "__main__":
    sys.exit(main())

