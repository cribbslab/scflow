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

    parser.add_option("-e", "expectedcells", dest="exp_cells", type ="int" ,
                      help="Expected number of cells")

    parser.add_option( "threads", dest="threads", type ="int" ,
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


'''
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
        
        if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
            del cell_gene[pbar]

    barcode_hist = collections.defaultdict(int)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        s = len([cg[g] for g in cg])
        barcode_hist[barcode] += s


    bcv = [x for b,x in barcode_hist.items() if x > 10 and x < 4000]
    plt.hist(bcv,bins=100)
    # Better wording needed for labels
    plt.xlabel('Number of cells with the same barcode')
    plt.ylabel('Frequency')
    plt.savefig('kallisto.dir/barcode_histogram.png')

    # Paramaterise based on how histogram looks. User can see plot and then run again based on plot
        
    bad_barcode = [x for x in barcode_hist if  barcode_hist[x] <= barcode_thresh]

    s = 0
    bad_s = 0
    bad_barcode_set = set(bad_barcode)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        cgs =  sum(cg[g] for g in cg)
        s += cgs
        if barcode in bad_barcode_set:
            bad_s += cgs

    gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
    barcodes_to_use = [b for b,x in barcode_hist.items() if x > 0]

    num_entries = 0
    for barcode in barcodes_to_use:
        num_entries += len([x for x in cell_gene[barcode].values() if round(x)>0])

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

'''
import mygene

    def g2n_dict():
        mg = mygene.MyGeneInfo()
        ginfo = mg.querymany(ENSGLIST, scopes='ensembl.gene',returnall=True)

        g2n = {}
        count_exept=0
        for g in ginfo['out']:
            try:
                gene_id=str(g['query'])
                gene_name=str(g['symbol'])
        
                g2n[gene_id] = g2n.get(gene_id, [])
                g2n[gene_id].append(str(g['symbol']))                
            except KeyError:
                count_exept+=1
                g2n[ str(g['query']) ] = [str(g['query'])]
        return(g2n)

    def t2g_dict(infile):
        d={}
        with open(infile) as f:
            next(f)
            for line in f:
                (key, value)=line.split()
                d[key]=value
        return(d)

    # Dictionaries for transcript to gene and for gene to gene symbol
    g2n = g2n_dict()
    tr2g = t2g_dict(t2gmap)

    # ec to gene names 
    ec2gn = {ec:frozenset([item for sublist in  [g2n[tr2g[trlist[t][:len_of_ens]]] for t in ecs[ec]]   for item in sublist ]) for ec in equivalence_classes}

    # load equivalence classes
    ecs = {}
    with open(matrix_ec) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs

    barcodes=np.array(pd.read_csv(busdir+'output.bus.sorted.txt',delimiter='\t',usecols=[0],header=None, dtype=str)).reshape(-1,)

    counts = Counter(barcodes)
    labels, values = zip(*counts.items())
    # sort the values in descending order
    indSort = np.argsort(values)[::-1]
    # rearrange the data
    labels = np.array(labels)[indSort]
    values = np.array(values)[indSort]

    indices = np.arange(len(labels))

    # User specifies how many cells they expect roughly in yml file

    t=int(math.floor(exp_cells * 0.1))
    exp_values = np.where(values >= values[t] / 10)

    NUM_OF_BARCODES = np.shape(exp_values)[1]

    # Plot diagram of UMIs per barcode

    fig, ax = plt.subplots()
    ax.plot(indices, (values))

    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylabel('Number of reads', color='k')
    ax.set_xlabel('barcode', color='k')
    ax.set_title('Number of umis per Barcode', color='b')
    ax.axvline(t, color='gray', linestyle='--',linewidth=.5)
    ax.axvline(NUM_OF_BARCODES, color='g', linestyle='--',linewidth=2.0)
    ax.axhline(values[t], color='gray', linestyle='--',linewidth=.5)
    ax.axhline(values[t]/10, color='red', linestyle='--',linewidth=.5)
    fig.savefig( bus_dir + '/Figures.dir/UMI_barcodes.png')

    # No whitelist info, otherwise would incorporate here
    codewords = labels[:NUM_OF_BARCODES]
    NUM_OF_UMIS_in_CELL_BARCODES=sum(values[:NUM_OF_BARCODES]

    ## Error correct barcodes
    dmin = 3 # Paramaterised in code, find out what it means?

    def hamdist(s1, s2):
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

    def is_far_enough(i):
        d=100000
        codi=codewords[i]
        j_range=list(np.arange(1,i))+list(np.arange(i+1,len(codewords)))
        while d>=dmin and len(j_range)>0:
            j=j_range.pop()
            d=np.min([d,hamdist(codi,codewords[j])])

        return i if d>=dmin else -1

    brc_idx_to_correct=[]
    p=Pool(parameter['NUM_THREADS'])
    ret_vec=p.map(is_far_enough, np.arange(len(codewords)) )
    p.close()
    p.join()

    brc_idx_to_correct = [i for i in ret_vec if i>=0]
    print("-- number of cell barcodes to error-correct:", len(brc_idx_to_correct), "( dmin >=", dmin,")" )




    def hamming_circle(s, n, alphabet='ATCG'):
        """Generate strings over alphabet whose Hamming distance from s is
        exactly n.
        """
        for positions in combinations(range(len(s)), n):
            for replacements in product(range(len(alphabet) - 1), repeat=n):
                cousin = list(s)
                for p, r in zip(positions, replacements):
                    if cousin[p] == alphabet[r]:
                        cousin[p] = alphabet[-1]
                    else:
                        cousin[p] = alphabet[r]
                yield ''.join(cousin)


def merge_barcodes(barcs):
    offset=barcs[0]
    barcs=barcs[1]
    retvec=[]
    for idd in range(len(codewords)):
        retvec+=[[]]
    for idx, barcode in enumerate(barcs):
        if barcode in codeword_set: retvec[cw[barcode]] +=[idx+offset]
        else:
            if barcode in brc_to_correct_neigbors:
                neighbors = hamming_circle(barcode,1)
                for neighbor in neighbors:
                    if neighbor in brc_to_correct: retvec[cw[neighbor]] +=[idx+offset]; break;
    return retvec

if __name__ == "__main__":
    sys.exit(main())
