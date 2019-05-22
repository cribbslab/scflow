''' 
ModuleBus2count.py - Helper functions for bus2count.py 
======================================================
:Tags: Python

Purpose
--------
Modules for bus2count, which processes bus files into gene count matrices.
Parallelisation parts in separate file using multiprocessing.

Usage
-----
Called by bus2count.py
Not to be ran using the command line

'''

import numpy as np
import pandas as pd
import math
from multiprocessing import Pool, Lock
from functools import partial
from itertools import chain, combinations, product, compress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from statsmodels.distributions.empirical_distribution import ECDF

# Get rec_vec, to detect number of CBs to error correct
def hamdist(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def is_far_enough(codewords, dmin, i):
    d=100000
    codi=codewords[i]
    j_range=list(np.arange(1,i))+list(np.arange(i+1,len(codewords)))
    while d>=dmin and len(j_range)>0:
        j=j_range.pop()
        d=np.min([d,hamdist(codi,codewords[j])])

    return i if d>=dmin else -1

def get_ret_vec(dmin, codewords, threads):
    p = Pool(threads)
    i = np.arange(len(codewords))
    func = partial(is_far_enough, codewords, dmin)
    ret_vec = p.map(func, i)

    p.close()
    p.join()
    return(ret_vec)

# ret_threads. Merge barcodes and barcode split
# ------------------------------------------------
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


def merge_barcodes(codewords, codeword_set, brc_to_correct_neigbors, brc_to_correct, cw, barcs):
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


def get_ret_threads(threads, barcodes, codewords, brc_idx_to_correct):

    chunksize=1+int(len(barcodes)/threads)
    
    cw={}
    for id in range(len(codewords)):
        cw[codewords[id]] = id
    
    barcode_split=[]
    for i in range(0, len(barcodes), chunksize):
        barcode_split+=[[i,barcodes[i:i+chunksize]]]
    
    codeword_set = set(codewords)
    codeword_list = list(codewords)
    brc_to_correct=set(codewords[brc_idx_to_correct])

    #### Generate the set of all dist-1 neighbors of brc_to_correct (for fast check in merge func)
    #### note: the number of barcodes in this set is len(brc_to_correct)*3*barcode_length
    brc_to_correct_neigbors=set()
    for brc in brc_to_correct:
        neighbors = hamming_circle(brc,1)
        for neighbor in neighbors:
            brc_to_correct_neigbors.add(neighbor)

    # Parallelise 
    p = Pool(threads)
    func = partial(merge_barcodes, codewords, codeword_set, brc_to_correct_neigbors, brc_to_correct, cw)
    ret_threads=p.map(func, barcode_split)

    p.close()
    p.join()
    return(ret_threads)

## Plots ##

# UMIs per barcode (knee plot)
def plot_UMI_per_barcode(indices, values, NUM_OF_BARCODES, t, bus_dir):
    
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
    fig.savefig( bus_dir + '/UMI_barcodes.png')

def plot_cell_before_after(codewords, cellsets, ecs, labels, new_cellsets, bus_dir):

    # figure of cells before and after
    fig, ax = plt.subplots(2,2,sharex=True)
    ax = ax.ravel()
    cnt=0
    for c in np.random.randint(0,len(codewords),4):
        labels, values = zip(*collections.Counter([i[0] for i in cellsets[c]]).items())
 
        ecdf = ECDF([len(ecs[i]) for i in labels])
        ax[cnt].plot(ecdf.x,ecdf.y)
        labels, values = zip(*collections.Counter([i[0] for i in new_cellsets[c]]).items())
 
        ecdf = ECDF([len(ecs[i]) for i in labels])
        ax[cnt].plot(ecdf.x,ecdf.y)
        plt.xlim([0,30])
        ax[cnt].set_title( 'cell barcode: {:1}'.format(codewords[c]))
        ax[cnt].set_xlabel('x: number of tx')
        ax[cnt].set_ylabel('Pr(ec_size < x)')
        ax[cnt].legend(['before','after'])
        cnt+=1

    fig.suptitle("Equivalence classes for 4 random CBs", fontsize=14)
    fig.savefig(bus_dir + '/ec_intersection_example_cells.png')


def plot_mean_gene_counts(B,t, bus_dir):

    fig, ax = plt.subplots()
    t=np.sum(np.array(B.mean(axis=1))>0.1)
    ax.grid()
    ax.plot(np.sort(np.array(B.mean(axis=1)),axis=0).T[0][::-1],color='r')
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylabel('average umi counts', color='k')
    ax.set_xlabel('genes', color='k')
    ax.set_title('reliably detected genes ('+str(t)+')', color='darkgreen')
    ax.axhline(0.1, color='firebrick', linestyle='--',linewidth=1)
    ax.axvline(t, color='firebrick', linestyle='--',linewidth=1)
    ax.plot(np.sort(np.array(B.mean(axis=1)),axis=0).T[0][::-1][:t],color='green',linewidth=2)
    fig.savefig(bus_dir+'/Mean_gene_counts.png')

if __name__ == "__main__":
    sys.exit(main())
