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

if __name__ == "__main__":
    sys.exit(main())
