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
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import mygene
import ModuleBus2count
import csv
import pickle

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
    
    # load equivalence classes
    ecs = {}
    with open(matrix_ec) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs
    
    barcodes=np.array(pd.read_csv(bus_dir+'/output.bus.sorted.txt',delimiter='\t',usecols=[0],header=None, dtype=str)).reshape(-1,)

    counts = collections.Counter(barcodes)
    labels, values = zip(*counts.items())
    # sort the values in descending order
    indSort = np.argsort(values)[::-1]
    # rearrange the data
    labels = np.array(labels)[indSort]
    values = np.array(values)[indSort]

    indices = np.arange(len(labels))

    f = open(bus_dir + "/bus_count.log", "a+")
    f.write("NUM_OF_DISTINCT_BARCODES \t %d \n"%len(indices))

    # User specifies how many cells they expect roughly in yml file

    t=int(math.floor(exp_cells * 0.1))
    exp_values = np.where(values >= values[t] / 10)

    NUM_OF_BARCODES = np.shape(exp_values)[1]

    ## PLOT, diagram of UMIs per barcode ##
    ModuleBus2count.plot_UMI_per_barcode(indices, values, NUM_OF_BARCODES, t, bus_dir)

    # No whitelist info, otherwise would incorporate here
    codewords = labels[:NUM_OF_BARCODES]
    f.write("CBs_detected \t %d \n"%len(codewords))

    NUM_OF_UMIS_in_CELL_BARCODES=sum(values[:NUM_OF_BARCODES])
    f.write("NUM_UMIs_in_CBs \t %d \n"%NUM_OF_UMIS_in_CELL_BARCODES)
    
    ## Error correct barcodes
    dmin = 3 # Paramaterised in code, find out what it means?

    brc_idx_to_correct=[]
    ret_vec = ModuleBus2count.get_ret_vec(dmin, codewords, threads)

    brc_idx_to_correct = [i for i in ret_vec if i>=0]
    f.write("DMIN \t %d \n"%dmin)
    f.write("NUM_CBs_TO_ERROR_CORRECT \t %d \n"%len(brc_idx_to_correct))
   
    ret_threads = ModuleBus2count.get_ret_threads(threads, barcodes, codewords, brc_idx_to_correct)

    ret_vec=[]
    for idd in range(len(codewords)):
        idx_list=[]
        for t in range(len(ret_threads)):
            idx_list+=ret_threads[t][idd]
        ret_vec+=[idx_list]


    reads_per_barcode=[]
    for i in range(len(codewords)):
        reads_per_barcode+=[len(ret_vec[i])]
    
    NUM_OF_UMIS_in_error_corrected_CELL_BARCODES=sum(reads_per_barcode)
    f.write("NUM_OF_UMIS_in_CBs_AFTER_ERROR_CORRECTION \t %d \n"%NUM_OF_UMIS_in_error_corrected_CELL_BARCODES)
    f.write("NUM_OF_UMIS_RESCUED \t %d \n"%(NUM_OF_UMIS_in_error_corrected_CELL_BARCODES - NUM_OF_UMIS_in_CELL_BARCODES))
    f.write("PERCENT_UMI_INCREASE \t  {:.3}% \n".format(100*(NUM_OF_UMIS_in_error_corrected_CELL_BARCODES - NUM_OF_UMIS_in_CELL_BARCODES)/NUM_OF_UMIS_in_CELL_BARCODES)) 


    b2c={} #barcode to codeword dict
    for i in range(len(ret_vec)):
        for j in ret_vec[i]:
            b2c[barcodes[j]]=i

    cellsets = {i:set() for i in range(len(codewords))} # this set will collapse duplicate ec-umis in each cell
    num_tx_aligned_reads = 0
    with open(bus_dir+'/output.bus.sorted.txt') as r:
        rdr = csv.reader(r, delimiter='\t')
        for bar,umi,ec,c in rdr:

            num_tx_aligned_reads+=int(c)
            try:
                cell = b2c[bar]
                cellsets[cell].add((int(ec),umi))
            
            except KeyError: pass  


    # from set back to sorted list
    #for c in range(len(cellsets)):
     #    cellsets[c]=sorted(list(cellsets[c]),key=lambda x: x[1])

    # Intersect ec/UMIs
    
    s2ec = {} #set to ec dict
    for ec in ecs.keys():
        s2ec[frozenset(ecs[ec])]=ec

    ecs_sets={i:set(ecs[i]) for i in ecs.keys()}

    new_cellsets={}

    for c in range(len(cellsets)):
        new_cellsets[c]=set()
        cell = cellsets[c]
        umi_dict=collections.defaultdict(list)
        for i in cell:
            umi_dict[i[1]].append(ecs_sets[i[0]])

        for umi in umi_dict.keys(): 
            ec_list=sorted(umi_dict[umi],key=len,reverse=1)
            k=len(ec_list)
            done = False
            new_ec = set.intersection(*ec_list)
            if len(new_ec) > 0:
                ec_list=[new_ec];
                done = True;     
            while(not done):
                done=True
                for i,j in combinations(range(k), 2): 

                    new_ec = ec_list[i].intersection(ec_list[j])
                    if len(new_ec) > 0:
                        msk=np.ones(k,dtype=bool);msk[[i,j]]=False
                        ec_list= [u for u in compress(ec_list, msk)]+[new_ec]
                        k-=1;done=False; break;
        
            for new_ec in ec_list:
                try: 
                    new_cellsets[c].add((s2ec[frozenset(new_ec)],umi))
                    #print((s2ec[frozenset(new_ec)],umi))
                except KeyError:
                    new_id=len(ecs)
                    listoftx=list(new_ec);listoftx.sort()
                    ecs[new_id] = listoftx
                    s2ec[frozenset(new_ec)] = new_id
                    new_cellsets[c].add((s2ec[frozenset(new_ec)],umi))
    
    ## PLOT, figure of cells before and after ##
    ModuleBus2count.plot_cell_before_after(codewords, cellsets, ecs, labels, new_cellsets, bus_dir)

    # Get TCC matrix
    ec_set=set()
    for c in range(len(codewords)):
        ec_set.update([i[0] for i in new_cellsets[c]])
    
    equivalence_classes=list(ec_set)
    equivalence_classes.sort()


    tmp={}
    for id in range(len(equivalence_classes)):
        tmp[equivalence_classes[id]] = id

    row=[]
    col=[]
    data=[]

    for c in range(len(codewords)):
        labels, values = zip(*collections.Counter([i[0] for i in new_cellsets[c]]).items())
        col+= [c] * len(labels)
        row+=[tmp[ec] for ec in labels]
        data+=values
    A=coo_matrix((data, (row, col)), shape=(len(equivalence_classes),len(codewords)))
    del row,col,data
    A_large = A.toarray()

    f.write("MEDIAN_UMI_COUNTS_TCC \t %d \n"%(int(np.median(np.array(A.sum(axis=0))[0]))))


    def g2n_dict(ENSGLIST, species_code):
        mg = mygene.MyGeneInfo()
        ginfo = mg.querymany(ENSGLIST, scopes='ensembl.gene',returnall=True, species = species_code )

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
   
 
    # load transcripts        
    trlist = []
    with open(bus_dir+'/transcripts.txt') as f:
        for line in f:
            trlist.append(line.rstrip('\n'))

    # Dictionaries for transcript to gene and for gene to gene symbol
    tr2g = t2g_dict(t2gmap)

    # 15 for human, 18 for mouse
    len_of_ens = len(list(tr2g.keys())[0])
    if len_of_ens == 15:
        species = 9606
    elif len_of_ens == 18:
        species = 10090
    
    ENSGLIST = list(np.unique(list(tr2g.values())))
    #g2n = g2n_dict(ENSGLIST, species) 

    # ec to gene names  
    #ec2gn = {ec:frozenset([item for sublist in  [g2n[tr2g[trlist[t][:len_of_ens]]] for t in ecs[ec]]   for item in sublist ]) for ec in equivalence_classes}
    # EC to ensembl id    
    ec2gn = {ec:frozenset([item for sublist in  [tr2g[trlist[t][:len_of_ens]] for t in ecs[ec]]   for item in sublist ]) for ec in equivalence_classes}
    ec2id = {ec:frozenset([tr2g[trlist[t]] for t in ecs[ec]]) for ec in equivalence_classes}

    # Stats
    ec_counts=np.array(A.sum(axis=1)).T[0]
    counts_sz =np.zeros(1000)
    for ec in equivalence_classes:
        cnt=ec_counts[tmp[ec]]
        sz=len(ec2id[ec])
        try:
            counts_sz[sz-1]+=cnt
        except:
            pass
    counts_sz=100*counts_sz/sum(counts_sz)

    f = open(bus_dir + "/bus_count.log", "a+")
    f.write('PERCENT_INCREASE_READS_KEPT_IN_TCCs \t {:.1f}% \n'.format(100*(100-counts_sz[0])/counts_sz[0]))
    f.write('NUM_OF_TCC_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS \t {:.1f}% \n'.format(100*A.sum()/num_tx_aligned_reads))

    f.close()

    # Write TCC matrix

    mmwrite(bus_dir+'/matrix.tcc.coord.mtx',A)
    
    with open(bus_dir+'/matrix_updated.ec','w') as of:
        for ec in equivalence_classes:
            transcripts = ",".join([str(i) for i in ecs[ec]])
            of.write("%s\t%s\n"%(str(ec),transcripts))

    with open(bus_dir+'/matrix.cells','w') as of:
        of.write('\n'.join(codewords))
        of.write('\n')

    # Gene counts
    gn2ec={}
    for ec in equivalence_classes:
        for gn in ec2id[ec]:
            try:
                gn2ec[gn].update([ec])
            except KeyError:
                gn2ec[gn] = set([ec])


    cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))

    for c in range(len(codewords)):
        for ec,umi in new_cellsets[c]:
            gs=ec2id[ec]
            if len(gs)==1:
                cell_gene[c][gs] += 1

    gene_set=set()
    for c in range(len(codewords)):
        gene_set.update(cell_gene[c])
    
    genes=list(gene_set)

    tmp2={}
    for id in range(len(genes)):
        tmp2[genes[id]] = id

    row=[]
    col=[]
    data=[]

    for c in range(len(codewords)):
        row+=[tmp2[gene] for gene,_ in cell_gene[c].items()]
        data+=[round(val) for _,val in cell_gene[c].items()]
        col+= [c] * len(cell_gene[c])
    
    
    B = coo_matrix((data, (row, col)), shape=(len(genes),len(codewords)))
    del row,col,data
    #B_large = B.toarray()
    #B_dense = B.todense()

    f = open(bus_dir + "/bus_count.log", "a+")
    f.write('MEDIAN_UMI_COUNTs_GENE \t %d \n'%int(np.median(np.array(B.sum(axis=0))[0])))

    f.write('CONSISTENCY_GENE_COUNTS/TCC \t {:.2f}% \n'.format(100*B.sum()/A.sum()))
    # UMI deduplication RATE
    f.write('NUM_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS_GENES \t{:.1f}% \n'.format(100*B.sum()/num_tx_aligned_reads))
    f.write('NUM_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS_TCCs  \t{:.1f}% \n'.format(100*A.sum()/num_tx_aligned_reads))
   
    t=np.sum(np.array(B.mean(axis=1))>0.1)

    ## PLOT mean_gene_counts ##
    ModuleBus2count.plot_mean_gene_counts(B,t, bus_dir)
 
    f.write('RELIABLY_DETECTED_GENES \t %d \n'%t)
    f.close()

    ## Write gene matrix ##

    # Coordinates
    mmwrite(outfile, B)
    
    regex_match = "frozenset\(\{\'(.*)\'\}\)"
    with open(bus_dir+'/GCmatrix.genes','w') as of:
        for g in genes:
            print(g)
            gname_re = re.search(regex_match, str(g))
            gene_name = gname_re.group(1)

            of.write("%s\n"%gene_name)
             
    with open(bus_dir +'/GCmatrix.cells','w') as of:
        of.write('\n'.join(codewords))
        
   # Save full array
   # np.save(outfile, B_large.int32)

if __name__ == "__main__":
    sys.exit(main())
