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
   
    # load transcripts        
    trlist = []
    with open(busdir+'transcripts.txt') as f:
        for line in f:
            trlist.append(line[:-3])

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

    f = open(bus_dir + "/bus_count.log", "a+")
    f.write("NUM_OF_DISTINCT_BARCODES \t %d \n"%len(indices))

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
    f.write("CBs_detected \t %d \n"%len(codewords))

    NUM_OF_UMIS_in_CELL_BARCODES=sum(values[:NUM_OF_BARCODES]
    f.write("NUM_UMIs_in_CBs \t %d \n"%NUM_OF_UMIS_in_CELL_BARCODES)
    
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
    f.write("DMIN \t %d \n"%dmin)
    f.write("NUM_CBS_TO_ERROR_CORRECT \t %d \n"%len(brc_idx_to_correct))




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

# Link barcodes to their closest codeword

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

    p = Pool(threads])
    ret_threads=p.map(merge_barcodes, barcode_split)
    p.close(); p.join()

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
    with open(busdir+'output.bus.sorted.txt') as f:
        rdr = csv.reader(f, delimiter='\t')
        for bar,umi,ec,c in rdr:
            num_tx_aligned_reads+=int(c)
            try:
                cell = b2c[bar]
                cellsets[cell].add((int(ec),umi))
            
            except KeyError: pass  
        
    # # from set back to sorted list
    # for c in range(len(cellsets)):
    #     cellsets[c]=sorted(list(cellsets[c]),key=lambda x: x[1])

    # Intersect ec/UMIs
    
    s2ec = {} #set to ec dict
    for ec in ecs.keys():
        s2ec[frozenset(ecs[ec])]=ec

    ecs_sets={i:set(ecs[i]) for i in ecs.keys()}

new_cellsets={}

    for c in range(len(cellsets)):
        new_cellsets[c]=set()
        cell = cellsets[c]
        umi_dict=defaultdict(list)
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
    

    # figure of cells before and after
    fig, ax = plt.subplots(2,2,sharex=True)
    ax = ax.ravel()
    cnt=0
    for c in np.random.randint(0,len(codewords),4):
        labels, values = zip(*Counter([i[0] for i in cellsets[c]]).items())
 
        ecdf = ECDF([len(ecs[i]) for i in labels])
        ax[cnt-1].plot(ecdf.x,ecdf.y)
        labels, values = zip(*Counter([i[0] for i in new_cellsets[c]]).items())
 
        ecdf = ECDF([len(ecs[i]) for i in labels])
        ax[cnt].plot(ecdf.x,ecdf.y)
        ax[cnt].xlim([0,30])
        ax[cnt].title( 'cell barcode: {:1}'.format(codewords[c]))
        ax[cnt].xlabel('x: number of tx')
        ax[cnt].ylabel('Pr(ec_size < x)')
        ax[cnt].legend(['before','after'])
        fig.savefig( bus_dir + '/Figures.dir/ec_intersection_example_cells.png')
        cnt+=1

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
        labels, values = zip(*Counter([i[0] for i in new_cellsets[c]]).items())
        col+= [c] * len(labels)
        row+=[tmp[ec] for ec in labels]
        data+=values
    A=coo_matrix((data, (row, col)), shape=(len(equivalence_classes),len(codewords)))
    del row,col,data

    f.write("MEDIAN_UMI_COUNTS_TCC \t %d \n"%(int(np.median(np.array(A.sum(axis=0))[0]))))

    ec_counts=np.array(A.sum(axis=1)).T[0]
    counts_sz =np.zeros(1000)
    for ec in equivalence_classes:
        cnt=ec_counts[tmp[ec]]
        sz=len(ec2gn[ec])
        try:
            counts_sz[sz-1]+=cnt
        except:
            pass
    counts_sz=100*counts_sz/sum(counts_sz)

    f.write('PERCENT_INCREASE_READS_KEPT_IN_TCCs \t {:.1f}% \n'.format(100*(100-counts_sz[0])/counts_sz[0]))
    f.write('NUM_OF_TCC_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS \t {:.1f}% \n'.format(100*A.sum()/num_tx_aligned_reads))

    f.close()

    # Write TCC matrix

    mmwrite(bus_dir+'/matrix.tcc.mtx',A)
    
    with open(save_dir+'matrix2.ec','w') as of:
        for ec in equivalence_classes:
            transcripts = ",".join([str(i) for i in ecs[ec]])
            of.write("%s\t%s\n"%(str(ec),transcripts))

    with open(save_dir+'matrix.cells','w') as of:
        of.write('\n'.join(x + '-' + str(parameter['LANE']) for x in codewords))
        of.write('\n')

    # Gene counts
    gn2ec={}
    for ec in equivalence_classes:
        for gn in ec2gn[ec]:
            try:
                gn2ec[gn].update([ec])
            except KeyError:
                gn2ec[gn] = set([ec])


    cell_gene = defaultdict(lambda: defaultdict(float))

    for c in range(len(codewords)):
        for ec,umi in new_cellsets[c]:
            gs=ec2gn[ec]
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
    
    
    B=coo_matrix((data, (row, col)), shape=(len(genes),len(codewords)))
    del row,col,data

    f = open(bus_dir + "/bus_count.log", "a+")
    f.write('MEDIAN_UMI_COUNTs_GENE \t %d \n'%int(np.median(np.array(B.sum(axis=0))[0])))

    f.write('CONSISTENCY_GENE_COUNTS/TCC \t {:.2f}% \n'.format(100*B.sum()/A.sum()))
    # UMI deduplication RATE
    f.write('NUM_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS_GENES \t{:.1f}% \n'.format(100*B.sum()/num_tx_aligned_reads))
    f.write('NUM_DEDUPLICATED_UMIs/TRANSCRIPTOME_ALIGNED_READS_TCCs  \t{:.1f}% \n'.format(100*A.sum()/num_tx_aligned_reads))

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
    fig.savefig(bus_dir+'/Figure.dir/Mean_gene_counts.png')


    f.write('RELIABLY_DETECTED_GENES \t %d \n'%t)
    f.close()

    # Write gene matrix

    mmwrite(bus_dir +'/GCmatrix.mtx',B)

    with open(bus_dir+'/GCmatrix.genes','w') as of:
        for g in genes:
            of.write("%s\n"%g)

    with open(+'GCmatrix.cells','w') as of:
        of.write('\n'.join(x + '-' + str(parameter['LANE']) for x in codewords))
        of.write('\n')


if __name__ == "__main__":
    sys.exit(main())
