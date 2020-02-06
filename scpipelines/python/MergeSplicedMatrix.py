#!/usr/bin/env python


from anndata import AnnData
import anndata
from scipy import sparse
import scipy
import anndata
import scipy.io
import os
import pandas as pd
import scanpy as sc
import cgatcore.experiment as E


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-d", "--dir-bus", dest="bus_path", type=str,
                        help="The path to the bus directory")

    parser.add_argument("-o", "--output", dest="output", type=str,
                        help="Output directory")

    parser.add_argument("-s", "--spliced", dest="spliced", type=str,
                        help="The spliced .mtx file")

    parser.add_argument("-c", "--sp-barcode", dest="spliced_barcode", type=str,
                        help="The unspliced .mtx file")

    parser.add_argument("-h", "--sp-genes", dest="spliced_genes", type=str,
                        help="The unspliced .mtx file")

    parser.add_argument("-u", "--unspliced", dest="unspliced", type=str,
                        help="The unspliced .mtx file")

    parser.add_argument("-b", "--unsp-barcode", dest="unspliced_barcode", type=str,
                        help="The unspliced .mtx file")

    parser.add_argument("-g", "--unsp-genes", dest="unspliced_genes", type=str,
                        help="The unspliced .mtx file")
    parser.set_defaults(
        bus_path=None,
        output=None,
    )

    (args) = E.start(parser)

## load unspliced data on anndata as sparse crs matrix
unspliced = anndata.AnnData(scipy.io.mmread(unspliced).tocsr())
unspliced.obs= pd.read_csv(unspliced_barcode, index_col = 0, header = None, names = ['barcode'])
unspliced.var = pd.read_csv(unspliced_genes, header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded unspliced count matrix.')

## load unspliced data on anndata as sparse crs matrix
spliced = anndata.AnnData(scipy.io.mmread(spliced).tocsr())
spliced.obs= pd.read_csv(spliced_barcode, index_col = 0, header = None, names = ['barcode'])
spliced.var = pd.read_csv(spliced_genes, header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded spliced count matrix')

# Now that we have spliced and unspliced matrices we can sum the counts of genes for barcodes common to both matrices We take the intersection of both matrices
# because presumably cells without a single count on either have very low counts anyway
idx = spliced.obs.index.intersection(unspliced.obs.index)
spliced_intersection = spliced[idx]
spliced_intersection = unspliced[idx]
spliced_intersection.X + unspliced_intersection.X

spliced_plus_unspliced = spliced_intersection.copy()
spliced_plus_unspliced.X = spliced_intersection.X + unspliced_intersection.X


# Use scipy to write the matrix to .mtx file
spliced_plus_unspliced.write(outfile)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
