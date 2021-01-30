from scipy import sparse
import anndata
import scipy
from scipy import io
import os
import pandas as pd
import argparse
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", dest="output", type=str,
                        help="Output directory")

parser.add_argument("-s", "--spliced", dest="spliced", type=str,
                        help="The spliced .mtx file")

parser.add_argument("-c", "--sp-barcode", dest="spliced_barcode", type=str,
                        help="The spliced barcodes")

parser.add_argument("-a", "--sp-genes", dest="spliced_genes", type=str,
                        help="The spliced genes")

parser.add_argument("-u", "--unspliced", dest="unspliced", type=str,
                        help="The unspliced .mtx file")

parser.add_argument("-b", "--unsp-barcode", dest="unspliced_barcode", type=str,
                        help="The unspliced barcodes")

parser.add_argument("-g", "--unsp-genes", dest="unspliced_genes", type=str,
                        help="The unspliced genes")
parser.set_defaults(
        output=None,
    )

args = parser.parse_args()

## load unspliced data on anndata as sparse crs matrix
unspliced = anndata.AnnData(scipy.io.mmread(args.unspliced).tocsr())
unspliced_csr = scipy.io.mmread(args.unspliced).tocsr()
unspliced.obs = pd.read_csv(args.unspliced_barcode, index_col = 0, header = None, names = ['barcode'])
unspliced.var = pd.read_csv(args.unspliced_genes, header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded unspliced count matrix.')

## load unspliced data on anndata as sparse crs matrix
spliced = anndata.AnnData(scipy.io.mmread(args.spliced).tocsr())
spliced_csr = scipy.io.mmread(args.spliced).tocsr()
spliced.obs= pd.read_csv(args.spliced_barcode, index_col = 0, header = None, names = ['barcode'])
spliced.var = pd.read_csv(args.spliced_genes, header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded spliced count matrix')

# Now that we have spliced and unspliced matrices we can sum the counts of genes for barcodes common to both matrices We take the intersection of both matrices
# because presumably cells without a single count on either have very low counts anyway
idx = spliced.obs.index.intersection(unspliced.obs.index)
spliced_intersection = spliced[idx]
unspliced_intersection = unspliced[idx]
spliced_intersection.X + unspliced_intersection.X

spliced_plus_unspliced = spliced_intersection.copy()
spliced_plus_unspliced.X = spliced_intersection.X + unspliced_intersection.X


# Use scipy to write the matrix to .mtx file

#spliced_plus_unspliced = spliced_plus_unspliced.to_df
io.mmwrite(args.output + "/genes.mtx", spliced_plus_unspliced.X)
barcodes_merged = spliced_plus_unspliced.obs_names
genes_merged = spliced_plus_unspliced.var_names

out_barcode = open(args.output + "/genes.barcodes.txt", "w")
for i in barcodes_merged:
    out_barcode.write("%s\n"% (i))

out_gene = open(args.output + "/genes.genes.txt", "w")

for i in genes_merged:
    out_gene.write("%s\n"% (i))

out_gene.close()
out_barcode.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
