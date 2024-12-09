#!/usr/bin/env python

import argparse
import os
import pandas as pd
from scipy import io
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description="Convert mtx files to h5ad")
    parser.add_argument('--matrix', required=True, help='Path to genes.mtx file')
    parser.add_argument('--genes', required=True, help='Path to genes.genes.txt file')
    parser.add_argument('--barcodes', required=True, help='Path to genes.barcodes.txt file')
    parser.add_argument('--output', required=True, help='Path to output h5ad file')

    args = parser.parse_args()

    # Read mtx file
    matrix = io.mmread(args.matrix).tocsr()

    # Read genes and barcodes
    genes = pd.read_csv(args.genes, header=None, names=["gene_id"])
    barcodes = pd.read_csv(args.barcodes, header=None, names=["barcode"])

    # Check dimensions
    if matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(genes):
        pass  # Correct orientation
    elif matrix.shape[0] == len(genes) and matrix.shape[1] == len(barcodes):
        matrix = matrix.transpose()  # Transpose if flipped
    else:
        raise ValueError(
            f"Dimension mismatch: matrix={matrix.shape}, barcodes={len(barcodes)}, genes={len(genes)}"
        )

    # Create AnnData object
    adata = sc.AnnData(X=matrix)
    adata.obs_names = barcodes["barcode"].values
    adata.var_names = genes["gene_id"].values

    # Save to h5ad
    adata.write_h5ad(args.output)

if __name__ == '__main__':
    main()
