## seurat qc-1

**Commands:**

Generate the pipeline.yml file

    scflow seurat qc-1 config

Run the pipeline

    scflow seurat qc-1 make full -v5

The pipeline pipeline_qc-1.py runs an R markdown file called QC.Rmd to assess the statistics of single-cell experiments run using pipeline_kb.py,  pipeline_quantnuclei.py or pipeline_quantcells.py.

### QC.Rmd

**Inputs:**

Count matrix generated from quantnuclei pipeline  

**Steps:**
1. Read in the count matrix
2. Create Seurat object, Normalize data, scale data and find variable features
3. Create some ggplot themes
4. Add columns to metadata: nUMI, nGene, log10GenesPerUMI
5. Identify mitochondrial genes and add to metadata
6. Map ensembl symbols to hgnc symbols using biomart
7. Filter out genes and cells with low number of counts
8. Create SingleCellExperiment object and save RDS files
9. Plot and save QC metrics
	- Cell counts per sample
	- UMI counts per cell
	- Genes detected per cell
	- UMIs vs genes detected
	- Mitochondrial counts ratio
	- Novelty

**Outputs:**

QC.Rmd knitted to html
SingleCellExperiment and Seurat Object RDS objects saved in RDS_objects.dir
QC plots saved as .eps files in QC_Figures.dir
