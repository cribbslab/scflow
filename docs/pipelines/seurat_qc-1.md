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
2. Create Seurat object
3. Generate additional QC metrics: percent mitochondrial genes, log10GenesPerUMI
4. Add gene symbols to the meta.features of the RNA assay
5. Save Seurat Object, Create SingleCellExperiment object and save RDS files
6. Plot and save QC metrics
	- Cell counts per sample
	- UMI counts per cell
	- Genes detected per cell
	- Mitochondrial counts ratio
	- UMIs vs genes detected
	- Novelty

**Outputs:**

QC.Rmd knitted to html  
SingleCellExperiment and Seurat Object RDS objects saved in RDS_objects.dir  
QC plots saved as .png files in QC_Figures.dir  
