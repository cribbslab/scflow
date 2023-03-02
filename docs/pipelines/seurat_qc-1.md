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
Patient metadata as tab-delimited file    


**Steps:**
1. Read in the count matrix
2. Create Seurat object
3. Generate additional QC metrics: percent mitochondrial genes, log10GenesPerUMI
4. Save the mapping reference for ensemble gene names to gene IDs
5. Add gene symbols to the meta.features of the RNA assay
6. Add patient metadata
7. Convert to SingleCellExperiment 
	- Identify and remove empty droplets (NB empty droplets must be removed in order to perform doublet detection).  
	- Identify ambient RNA
	- Identify doublets
8. Save Seurat Objects and SingleCellExperiment objects as RDS files (both complete and with emptry droplets filtered out)
9. Plot and save QC metrics
	- Cell counts per sample
	- UMI counts per cell
	- Genes detected per cell
	- Mitochondrial counts ratio
	- UMIs vs genes detected
	- Novelty

**Outputs:**

QC.Rmd knitted to html  

SingleCellExperiment and Seurat Object RDS objects saved in RDS_objects.dir/unfiltered    
QC plots saved as .png files in QC_Figures.dir  
mapping.txt saved in Files.dir  
