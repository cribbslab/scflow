## seurat doublet-4

pipeline_doublet-4.py overview:  
- Make directory called Doublet_Figures.dir
- Run the R script doublet_finder.R
- Run the Rmd Doublets.Rmd to generate Doublets.html

**Commands:**
> scflow seurat doublet-4 config  
> nohup scflow seurat doublet-4 make full -v5


### doublet_finder.R

**Inputs:**  
Clustered filtered Seurat object created by seurat_cluster.R

**Options:**
Defined in pipeline.yml

| Option | Description
|---------|-------------|
|	-i -input |	input files
|	-s --sample	|sample name

**Steps:**
- Read in the seurat object
- Convert to SingleCellExperiment object
- Call doublets with scDblfinder
- Save a table describing the number of doublets identified
- Save RDS object

**Outputs:**
- Table (.csv) describing number of doublets identified
- RDS object including doublet information

### Doublets.Rmd

**Inputs:**  
Clustered filtered doublets Seurat object created by doublet_finder.R

**Steps:**
- Read in the seurat object
- Normalise the data
- Find Variable features
- Scale the data
- Perform PCA
- Find Neighbours
- Find Clusters
- Generat UMAP
- Plot UMAP displaying doublets and singlets

**Outputs:**
- UMAP (.eps) displaying doublets and singlets
- Doublets.html
