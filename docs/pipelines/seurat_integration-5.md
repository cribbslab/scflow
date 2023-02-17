## seurat integration-5

The pipeline follows pipeline quantnuclei, qc-1, filter-2 and optionally cluster-3.

It takes filtered Seurat objects and performs integrations by two methods, Seurat and Harmony.
Each dataset is first processed individually by SCTransform normalization and then perform integration and clustering by Seurat.
The integrated data could be furthered accepted by Harmony with integration/clustering.
The integrations performed by both the methods could be visualized as tSNE and UMAP dimensional reduction.

**Overview**  
Runs the R script "seurat_integrate.R"  
Runs the R script "harmony_integrate.R"  
Runs Rmd file "Integration.Rmd" to visualise the results.

**Commands:**  

Configure the pipeline.yml file

> scflow seurat doublet-4 config  

Run the pipeline
> nohup scflow seurat doublet-4 make full -v5

### seurat_integrate.R

**Inputs:**  
Filtered clustered seurat object.  
(Or would work using filtered seurat object from output of filter-2.

**Steps:**  
- Read in the .yml file
- Extract parameters from the .yml file
- Read in the samples
- Perform integration steps (creates a new seurat object which contains the integrated data):
    - Normalise the data by SCTransform
    - Select variable features that are common across the Seurat Objects
    - Use PrepSCTIntegration to get residuals for all features
    - Identify anchors
    - Integrate the data  
- Perform PCA
- Find Neighbours
- Find Clusters
- Perform UMAP
- Plot and save UMAPs
- Perform tSNE
- Plot and save tSNEs
- Save integrated object as RDS file

**Outputs**
UMAP and tSNE Plots  
Integrated RDS object  
