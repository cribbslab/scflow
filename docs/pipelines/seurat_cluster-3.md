## seurat cluster-3

pipeline_cluster-3.py overview:
- runs an R script called seurat_cluster.R.
- runs a second script called find_markers.R is run to find markers for each cluster.
- runs an R markdown called Cluster.Rmd to visualise clustering and dimensional reduction

**Pipeline commands**

Generate the .yml file
> scflow seurat cluster-3 config

Run the pipeline
> nohup scflow seurat cluster-3 make full -v5

### seurat_cluster.R

**Inputs:**
Seurat objects filtered and unfiltered RDS files created from qc-1 and filter-2 pipelines.

**Options**
Defined in pipeline.yml

| Option | Description | Default |
|---------|-------------|--------|
|	-i -input |	input files |
|	-s --sample	|sample name |
|	-v --variableFeatures	|Number of variable features | 2000 |
|	--reddim	|What dimensionality reduction to use. PCA, harmony, zinbwave, CCA |pca|
|	-d --numdim	|How many dimensions to use. Default to all calculated ones, but can be user defined according to the Elbow Plot and Jack Straw plots|
|	--resolution	|Resolution for finding clusters. Sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. Between 0.4-1.2 good for ~3K cells. |0.5|
	|-m --metadata	|The input metadata file from doublet pipeline|

**Steps:**
- Define some arguments that can be included when the Rscript command is called, using the package "optparse".
- Read in some metadata about doublets if this has been performed (not yet, maybe we come back to it?)
- Read in the RDS files (Seurat objects, filtered and unfiltered).
- Add doublet metadata to the Seurat object if relevant.

- Perform usual steps in Seurat pipeline (see https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
	- Normalise the data
	- Find variable features
	- list top 10 variable genes (top10)
	- plot variable features (labelled_variable_feature_plot)
	- Scale the data
	- Perform PCA
	- Perform clustering (UMAP and tSNE)

**Outputs:**
RDS files filtered and clustered in RDS_objects.dir.


### find_markers.R

**Inputs:**
Clustered filtered Seurat object created by seurat_cluster.R

**Options**
Defined in pipeline.yml

| Option | Description | Default |
|---------|-------------|--------|
|-i --input	|input rds path of filtered clustered Seurat object
|-s --sample|	sample name
|-m --minPct|	Genes only tested if found in minimum percentage of cells in either population. | 0.1
|-l --logfc	|"Limit testing to genes which have (on average) a log fold change greater than this threshold |0.25
|-t --testuse|	Test to use. Options: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2. |wilco
|-c --maxClusters	|If you want to set a maximum number of clusters to find markers for. 0 = Do it for all clusters.| 0

**Steps:**
- Define some arguments that can be included when the Rscript command is called, using the package "optparse".
- Read in the RDS files of the clustered filtered object that was created in seurat_cluster.R
- Find the total number of clusters and choose whether to use this number or assign it manually with -c parameter.

- For each sample
- Generate some statistics for each gene in each cluster
	- First define which cells are in the active cluster
	- Calculate the mean and experimental mean of each gene expression within the active cluster and all the other clusters, then compare the two means.
	- Identify markers for each cluster with FindMarkers
	- Generate a summary table of markers (combined)
	- Generate a summary table of markers and means (markers_filter_stats_combined)
	- Add the BH corrected p value to the summary tables (combined and markers_filter_stats_combined)
	- Select appropriate columns and order by padj/log2fc
	- Save as csv and tsv files


**Outputs:**
Summary tables of cluster markers

### Cluster.Rmd

**Inputs:**
Filtered clustered Seurat Objects created by seurat_cluster.R

**Steps:**
Read in seurat objects
Plot and save
- top 2000 variable features, labelled with EnsemblID
	(could be good to include gene name?)
- JackStraw plot and ElbowPlot
- PCA showing cluster membership
- PCA loading for PCs 1 and 2
- PC heatmaps for PCs 1-9
- tSNE
- UMAP

**Outputs:**

Plots in Clustering_Figures.dir for each sample
- Variable features
- JackStraw plot
- Elbow plot
- PCA showing cluster membership
- PCA loading for PCs 1 and 2
- PC heatmaps for PCs 1-9
- tSNE
- UMAP

Cluster.html markdown
