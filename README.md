# scflow

<p align="left">
	<a href='https://single-cell.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/single-cell/badge/?version=latest' alt='Documentation Status' /></a>
	<a href="https://github.com/Acribbs/scflow/actions/workflows/aattggcc_python.yml/badge.svg", alt="Actions">
		<img src="https://github.com/Acribbs/scflow/actions/workflows/aattggcc_python.yml/badge.svg" /></a>
	<a href="https://twitter.com/CribbsP?lang=en", alt="Twitter followers">
		<img src="https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter" /></a>
</p>


This repository contains a collection of pipelines that aid the analysis of single cell sequencing experiments. Currently there is one pipeline implimented that allows the analysis of drop-seq and 10X sequencing analysis. Current pipelines in development: 1) pseudoalignment scpipeline 2) velocyto pipeline 2) kallisto bustools pipeline.

## Installation

### pip install

You can install scflow using pip, this will only install the package without any dependancies, which will have to be installed seperately.::

	pip install scflow

### Conda installation - **in progress**

The preferred method for installation is through conda. Currently this installation is still in working progress. Preferably the
installation should be in a seperate environment::

    conda create -n scflow -c cgat scflow
    conda activate scflow
    scflow --help

### Manual installation

The repository can also be installed manually, but dependencies will need to be installed seperately::

    python setup.py install
    scflow --help

## Usage

Run the ``scflow --help`` command view the help documentation for how to run the single-cell repository.

To run the main single_cell droplet based pipeline run first generate a configuration file::

    scflow singlecell config

Then run the pipeline::

    scflow singlecell make full -v5

Then to run the report::

    scflow singlecell make build_report

## Documentation

Further help that introduces single-cell and provides a tutorial of how to run example
code can be found at [read the docs](http://single-cell.readthedocs.io/)

# Pipelines overview

## scflow main quantnuclei


- [ ] [Introduction to the quantnuclei pipeline](docs/pipelines/Singlenuclei.rst)

**Commands**

Generate the pipeline.yml file

    scflow main quantnuclei config

Run the pipeline

    scflow main quantnuclei make full -v5


**Inputs:**  

Genome reference files  
Fastq files from 10X experiment  
pipeline.yml

**Steps:**
1. Builds kallisto index using kb ref
2. Performs read quality steps with fastqc
3. Performs pseudoalignment using kb count
4. Merges the spliced and unspliced matrix using custom python script

**Outputs**  

Kallisto index  
Fastqc html files  
Count matrix  

## seurat qc-1  

**Commands**

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


## seurat filter-2

**Commands**

Generate the pipeline.yml file

    scflow seurat filter-2 config

Run the pipeline

    scflow seurat filter-2 make full -v5

The first step in pipeline_filter-2.py runs an R markdown file called Filter.Rmd

### Filter.Rmd


**Inputs:**  

Seurat Object generated in qc-1

**Steps:**
1. Read in Seurat Object from qc-1
2. Set publication themes for ggplot
3. Filtering - more stringent than in qc-1. nGene 300-6000, nUMI > 100 and mitoRatio <0.1 are defaults, these can be updated in the pipeline.yml file.
4. Save Seurat Object as RDS file.
5. Make QC plots post-filtering as in qc-1 pipeline.


**Outputs:**  

Filter.Rmd knitted to html
Filtered SingleCellExperiment and Seurat Object RDS objects saved in RDS_objects.dir
QC plots saved as .eps files in Filtered_Figures.dir

At this point it is good to compare the outputs of QC.html and Filter.html to check that the filtering steps have tidied up the dataset as you would expect (removed outliers and so on).

pipeline_filter-2.py then converts the Seurat Objects (.rds files) into AnnData files stored in Andata.dir. This is another format of annotated data matrices that is used in python scripts.

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
