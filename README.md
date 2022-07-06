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

The respository can also be installed manually, but dependancies will need to be installed seperately::

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




