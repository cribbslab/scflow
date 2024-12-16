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

### Conda installation - **in progress**

The preferred method for installation is through conda/mamba as a separate environment:

    (Optional) # Create a new folder in your desired location and cd into it:
    mkdir scflow_install && cd scflow_install
    
    # Install the scflow package:
    git clone git@github.com:cribbslab/scflow.git
    mamba env create -f scflow/conda/environments/scflow.yml
    conda activate scflow
    cd scflow
    python setup.py develop

    # Install a modified version of kb-tools from Adam's cloned repo:
    git clone git@github.com:Acribbs/kb_python.git
    cd kb_python
    python setup.py develop

    # Check that the installation has worked:
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

- [ ] [Overview of the quantnuclei pipeline](docs/pipelines/quantnuclei.md)

## scflow main quantcells

- [ ] [Introduction to the quantcells pipeline](docs/pipelines/Singlecell.rst)

## seurat qc-1  

- [ ] [Overview of the seurat qc-1 pipeline](docs/pipelines/seurat_qc-1.md)

## seurat filter-2

- [ ] [Overview of the seurat filter-2 pipeline](docs/pipelines/seurat_filter-2.md)

## seurat cluster-3

- [ ] [Overview of the seurat cluster-3 pipeline](docs/pipelines/seurat_cluster-3.md)

## seurat doublet-4

- [ ] [Overview of the seurat doublet-4 pipeline](docs/pipelines/seurat_doublet-4.md)

## seurat integration-5

- [ ] [Overview of the seurat doublet-4 pipeline](docs/pipelines/seurat_doublet-4.md)

# Project Info

- [ ] [Contributors](docs/project_info/Contributing.rst)

- [ ] [License](docs/project_info/License.rst)

- [ ] [How to contribute](docs/project_info/how_to_contribute.rst)
