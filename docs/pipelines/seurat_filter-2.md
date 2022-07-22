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
