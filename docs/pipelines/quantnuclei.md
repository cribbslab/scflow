## scflow main quantnuclei overview

**Commands:**

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
