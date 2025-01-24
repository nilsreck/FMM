# Reproducibility Study
This repository contains the Snakemake project of our reproducibility study of the Bioinformatics paper
"A functional analysis of omic network embedding spaces reveals key altered functions in cancer" by Doria-Belenguer et al. The original
paper can be found [here](https://academic.oup.com/bioinformatics/article/39/5/btad281/7135836?login=false). Snakemake is a workflow tool developed by [Felix Mölder et al.](https://f1000research.com/articles/10-33/v1). 

## Data Availability 
The following files: `Human_Biogrid_Adj_PPI_.npy`, `_Matrix_Genes_GO_BP_PPI.csv`, `enriched_in_cancer_gobp_terms_cosmic.txt`, 
corresponding to the human PPI network, the GO BP annotations, and the semantic similarity of the GO BP terms, respectively, used in this project 
can be downloaded [here](https://drive.google.com/drive/folders/1SlZ1QixgQu0DoCJabR_cMzjiJv6aM7Pr).
For reasons of reproducibility, we have saved the data from the original paper in a separate Google drive folder. 
Please donwload these files first and move them to the `Data` folder.

## Parameters
The pipeline takes a few parameters, most of which are set to replicate the results Doria Belenguer et al. . These Parameters can be found in a configuration file at "workflow/config.yaml" . Because this pipeline was created in the scope of a reproducibility study, some parameters will not function, as the original code was erroneous and we only fixed those parameter choices that copy the results by the original authors.
Note that the necessary parameter "NCBI_account" is empty. To successfully execute the script, this has to be filled with an email linked to an NCBI account. \\

## Dependencies
This pipeline requires four different environments. Yaml files can be found in the "envs" folder. However, only the "snakemake.yaml" needs to be created manually to execute this pipeline. The other two will be created automatically during runtime upon the first execution.

## Computational Requirements
This pipeline was designed for cluster execution. To configure cluster specific job allocation methods, open the config file found at "workflow/profiles/default/config.yaml". See the [snakemake cluster execution](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html) and [snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) for further details on how to set this up. Each step of the workflow has different computational requirements that can be found in workflow/rules/steps.smk. The range of required resources reaches from 1 to 8 parallel processes with up to 32 GB of memory. Additionally, the results will take up about 8 additional GB of disc space, and one process requires access to the NCBI website using this URL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi". 

## Execution
To exucute the pipeline, use conda activate snakemake (the environment found in snakemake.yaml) and enter the workflow directory. The pipeline can then simply be activated by typing "snakemake" to the command line.