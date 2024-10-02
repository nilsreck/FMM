# Reproducibility Study
This repository contains the Snakemake project of our reproducibility study of the Bioinformatics paper
"A functional analysis of omic network embedding spaces reveals key altered functions in cancer" by Doria-Belenguer et al. The original
paper can be found [here](https://academic.oup.com/bioinformatics/article/39/5/btad281/7135836?login=false). Snakemake is a workflow tool developed by [[Felix Mölder et al.](https://f1000research.com/articles/10-33/v1)].

## Data Availability 
The following files: `Human_Biogrid_Adj_PPI_.npy`, `_Matrix_Genes_GO_BP_PPI.csv`, `Semantic_Human.csv`, `enriched_in_cancer_gobp_terms_cosmic.txt`, 
corresponding to the human PPI network, the GO BP annotations, and the semantic similarity of the GO BP terms, respectively, used in this project 
can be downloaded [here](https://drive.google.com/drive/folders/1SlZ1QixgQu0DoCJabR_cMzjiJv6aM7Pr).
For reasons of reproducibility, we have saved the data from the original paper in a separate Google drive folder. 
Please donwload these files first and move them to the `Data` folder.

## Parameters
The pipeline takes a few parameters, most of which are already set to replicate the results Doria Belenguer et al. . These Parameters can be found in a configuration file at workflow/config.yaml .
Note that the necessary parameter "NCBI_account" is empty. To successfully execute the script, This has to be filled with an email linked to an NCBI account.

## Dependencies
This pipeline requires three different environments. Yaml files can be found in the "envs" folder. However, only the "snakemake.yaml" needs to be created manually to execute this pipeline. The other two will be created automatically during runtime upon the first execution.

## Computational Requirements
Each step of the workflow has different computational requirements and can be found in workflow/rules/steps.smk. The range of required resources reaches from 1 to 6 parallel processes with up to 4 GB of memory for each of these 6 threads. Additionally, the results will take up about 8 additional GB of disc space, and one preocess requires access to the NCBI website using this URL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi". 

