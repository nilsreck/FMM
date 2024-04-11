# Reproducibility Study
This repository contains the customized code of our reproducibility study of the Bioinformatics paper
"A functional analysis of omic network embedding spaces reveals key altered functions in cancer" by Doria-Belenguer et al. The original
paper can be found [here](https://academic.oup.com/bioinformatics/article/39/5/btad281/7135836?login=false).

## Data Availability 
The following files: `Human_Biogrid_Adj_PPI_.npy`, `Matrix_Genes_GO_BP_PPI.csv`, `Semantic_Human.csv`, `enriched_in_cancer_gobp_terms_cosmic.txt`, 
corresponding to the human PPI network, the GO BP annotations, and the semantic similarity of the GO BP terms, respectively, used in this project 
can be downloaded [here](https://drive.google.com/drive/folders/15byedX0FchunBx5Y8Y4lztrZvODrPjO1?usp=drive_link).
Please donwload these files first and move them to the `./Cleaned_version/Data` folder.

## Code Execution
The file `Cancer_Functional_Organization.py` contains the workflow and is the entry point to the analysis.

## Computational Requirements
Each step of the workflow has different computational requirements, listed in the table below:

| Step                    | CPUs | Memory | Time     |
|-------------------------|------|--------|----------|
| NMTF                    | 64   | 32GB   | 12 hours |
| Functional Organisation | 1    | 1GB    | <1 hour  |
| Dimensionality          | 4    | 64GB   | 4 hours  |
| Movements               | 1    | 6GB    | <1 hour  |
| Literature Search       | 1    | 1GB    | 12 hours |

