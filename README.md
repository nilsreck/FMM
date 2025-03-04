# Reproducibility Study of "A functional analysis of omic network embedding spaces reveals key altered functions in cancer" by Doria-Belenguer et al.

**Authors:** Nils Reck, Henning Mitzinger, Arne Kugel <br>
**Affiliation:** Chair of Algorithmic Bioinformatics, Heinrich Heine University Düsseldorf

## Summary of the Reproducibility Study
This repository contains a reproducibility study of [*"A functional analysis of omic network embedding spaces reveals key altered functions in cancer"*](https://academic.oup.com/bioinformatics/article/39/5/btad281/7135836?login=false) by Doria-Belenguer et al., published in Bioinformatics 39.5 (2023).
The original paper introduces a Functional Mapping Matrix (FMM) to embed protein-protein interaction (PPI) networks while incorporating functional annotations. 
However, significant challenges were encountered during our reproducibility efforts.

## Key Findings
- **Reproducibility Issues:** The provided code and datasets were insufficient for reproducing the results as described in the paper, requiring major adjustments.
- **Codebase Limitations:** The software is highly dataset-specific, making it difficult to apply to different biological domains.
- **Contradictions in Results:** Contrary to the original claims, the method did not significantly enrich lung cancer-related annotations.
- **Dimensionality Selection Unclear:** While we could reproduce the relative squared errors for different embedding dimensions, the optimal dimensionality remains uncertain.
- **Software Improvements:** To enhance usability and reproducibility, we implemented a Snakemake workflow that improves scalability, resource efficiency, and ease of use.

For a detailed analysis of our methodology, challenges, and findings, read our full [ReScience report](insert link here).

## Improvements in This Repo
- Fix software issues (incompatible environments, missing data, hardcoded hyperparameters).
- Introduce a [Snakemake](https://f1000research.com/articles/10-33/v1) workflow for modular, scalable execution of the workflow.
- Calculate resource requirements for every step of the pipeline.

## Full Reproducibility Report

## How to Reproduce Our Results

### Data Availability 
The files `Human_Biogrid_Adj_PPI_.npy`, `_Matrix_Genes_GO_BP_PPI.csv` and  `enriched_in_cancer_gobp_terms_cosmic.txt` correspond to the human PPI network, the GO BP annotations, and the semantic similarity of the GO BP terms, respectively, used in this project and can be downloaded [here](https://drive.google.com/drive/folders/1SlZ1QixgQu0DoCJabR_cMzjiJv6aM7Pr).
For reasons of reproducibility, we have saved the data from the original paper in a separate Google drive folder. 
Please donwload these files first and move them to the `Data` folder.

### Parameters
The pipeline takes a few parameters, most of which are set to replicate the results Doria Belenguer et al. These Parameters can be found in a configuration file at `workflow/config.yaml`.
Because this pipeline was created in the scope of a reproducibility study, some parameters will not function, as the original code was erroneous and we only fixed those parameter choices that copy the results by the original authors.
Note that the necessary parameter `NCBI_account` is empty. To successfully execute the script, this has to be filled with an email linked to an NCBI account for the automated literature search. \\

### Dependencies
This pipeline requires four different environments. Yaml files can be found in the `envs` folder.
However, only the `snakemake.yaml` needs to be created manually to execute this pipeline. 
The other two will be created automatically during runtime upon the first execution.

### Computational Requirements
This pipeline was designed for cluster execution.
To configure cluster-specific job allocation methods, open the config file found at `workflow/profiles/default/config.yaml`.
See the [snakemake cluster execution](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html) and [snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) for further details on how to set this up.
Each step of the workflow has different computational requirements that can be found in `workflow/rules/steps.smk`.
The range of required resources reaches from 1 to 8 parallel processes with up to 32 GB of memory.
Additionally, the results will take up about 8 additional GB of disc space, and one process requires access to the NCBI website using [this URL](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi). 

### Execution
To exucute the pipeline, use `conda activate snakemake` (the environment found in `snakemake.yaml`) and enter the workflow directory. 
The pipeline can then simply be activated by typing `snakemake` to the command line.

### Results
The results will be copied into the `Results` directory. 
Those consist of the prediction table and results from different statistical analyses:
+ The prediction tables $`\text{(Predictions\_Rank\_{tissue}.csv)}`$ are sorted by predictive confidence and entail the amount of literature validated occurences of the corresponding gene. These are produced by rule "eval_predictions".
+ $`\text{Fold\_Rank\_Table\_Dos.txt}`$ contains the results of a statistical test using the hypergeometric distribution where the H0-hypothesis is "those genes, that change their distance to moving/stable annotations the most/least, are not enriched/depleted in literature validated cancer-genes". The resulting table contains the number of such genes and the p-value in brackets. This is produced by rule "eval_predictions". 
+ $`\text{Network\_Statistics.txt}`$ contains details for each gene network. This is produced by rule "calculate_networks".
+ $`\text{Venn\_Diagrams\_Networks.svg}`$ shows the amount of overlap between the different gene networks for each tissue. This is produced by rule "calculate_networks".
+ $`\text{Relative\_Error\_Leaf.svg}`$ shows the effects of increasing dimensionality for the embedding spaces on the FMM Matrices. The original authors use this to deduce the optimal dimensionality where the relative error stops declining. This is produced by rule "eval_optimal_dimensionality".
+ $`\text{cancer\_enrichments\_2std.svg}`$ shows cancer-relatedness of annotations grouped by their movement between the embedding spaces for cancer and control tissues. This shows that this movement is indeed correlated to annotations adopting cancer functions. This is produced by rule "eval_movements".
+ $`\text{movement\_evaluation.txt}`$ shows the results of two statistical evaluations of movement. The first is a Mann–Whitney U test using movement of known cancer-related annotations (calculated by the original authors) and annotations that are not. The file used for these is "enriched_in_cancer_gobp_terms_cosmic.txt". The second Test is using a Hypergeometric distribution to do the same evaluation; H0-hypothesis being "High annotation movement is un-correlated to cancer-relatedness". This file is produced by rule "eval_movemants".
+ $`\text{Functional\_Organization\_Leaf\_Common\_Set.svg}`$ shows semantic similarity between annotations closely together (Similar) and annotations far apart (Dissimilar) within the embedding spaces. This portrays that biological functionality is captured by distances in this embedding space. This is produced by rule "eval_functional_organization_1".
+ $`\text{Organization\_Common\_Space.txt}`$ contains a table with the results of another Mann–Whitney U test. The two groups here are distances for each gene to belonging annotations and non-belonging annotations. This is produced by rule "eval_functional_organization_2".


## Brief Rule Descriptions

### `prepare_resources`
**Creates auxiliary files**

---

### `calculate_common_gene`
**Filters genes that are common in all input files**

#### Steps:
- **Load annotations** → `_Matrix_Genes_GO_BP_PPI.csv`
- **Load genes** → `{Group}_{tissue}_Genes.csv`
- **Intersect the two files** for cancer and control  
  **→ Common genes** → `Common_Set_{tissue}_Leaf.csv`

---

### `calculate_networks`
**Generate the networks and their statistics**

#### Steps:
- **Load genes** → `{group}_{tissue}_Genes.csv`
- **Load PPI** → `Human_Biogrid_Genes_PPI_`
- **Load adjacencies** → `Human_Biogrid_Adj_PPI_.npy`
- **Generate network** → `{group}_{tissue}_PPI.npy`
- **Generate Venn Diagram** → `Venn_Diagrams_Networks.svg`
- **Compute Network Statistics** → `Network_Statistics.txt`

---

### `calculate_PPMI`
**Generate the PPMI matrix**

#### Steps:
- **Load network** → `{group}_{tissue}_PPI.npy`
- **Perform a deep-walk**
- **Generate PPMI** → `{group}_PPMI_{tissue}.npy`

---

### `calculate_embeddings`
**Generate the embedding coordinates for the genes**

#### Steps:
- **Load network** → `{group}_{tissue}_PPI.npy`
- **Load PPMI** → `{group}_PPMI_{tissue}.npy`
- **Perform NMTF**:
  - **Embedding space (G):** `Embeddings/G_Matrix_{dim}_PPI_{tissue}_{matrix}_{group}.npy`
  - **Gene coordinates (P):** `Embeddings/P_Matrix_{dim}_PPI_{tissue}_{matrix}_{group}.npy`
  - **Gene coordinates (S):** `Embeddings/U_Matrix_{dim}_PPI_{tissue}_{matrix}_{group}.npy`

---

### `calculate_annotation_coordinates`
**Calculate the embedding coordinates for annotations (U-matrix)**

#### Steps:
- **Load annotations** → `_Matrix_Genes_GO_BP_PPI.csv`
- **Load genes** → `{group}_{tissue}_Genes.csv`
- **Intersect annotations with genes and clean data**
- **Load embedding space (G)** → `_G_Matrix_{dim}_PPI_{tissue}_PPMI_{group}.npy`
- **Generate annotation coordinates (U)** → `_GO_Embeddings_Leaf_PPI_{tissue}_{dim}_{matrix}_{group}.csv`

---

### `calculate_FMMs`
**Generate the FMMs**

#### Steps:
- **Load annotation coordinates (U)** → `_GO_Embeddings_Leaf_PPI_{tissue}_{dim}_{matrix}_{group}.csv`
- **Calculate cosine annotation distances**
- **Generate FMM** → `Cosine_{group}_{tissue}_{dim}_PPMI_Leaf.csv`

---

### `calculate_movements`
**Perform enrichment analysis (self-made cancer relatedness measure)**

#### Steps:
- **Load cosine annotation distances (FMM)** → `Cosine_{group}_{tissue}_{dim}_PPMI_Leaf.csv`
- **Load common set** → `Common_Set_{tissue}_Leaf.csv`
- **Intersect the two files**
- **Calculate movement between cancer and control**
- **Generate annotation movement** → `Rank_movement_{tissue}_PPMI_Leaf.csv`

---

### `literature_search`
**Perform literature validation**

#### Steps:
- **Load top moving annotations** → `top_100_moving_{tissue}.csv`
- **Query top moving to NCBI** → `Top_moving_{tissue}_Table.csv`
- **Load gene descriptions (arbitrary genes)** → `Transformed_Common_{tissue}.csv`
- **Query gene descriptions to NCBI**
- **Generate gene literature results** → `Cancer_Count_{tissue}.csv`

---

### `calculate_annotation_gene_distances`
**Calculate the distance between genes and annotations in the embedding space**

#### Steps:
- **Load gene coordinates (P, S)**
- **Compute gene coordinates by multiplying P and S**
- **Load annotation coordinates (U)**
- **Load annotation movement**
- **Subset the top 100 moving annotations**
- **Compute distances between top 100 annotations and genes**
- **Generate annotation gene distances** → `{group}_Gene_GO_dist_{tissue}.csv`

---

### `eval_functional_organization_2`
**Test whether embedding spaces capture annotation belongingness for genes**

#### Steps:
- **Load annotation gene distances** → `{group}_Gene_GO_dist_{tissue}.csv`
- **Split distances into belonging and non-belonging annotations**
- **Perform Mann-Whitney U test**
- **Generate evaluation results** → `Organization_Common_Space.txt`

---

### `calculate_optimal_dimensionality`
**Calculate relative error between spaces and deduce the optimal dimensionality**

#### Steps:
- **Load cosine annotation distances (FMM)**
- **Load annotations** → `gene2go_Human_PPIGO_Specific_BP.json`
- **Intersect the two files**
- **Calculate relative error between embedding spaces**
- **Generate relative error results** → `FMM/Relative_{group}_{tissue}_PPMI_Leaf.txt`

---

### `eval_optimal_dimensionality`
**Plot relative error between embedding spaces of increasing dimensionality**

#### Steps:
- **Load relative error** → `FMM/Relative_{Group}_{tissue}_PPMI_Leaf.txt`
- **Generate error plot** → `FMM/Relative_Error_Leaf.svg`

---

### `eval_movements`
**Statistical tests to evaluate correlation between annotation movement and cancer**

#### Steps:
- **Load cancer annotations** → `enriched_in_cancer_gobp_terms_cosmic.txt`
- **Load annotation movement**
- **Perform Mann-Whitney U and Hypergeometric tests**
- **Generate movement evaluation** → `movement_evaluation.txt`
- **Generate movement evaluation plot** → `cancer_enrichments_2std.svg`

---

### `calculate_semantic_similarity`
**Assess functional organization of genes and functions in the embedding space**

#### Steps:
- **Load cosine annotation distances (FMM)**
- **Load common genes**
- **Intersect the two**
- **Subset most similar/dissimilar annotation pairs**
- **Load annotation descriptions** → `go-basic.obo`
- **Calculate Lin similarity & Jaccard distance**
- **Generate similarities** → `Similarity_{tissue}_Common_Set_500.json`

---

### `eval_functional_organization_1`
**Plot how well the embedding space captures semantic similarity between annotations**

#### Steps:
- **Load similarities** → `Similarity_{tissue}_Common_Set_500.json`
- **Generate similarity plot** → `Functional_Organization_Leaf_Common_Set.svg`

---

### `eval_predictions`
**Perform gene predictions and validate them with enrichment analysis**

#### Steps:
- **Load annotation movement**
- **Subset most moving annotations (2*std)**
- **Load annotation gene distances**
- **Compute change of distance for genes**
- **Normalize distribution & subset top 5%**
- **Load gene literature**
- **Perform enrichment analysis**
- **Generate predictions** → `Predictions_Rank_{tissue}.csv`

---

### `gather_results`
**Copy results into result folder**




