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

## Results
The results will be copied into the "Results/" directory. 
Those consist of the prediction table and results from different statistical analyses:
+ The prediction tables $\text{(Predictions_Rank_{tissue}.csv)}$ are sorted by predictive confidence and entail the amount of literature validated occurences of the corresponding gene. These are produced by rule "eval_predictions".
+ $\text{Fold_Rank_Table_Dos.txt}$ contains the results of a statistical test using the hypergeometric distribution where the H0-hypothesis is "those genes, that change their distance to moving/stable annotations the most/least, are not enriched/depleted in literature validated cancer-genes". The resulting table contains the number of such genes and the p-value in brackets. This is produced by rule "eval_predictions". 
+ $\text{Network_Statistics.txt}$ contains details for each gene network. This is produced by rule "calculate_networks".
+ $\text{Venn_Diagrams_Networks.svg}$ shows the amount of overlap between the different gene networks for each tissue. This is produced by rule "calculate_networks".
+ $\text{Relative_Error_Leaf.svg}$ shows the effects of increasing dimensionality for the embedding spaces on the FMM Matrices. The original authors use this to deduce the optimal dimensionality where the relative error stops declining. This is produced by rule "eval_optimal_dimensionality".
+ $\text{cancer_enrichments\_2std.svg}$ shows cancer-relatedness of annotations grouped by their movement between the embedding spaces for cancer and control tissues. This shows that this movement is indeed correlated to annotations adopting cancer functions. This is produced by rule "eval_movements".
+ $\text{movement_evaluation.txt}$ shows the results of two statistical evaluations of movement. The first is a Mann–Whitney U test using movement of known cancer-related annotations (calculated by the original authors) and annotations that are not. The file used for these is "enriched_in_cancer_gobp_terms_cosmic.txt". The second Test is using a Hypergeometric distribution to do the same evaluation; H0-hypothesis being "High annotation movement is un-correlated to cancer-relatedness". This file is produced by rule "eval_movemants".
+ $\text{Functional_Organization_Leaf_Common_Set.svg}$ shows semantic similarity between annotations closely together (Similar) and annotations far apart (Dissimilar) within the embedding spaces. This portrays that biological functionality is captured by distances in this embedding space. This is produced by rule "eval_functional_organization_1".
+ $\text{Organization_Common_Space.txt}$ contains a table with the results of another Mann–Whitney U test. The two groups here are distances for each gene to belonging annotations and non-belonging annotations. This is produced by rule "eval_functional_organization".


## Brief Rule Descriptions

$ \text{prepare\_resources} $ \
\## _creates auxiliary files_

$ \text{calculate\_common\_gene} $ \
\## _filters genes that are common in all input files_ \
\# load annotations -------------------------------------------------	"\_Matrix\_Genes\_GO\_BP\_PPI.csv" \
\# load genes -------------------------------------------------------- "{Group}\_{tissue}\_Genes.csv" \
\# intersect the two files for cancer and control \
\# -> common genes ------------------------------------------------ "Common\_Set\_{tissue}\_Leaf.csv"

$ \text{calculate\_networks} $ \
\## _generate the networks and their statistics_ \
\# load genes -------------------------------------------------------- "{group}\_{tissue}\_Genes.csv" \
\# load PPI ----------------------------------------------------------- "Human\_Biogrid\_Genes\_PPI\_" \
\# load adjacencies -------------------------------------------------- "Human\_Biogrid\_Adj\_PPI\_.npy" \
\# -> network -------------------------------------------------------- "{group}\_{tissue}\_PPI.npy" \
\# -> Venn-Diagram ------------------------------------------------- "Venn\_Diagrams\_Networks.svg" \
\# -> Network-Statistics --------------------------------------------- "Network\_Statistics.txt"

$ \text{calculate\_PPMI} $ \
\## _generate the ppmi matrix_ \
\# load network ------------------------------------------------------ "{group}\_{tissue}\_PPI.npy" \
\# perform a deep-walk \
\# -> PPMI ----------------------------------------------------------- "{group}\_PPMI\_{tissue}.npy"

$ \text{calculate\_embeddings} $ \
\## _generate the embedding coordinates for the genes_ \
\# load network ------------------------------------------------------ "{group}\_{tissue}\_PPI.npy" \
\# load PPMI --------------------------------------------------------- "{group}\_PPMI\_{tissue}.npy" \
\# do the NMTF \
\# -> embedding space (G) ------------------------------------------ "Embeddings/\_G\_Matrix\_{dim}\_PPI\_{tissue}\_{matrix}\_{group}.npy" \
\# -> gene coordinates (P) ------------------------------------------- "Embeddings/\_P\_Matrix\_{dim}\_PPI\_{tissue}\_{matrix}\_{group}.npy" \
\# -> gene coordinates (S) ------------------------------------------- "Embeddings/\_U\_Matrix\_{dim}\_PPI\_{tissue}\_{matrix}\_{group}.npy"

$ \text{calculate\_annotation\_coordinates} $ \
\## _calculate the embedding coordinates for the annotations (U-matrix)_ \
\# load annotations -------------------------------------------------- "\_Matrix\_Genes\_GO\_BP\_PPI.csv" \
\# load genes --------------------------------------------------------- "{group}\_{tissue}\_Genes.csv" \
\# intersect the annotations with genes and clean the data \
\# load embedding space (G) ---------------------------------------- "\_G\_Matrix\_{dim}\_PPI\_{tissue}\_PPMI\_{group}.npy" \
\# -> annotation coordinates (U) ------------------------------------ "\_GO\_Embeddings\_Leaf\_PPI\_{tissue}\_{dim}\_{matrix}\_{group}.csv"

$ \text{calculate\_FMMs} $ \
\## _generate the FMMs_ \
\# load annotation coordinates (U) ---------------------------------- "\_GO\_Embeddings\_Leaf\_PPI\_{tissue}\_{dim}\_{matrix}\_{group}.csv" \
\# calculate cosine annotation distances \
\# -> cosine annotation distances (FMM) ---------------------------- "Cosine\_{group}\_{tissue}\_{dim}\_PPMI\_Leaf.csv" 

$ \text{calculate\_movements} $ \
\## _enrichment analysis (with self made cancer relatedness)_ \
\# load cosine annotation distances (FMM) --------------------------"Cosine\_{group}\_{tissue}\_{dim}\_PPMI\_Leaf.csv" \
\# load common set -------------------------------------------------- "Common\_Set\_{tissue}\_Leaf.csv" \
\# intersect the two files \
\# calculate movement between cancer and control \
\# -> annotation movement ------------------------------------------ "Rank\_movement\_{tissue}\_PPMI\_Leaf.csv"

$ \text{literature\_search} $ \
\## _Literature validation_ \
\# load top moving --------------------------------------------------- "top\_100\_moving\_{tissue}.csv" \
\# query top moving to ncbi  \
\# -> ------------------------------------------------------------------- "Top\_moving\_{tissue}\_Table.csv" \
\# load gene descriptions (arbitrary genes) -------------------------- "Transformed\_Common\_{tissue}.csv" \
\# query gene descriptions to ncbi \
\# -> gene literature -------------------------------------------------- "Cancer\_Count\_{tissue}.csv"

$ \text{calculate\_annotation\_gene\_distances} $ \
\## _calculate the distance between genes and annotations withn the embedding space_ \
\# load gene coordinates (P) ------------------------------------------ "Embeddings/\_P\_Matrix\_{dim}\_PPI\_{tissue}\_{matrix}\_{group}.npy" \
\# load gene coordinates (S) ------------------------------------------ "Embeddings/\_U\_Matrix\_{dim}\_PPI\_{tissue}\_{matrix}\_{group}.npy" \
\# calculate gene coordinates by multiplying P and S \
\# load annotation coordinates (U) ----------------------------------- "\_GO\_Embeddings\_Leaf\_PPI\_{tissue}\_{dim}\_PPMI\_{group}.csv" \
\# load annotation movement ---------------------------------------- "Rank\_movement\_{tissue}\_PPMI\_Leaf.csv" \
\# subset the top 100 moving from the annotation coordinates \
\# calculate distances between top 100 annotation coordinates and gene coordinates \
\# -> annotation gene distances -------------------------------------- "{group}\_Gene\_GO\_dist\_{tissue}.csv"

$ \text{eval\_functional\_organization\_2} $ \
\## _test whether embedding spaces capture annotation beloningness for genes_ \
\# load annotation gene distances ------------------------------------ "{group}\_Gene\_GO\_dist\_{tissue}.csv" \
\# split in distances to belonging annotations and non belonging annotations \
\# mannwhitneyu test: H0 = distance to belonging annotations does not differ for non-belonging annotations \
\# -> organization evaluation ------------------------------------------ "Organization\_Common\_Space.txt"

$ \text{calculate\_optimal\_dimensionality} $ \
\## _Calculate relative error between spaces and deduce the optimal dimensionality_ \
\# load cosine annotation distances (FMM) --------------------------- "FMM/Cosine\_{group}\_{tissue}\_{dim}\_PPMI\_Leaf.csv" \
\# load annotations ----------------------------------------------------- "gene2go\_Human\_PPIGO\_Specific\_BP.json" \
\# intersect the two files \
\# calculate relative error between embedding spaces of different dimensionalities \
\# -> relative error ------------------------------------------------------ "FMM/Relative\_{group}\_{tissue}\_PPMI\_Leaf.txt"

$ \text{eval\_optimal\_dimensionality} $ \
\## _plot the relative error between embedding spaces of increasing dimensionality_ \
\# load relative error --------------------------------------------------- "FMM/Relative\_{Group}\_{tissue}\_PPMI\_Leaf.txt" \
\# -> Error between dimensionalities ---------------------------------- "FMM/Relative\_Error\_Leaf.svg"

$ \text{eval\_movements} $ \
\## _statistical tests to evaluate the corrlatio between annotation movement and cancer using author defined cancer hallmarks_ \
\# load cancer annotations (by the authors) -------------------------- "enriched\_in\_cancer\_gobp\_terms\_cosmic.txt" \
\# load annotation movement ----------------------------------------- "Rank\_movement\_{tissue}\_PPMI\_Leaf.csv" \
\# mannwhitneyu test: H0 = high movement does not correlate with cancer \
\# hypergeom test: H0 = high movement does not correlate with cancer \
\# -> movement evaluation --------------------------------------------"movement\_evaluation.txt" \
\# plot movement evaluation \
\# -> movement evaluation plot -------------------------------------- "cancer\_enrichments\_2std.svg"

$ \text{calculate\_semantic\_similarity} $ \
\## _Assess the functional organization of genes and functions in the embedding space_ \
\# load cosine annotation distances (FMM) -------------------------- "FMM/Cosine\_{Group}\_{tissue}\_{dim}\_PPMI_Leaf.csv" \
\# load common genes ----------------------------------------------- "Common\_Set\_{tissue}\_Leaf.csv" \
\# intersect the two \
\# subset the most similar and dissimilar annotation pairs \
\# load annotation descriptions ----------------------------------------"go-basic.obo" \
\# calulate Lin similarity for these annotation pairs based on their description \
\# calculate Jaccard distance between cancer and control \
\# -> similarities -------------------------------------------------------- "Similarity\_{tissue}\_Common\_Set\_500.json"

$ \text{eval\_functional\_organization\_1} $ \
\## _plot how well the embedding space captures semantic similarity between annotations via distance_ \
\# load similarities ------------------------------------------------------ "Similarity\_{tissue}\_Common\_Set\_500.json" \
\# -> similarity plot ----------------------------------------------------- "Functional\_Organization\_Leaf\_Common\_Set.svg"

$ \text{eval\_predictions} $ \
\## _Do the gene predictions and validate them with an enrichment analyses_ \
\# load annotation movement -------------------------------------------- "Rank\_movement\_{tissue}\_PPMI\_Leaf.csv" \
\# subset the most moving annotations (2*std) \
\# load annotation gene distances --------------------------------------- "{Group}\_Gene\_GO\_dist\_{tissue}.csv" \
\# intersect the two \
\# calculate the change of distance of genes to the most moving annotations between cancer and control \	
\# normalize the distribution \
\# subset the the top 5% \
\# load gene literature ---------------------------------------------------- "Cancer\_Count\_{tissue}.csv" \
\# intersect the two \
\# -> gene predictions ---------------------------------------------------- "Predictions\_Rank\_{tissue}.csv" \
\# hypergeom test: H0 = top 5% are not enriched in literature validated genes \
\# hypergeom test: H0 = bottom 5% are not depleted in literature validated genes \
\# -> prediction evaluation ----------------------------------------------- "Fold\_Rank\_Table\_Dos.txt"

$ \text{gather\_results} $ \
\## _copy results into result folder_


