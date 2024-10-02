import os
import sys
import json
import subprocess
import numpy as np
import pandas as pd
from Bio import Entrez
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

# 9.4 Literature validation
# 9.5 Define cancer-related genes (based on the literature)

def Query_GO_terms(Normal_tissue_list, data_dir, annotation, NCBI_account):
    """
    To run this function is mandatory to manually obtain the top 100 moving file is needed. This file can be recomputed manually
    for new runs. The file contains the description of the GO term and its total movement. Since this ranking can variate depending
    on the annotation file that is used (and also on the embeddings) it should be recomputed with the new ranking if needed.
    """

    movement_path = f"{data_dir}/"
    path_rank = f"{data_dir}/"

    for tissue in Normal_tissue_list:
        file = pd.read_csv(f"{movement_path}top_100_moving_{tissue}.csv")

        rank = pd.read_csv(f"{path_rank}Rank_movement_{tissue}_PPMI_{annotation}.csv")
        moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2 * np.std(rank["0"]))][
            "Unnamed: 0"
        ]

        file = file.head(len(moving_GO))

        cancer_q = f"{tissue} cancer"

        counts = []

        for row in range(len(file)):
            query = (
                f'({file.iloc[row]["GO name"]}[Text Word]) AND ({cancer_q}[Text Word])'
            )
            query_file = search(query, NCBI_account)
            citation_list = query_file["IdList"]

            counts.append(len(citation_list))

        Result = pd.DataFrame()
        Result["Annotation"] = list(file["GO name"])
        Result["norm"] = list(file["norm"])
        Result["Cancer Related"] = list(file["Cancer Related"])
        Result["Bibliography"] = counts

        Result.to_csv(f"{movement_path}Top_moving_{tissue}_Table.csv")


def search(query, NCBI_account):
    Entrez.email = NCBI_account
    print(query, flush=True)
    handle = Entrez.esearch(
        db="pubmed", sort="relevance", retmax="100000", retmode="xml", term=query
    )
    results = Entrez.read(handle)
    return results


def Query_Common_gene_set(Normal_tissue_list, data_dir, NCBI_account):
    """
    To run this function the file Transformed_Common_{tissue}.csv is needed. This file corresponds to a query in
    bioMart. For each gene, we retrieve the description of the gene and different name codes. The file is attached
    but it can be manually recomputed. By using this file, the function performs a literature search for assessing if a
    gene is related to a particular cancer type.
    """

    # Paths:

    gene_GO_path = f"{data_dir}/"

    # For each tissue:

    for tissue in Normal_tissue_list:
        common_genes = pd.read_csv(f"{gene_GO_path}Transformed_Common_{tissue}.csv")
        common_genes = common_genes.drop_duplicates(subset=["initial_alias"])
        common_genes = common_genes[["initial_alias", "name"]]
        common_genes = common_genes.dropna()
        common_genes = common_genes.reset_index(drop=True)

        cancer_q = f"{tissue} cancer"
        counts = []
        
        for row in range(len(common_genes)):
            query = f'({common_genes.iloc[row]["name"]}[Text Word]) AND ({cancer_q}[Text Word])'
            query_file = search(query, NCBI_account)
            citation_list = query_file["IdList"]

            counts.append(len(citation_list))
        

        # Prepare the data frame:

        common_genes["Counts"] = counts

        common_genes.to_csv(
            f"{gene_GO_path}Cancer_Count_{tissue}.csv", header=True, index=True
        )

data_path = sys.argv[1]
annotation = sys.argv[2]
NCBI_account = sys.argv[3]
log = sys.argv[4]
Normal_Tissue_list = sys.argv[5:]

f = open(log, "a")
sys.stdout = f

Query_GO_terms(Normal_Tissue_list, data_path, annotation, NCBI_account)

Query_Common_gene_set(Normal_Tissue_list, data_path, NCBI_account) 

f.close()

