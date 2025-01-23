import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.spatial import distance


def solve(process):
    Gen_GO_Common_Space(*process)


def Generate_Final_DataFrame(GO_Gene_Matrix, distances):
    # Lists to fill:

    group = []  # Cluster based on the distance to the GO term
    shape = (
        []
    )  # Shape for the plot (if the gene clusteres is annotated with the corresponding GO)
    size = (
        []
    )  # Size for the plot (if the gene clusteres is annotated with the corresponding GO)
    dist = []  # Minimum distance (the minimim distance to a GO term)
    GO_min = []  # The closest GO term in the space

    # Iteration per gene:

    for gene in list(distances.index):
        it_gene = distances.loc[gene]
        it_annotation = GO_Gene_Matrix.loc[str(gene)]

        if min(it_gene) <= 0.5:
            group.append(it_gene.idxmin())

            if it_annotation.loc[it_gene.idxmin()] == 1:
                shape.append("Annotated")
                size.append(75)
                dist.append(min(it_gene))
                GO_min.append(it_gene.idxmin())
            else:
                shape.append("NO Annotated")
                size.append(55)
                dist.append(min(it_gene))
                GO_min.append(it_gene.idxmin())
        else:
            group.append("NO")
            shape.append("NO Annotated")
            size.append(15)
            dist.append(min(it_gene))
            GO_min.append(it_gene.idxmin())

    # Fill the final dataframe:

    Cosine_it = distances.copy()

    Cosine_it["Group"] = group
    Cosine_it["Shape"] = shape
    Cosine_it["size"] = size
    Cosine_it["Min_Dist"] = dist
    Cosine_it["Min_GO"] = GO_min

    # Return the dataFrame:

    return Cosine_it


def Distance_GO_Gene(vector_GO, matrix_genes):
    vector_GO = vector_GO.loc[:, ~vector_GO.columns.duplicated()]

    # Cosine

    result_pd_cos = pd.DataFrame(
        0.0, index=np.arange(len(matrix_genes.index)), columns=list(vector_GO.columns)
    )

    result_pd_cos.columns = vector_GO.columns
    result_pd_cos.index = matrix_genes.index

    for gene in list(matrix_genes.index):
        gene_vector = matrix_genes.loc[gene]

        for GO in list(vector_GO.columns):
            GO_vector = vector_GO[GO]

            result_pd_cos.loc[gene, GO] = distance.cosine(GO_vector, gene_vector)

    return result_pd_cos


def Gen_GO_Common_Space(cancer, tissue, cell, dim, data_dir, annotation):
    """
    This functions calculates the cosine distance between the gene/GO common space
    """

    # Paths:

    network_path = f"{data_dir}/Embeddings/"
    gene_GO_path = f"{data_dir}/"
    moving_path = f"{data_dir}/"
    FMM_path = f"{data_dir}/FMM/"

    # For each combination:

    # for cancer, tissue, cell in zip(Cancer_list, Normal_tissue_list, Cell_type_list):
    #     # Prepare the paths:

    # Genes:

    genes_Control = f"{gene_GO_path}Control_{tissue}_{cell}_Genes.csv"
    genes_Cancer = f"{gene_GO_path}Cancer_{cancer}_Genes.csv"

    # Gene embedded coordinates:

    Gene_embedd_Control = (
        f"{network_path}_P_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy"
    )
    Gene_embedd_Cancer = (
        f"{network_path}_P_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy"
    )

    # Middle factor:

    midle_Control = (
        f"{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy"
    )
    midle_Cancer = (
        f"{network_path}_U_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy"
    )

    Gene_midle_Control_db = pd.DataFrame(np.load(midle_Control, allow_pickle=True))
    Gene_midle_Cancer_db = pd.DataFrame(np.load(midle_Cancer, allow_pickle=True))

    # GO embedded coordinates:

    GO_embedd_Control = (
        f"{network_path}_GO_Embeddings_{annotation}_PPI_{tissue}_{cell}_{dim}_PPMI_Control.csv"
    )
    GO_embedd_Cancer = (
        f"{network_path}_GO_Embeddings_{annotation}_PPI_{cancer}_{dim}_PPMI_Cancer.csv"
    )

    # Top moving GO terms:

    ranking_p = f"{moving_path}Rank_movement_{tissue}_PPMI_{annotation}.csv"
    Ranking_db = pd.read_csv(ranking_p, index_col=0)

    # Load and preprocess the data:

    # Cancer:

    Gene_embedd_Cancer_db = pd.DataFrame(
        np.load(Gene_embedd_Cancer, allow_pickle=True)
    )
    genes_Cancer_db = list(pd.read_csv(genes_Cancer)["0"])
    genes_Cancer_db = [str(i) for i in genes_Cancer_db]
    Gene_embedd_Cancer_db.index = genes_Cancer_db

    GO_embedd_Cancer_db = pd.read_csv(GO_embedd_Cancer, index_col=0)
    GO_embedd_Cancer_db[GO_embedd_Cancer_db < 0] = 0

    # Control:

    Gene_embedd_Control_db = pd.DataFrame(
        np.load(Gene_embedd_Control, allow_pickle=True)
    )
    genes_Control_db = list(pd.read_csv(genes_Control)["0"])
    genes_Control_db = [str(i) for i in genes_Control_db]
    Gene_embedd_Control_db.index = genes_Control_db

    GO_embedd_Control_db = pd.read_csv(GO_embedd_Control, index_col=0)
    GO_embedd_Control_db[GO_embedd_Control_db < 0] = 0

    # Coordinates: 

    Gene_Cancer_Coordinates = Gene_embedd_Cancer_db.dot(Gene_midle_Cancer_db)
    Gene_Control_Coordinates = Gene_embedd_Control_db.dot(Gene_midle_Control_db)

    # Subset GO by the top 100 moving:

    vector_Ranking_Cancer = GO_embedd_Cancer_db[Ranking_db.index]
    vector_Ranking_Control = GO_embedd_Control_db[Ranking_db.index]

    # Calculate distances:

    Cancer_distances = Distance_GO_Gene( vector_Ranking_Cancer, Gene_Cancer_Coordinates)
    Control_distances = Distance_GO_Gene(vector_Ranking_Control, Gene_Control_Coordinates)

    # Add more information to the final data frame:

    Matrix_GO_Gene = pd.read_csv(
        f"{data_dir}/_Matrix_Genes_GO_BP_PPI.csv",
        index_col=0,
        dtype={0: str},
    )

    Matrix_GO_Gene.index = [str(i) for i in Matrix_GO_Gene.index] # added line to fix access type mismatch when comparing

    # Cancer:

    GO_Gene_Matrix_Cancer = Matrix_GO_Gene[
        Matrix_GO_Gene.index.isin(genes_Cancer_db)
    ]
    GO_Gene_Matrix_Cancer.index = [str(i) for i in GO_Gene_Matrix_Cancer.index]

    # Control:        
    
    GO_Gene_Matrix_Control = Matrix_GO_Gene[
        Matrix_GO_Gene.index.isin(genes_Control_db)
    ]
    GO_Gene_Matrix_Control.index = [str(i) for i in GO_Gene_Matrix_Control.index]

    # Compute the information

    Cancer_Cosine_it = Generate_Final_DataFrame(
        GO_Gene_Matrix_Cancer, Cancer_distances
    )
    Control_Cosine_it = Generate_Final_DataFrame(
        GO_Gene_Matrix_Control, Control_distances
    )

    # Save the data frame:

    Cancer_Cosine_it.to_csv(
        f"{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv", header=True, index=True
    )
    Control_Cosine_it.to_csv(
        f"{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv", header=True, index=True
    )


log_file = open(snakemake.log[0], "w")
sys.stdout = log_file

processes = []
for cancer, tissue, cell in zip(snakemake.params.Cancer_list, snakemake.params.Normal_tissue_list, snakemake.params.Control_list):
    processes.append((
            cancer, 
            tissue, 
            cell, 
            snakemake.params.optimal_dim, 
            snakemake.params.data_path,
            snakemake.params.annotation))

with Pool(len(processes)) as pool:
    pool.map(solve, processes)
    pool.close()
    pool.join()

log_file.close()