import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from sklearn.metrics import pairwise_distances

# 5 Generate the FMMs

def Embedding_Structure(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    data_dir,
    Gene_Matrix=False,
    matrix="PPMI",
    annotation="GO",
):
    """
    This function calculates the relative position of any embedded entity (GO terms by defaut). It producess a csv with the corresponging
    results.

        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - Gene_Matrix     : bool, if the entities are genes.
                - matrix          : string, network matrix representation.
    """

    # Path:

    network_path = f"{data_dir}/Embeddings/"
    save_cosine = f"{data_dir}/FMM/"

    control = 0

    # Get the structure of the spaces:

    path_Cancer = None
    path_Control = None
    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):


        for dim in dimension_list:
            if Gene_Matrix == False:
                if annotation == "GO":
                    path_Cancer = f"{network_path}_GO_Embeddings_GO_PPI_{cancer}_{dim}_{matrix}_Cancer.csv"
                    path_Control = f"{network_path}_GO_Embeddings_GO_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv"

                elif annotation == "Reactome":
                    path_Cancer = f"{network_path}_GO_Embeddings_Reactome_PPI_{cancer}_{dim}_{matrix}_Cancer.csv"
                    path_Control = f"{network_path}_GO_Embeddings_Reactome_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv"

                elif annotation == "Leaf":
                    path_Cancer = f"{network_path}_GO_Embeddings_Leaf_PPI_{cancer}_{dim}_{matrix}_Cancer.csv"
                    path_Control = f"{network_path}_GO_Embeddings_Leaf_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv"

                embeddings_Cancer = pd.read_csv(path_Cancer, index_col=0).T
                embeddings_Control = pd.read_csv(path_Control, index_col=0).T

                embeddings_Cancer[embeddings_Cancer < 0] = 0
                embeddings_Control[embeddings_Control < 0] = 0

                embeddings_Cancer_index = embeddings_Cancer.index
                embeddings_Control_index = embeddings_Control.index

                embeddings_Cancer = np.array(embeddings_Cancer)
                embeddings_Control = np.array(embeddings_Control)

                cosine_Cancer = pairwise_distances(embeddings_Cancer, metric="cosine")
                cosine_Control = pairwise_distances(embeddings_Control, metric="cosine")

                cosine_Cancer = pd.DataFrame(
                    cosine_Cancer, embeddings_Cancer_index, embeddings_Cancer_index
                )
                cosine_Control = pd.DataFrame(
                    cosine_Control, embeddings_Control_index, embeddings_Control_index
                )

                if annotation == "GO":
                    cosine_Cancer.to_csv(
                        f"{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}.csv",
                        header=True,
                        index=True,
                    )
                    cosine_Control.to_csv(
                        f"{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv",
                        header=True,
                        index=True,
                    )

                    control = control + 2

                elif annotation == "Reactome":
                    cosine_Cancer.to_csv(
                        f"{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Reactome.csv",
                        header=True,
                        index=True,
                    )
                    cosine_Control.to_csv(
                        f"{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Reactome.csv",
                        header=True,
                        index=True,
                    )

                    control = control + 2

                elif annotation == "Leaf":
                    cosine_Cancer.to_csv(
                        f"{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Leaf.csv",
                        header=True,
                        index=True,
                    )
                    cosine_Control.to_csv(
                        f"{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Leaf.csv",
                        header=True,
                        index=True,
                    )

                    control = control + 2

                print(f"{control}/{len(Cancer_type_list) * 2 * len(dimension_list)}")

            else:
                path_Cancer = (
                    f"{network_path}_G_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy"
                )
                path_Control = (
                    f"{network_path}_G_Matrix_{dim}_PPI_{cancer}_Adj_Cancer.npy"
                )

                embeddings_Cancer = np.load(path_Cancer, allow_pickle=True)
                embeddings_Control = np.load(path_Control, allow_pickle=True)

                Cancer_cosine = pairwise_distances(embeddings_Cancer, metric="cosine")
                Control_cosine = pairwise_distances(embeddings_Control, metric="cosine")

                Cancer_cosine.to_csv(
                    f"{save_cosine}Gene_Cosine_Cancer_{cancer}_{dim}_{matrix}.csv",
                    header=True,
                    index=True,
                )
                Control_cosine.to_csv(
                    f"{save_cosine}Gene_Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv",
                    header=True,
                    index=True,
                )

def solve(process):
    Embedding_Structure(*process)

processes = []
for cancer, tissue, control in zip(
        snakemake.params.Cancer_list, 
        snakemake.params.Normal_tissue_list, 
        snakemake.params.Control_list
    ):
    processes.append((
            [cancer],
            [tissue],
            [control],
            [int(open(snakemake.input.dimension, "w").name.split("/")[-1].split("_")[-1])],
            snakemake.params.data_path,
            False,
            "PPMI",
            snakemake.params.annotation)
            )
with Pool(4) as pool:
    pool.map(solve, processes)
    pool.close()
    pool.join()

# Embedding_Structure(
#     snakemake.params.Cancer_list,
#     snakemake.params.Normal_tissue_list,
#     snakemake.params.Control_list,
#     [int(open(snakemake.input.dimension, "w").name.split("/")[-1])],
#     snakemake.params.data_path,
#     Gene_Matrix=False,
#     matrix="PPMI",
#     annotation=snakemake.params.annotation,
# )