import sys
import numpy as np
import pandas as pd

# 2. Generate the PPMI Matrices:

def deep_walk_ppmi(adj_matrix, context_window=10):
    """
    Input: Adj Matrix as numpy array
    """
    degrees = np.sum(adj_matrix, axis=0)
    volume_of_graph = sum(degrees)
    diag_degrees_inv = np.diag([1 / float(i) for i in degrees])

    pmi = np.zeros([len(adj_matrix), len(adj_matrix)])

    transition_matrix = np.matmul(diag_degrees_inv, adj_matrix)
    pmi_temp = transition_matrix
    pmi += transition_matrix

    for r in range(1, context_window):
        print("Iteration:")
        pmi_temp = np.matmul(pmi_temp, transition_matrix)
        pmi += pmi_temp

    pmi = pmi / float(context_window)
    pmi = volume_of_graph * np.matmul(pmi, diag_degrees_inv)

    pmi_matrix = np.log(pmi, out=np.zeros_like(pmi), where=(pmi != 0))

    for i in pmi_matrix:
        i[i < 0] = 0

    return pmi_matrix

def Generate_PPMI(
    Cancer_type_list, Normal_Tissue_list, Cell_type_list, data_dir, context_window=10
):
    """
    This function generates the PPMI matrix from the Adj Matrices and save it as a np:

        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - context_window  : int, contex window for the PPMI matrix.
    """

    # Networks path:

    network_path = f"{data_dir}/"
    counter = 1

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Load the Adj:

        Cancer_adj = np.load(
            f"{network_path}Cancer_{cancer}_PPI.npy", allow_pickle=True
        )
        Control_adj = np.load(
            f"{network_path}Control_{tissue}_{cell}_PPI.npy", allow_pickle=True
        )

        # Generate the PPMI:

        PPMI_Cancer = deep_walk_ppmi(
            Cancer_adj, context_window
        )
        PPMI_Control = deep_walk_ppmi(
            Control_adj, context_window
        )

        # Save the PPMI:

        np.save(f"{network_path}Cancer_PPMI_{cancer}", PPMI_Cancer)
        np.save(f"{network_path}Control_PPMI_{tissue}_{cell}", PPMI_Control)

        print(f"{counter}/{len(Cancer_type_list)}")

        counter = counter + 1


Generate_PPMI(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.data_path,
    context_window=10
)