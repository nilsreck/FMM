import sys
import numpy as np
import pandas as pd

# 9.1 Movement of the GO terms:
# 9.2 Link movement with cancer (globally):

def Movement_Ranking(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    optimal_dim,
    data_dir,
    matrix="PPMI",
    annotation="Leaf",
):
    # Paths:

    cosine_path = f"{data_dir}/FMM/"
    Filter_path = f"{data_dir}/"
    movement_path = f"{data_dir}/"

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        Cancer_structure = pd.read_csv(
            f"{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_Leaf.csv",
            index_col=0,
        )
        Control_structure = pd.read_csv(
            f"{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_Leaf.csv",
            index_col=0,
        )

        # Filter the matrices with the common set:

        common_set = list(
            pd.read_csv(f"{Filter_path}Common_Set_{tissue}_Leaf.csv")["0"]
        )

        Cancer_structure = Cancer_structure.loc[common_set, common_set]
        Control_structure = Control_structure.loc[common_set, common_set]

        # Substract the cancer to the control:

        Substract = Control_structure - Cancer_structure

        # Calculate the total movement of the GO terms in the space by the 1-norm of the vectors:

        Movement = Substract.apply(np.linalg.norm, axis=1)

        # Rank the movement:

        Movement = Movement.sort_values(ascending=False)

        # Save the rank:

        Movement.to_csv(
            f"{movement_path}Rank_movement_{tissue}_{matrix}_{annotation}.csv"
        )


Movement_Ranking(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.optimal_dim, 
    snakemake.params.data_path,
    matrix="PPMI",
    annotation=snakemake.params.annotation
)