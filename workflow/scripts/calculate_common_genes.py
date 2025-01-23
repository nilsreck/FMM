import sys
import pandas as pd
import matplotlib.pyplot as plt


# 6. Filter FMMs to contain only annotations shared by cancer and control.
# function was taken from Cancer_plot_Scripts.py

def Common_GO_Terms(
    Cancer_type_list, Normal_Tissue_list, Cell_type_list, data_dir, annotation="GO"
):
    """
    This gets the common set of GO terms between cancer/control. Then, the annotated space can be filtered for comparisons.

        Inputs:
                - Cancer_type_list     : string list, Name of the Cancer.
                - Normal_Tissue_list   : string list, Name of the normal tissue.
                - Cell_type_list       : string list, Name of the normal cell in the tissue.
    """

    path = f"{data_dir}/"
    save = f"{data_dir}/"

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Load GO Matrix:

        if annotation == "GO":
            GO_Matrix = pd.read_csv(
                f"{path}_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv",
                index_col=0,
                dtype={0: str},
            )
            GO_Matrix.index = GO_Matrix.index.astype(str)

        else:
            GO_Matrix = pd.read_csv(
                f"{path}_Matrix_Genes_GO_BP_PPI.csv", index_col=0, dtype={0: str}
            )
            GO_Matrix.index = GO_Matrix.index.astype(str)

        # Load the GO/gene matrix:

        Cancer_genes = pd.read_table(
            f"{path}Cancer_{cancer}_Genes.csv", header=0, dtype={0: str}
        )
        Cancer_genes_list = list(Cancer_genes["0"])

        Control_genes = pd.read_table(
            f"{path}Control_{tissue}_{cell}_Genes.csv", header=0, dtype={0: str}
        )
        Control_genes_list = list(Control_genes["0"])

        # Subset the GO embeddings to keep the same genes:

        GO_Gene_Matrix_Cancer = GO_Matrix[GO_Matrix.index.isin(Cancer_genes_list)]
        GO_Gene_Matrix_Control = GO_Matrix[GO_Matrix.index.isin(Control_genes_list)]

        # Reorder the Matrices:

        GO_Gene_Matrix_Cancer = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
        GO_Gene_Matrix_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]

        # Cancer:

        GO_Cancer = GO_Gene_Matrix_Cancer.sum(axis=0)
        GO_terms_filtered_Cancer = set(GO_Cancer[GO_Cancer >= 3].index)

        # Control:

        GO_Control = GO_Gene_Matrix_Control.sum(axis=0)
        GO_terms_filtered_Control = set(GO_Control[GO_Control >= 3].index)

        # Intersecions:

        Common_Annotations = GO_terms_filtered_Cancer.intersection(
            GO_terms_filtered_Control
        )

        if annotation == "GO":
            pd.DataFrame(Common_Annotations).to_csv(f"{save}Common_Set_{tissue}.csv")

        else:
            pd.DataFrame(Common_Annotations).to_csv(
                f"{save}Common_Set_{tissue}_{annotation}.csv"
            )

Common_GO_Terms(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.data_path,
    annotation=snakemake.params.annotation
)
