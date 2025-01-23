import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2_unweighted


# 1.1 Generate the networks:
# 1.2 Plot the Venn Diagrams:
# 1.3 Get the statistics for the networks:


def edges(A):
    m = A.shape[0]
    r, c = np.triu_indices(m, 1)
    A = A[r, c]
    np.count_nonzero(A == 1.0)
    return np.count_nonzero(A == 1.0)


def Density(Cancer_nodes, Cancer_edges):
    PC = (Cancer_nodes * (Cancer_nodes - 1)) / 2
    return Cancer_edges / PC * 100


def control_infected_venn(list_cotrol, list_infected, Cancer_type):
    cancer_genes = set(list_infected[0])
    control_genes = set(list_cotrol[0])
    common_genes = len(cancer_genes.intersection(control_genes))
    total_genes = len(cancer_genes.union(control_genes))
    unique_control = len(control_genes) - common_genes
    unique_cancer = len(cancer_genes) - common_genes

    label_control = "Control"
    label_infected = "Cancer"

    out = venn2_unweighted(
        subsets=(unique_control, unique_cancer, common_genes),
        set_labels=(label_control, label_infected),
        alpha=0.5,
        subset_label_formatter=lambda x: str(x)
        + "\n("
        + f"{(x/total_genes):1.1%}"
        + ")",
    )
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)

    plt.title(str(Cancer_type).title(), fontsize=16, fontweight="bold")



def Create_Tissue_Specific_Network(
    Cancer_type, Normal_Tissue, Cell_type, data_dir, plot=True
):
    # Load the files:

    # Load genes that are expressed in each tissue and cancer type:

    save_path = f"{data_dir}/"

    cancer_genes = list(
        pd.read_csv(f"{save_path}Cancer_{Cancer_type}_Genes.csv", dtype={0: str})["0"]
    )
    control_genes = list(
        pd.read_csv(
            f"{save_path}Control_{Normal_Tissue}_{Cell_type}_Genes.csv", dtype={0: str}
        )["0"]
    )

    # Load the human PPI network:

    file_name_array = "Human_Biogrid_Adj"
    file_name_genelist = "Human_Biogrid_Genes"
    data_path = save_path

    PPI_adj = np.load(data_path + file_name_array + "_PPI_.npy", allow_pickle=True)
    PPI_genes = list(
        pd.read_table(
            data_path + file_name_genelist + "_PPI_", header=None, dtype={0: str}
        )[0]
    )

    PPI_df = pd.DataFrame(PPI_adj, PPI_genes, PPI_genes)

    # Create both networks from the original human PPI networks:

    PPI_Cancer = PPI_df.loc[cancer_genes, cancer_genes]
    PPI_Control = PPI_df.loc[control_genes, control_genes]

    Cancer_np = np.array(PPI_Cancer)
    Control_np = np.array(PPI_Control)

    # Save the networks:

    np.save(save_path + "Cancer_" + str(Cancer_type) + "_PPI", Cancer_np)
    np.save(
        save_path + "Control_" + str(Normal_Tissue) + "_" + str(Cell_type) + "_PPI",
        Control_np,
    )

    # Print Venn:

    if plot == True:
        control_infected_venn(
            pd.DataFrame(cancer_genes), pd.DataFrame(control_genes), Cancer_type
        )


def Venn_Diagram_Cancer_Networks(Cancer_type_list, Normal_Tissue_list, Cell_type_list, data_dir):
    """
    This fucntions plots the Venn diagram of the tissue specific Networks:

        Inputs:
                - Cancer_type   : string list, Name of the Cancer.
                - Normal_Tissue : string list, Name of the normal tissue.
                - Cell_type     : string list, Name of the normal cell in the tissue.
    """

    # Networks path:

    network_path = f"{data_dir}/"

    # Plot Proprieties:

    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({"font.size": 17})
    plt.rc("font", weight="bold")

    # Set the grid:

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.tight_layout(pad=0.5)

    row_count = 0
    column_count = 0

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Load the data:

        path_Cancer = f"{network_path}Cancer_{cancer}_Genes.csv"
        path_Control = f"{network_path}Control_{tissue}_{cell}_Genes.csv"

        Cancer_net = pd.read_csv(path_Cancer)
        Control_net = pd.read_csv(path_Control)

        # Get the counts for the Venn Diagram:

        cancer_genes = set(Cancer_net["0"])
        control_genes = set(Control_net["0"])
        common_genes = len(cancer_genes.intersection(control_genes))
        total_genes = len(cancer_genes.union(control_genes))
        unique_control = len(control_genes) - common_genes
        unique_cancer = len(cancer_genes) - common_genes

        # Plot the Venn:

        venn2_unweighted(
            subsets=(unique_control, unique_cancer, common_genes),
            set_labels=("", ""),
            alpha=0.5,
            subset_label_formatter=lambda x: str(x)
            + "\n("
            + f"{(x / total_genes):1.1%}"
            + ")",
            ax=axes[row_count, column_count],
        )

        axes[row_count, column_count].set_title(
            f"{tissue}".capitalize(), fontweight="bold", fontsize=23, y=1
        )

        # Move in the grid:

        if column_count != 1:
            column_count = column_count + 1
        else:
            column_count = 0
            row_count = row_count + 1

    # Set the legend for all the plots:

    fig.legend(
        labels=["Control", "Cancer"],
        borderaxespad=0.1,
        bbox_to_anchor=(0.5, 0.54),
        loc="upper center",
        frameon=True,
        ncol=3,
    )

    # Save the plot:

    fig.savefig(
        f"{network_path}Venn_Diagrams_Networks.svg",
        format="svg",
        dpi=600,
        bbox_inches="tight",
    )


def Network_Statistics(Cancer_type_list, Normal_Tissue_list, Cell_type_list, data_dir):
    """
    This function generate a csv with the network statistics (nodes, edges, density)

        Inputs:
                - Cancer_type   : string list, Name of the Cancer.
                - Normal_Tissue : string list, Name of the normal tissue.
                - Cell_type     : string list, Name of the normal cell in the tissue.
    """
    # Networks path:

    network_path = f"{data_dir}/"

    # Open the file to write:

    with open(network_path + "Network_Statistics.txt", "a") as the_file:
        # Write the titles per each column:

        the_file.write("# Name" + "\t")
        the_file.write("# Nodes" + "\t")
        the_file.write("# Edges" + "\t")
        the_file.write("% Density" + "\n")

        # Start writting the file:

        for cancer, tissue, cell in zip(
            Cancer_type_list, Normal_Tissue_list, Cell_type_list
        ):
            # Load the Adj matrices:

            Cancer_adj = np.load(
                f"{network_path}Cancer_{cancer}_PPI.npy", allow_pickle=True
            )
            Control_adj = np.load(
                f"{network_path}Control_{tissue}_{cell}_PPI.npy", allow_pickle=True
            )

            # Get info of statistics:

            Cancer_nodes = len(Cancer_adj)
            Cancer_edges = edges(Cancer_adj)
            Cancer_density = round(Density(Cancer_nodes, Cancer_edges), 3)

            Control_nodes = len(Control_adj)
            Control_edges = edges(Control_adj)
            Control_density = round(Density(Control_nodes, Control_edges), 3)

            # Write the information:

            the_file.write(f"{cancer}" + "\t")
            the_file.write(f"{str(Cancer_nodes)}" + "\t")
            the_file.write(f"{str(Cancer_edges)}" + "\t")
            the_file.write(f"{str(Cancer_density)}" + "\n")

            the_file.write(f"{tissue}" + "\t")
            the_file.write(f"{str(Control_nodes)}" + "\t")
            the_file.write(f"{str(Control_edges)}" + "\t")
            the_file.write(f"{str(Control_density)}" + "\n")

        # Close the file:

        the_file.close()


for control, tissue, cancer in zip(snakemake.params.Control_list, snakemake.params.Normal_tissue_list, snakemake.params.Cancer_list):
    Create_Tissue_Specific_Network(
        cancer, tissue, control, snakemake.params.data_path, plot=False
    )
Venn_Diagram_Cancer_Networks(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.data_path
)
Network_Statistics(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.data_path
)