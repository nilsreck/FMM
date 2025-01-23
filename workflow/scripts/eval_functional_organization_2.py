import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.stats import mannwhitneyu


# 9.6 Assess the functional organization of genes and functions in the embedding space

def solve(process):
    return Demonstrate_Gene_GO_Org(*process)


def Distance_Check(GO_Gene_Matrix, distances):
    # Lists to fill:

    x = []  # Annotated distance
    y = []  # Non-Annotated distance

    # Iteration per gene:

    for gene in list(distances.index):
        it_gene = distances.loc[gene]
        it_annotation = GO_Gene_Matrix.loc[str(gene)]

        for annotation in list(it_gene.index):
            if it_annotation[annotation] == 1:
                x.append(it_gene[annotation])
            else:
                y.append(it_gene[annotation])

    # Return the dataFrame:

    return (x, y)


def Demonstrate_Gene_GO_Org(cancer, tissue, cell, dim, data_dir, annotation):
    # Paths:

    network_path = f"{data_dir}/"
    gene_GO_path = f"{data_dir}/"
    matrix_path = f"{data_dir}/"

    # Open the file:

    # with open(gene_GO_path + "Organization_Common_Space.txt", "a") as the_file:
    #     # Write the columns names in the file:

    #     the_file.write("# Sample" + "\t")
    #     the_file.write("# Annotated_NonAnnotated(less distance)" + "\t")
    #     the_file.write("# Distance Annotated" + "\t")
    #     the_file.write("# Distance Not-Annotated" + "\n")

        # For each combination:


    # for cancer, tissue, cell in zip(
    #     Cancer_list, Normal_tissue_list, Cell_type_list
    # ):

    # Load the data:

    # Genes:

    genes_Control = f"{network_path}Control_{tissue}_{cell}_Genes.csv"
    genes_Cancer = f"{network_path}Cancer_{cancer}_Genes.csv"


    genes_Cancer_db = list(pd.read_csv(genes_Cancer)["0"])
    genes_Cancer_db = [str(i) for i in genes_Cancer_db]
    genes_Control_db = list(pd.read_csv(genes_Control)["0"])
    genes_Control_db = [str(i) for i in genes_Control_db]

    # Annotations:

    Matrix_GO_Gene = pd.read_csv(
        f"{matrix_path}_Matrix_Genes_GO_BP_PPI.csv", index_col=0, dtype={0: str}
    )

    Matrix_GO_Gene.index = [str(i) for i in Matrix_GO_Gene.index] # added line to fix access type mismatch when comparing

    # Cancer annotations:

    GO_Gene_Matrix_Cancer = Matrix_GO_Gene[
        Matrix_GO_Gene.index.isin(genes_Cancer_db)
    ]
    GO_Gene_Matrix_Cancer.index = [str(i) for i in GO_Gene_Matrix_Cancer.index]

    # Control annotations:

    GO_Gene_Matrix_Control = Matrix_GO_Gene[
        Matrix_GO_Gene.index.isin(genes_Control_db)
    ]
    GO_Gene_Matrix_Control.index = [
        str(i) for i in GO_Gene_Matrix_Control.index
    ]

    # Distances:

    Cancer_Cosine_it = pd.read_csv(
        f"{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv",
        index_col=0,
        dtype={0: str},
    )
    Control_Cosine_it = pd.read_csv(
        f"{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv",
        index_col=0,
        dtype={0: str},
    )

    Cancer_Cosine_it.index = [str(i) for i in Cancer_Cosine_it.index]
    Control_Cosine_it.index = [str(i) for i in Control_Cosine_it.index]

    Cancer_Cosine_it = Cancer_Cosine_it.drop(
        ["Group", "Shape", "size", "Min_Dist", "Min_GO"], axis=1
    )
    Control_Cosine_it = Control_Cosine_it.drop(
        ["Group", "Shape", "size", "Min_Dist", "Min_GO"], axis=1
    )

    # Compute the sets:

    x_cancer, y_cancer = Distance_Check(GO_Gene_Matrix_Cancer, Cancer_Cosine_it)
    x_control, y_control = Distance_Check(
        GO_Gene_Matrix_Control, Control_Cosine_it
    )

    # Compare distance distributions:

    p_value_cancer = mannwhitneyu(x_cancer, y_cancer, alternative="less").pvalue
    p_value_control = mannwhitneyu(
        x_control, y_control, alternative="less"
    ).pvalue

    #print([tissue, p_value_cancer, p_value_control], flush=True)
    return [tissue, p_value_cancer, x_cancer, y_cancer, p_value_control, x_control, y_control]

            # # Write into the file:

            # the_file.write(f"Cancer {tissue}\t")
            # the_file.write(f"{p_value_cancer}\t")
            # the_file.write(f"{np.mean(x_cancer)} ({np.std(x_cancer)})\t")
            # the_file.write(f"{np.mean(y_cancer)} ({np.std(y_cancer)})\n")

            # the_file.write(f"Control {tissue}\t")
            # the_file.write(f"{p_value_control}\t")
            # the_file.write(f"{np.mean(x_control)} ({np.std(x_control)})\t")
            # the_file.write(f"{np.mean(y_control)} ({np.std(y_control)})\n")

        # Close the file:

        # the_file.close()
        ###################################


# Examples of ways for evaluating the organization:
# Semantic similarity of Top 500 closest/farthest functional annotation embedding vectors:

log_file = open(snakemake.log[0], "w")
sys.stdout = log_file

processes = []
for cancer, tissue, cell in zip(snakemake.params.Cancer_list, snakemake.params.Normal_tissue_list, snakemake.params.Control_list):
    processes.append([
            cancer, 
            tissue, 
            cell, 
            snakemake.params.optimal_dim, 
            snakemake.params.data_path,
            snakemake.params.annotation])


with Pool(len(processes)) as pool:
    results = pool.map(solve, processes)
    pool.close()
    pool.join()


with open(snakemake.params.data_path + "/Organization_Common_Space.txt", "a") as the_file:
    # Write the columns names in the file:

    the_file.write("# Sample" + "\t")
    the_file.write("# p-value" + "\t")
    the_file.write("# Distance Annotated mean (std)" + "\t")
    the_file.write("# Distance Not-Annotated mean (std)" + "\n")

    for result in results:

        tissue = result[0]
        p_value_cancer = result[1]
        x_cancer = result[2]
        y_cancer = result[3]
        p_value_control = result[4]
        x_control = result[5]
        y_control = result[6]

        the_file.write(f"Cancer {tissue}\t")
        the_file.write(f"{p_value_cancer}\t")
        the_file.write(f"{np.mean(x_cancer)} ({np.std(x_cancer)})\t")
        the_file.write(f"{np.mean(y_cancer)} ({np.std(y_cancer)})\n")

        the_file.write(f"Control {tissue}\t")
        the_file.write(f"{p_value_control}\t")
        the_file.write(f"{np.mean(x_control)} ({np.std(x_control)})\t")
        the_file.write(f"{np.mean(y_control)} ({np.std(y_control)})\n")

    the_file.close()

log_file.close()

# # Demonstrate that the space is functionally organized (between genes and annotations):
# Demonstrate_Gene_GO_Org(
#     snakemake.params.Cancer_list, 
#     snakemake.params.Normal_tissue_list, 
#     snakemake.params.Control_list, 
#     snakemake.params.optimal_dim, 
#     snakemake.params.data_path
# )