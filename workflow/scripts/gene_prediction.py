import sys
import numpy as np
import pandas as pd

# 16. Do the genes predictions and validate them with an enrichment analyses:

def Do_Fold_Enrichment_Complete(Normal_tissue_list, data_dir):
    path_rank = f"{data_dir}/"
    gene_GO_path = f"{data_dir}/"

    # Open the file:

    with open(f"{gene_GO_path}Fold_Rank_Table_Dos.txt", "a") as the_file:
        # Write the name of the columns:

        the_file.write(f"Cancer Type\t")
        the_file.write(f"Moving Genes\t")
        the_file.write(f"Stable Genes\n")

        for tissue in Normal_tissue_list:
            # Load movement and choose the most moving:

            rank = pd.read_csv(f"{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv")
            moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2 * np.std(rank["0"]))][
                "Unnamed: 0"
            ]

            # Load the bibliography information:

            counter_list = pd.read_csv(f"{gene_GO_path}Cancer_Count_{tissue}.csv")
            counter_list = counter_list.drop_duplicates("name")

            # Load cosines:

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

            # For all the GO terms.

            # Subset the data:

            gene_common = Cancer_Cosine_it.index[
                Cancer_Cosine_it.index.isin(Control_Cosine_it.index)
            ]
            Cancer_Cosine_it_com = Cancer_Cosine_it.loc[gene_common]
            Control_Cosine_it_com = Control_Cosine_it.loc[gene_common]

            moving_GO = list(
                set(Cancer_Cosine_it_com.columns).intersection(set(moving_GO))
            )

            # Filter the moving GO terms (based on the distribution)

            Cancer_structure_GO1 = Cancer_Cosine_it_com[moving_GO]
            Control_structure_GO1 = Control_Cosine_it_com[moving_GO]

            # Analyze the movement of the genes together:

            move = Control_structure_GO1 - Cancer_structure_GO1

            # Get the sets to analyze:

            # A) Moving Vs Not moving (no matter the direction):

            move_abs = abs(move)
            Restuls_abs = pd.DataFrame(move_abs.T.max())
            Restuls_abs["GO"] = move_abs.T.idxmax()
            Restuls_abs["Gene"] = list(Restuls_abs.index)
            Restuls_abs = Restuls_abs.sort_values(by=0, ascending=False)
            Restuls_abs = Restuls_abs.reset_index(drop=True)

            # Normalize the distribution by using a sqrt transformation:

            Restuls_abs["Original"] = Restuls_abs[0]
            z = [np.log(i) for i in list(Restuls_abs[0])]
            Restuls_abs[0] = z

            # Choose the tails:

            Restuls_abs_top = Restuls_abs.head(round(len(Restuls_abs) * 5 / 100))
            Restuls_abs_tail = Restuls_abs.tail(round(len(Restuls_abs) * 5 / 100))

            # To save the predictions:

            gene_names = counter_list[
                counter_list.initial_alias.isin(Restuls_abs_top.Gene)
            ]
            gene_names = gene_names.drop(["Unnamed: 0", "initial_alias"], axis=1)
            gene_names = gene_names.drop_duplicates("name")
            gene_names = gene_names.reset_index(drop=True)

            gene_names.to_csv(
                f"{path_rank}Predictions_Rank_{tissue}.csv", header=True, index=True
            )

            # Get the total set

            genes_succes = list(counter_list[counter_list.Counts > 0]["initial_alias"])
            genes_succes = [str(i) for i in genes_succes]

            Total_succes = len(genes_succes)
            Total = len(counter_list)

            # Do the enrichment analyses:

            set_stable_succes = len(
                set(Restuls_abs_tail.Gene).intersection(set(genes_succes))
            )
            set_stable_total = len(Restuls_abs_tail)

            set_moving_succes = len(
                set(Restuls_abs_top.Gene).intersection(set(genes_succes))
            )
            set_moving_total = len(Restuls_abs_top)

            fold_stb, p_value_stb = Fold_enriched(
                Total, Total_succes, set_stable_total, set_stable_succes
            )
            fold_mv, p_value_mv = Fold_enriched(
                Total, Total_succes, set_moving_total, set_moving_succes
            )

            the_file.write(f"{tissue}\t")
            the_file.write(f"{round(fold_mv, 3)} ({p_value_mv})\t")
            the_file.write(f"{round(fold_stb, 3)} ({p_value_stb})\n")

        # Close the file:

        the_file.close()


Do_Fold_Enrichment_Complete(snakemake.params.Normal_tissue_list, snakemake.params.data_path)