import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, hypergeom


# 9.1 Movement of the GO terms:
# 9.2 Link movement with cancer (globally):
# 9.3 Enrichment analyses validation (shifted annotations)

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



def Global_moving_in_the_space(data_dir, annotation):
    """
    This function test if the cancer-related annotations are moving statistically more than the rest.
    """

    # Load our set of cancer-related annotations:

    with open(
        f"{data_dir}/enriched_in_cancer_gobp_terms_cosmic.txt", "r"
    ) as fp:
        cancer_related_go_terms = json.load(fp)

    # Do the test

    cancers = ["breast", "prostate", "lung", "colon"]

    for i in cancers:
        Ranking = pd.read_csv(
            f"{data_dir}/Rank_movement_{i}_PPMI_{annotation}.csv",
            index_col=0,
            names=["norm"],
            header=0,
        )

        # Globally validate the movement is connected with cancer:

        a = Ranking.loc[~Ranking.index.isin(cancer_related_go_terms)]["norm"]
        b = Ranking.loc[Ranking.index.isin(cancer_related_go_terms)]["norm"]
        print(i, mannwhitneyu(a, b, alternative="less"))


def moving_in_the_space(data_dir, annotation):
    """
    This function performs the enrichment analyses of cancer-related annotations.
    """

    # Load our set of cancer-related annotations:

    with open(
        f"{data_dir}/enriched_in_cancer_gobp_terms_cosmic.txt", "r"
    ) as fp:
        cancer_related_go_terms = json.load(fp)

    # Set the number of shifted or stable annotations based on 2std:

    percentile_moving = {"lung": 75, "breast": 58, "colon": 68, "prostate": 49}
    percentile_stable = {"lung": 15, "breast": 29, "colon": 68, "prostate": 26}

    # To save the values of the enrichment analyses

    percentage_most_moving = []
    percentage_less_moving = []
    percentage_by_random = []

    success_in_the_population = len(cancer_related_go_terms)

    # Do the enrichment per each cancer type:

    cancers = ["breast", "prostate", "lung", "colon"]
    for i in cancers:
        number_of_draws = percentile_moving[i]
        number_of_draws_stable = percentile_stable[i]

        print("\n")
        Ranking = pd.read_csv(
            f"{data_dir}/Rank_movement_{i}_PPMI_{annotation}.csv",
            index_col=0,
            names=["norm"],
            header=0,
        )

        # Observed by random
        population_size = len(Ranking)
        probability_of_success = success_in_the_population / population_size
        expected_random = number_of_draws * probability_of_success
        print(expected_random)

        top_moving_norm_df = Ranking.sort_values(by="norm", ascending=False)[
            0:number_of_draws
        ]
        top_moving_norm = set(top_moving_norm_df.index)
        top_success = len(top_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_head = round(
            hypergeom.sf(
                top_success, population_size, success_in_the_population, number_of_draws
            ),
            5,
        )

        less_moving_norm_df = Ranking.sort_values(by="norm", ascending=False)[
            -number_of_draws_stable:
        ]
        less_moving_norm = set(less_moving_norm_df.index)
        less_success = len(less_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_tail = round(
            hypergeom.sf(
                less_success,
                population_size,
                success_in_the_population,
                number_of_draws_stable,
            ),
            5,
        )

        print(
            f"% cancer related genes in the top {number_of_draws}:",
            top_success,
            "and the corresponding p-value",
            pvalue_of_head,
        )
        print(
            f"% cancer related genes in the less {number_of_draws_stable}:",
            less_success,
            "and the corresponding p-value",
            pvalue_of_tail,
        )

        percentage_most_moving.append(top_success * 100 / number_of_draws)
        percentage_less_moving.append(less_success * 100 / number_of_draws_stable)
        percentage_by_random.append(expected_random * 100 / number_of_draws)

    # Do the plot:

    plt.rcParams.update({"font.size": 15})
    categories = snakemake.params.Normal_tissue_list

    if not snakemake.params.plot_bar:

        plt.figure(figsize=(8, 5))
        plt.plot(percentage_most_moving, "--o", label=f" shifted annotations", linewidth=2)
        plt.plot(percentage_less_moving, "--o", label=f" stable annotations", linewidth=2)
        plt.plot(percentage_by_random, "--o", label=f" expected by random", linewidth=2)
        plt.ylabel("% cancer-related annotations", fontsize=16, fontweight="bold")
        plt.xticks(range(4), categories, fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=16)
        plt.title("A", fontdict=None, loc='center', pad=None)
        plt.tight_layout()
        plt.savefig(f"{data_dir}/cancer_enrichments_2std.svg", dpi=500)

    else:

        bar_width = 0.25  # Breite der Säulen
        x_indices = np.arange(len(categories))
        plt.figure(figsize=(8, 8))
        plt.bar(x_indices - bar_width, percentage_most_moving, width=bar_width, label='shifted annotations', edgecolor=None)
        plt.bar(x_indices + bar_width, percentage_less_moving, width=bar_width, label='stable annotations', edgecolor=None)
        plt.bar(x_indices , percentage_by_random, width=bar_width, label='expected by random', edgecolor=None)
        for i in range(0, 13, 2):
            plt.axhline(y=i, color='gray', linestyle='-', linewidth=1)
        plt.ylabel('% cancer-related annotations', fontsize=24 , fontweight='bold')
        plt.xticks(x_indices, categories, fontsize=24)
        plt.yticks(fontsize=24)
        #plt.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=16)
        plt.legend(fontsize=16)
        plt.title("A", fontdict=None, loc='center', pad=None)
        #plt.tight_layout()
        plt.savefig(f"{data_dir}/cancer_enrichments_2std.svg", dpi=500)

    #plt.show()

Movement_Ranking(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.optimal_dim, 
    snakemake.params.data_path,
    matrix="PPMI",
    annotation=snakemake.params.annotation
)

Global_moving_in_the_space(snakemake.params.data_path, snakemake.params.annotation)

moving_in_the_space(snakemake.params.data_path, snakemake.params.annotation)