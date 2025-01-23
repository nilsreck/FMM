import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, hypergeom


# 9.3 Enrichment analyses validation (shifted annotations)


def Global_moving_in_the_space(data_dir, annotation):
    """
    This function test if the cancer-related annotations are moving statistically more than the rest.
    """

    # Load our set of cancer-related annotations:

    with open(f"{data_dir}/enriched_in_cancer_gobp_terms_cosmic.txt", "r") as fp:
        cancer_related_go_terms = json.load(fp)

    # Do the test

    cancers = snakemake.params.Normal_tissue_list

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
        print(f"validating movement to cancer relatedness using \"enriched_in_cancer_gobp_terms_cosmic\" for {i}: \n", 
            mannwhitneyu(a, b, alternative="less"), flush=True)


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
    
    #     percentile_moving = {"lung": 75, "breast": 58, "colon": 68, "prostate": 49}
    #     percentile_stable = {"lung": 15, "breast": 29, "colon": 68, "prostate": 26}

    # To save the values of the enrichment analyses

    percentage_most_moving = []
    percentage_less_moving = []
    percentage_by_random = []

    success_in_the_population = len(cancer_related_go_terms)

    # Do the enrichment per each cancer type:

    cancers = snakemake.params.Normal_tissue_list
    for i in cancers:

        Ranking = pd.read_csv(f"{data_dir}/Rank_movement_{i}_PPMI_{annotation}.csv",
            index_col=0,
            names=["norm"],
            header=0,
        )

        # Set the number of shifted or stable annotations based on 2std:
    
        if snakemake.params.shifted[i]: number_of_draws = int(snakemake.params.shifted[i])
        else: number_of_draws = len(Ranking[Ranking["norm"] > (np.mean(Ranking["norm"]) + 2 * np.std(Ranking["norm"]))])

        if snakemake.params.stable[i]: number_of_draws_stable = int(snakemake.params.stable[i])
        else: number_of_draws_stable = len(Ranking[Ranking["norm"] < (np.mean(Ranking["norm"]) - 2 * np.std(Ranking["norm"]))])
        

        print("\n", flush=True)
        print(f"shifted {i}: ", number_of_draws, flush=True)
        print(f"stable {i}: ", number_of_draws_stable, flush=True)


        # Observed by random
        population_size = len(Ranking)
        probability_of_success = success_in_the_population / population_size
        expected_random = number_of_draws * probability_of_success
        print("shifted group size expected at random: ", expected_random, flush=True)

        top_moving_norm_df = Ranking.sort_values(by="norm", ascending=False)[0:number_of_draws]
        top_moving_norm = set(top_moving_norm_df.index)
        top_success = len(top_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_head = round(hypergeom.sf(top_success, population_size, success_in_the_population, number_of_draws), 5,)

        less_moving_norm_df = Ranking.sort_values(by="norm", ascending=False)[-number_of_draws_stable:]
        less_moving_norm = set(less_moving_norm_df.index)
        less_success = len(less_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_tail = round(hypergeom.sf(less_success, population_size, success_in_the_population, number_of_draws_stable,),5,)

        print(
            f"% cancer related genes in the top {number_of_draws}:",
            top_success,
            f" p-value ({pvalue_of_head}) of this cancer relatedness with H0:=\"unrelated with cancer\" using hypergeom distribution for {i}"
        )
        print(
            f"% cancer related genes in the less {number_of_draws_stable}:",
            less_success,
            f" p-value ({pvalue_of_tail}) of this cancer relatedness with H0:=\"unrelated with cancer\" using hypergeom distribution for {i}"
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

        bar_width = 0.25
        x_indices = np.arange(len(categories))
        
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        ax.bar(x_indices - bar_width, percentage_most_moving, width=bar_width, label='shifted annotations', edgecolor=None)
        ax.bar(x_indices + bar_width, percentage_less_moving, width=bar_width, label='stable annotations', edgecolor=None)
        ax.bar(x_indices , percentage_by_random, width=bar_width, label='expected by random', edgecolor=None)
        for i in range(0, 13, 2):
            ax.axhline(y=i, color='gray', linestyle='-', linewidth=1)
        plt.ylabel('% cancer-related annotations', fontsize=24 , fontweight='bold')
        plt.xticks(x_indices, categories, fontsize=24)
        plt.yticks(fontsize=24)
        handles, labels = ax.get_legend_handles_labels()
        text = ax.text(-0.2,1.05, "", transform=ax.transAxes)
        lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
        plt.title("movement to cancer relatedness:", fontsize=24, fontdict=None, loc='center', pad=None)
        fig.savefig(f"{data_dir}/cancer_enrichments_2std.svg", dpi=500, bbox_extra_artists=(lgd,text), bbox_inches='tight')


f = open(snakemake.output.statistics, "w")
sys.stdout = f

Global_moving_in_the_space(snakemake.params.data_path, snakemake.params.annotation)

moving_in_the_space(snakemake.params.data_path, snakemake.params.annotation)

f.close()