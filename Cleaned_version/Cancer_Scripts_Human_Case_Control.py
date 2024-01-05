# General Packages:

from sklearn.metrics import pairwise_distances, silhouette_score
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import euclidean
from matplotlib_venn import venn2_unweighted
from scipy.stats import mannwhitneyu
from collections import Counter
from goatools import obo_parser
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import networkx as nx
import pandas as pd
import numpy as np
import itertools
import json
import os

# Functions:


def Create_Tissue_Specific_Network(
    Cancer_type, Normal_Tissue, Cell_type, save_path, plot=True
):
    # Load the files:

    # Load genes that are expressed in each tissue and cancer type:

    cancer_genes = list(
        pd.read_csv(f"{save_path}Cancer_{Cancer_type}_Genes.csv", dtype={0: str})["0"]
    )
    control_genes = list(
        pd.read_csv(f'{save_path}Control_{Normal_Tissue}_{Cell_type}_Genes.csv', dtype={0: str})["0"]
    )

    # Load the human PPI network:

    file_name_array = "Human_Biogrid_Adj"
    file_name_genelist = "Human_Biogrid_Genes"
    data_path = "./Data/"

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


def Prepare_GO_Gene_Matrix(
    GO_Gene_Matrix, Cancer_type, Normal_Tissue, Cell_type, save_path
):
    # Read the gene files:

    Cancer_genes = pd.read_table(
        save_path + "Cancer_" + Cancer_type + "_Genes.csv", header=0, dtype={0: str}
    )
    Cancer_genes_list = list(Cancer_genes["0"])

    Control_genes = pd.read_table(
        save_path + "Control_" + Normal_Tissue + "_" + Cell_type + "_Genes.csv",
        header=0,
        dtype={0: str},
    )
    Control_genes_list = list(Control_genes["0"])

    # Subser the GO embeddings:

    GO_Gene_Matrix_Cancer = GO_Gene_Matrix[GO_Gene_Matrix.index.isin(Cancer_genes_list)]
    GO_Gene_Matrix_Control = GO_Gene_Matrix[
        GO_Gene_Matrix.index.isin(Control_genes_list)
    ]

    # Reorder the Matrices:

    GO_Gene_Matrix_Cancer = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
    GO_Gene_Matrix_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]

    # Return the matrices:

    return (GO_Gene_Matrix_Cancer, GO_Gene_Matrix_Control)


def get_cosine_matrices_of_embeddings(
    dimension,
    Tissue,
    Sample,
    data_path,
    distance_measure="cosine",
    annnotation="BP_Back_Propagation",
    Gene_Matrix=False,
):
    if Gene_Matrix == False:
        path = (
            data_path
            + "_GO_Embeddings_"
            + annnotation
            + "_PPI_"
            + str(Tissue)
            + "_"
            + str(dimension)
            + "_PPMI_"
            + str(Sample)
            + ".csv"
        )
        embeddings = pd.read_csv(path, index_col=0).T
        embeddings[embeddings < 0] = 0
        embeddings_index = embeddings.index
        embeddings = np.array(embeddings)

        cosine_matrix_ppmi = pairwise_distances(embeddings, metric=distance_measure)
        cosine_matrix_ppmi_df = pd.DataFrame(
            cosine_matrix_ppmi, embeddings_index, embeddings_index
        )

        return cosine_matrix_ppmi_df

    else:
        path = (
            data_path
            + "_G_Matrix_"
            + str(dimension)
            + "_PPI_"
            + str(Tissue.title())
            + "_PPMI_"
            + str(Sample)
            + ".npy"
        )
        embeddings = np.load(path, allow_pickle=True)

        cosine_matrix_ppmi = pairwise_distances(embeddings, metric=distance_measure)

        return cosine_matrix_ppmi


def diff_Cosine_Distances_case_vs_control(
    k, Tissue, save_path, distance_measure="cosine", filtered=True, GO_filt_common=[]
):
    # Calculate the cosine distance:

    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Cancer", save_path, distance_measure
    )
    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Control", save_path, distance_measure
    )

    # Claculate the distances:

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        number_GO = 3

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]

        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_filt_common, GO_filt_common]
        cos_matrix_2 = cos_matrix_2.loc[GO_filt_common, GO_filt_common]

    # Get the distances:

    cos_matrix_1_distance = np.array(
        cos_matrix_1.values[np.triu_indices(len(cos_matrix_1), k=1)]
    )
    cos_matrix_2_distance = np.array(
        cos_matrix_2.values[np.triu_indices(len(cos_matrix_2), k=1)]
    )

    temp_eucl = euclidean(cos_matrix_1_distance, cos_matrix_2_distance)

    # At this point only returns the non normalized version.

    return temp_eucl


def Common_GO_Terms(Gene_GO_matrix_Cancer, Gene_GO_matrix_Control, number_GO):
    # Cancer:

    GO_Cancer = Gene_GO_matrix_Cancer.sum(axis=0)
    GO_terms_filtered_Cancer = set(GO_Cancer[GO_Cancer >= number_GO].index)

    # Control:

    GO_Control = Gene_GO_matrix_Control.sum(axis=0)
    GO_terms_filtered_Control = set(GO_Control[GO_Control >= number_GO].index)

    # Intersecions:

    Common_Annotations = GO_terms_filtered_Cancer.intersection(
        GO_terms_filtered_Control
    )

    # Return:

    return Common_Annotations


def Plot_Distances(
    k,
    Tissue,
    Sample,
    save_path,
    distance_measure="cosine",
    filtered=True,
    number_GO=3,
    GO_filt_common=[],
):
    # Calculate the cosine distance:

    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, Sample, save_path, distance_measure
    )

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]
        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_filt_common, GO_filt_common]

    # Print the average of the distances:

    print(np.mean(list(cos_matrix_1.values[np.triu_indices(len(cos_matrix_1), k=1)])))

    # Plot the distance:

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    cax = ax.matshow(cos_matrix_1.values, interpolation="nearest", cmap="jet")
    fig.colorbar(cax)

    # plt.savefig(save_path + "Color_map" + str(Tissue) + "_" + str(Sample) + ".png", dpi = 500)


def Compare_The_Similar_GO_Control_Infected(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered=True,
    GO_filt_common=[],
    number_similar=500,
    number_GO=3,
):
    # Calculate the distances:

    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Cancer", save_path, distance_measure
    )
    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Control", save_path, distance_measure
    )

    # Filter if needed:

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]

        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_filt_common, GO_filt_common]
        cos_matrix_2 = cos_matrix_2.loc[GO_filt_common, GO_filt_common]

    # Calculate the most similar and dissimilar:

    all_values_1 = (
        cos_matrix_1.mask(np.triu(np.ones(cos_matrix_1.shape)).astype(bool))
        .stack()
        .sort_values(ascending=True)
    )
    all_values_2 = (
        cos_matrix_2.mask(np.triu(np.ones(cos_matrix_2.shape)).astype(bool))
        .stack()
        .sort_values(ascending=True)
    )

    # Calculate the similar and disimilar:

    top_100_similar_values_Cancer = all_values_1[0:number_similar]
    top_100_similar_values_Control = all_values_2[0:number_similar]

    # For Cancer:

    # Plot Venn Diagram:

    control_infected_GO_Venn(
        top_100_similar_values_Control.index,
        top_100_similar_values_Cancer.index,
        Tissue,
        all_values_1,
        all_values_2,
    )

    # Return information:

    # return top_100_similar_values_Cancer, top_100_dissimilar_values_Cancer, top_100_similar_values_Control, top_100_dissimilar_values_Control


def control_infected_GO_Venn(
    list_cotrol, list_infected, Cancer_type, all_values_1, all_values_2
):
    cancer_genes = set(list_infected)
    control_genes = set(list_cotrol)
    common_genes = len(cancer_genes.intersection(control_genes))
    total_genes = len(cancer_genes.union(control_genes))
    unique_control = len(control_genes) - common_genes
    unique_cancer = len(cancer_genes) - common_genes

    # Plot statistics:

    common_genes_list = cancer_genes.intersection(control_genes)
    unique_control_list = control_genes - common_genes_list
    unique_cancer_list = cancer_genes - common_genes_list

    # Comparisons of distances in the corresponding space:

    print(
        "Common GO in the Control Space: ",
        str(all_values_2[list(common_genes_list)].mean()),
        "\n",
    )
    print(
        "Common GO in the Cancer  Space: ",
        str(all_values_1[list(common_genes_list)].mean()),
        "\n",
    )

    print(
        "Cancer GO in the Cancer  Space: ",
        str(all_values_1[list(unique_cancer_list)].mean()),
        "\n",
    )
    print(
        "Cancer GO in the Control  Space: ",
        str(all_values_2[list(unique_cancer_list)].mean()),
        "\n",
    )

    print(
        "Control GO in the Cancer  Space: ",
        str(all_values_1[list(unique_control_list)].mean()),
        "\n",
    )
    print(
        "Control GO in the Control  Space: ",
        str(all_values_2[list(unique_control_list)].mean()),
        "\n",
    )

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


def diff_Cosine_Distances_Optimal_Dim(
    k1,
    k2,
    Tissue,
    Sample,
    save_path,
    distance_measure="cosine",
    filtered=True,
    GO_filt_common=[],
    error=False,
    Gene_Embeddings=False,
):
    # Calculate the cosine distance:

    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k1, Tissue, Sample, save_path, distance_measure, Gene_Matrix=Gene_Embeddings
    )
    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k2, Tissue, Sample, save_path, distance_measure, Gene_Matrix=Gene_Embeddings
    )

    # Claculate the distances:

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        number_GO = 3

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]

        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_filt_common, GO_filt_common]
        cos_matrix_2 = cos_matrix_2.loc[GO_filt_common, GO_filt_common]

    if error == True:
        # Get the relative error:
        norm_R1 = np.linalg.norm(cos_matrix_1, ord="fro")
        norm_R2 = np.linalg.norm(cos_matrix_2, ord="fro")
        ErR1 = cos_matrix_1 - cos_matrix_2
        norm_erR1 = np.linalg.norm(ErR1, ord="fro")
        rel_erR1 = norm_erR1 / max(norm_R1, norm_R2)
        return rel_erR1
    else:
        # Get the distances:
        cos_matrix_1_distance = np.array(
            cos_matrix_1.values[np.triu_indices(len(cos_matrix_1), k=1)]
        )
        cos_matrix_2_distance = np.array(
            cos_matrix_2.values[np.triu_indices(len(cos_matrix_2), k=1)]
        )
        temp_eucl = euclidean(cos_matrix_1_distance, cos_matrix_2_distance)
        return temp_eucl


def Line_Plot_Distances(temp_dist, labels, error=False):
    # Calculate the distances

    temp_dist_euclidean = temp_dist

    # Plot them:

    plt.figure(figsize=(16, 6))
    plt.plot(labels, temp_dist_euclidean, marker="o")
    plt.xticks(labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Dimensions", fontsize=18)

    if error == False:
        plt.ylabel("Euclidean Distance", fontsize=18)
    else:
        plt.ylabel("Relative Error", fontsize=18)


def Euclidean_Heatmap(temp_dist, labels):
    pearson_corr_df = pd.DataFrame(temp_dist, index=labels, columns=labels)
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(
        pearson_corr_df, annot=False, ax=ax, cmap="jet"
    )  # Headmap change parametres if needed


def Compare_Cluster_Maps(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered=True,
    number_GO=3,
    GO_filt_common=[],
):
    # Calculate the cosine distance:

    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Cancer", save_path, distance_measure
    )
    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Control", save_path, distance_measure
    )

    # Filtering if needed

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]
        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_filt_common, GO_filt_common]
        cos_matrix_2 = cos_matrix_2.loc[GO_filt_common, GO_filt_common]

    # Create and plot the first cluster map (cluster map on Control and keep the order for cancer):

    linkage = hc.linkage(sp.distance.squareform(cos_matrix_2), method="average")
    control_clust = sns.clustermap(
        cos_matrix_2, cmap="jet", row_linkage=linkage, col_linkage=linkage
    )
    order = control_clust.dendrogram_row.reordered_ind

    order_cancer = cos_matrix_1.index[order]
    cancer_final = cos_matrix_1.reindex(order_cancer)
    cancer_final = cancer_final[cancer_final.index]

    order_control = cos_matrix_2.index[order]
    control_final = cos_matrix_2.reindex(order_control)
    control_final = control_final[control_final.index]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    cax = ax.matshow(
        cancer_final, interpolation="nearest", cmap="jet"
    )  # jet or rocket seismic
    fig.colorbar(cax)

    difference = control_final - cancer_final

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    cax = ax.matshow(difference, interpolation="nearest", cmap="jet")  # jet or rocket
    fig.colorbar(cax)

    # Difference between them:

    cos_matrix_3_diff = cos_matrix_2 - cos_matrix_1

    sns.clustermap(cos_matrix_3_diff, metric="euclidean", cmap="seismic")


def clusters_Analyses_Generate(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered=True,
    number_GO=3,
    GO_common_filt=[],
    top=5,
    threshold=0.80,
    Report=False,
    Summarize=False,
):
    # Calculate the cosine distance:

    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Cancer", save_path, distance_measure
    )
    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Control", save_path, distance_measure
    )

    # Filtering if needed

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]
        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_common_filt, GO_common_filt]
        cos_matrix_2 = cos_matrix_2.loc[GO_common_filt, GO_common_filt]

    # Difference between the cancer and the control:

    cos_matrix_3 = cos_matrix_1 - cos_matrix_2

    # Get the thresholds:

    upper_whisker, lower_whisker = Define_Threshold_whiskers(cos_matrix_3)

    # Get the GO terms using the thresholds:

    all_pairs = (
        cos_matrix_3.mask(np.triu(np.ones(cos_matrix_3.shape)).astype(bool))
        .stack()
        .sort_values(ascending=True)
    )

    # Here we have what happends if we go from control to cancer:

    # GO_Separating = all_pairs[all_pairs > upper_whisker]
    # GO_Clustering = all_pairs[all_pairs < lower_whisker]

    GO_Clustering = all_pairs[all_pairs >= threshold]
    GO_Separating = all_pairs[all_pairs <= -threshold]

    # Check the one that appears in more changes in each group:

    Separating_list = []
    for i in GO_Separating.index:
        Separating_list.append(i[0])
        Separating_list.append(i[1])

    Separating_list_count = Counter(Separating_list).most_common()

    Clustering_list = []
    for i in GO_Clustering.index:
        Clustering_list.append(i[0])
        Clustering_list.append(i[1])

    Clustering_list_count = Counter(Clustering_list).most_common()

    # Check which ones they are puting together:

    dic_cluster = {el[0]: [] for el in Clustering_list_count[:top]}
    dic_separet = {el[0]: [] for el in Separating_list_count[:top]}

    keys_clust = list(dic_cluster.keys())

    # GO terms that are closer in the space:

    for keys in keys_clust:
        for GO_clust in GO_Clustering.index:
            if GO_clust[0] == keys:
                dic_cluster[keys].append(GO_clust[1])
            elif GO_clust[1] == keys:
                dic_cluster[keys].append(GO_clust[0])
            else:
                continue

    # GO terms that are separating in the space:

    keys_separet = list(dic_separet.keys())

    for keys in keys_separet:
        for GO_clust in GO_Separating.index:
            if GO_clust[0] == keys:
                dic_separet[keys].append(GO_clust[1])
            elif GO_clust[1] == keys:
                dic_separet[keys].append(GO_clust[0])
            else:
                continue

    # Get the corresponding output depending on the user preferences:

    if Report == False:
        if Summarize == False:
            return (dic_cluster, dic_separet)
        else:
            dic_cluster = Summarize_Clusters_Semantic(dic_cluster)
            dic_separet = Summarize_Clusters_Semantic(dic_separet)
            return (dic_cluster, dic_separet)
    else:
        Analyze_GO_Clusters(
            cos_matrix_2,
            cos_matrix_1,
            k,
            Tissue,
            "Closer_Cancer",
            dic_cluster,
            save_path,
        )
        Analyze_GO_Clusters(
            cos_matrix_2,
            cos_matrix_1,
            k,
            Tissue,
            "Farther_Cancer",
            dic_separet,
            save_path,
        )
        print("Your report has beem computed into a txt file")


def Summarize_Clusters_Semantic(Similarity, clusters, all_comparison=False):
    # Result File:

    Result = {el: [] for el in range(len(clusters.keys()))}

    # For clusters

    for group in range(len(clusters.keys())):
        # GO terms in cluster:

        GO_list = [list(clusters.keys())[group]]
        GO_list.extend(clusters[GO_list[0]])

        # Terms that we have information in the semantic file:

        # Summary terms:

        GO_intersection = Similarity.index[Similarity.index.isin(GO_list)]
        Semantic_sub = Similarity.loc[GO_intersection, GO_intersection]

        # Calculat the mean of similarity taking out the 1:

        Semantic_sub["Semantic"] = round(
            (Semantic_sub.mean(axis=1) * len(Semantic_sub.mean(axis=1)) - 1)
            / (len(Semantic_sub.mean(axis=1) - 1)),
            2,
        )

        # Add the information:

        # if all comparisons is True the function gives back the mean semantic of each terms

        if all_comparison == True:
            Result[group] = {GO: [] for GO in Semantic_sub.index}

            for semantic_GO in Semantic_sub.index:
                Result[group][semantic_GO] = Semantic_sub.loc[semantic_GO, "Semantic"]

        # if not, the fucntion only return the GO terms with the maximoum semantic similarity between the rest:

        else:
            term = Semantic_sub[
                Semantic_sub.Semantic == Semantic_sub["Semantic"].max()
            ].index
            Result[group] = list(term)

    # Give back the information:

    return Result


def Get_Cluster_Score(clusters, cosine_distance_cancer, cosine_distance_control):
    # Per each cluster:

    mean_dist_cancer = []
    std_dist_cancer = []

    mean_dist_control = []
    std_dist_control = []

    for cluster in clusters.keys():
        # Get the list of GO terms to compare:

        GO_list = [cluster]
        GO_list.extend(clusters[cluster])

        # Calculate the mean distance:

        cosine_distance_cancer_filtered = cosine_distance_cancer.loc[GO_list, GO_list]
        cosine_distance_control_filtered = cosine_distance_control.loc[GO_list, GO_list]

        mean_dist_cancer.append(
            np.mean(
                np.array(
                    cosine_distance_cancer_filtered.values[
                        np.triu_indices(len(cosine_distance_cancer_filtered), k=1)
                    ]
                )
            )
        )
        std_dist_cancer.append(
            np.std(
                np.array(
                    cosine_distance_cancer_filtered.values[
                        np.triu_indices(len(cosine_distance_cancer_filtered), k=1)
                    ]
                )
            )
        )

        mean_dist_control.append(
            np.mean(
                np.array(
                    cosine_distance_control_filtered.values[
                        np.triu_indices(len(cosine_distance_control_filtered), k=1)
                    ]
                )
            )
        )
        std_dist_control.append(
            np.std(
                np.array(
                    cosine_distance_control_filtered.values[
                        np.triu_indices(len(cosine_distance_control_filtered), k=1)
                    ]
                )
            )
        )

    return (mean_dist_cancer, std_dist_cancer, mean_dist_control, std_dist_control)


def Analyze_GO_Clusters(
    cosine_distance_cancer,
    cosine_distance_control,
    k,
    Tissue,
    Type,
    clusters,
    save_path,
):
    # Load the ontology:

    GO_file = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/go-basic.obo"
    go = obo_parser.GODag(GO_file)

    # For Clusters:

    with open(
        save_path
        + "Description_File_"
        + str(Tissue)
        + "_"
        + str(k)
        + "_"
        + Type
        + ".txt",
        "a",
    ) as the_file:
        mean_can, std_can, mean_cont, std_cont = Get_Cluster_Score(
            clusters, cosine_distance_cancer, cosine_distance_control
        )

        statistics_control = 0

        for dimension in clusters.keys():
            # General Information about the cluster:

            if statistics_control <= len(clusters.keys()):
                the_file.write(
                    "Cancer_Mean: " + str(mean_can[statistics_control]) + "\n"
                )
                the_file.write("Cancer_Std: " + str(std_can[statistics_control]) + "\n")
                the_file.write(
                    "Control_Mean: " + str(mean_cont[statistics_control]) + "\n"
                )
                the_file.write(
                    "Control_Std: " + str(std_cont[statistics_control]) + "\n"
                )
                the_file.write("\n")

                statistics_control += 1

            # General information about the centroid:

            the_file.write(str(dimension) + "\t")
            the_file.write(str(dimension) + "\t")
            definition = go[dimension].name
            the_file.write(str(definition) + "\n")

            # Information about the GO terms that are clustered with the centroid:

            for GO in range(len(clusters[dimension])):
                the_file.write(str(dimension) + "\t")
                GO_term = clusters[dimension][GO]
                the_file.write(str(GO_term) + "\t")
                definition = go[GO_term].name
                the_file.write(str(definition) + "\n")

            # New line for the new cluster:

            the_file.write("\n")
            the_file.write("\n")

        # Close the file:

        the_file.close()


def Define_Threshold_whiskers(cos_matrix_3):
    cos_matrix_3_threshold = np.array(
        cos_matrix_3.values[np.triu_indices(len(cos_matrix_3), k=1)]
    )

    upper_quartile = np.percentile(cos_matrix_3_threshold, 75)
    lower_quartile = np.percentile(cos_matrix_3_threshold, 25)

    iqr = upper_quartile - lower_quartile
    upper_whisker = cos_matrix_3_threshold[
        cos_matrix_3_threshold <= upper_quartile + 1.5 * iqr
    ].max()
    lower_whisker = cos_matrix_3_threshold[
        cos_matrix_3_threshold >= lower_quartile - 1.5 * iqr
    ].min()

    return upper_whisker, lower_whisker


def Get_GO_Cosine_Differences(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered=True,
    number_GO=3,
    GO_common_filt=[],
):
    # Calculate the cosine distance:

    cos_matrix_2 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Cancer", save_path, distance_measure
    )
    cos_matrix_1 = get_cosine_matrices_of_embeddings(
        k, Tissue, "Control", save_path, distance_measure
    )

    # Filtering if needed

    if filtered == True:
        # Read the annotations that are being used:

        BP = json.load(
            open(
                "/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"
            )
        )

        # Fitler the GO terms:

        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [
            key
            for key, values in occurance_of_each_annotation_in_network.items()
            if values >= number_GO
        ]
        cos_matrix_1 = cos_matrix_1.loc[terms_filtering, terms_filtering]
        cos_matrix_2 = cos_matrix_2.loc[terms_filtering, terms_filtering]

    elif filtered == "Specific":
        cos_matrix_1 = cos_matrix_1.loc[GO_common_filt, GO_common_filt]
        cos_matrix_2 = cos_matrix_2.loc[GO_common_filt, GO_common_filt]

    cos_matrix_3 = cos_matrix_1 - cos_matrix_2

    return cos_matrix_3


def Get_GO_Difference_Ranking(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered="Specific",
    number_GO=3,
    GO_common_filt=[],
):
    # Get the difference:

    cos_matrix_3 = Get_GO_Cosine_Differences(
        k, Tissue, save_path, "cosine", filtered, 3, GO_common_filt
    )

    # Get the mean per row without taking into account the diagonal:

    cos_matrix_3["mean"] = cos_matrix_3.mean(axis=1)

    # Rank the values:

    cos_matrix_3 = cos_matrix_3.sort_values(by=["mean"], ascending=False)
    Rank_db = pd.DataFrame(cos_matrix_3["mean"])

    # Return the values:

    return Rank_db


def Define_Disrupted_Functional(Rank_db):
    # Define the threshold (two times the std):

    sns.boxplot(Rank_db)
    threshold = Rank_db.std() * 2

    # Generate the two groups:

    disrupeted_postive = Rank_db[(Rank_db > threshold)].dropna()
    disrupeted_negative = Rank_db[(Rank_db < -threshold)].dropna()

    disrupted = disrupeted_postive.append(disrupeted_negative)
    functional = Rank_db[~Rank_db.index.isin(disrupted.index)]

    # return the info:

    return (disrupted, functional)


def Taxons_Ranking(Rank_db, resolution=False):
    # Absolute values:

    Rank_db_abs = abs(Rank_db)

    # Load the Taxons

    Taxons_Dic = json.load(
        open("/media/sergio/sershiosdisk/Human/Axes/GO_Taxons_Whole_Annotations.json")
    )

    threshold = Rank_db.std() * 2

    disrupted = Rank_db_abs[(Rank_db_abs > threshold)].dropna()
    not_disrupted = Rank_db_abs[~Rank_db_abs.index.isin(disrupted.index)].dropna()

    Taxons_disrupted = {
        key: values for key, values in Taxons_Dic.items() if key in disrupted.index
    }
    Taxons_not = {
        key: values for key, values in Taxons_Dic.items() if key in not_disrupted.index
    }

    n_taxons_disrupted = [len(values) for values in Taxons_disrupted.values()]
    n_taxons_not = [len(values) for values in Taxons_not.values()]

    print(
        len(n_taxons_disrupted), np.mean(n_taxons_disrupted), np.std(n_taxons_disrupted)
    )
    print(len(n_taxons_not), np.mean(n_taxons_not), np.std(n_taxons_not))
    print(
        round(stats.ttest_ind(n_taxons_disrupted, n_taxons_not, equal_var=True).pvalue)
    )

    # For more resolution:

    # Get the bins

    if resolution == True:
        print("You choose more resolution\n")

        for threshold_it in [
            (1, 0.07),
            (0.07, 0.05),
            (0.05, 0.03),
            (0.03, 0.01),
            (0.01, 0),
        ]:
            itera = Rank_db_abs[
                (Rank_db_abs > threshold_it[1]) & (Rank_db_abs < threshold_it[0])
            ].dropna()
            Taxons_it = {
                key: values for key, values in Taxons_Dic.items() if key in itera.index
            }
            n_taxons = [len(values) for values in Taxons_it.values()]
            print(len(itera), np.mean(n_taxons), np.std(n_taxons))


def Leaf_Ranking(Rank_db, resolution=False):
    # Absolute values:

    Rank_db_abs = abs(Rank_db)

    # Load the Leafs:

    leaf_list = Define_Leaf_Annotations()
    leaf_list = Rank_db.index[Rank_db.index.isin(leaf_list)]

    # Separate groups:

    threshold = Rank_db.std() * 2

    disrupted = Rank_db_abs[(Rank_db_abs > threshold)].dropna()
    not_disrupted = Rank_db_abs[~Rank_db_abs.index.isin(disrupted.index)].dropna()

    # Calculate percentages:

    number_Leaf_dis = disrupted.index[disrupted.index.isin(leaf_list)]
    percentage_it_dis = len(set(number_Leaf_dis)) * 100 / len(disrupted)

    number_Leaf_not = not_disrupted.index[not_disrupted.index.isin(leaf_list)]
    percentage_it_not = len(set(number_Leaf_not)) * 100 / len(not_disrupted)

    print(round(percentage_it_dis, 1), round(percentage_it_not, 1))

    # Get the bins

    if resolution == True:
        for threshold_it in [
            (1, 0.07),
            (0.07, 0.05),
            (0.05, 0.03),
            (0.03, 0.01),
            (0.01, 0),
        ]:
            itera = Rank_db_abs[
                (Rank_db_abs > threshold_it[1]) & (Rank_db_abs < threshold_it[0])
            ].dropna()
            number_Leaf = itera.index[itera.index.isin(leaf_list)]
            percentage_it = len(set(number_Leaf)) * 100 / len(itera)

            print(len(itera), percentage_it)


def Define_Leaf_Annotations():
    save_path = "/media/sergio/sershiosdisk/Human/Cancer/"
    GO_Matrix = pd.read_csv(
        save_path + "_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv",
        index_col=0,
        dtype={0: str},
    )
    GO_Matrix.index = GO_Matrix.index.astype(str)

    # Prepare the Goatools:

    GO_file = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/go-basic.obo"
    go = obo_parser.GODag(GO_file)

    # Get the list of the GO terms:

    GO_terms = GO_Matrix.columns

    # Get the GO terms:

    GO_list = []

    for GO in GO_terms:
        try:
            if len(go[GO].children) == 0:
                GO_list.append(GO)
        except KeyError:
            continue

    # Give Back the list:

    return GO_list


def Extract_Clusters_From_Cluster_Map(
    k,
    Tissue,
    save_path,
    distance_measure="cosine",
    filtered=True,
    number_GO=3,
    GO_common_filt=[],
):
    # Get the difference between he cosine distances:

    cos_matrix_3 = Get_GO_Cosine_Differences(
        k, Tissue, save_path, distance_measure, filtered, number_GO, GO_common_filt
    )

    # Split the matrices in two parts (positive and negatives):

    # Positives:

    cos_matrix_3_positive = cos_matrix_3.copy()
    cos_matrix_3_positive[cos_matrix_3_positive < 0] = 0.001
    cos_matrix_3_positive = 1 - cos_matrix_3_positive
    np.fill_diagonal(cos_matrix_3_positive.values, 0)

    # Plot the cluster:

    linkage = hc.linkage(
        sp.distance.squareform(cos_matrix_3_positive), method="average"
    )
    sns.clustermap(
        cos_matrix_3_positive, cmap="rocket", row_linkage=linkage, col_linkage=linkage
    )

    # Get the number of clusters:

    optimal_dimension = Find_Optimal_Threshold(cos_matrix_3_positive)


def Find_Optimal_Threshold(cos_matrix_3):
    Optimal = {}

    for clusters in np.arange(100, 4000, 150):
        print(clusters)

        # Try with a concrete number of clusters:

        cluster = AgglomerativeClustering(
            n_clusters=clusters, affinity="precomputed", linkage="average"
        )

        # Get the clusters:

        Semantic_Cluster = cluster.fit_predict(cos_matrix_3)

        # Calculate the Silhouette score:

        Siluette = silhouette_score(
            cos_matrix_3, Semantic_Cluster, metric="precomputed"
        )

        # Save the score:

        Optimal[clusters] = round(Siluette, 2)

        print(round(Siluette, 2))

    # Calculate the optimal one:

    # First we keep the those dimensions that have a maximum value of siluette:

    max_value = max(Optimal.values())

    Optimal_dimensions = [k for k, v in Optimal.items() if v == max_value]

    # For the maximum Siluette values, return the one that corresponds to the lower number of clusters:

    return min(Optimal_dimensions)


def run_process(process):
    os.system(format(process))


def Run_Permutations():
    print("Running Permutations")

    command = "/home/sergio/Project_2/Scripts/python3 Terminal_Permutations_Cancer.py"
    run_process(command)

    print("Finished")


def Permutation_pvalues(
    Cancer_type, Normal_Tissue, Cell_type, k, Tissue, Sample, times=60000
):
    # Read the files.

    # GO Embeddings from the same k (at this point for back propagation).

    print("Cancer and Control Permutations, Cancer is not good be careful")

    # Gene and GO embeddings:

    save_path = "/media/sergio/sershiosdisk/Human/Cancer/"
    Gene_embeddings = np.load(
        save_path
        + "_G_Matrix_"
        + str(k)
        + "_PPI_"
        + str(Tissue)
        + "_PPMI_"
        + str(Sample)
        + ".npy",
        allow_pickle=True,
    )
    GO_Embeddings_pd = pd.read_csv(
        save_path
        + "_GO_Embeddings_BP_Back_Propagation_PPI_"
        + str(Tissue)
        + "_"
        + str(k)
        + "_PPMI_"
        + str(Sample)
        + ".csv",
        index_col=0,
        dtype={0: str},
    )
    # Matrix Gene/GO:

    GO_Matrix = pd.read_csv(
        save_path + "_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv",
        index_col=0,
        dtype={0: str},
    )
    GO_Matrix.index = GO_Matrix.index.astype(str)

    GO_Gene_Matrix_Cancer, GO_Gene_Matrix_Control = Prepare_GO_Gene_Matrix(
        GO_Matrix, Cancer_type, Normal_Tissue, Cell_type, save_path
    )

    if Sample == "Cancer":
        Gene_GO_matrix = GO_Gene_Matrix_Cancer
    else:
        Gene_GO_matrix = GO_Gene_Matrix_Control

    # Start the permutations:

    # Method to get the p-values using permutations:

    # Init the count matrix:

    Count_db = pd.DataFrame(
        0, columns=GO_Embeddings_pd.columns, index=GO_Embeddings_pd.index
    )

    # Keep the names:

    index_orig = Gene_GO_matrix.index
    columns_orig = Gene_GO_matrix.columns

    # Init the permutation test for N times:

    for rep in range(times):
        # Control de iteration:

        if rep % 100 == 0:
            print(rep)

        # Randomly suffle the Gene/GO matrix:

        Gene_GO_matrix_it = Gene_GO_matrix.sample(frac=1).reset_index(drop=True)
        Gene_GO_matrix_it.index = index_orig
        Gene_GO_matrix_it.columns = columns_orig

        # Do the Embeddings using the suffled Gene/GO matrix (direct calculation):

        gene_embeddings_db_inverse = pd.DataFrame(
            np.linalg.pinv(Gene_embeddings), columns=index_orig
        )
        GO_emebddings_it = gene_embeddings_db_inverse.dot(Gene_GO_matrix_it)

        # Compare the new scores vs the original scores:

        comparison = GO_emebddings_it >= GO_Embeddings_pd

        # Update the counts:

        Count_db = Count_db + comparison

    # Finished:

    print("the " + str(rep + 1) + " iterations are finished")

    # Calculate the p-values:

    Count_db_Final = (Count_db + 1) / (times + 1)

    # Save the matrix.

    Count_db_Final.to_csv(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "P_value_PPMI_PPI_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv",
        header=True,
        index=True,
    )


def Get_Adjusted_P_values_Filtered(
    k, Tissue, Sample, alpha=0.05, plot=True, cut=11, bins=10
):
    # Load the data (GO embeddings and the permutation matrix for the corresponding GO embedding matrix):

    save_path = "/media/sergio/sershiosdisk/Human/Cancer/"

    # Embeddings:

    GO_Embeddings_pd = pd.read_csv(
        save_path
        + "_GO_Embeddings_BP_Back_Propagation_PPI_"
        + str(Tissue)
        + "_"
        + str(k)
        + "_PPMI_"
        + str(Sample)
        + ".csv",
        index_col=0,
        dtype={0: str},
    )
    # P_values:

    Permutation_matrix_trans = pd.read_csv(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "P_value_PPMI_PPI_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv",
        index_col=0,
        dtype={0: str},
    )

    # Per each dimension we have performed a total test equals to the number of GO terms:

    NMTF_matrix = GO_Embeddings_pd.T

    # Prepare an empty matrix:

    P_adj_Matrix = pd.DataFrame(
        0,
        columns=Permutation_matrix_trans.columns,
        index=Permutation_matrix_trans.index,
    )

    # Get the information for the plots:

    if plot == True:
        dimensions_list = []
        go_scores_sign = []
        go_scores_list_nosign = []

    # We consider that per each GO term we did a total of N dimension test:

    for GO in Permutation_matrix_trans.columns:
        # Get the p-values:

        p_values_list = Permutation_matrix_trans[GO]

        # Correct the p-values:

        p_values_list_corrected = multipletests(
            p_values_list.values,
            alpha=alpha,
            method="fdr_bh",
            is_sorted=False,
            returnsorted=False,
        )

        # Save the p-values adjusted:

        P_adj_Matrix[GO] = p_values_list_corrected[1]

        # If the user want the plot:

        if plot == True:
            # Get the p-values that are significant:

            corrected_list = [i for i in p_values_list_corrected[1] if i <= alpha]

            # Get the information using these values for lists that at least have one significative value:

            if len(corrected_list) > 0:
                # Get the info and the scores:

                corrected_list_idx = [
                    number
                    for number, i in enumerate(p_values_list_corrected[1])
                    if i <= alpha
                ]
                GO_scores = [NMTF_matrix.loc[GO, i] for i in corrected_list_idx]
                GO_scores_nosign = [
                    NMTF_matrix.loc[GO, i]
                    for i in range(0, len(NMTF_matrix.columns))
                    if i not in corrected_list_idx
                ]

                dimensions_list.append(corrected_list_idx)
                go_scores_sign.append(GO_scores)
                go_scores_list_nosign.append(GO_scores_nosign)

            else:
                GO_scores_nosign = [
                    NMTF_matrix.loc[GO, i] for i in range(0, len(NMTF_matrix.columns))
                ]
                go_scores_list_nosign.append(GO_scores_nosign)

    # Do the statistical analyses

    if plot == True:
        General_Statistic_Analyses(
            Tissue,
            Sample,
            "P_Value_Corrected",
            go_scores_sign,
            go_scores_list_nosign,
            dimensions_list,
            alpha,
            k,
            cut,
            bins,
        )

    # Save the p-adjusted values:

    P_adj_Matrix.to_csv(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv",
        header=True,
        index=True,
    )


def General_Statistic_Analyses(
    Tissue,
    Sample,
    Name,
    go_scores_sign,
    go_scores_list_nosign,
    dimensions_list,
    k,
    alpha,
    cut,
    bins,
):
    # GO terms assigned:

    Total_GO_Terms_Assigned = len(go_scores_sign)

    # Mean Dimensions per GO term:

    Average_Dimensions = (
        sum([len(i) for i in go_scores_sign])
    ) / Total_GO_Terms_Assigned

    # Average GO terms per dimension:

    count_list = [item for sublist in dimensions_list for item in sublist]
    count_list = Counter(count_list)

    # Number of dimensions with a GO terms associated:

    # To control the titles of the plots:

    plt.figure(figsize=(10, 4))
    plt.title(
        "Dimensions with association: "
        + str(len(count_list.keys()))
        + "\n"
        + "Total_GO: "
        + str(Total_GO_Terms_Assigned)
        + "\n"
        + "Mean_Dimensions: "
        + str(round(Average_Dimensions))
    )
    plt.bar(count_list.keys(), count_list.values())
    plt.xlabel("# Dimensions")
    plt.ylabel("# GO terms")
    plt.savefig(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Dimensions_Association_P_Adj_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".png"
    )

    # Scores for those values that are significant or not:

    Scores_sign = [item for sublist in go_scores_sign for item in sublist]
    Scores_no_sign = [item for sublist in go_scores_list_nosign for item in sublist]

    # Filter the scores:

    Scores_sign_filt = [i for i in Scores_sign if i < cut]
    Scores_no_sign_filt = [i for i in Scores_no_sign if i < cut]

    print(
        "With the cut selected you lost: "
        + "\n"
        + "   - Significative: "
        + str(len(Scores_sign) - len(Scores_sign_filt))
        + "\n"
        + "   - Non Sign: ",
        len(Scores_no_sign) - len(Scores_no_sign_filt),
    )

    # Plot the distributions:

    print(
        "\n" + "The P-value of the Mann Whitney test is : ",
        mannwhitneyu(Scores_sign, Scores_no_sign).pvalue,
    )

    f = plt.figure(figsize=(12, 4))
    ax1 = f.add_subplot(121)
    # plt.subplot(1, 2, 1)
    plt.title("Non-Significant")
    plt.hist(Scores_no_sign_filt, bins, color="orange")
    plt.text(
        0.55,
        0.8,
        "Mean: %s" % round(np.mean(Scores_no_sign_filt), 2),
        transform=ax1.transAxes,
        fontsize=16,
    )
    plt.xlabel("# Scores")
    plt.ylabel("# Count")

    ax2 = f.add_subplot(122)
    plt.title("Significant " + str(Sample))
    plt.hist(Scores_sign_filt, bins)
    plt.text(
        0.55,
        0.8,
        "Mean: %s" % round(np.mean(Scores_sign_filt), 2),
        transform=ax2.transAxes,
        fontsize=16,
    )
    plt.xlabel("# Scores")
    plt.ylabel("# Count")
    plt.savefig(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Dimensions_Scores_P_value_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".png"
    )


def adj_pvalues_distr(Tissue, Sample, k):
    # Load Matrix

    pvalues_path = (
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv"
    )

    adjusted_pvalues_df = pd.read_csv(pvalues_path, index_col=0, header=0)

    # Select only p-values that are < 0.1:

    values = [
        item
        for sublist in adjusted_pvalues_df[adjusted_pvalues_df < 0.1].values
        for item in sublist
    ]
    cleanedList = [x for x in values if str(x) != "nan"]

    # Plot the distributions:

    plt.hist(cleanedList, bins=20)
    plt.xticks(
        [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
        [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
    )
    plt.title(
        "P-values distribution (" + str(Tissue) + " " + str(Sample) + " " + str(k) + ")"
    )
    plt.xlabel("p-value", fontsize=16)

    # Save the plot:

    plt.savefig(
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/P_Distributions_"
        + str(Sample)
        + "_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".png"
    )


def Compare_GO_Associations(
    k, Tissue, filtered="specific", number_GO=3, GO_common_filt=[]
):
    # Load the Corrected p-values:

    pvalues_path_Cancer = (
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_Cancer_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv"
    )
    pvalues_path_Control = (
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_Control_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv"
    )

    Cancer = pd.read_csv(pvalues_path_Cancer, index_col=0, dtype={0: str})
    Control = pd.read_csv(pvalues_path_Control, index_col=0, dtype={0: str})

    # Filter the associations:

    if filtered == "specific":
        Cancer = Cancer[GO_common_filt]
        Control = Control[GO_common_filt]

    # Associate the GO terms to dumensions:

    # Cancer:

    Associated_Cancer = []

    for GO in Cancer.columns:
        Cancer_it = Cancer[GO]
        if sum(Cancer_it <= 0.05) > 0:
            Associated_Cancer.append(GO)

    # Control:

    Associated_Control = []

    for GO in Cancer.columns:
        Control_it = Control[GO]
        if sum(Control_it <= 0.05) > 0:
            Associated_Control.append(GO)

    # Compare the GO terms:

    Cancer_Set = set(Associated_Cancer)
    Control_Set = set(Associated_Control)

    intersection = Cancer_Set.intersection(Control_Set)
    only_cancer = Cancer_Set - Control_Set
    only_control = Control_Set - Cancer_Set

    return (only_cancer, only_control, intersection)


def Permutations_Taxons_Specific(only_cancer, only_control, intersection):
    Taxons_Dic = json.load(
        open("/media/sergio/sershiosdisk/Human/Axes/GO_Taxons_Whole_Annotations.json")
    )

    Taxons_cancer = {
        key: values for key, values in Taxons_Dic.items() if key in only_cancer
    }
    Taxons_control = {
        key: values for key, values in Taxons_Dic.items() if key in only_control
    }
    Taxons_intersection = {
        key: values for key, values in Taxons_Dic.items() if key in intersection
    }

    Taxons_cancer = [len(values) for values in Taxons_cancer.values()]
    Taxons_control = [len(values) for values in Taxons_control.values()]
    Taxons_intersect = [len(values) for values in Taxons_intersection.values()]

    # Print results:

    print(len(Taxons_cancer), np.mean(Taxons_cancer), np.std(Taxons_cancer))
    print(len(Taxons_control), np.mean(Taxons_control), np.std(Taxons_control))
    print(len(Taxons_intersect), np.mean(Taxons_intersect), np.std(Taxons_intersect))


def Permutations_Leaf_Specific(only_cancer, only_control, intersection):
    # Load the Leafs:

    leaf_list = Define_Leaf_Annotations()

    # Calculate percentages:

    Leaf_cancer = sum(pd.DataFrame(only_cancer).isin(leaf_list)[0]) / len(only_cancer)
    Leaf_control = sum(pd.DataFrame(only_control).isin(leaf_list)[0]) / len(
        only_control
    )
    Leaf_intersection = sum(pd.DataFrame(intersection).isin(leaf_list)[0]) / len(
        intersection
    )

    # Print results:

    print("Cancer: ", str(Leaf_cancer))
    print("Control: ", str(Leaf_control))
    print("INtersectio: ", str(Leaf_intersection))


def Semantic_Similarity_Clustering(Similarity, GO_list):
    # Filter the Similarity based on the GO list:

    common_GO = Similarity.index[Similarity.index.isin(list(GO_list))]

    Similarity_filt = Similarity.loc[common_GO, common_GO]
    Distance = 1 - Similarity_filt

    # Calculate the Optimal number of clusters:

    n_clusters = Optimal_Clusters(Distance)

    # Apply the hierarchycal clustering:

    cluster = AgglomerativeClustering(
        n_clusters=n_clusters, affinity="precomputed", linkage="complete"
    )
    Semantic_Cluster = cluster.fit_predict(Distance)

    # Create the clusters:

    Dict_clust = {k: [] for k in list(set(Semantic_Cluster))}

    # Fill the clusters:

    i = 0
    for GO_terms in Similarity_filt.index:
        Dict_clust[Semantic_Cluster[i]].append(GO_terms)
        i += 1

    # Return the clusters:

    return Dict_clust


def Optimal_Clusters(distance_dim):
    Optimal = {}

    for clusters in range(2, len(distance_dim)):
        # Try with a concrete number of clusters:

        cluster = AgglomerativeClustering(
            n_clusters=clusters, affinity="precomputed", linkage="complete"
        )

        # Get the clusters:

        Semantic_Cluster = cluster.fit_predict(distance_dim)

        # Calculate the Silhouette score:

        Siluette = silhouette_score(
            distance_dim, Semantic_Cluster, metric="precomputed"
        )

        # Save the score:

        Optimal[clusters] = round(Siluette, 2)

    # Calculate the optimal one:

    # First we keep the those dimensions that have a maximum value of siluette:

    max_value = max(Optimal.values())

    Optimal_dimensions = [k for k, v in Optimal.items() if v == max_value]

    # For the maximum Siluette values, return the one that corresponds to the lower number of clusters:

    return min(Optimal_dimensions)


def File_Definitions_Permutations_Semantic_Clusters(
    k, Cluster_Summary, Cluster, Tissue, Sample, save_path
):
    # Load the ontology:

    GO_file = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/go-basic.obo"
    go = obo_parser.GODag(GO_file)

    # For Clusters:

    with open(
        save_path
        + "Permutation_Summary_file"
        + str(Tissue)
        + "_"
        + str(k)
        + "_"
        + Sample
        + ".txt",
        "a",
    ) as the_file:
        for dimension in list(Cluster_Summary.keys()):
            # Write the informations about the GO terms:

            for GO in Cluster_Summary[dimension]:
                # Write the cluster :

                the_file.write(str(dimension) + "\t")

                # Write number of GO in the cluster:

                the_file.write(str(len(Cluster[dimension])) + "\t")

                # Get the definition:

                definition = go[GO].name

                # Write the name:

                the_file.write(str(GO) + "\t")

                # Write the information:

                the_file.write(str(definition) + "\n")

            # Close the file:

        the_file.close()


def Get_Most_Important_Clusters(k, Tissue, Sample, summary_path, number):
    # Read the file:

    Clusters = pd.read_csv(
        summary_path
        + "Permutation_Summary_file"
        + str(Tissue)
        + "_"
        + str(k)
        + "_"
        + Sample
        + ".txt",
        sep="\t",
        header=None,
    )

    Clusters.columns = ["Cluster", "n_GO", "Representative", "Description"]

    # Order the file by number of GO terms in the cluster:

    Clusters = Clusters.sort_values("n_GO", ascending=False)

    # Subset the file for the most important ones:

    Clusters_select = Clusters.drop_duplicates("Cluster").Cluster[:number]
    Clusters_selection = Clusters[Clusters.Cluster.isin(Clusters_select)]

    return Clusters_selection


def Agreement_Permutation_Disrupted(
    k,
    Tissue,
    Sample,
    summary_path,
    rank_clusters,
    cancer_set,
    control_set,
    intersection_set,
    number,
):
    # Get the  disrupted clusters:

    Clusters = pd.read_csv(
        summary_path
        + "Permutation_Summary_file"
        + str(Tissue)
        + "_"
        + str(k)
        + "_"
        + Sample
        + ".txt",
        sep="\t",
        header=None,
    )

    Clusters.columns = ["Cluster", "n_GO", "Representative", "Description"]

    # Order the file by number of GO terms in the cluster:

    Clusters = Clusters.sort_values("n_GO", ascending=False)

    # Subset the file for the most important ones:

    Clusters_select = Clusters.drop_duplicates("Cluster").Cluster[:number]
    Clusters_selection = Clusters[Clusters.Cluster.isin(Clusters_select)]

    # Get the GO terms:

    in_cancer_list = []
    in_control_list = []
    in_intersection_list = []

    # Compare the top with the sets:

    for cluster in list(set(Clusters_selection.Cluster)):
        rank_clusters[cluster]

        GO_list = rank_clusters[cluster]

        # Compare with cancer set:

        in_cancer = set(GO_list).intersection(set(cancer_set))
        in_control = set(GO_list).intersection(set(control_set))
        in_intersection = set(GO_list).intersection(set(intersection_set))

        in_cancer_list.extend(list(in_cancer))
        in_control_list.extend(list(in_control))
        in_intersection_list.extend(list(in_intersection))

        # Print:

        print(
            "Cancer: ",
            str(len(in_cancer)),
            " Control: ",
            str(len(in_control)),
            " Intersection: ",
            str(len(in_intersection)),
            " Cluster: ",
            str(cluster),
        )

    # return the information:

    return (in_cancer_list, in_control_list, in_intersection_list)


def Compare_Distance_Sets(
    disrupted, in_cancer_list, in_control_list, in_intersection_list
):
    # Compare distances:

    distance_cancer = disrupted.loc[in_cancer_list].dropna().mean()
    distance_control = disrupted.loc[in_control_list].dropna().mean()
    distance_inter = disrupted.loc[in_intersection_list].dropna().mean()

    distance_cancer_std = disrupted.loc[in_cancer_list].dropna().std()
    distance_control_std = disrupted.loc[in_control_list].dropna().std()
    distance_inter_std = disrupted.loc[in_intersection_list].dropna().std()

    print(
        "Cancer: ",
        str(distance_cancer[0]),
        " ",
        str(distance_cancer_std[0]),
        "\n" "Control: ",
        str(distance_control[0]),
        " ",
        str(distance_control_std[0]),
        "\n",
        "Inter: ",
        str(distance_inter[0]),
        " ",
        str(distance_inter_std[0]),
    )


def Associate_GO_Axes(k, Tissue, filtered="specific", number_GO=3, GO_common_filt=[]):
    # Load the Corrected p-values:

    pvalues_path_Cancer = (
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_Cancer_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv"
    )
    pvalues_path_Control = (
        "/media/sergio/sershiosdisk/Human/Cancer/P_values/"
        + "Adjusted_P_value_Control_"
        + str(Tissue)
        + "_"
        + str(k)
        + ".csv"
    )

    Cancer = pd.read_csv(pvalues_path_Cancer, index_col=0, dtype={0: str})
    Control = pd.read_csv(pvalues_path_Control, index_col=0, dtype={0: str})

    # If the filtering is needed:

    if filtered == "specific":
        Cancer = Cancer[GO_common_filt]
        Control = Control[GO_common_filt]

    # Create the associations:

    Cluster_Cancer = []
    Cluster_Control = []

    # Cancer:

    for dimension in Cancer.index:
        list_Assoc = list(Cancer.columns[Cancer.loc[dimension] <= 0.05])

        if not list_Assoc:
            continue
        else:
            Cluster_Cancer.append(list_Assoc)

    # Control:

    for dimension in Control.index:
        list_Assoc = list(Control.columns[Control.loc[dimension] <= 0.05])

        if not list_Assoc:
            continue
        else:
            Cluster_Control.append(list_Assoc)

    # Return the clusters:

    return (Cluster_Cancer, Cluster_Control)


def Jaccar_index(cluster1, cluster2):
    intersection = len(list(set(cluster1).intersection(cluster2)))
    union = (len(cluster1) + len(cluster2)) - intersection
    return float(intersection) / union


def overlap_coefficient(cluster1, cluster2):
    intersection = len(list(set(cluster1).intersection(cluster2)))
    return float(intersection) / min(len(cluster1), len(cluster2))


def Get_Matrix_Jaccard_Overlap(Cluster_Cancer, Cluster_Control):
    # DataFrame for the results:

    Results_Jaccard = pd.DataFrame(
        0, index=np.arange(len(Cluster_Cancer)), columns=np.arange(len(Cluster_Control))
    )
    Results_Overlap = pd.DataFrame(
        0, index=np.arange(len(Cluster_Cancer)), columns=np.arange(len(Cluster_Control))
    )

    # Iterate per each cluster:

    for cancer_cluster in range(len(Cluster_Cancer)):
        cluster1 = Cluster_Cancer[cancer_cluster]

        for control_cluster in range(len(Cluster_Control)):
            cluster2 = Cluster_Control[control_cluster]

            jaccard = Jaccar_index(cluster1, cluster2)
            coeff = overlap_coefficient(cluster1, cluster2)

            Results_Jaccard.loc[cancer_cluster, control_cluster] = jaccard
            Results_Overlap.loc[cancer_cluster, control_cluster] = coeff

    # Result

    return (Results_Jaccard, Results_Overlap)


def Plot_Jaccard_Overlap(Results):
    # Plot the clustermap using a precomputed distance:

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    cax = ax.imshow(Results.values, interpolation="nearest", cmap="jet")
    fig.colorbar(cax)


def Get_Similar_Axes_Distribution(Results):
    result = Results.max()

    # Get the distributions:

    sns.boxplot(result)


def Get_Similar_Axes_Pair_list(Results):
    Pairs = pd.DataFrame(Results.idxmax()).reset_index()
    values = []

    # Get its value:

    for pair in range(len(Pairs)):
        position = Pairs.loc[pair]
        values.append(Results.loc[position[0], position["index"]])

    # Create the DataFrame:

    Pairs["Value"] = values
    Pairs.columns = ["Control", "Cancer", "Value"]

    # Order the information:

    Pairs = Pairs.sort_values("Value", ascending=False)

    # Return the information:

    return Pairs


def Divide_groups_Axes(Pairs, Cluster_Cancer, Cluster_Control):
    # Distribution of the pairs to binarize them:

    Thresholds = Pairs.Value.quantile([0.25, 0.75])

    not_moving_threshold = round(Thresholds.loc[0.75], 2)
    moving_threshold = round(Thresholds.loc[0.25], 2)

    Pairs["Len Control"] = [len(Cluster_Control[i]) for i in Pairs.Control]
    Pairs["Len Cancer"] = [len(Cluster_Cancer[i]) for i in Pairs.Cancer]

    # Get the three groups:

    stable = Pairs[Pairs["Value"] >= not_moving_threshold]
    moving = Pairs[Pairs["Value"] <= moving_threshold]
    between = Pairs[
        (Pairs["Value"] >= moving_threshold) & (Pairs["Value"] <= not_moving_threshold)
    ]

    # Information of the three sets:

    print(
        "Stable Axes: ",
        str(len(stable)),
        "\n",
        "Moving Axes: ",
        str(len(moving)),
        "\n",
        "In between: ",
        str(len(between)),
    )

    # Information of extreme cases:

    print("Equal Axes: ", str(len(stable[stable.Value == 1])))
    print("Different Axes: ", str(len(moving[moving.Value == 0])))

    # Return the options:

    return (stable, moving, between)


def Get_Statistics_group_Axes(equals_list, only_cancer, only_control, dic):
    # For axes that are changing:

    # For each cluster:

    intersection_n = []
    cancer_n = []
    control_n = []

    for key in dic.keys():
        intersection_n.append(len(dic[key]["Intersect"]))
        cancer_n.append(len(dic[key]["Cancer"]))
        control_n.append(len(dic[key]["Control"]))

    # Print results:

    print("Axes Changing: ", len(dic.keys()), "\n")
    print("Statistics\n")

    print(
        "Intersect mean: ",
        str(np.mean(intersection_n)),
        " std: ",
        str(np.std(intersection_n)),
    )
    print("Cancer mean: ", str(np.mean(cancer_n)), " std: ", str(np.std(cancer_n)))
    print("Control mean: ", str(np.mean(control_n)), " std: ", str(np.std(control_n)))

    # Statistics:

    print(
        "P_value Cancer Vs Control: ", str(stats.ttest_ind(cancer_n, control_n).pvalue)
    )

    # For equal lists:

    print("Number of Axes that are the same: ", str(len(equals_list)), "\n")

    # For only control (not repeated):

    only_control.sort()
    only_control_not_repeated = list(k for k, _ in itertools.groupby(only_control))
    print("Number of Axes only in control: ", str(len(only_control_not_repeated)), "\n")

    # For only in cancer (not repeated):

    only_cancer.sort()
    only_cancer_not_repeated = list(k for k, _ in itertools.groupby(only_cancer))
    print("Number of Axes only in cancer: ", str(len(only_cancer_not_repeated)), "\n")


def Summarize_The_Changing_Axes(dic, Similarity):
    # For equal_list:

    dic_equals = {k: [] for k in range(len(dic.keys()))}

    # For each dimension:

    for dimension in dic.keys():
        intersection_go = dic[dimension]["Intersect"]
        cancer_go = dic[dimension]["Cancer"]
        control_go = dic[dimension]["Control"]

        dic_equals[dimension] = {}

        # For intersection:

        if len(intersection_go) > 2:
            dic_equals[dimension]["Intersect"] = Semantic_Similarity_Clustering_Axes(
                Similarity, intersection_go
            )
        else:
            dic_equals[dimension]["Intersect"] = intersection_go

        # For cancer:

        if len(cancer_go) > 2:
            dic_equals[dimension]["Cancer"] = Semantic_Similarity_Clustering_Axes(
                Similarity, cancer_go
            )
        else:
            dic_equals[dimension]["Cancer"] = cancer_go

        # For control:

        if len(control_go) > 2:
            dic_equals[dimension]["Control"] = Semantic_Similarity_Clustering_Axes(
                Similarity, control_go
            )
        else:
            dic_equals[dimension]["Control"] = control_go

    # Return the clusters:

    return dic_equals


def Semantic_Similarity_Clustering_Axes(Similarity, GO_list):
    # Filter the Similarity based on the GO list:

    common_GO = Similarity.index[Similarity.index.isin(list(GO_list))]

    # Control the GO that we dont have information:

    Common_set = set(common_GO)
    GO_list_set = set(GO_list)
    no_info = GO_list_set - Common_set

    # If there is more than 2 that we can cluster:

    if len(Common_set) > 2:
        Similarity_filt = Similarity.loc[common_GO, common_GO]
        Distance = 1 - Similarity_filt

        # Calculate the Optimal number of clusters:

        n_clusters = Optimal_Clusters(Distance)

        # Apply the hierarchycal clustering:

        cluster = AgglomerativeClustering(
            n_clusters=n_clusters, affinity="precomputed", linkage="complete"
        )
        Semantic_Cluster = cluster.fit_predict(Distance)

        # Create the clusters:

        Dict_clust = {k: [] for k in list(set(Semantic_Cluster))}

        # Fill the clusters:

        i = 0
        for GO_terms in Similarity_filt.index:
            Dict_clust[Semantic_Cluster[i]].append(GO_terms)
            i += 1

        # Add those that we dont have information:

        if len(no_info) > 0:
            Dict_clust["No Info"] = {}
            Dict_clust["No Info"] = list(no_info)

        # Return the clusters:

        return Dict_clust

    # If we can not perform clusters then return all

    else:
        return GO_list


def Get_The_Representative_Axes(Similarity, dic_equals):
    dic_summary = {k: [] for k in range(len(dic_equals.keys()))}

    for dimensions in dic_equals.keys():
        dimension_clusters_intersection = dic_equals[dimensions]["Intersect"]
        dimension_clusters_cancer = dic_equals[dimensions]["Cancer"]
        dimension_clusters_control = dic_equals[dimensions]["Control"]

        dic_summary[dimensions] = {}

        # Intersection:

        if type(dimension_clusters_intersection) == dict:
            dic_summary[dimensions][
                "Intersect"
            ] = Summarize_Clusters_Semantic_Changing_Axes(
                dimension_clusters_intersection, Similarity
            )
        else:
            dic_summary[dimensions]["Intersect"] = dimension_clusters_intersection

        # Cancer:

        if type(dimension_clusters_cancer) == dict:
            dic_summary[dimensions][
                "Cancer"
            ] = Summarize_Clusters_Semantic_Changing_Axes(
                dimension_clusters_cancer, Similarity
            )
        else:
            dic_summary[dimensions]["Cancer"] = dimension_clusters_cancer

        # Control:

        if type(dimension_clusters_control) == dict:
            dic_summary[dimensions][
                "Control"
            ] = Summarize_Clusters_Semantic_Changing_Axes(
                dimension_clusters_control, Similarity
            )
        else:
            dic_summary[dimensions]["Control"] = dimension_clusters_control

    # Return the information:

    return dic_summary


def Summarize_Clusters_Semantic_Changing_Axes(dictionary, Similarity):
    # Create the result:

    Result = {el: [] for el in dictionary.keys()}

    # For each cluster:

    for cluster in dictionary.keys():
        if cluster != "No Info":
            GO_list = dictionary[cluster]
            GO_intersection = Similarity.index[Similarity.index.isin(GO_list)]
            Semantic_sub = Similarity.loc[GO_intersection, GO_intersection]

            Semantic_sub["Semantic"] = round(
                (Semantic_sub.mean(axis=1) * len(Semantic_sub.mean(axis=1)) - 1)
                / (len(Semantic_sub.mean(axis=1) - 1)),
                2,
            )
            term = Semantic_sub[
                Semantic_sub.Semantic == Semantic_sub["Semantic"].max()
            ].index

            Result[cluster] = list(term)

        else:
            GO_list = dictionary[cluster]
            Result[cluster] = GO_list

    return Result


def write_the_Axes_Changes_information(k, Tissue, dic_total, dic_summ, save_path):
    # Load the ontology:

    GO_file = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/go-basic.obo"
    go = obo_parser.GODag(GO_file)

    # For Clusters:

    with open(
        save_path
        + "Permutation_Summary_file_Axes_Changing"
        + str(Tissue)
        + "_"
        + str(k)
        + ".txt",
        "a",
    ) as the_file:
        # For axes:

        for axes in dic_summ.keys():
            # For the group:

            for comparison in ["Cancer", "Control", "Intersect"]:
                # If something was clustered:

                if type(dic_summ[axes][comparison]) == dict:
                    for cluster in dic_summ[axes][comparison].keys():
                        # Get more information:

                        GO_list = dic_summ[axes][comparison][cluster]
                        number = len(dic_total[axes][comparison][cluster])

                        for GO in GO_list:
                            the_file.write(str(axes) + "\t")
                            the_file.write(str(comparison) + "\t")
                            the_file.write(str(cluster) + "\t")

                            # write the number of GO in the cluster:

                            the_file.write(str(number) + "\t")

                            # write the name of the GO:

                            the_file.write(str(GO) + "\t")

                            # write the description:

                            definition = go[GO].name
                            the_file.write(str(definition) + "\n")

                else:
                    GO_list = dic_summ[axes][comparison]
                    number = len(dic_total[axes][comparison])

                    for GO in GO_list:
                        the_file.write(str(axes) + "\t")
                        the_file.write(str(comparison) + "\t")
                        the_file.write(str("NO_cluster") + "\t")
                        the_file.write(str(number) + "\t")

                        # write the name of the GO:

                        the_file.write(str(GO) + "\t")

                        # write the description:

                        definition = go[GO].name
                        the_file.write(str(definition) + "\n")

    the_file.close()


def Get_Functional_Domains(Pairs, Cluster_Cancer, Cluster_Control, Similarity):
    # Get axes that could have different Functional domains:

    changes = Pairs[(Pairs["Value"] != 1) & (Pairs["Value"] != 0)]
    changes = changes.reset_index(drop=True)

    # Dictinaries with the results:

    cancer_domains = {k: [] for k in changes["Cancer"]}
    control_domains = {k: [] for k in changes["Control"]}

    # Obtain the functional domains per each axe:

    for axe in range(len(changes)):
        # Take the axe:

        Cancer_axe = changes.loc[axe, "Cancer"]
        Control_axe = changes.loc[axe, "Control"]

        # GO terms in the axes:

        Cancer_GO = Cluster_Cancer[Cancer_axe]
        Control_GO = Cluster_Control[Control_axe]

        # Get the functional domains:

        cancer_domains[Cancer_axe] = {}
        control_domains[Control_axe] = {}

        if len(Cancer_GO) > 3:
            cancer_domains[Cancer_axe] = Semantic_Similarity_Clustering(
                Similarity, Cancer_GO
            )
        else:
            cancer_domains[Cancer_axe] = Cancer_GO
        if len(Control_GO) > 3:
            control_domains[Control_axe] = Semantic_Similarity_Clustering(
                Similarity, Control_GO
            )
        else:
            cancer_domains[Control_axe] = Control_GO
