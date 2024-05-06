from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom
import operator
from collections import Counter
from functools import reduce
from sklearn.cluster import KMeans
import numpy as np
import json
import pandas as pd


def flatten(l):
    """flattens a nested list"""
    return reduce(operator.concat, l)


def enrichment_analysis(
    clusters,
    annotation_mapping,
    dim,
    data_path,
    alpha=0.05,
    enrichment="BP",
    network="PPI",
    matrix="Adj",
    enriched="statistics",
    get_p_values=True,
):
    """
    Input: clusters --> list of lists
           annotation_mapping --> dictionary (key --> genes, values --> annotation)
    """

    """
    N: Number of annotated genes in the cluster.
    X: Number of genes annotated with the given annotation in the cluster.
    M: Number of annotated genes in the dataset.
    K: Number of genes annotated with the given annotation in the dataset/network.
    """

    nodes_in_the_ppi = flatten(clusters)
    # nodes_in_the_ppi = [str(i) for i in nodes_in_the_ppi]
    common_nodes = set(nodes_in_the_ppi).intersection(set(annotation_mapping.keys()))
    filtered_mapping = {i: annotation_mapping[i] for i in common_nodes}

    M = len(common_nodes)
    annotation_list = flatten(filtered_mapping.values())
    occurance_of_each_annotation_in_network = Counter(annotation_list)

    enriched_term_in_each_cluster = {}
    enriched_genes_in_network = []

    save_path = (
        data_path
        + "Enrichments_"
        + str(enrichment)
        + "_"
        + str(network)
        + "_"
        + str(matrix)
        + "_"
        + str(dim)
    )

    # For the p_values:

    if get_p_values == True:
        p_values_df = pd.DataFrame(
            np.nan, index=set(annotation_list), columns=list(range(dim))
        )

    # Enrichment:

    for number, cluster in enumerate(clusters):
        N = len(cluster)

        if N < 5:
            print("Length of cluster ", N)
            enriched_term_in_each_cluster[number] = []
        else:
            cluster_annotation = []
            for element in cluster:
                try:
                    cluster_annotation.extend(filtered_mapping[element])
                except KeyError:
                    N -= 1
            occurance_of_each_annotation_in_cluster = Counter(cluster_annotation)

            p_vals = []
            unique_annotations_in_cluster = list(set(cluster_annotation))
            terms_enriched = []
            for i in unique_annotations_in_cluster:
                try:
                    X = occurance_of_each_annotation_in_cluster[i]
                    K = occurance_of_each_annotation_in_network[i]
                    p_vals.append(hypergeom.sf(X - 1, M, K, N))
                    terms_enriched.append(i)
                except KeyError:
                    pass

            if len(p_vals) == 0:
                enriched_term_in_each_cluster[number] = []
            else:
                corrected_pvalues = multipletests(
                    p_vals,
                    alpha=alpha,
                    method="fdr_bh",
                    is_sorted=False,
                    returnsorted=False,
                )[1]
                enriched_terms = [
                    terms_enriched[position]
                    for position, i in enumerate(corrected_pvalues)
                    if i <= alpha
                ]

                if get_p_values == True:
                    Get_P_Values(p_values_df, corrected_pvalues, number, terms_enriched)

                if "GO:0008150" in enriched_terms:
                    enriched_terms.remove("GO:0008150")
                    print(
                        "Cluster %s enriched in the root of the ontology and %s more terms"
                        % (number, len(enriched_terms))
                    )

                enriched_term_in_each_cluster[number] = list(set(enriched_terms))

                for node in cluster:
                    try:
                        if (
                            len(
                                set(filtered_mapping[node]).intersection(
                                    set(enriched_terms)
                                )
                            )
                            != 0
                        ):
                            enriched_genes_in_network.append(node)
                    except KeyError:
                        pass

    if enriched == "terms":
        # enriched_terms = set(flatten(enriched_term_in_each_cluster.values()))
        enriched_go_terms = set(
            [
                go_term
                for cluster in enriched_term_in_each_cluster.values()
                for go_term in cluster
            ]
        )

        return enriched_term_in_each_cluster.values()

    elif enriched == "terms in clusters":
        # Save it:
        save_path_1 = save_path + "Per_Cluster"
        with open(save_path_1 + ".json", "w") as fp:
            json.dump(enriched_term_in_each_cluster, fp)

        if get_p_values == True:
            save_path_P_value = save_path + "P_values.csv"
            p_values_df.to_csv(save_path_P_value)

    elif enriched == "genes":
        # Save it:
        save_path_2 = save_path + "Genes_Network"
        with open(save_path_2 + ".json", "w") as fp:
            json.dump(enriched_term_in_each_cluster, fp)
    elif enriched == "both":
        return enriched_terms, enriched_genes_in_network

    else:  # Statistics:
        total_anotations_in_network = len(occurance_of_each_annotation_in_network)

        annotated_clusters = [
            i for i in enriched_term_in_each_cluster.values() if len(i) != 0
        ]
        percentage_of_enriched_clusters = round(
            len(annotated_clusters) / float(len(clusters)) * 100, 2
        )
        percentage_of_enriched_genes = round(
            len(enriched_genes_in_network) / float(len(filtered_mapping)) * 100, 2
        )

        try:
            percentage_of_enriched_annotations = round(
                len(list(set(flatten(annotated_clusters))))
                / float(total_anotations_in_network)
                * 100,
                2,
            )
        except TypeError:
            percentage_of_enriched_annotations = 0.0
            print("No annotations are enriched in the network")

        # Save tables: #

        table_Statistics = pd.DataFrame(
            {
                "Enriched_Clusters": [percentage_of_enriched_clusters],
                "%GO": [percentage_of_enriched_annotations],
                "%Genes": [percentage_of_enriched_genes],
            }
        )

        save_path_3 = save_path + "_Statistics.csv"
        table_Statistics.to_csv(save_path_3)

        # return percentage_of_enriched_clusters, percentage_of_enriched_annotations ,percentage_of_enriched_genes


def kmeans_clusters(matrix, labels, number_of_clusters=95):
    kmeans = KMeans(n_clusters=number_of_clusters, random_state=0).fit(matrix)
    kmeans_labels = kmeans.labels_
    clusters = [[] for i in range(number_of_clusters)]
    n = len(kmeans_labels)
    for i in range(n):
        clusters[kmeans_labels[i]].append(labels[i])
    return clusters


def embedding_space_enrichments(embedding_space, labels, annotation_mapping):
    clusters_embedding_space = kmeans_clusters(embedding_space, labels)
    (
        percentage_of_enriched_clusters,
        percentage_of_enriched_annotations,
        percentage_of_enriched_genes,
    ) = enrichment_analysis(
        clusters_embedding_space, annotation_mapping, alpha=0.05, enriched="statistics"
    )
    return (
        percentage_of_enriched_clusters,
        percentage_of_enriched_annotations,
        percentage_of_enriched_genes,
    )


# New Functions:


def Get_Hard_Clusters(M, nodes):
    n, k = np.shape(M)
    Clusters = [[] for i in range(k)]
    for i in range(n):
        idx = np.argmax(M[i])
        Clusters[idx].append(nodes[i])
    return Clusters


"""
def Get_P_Values(p_values_dic, corrected_pvalues, number,terms_enriched):
    p_values_dic[str(number)] = {}
    GO_terms = [terms_enriched[position] for position,i in enumerate(corrected_pvalues)]
    for GO in range(len(GO_terms)):
        p_values_dic[str(number)][GO_terms[GO]] = corrected_pvalues[GO] 
"""


def Get_P_Values(p_values_df, corrected_pvalues, number, terms_enriched):
    GO_terms = [
        terms_enriched[position] for position, i in enumerate(corrected_pvalues)
    ]
    for GO in range(len(GO_terms)):
        p_values_df.loc[GO_terms[GO]][number] = corrected_pvalues[GO]
