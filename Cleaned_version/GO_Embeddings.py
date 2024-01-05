# Script with functions to get the GO terms in the genes space.

# Libraries:

import pandas as pd
import numpy as np
import json
import NMF_Implementation
import NMF_Implementation_Human
from collections import Counter
from sklearn.metrics import roc_auc_score
from sklearn import metrics
from scipy.stats import mannwhitneyu


def Gene_GO_Matrix(Gene_Annotations, gene_names, save_path, GO="BP", network="PPI"):
    """
    Imputs:

        Gene_Annotations: Dictionary of Genes as index and GO terms per gene
        gene_names : list of genes with the same order as the genes embedding Matrix

    """
    # Matrix Genes/GO:

    go_terms = list(
        set([item for sublist in Gene_Annotations.values() for item in sublist])
    )
    array_np = np.zeros((len(gene_names), len(go_terms)))

    # Binarize the matrix:

    # Not all the genes have annotation:

    filtered_gene_annotations = [i for i in Gene_Annotations.keys() if i in gene_names]

    for key in filtered_gene_annotations:
        for value in Gene_Annotations[key]:
            array_np[gene_names.index(key)][go_terms.index(value)] = 1

    # Generate the Binary pandas:

    Binary_pd = pd.DataFrame(array_np, columns=go_terms, index=gene_names)
    Binary_pd = Binary_pd.astype(int)
    Binary_pd = Binary_pd[sorted(Binary_pd.columns)]

    # Save Matrix Genes/GO

    save_path_1 = save_path + "_Matrix_Genes_GO_" + str(GO) + "_" + str(network)
    Binary_pd.to_csv(save_path_1 + ".csv", index=True)


def GO_Embeddings(
    Matrix_Gene_GO,
    gene_embeddings,
    save_path,
    GO="BP",
    network="PPI",
    matrix="PPMI",
    ortogonal=False,
    init="SVD",
    repetitions=200,
    save=True,
):
    """
    Imputs:

        Matrix_Gene_GO : DataFrame with genes and GO terms (0 : Not annotated, 1 : Annotated)
        gene_embeddings: Array with the embeddings of the genes.

    """

    # Matrix with gene embeddings (if the previous function was used to create the gene/GO the order of the genes is the same):

    gene_embeddings_db = pd.DataFrame(gene_embeddings)
    gene_embeddings_db.index = Matrix_Gene_GO.index

    # Calculate the GO embeddings:

    # If is ortogonal V transpose is equal to V^-1 then we can use it we want to get the GO embeddings.
    # This is equals to (V^-1) * R = W

    if ortogonal == True:
        gene_embeddings_db_inverse = pd.DataFrame(
            np.linalg.pinv(gene_embeddings_db.values),
            gene_embeddings_db.columns,
            gene_embeddings_db.index,
        )
        GO_embeddings = gene_embeddings_db_inverse.dot(Matrix_Gene_GO)

        if save == True:
            save_path_1 = (
                save_path
                + "_GO_Embeddings_"
                + str(GO)
                + "_"
                + str(network)
                + "_"
                + str(len(gene_embeddings_db.columns))
                + "_"
                + str(matrix)
                + ".csv"
            )
            GO_embeddings.to_csv(save_path_1)
        else:
            return GO_embeddings

    # If not we have to use the NMF.
    # In this case by default we used 500 reps, fix_U and initialitation of SVD.

    elif ortogonal == False:
        GO_embeddings = NMF_Implementation.NMF(
            Matrix_Gene_GO,
            len(gene_embeddings_db.columns),
            repetitions,
            init=init,
            fix_U=True,
            U=gene_embeddings,
        )

        save_path_1 = (
            save_path
            + "_GO_Embeddings_"
            + str(GO)
            + "_"
            + str(network)
            + "_"
            + str(len(gene_embeddings_db.columns))
            + "_"
            + str(matrix)
            + ".csv"
        )
        GO_embeddings.to_csv(save_path_1)


def GO_Embeddings_Human(
    Matrix_Gene_GO,
    gene_embeddings,
    save_path,
    GO="BP",
    network="PPI",
    matrix="PPMI",
    ortogonal=False,
    init="SVD",
    repetitions=200,
    save=True,
):
    """
    Imputs:

        Matrix_Gene_GO : DataFrame with genes and GO terms (0 : Not annotated, 1 : Annotated)
        gene_embeddings: Array with the embeddings of the genes.

    """

    # Matrix with gene embeddings (if the previous function was used to create the gene/GO the order of the genes is the same):

    gene_embeddings_db = pd.DataFrame(gene_embeddings)
    gene_embeddings_db.index = Matrix_Gene_GO.index

    # Calculate the GO embeddings:

    # If is ortogonal V transpose is equal to V^-1 then we can use it we want to get the GO embeddings.
    # This is equals to (V^-1) * R = W

    if ortogonal == True:
        gene_embeddings_db_inverse = pd.DataFrame(
            np.linalg.pinv(gene_embeddings_db.values),
            gene_embeddings_db.columns,
            gene_embeddings_db.index,
        )
        GO_embeddings = gene_embeddings_db_inverse.dot(Matrix_Gene_GO)

        if save == True:
            save_path_1 = (
                save_path
                + "_GO_Embeddings_"
                + str(GO)
                + "_"
                + str(network)
                + "_"
                + str(len(gene_embeddings_db.columns))
                + "_"
                + str(matrix)
                + ".csv"
            )
            GO_embeddings.to_csv(save_path_1)
        else:
            return GO_embeddings

    # If not we have to use the NMF.
    # In this case by default we used 500 reps, fix_U and initialitation of SVD.

    elif ortogonal == False:
        GO_embeddings = NMF_Implementation_Human.NMF(
            Matrix_Gene_GO,
            len(gene_embeddings_db.columns),
            repetitions,
            init=init,
            fix_U=True,
            U=gene_embeddings,
        )

        save_path_1 = (
            save_path
            + "_GO_Embeddings_"
            + str(GO)
            + "_"
            + str(network)
            + "_"
            + str(len(gene_embeddings_db.columns))
            + "_"
            + str(matrix)
            + ".csv"
        )
        GO_embeddings.to_csv(save_path_1)


def Get_Hard_Clusters(M, nodes):
    n, k = np.shape(M)
    Clusters = [[] for i in range(k)]
    for i in range(n):
        idx = np.argmax(M[i])
        Clusters[idx].append(nodes[i])
    return Clusters


def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    if union == 0:
        return 0
    else:
        return float(intersection) / union


def Test_GO_Embeddings_Jaccard(
    data_path,
    dimensions,
    network="PPI",
    matrix="PPMI",
    enrichment="BP",
    intersection=False,
    soft_clustering=False,
):
    """
    Imputs:
            data_path    : path to the data folder.
            dimensions   : number of dimensions.
            network      : network to analyse.
            matrix       : matrix that was factorized.
            enrichments  : annotations.
    """

    # Hard clustering of the GO terms:

    # load the embeddings:

    # If soft clustering is or not choosen:

    if soft_clustering == False:
        GO_embeddings = pd.read_csv(
            data_path
            + "_GO_Embeddings_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(dimensions)
            + "_"
            + str(matrix)
            + ".csv",
            index_col=0,
            dtype={0: str},
        )

    elif soft_clustering == True:
        GO_embeddings = np.load(
            data_path
            + "_Soft_Clustering_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "_"
            + str(enrichment)
            + ".npy",
            allow_pickle=True,
        )

    # Load dictionary:

    GO_enrichments = json.load(
        open(
            data_path
            + "Enrichments_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "Per_Cluster.json"
        )
    )

    if intersection == False:
        # Prepare the file:

        GO_terms = GO_embeddings.columns.tolist()
        GO_embeddings_t = GO_embeddings.T
        GO_embeddings_np = GO_embeddings_t.to_numpy()

        # Get the clustering For hard clustering:

        if soft_clustering == False:
            GO_Clusters = Get_Hard_Clusters(GO_embeddings_np, GO_terms)

        # Compare the enrichments per cluster by applying jaccard index:

    elif intersection == True:
        # Get the GO terms that are enriched:

        GO_enrichments_filt = GO_enrichments.values()
        GO_enrichments_filt = set(
            [item for sublist in GO_enrichments_filt for item in sublist]
        )
        GO_emb_filt = set([item for sublist in GO_embeddings for item in sublist])

        # Filters (intersection between two sets):

        filters = GO_enrichments_filt.intersection(GO_emb_filt)

        if soft_clustering == False:
            GO_embeddings = GO_embeddings.T
            GO_embeddings = GO_embeddings[GO_embeddings.index.isin(GO_enrichments_filt)]

            # Prepare the file:

            GO_terms = GO_embeddings.T.columns.tolist()
            GO_embeddings_t = GO_embeddings
            GO_embeddings_np = GO_embeddings_t.to_numpy()

            # Get the clustering:

            GO_Clusters = Get_Hard_Clusters(GO_embeddings_np, GO_terms)

        elif soft_clustering == True:
            # Generate the final result:

            GO_Clusters = [[] for i in range(dimensions)]

            # Fill the clusters only with the intersection of GO terms:

            # For the embeddings:

            for i in range(len(GO_embeddings)):
                GOs = pd.Series(GO_embeddings[i])
                Selection = GOs[GOs.isin(filters)]

                GO_Clusters[i] = list(Selection)

            # For the Enrichments:

            GO_emb_filt = set(GO_emb_filt)
            Cluster_Enrich = [[] for i in range(dimensions)]

            for n in range(len(GO_enrichments)):
                GOs_Enr = pd.Series(GO_enrichments[str(n)])
                Selection_Enr = GOs_Enr[GOs_Enr.isin(filters)]
                Cluster_Enrich[n] = list(Selection_Enr)

            GO_enrichments = Cluster_Enrich

    jaccard_Result = []

    for cluster in range(dimensions):
        # Get the GO in the cluster n:

        GO_in_cluster = GO_Clusters[cluster]  # From GO embeddings.

        if soft_clustering == False:
            GO_enriched_in_cluster = GO_enrichments[str(cluster)]  # From enrichments.
        elif soft_clustering == True:
            GO_enriched_in_cluster = GO_enrichments[cluster]

        # Calculating the jaccard index:

        if len(GO_enriched_in_cluster) > 0 or len(GO_in_cluster) > 0:
            jaccard = jaccard_similarity(GO_in_cluster, GO_enriched_in_cluster)
            jaccard_Result.append(jaccard)

    # Calculating the average of the jaccard index:

    average_jaccard = sum(jaccard_Result) / len(jaccard_Result)

    # Return the information:

    print(average_jaccard)


def Get_Distribution_Nans(P_values):
    """
    Function to test the NaN distributions of the P-values
    Imput:

            - P_values: DataFrame with p-values per eachGO/k
    """
    Annotated = []
    for i in range(len(P_values.index)):
        Annotated.append(len(P_values.columns) - P_values.iloc[i].isnull().sum())

    return Annotated


def Mean_GO_Clusters(
    data_path,
    dimensions,
    network="PPI",
    matrix="PPMI",
    enrichment="BP",
    intersection=False,
):
    GO_enrichments = json.load(
        open(
            data_path
            + "Enrichments_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "Per_Cluster.json"
        )
    )
    GO_enrichments_filt = GO_enrichments.values()
    GO_enrichments_filt = [item for sublist in GO_enrichments_filt for item in sublist]
    print(np.mean(list(Counter(GO_enrichments_filt).values())))


def Get_AUC_ROC_Curve(
    data_path, dimensions, network="PPI", matrix="PPMI", enrichment="BP", p_value=False
):
    # load the embeddings:

    GO_embeddings = pd.read_csv(
        data_path
        + "_GO_Embeddings_"
        + str(enrichment)
        + "_"
        + str(network)
        + "_"
        + str(dimensions)
        + "_"
        + str(matrix)
        + ".csv",
        index_col=0,
    )

    # Load dictionary:

    GO_enrichments = json.load(
        open(
            data_path
            + "Enrichments_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "Per_Cluster.json"
        )
    )

    # Filter the GO embeddings by the enrichments:

    GO_enrichments_filt = GO_enrichments.values()
    GO_enrichments_filt = [item for sublist in GO_enrichments_filt for item in sublist]
    GO_embeddings = GO_embeddings.T
    GO_embeddings = GO_embeddings[GO_embeddings.index.isin(GO_enrichments_filt)]

    # Create the binary version of the GO enrichments:

    Matrix_Gene_GO = pd.concat(
        [pd.Series(v, name=k).astype(str) for k, v in GO_enrichments.items()], axis=1
    )
    Matrix_Gene_GO = pd.get_dummies(Matrix_Gene_GO.stack()).sum(level=1).clip_upper(1)

    # P-value of the ROC Curve (printed):

    if p_value == True:
        Get_P_Value_ROC_Curve(Matrix_Gene_GO, GO_embeddings)

    # Calculate the Area Under the ROC Curve per each dimension:

    roc_scores = [
        roc_auc_score(
            Matrix_Gene_GO.T.loc[:, column], GO_embeddings.loc[:, int(column)]
        )
        for column in Matrix_Gene_GO.T.columns
    ]
    np.mean(roc_scores)

    # Get the mean:

    roc_scores_mean = np.mean(roc_scores)

    # Return the mean:

    print("The auc of ROC is: ", roc_scores_mean)

    return roc_scores_mean


def Get_AUC_PreRec_Curve(
    data_path, dimensions, network="PPI", matrix="PPMI", enrichment="BP", random=True
):
    # load the embeddings:

    GO_embeddings = pd.read_csv(
        data_path
        + "_GO_Embeddings_"
        + str(enrichment)
        + "_"
        + str(network)
        + "_"
        + str(dimensions)
        + "_"
        + str(matrix)
        + ".csv",
        index_col=0,
    )

    # Load dictionary:

    GO_enrichments = json.load(
        open(
            data_path
            + "Enrichments_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "Per_Cluster.json"
        )
    )

    # Filter the GO embeddings by the enrichments:

    GO_enrichments_filt = GO_enrichments.values()
    GO_enrichments_filt = [item for sublist in GO_enrichments_filt for item in sublist]
    GO_embeddings = GO_embeddings.T
    GO_embeddings = GO_embeddings[GO_embeddings.index.isin(GO_enrichments_filt)]

    # Create the binary version of the GO enrichments:

    Matrix_Gene_GO = pd.concat(
        [pd.Series(v, name=k).astype(str) for k, v in GO_enrichments.items()], axis=1
    )
    Matrix_Gene_GO = pd.get_dummies(Matrix_Gene_GO.stack()).sum(level=1).clip_upper(1)

    # Calculate the Area Under the Precission/Recall Curve per each dimension:

    auc_list = []
    for column in Matrix_Gene_GO.T.columns:
        precision, recall, thresholds = metrics.precision_recall_curve(
            Matrix_Gene_GO.T.loc[:, column],
            GO_embeddings.loc[:, int(column)],
            pos_label=1,
        )
        auc_score = metrics.auc(recall, precision)
        auc_list.append(auc_score)

    # Get the mean:

    auc_list_mean = np.mean(auc_list)

    # Random Experiment:

    if random == True:
        auc_list_random_final = []

        for rep in range(100):  # 100 repetitions as a first step
            # Randomly change the rows

            df = Matrix_Gene_GO.T.sample(frac=1).reset_index(drop=True)

            # Calculate the auc:

            auc_list_random = []

            for column in Matrix_Gene_GO.T.columns:
                precision, recall, thresholds = metrics.precision_recall_curve(
                    df.loc[:, column], GO_embeddings.loc[:, int(column)], pos_label=1
                )
                auc_score = metrics.auc(recall, precision)
                auc_list_random.append(auc_score)

            # Get the mean:

            auc_list_random_mean = np.mean(auc_list_random)

            # Save the mean for repetition rep:

            auc_list_random_final.append(auc_list_random_mean)

        # Get the final mean of the random experiment:

        auc_list_random_final_mean = np.mean(auc_list_random_final)

        print(
            f"Auc is:",
            str(auc_list_mean),
            " and the auc of the random (10 rep is) ",
            str(auc_list_random_final_mean),
        )

        return auc_list_mean, auc_list_random_final_mean


def Get_P_Value_ROC_Curve(Binary_Enrichments, GO_embeddings):
    # Prepare the files:

    Binary_Enrichments = Binary_Enrichments.T
    cero_matrix = Binary_Enrichments == 0
    one_matrix = Binary_Enrichments == 1

    # Prepare the lists:

    values_zero = []
    values_one = []

    # Divide the matrices of embeddings:

    for column in Binary_Enrichments.columns:
        GO_embeddings_column = GO_embeddings.loc[:, int(column)]
        cero_matrix_columns = cero_matrix.loc[:, str(column)]
        one_matrix_columns = one_matrix.loc[:, str(column)]

        # For Positives (1):

        GO_embeddings_column_one = GO_embeddings_column[one_matrix_columns]
        values_one.append(GO_embeddings_column_one.tolist())

        # For Negatives (0):

        GO_embeddings_column_cero = GO_embeddings_column[cero_matrix_columns]
        values_zero.append(GO_embeddings_column_cero.tolist())

    # Do the Mann-whitney test

    x, y = mannwhitneyu(values_zero, values_one)

    # Tell the p-value

    print("The p-value of ROC Curve is: ", y)


def Get_Soft_Clustering(
    data_path, dimensions, network="PPI", matrix="PPMI", enrichment="BP"
):
    # load the embeddings:

    GO_embeddings = pd.read_csv(
        data_path
        + "_GO_Embeddings_"
        + str(enrichment)
        + "_"
        + str(network)
        + "_"
        + str(dimensions)
        + "_"
        + str(matrix)
        + ".csv",
        index_col=0,
    )

    GO_embeddings_original = GO_embeddings

    # Load dictionary:

    GO_enrichments = json.load(
        open(
            data_path
            + "Enrichments_"
            + str(enrichment)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
            + "_"
            + str(dimensions)
            + "Per_Cluster.json"
        )
    )

    # Filter the GO embeddings by the enrichments for the threshold:

    GO_enrichments_filt = GO_enrichments.values()
    GO_enrichments_filt = [item for sublist in GO_enrichments_filt for item in sublist]
    GO_embeddings = GO_embeddings.T
    GO_embeddings = GO_embeddings[GO_embeddings.index.isin(GO_enrichments_filt)]

    # Create the binary version of the GO enrichments:

    Matrix_Gene_GO = pd.concat(
        [pd.Series(v, name=k).astype(str) for k, v in GO_enrichments.items()], axis=1
    )
    Matrix_Gene_GO = pd.get_dummies(Matrix_Gene_GO.stack()).sum(level=1).clip(upper=1)

    # Calculate the Area Under the Precission/Recall Curve per each dimension:

    Precision = []
    Recall = []
    Threshold = []

    for column in Matrix_Gene_GO.T.columns:
        precision, recall, thresholds = metrics.precision_recall_curve(
            Matrix_Gene_GO.T.loc[:, column],
            GO_embeddings.loc[:, int(column)],
            pos_label=1,
        )

        # Get thte values :

        Precision.append(list(precision))
        Recall.append(list(recall))
        Threshold.append(list(thresholds))

    # Get the threshold:

    Precision = [item for sublist in Precision for item in sublist]
    Recall = [item for sublist in Recall for item in sublist]
    Threshold = [item for sublist in Threshold for item in sublist]

    f1 = [2 * (p * r) / (p + r + 0.0001) for p, r in zip(Precision, Recall)]
    pMaxF1 = np.argmax(np.array(f1))

    # Get Threshold:

    Final_Threshold = Threshold[pMaxF1]

    # Apply the threshold to the matrix of embeddings (the full matrix, not the filtered one):

    GO = GO_embeddings_original.columns.tolist()
    GO_embeddings_numpy = np.array(GO_embeddings_original.T)
    n, k = np.shape(GO_embeddings_original.T)

    Clusters = [[] for i in range(k)]
    for i in range(n):
        idx = np.argwhere(GO_embeddings_numpy[i] >= Final_Threshold)
        for clust in range(len(idx)):
            Clusters[idx[clust][0]].append(GO[i])

    # Save the clusters:

    np.save(
        data_path
        + "_Soft_Clustering_"
        + str(network)
        + "_"
        + str(matrix)
        + "_"
        + str(dimensions)
        + "_"
        + str(enrichment),
        Clusters,
    )


# Compute PPMI matrix from adjacency matrix A
def Make_PPMI(A, context_window=10, Bin=False):
    print("Computing PPMI matrix")
    degrees = np.sum(A, axis=0)
    volume_of_graph = sum(degrees)

    # fix for disconnected nodes
    degree_corrected = np.array([float(i) for i in degrees])

    for i in range(len(degree_corrected)):
        if degree_corrected[i] == 0.0:
            degree_corrected[i] = 1.0

    diag_degrees_inv = np.diag([1.0 / float(i) for i in degree_corrected])

    tpm = 0.5 * (np.eye(len(A)) + np.matmul(diag_degrees_inv, A))
    power_tpm = np.zeros([len(A), len(A)]) + tpm
    pmi = np.zeros([len(A), len(A)]) + power_tpm

    for r in range(2, context_window + 1):
        # print 'Iteration: ',r
        power_tpm = np.matmul(power_tpm, tpm)
        pmi += power_tpm

    pmi = pmi / float(context_window)
    pmi = volume_of_graph * np.matmul(pmi, diag_degrees_inv)
    pmi_matrix = np.log(pmi, out=np.zeros_like(pmi), where=(pmi != 0))

    for i in pmi_matrix:
        i[i < 0] = 0.0

    if Bin == True:
        for i in pmi_matrix:
            i[i > 0] = 1.0

    return pmi_matrix


# Compute the PPMI matrix of a bipartite adjacency matrix between x and y nodes.
def Make_PPMI_Bipartite(Mat, context_window=10, Bin=False):
    nb_x, nb_y = Mat.shape
    Full_Adj = np.zeros((nb_x + nb_y, nb_x + nb_y))
    for i in range(nb_x):
        # Full_Adj[i][i] = 1.
        for j in range(nb_y):
            Full_Adj[i][nb_x + j] = Mat[i][j]
            Full_Adj[nb_x + j][i] = Mat[i][j]

    # for j in range(nb_y):
    # 	Full_Adj[nb_x+j][nb_x+j] = 1.

    Full_PPMI = Make_PPMI(Full_Adj, context_window, Bin)
    Sub_PPMI = np.zeros((nb_x, nb_y))
    for i in range(nb_x):
        for j in range(nb_y):
            Sub_PPMI[i][j] = Full_PPMI[i][nb_x + j]
    return Sub_PPMI
