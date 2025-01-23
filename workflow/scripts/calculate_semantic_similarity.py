import sys
import json
import numpy as np
import pandas as pd
from multiprocessing import Pool
from pygosemsim import graph, similarity


# 7. Evaluate if the embedding spaces are organized using the FMM

def solve(process):
    Test_Functional_Organization(*process)


def Calculate_Semantic_Similarity(
    dim, save_cosine, common_list=[], filtered=True, Common=False, number_similar=500
):
    """
    This function calculates the semantic similarity between the top500 GO terms pairs that are closer in an embedding space.

         Inputs:
                 - similarity     : pandas DataFrame, with the similarity between the pairs.
                 - dim            : int, dimension of the embedding space.
                 - save_cosine    : string, path to the cosine/cosine matrix of the annotations.
                 - common_list    : list, list to filter the GO terms in the cosine/cosine matrix.
                 - filtered       : booleran, if true the filtering of the matrix is based on its individual list
                 - Common         : booleran, if true the filtering of the matrix is based on the common cancer/control list.
                 - number_similar : int, top pairs to calculate its semantic similarity.
    """

    # Load the cosine matrix:

    cos_matrix_1 = pd.read_csv(save_cosine, index_col=0)

    if filtered == True:
        # Filtered by its individual set:

        cos_matrix_1 = cos_matrix_1.loc[common_list, common_list]

    elif Common == True:
        # Filtered by the comparison set:

        cos_matrix_1 = cos_matrix_1.loc[common_list, common_list]

    else:
        # Without filtering:

        cos_matrix_1 = cos_matrix_1.loc[similarity.index, similarity.index]

    # Get the top similar and dissimilar GO terms based on the corresponding distance measure:

    disimilar_values = (
        cos_matrix_1.mask(np.triu(np.ones(cos_matrix_1.shape)).astype(bool))
        .stack()
        .sort_values(ascending=True)
    )

    top_100_similar_values = disimilar_values[0:number_similar]
    top_100_dissimilar_values = disimilar_values[-number_similar:]

    a = []
    b = []

    G = graph.from_resource("go-basic", snakemake.params.data_path)  # ./Data/go-basic.obo"
    similarity.precalc_lower_bounds(G)

    for i in range(len(top_100_similar_values)):
        GO1 = top_100_similar_values.index.get_level_values(0)[i]
        GO2 = top_100_similar_values.index.get_level_values(1)[i]

        try:
            a.append(similarity.lin(G, GO1, GO2))
        except Exception as PGSSLookupError:
            continue

    for i in range(len(top_100_dissimilar_values)):
        GO1 = top_100_dissimilar_values.index.get_level_values(0)[i]
        GO2 = top_100_dissimilar_values.index.get_level_values(1)[i]

        try:
            b.append(
                similarity.lin(G, GO1, GO2)
            )  # Computes the Lim's semantic similarity, other similarities can be also computed changing this function.
        except Exception as PGSSLookupError:
            continue

    # Relate the distance with the semantic similarity:

    # c = [similarity.values[np.triu_indices(len(similarity), k = 1)]]

    c = []

    # Give back the results:

    if len(b) > 0:
        return (
            np.mean(a),
            np.std(a),
            np.mean(b),
            np.std(b),
            np.mean(c),
            np.std(c),
            top_100_similar_values,
            top_100_dissimilar_values,
        )
    else:
        print("error")


def Compute_Jaccard_Top_Pairs(Cancer, Control):
    Jacc_list = []

    for pair in range(len(Control)):
        # Get the pair:

        pair_control = Control.index[pair]
        reverse_pair_control = (Control.index[pair][1], Control.index[pair][0])

        if pair_control in Cancer.index:
            Jacc_list.append(1)
        elif reverse_pair_control in Cancer.index:
            Jacc_list.append(1)
        else:
            Jacc_list.append(0)

    # Return Jaccard:

    return sum(Jacc_list) / (len(Jacc_list) + len(Jacc_list) - sum(Jacc_list))


def Test_Functional_Organization(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    data_dir,
    matrix="PPMI",
    filtered=True,
    Common=False,
    Jaccard=False,
    number_similar=500,
    annotation="GO",
):
    """
    This fucntion produces two dictionaries that contains the semantic similarity results and the top500 GO terms pairs for Cancer
    and control.

        Inputs:

                 - similarity     : pandas DataFrame, with the similarity between the pairs.
                 - Cancer_type_list      : string list, Name of the Cancer.
                 - Normal_Tissue_list   : string list, Name of the normal tissue.
                 - Cell_type_list       : string list, Name of the normal cell in the tissue.
                 - dimension_list       : int list, list with the dimensions of the embedding spaces.
                 - matrix               : string, network matrix representation.
                 - filtered             : booleran, if true the filtering of the matrix is based on its individual list
                 - Common               : booleran, if true the filtering of the matrix is based on the common cancer/control list.
    """

    # Paths:

    save_cosine = f"{data_dir}/FMM/"
    functional_path = f"{data_dir}/"
    Filter_path = f"{data_dir}/"

    # Calculate the Semantic Similarity:

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        #print("Start Computations")

        # List to save the information:

        Similar_Cancer_mean = []
        Similar_Cancer_std = []

        Disimilar_Cancer_mean = []
        Disimilar_Cancer_std = []

        All_Cancer_mean = []
        All_Cancer_std = []

        Similar_Control_mean = []
        Similar_Control_std = []

        Disimilar_Control_mean = []
        Disimilar_Control_std = []

        All_Control_mean = []
        All_Control_std = []

        Jaccard_Sim = []
        Jaccard_Dis = []

        for dim in dimension_list:
            if annotation == "GO":
                path_Cancer = f"{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}.csv"
                path_Control = (
                    f"{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv"
                )
            else:
                path_Cancer = f"{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_{annotation}.csv"
                path_Control = f"{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_{annotation}.csv"

            if Common == False:
                # Individual list for them:

                if annotation == "GO":
                    Set_Cancer = f"{Filter_path}Filter_{cancer}_Set.csv"
                    Set_Control = f"{Filter_path}Filter_{tissue}_{cell}_Set.csv"
                else:
                    Set_Cancer = f"{Filter_path}Filter_{cancer}_Set_{annotation}.csv"
                    Set_Control = (
                        f"{Filter_path}Filter_{tissue}_{cell}_Set_{annotation}.csv"
                    )

                filt_Cancer = list(pd.read_csv(Set_Cancer)["0"])
                filt_Control = list(pd.read_csv(Set_Control)["0"])

                Cancer_Result = Calculate_Semantic_Similarity(
                    dim, path_Cancer, filt_Cancer, filtered, Common, number_similar
                )
                Control_Result = Calculate_Semantic_Similarity(
                    dim, path_Control, filt_Control, filtered, Common, number_similar
                )
            else:
                # Unique common list for both networks:
                if annotation == "GO":
                    common_list = list(pd.read_csv(f"{Filter_path}Common_Set_{tissue}.csv")["0"])
                else:
                    common_list = list(
                        pd.read_csv(f"{Filter_path}Common_Set_{tissue}_{annotation}.csv")["0"])

                Cancer_Result = Calculate_Semantic_Similarity(
                    dim, path_Cancer, common_list, filtered, Common, number_similar
                )
                Control_Result = Calculate_Semantic_Similarity(
                    dim, path_Control, common_list, filtered, Common, number_similar
                )

            # Keep the results:

            Similar_Cancer_mean.append(Cancer_Result[0])
            Similar_Cancer_std.append(Cancer_Result[1])
            Disimilar_Cancer_mean.append(Cancer_Result[2])
            Disimilar_Cancer_std.append(Cancer_Result[3])
            All_Cancer_mean.append(Cancer_Result[4])
            All_Cancer_std.append(Cancer_Result[5])

            Similar_Control_mean.append(Control_Result[0])
            Similar_Control_std.append(Control_Result[1])
            Disimilar_Control_mean.append(Control_Result[2])
            Disimilar_Control_std.append(Control_Result[3])
            All_Control_mean.append(Control_Result[4])
            All_Control_std.append(Control_Result[5])

            if Jaccard == True:
                Jaccard_Sim.append(
                    Compute_Jaccard_Top_Pairs(Cancer_Result[6], Control_Result[6])
                )
                Jaccard_Dis.append(
                    Compute_Jaccard_Top_Pairs(Cancer_Result[7], Control_Result[7])
                )

        # Save the information:

        save_dic = {
            "Sim_Cancer": Similar_Cancer_mean,
            "Sim_Cancer_std": Similar_Cancer_std,
            "Diss_Cancer": Disimilar_Cancer_mean,
            "Diss_Cancer_std": Disimilar_Cancer_std,
            "All_Cancer": All_Cancer_mean,
            "All_Cancer_std": All_Cancer_std,
            "Sim_Control": Similar_Control_mean,
            "Sim_Control_std": Similar_Control_std,
            "Diss_Control": Disimilar_Control_mean,
            "Diss_Control_std": Disimilar_Control_std,
            "All_Control": All_Control_mean,
            "All_Control_std": All_Control_std,
            "Jaccard_Similar": Jaccard_Sim,
            "Jaccard_Diss": Jaccard_Dis,
        }

        # Save to a file:

        if Common == False:
            note = "Individual_Set"
        else:
            note = "Common_Set"

        sa = json.dumps(save_dic)
        if annotation == "GO":
            f = open(
                f"{functional_path}Similarity_{tissue}_{note}_{number_similar}.json",
                "w",
            )
        else:
            f = open(
                f"{functional_path}Similarity_{tissue}_{note}_{number_similar}_{annotation}.json",
                "w",
            )
        f.write(sa)
        f.close()

        #print("Finish Computations")

# Examples of ways for evaluating the organization:
# Semantic similarity of Top 500 closest/farthest functional annotation embedding vectors:

log_file = open(snakemake.log[0], "w")
sys.stdout = log_file

processes = []
for cancer, tissue, cell in zip(snakemake.params.Cancer_list, snakemake.params.Normal_tissue_list, snakemake.params.Control_list):
    processes.append([
            [cancer], 
            [tissue], 
            [cell], 
            snakemake.params.dimension_list, 
            snakemake.params.data_path,
            "PPMI",
            snakemake.params.filtered,
            snakemake.params.use_common,
            True,
            500,
            snakemake.params.annotation])


with Pool(len(processes)) as pool:
    results = pool.map(solve, processes)
    pool.close()
    pool.join()


log_file.close()
