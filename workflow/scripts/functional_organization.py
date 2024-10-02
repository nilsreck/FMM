import sys
import json
import pandas as pd
import multiprocessing

# 7. Evaluate if the embedding spaces are organized using the FMM

def Parallel_Functional_Organization(
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
    This function paralelize the computation of Test_Functional_Organization.

        Inputs:

                 - Cancer_type_list      : string list, Name of the Cancer.
                 - Normal_Tissue_list   : string list, Name of the normal tissue.
                 - Cell_type_list       : string list, Name of the normal cell in the tissue.
                 - dimension_list       : int list, list with the dimensions of the embedding spaces.
                 - matrix               : string, network matrix representation.
                 - filtered             : booleran, if true the filtering of the matrix is based on its individual list
                 - Common               : booleran, if true the filtering of the matrix is based on the common cancer/control list.
    """

    # Similarity:

    Process = []

    print(str(number_similar))

    for i in range(len(Cancer_type_list)):
        Process.append(
            (
                [Cancer_type_list[i]],
                [Normal_Tissue_list[i]],
                [Cell_type_list[i]],
                dimension_list,
                data_dir,
                matrix,
                filtered,
                Common,
                Jaccard,
                number_similar,
                annotation,
            )
        )

    n_cores = multiprocessing.cpu_count()
    pool1 = Pool(processes=n_cores)
    pool1.starmap(Test_Functional_Organization, Process)
    pool1.close()
    pool1.join()


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
    functional_path = f"{data_dir}/Data/"
    Filter_path = f"{data_dir}/Data/"

    # Calculate the Semantic Similarity:

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        print("Start Computations")

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
                    common_list = list(
                        pd.read_csv(f"{Filter_path}Common_Set_{tissue}.csv")["0"]
                    )
                else:
                    common_list = list(
                        pd.read_csv(
                            f"{Filter_path}Common_Set_{tissue}_{annotation}.csv"
                        )["0"]
                    )

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

        print("Finish Computations")


# Examples of ways for evaluating the organization:
# Semantic similarity of Top 500 closest/farthest functional annotation embedding vectors:
Parallel_Functional_Organization(
    snakemake.params.Cancer_list,
    snakemake.params.Normal_tissue_list,
    snakemake.params.Control_list,
    snakemake.params.dimension_list, 
    snakemake.params.data_path,
    matrix="PPMI",
    filtered=True,
    Common=True,
    Jaccard=True,
    number_similar=500,
    annotation=snakemake.params.annotation,
)