import sys
import json
import numpy as np
import pandas as pd
import multiprocessing
from collections import Counter

# 8. Obtain the optimal dimensionlity of the space
# 8.1 Calculate the relative error between the spaces
# 8.2 Plot these errors

   
def Relative_Error(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    comparison_list,
    common_list,
    data_dir,
    matrix="PPMI",
    filtered=True,
    Common=False,
    annotation="GO",
):
    """
    This fucntiona calculates the relative error between two dimensional embedding spaces. To do it, it takes the relative position
    of the embedded entities of two different embedding spaces and calculates its relative error.

        Inputs:
                - Cancer_type_list     : string list, Name of the Cancer.
                - Normal_Tissue_list   : string list, Name of the normal tissue.
                - Cell_type_list       : string list, Name of the normal cell in the tissue.
                - dimension_list       : int list, list with the dimensions of the embedding spaces.
                - comparison_list      : string list in the format "d1-d2", pairs of comparisons.
                - common_list          : string list, if Common equal True this list contains the common set of GO terms between Cancer and Control.
                - matrix               : string, network matrix representation.
                - filtered             : bool, if true the relative error is done only for GO terms that has a minimum of 3 genes annotated.
                - Common               : bool, if true the fucntion only takes into account the common set of Cancer and Control annotations.
    """

    # Path:
    

    save_cosine = f"{data_dir}/FMM/"

    control = 0
    bp = None

    # Calculate the Relative error:

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        Cancer_list = []
        Control_list = []

        for dim in comparison_list:
                        
            comparison = dim.split("-")
            d1 = int(comparison[0])
            d2 = int(comparison[1])

            # Load the structure:

            if annotation == "GO":
                Cancer_Structure_1 = pd.read_csv(
                    f"{save_cosine}Cosine_Cancer_{cancer}_{d1}_{matrix}.csv",
                    index_col=0,
                )
                Cancer_Structure_2 = pd.read_csv(
                    f"{save_cosine}Cosine_Cancer_{cancer}_{d2}_{matrix}.csv",
                    index_col=0,
                )

                Control_Structure_1 = pd.read_csv(
                    f"{save_cosine}Cosine_Control_{tissue}_{cell}_{d1}_{matrix}.csv",
                    index_col=0,
                )
                Control_Structure_2 = pd.read_csv(
                    f"{save_cosine}Cosine_Control_{tissue}_{cell}_{d2}_{matrix}.csv",
                    index_col=0,
                )

            else:
                Cancer_Structure_1 = pd.read_csv(
                    f"{save_cosine}Cosine_Cancer_{cancer}_{d1}_{matrix}_{annotation}.csv",
                    index_col=0,
                )
                Cancer_Structure_2 = pd.read_csv(
                    f"{save_cosine}Cosine_Cancer_{cancer}_{d2}_{matrix}_{annotation}.csv",
                    index_col=0,
                )

                Control_Structure_1 = pd.read_csv(
                    f"{save_cosine}Cosine_Control_{tissue}_{cell}_{d1}_{matrix}_{annotation}.csv",
                    index_col=0,
                )
                Control_Structure_2 = pd.read_csv(
                    f"{save_cosine}Cosine_Control_{tissue}_{cell}_{d2}_{matrix}_{annotation}.csv",
                    index_col=0,
                )

            # If we calculate it not only with the common set:

            if filtered == True:
                if annotation == "GO":
                    bp = json.load(open(f"{data_dir}/gene2go_Human_PPIGO_Specific_BP.json"))
                else:
                    bp = json.load(open(f"{data_dir}/gene2go_Human_PPIGO_Specific_BP.json"))

                number_GO = 3
                annotation_list = [name for sublist in bp.values() for name in sublist]
                occurance_of_each_annotation_in_network = Counter(annotation_list)
                terms_filtering = [
                    key
                    for key, values in occurance_of_each_annotation_in_network.items()
                    if values >= number_GO
                ]
                terms_filtering_1 = list(
                    set(terms_filtering).intersection(set(Cancer_Structure_1.index))
                )

                terms_filtering_2 = list(
                    set(terms_filtering).intersection(set(Control_Structure_1.index))
                )

                Cancer_Structure_1 = Cancer_Structure_1.loc[
                    terms_filtering_1, terms_filtering_1
                ]
                Cancer_Structure_2 = Cancer_Structure_2.loc[
                    terms_filtering_1, terms_filtering_1
                ]
                Control_Structure_1 = Control_Structure_1.loc[
                    terms_filtering_2, terms_filtering_2
                ]
                Control_Structure_2 = Control_Structure_2.loc[
                    terms_filtering_2, terms_filtering_2
                ]

            # If we calculate it for the common set of GO terms:

            if Common == True:
                Cancer_Structure_1 = Cancer_Structure_1.loc[common_list, common_list]
                Cancer_Structure_2 = Cancer_Structure_2.loc[common_list, common_list]

                Control_Structure_1 = Control_Structure_1.loc[common_list, common_list]
                Control_Structure_2 = Control_Structure_2.loc[common_list, common_list]

            # Calculate the relative error:

            norm_R1_Cancer = np.linalg.norm(Cancer_Structure_1, ord="fro")
            norm_R2_Cancer = np.linalg.norm(Cancer_Structure_2, ord="fro")

            norm_R1_Control = np.linalg.norm(Control_Structure_1, ord="fro")
            norm_R2_Control = np.linalg.norm(Control_Structure_2, ord="fro")

            Error_Cancer = Cancer_Structure_1 - Cancer_Structure_2
            Error_Control = Control_Structure_1 - Control_Structure_2

            norm_Cancer = np.linalg.norm(Error_Cancer, ord="fro")
            rel_erR1_Cancer = norm_Cancer / max(norm_R1_Cancer, norm_R2_Cancer)

            norm_Control = np.linalg.norm(Error_Control, ord="fro")
            rel_erR1_Control = norm_Control / max(norm_R1_Control, norm_R2_Control)

            Cancer_list.append(rel_erR1_Cancer) 
            Control_list.append(rel_erR1_Control) 

        # Save the corresponding Errors:

        with open(
            f"{save_cosine}Relative_Cancer_{cancer}_{matrix}_{annotation}.txt", "w"
        ) as fp:
            json.dump(Cancer_list, fp)
        with open(
            f"{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}_{annotation}.txt",
            "w",
        ) as fp:
            json.dump(Control_list, fp)

        control = control + 2

        print(f"{control}/{len(Cancer_type_list) * 2}", flush=True)


def Parallel_Error(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    comparison_list,
    data_dir,
    common_list,
    matrix="PPMI",
    filtered=True,
    Common=False,
    annotation="GO",
):
    # Paralelize the tasks:

    Process = []

    for i in range(len(Cancer_type_list)):
        Process.append(
            (
                [Cancer_type_list[i]],
                [Normal_Tissue_list[i]],
                [Cell_type_list[i]],
                dimension_list,
                comparison_list,
                common_list,
                data_dir,
                matrix,
                filtered,
                Common,
                annotation,
            )
        )


    f = open(snakemake.log[0], "w")
    sys.stdout = f

    print(Process, flush=True)
    print(multiprocessing.cpu_count(), flush=True)
    n_cores = snakemake.resources.nodes
    pool1 = multiprocessing.Pool(processes=n_cores)
    pool1.daemon = False
    pool1.starmap(Relative_Error, Process)
    pool1.close()
    pool1.join()
    
    f.close()


comparison_list = []
dim_list=snakemake.params.dimension_list
for i in range(len(dim_list)-1):
    comparison_list.append(f"{dim_list[i]}-{dim_list[i+1]}")
Parallel_Error(
    snakemake.params.Cancer_list,
    snakemake.params.Normal_tissue_list,
    snakemake.params.Control_list,
    snakemake.params.dimension_list,
    comparison_list,
    snakemake.params.data_path,
    common_list=None,
    matrix="PPMI",
    filtered=True,
    Common=False,
    annotation=snakemake.params.annotation,
)


