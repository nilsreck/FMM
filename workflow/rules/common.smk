
import os
import subprocess

# configurations = config.get("configurations", dict())

# Control_list = configurations["Control_list"]
# Normal_tissue_list = configurations["Normal_tissue_list"]
# Cancer_list = configurations["Cancer_list"]
# plot_network = configurations["plot_network"]
# dimension_list = configurations["dimension_list"]
# max_iter = configurations["max_iter"]
# annotation = configurations["Annotation"]
# optimal_dim = configurations["optimal_dim"]

# data_path = os.path.abspath("../Data")
# resource_path = os.path.abspath("./resources")


def get_genes(Control_list, Normal_tissue_list, Cancer_list, data_path):
    inputs = []
    for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
        inputs.append(data_path+"/Cancer_"+cancer+"_Genes.csv")
        inputs.append(data_path+"/Control_"+tissue+"_"+control+"_Genes.csv")
    return inputs

def get_Normal_tissue_string(Normal_tissue_list):
    return ' '.join(Normal_tissue_list)

def get_enrichment_script():
    return str(os.path.abspath("scripts/cancer_enrichment.py"))

def get_enrichment_log():
    return os.path.abspath("logs/cancer_enrichment.log")

def get_enrichment_cmd(network_wrapper, data_path, annotation, NCBI_account, Normal_tissue_list):
    script = get_enrichment_script()
    log = get_enrichment_log()
    tissues = get_Normal_tissue_string(Normal_tissue_list)

    script = f"python {script} {data_path} {annotation} {NCBI_account} {log} {tissues}"
    work_dir = str(os.path.abspath("./"))

    cmd = network_wrapper
    cmd = cmd.replace("{script}", script)
    cmd = cmd.replace("{log}", log)
    cmd = cmd.replace("{workdir}", work_dir)
    return cmd

def get_networks(Control_list, Normal_tissue_list, Cancer_list, data_path):
    networks = []
    for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
        networks.append(data_path+"/Cancer_"+cancer+"_PPI.npy")
        networks.append(data_path+"/Control_"+tissue+"_"+control+"_PPI.npy")
    return networks

    
def get_PPMI_matrices(Control_list, Normal_tissue_list, Cancer_list, data_path):
    matrices = []
    for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
        matrices.append(data_path+"/Cancer_PPMI_"+cancer+".npy")
        matrices.append(data_path+"/Control_PPMI_"+tissue+"_"+control+".npy")
    return matrices


def get_dimensions(dimension_list, resource_path):
    dimensions = []
    for dim in dimension_list:
        dimensions.append(resource_path+"/dimension_"+str(dim))
    return dimensions


def get_comparison_list(dimension_list):
    comparison_list = []
    for i in range(len(dimension_list)-1):
        comparison_list.append(f"{dimension_list[i]}-{dimension_list[i+1]}")
    return comparison_list


def get_comparisons(dimension_list, resource_path):
    comparisons = []
    for i in range(len(dimension_list)-1):
        comparisons.append(f"{resource_path}/comparison_{dimension_list[i]}-{dimension_list[i+1]}")
    return comparisons


def get_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, data_path):
    embeddings = []
    for dim in dimension_list:
        for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
            embeddings.append(data_path+"/Embeddings/_P_Matrix_"+str(dim)+"_PPI_"+cancer+"_Adj_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_P_Matrix_"+str(dim)+"_PPI_"+cancer+"_PPMI_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_P_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_Adj_Control.npy")
            embeddings.append(data_path+"/Embeddings/_P_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_PPMI_Control.npy")
            
            embeddings.append(data_path+"/Embeddings/_U_Matrix_"+str(dim)+"_PPI_"+cancer+"_Adj_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_U_Matrix_"+str(dim)+"_PPI_"+cancer+"_PPMI_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_U_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_Adj_Control.npy")
            embeddings.append(data_path+"/Embeddings/_U_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_PPMI_Control.npy")
            
            embeddings.append(data_path+"/Embeddings/_G_Matrix_"+str(dim)+"_PPI_"+cancer+"_Adj_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_G_Matrix_"+str(dim)+"_PPI_"+cancer+"_PPMI_Cancer.npy")
            embeddings.append(data_path+"/Embeddings/_G_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_Adj_Control.npy")
            embeddings.append(data_path+"/Embeddings/_G_Matrix_"+str(dim)+"_PPI_"+tissue+"_"+control+"_PPMI_Control.npy")
    return embeddings


def get_annotations(annotation, data_path):
    if annotation == "GO":
        return data_path+"/_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv"
    elif annotation == "Leaf":
        return data_path+"/_Matrix_Genes_GO_BP_PPI.csv"


def get_mapped_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path):
    embeddings = []
    for dim in dimension_list:
        for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
            embeddings.append(data_path+"/Embeddings/_GO_Embeddings_"+annotation+"_PPI_"+cancer+"_"+str(dim)+"_Adj_Cancer.csv")
            embeddings.append(data_path+"/Embeddings/_GO_Embeddings_"+annotation+"_PPI_"+cancer+"_"+str(dim)+"_PPMI_Cancer.csv")
            embeddings.append(data_path+"/Embeddings/_GO_Embeddings_"+annotation+"_PPI_"+tissue+"_"+control+"_"+str(dim)+"_Adj_Control.csv")
            embeddings.append(data_path+"/Embeddings/_GO_Embeddings_"+annotation+"_PPI_"+tissue+"_"+control+"_"+str(dim)+"_PPMI_Control.csv")
    return embeddings


def get_FMMs(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path):
    fmms = []
    for dim in dimension_list:
        for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
            fmms.append(data_path+"/FMM/Cosine_Cancer_"+cancer+"_"+str(dim)+"_PPMI_"+annotation+".csv")
            fmms.append(data_path+"/FMM/Cosine_Control_"+tissue+"_"+control+"_"+str(dim)+"_PPMI_"+annotation+".csv")
    return fmms


def get_common(Normal_tissue_list, annotation, data_path):
    common = []
    for tissue in Normal_tissue_list:
        if annotation == "GO":
            common.append(data_path+"/Common_Set_"+tissue+".csv")
        else:
            common.append(data_path+"/Common_Set_"+tissue+"_"+annotation+".csv")
    return common


def get_similarity(Normal_tissue_list, annotation, data_path): 
    similarities = []
    for tissue in Normal_tissue_list:         
        if annotation == "GO":
            similarities.append(data_path+"/Similarity_"+tissue+"_Common_Set_500.json")
        else:
            similarities.append(data_path+"/Similarity_"+tissue+"_Common_Set_500_"+annotation+".json")
    return similarities


def get_error(Control_list, Normal_tissue_list, Cancer_list, annotation, data_path):
    errors = []
    for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
        errors.append(data_path+"/Relative_Cancer_"+cancer+"_PPMI_"+annotation+".txt")
        errors.append(data_path+"/Relative_Control_"+tissue+"_"+control+"_PPMI_"+annotation+".txt")
    return errors


def get_movements(Normal_tissue_list, annotation, data_path):
    movements = []
    for tissue in Normal_tissue_list: 
        movements.append(data_path+"/Rank_movement_"+tissue+"_PPMI_"+annotation+".csv")
    return movements


def get_top_movements(Normal_tissue_list, data_path):
    movements = []
    for tissue in Normal_tissue_list: 
        movements.append(data_path+"/top_100_moving_"+tissue+".csv")
    return movements
    

def get_transformed(Normal_tissue_list, data_path):
    transformed = []
    for tissue in Normal_tissue_list: 
        transformed.append(data_path+"/Transformed_Common_"+tissue+".csv")
    return transformed


def get_movement_table(Normal_tissue_list, data_path):
    tables = []
    for tissue in Normal_tissue_list: 
        tables.append(data_path+"/Top_moving_"+tissue+"_Table.csv")
    return tables


def get_counts(Normal_tissue_list, data_path):
    counts = []
    for tissue in Normal_tissue_list: 
        counts.append(data_path+"/Cancer_Count_"+tissue+".csv")
    return counts


def get_distances(Normal_tissue_list, data_path):
    distances = []
    for tissue in Normal_tissue_list: 
        distances.append(data_path+"/Cancer_Gene_GO_dist_"+tissue+".csv")
        distances.append(data_path+"/Control_Gene_GO_dist_"+tissue+".csv")
    return distances


def get_error_plot(annotation, data_path):
    return data_path+"/Relative_Error_"+annotation+".png"


def get_predictions(Normal_tissue_list, data_path):
    predictions = []
    for tissue in Normal_tissue_list: 
        predictions.append(data_path+"/Predictions_Rank_"+tissue+".csv")
        predictions.append(data_path+"/Predictions_Rank_"+tissue+".csv")
    return predictions