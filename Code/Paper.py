# This Script Contains all the information to Analyze the fucntional organization of Cancer/Control.

# General Packages:

import os

# Our Packages:

os.chdir("./Code/")

import Script_Paper

# Paths:

save_path = "./Data/"

# 1. Generate the different Cancer/Control Networks:

Control_list       = ["glandular cells", "glandular cells", "glandular cells", "alveolar cells"]
Normal_tissue_list = ["breast", "prostate", "colon", "lung"]
Cancer_list        = ["breast cancer", "prostate cancer", "colorectal cancer", "lung cancer"]

# 1.1 Generate the tissues-specific networks:

for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
    
    Script_Paper.Create_Tissues_specific_Networks(cancer, tissue, control, save_path, 
                                                                     plot = False, Main_Component = True)

# 2. Generate the PPMI matrix representation of the tissues-specific PPI networks:

Script_Paper.Generate_PPMI(Cancer_list, Normal_tissue_list, Control_list, context_window = 10)

# 3. Generate the gene embeddings using NMTF algorithm:

dimension_list = [50, 100, 150, 200, 250, 300]

Script_Paper.Generate_Gene_Embeddings(Cancer_list, Normal_tissue_list, 
                                             Control_list, dimension_list, max_iter=400, verbose = 50)

# 4. Embed the functional annotations in the tissues-specific embedding spaces:

Script_Paper.Annotate_Gene_Space(Cancer_list, Normal_tissue_list, Control_list, 
                                        dimension_list)

# 5. Compute the FMMs:

Script_Paper.Embedding_Structure(Cancer_list, Normal_tissue_list, Control_list,
                                        dimension_list, Gene_Matrix = False, matrix = "PPMI")

# 6. Filter FMMs to contain only annotations shared by cancer and control.

Script_Paper.Cancer_Plot_Scripts.Filter_Set(Cancer_list, Normal_tissue_list, Control_list)

# 7. Evaluate if the gene embedding space is organized using the FMMs:

# 7.1 Semantic similarity of Top 500 closest/farthest functional annotation embedding vectors:

Script_Paper.Parallel_Functional_Organization(Cancer_list, Normal_tissue_list, Control_list, dimension_list,
                    matrix = "PPMI", filtered = True, Common = False, Jaccard = True, number_similar = 500,annotation = "Leaf")

# 7.2 Intra/Inter cluster distance:

Script_Paper.Calculate_Clusters_FMM(Cancer_list, Normal_tissue_list, Control_list,dimension_list)

# 8. Obtain the optimal dimensionlity of the space:

comparison_list = ["50-100", "100-150", "150-200", "200-250","250-300"]

Script_Paper.Parallel_Error(Cancer_list, Normal_tissue_list, Control_list, dimension_list,
                   comparison_list, common_list=None, matrix = "PPMI", filtered = True, Common = False)

# 9. Rank the annotations based on their movement:

optimal_dim = 200 # This is the optimal dimensionality (based on the previous analyses).
Script_Paper.Movement_Ranking(Cancer_list, Normal_tissue_list, 
                                     Control_list, optimal_dim, matrix = "PPMI")

# 10. Link movement with cancer (globally):

Script_Paper.Global_moving_in_the_space()

# 11. Enrichment analyses validation (shifted annotations):

Script_Paper.moving_in_the_space()

# 12. Authomatic literature validation (shifted annotations):

Script_Paper.Query_GO_terms(Normal_tissue_list)

# 13. Define cancer-related genes (based on the literature):

Script_Paper.Query_Common_gene_set(Normal_tissue_list)

# 14. Calculate distances between the genes and the shifted annotations:

Script_Paper.Gen_GO_Common_Space(Cancer_list, Normal_tissue_list, Control_list, optimal_dim)

# 15. Demonstrate that the space is functionally organized (between genes and annotations):

Script_Paper.Demonstrate_Gene_GO_Org(Cancer_list, Normal_tissue_list, Control_list, optimal_dim)

# 16. Do the genes predictions and validate them with an enrichment analyses:

Script_Paper.Do_Fold_Enrichment_Complete(Normal_tissue_list)













































    
    
    
    
    











