# Main

"""
The code showed below represents an expample of how to use our new FMM-based methodology to
analyze the functional organization of four different cancer-types and their corresponding control
tissues.

Note that the code can be easily generalized to be applied with other networks.
"""

import Cancer_Scripts_Human_Case_Control
import Cancer_Plot_Scripts
import os

# Set the base directory dynamically
base_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Use os.path.join to create paths in a platform-independent way
data_directory = os.path.join(base_directory, "Cleaned_version", "Data", "")

# 1. Generate the different Cancer/Control Networks:

Control_list = [
    "glandular cells",
    "glandular cells",
    "glandular cells",
    "alveolar cells",
]
Normal_tissue_list = ["breast", "prostate", "colon", "lung"]
Cancer_list = ["breast cancer", "prostate cancer", "colorectal cancer", "lung cancer"]

# 1.1 Generate the networks:

#for control, tissue, cancer in zip(Control_list, Normal_tissue_list, Cancer_list):
#    Cancer_Scripts_Human_Case_Control.Create_Tissue_Specific_Network(
#        cancer, tissue, control, data_directory, plot=False
#    )

# 1.2 Plot the Venn Diagrams:

#Cancer_Plot_Scripts.Venn_Diagram_Cancer_Networks(
#    Cancer_list, Normal_tissue_list, Control_list
#)

# 1.3 Get the statistics for the networks:

#Cancer_Plot_Scripts.Network_Statistics(Cancer_list, Normal_tissue_list, Control_list)

# 2. Generate the PPMI Matrices:

#Cancer_Plot_Scripts.Generate_PPMI(
#    Cancer_list, Normal_tissue_list, Control_list, context_window=10
#)

# 3. Generate the Gene embeddings:

dimension_list = [50, 100, 150, 200, 250, 300]

#Cancer_Plot_Scripts.Generate_Gene_Embeddings(
#    Cancer_list,
#    Normal_tissue_list,
#    Control_list,
#    dimension_list,
#    max_iter=1000,
#    verbose=50,
#)

# 4 Map Functional annotations (GO terms) to the embedding spaces:

#Cancer_Plot_Scripts.Annotate_Gene_Space(
#    Cancer_list, Normal_tissue_list, Control_list, dimension_list, Annotation="Leaf"
#)

# 5 Generate the FMMs:

#Cancer_Plot_Scripts.Embedding_Structure(
#    Cancer_list,
#    Normal_tissue_list,
#    Control_list,
#    dimension_list,
#    Gene_Matrix=False,
#    matrix="PPMI",
#    annotation="Leaf",
#)

# 6. Filter FMMs to contain only annotations shared by cancer and control.

#Cancer_Plot_Scripts.Common_GO_Terms(
#    Cancer_list, Normal_tissue_list, Control_list, annotation="Leaf"
#)

# 7. Evaluate if the embedding spaces are organized using the FMM:

# Examples of ways for evaluating the organization:

# Semantic similarity of Top 500 closest/farthest functional annotation embedding vectors:

#Cancer_Plot_Scripts.Parallel_Functional_Organization(
#    Cancer_list,
#    Normal_tissue_list,
#    Control_list,
#    dimension_list,
#    matrix="PPMI",
#    filtered=True,
#    Common=True,
#    Jaccard=True,
#    number_similar=500,
#    annotation="Leaf",
#)

# 8. Obtain the optimal dimensionlity of the space:

# 8.1 Calculate the relative error between the spaces:

comparison_list = ["50-100", "100-150", "150-200", "200-250", "250-300"]

#Cancer_Plot_Scripts.Parallel_Error(
#    Cancer_list,
#    Normal_tissue_list,
#    Control_list,
#    dimension_list,
#    comparison_list,
#    common_list=None,
#    matrix="PPMI",
#    filtered=True,
#    Common=False,
#    annotation="Leaf",
#)

# 8.2 Plot these errors:

#Cancer_Plot_Scripts.Plot_Relative_Error(
#    Cancer_list,
#    Normal_tissue_list,
#    Control_list,
#    dimension_list,
#    comparison_list,
#    matrix="PPMI",
#    annotation="GO",
#)

# Based on these results we consider the following dimensionaity as the optimal:

optimal_dim = 200

# 9. Application to cancer:

# 9,1 Movement of the GO terms:

Cancer_Plot_Scripts.Movement_Ranking(
    Cancer_list, Normal_tissue_list, Control_list, optimal_dim, matrix="PPMI"
)

# 9.2 Link movement with cancer (globally):

Cancer_Plot_Scripts.Global_moving_in_the_space()

# 9.3 Enrichment analyses validation (shifted annotations):

Cancer_Plot_Scripts.moving_in_the_space()

# 9.4 Literature validation:

Cancer_Plot_Scripts.Query_GO_terms(
    Normal_tissue_list
)  # Read the description of the funciton for details about the needed files.

# 9.5 Define cancer-related genes (based on the literature):

Cancer_Plot_Scripts.Query_Common_gene_set(
    Normal_tissue_list
)  # Read the description of the funciton for details about the needed files.

# 9.6 Assess the functional organization of genes and functions in the embedding space:

# Distance between genes and functions:

Cancer_Plot_Scripts.Gen_GO_Common_Space(
    Cancer_list, Normal_tissue_list, Control_list, optimal_dim
)

# Demonstrate that the space is functionally organized (between genes and annotations):

Cancer_Plot_Scripts.Demonstrate_Gene_GO_Org(
    Cancer_list, Normal_tissue_list, Control_list, optimal_dim
)

# 16. Do the genes predictions and validate them with an enrichment analyses:

Cancer_Plot_Scripts.Do_Fold_Enrichment_Complete(Normal_tissue_list)
