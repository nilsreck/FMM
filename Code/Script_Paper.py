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
from scipy.linalg import svd
import itertools
import json
import os
import os
import re
import time
import json
import math
import random
import requests
import Ontology
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Entrez
import multiprocessing
import scipy.stats as ss   
from ast import literal_eval
from matplotlib import lines
from selenium import webdriver
from collections import Counter
import matplotlib.pyplot as plt
from multiprocessing import Pool
from pySankey.sankey import sankey
from itertools import combinations
from scipy.spatial import distance
from string import ascii_lowercase
import matplotlib.patches as mpatches
from pygosemsim import graph, similarity
from sklearn_extra.cluster import KMedoids
from matplotlib.ticker import NullFormatter
from matplotlib_venn import venn2_unweighted
from sklearn.metrics import pairwise_distances
from statsmodels.stats.multitest import multipletests
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.stats import mannwhitneyu, pearsonr, gaussian_kde
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
import random
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu,hypergeom


def Create_Tissues_specific_Networks(save_path, path_pathology, Main_Component = True):
    
    '''
    This function generate the tissues-specific PPI network:
        
        - save_path      : string, path to save the Adj matrix.
        - Main_Component : booleran, if True the Adj matrix produced contains the main connected part of the network.
        - path_pathology : string, path to the folder in which the file pathology.tsv is. 
    '''

    save_path = "./Data/" 
    
    cancer_pb = pd.read_csv(f'{path_pathology}', sep='\t')
    
    # Cancer list (the ones we do not have, 16 more):
    
    cancers = list(set(cancer_pb.Cancer))    
    
    # Generate the PPIs:
    
    for tiss in cancers:
        
        cancer = cancer_pb[cancer_pb.Cancer == tiss] # Cancer type
    
        # Cancer Expressed genes:
        
        cancer["Total appearances"] = cancer["High"] + cancer["Medium"] + cancer["Low"]
        cancer["Percentage of Appear"] = cancer["Total appearances"]/(cancer["Total appearances"]+ cancer["Not detected"])
        gene_cancer_expresed = cancer[cancer["Total appearances"].isna() == False]
        gene_cancer_expresed = gene_cancer_expresed[gene_cancer_expresed["Percentage of Appear"] >= 0.5].Gene    

        # Load the PPI network:
        
        file_name_array    = "Human_Biogrid_Adj"
        file_name_genelist = "Human_Biogrid_Genes"
        data_path          = "./Data/" 
        
        PPI_adj   = np.load(data_path + file_name_array + "_PPI_.npy", allow_pickle= True ) 
        PPI_genes = pd.read_table(data_path + file_name_genelist + "_PPI_", header = None, dtype={0: str})
        
        PPI_df = pd.DataFrame(PPI_adj, PPI_genes[0], PPI_genes[0])   
        
        # Change the name of the genes:
        
        ensmbl_code = pd.read_csv("./Data/" + "mart_export.txt", sep='\t', dtype={1: str})
        ensmbl_code = ensmbl_code.dropna()     
        ensmbl_code_cancer  = ensmbl_code[ensmbl_code["Gene stable ID"].isin(gene_cancer_expresed)]["NCBI gene (formerly Entrezgene) ID"]      
        common_ensmbl_PPI         = set(ensmbl_code_cancer).intersection(set(PPI_genes[0]))
    
        PPI_Cancer  = PPI_df.loc[common_ensmbl_PPI,common_ensmbl_PPI]
         
        if Main_Component == True:
        
        # Get the biggest connected components:
        
            PPI_Cancer_Network  = nx.from_pandas_adjacency(PPI_Cancer)           
            connected_components_Cancer  = sorted((PPI_Cancer_Network.subgraph(c) for c in nx.connected_components(PPI_Cancer_Network)), key=len, reverse=True)
            main_component_PPI_Cancer  = connected_components_Cancer[0]           
            node_list_Cancer  =  list(main_component_PPI_Cancer.nodes())
            Cancer_np  = nx.to_numpy_array(main_component_PPI_Cancer)
            
            print("The len of the main component in Cancer is: ", str(len(node_list_Cancer)), " / ", str(len(PPI_Cancer.index)), "\n")

         # Save the information:
        
        np.save(save_path + "Pan_Cancer_" + str(tiss) + "_PPI", Cancer_np)        
        pd.DataFrame(node_list_Cancer).to_csv(save_path + "Pan_Cancer_" + str(tiss) + "_Genes.csv", index=False)
 

def deep_walk_ppmi(adj_matrix, context_window=10):
    '''
    Input  : Adj Matrix as numpy array
    Output :  Numpy array with the PPMI matrix of the Adj.
    '''
    degrees = np.sum(adj_matrix,axis = 0)
    volume_of_graph = sum(degrees)
    diag_degrees_inv = np.diag([1/float(i) for i in degrees])

    pmi = np.zeros([len(adj_matrix),len(adj_matrix)])

    transition_matrix = np.matmul(diag_degrees_inv,adj_matrix)
    pmi_temp = transition_matrix
    pmi += transition_matrix

    for r in range(1,context_window):
        print("Iteration:")
        pmi_temp = np.matmul(pmi_temp,transition_matrix)
        pmi += pmi_temp

    pmi = pmi/float(context_window)
    pmi = volume_of_graph*np.matmul(pmi,diag_degrees_inv)

    pmi_matrix = np.log(pmi, out=np.zeros_like(pmi), where=(pmi!=0))

    for i in pmi_matrix:
        i[i<0]=0

    return pmi_matrix


def Pann_Generate_PPMI(context_window = 10): 
    
    '''
    This function generates the PPMI matrix from the Adj Matrices and save it as a np:
        
        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - context_window  : int, contex window for the PPMI matrix.
    '''    
    
    # Networks path:
    
    network_path = "./Data/" 
    
    cancer_pb = pd.read_csv("./Data/" + "pathology.tsv", sep='\t')
    
    # Cancer list (the ones we do not have, 16 more):
    
    cancers = list(set(cancer_pb.Cancer)) 
    
    for cancer  in cancers:
        
        # Load the Adj: 
        
        Cancer_adj  = np.load(f'{network_path}Pan_Cancer_{cancer}_PPI.npy', allow_pickle= True)
        
        # Generate the PPMI:
        
        PPMI_Cancer  = deep_walk_ppmi(Cancer_adj,  context_window)
        
        # Save the PPMI:
        
        np.save(f'{network_path}Pan_Cancer_PPMI_{cancer}', PPMI_Cancer)


class SNMTF(object):
    """Compute Partial Orthonormal Non-negative Matrix Tri-Factorization (NMTF) Embeddings"""
    def __init__(self, max_iter=1000, verbose = 10):

        super(SNMTF,self).__init__()
        self.max_iter = max_iter
        self.verbose = verbose


    def Score(self, R1, P, U, G, norm_R1):
        GT = np.transpose(G)

        ErR1 = R1 - np.matmul(P, np.matmul(U, GT))
        norm_erR1 = np.linalg.norm(ErR1, ord='fro')

        rel_erR1 = norm_erR1 / norm_R1

        return norm_erR1, rel_erR1


    def Get_Clusters(self, M, nodes):
        n,k = np.shape(M)
        Clusters = [[] for i in range(k)]
        for i in range(n):
            idx = np.argmax(M[i])
            Clusters[idx].append(nodes[i])
        return Clusters
    
    
    def SVD_Matrix(self, R1, matrix,file_name_array_path):
        
        print("Generating the SVD Matrix")
        
        J,K,L = svd(R1)
        
        file_name_array_path_1 = file_name_array_path + "_J_Matrix_" + str(matrix) + "_Human"
        file_name_array_path_2 = file_name_array_path + "_K_Matrix_" + str(matrix) + "_Human"
        file_name_array_path_3 = file_name_array_path + "_L_Matrix_" + str(matrix) + "_Human"
        
        np.save(file_name_array_path_1, J)
        np.save(file_name_array_path_2, K)
        np.save(file_name_array_path_3, L)
        

    def Solve_MUR_Hum(self, R1, k1, file_name_array_path, network = "PPI", matrix = "PPMI", init="rand"):
        
        print("Starting NMTF")
        n,n=np.shape(R1)
        norm_R1 = np.linalg.norm(R1, ord='fro')

        if init == "rand":
            print(" * Using random initialization")
            P = np.random.rand(n,k1) + 1e-5
            U = np.random.rand(k1,k1) + 1e-5
            G = np.random.rand(n,k1) + 1e-5
        elif init == 'SVD':
            # Initialize matrix factor P, U, and G using SVD decomposition

            print(" * -- Eig decomposition on R1 to get PxUxG^T")
            
            # Load the SVD to be faster:
            
            print(" * -- Loading the files for the SVD:")
            
            J = np.load(file_name_array_path + "_J_Matrix_" + str(matrix) + "_Human.npy")
            K = np.load(file_name_array_path + "_K_Matrix_" + str(matrix) + "_Human.npy")
            L = np.load(file_name_array_path + "_L_Matrix_" + str(matrix) + "_Human.npy")
            
            # The rest depends on the number of dimensions:
            
            P = np.zeros((n,k1)) + 1e-5
            G = np.zeros((n,k1)) + 1e-5
            U = np.zeros((k1,k1)) + 1e-5

            for i in range(n):
                for j in range(k1):
                    if J[i][j] > 0.:
                        P[i][j] = J[i][j]

            for i in range(k1):
                U[i][i] = abs(K[i])

            for i in range(n):
                for j in range(k1):
                    if L[i][j] > 0.:
                        G[i][j] = L[i][j]

        else :
            print("Unknown initializer: %s"%(init))
            exit(0)
        OBJ, REL2 = self.Score(R1,P, U, G, norm_R1)
        print(" - Init:\t OBJ:%.4f\t REL2:%.4f"%(OBJ, REL2))


        #Begining M.U.R.
        for it in range(1,self.max_iter+1):

            R1T = np.transpose(R1)
            UT = np.transpose(U)

            PT = np.transpose(P)
            PT_P = np.matmul(PT,P)

            GT = np.transpose(G)
            GT_G = np.matmul(GT,G)

            #update rule for G

            #R1 matrix

            #update rule for P

            R1_G_UT = np.matmul(np.matmul(R1, G), UT)
            U_GT_G_UT = np.matmul(np.matmul(U,GT_G),UT)

            P_U_GT_G_UT = np.matmul(P,U_GT_G_UT)

            P_mult = np.sqrt( np.divide(R1_G_UT,P_U_GT_G_UT + 1e-14) )

            #update rule for U

            PT_R1_G = np.matmul(PT, np.matmul(R1, G))
            PT_P_U_GT_G = np.matmul(PT_P, np.matmul(U, GT_G))

            U_mult = np.sqrt( np.divide(PT_R1_G,PT_P_U_GT_G + 1e-14))

            #update rule for G

            R1T_P_U = np.matmul(np.matmul(R1T, P), U)
            G_GT_R1T_P_U = np.matmul(G, np.matmul(GT, R1T_P_U))

            G_mult = np.sqrt( np.divide(R1T_P_U, G_GT_R1T_P_U + 1e-14) )

            # Applying M.U.R.
            P = np.multiply(P, P_mult) + 1e-14
            U = np.multiply(U, U_mult) + 1e-14
            G = np.multiply(G, G_mult) + 1e-14


            if (it%self.verbose == 0) or (it==1):
                OBJ, REL2 = self.Score(R1, P, U, G, norm_R1)
                print(" - It %i:\t OBJ:%.4f\t REL2:%.4f"%(it, OBJ, REL2))
         
        # Save the matrices:
        
        file_name_array_path_1 = file_name_array_path + "_P_Matrix_" + str(k1) + "_" + str(network) + "_" + str(matrix)
        file_name_array_path_2 = file_name_array_path + "_U_Matrix_" + str(k1) + "_" + str(network) + "_" + str(matrix)
        file_name_array_path_3 = file_name_array_path + "_G_Matrix_" + str(k1) + "_" + str(network) + "_" + str(matrix)
        
        np.save(file_name_array_path_1, P)
        np.save(file_name_array_path_2, U)
        np.save(file_name_array_path_3, G)
        
        return P, U, G        


def Generate_Gene_Embeddings(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, max_iter=400, verbose = 50):
    
    '''
    This function generates the Gene Embedding from the PPMI and Adj:
    
        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - context_window  : int, contex window for the PPMI matrix.
                - max_iter        : int, iterations for the NMTF.
                - verbose         : int.
    '''    
    # Networks path:
    
    network_path = "./Data/" 
    
    counter = 0
    
    # Generate the embeddings:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        Solver = SNMTF(max_iter, verbose) 
        
        # Load the information:
        
        
        Cancer_PPMI  = np.load(f'{network_path}Cancer_PPMI_{cancer}.npy', allow_pickle= True)
        Control_PPMI = np.load(f'{network_path}Control_PPMI_{tissue}_{cell}.npy', allow_pickle= True) 
        
        # Do the embeddings:
        
        for dim in dimension_list:
        
            Solver.Solve_MUR(Cancer_PPMI, dim, network_path, network = f'PPI_{cancer}', matrix = "PPMI_Cancer", init="SVD")
            Solver.Solve_MUR(Control_PPMI, dim, network_path, network = f'PPI_{tissue}_{cell}', matrix = "PPMI_Control", init="SVD")

            counter = counter + 4
            
            print(f'{counter}/{len(Cancer_type_list)*4*(len(dimension_list))}')        

        
def GO_Embeddings_Human(Matrix_Gene_GO, gene_embeddings, save_path, GO = "BP", network = "PPI", matrix = "PPMI",
                  ortogonal = False, init = "SVD", repetitions = 200, save = True):
    
    '''
    This function embeddes the functional annotations into a given gene embedding space.
       
    Imputs:
            
        Matrix_Gene_GO : DataFrame with genes and GO terms (0 : Not annotated, 1 : Annotated)
        gene_embeddings: Array with the embeddings of the genes.
                  
    '''
    
    # Matrix with gene embeddings (if the previous function was used to create the gene/GO the order of the genes is the same):
    
    gene_embeddings_db = pd.DataFrame(gene_embeddings)
    gene_embeddings_db.index = Matrix_Gene_GO.index
    
    # Calculate the GO embeddings:
    
    # If is ortogonal V transpose is equal to V^-1 then we can use it we want to get the GO embeddings.
    # This is equals to (V^-1) * R = W
    
    if ortogonal == True:
        
        gene_embeddings_db_inverse = pd.DataFrame(np.linalg.pinv(gene_embeddings_db.values), gene_embeddings_db.columns, gene_embeddings_db.index)
        GO_embeddings = gene_embeddings_db_inverse.dot(Matrix_Gene_GO)
        
        if save == True:
            save_path_1 = save_path + "_GO_Embeddings_" + str(GO) + "_" + str(network) + "_"+ str(len(gene_embeddings_db.columns))+ "_" + str(matrix) + ".csv"
            GO_embeddings.to_csv(save_path_1)
        else:
            return(GO_embeddings)
        
def Annotate_Gene_Space(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list):
    
    '''
    This function generates the Annotation of the gene embedding space:
    
        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - context_window  : int, contex window for the PPMI matrix.
                - Annotation      : string, annotation for the space.
    ''' 
    
    # Networks path:
    
    network_path = "./Data/" 
    save_path    = "./Data/" 
    
    counter = 0
    
    # Get the annotations:
          
    GO_Matrix = pd.read_csv(f'"./Data/_Matrix_Genes_GO_BP_PPI.csv', 
                                index_col = 0, dtype={0: str}) 
    GO_Matrix.index = GO_Matrix.index.astype(str)  
               
    # Start the embeddings:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
               
        # Read the gene files: 
            
            Cancer_genes      = pd.read_table(f'{network_path}Cancer_{cancer}_Genes.csv', header = 0, dtype={0: str})
            Cancer_genes_list = list(Cancer_genes["0"]) 
        
            Control_genes      = pd.read_table(f'{network_path}Control_{tissue}_{cell}_Genes.csv', header = 0, dtype={0: str})
            Control_genes_list = list(Control_genes["0"])
            
            # Subset the annotation to keep ontly the genes:
        
            GO_Gene_Matrix_Cancer  = GO_Matrix[GO_Matrix.index.isin(Cancer_genes_list)]
            GO_Gene_Matrix_Control = GO_Matrix[GO_Matrix.index.isin(Control_genes_list)]
            
            # Delete annotations without genes (> 0):
            
            GO_Cancer          = GO_Gene_Matrix_Cancer.sum(axis=0)
            filter_Cancer      = set(GO_Cancer[GO_Cancer > 0].index)
        
            GO_Control          = GO_Gene_Matrix_Control.sum(axis=0)
            filter_Control      = set(GO_Control[GO_Control > 0].index)
            
            GO_Gene_Matrix_Cancer  = GO_Gene_Matrix_Cancer[filter_Cancer]
            GO_Gene_Matrix_Control = GO_Gene_Matrix_Control[filter_Control]
            
            # Reorder the Matrices:
    
            Annotation_Cancer  = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
            Annotation_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]
            
            for dim in dimension_list:
            
                # Fix the gene vectorial representation:
                
                G_Cancer_PPMI  = np.load(f'{network_path}_G_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)
                G_Control_PPMI = np.load(f'{network_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True)

                # Do the embeddings:
                
                GO_Embeddings_Human(Annotation_Cancer,  G_Cancer_PPMI, 
                                              save_path,  GO = "Leaf", network = f'PPI_{cancer}', 
                                              matrix = "PPMI_Cancer",  ortogonal = True)
                
                GO_Embeddings_Human(Annotation_Control,  G_Control_PPMI, 
                                              save_path,  GO = "Leaf", network = f'PPI_{tissue}_{cell}', 
                                              matrix = "PPMI_Control",  ortogonal = True)
            
                counter = counter + 4
            
                print(f'{counter}/{len(Cancer_type_list)*4*(len(dimension_list))}')        
        
    
def Embedding_Structure(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, Gene_Matrix = False, matrix = "PPMI"
                    ):
    
    '''
    This function calculates the relative position of any embedded entity (GO terms by defaut). It producess a csv with the corresponging
    results.
    
        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - Gene_Matrix     : bool, if the entities are genes.
                - matrix          : string, network matrix representation.
    ''' 
    
    # Path:
    
    network_path = "./Data/" 
    save_cosine  = "./Data/" 
    
    control = 0
    
    # Get the structure of the spaces:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        for dim in dimension_list:
            
            path_Cancer  = f'{network_path}_GO_Embeddings_Leaf_PPI_{cancer}_{dim}_{matrix}_Cancer.csv'
            path_Control = f'{network_path}_GO_Embeddings_Leaf_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv'                    
                
            embeddings_Cancer  = pd.read_csv(path_Cancer, index_col = 0).T
            embeddings_Control = pd.read_csv(path_Control, index_col = 0).T
                
            embeddings_Cancer[embeddings_Cancer < 0] = 0
            embeddings_Control[embeddings_Control < 0] = 0
                
            embeddings_Cancer_index  = embeddings_Cancer.index
            embeddings_Control_index = embeddings_Control.index
                
            embeddings_Cancer  = np.array(embeddings_Cancer)
            embeddings_Control = np.array(embeddings_Control)
                
            cosine_Cancer  = pairwise_distances(embeddings_Cancer, metric="cosine")
            cosine_Control = pairwise_distances(embeddings_Control, metric="cosine")
        
                
            cosine_Cancer  = pd.DataFrame(cosine_Cancer,embeddings_Cancer_index,embeddings_Cancer_index)
            cosine_Control = pd.DataFrame(cosine_Control,embeddings_Control_index,embeddings_Control_index)
                                   
            cosine_Cancer.to_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Leaf.csv', header = True, index=True) 
            cosine_Control.to_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Leaf.csv', header = True, index=True)
                    
            control = control + 2
                
            print(f'{control}/{len(Cancer_type_list) * 2 * len(dimension_list)}')
                
    
def Filter_Set(Cancer_type_list, Normal_Tissue_list, Cell_type_list):

    '''
    For each Control/Cancer filter GO terms with less of 3 annotations. It produces a panda Series with the list of GO terms
    per each of the samples.
    
        Inputs:
                - Cancer_type_list     : string list, Name of the Cancer.
                - Normal_Tissue_list   : string list, Name of the normal tissue.
                - Cell_type_list       : string list, Name of the normal cell in the tissue.
    ''' 
    
    path = "./Data/" 
    save = "./Data/" 
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load GO Matrix:

            
        GO_Matrix = pd.read_csv(f'./Data/_Matrix_Genes_GO_BP_PPI.csv', 
                                    index_col = 0, dtype={0: str}) 
        GO_Matrix.index = GO_Matrix.index.astype(str) 
        
        # Load the GO/gene matrix:
        
        Cancer_genes      = pd.read_table(f'{path}Cancer_{cancer}_Genes.csv', header = 0, dtype={0: str})
        Cancer_genes_list = list(Cancer_genes["0"]) 
        
        Control_genes      = pd.read_table(f'{path}Control_{tissue}_{cell}_Genes.csv', header = 0, dtype={0: str})
        Control_genes_list = list(Control_genes["0"])        
        
        # Subset the GO embeddings to keep the same genes:
    
        GO_Gene_Matrix_Cancer  = GO_Matrix[GO_Matrix.index.isin(Cancer_genes_list)]
        GO_Gene_Matrix_Control = GO_Matrix[GO_Matrix.index.isin(Control_genes_list)]
        
        # Reorder the Matrices:

        GO_Gene_Matrix_Cancer  = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
        GO_Gene_Matrix_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]
    
        # Cancer:
    
        GO_Cancer                = GO_Gene_Matrix_Cancer.sum(axis=0)
        GO_terms_filtered_Cancer = set(GO_Cancer[GO_Cancer >= 3].index)
        
        # Control:
        
        GO_Control                = GO_Gene_Matrix_Control.sum(axis=0)
        GO_terms_filtered_Control = set(GO_Control[GO_Control >= 3].index)
        
        # Save them:

    pd.DataFrame(GO_terms_filtered_Cancer).to_csv(f'{save}Filter_{cancer}_Set_Leaf.csv') 
    pd.DataFrame(GO_terms_filtered_Control).to_csv(f'{save}Filter_{tissue}_{cell}_Set_Leaf.csv')     
    
        

def Parallel_Functional_Organization(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                    matrix = "PPMI", filtered = True, Common = False, Jaccard = False, number_similar = 500, annotation = "Leaf"):
    
    '''
   This function paralelize the computation of Test_Functional_Organization.
   
       Inputs:
    
                - Cancer_type_list      : string list, Name of the Cancer.
                - Normal_Tissue_list   : string list, Name of the normal tissue.
                - Cell_type_list       : string list, Name of the normal cell in the tissue.
                - dimension_list       : int list, list with the dimensions of the embedding spaces.
                - matrix               : string, network matrix representation.
                - filtered             : booleran, if true the filtering of the matrix is based on its individual list
                - Common               : booleran, if true the filtering of the matrix is based on the common cancer/control list.
    ''' 
    
    # Similarity:
    
    Process = []
    
    print(str(number_similar))
    
    for i in range(len(Cancer_type_list)):
        
        Process.append(([Cancer_type_list[i]],  [Normal_Tissue_list[i]], 
                        [Cell_type_list[i]], dimension_list, matrix, filtered, Common, Jaccard,number_similar,annotation))
                        
    n_cores = multiprocessing.cpu_count()
    pool1   = Pool(processes  = n_cores)
    pool1.starmap(Test_Functional_Organization, Process)
    pool1.close()
    pool1.join()        

def Compute_Jaccard_Top_Pairs(Cancer, Control):
    
    '''
    This function produce the JI between a set of pairs of annotations:
        - Cancer  : list, contain the top 500 closest and top 500 farthes GO terms in the cancer space.
        - Control : list, contain the top 500 closest and top 500 farthes GO terms in the control space.
    '''       
    Jacc_list = []
       
    for pair in range(len(Control)):
        
        # Get the pair:
        
        pair_control         = Control.index[pair]
        reverse_pair_control = (Control.index[pair][1], Control.index[pair][0])
        
        if pair_control in Cancer.index:
            Jacc_list.append(1)
        elif reverse_pair_control in Cancer.index:
            Jacc_list.append(1)
        else:
            Jacc_list.append(0)
           
    # Return Jaccard:
        
    return(sum(Jacc_list)/(len(Jacc_list) + len(Jacc_list) - sum(Jacc_list)))
          
def Test_Functional_Organization(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                   matrix = "PPMI", filtered = True, Common = False, Jaccard = False, number_similar = 500, annotation = "GO"):
    
    '''
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
    ''' 
    
    # Paths:
    
    save_cosine     = "./Data/" 
    functional_path = "./Data/" 
    Filter_path     = "./Data/" 
    
    # Calculate the Semantic Similarity:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        print("Start Computations")
        
        # List to save the information:
        
        Similar_Cancer_mean = []
        Similar_Cancer_std  = []
        
        Disimilar_Cancer_mean = []
        Disimilar_Cancer_std  = []
        
        All_Cancer_mean = []
        All_Cancer_std  = []        
        
        Similar_Control_mean = []
        Similar_Control_std  = []
        
        Disimilar_Control_mean = []
        Disimilar_Control_std  = []
    
        All_Control_mean = []
        All_Control_std  = []  
        
        Jaccard_Sim = []
        Jaccard_Dis = []

        for dim in dimension_list:
            
            if annotation == "GO":
                path_Cancer  = f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}.csv'
                path_Control = f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv'
            else:
                path_Cancer  = f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_{annotation}.csv'
                path_Control = f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_{annotation}.csv'
            
            if Common == False: 
                
                # Individual list for them:
                
                if annotation == "GO":
                    Set_Cancer  = f'{Filter_path}Filter_{cancer}_Set.csv' 
                    Set_Control = f'{Filter_path}Filter_{tissue}_{cell}_Set.csv'
                else:
                    Set_Cancer  = f'{Filter_path}Filter_{cancer}_Set_{annotation}.csv' 
                    Set_Control = f'{Filter_path}Filter_{tissue}_{cell}_Set_{annotation}.csv'                    
                
                filt_Cancer  = list(pd.read_csv(Set_Cancer)["0"])
                filt_Control = list(pd.read_csv(Set_Control)["0"])
                
                Cancer_Result  = Calculate_Semantic_Similarity(dim, path_Cancer, filt_Cancer, filtered, Common, number_similar)
                Control_Result = Calculate_Semantic_Similarity(dim, path_Control, filt_Control, filtered, Common, number_similar) 
            else:
                
                # Unique common list for both networks:
                if annotation == "GO":
                    common_list = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}.csv')["0"])
                else: 
                    common_list = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_{annotation}.csv')["0"])
                    
                Cancer_Result  = Calculate_Semantic_Similarity(dim, path_Cancer, common_list, filtered, Common, number_similar)
                Control_Result = Calculate_Semantic_Similarity(dim, path_Control, common_list, filtered, Common, number_similar)     
                
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
                
               Jaccard_Sim.append(Compute_Jaccard_Top_Pairs(Cancer_Result[6], Control_Result[6]))
               Jaccard_Dis.append(Compute_Jaccard_Top_Pairs(Cancer_Result[7], Control_Result[7]))
                
        # Save the information:
        
        save_dic = {"Sim_Cancer"       : Similar_Cancer_mean,
                    "Sim_Cancer_std"   : Similar_Cancer_std,
                    "Diss_Cancer"      : Disimilar_Cancer_mean,
                    "Diss_Cancer_std"  : Disimilar_Cancer_std,
                    "All_Cancer"       : All_Cancer_mean,
                    "All_Cancer_std"   : All_Cancer_std,
                    "Sim_Control"      : Similar_Control_mean,
                    "Sim_Control_std"  : Similar_Control_std,
                    "Diss_Control"     : Disimilar_Control_mean,
                    "Diss_Control_std" : Disimilar_Control_std,
                    "All_Control"      : All_Control_mean,
                    "All_Control_std"  : All_Control_std,
                    "Jaccard_Similar"  : Jaccard_Sim,
                    "Jaccard_Diss"     : Jaccard_Dis}
        
        # Save to a file:
        
        if Common == False:
            note = "Individual_Set"
        else:
            note = "Common_Set"
        
        sa = json.dumps(save_dic)
        if annotation == "GO":
            f = open(f'{functional_path}Similarity_{tissue}_{note}_{number_similar}.json' ,"w")
        else:
            f = open(f'{functional_path}Similarity_{tissue}_{note}_{number_similar}_{annotation}.json' ,"w")
        f.write(sa)
        f.close()
        
        print("Finish Computations")        

def Calculate_Semantic_Similarity(dim, save_cosine, common_list = [], filtered = True, Common = False, number_similar= 500):
    
    '''
   This function calculates the semantic similarity between the top500 GO terms pairs that are closer in an embedding space.
    
        Inputs:
                - similarity     : pandas DataFrame, with the similarity between the pairs.
                - dim            : int, dimension of the embedding space.
                - save_cosine    : string, path to the cosine/cosine matrix of the annotations.
                - common_list    : list, list to filter the GO terms in the cosine/cosine matrix.
                - filtered       : booleran, if true the filtering of the matrix is based on its individual list
                - Common         : booleran, if true the filtering of the matrix is based on the common cancer/control list.
                - number_similar : int, top pairs to calculate its semantic similarity.
    ''' 
    
    # Load the cosine matrix:

    cos_matrix_1 = pd.read_csv(save_cosine, index_col = 0)
    
    if filtered == True:
        
        # Filtered by its individual set:
        
        cos_matrix_1 =  cos_matrix_1.loc[common_list, common_list]
        
    elif Common == True:
        
        # Filtered by the comparison set:
        
        cos_matrix_1 = cos_matrix_1.loc[common_list, common_list]
        
    else:
        
        # Without filtering:
        
        cos_matrix_1 = cos_matrix_1.loc[similarity.index, similarity.index ]
        
    # Get the top similar and dissimilar GO terms based on the corresponding distance measure:
    
    disimilar_values = cos_matrix_1.mask(np.triu(np.ones(cos_matrix_1.shape)).astype(bool)).stack().sort_values(ascending=True)

    top_100_similar_values    = disimilar_values[0:number_similar]
    top_100_dissimilar_values = disimilar_values[-number_similar:]
    
    a = []
    b = []
    
    G = graph.from_resource("go-basic")
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
            b.append(similarity.lin(G, GO1, GO2))
        except Exception as PGSSLookupError:
            continue   

    # Give back the results:
    c = []
    if len(b) > 0:
        return(np.mean(a), np.std(a), np.mean(b), np.std(b), np.mean(c), np.std(c), top_100_similar_values,top_100_dissimilar_values)
    else:
        print("error") 
        

def Calculater_Intra_Inter_SS(Semantic, Cluster_db_1):
    
    intra_distance = []
    inter_distance = []
                    
    # Intra cluster semantic similarity:
                    
    for cluster in list(set(Cluster_db_1.Cluster)):
                        
        # Get the intra cluster distance:
                        
        GO_cluster_1 = Cluster_db_1[Cluster_db_1.Cluster == cluster].GO                    
        Semantic_1   = Semantic.loc[GO_cluster_1, GO_cluster_1]
        Semantic_1   = np.array(Semantic_1)[np.triu_indices(np.array(Semantic_1).shape[0], k = 1)]
        intra_distance.append(round(np.mean(Semantic_1),3))
                    
        # Inter cluster semantic similarity:
                    
        lista = list(set(Cluster_db_1.Cluster))
                    
        for i in range(len(lista)):
                        
            cluster_1    = lista[i]
            GO_cluster_1 = Cluster_db_1[Cluster_db_1.Cluster == cluster_1].GO 
                        
            for j in range(i+1, len(lista)):
                            
                cluster_2    = lista[j]
                GO_cluster_2 = Cluster_db_1[Cluster_db_1.Cluster == cluster_2].GO 
                            
                # Take the semantic similarity:
                            
                Semantic_1 = Semantic.loc[GO_cluster_1, GO_cluster_2]
                inter_distance.append(round(Semantic_1.mean().mean(),3))
        
    return(intra_distance,inter_distance)
        

def Calculate_Clusters_FMM(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, Common = True): 
    
    # Paths:

    save_cosine    = "./Data/" 
    Path_Semantic  = "./Data/" 
    Filter_path    = "./Data/" 
    path_organit   = "./Data/" 
       
    # Load Semantic Similarity between the GO terms in human:
    
    Semantic = pd.read_csv(f'{Path_Semantic}Semantic_Human.csv', index_col = 0)
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Set the dictionary:
            
        keys = []
            
        for dm in dimension_list:
            keys.append(f'{dm}_Inter')
            keys.append(f'{dm}_Intra')
            
        Result_Cancer  = dict(([(key, []) for key in keys]))
        Result_Control = dict(([(key, []) for key in keys]))
        
        for dim in dimension_list: 
            
            print(dim)
        
            # Load the corresponding cosine:
            
            path_Cancer  = f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_PPMI_Leaf.csv'
            path_Control = f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_PPMI_Leaf.csv'
        
            cos_matrix_1 = pd.read_csv(path_Cancer, index_col  = 0)
            cos_matrix_2 = pd.read_csv(path_Control, index_col = 0)
            
            # Filter the cosines:
            
            common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Leaf.csv')["0"])
            
            cos_matrix_1  = cos_matrix_1.loc[common_set, common_set]
            cos_matrix_2  = cos_matrix_2.loc[common_set, common_set]
            
            # K medoids:
            
            n_clusters_1 = round(math.sqrt(len(cos_matrix_1.index)/2))
            n_clusters_2 = round(math.sqrt(len(cos_matrix_2.index)/2))
            
            grouping_1 = KMedoids(n_clusters=n_clusters_1, metric='precomputed',).fit(cos_matrix_1)
            grouping_2 = KMedoids(n_clusters=n_clusters_2, metric='precomputed',).fit(cos_matrix_2)

            Cluster_db_1 = pd.DataFrame(cos_matrix_1.index)  
            Cluster_db_2 = pd.DataFrame(cos_matrix_2.index) 

            Cluster_db_1.columns = ["GO"]    
            Cluster_db_2.columns = ["GO"]     
            
            Cluster_db_1["Cluster"] = grouping_1.labels_
            Cluster_db_2["Cluster"] = grouping_2.labels_ 
            
            # For Cancer:
            
            intra_distance_1, inter_distance_1 = Calculater_Intra_Inter_SS(Semantic, Cluster_db_1)
            intra_distance_2, inter_distance_2 = Calculater_Intra_Inter_SS(Semantic, Cluster_db_2)
                   
            # Add the info to the list:
                    
            Result_Cancer[f'{dim}_Inter'] = inter_distance_1
            Result_Cancer[f'{dim}_Intra'] = intra_distance_1
            
            Result_Control[f'{dim}_Inter'] = inter_distance_2
            Result_Control[f'{dim}_Intra'] = intra_distance_2
            
            x = mannwhitneyu(intra_distance_1,inter_distance_1, alternative = "greater")
            y = mannwhitneyu(intra_distance_2,inter_distance_2, alternative = "greater")
            
            print(x)
            print(y)
            
            # Save the information:
                    
        sa = json.dumps(Result_Cancer)
        f = open(f'{path_organit}Cancer_{tissue}_Cluster_SS.json' ,"w")
        f.write(sa)
        f.close() 

        sa = json.dumps(Result_Control)
        f = open(f'{path_organit}Control_{tissue}_Cluster_SS.json' ,"w")
        f.write(sa)
        f.close()      
                
def Parallel_Error(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                   comparison_list, common_list, matrix = "PPMI", filtered = True, Common = False, annotation = "GO"):
    
    # Paralelize the tasks:
    
    Process = []
    
    for i in range(len(Cancer_type_list)):
        Process.append(([Cancer_type_list[i]],  [Normal_Tissue_list[i]], 
                        [Cell_type_list[i]], dimension_list, comparison_list, 
                        common_list, matrix, filtered, Common,annotation))
    
    n_cores = multiprocessing.cpu_count()
    pool1   = Pool(processes  = n_cores)
    pool1.starmap(Relative_Error, Process)
    pool1.close()
    pool1.join()        
        



def Relative_Error(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                   comparison_list, common_list, matrix = "PPMI", filtered = True, Common = False, annotation = "GO"):
    '''
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
    ''' 
    
    # Path:
    
    save_cosine  = "./Data/" 
    
    control = 0
    
    # Calculate the Relative error:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        Cancer_list  = []
        Control_list = []
        
        for dim in comparison_list:
            
            # Break the comparions:
            
            comparison = dim.split("-")
            d1 = int(comparison[0])
            d2 = int(comparison[1])
            
            # Load the structure:

                
            Cancer_Structure_1  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d1}_{matrix}_{annotation}.csv', index_col = 0)
            Cancer_Structure_2  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d2}_{matrix}_{annotation}.csv', index_col = 0)
                
            Control_Structure_1 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d1}_{matrix}_{annotation}.csv', index_col = 0)
            Control_Structure_2 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d2}_{matrix}_{annotation}.csv', index_col = 0)
            
            # If we calculate it not only with the common set:
            
            if filtered == True:

                BP = json.load(open("./Data/gene2go_Human_PPIGO_Specific_BP.json"))
                
                number_GO = 3
                annotation_list = [name for sublist in BP.values() for name in sublist]
                occurance_of_each_annotation_in_network = Counter(annotation_list)
                terms_filtering = [key for key,values in occurance_of_each_annotation_in_network.items() if values >= number_GO]
                
                terms_filtering_1 = list(set(terms_filtering).intersection(set(Cancer_Structure_1.index)))
                
                terms_filtering_2 = list(set(terms_filtering).intersection(set(Control_Structure_1.index)))
                
                Cancer_Structure_1 = Cancer_Structure_1.loc[terms_filtering_1,terms_filtering_1]
                Cancer_Structure_2 = Cancer_Structure_2.loc[terms_filtering_1,terms_filtering_1]
                Control_Structure_1 = Control_Structure_1.loc[terms_filtering_2,terms_filtering_2]
                Control_Structure_2 = Control_Structure_2.loc[terms_filtering_2,terms_filtering_2] 
            
            # If we calculate it for the common set of GO terms:
            
            if Common == True:
                
                Cancer_Structure_1 = Cancer_Structure_1.loc[common_list, common_list]
                Cancer_Structure_2 = Cancer_Structure_2.loc[common_list, common_list]
                
                Control_Structure_1 = Control_Structure_1.loc[common_list, common_list]
                Control_Structure_2 = Control_Structure_2.loc[common_list, common_list]
            
            # Calculate the relative error:
            
            norm_R1_Cancer = np.linalg.norm(Cancer_Structure_1, ord='fro')
            norm_R2_Cancer = np.linalg.norm(Cancer_Structure_2, ord='fro')
            
            norm_R1_Control = np.linalg.norm(Control_Structure_1, ord='fro')
            norm_R2_Control = np.linalg.norm(Control_Structure_2, ord='fro')
            
            Error_Cancer  = Cancer_Structure_1  - Cancer_Structure_2
            Error_Control = Control_Structure_1 - Control_Structure_2

            norm_Cancer     = np.linalg.norm(Error_Cancer, ord='fro')
            rel_erR1_Cancer = norm_Cancer/max(norm_R1_Cancer,norm_R2_Cancer)
            
            norm_Control     = np.linalg.norm(Error_Control, ord='fro')
            rel_erR1_Control = norm_Control/max(norm_R1_Control, norm_R2_Control) 
            
            Cancer_list.append(rel_erR1_Cancer)
            Control_list.append(rel_erR1_Control)
        
        # Save the corresponding Errors:
        
        with open(f'{save_cosine}Relative_Cancer_{cancer}_{matrix}.txt', "w") as fp:
            json.dump(Cancer_list, fp)
        with open(f'{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}.txt', "w") as fp:
            json.dump(Control_list, fp)
        
        control = control + 2
        
        print(f'{control}/{len(Cancer_type_list) * 2}')        
        
def Movement_Ranking(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI", annotation = "Leaf"):
    
    # Paths:
    
    cosine_path   = "./Data/" 
    Filter_path   = "./Data/"   
    movement_path = "./Data/"  
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_Leaf.csv',         index_col = 0)
        Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_Leaf.csv', index_col = 0)
            
        # Filter the matrices with the common set:
            
        common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Leaf.csv')["0"])
            
        Cancer_structure  = Cancer_structure.loc[common_set, common_set]
        Control_structure = Control_structure.loc[common_set, common_set]    
        
        # Substract the cancer to the control:
        
        Substract = Control_structure - Cancer_structure
        
        # Calculate the total movement of the GO terms in the space by the 1-norm of the vectors:
        
        Movement = Substract.apply(np.linalg.norm, axis=1) 
        
        # Rank the movement:
        
        Movement = Movement.sort_values(ascending=False)
        
        # Save the rank:
        
        Movement.to_csv(f'{movement_path}Rank_movement_{tissue}_{matrix}_{annotation}.csv') 
                        
def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='100000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results   

def Query_GO_terms(Normal_tissue_list):
    
    '''
    This functions query the shifted annotations in pubmed. The output is a pandas dataframe with
    the name of the annotation (its definition) and the number of publications in pubmed.
    
    Input:
            - Normal_tissue_list : list of strings, name of the tissue (breast, lung, prostate, and colon)       
    '''
    
    movement_path = "./Data/" 
    path_rank     = "./Data/" 
    
    for tissue in Normal_tissue_list:
        
        rank      = pd.read_csv(f'{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv')
        moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"]
        
        file     = rank.head(len(moving_GO))       
        cancer_q = f'{tissue} cancer'        
        counts   = []
        
        for row in range(len(file)): 
    
            query         = f'({file.iloc[row]["GO name"]}[Text Word]) AND ({cancer_q}[Text Word])'
            query_file    = search(query)
            citation_list = query_file["IdList"]
            
            counts.append(len(citation_list))
        
        Result                   = pd.DataFrame()
        Result["Annotation"]     = list(file["GO name"])
        Result["norm"]           = list(file["norm"])
        Result["Cancer Related"] = list(file["Cancer Related"])
        Result["Bibliography"]   = counts
        
        Result.to_csv(f'{movement_path}Top_moving_{tissue}_Table.csv')        


def Global_moving_in_the_space():
    
    '''
    This function test if the cancer-related annotations are moving statistically more than the rest.
    '''
    
    # Load our set of cancer-related annotations:
    
    with open(f'./Data/enriched_in_cancer_gobp_terms_cosmic.txt', 'r') as fp:
        cancer_related_go_terms =json.load(fp)
    
    # Do the test
        
    cancers = ["breast","prostate","lung","colon"]
    
    for i in cancers:

        Ranking = pd.read_csv(f'.\Data\Rank_movement_{i}_PPMI_Leaf.csv',index_col=0,names=['norm'],header=0)

        # Globally validate the movement is connected with cancer:

        a = Ranking.loc[~Ranking.index.isin(cancer_related_go_terms)]['norm']
        b = Ranking.loc[Ranking.index.isin(cancer_related_go_terms)]['norm']
        print(i,mannwhitneyu(a,b,alternative='less'))        
       

def moving_in_the_space():
    
    '''
    This function performs the enrichment analyses of cancer-related annotations.
    '''
    
    # Load our set of cancer-related annotations:
    
    with open(f'./Data/enriched_in_cancer_gobp_terms_cosmic.txt', 'r') as fp:
        cancer_related_go_terms =json.load(fp)
        
    # Set the number of shifted or stable annotations based on 2std:
    
    percentile_moving = {'lung' : 53 , 'breast' : 58, 'colon' : 68, 'prostate' : 49}
    percentile_stable = {'lung' : 15 , 'breast' : 29, 'colon' : 22, 'prostate' : 26}  
    
    # To save the values of the enrichment analyses

    percentage_most_moving = []
    percentage_less_moving = []
    percentage_by_random = []

    success_in_the_population = len(cancer_related_go_terms)
    
    # Do the enrichment per each cancer type:

    cancers = ["breast","prostate","lung","colon"]
    for i in cancers:
        number_of_draws = percentile_moving[i]
        number_of_draws_stable = percentile_stable[i]

        print('\n')
        Ranking = pd.read_csv(f'.\Data\Rank_movement_{i}_PPMI_Leaf.csv',index_col=0,names=['norm'],header=0)

        # Observed by random
        population_size = len(Ranking)
        probability_of_success = success_in_the_population/population_size
        expected_random = number_of_draws * probability_of_success
        print(expected_random)
        
        top_moving_norm_df = Ranking.sort_values(by='norm',ascending=False)[0:number_of_draws]
        top_moving_norm = set(top_moving_norm_df.index)
        top_success = len(top_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_head = round(hypergeom.sf(top_success, population_size, success_in_the_population, number_of_draws),5)

        less_moving_norm_df = Ranking.sort_values(by='norm',ascending=False)[-number_of_draws_stable:]
        less_moving_norm = set(less_moving_norm_df.index)
        less_success = len(less_moving_norm.intersection(cancer_related_go_terms))
        pvalue_of_tail = round(hypergeom.sf(less_success, population_size, success_in_the_population, number_of_draws_stable),5)

        print(f"% cancer related genes in the top {number_of_draws}:",top_success,"and the corresponding p-value",pvalue_of_head)
        print(f"% cancer related genes in the less {number_of_draws_stable}:",less_success,"and the corresponding p-value",pvalue_of_tail)

        percentage_most_moving.append(top_success*100/number_of_draws)
        percentage_less_moving.append(less_success*100/number_of_draws_stable)
        percentage_by_random.append(expected_random*100/number_of_draws)
    
    # Do the plot:
    
    plt.rcParams.update({'font.size': 15})

    plt.figure(figsize=(8,5))
    plt.plot(percentage_most_moving,'--o',label=f' shifted annotations',linewidth=2)
    plt.plot(percentage_less_moving,'--o',label=f' stable annotations',linewidth=2)
    plt.plot(percentage_by_random,'--o',label=f' expected by random',linewidth=2)
    plt.ylabel('% cancer-related annotations',fontsize=16,fontweight='bold')
    plt.xticks(range(4),["Breast","Prostate","Lung","Colon"],fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(bbox_to_anchor=(1,1), loc="upper left",fontsize=16)
    plt.tight_layout()
    plt.savefig('cancer_enrichments_2std.svg',dpi=500)
    plt.show()       
        
def Query_Common_gene_set(Normal_tissue_list):
    
    '''
    This function takes the common set of genes between cancer and control and performs a systematic
    literature search. We use this functions to define our set of cancer-related genes.
    '''   
    # Paths:
    
    gene_GO_path = "./Data/"
    
    # For each tissue:
    
    for tissue in Normal_tissue_list: 
        
        common_genes = pd.read_csv(f'{gene_GO_path}Transformed_Common_{tissue}.csv')
        common_genes = common_genes.drop_duplicates(subset=['initial_alias'])
        common_genes = common_genes[["initial_alias", "name"]]
        common_genes = common_genes.dropna()
        common_genes = common_genes.reset_index(drop = True)

        cancer_q = f'{tissue} cancer'
        counts = []
        
        for row in range(len(common_genes)):
            
            query         = f'({common_genes.iloc[row]["name"]}[Text Word]) AND ({cancer_q}[Text Word])'
            query_file    = search(query)
            citation_list = query_file["IdList"] 
            
            counts.append(len(citation_list))

        # Prepare the data frame:
    
        common_genes["Counts"] = counts 
        
        common_genes.to_csv(f'{gene_GO_path}Cancer_Count_{tissue}.csv',   
                            header = True, index=True)        


def Distance_GO_Gene(vector_GO, matrix_genes):
    
    vector_GO = vector_GO.loc[:,~vector_GO.columns.duplicated()]
    
    # Cosine
    
    result_pd_cos = pd.DataFrame(0., index=np.arange(len(matrix_genes.index)), 
                             columns=list(vector_GO.columns))
    
    result_pd_cos.column = vector_GO.columns
    result_pd_cos.index = matrix_genes.index
    
    for gene in list(matrix_genes.index):
        
        gene_vector = matrix_genes.loc[gene]
        
        for GO in list(vector_GO.columns):
            
            GO_vector = vector_GO[GO]
            
            result_pd_cos.loc[gene, GO] = distance.cosine(GO_vector, gene_vector)
            
    return(result_pd_cos)

def Generate_Final_DataFrame(GO_Gene_Matrix, distances):
    
    # Lists to fill:
    
    group  = [] # Cluster based on the distance to the GO term
    shape  = [] # Shape for the plot (if the gene clusteres is annotated with the corresponding GO)
    size   = [] # Size for the plot (if the gene clusteres is annotated with the corresponding GO)
    dist   = [] # Minimum distance (the minimim distance to a GO term)
    GO_min = [] # The closest GO term in the space
    
    # Iteration per gene:

    for gene in list(distances.index):
        
        it_gene       = distances.loc[gene]
        it_annotation = GO_Gene_Matrix.loc[gene]
    
        if min(it_gene) <= 0.5:
            group.append(it_gene.idxmin())
            
            if it_annotation.loc[it_gene.idxmin()] == 1:
                shape.append("Annotated")
                size.append(75)
                dist.append(min(it_gene))
                GO_min.append(it_gene.idxmin())
            else:
                shape.append("NO Annotated")
                size.append(55) 
                dist.append(min(it_gene))
                GO_min.append(it_gene.idxmin())
        else:
            group.append("NO")
            shape.append("NO Annotated")
            size.append(15)
            dist.append(min(it_gene))
            GO_min.append(it_gene.idxmin())
    
    # Fill the final dataframe:
    
    Cosine_it = distances.copy()
    
    Cosine_it["Group"]     = group
    Cosine_it["Shape"]     = shape
    Cosine_it["size"]      = size
    Cosine_it["Min_Dist"]  = dist
    Cosine_it["Min_GO"]    = GO_min
    
    # Return the dataFrame:
    
    return(Cosine_it)

         
def Gen_GO_Common_Space(Cancer_list, Normal_tissue_list, Cell_type_list, dim):
    
    '''
    This functions calculates the cosine distance between the gene/GO common space. This functions
    needs the function "Distance_GO_Gene" to calculate the cosine distances and the function "Generate_Final_DataFrame"
    to generate the final dataframe.
    '''
    
    # Paths:
    
    network_path = "./Data/"
    gene_GO_path = "./Data/"
    moving_path  = "./Data/"
    
    # For each combination:
    
    for cancer, tissue, cell in zip(Cancer_list,Normal_tissue_list, Cell_type_list):
        
        # Prepare the paths:
        
        # Genes:
        
        genes_Control =  f'{network_path}Control_{tissue}_{cell}_Genes.csv'
        genes_Cancer  =  f'{network_path}Cancer_{cancer}_Genes.csv'
        
        # Gene embedded coordinates:
        
        Gene_embedd_Control = f'{network_path}_P_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        Gene_embedd_Cancer  =  f'{network_path}_P_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy'
        
        # Middle factor:
        
        midle_Control = f'{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        midle_Cancer  = f'{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        
        Gene_midle_Control_db  = pd.DataFrame(np.load(midle_Control, allow_pickle=True))
        Gene_midle_Cancer_db   = pd.DataFrame(np.load(midle_Cancer,  allow_pickle=True))
        
        # GO embedded coordinates:
        
        GO_embedd_Control = f'{network_path}_GO_Embeddings_GO_PPI_{tissue}_{cell}_{dim}_PPMI_Control.csv'
        GO_embedd_Cancer  = f'./Data/_GO_Embeddings_GO_PPI_{cancer}_{dim}_PPMI_Cancer.csv'
        
        # Top moving GO terms:
        
        ranking_p    = f'{moving_path}Rank_movement_{tissue}_PPMI.csv'
        Ranking_db = pd.read_csv(ranking_p,  index_col = 0)
        
        # Load and preprocess the data:
        
        # Cancer:
        
        Gene_embedd_Cancer_db       = pd.DataFrame(np.load(Gene_embedd_Cancer, allow_pickle=True))
        genes_Cancer_db             = list(pd.read_csv(genes_Cancer)["0"])
        genes_Cancer_db             = [str(i) for i in genes_Cancer_db]
        Gene_embedd_Cancer_db.index = genes_Cancer_db
        GO_embedd_Cancer_db         = pd.read_csv(GO_embedd_Cancer,  index_col = 0)
        GO_embedd_Cancer_db[GO_embedd_Cancer_db < 0] = 0
        
        # Coordinates:
        
        Gene_Cancer_Coordinates = Gene_embedd_Cancer_db.dot(Gene_midle_Cancer_db)
        
        # Control:
        
        Gene_embedd_Control_db       = pd.DataFrame(np.load(Gene_embedd_Control, allow_pickle=True))
        genes_Control_db             = list(pd.read_csv(genes_Control)["0"])
        genes_Control_db             = [str(i) for i in genes_Control_db]
        Gene_embedd_Control_db.index = genes_Control_db
        GO_embedd_Control_db         = pd.read_csv(GO_embedd_Control,  index_col = 0)
        GO_embedd_Control_db[GO_embedd_Control_db < 0] = 0   
        
        # Coordinates:
        
        Gene_Control_Coordinates = Gene_embedd_Control_db.dot(Gene_midle_Control_db)
        
        # Subset GO by the top 100 moving:
        
        vector_Ranking_Cancer  = GO_embedd_Cancer_db[Ranking_db.index]
        vector_Ranking_Control = GO_embedd_Control_db[Ranking_db.index]
        
        # Calculate distances:
        
        Cancer_distances  = Distance_GO_Gene(vector_Ranking_Cancer,  Gene_Cancer_Coordinates)
        Control_distances = Distance_GO_Gene(vector_Ranking_Control, Gene_Control_Coordinates)
        
        # Add more information to the final data frame:
        
        Matrix_GO_Gene = pd.read_csv(f'./Data/_Matrix_Genes_GO_BP_PPI.csv', 
                                        index_col = 0, dtype={0: str}) 
        
        # Cancer:
        
        GO_Gene_Matrix_Cancer        = Matrix_GO_Gene[Matrix_GO_Gene.index.isin(genes_Cancer_db)]
        GO_Gene_Matrix_Cancer.index  = [str(i) for i in GO_Gene_Matrix_Cancer.index]   
        
        # Control: 
        
        GO_Gene_Matrix_Control        = Matrix_GO_Gene[Matrix_GO_Gene.index.isin(genes_Control_db)]
        GO_Gene_Matrix_Control.index  = [str(i) for i in GO_Gene_Matrix_Control.index]  
        
        # Compute the information
        
        Cancer_Cosine_it  = Generate_Final_DataFrame(GO_Gene_Matrix_Cancer,  Cancer_distances)
        Control_Cosine_it = Generate_Final_DataFrame(GO_Gene_Matrix_Control, Control_distances)
        
        # Save the data frame:
        
        Cancer_Cosine_it.to_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',   header = True, index=True)
        Control_Cosine_it.to_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', header = True, index=True)

def Distance_Check(GO_Gene_Matrix, distances):
    
    # Lists to fill:
    
    x  = [] # Annotated distance
    y  = [] # Non-Annotated distance

    # Iteration per gene:

    for gene in list(distances.index):
  
        it_gene       = distances.loc[gene]
        it_annotation = GO_Gene_Matrix.loc[gene]
        
        for annotation in list(it_gene.index):
            
            if it_annotation[annotation] == 1:
                x.append(it_gene[annotation])
            else:
                y.append(it_gene[annotation])
    
    # Return the dataFrame:
    
    return(x,y)


def Demonstrate_Gene_GO_Org(Cancer_list, Normal_tissue_list, Cell_type_list, dim):
    
    '''
    This function test if the embedding vectors of the genes are close to the embedding vectors
    of those annotations that annotate them.
    '''
    
    # Paths:

    network_path = "./Data/"
    gene_GO_path = "./Data/"
    matrix_path  = "./Data/"
    
    # Open the file:
    
    with open(gene_GO_path + "Organization_Common_Space.txt", 'a') as the_file:
        
        # Write the columns names in the file:
        
        the_file.write("# Sample"   + "\t")
        the_file.write("# Annotated_NonAnnotated(less distance)"   + "\t")
        the_file.write("# Distance Annotated"   + "\t")
        the_file.write("# Distance Not-Annotated"   + "\n")


        # For each combination:
        
        for cancer, tissue, cell in zip(Cancer_list,Normal_tissue_list, Cell_type_list):
            
            # Load the data:
            
            # Genes:
            
            genes_Control =  f'{network_path}Control_{tissue}_{cell}_Genes.csv'
            genes_Cancer  =  f'{network_path}Cancer_{cancer}_Genes.csv'
            
            genes_Cancer_db  = list(pd.read_csv(genes_Cancer)["0"])
            genes_Cancer_db  = [str(i) for i in genes_Cancer_db]
            genes_Control_db = list(pd.read_csv(genes_Control)["0"])
            genes_Control_db = [str(i) for i in genes_Control_db]
            
            # Annotations:
            
            Matrix_GO_Gene = pd.read_csv(f'{matrix_path}_Matrix_Genes_GO_BP_PPI.csv', 
                                            index_col = 0, dtype={0: str}) 
            # Cancer annotations:
            
            GO_Gene_Matrix_Cancer        = Matrix_GO_Gene[Matrix_GO_Gene.index.isin(genes_Cancer_db)]
            GO_Gene_Matrix_Cancer.index  = [str(i) for i in GO_Gene_Matrix_Cancer.index]   
            
            # Control annotations: 
            
            GO_Gene_Matrix_Control        = Matrix_GO_Gene[Matrix_GO_Gene.index.isin(genes_Control_db)]
            GO_Gene_Matrix_Control.index  = [str(i) for i in GO_Gene_Matrix_Control.index] 
            
            # Distances:
             
            Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
            Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
            
            Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
            Control_Cosine_it.index  = [str(i) for i in Control_Cosine_it.index] 
            
            Cancer_Cosine_it  = Cancer_Cosine_it.drop( ['Group', 'Shape', 'size', 'Min_Dist', 'Min_GO'], axis=1)
            Control_Cosine_it = Control_Cosine_it.drop(['Group', 'Shape', 'size', 'Min_Dist', 'Min_GO'], axis=1)
            
            # Compute the sets:
            
            x_cancer,  y_cancer  = Distance_Check(GO_Gene_Matrix_Cancer,  Cancer_Cosine_it)
            x_control, y_control = Distance_Check(GO_Gene_Matrix_Control, Control_Cosine_it)
        
            # Compare distance distributions:
        
            p_value_cancer  = mannwhitneyu(x_cancer,  y_cancer,   alternative = "less").pvalue
            p_value_control = mannwhitneyu(x_control, y_control,  alternative = "less").pvalue
            
            # Write into the file:
            
            the_file.write(f'Cancer {tissue}\t')
            the_file.write(f'{p_value_cancer}\t')
            the_file.write(f'{np.mean(x_cancer)} ({np.std(x_cancer)})\t')
            the_file.write(f'{np.mean(y_cancer)} ({np.std(y_cancer)})\n')
        
            the_file.write(f'Control {tissue}\t')
            the_file.write(f'{p_value_control}\t')  
            the_file.write(f'{np.mean(x_control)} ({np.std(x_control)})\t')
            the_file.write(f'{np.mean(y_control)} ({np.std(y_control)})\n')
            
        # Close the file:
        
        the_file.close()

def Fold_enriched(Total, Total_succes, sample, sample_succes):
    
    expected = (sample*Total_succes)/Total
    
    if sample_succes >= expected:
        print("Over")
        fold     = sample_succes/expected
        p_value  = 1 - ss.hypergeom.cdf(sample_succes-1, Total, Total_succes, sample)
        return(fold, p_value)
    else:
        print("Under")
        fold     = expected/sample_succes
        p_value  = 1 - ss.hypergeom.cdf(sample_succes-1, Total, Total_succes, sample)
        return(fold, 1-p_value)


def Do_Fold_Enrichment_Complete(Normal_tissue_list):
    
    '''
    This function generate the gene predictions based on the movement of the genes to the
    set of shifted annotations. In addition, it also validates the set by performing a fold
    enrichment analyses.
    '''
     
    path_rank    = "./Data/"
    gene_GO_path = "./Data/"
    
    # Open the file:
    
    with open(f'{gene_GO_path}Fold_Rank_Table_Dos.txt', 'a') as the_file:
        
        # Write the name of the columns:
        
        the_file.write(f'Cancer Type\t')
        the_file.write(f'Moving Genes\t')
        the_file.write(f'Stable Genes\n')
        
        for tissue in Normal_tissue_list:
            
            # Load movement and choose the shifted annotations:
            
            rank      = pd.read_csv(f'{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv')      
            moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"]

            # Load the bibliography information:
                
            counter_list = pd.read_csv(f'{gene_GO_path}Cancer_Count_{tissue}.csv')
            counter_list = counter_list.drop_duplicates("name")
                
            # Load cosines:
                
            Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
            Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
              
            Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
            Control_Cosine_it.index = [str(i) for i in Control_Cosine_it.index]
            
            # For all the GO terms.
            
            # Subset the data:
            
            gene_common           = Cancer_Cosine_it.index[Cancer_Cosine_it.index.isin(Control_Cosine_it.index)]
            Cancer_Cosine_it_com  = Cancer_Cosine_it.loc[gene_common]
            Control_Cosine_it_com = Control_Cosine_it.loc[gene_common]
            
            moving_GO = list(set(Cancer_Cosine_it_com.columns).intersection(set(moving_GO)))
            
            # Filter the moving GO terms (based on the distribution)
            
            Cancer_structure_GO1  = Cancer_Cosine_it_com[moving_GO]
            Control_structure_GO1 = Control_Cosine_it_com[moving_GO]
            
            # Analyze the movement of the genes together:
            
            move = Control_structure_GO1 - Cancer_structure_GO1
                   
            # Get the sets to analyze:
            
            # A) Moving Vs Not moving (maximum movement):
            
            move_abs            = abs(move)
            Restuls_abs         = pd.DataFrame(move_abs.T.max())
            Restuls_abs["GO"]   = move_abs.T.idxmax()
            Restuls_abs["Gene"] = list(Restuls_abs.index)
            Restuls_abs         = Restuls_abs.sort_values(by = 0,ascending = False)
            Restuls_abs = Restuls_abs.reset_index(drop = True)
            
            # Normalize the distribution by using a sqrt transformation:        
            
            Restuls_abs["Original"] = Restuls_abs[0]
            z = [np.log(i) for i in list(Restuls_abs[0])]
            Restuls_abs[0] = z 
            
            # Choose the tails:
        
            Restuls_abs_top  = Restuls_abs.head(round(len(Restuls_abs)*5/100)) # 5th and 95th percentiles
            Restuls_abs_tail = Restuls_abs.tail(round(len(Restuls_abs)*5/100))
             
            # To save the predictions:
            
            gene_names = counter_list[counter_list.initial_alias.isin(Restuls_abs_top.Gene)] 
            gene_names = gene_names.drop(["Unnamed: 0", "initial_alias"], axis=1) 
            gene_names = gene_names.drop_duplicates("name")
            gene_names = gene_names.reset_index(drop = True)
            
            gene_names.to_csv(f'{path_rank}Predictions_Rank_{tissue}.csv', header = True, index=True)
            
            # Get the total set:

            genes_succes = list(counter_list[counter_list.Counts > 0]["initial_alias"])
            genes_succes = [str(i) for i in genes_succes]
                                
            Total_succes = len(genes_succes)
            Total        = len(counter_list)
            
            # Do the enrichment analyses:
            
            set_stable_succes = len(set(Restuls_abs_tail.Gene).intersection(set(genes_succes)))
            set_stable_total  = len(Restuls_abs_tail)
            
            set_moving_succes = len(set(Restuls_abs_top.Gene).intersection(set(genes_succes)))
            set_moving_total  = len(Restuls_abs_top)        
            
            fold_stb, p_value_stb = Fold_enriched(Total, Total_succes, set_stable_total, set_stable_succes) 
            fold_mv, p_value_mv   = Fold_enriched(Total, Total_succes, set_moving_total, set_moving_succes) 
                        
            the_file.write(f'{tissue}\t')
            the_file.write(f'{round(fold_mv, 3)} ({p_value_mv})\t')
            the_file.write(f'{round(fold_stb, 3)} ({p_value_stb})\n')
        
        # Close the file:
        
        the_file.close()                
            


















        