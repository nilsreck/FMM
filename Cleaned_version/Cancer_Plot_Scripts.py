# This Script Contains the fucntions for the Cancer Study:

# Packages:

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
from sklearn import metrics 
from ast import literal_eval
from matplotlib import lines
from selenium import webdriver
from collections import Counter
import matplotlib.pyplot as plt
from multiprocessing import Pool
from pySankey.sankey import sankey
from itertools import combinations, cycle
from scipy.spatial import distance
from string import ascii_lowercase
import matplotlib.patches as mpatches
from pygosemsim import graph, similarity
from sklearn_extra.cluster import KMedoids
from matplotlib.ticker import NullFormatter
from matplotlib_venn import venn2_unweighted
from statsmodels.stats.multitest import multipletests
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.stats import mannwhitneyu, pearsonr, gaussian_kde
from sklearn.metrics import pairwise_distances, silhouette_score
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from scipy.stats import hypergeom

# Our packages:

os.chdir("./Cleaned_version/")

import generate_embedding_spaces_Human 
import generate_embedding_spaces
import enrichment_analysis
import GO_Embeddings

# External package:

import clusterp


# Functions:

def Analysis_Dendograms(Normal_Tissue_list):
    
    network_path = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    Filter_path  = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    save_path    = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    for tissue in Normal_Tissue_list :
        
        # Load the FMM:
        
        Cancer_FMM  = pd.read_csv(f'{network_path}_GO_Embeddings_GO_PPI_Cancer_{tissue}.csv', index_col = 0)
        Control_FMM = pd.read_csv(f'{network_path}_GO_Embeddings_GO_PPI_Control_{tissue}.csv', index_col = 0)
        
        # Load the common set of annotations (to speed up the process):
               
        common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}.csv')["0"])
        
        # Filter the FMMs:
        
        Cancer_FMM  = Cancer_FMM[common_set]
        Control_FMM = Control_FMM[common_set]
        
        Cancer_FMM[Cancer_FMM < 0] = 0
        Control_FMM[Control_FMM < 0] = 0
        
        # Start the computations:
        
        print("Starting with control")
        pv_control = clusterp.PvClust(Control_FMM, method="average", metric="cosine", nboot=1000, parallel=True)
        print("Starting with cancer")
        pv_cancer  = clusterp.PvClust(Cancer_FMM,  method="average", metric="cosine", nboot=1000, parallel=True)
        
        # Save the plots:
        
        pv_control.seplot(save_path, f'Control_{tissue}')
        pv_cancer.seplot(save_path, f'Cancer_{tissue}') 
        
        # Save the information of the p-values in a table:
        
        pv_control = pv_control.print_result()
        pv_cancer  = pv_cancer.print_result()
        
        pv_control.to_csv(f'{save_path}Table_Dendogram_{tissue}_Control.csv')
        pv_cancer.to_csv(f'{save_path}Table_Dendogram_{tissue}_Cancer.csv')
        
        print(f'Finished with {tissue}')
      
        
def Analysis_Dendograms_Results(Normal_Tissue_list):        
        
    save_path    = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    for tissue in Normal_Tissue_list:   
        
        # Load the information:
        
        control = pd.read_csv(f'{save_path}Table_Dendogram_{tissue}_Control.csv')
        cancer  = pd.read_csv(f'{save_path}Table_Dendogram_{tissue}_Cancer.csv')
        
        # Show results:
        
        print(f'Control {tissue}: {sum(control["AU"] > 0.95)/len(control)}% Axes: {sum(control["AU"] > 0.95)}')
        print(f'Cancer {tissue}: {sum(cancer["AU"] > 0.95)/len(cancer)}% Axes: {sum(cancer["AU"] > 0.95)}')
    

def Comparison_Gene_GO_Enrichments(Cancer_type_list, Normal_Tissue_list, Cell_type_list, annotation = "terms"):
    
   # Networks path:
    
    network_path = "/media/sergio/sershiosdisk/Cancer/Networks/"
    
    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    # Annotations:

    BP = json.load(open("/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPI_GO_BackProp_biological_process.json"))
    
    # Set a dictionary:
    
    list_keys = []
    
    for i in Normal_Tissue_list:      
        list_keys.append(f'Cancer_{i}')
        list_keys.append(f'Control_{i}')
    
    Final = dict.fromkeys(list_keys)
    
    # Do the enrichments and save the value: 
   
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):

        for dim in [200]:
            
            # Load the name of the genes:
            
            path_Cancer  = f'{network_path}Cancer_{cancer}_Genes.csv'
            path_Control = f'{network_path}Control_{tissue}_{cell}_Genes.csv'
        
            Cancer_genes  = list(pd.read_csv(path_Cancer)["0"].astype(str))
            Control_genes = list(pd.read_csv(path_Control)["0"].astype(str))
            
            # Load the G Matrices:
            
            # Cancer:
            
            P_PPMI_Cancer = np.load(f'{network_path}_P_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)            
            U_PPMI_Cancer = np.load(f'{network_path}_U_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True) 
            
            P_PPMI_Control = np.load(f'{network_path}_P_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True)
            U_PPMI_Control = np.load(f'{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True) 
            
            # Compute the coordinates of the genes:
            
            G_PPMI_Cancer = P_PPMI_Cancer.dot(U_PPMI_Cancer)
            G_PPMI_Control = P_PPMI_Control.dot(U_PPMI_Control) 
            
            # Do the clustering:
            
            G_Cancer_PPMI_Cluster  = KMedoids_format(G_PPMI_Cancer, Cancer_genes) 
            G_Control_PPMI_Cluster = KMedoids_format(G_PPMI_Control, Control_genes) 
          
            # Compute the Enrichments:
            
            Enrichments_Control  = enrichment_analysis.enrichment_analysis(G_Control_PPMI_Cluster, BP, dim, save, 0.05, annotation, f'PPI_{cancer}', "PPMI", "terms")
            Enrichments_Cancer   = enrichment_analysis.enrichment_analysis(G_Cancer_PPMI_Cluster,  BP, dim, save, 0.05, annotation, f'PPI_{cancer}', "PPMI", "terms")
            
            # Append the results to the dictionary:
                       
            Final[f'Control_{tissue}'] = list(Enrichments_Control)
            Final[f'Cancer_{tissue}']  = list(Enrichments_Cancer) 
            
    # Saver the results:
    
    json_id = json.dumps(Final)
    f       = open(f'{save}GO_Enrichments.json',"w")
    f.write(json_id)

    # close file
    f.close()


def Fit_Data_Frame(clusters):

    # Do the Final DataFrame:
        
    set_GO = list(set([item for sublist in clusters for item in sublist]))
    Final  = pd.DataFrame(0, index=set_GO,  columns=set_GO)
    
    for clust in clusters:        
        Final.loc[clust, clust] = 1   
    return(Final)
        

def Generate_Matrix_Enrichmetns(Normal_tissue_list):
    
    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/" 
    
    # Load the enrichments:
    
    f           = open (f'{save}GO_Enrichments.json', "r")
    enrichments = json.loads(f.read())
    f.close()
    
    # Generate the matrices:
            
    for tissue in Normal_tissue_list:
        
        # Choose the corresponding tissue:
        
        Control = enrichments[f'Control_{tissue}']
        Cancer  = enrichments[f'Cancer_{tissue}']
        
        Control = [x for x in Control if x]
        Cancer  = [x for x in Cancer if x]
              
        # Fit the dataframe
        
        Final_Control  = Fit_Data_Frame(Control)
        Final_Cancer   = Fit_Data_Frame(Cancer)
        
        # Save the information:
        
        Final_Control.to_csv(f'{save}Matrix_Enrichments_Control_{tissue}.csv')
        Final_Cancer.to_csv(f'{save}Matrix_Enrichments_Cancer_{tissue}.csv')
 

def Plot_Differences(Normal_tissue_list):
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    for tissue in Normal_tissue_list:
        
        # Load the enrichments and the FMMs:
        
        Cancer_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Cancer_{tissue}.csv', index_col = 0)
        Cancer_Embeddings  = pd.read_csv(f'{save}Cosine_Cancer_{tissue}.csv', index_col = 0)
        
        Control_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Control_{tissue}.csv', index_col = 0)
        Control_Embeddings  = pd.read_csv(f'{save}Cosine_Control_{tissue}.csv', index_col = 0)
        
        # Keep only those GO terms that are enriched for the comparison:
        
        FMM_Cancer  = Cancer_Embeddings.loc[Cancer_Enrichments.index, Cancer_Enrichments.index]
        FMM_Control = Control_Embeddings.loc[Control_Enrichments.index, Control_Enrichments.index]
        
        # Transform the FMM to a similarity:
        
        FMM_Cancer  = 1 - FMM_Cancer
        FMM_Control = 1 - FMM_Control
        
        # Do the plot:
        
        diff_Cancer  = FMM_Cancer  - Cancer_Enrichments
        diff_Control = FMM_Control - Control_Enrichments 
        
        diff_Control = diff_Control.mask(np.triu(np.ones(diff_Control.shape)).astype(bool)).stack().sort_values(ascending=True)
        diff_Cancer  = diff_Cancer.mask(np.triu(np.ones(diff_Cancer.shape)).astype(bool)).stack().sort_values(ascending=True)
                
        plt.style.use("seaborn-whitegrid")
        plt.rcParams.update({'font.size': 17})
        plt.rc('font', weight='bold')
        
        fig, axes = plt.subplots(1, 1, figsize=(5,2.5))
        fig.tight_layout(pad = 0.5)

        axes.hist(diff_Cancer, bins = 20) 
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['left'].set_linewidth(3)
        axes.spines['bottom'].set_linewidth(3)
        
        axes.set_ylabel("Pairs", fontsize=20, fontweight='bold')
        axes.set_xlabel("Difference", fontsize=20, fontweight='bold')
        
        axes.spines['left'].set_color("grey")
        axes.spines['bottom'].set_color("grey") 
        
        fig.savefig(f'{save}Distribution_Cancer_{tissue}.svg', format="svg", dpi=600,
                    bbox_inches='tight')
        
        # Control:
               
        plt.style.use("seaborn-whitegrid")
        plt.rcParams.update({'font.size': 17})
        plt.rc('font', weight='bold')
        
        fig, axes = plt.subplots(1, 1, figsize=(5,2.5))
        fig.tight_layout(pad = 0.5)

        axes.hist(diff_Control, bins = 20) 
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['left'].set_linewidth(3)
        axes.spines['bottom'].set_linewidth(3)
        
        axes.set_ylabel("Pairs", fontsize=20, fontweight='bold')
        axes.set_xlabel("Difference", fontsize=20, fontweight='bold')
        
        axes.spines['left'].set_color("grey")
        axes.spines['bottom'].set_color("grey") 
        
        fig.savefig(f'{save}Distribution_Control_{tissue}.svg', format="svg", dpi=600,
                    bbox_inches='tight')                 


def Show_Functional_Interactions_Enrichments(Normal_tissue_list):
    
    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    # For each tissue:
    
    for tissue in Normal_tissue_list:
        
        # Load the enrichments and the FMMs:
        
        Cancer_Enrichments  = pd.read_csv(f'{save}Matrix_Enrichments_Cancer_{tissue}.csv', index_col = 0)
        Control_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Control_{tissue}.csv', index_col = 0)
        
        # Take the upper diagonal matrix of each:
        
        # Fit diagonal with 0s:
        
        Cancer_Enrichments.values[[np.arange(Cancer_Enrichments.shape[0])]*2]   = 0
        Control_Enrichments.values[[np.arange(Control_Enrichments.shape[0])]*2] = 0

        # Take the diagonal:        

        diag_Enrch_Cancer  = list(Cancer_Enrichments.values[np.triu_indices(len(Cancer_Enrichments),   k = 1)])
        diag_Enrch_Control = list(Control_Enrichments.values[np.triu_indices(len(Control_Enrichments),   k = 1)])
        
        # Show the number of interactions:
        
        print(f'Functional_Interactions_Cancer_{tissue}  = {sum(diag_Enrch_Cancer)}')
        print(f'Functional_Interactions_Control_{tissue} = {sum(diag_Enrch_Control)}')        

def GO_Clusters_Matrix(clusters, distance):
    
    # Init empty matrix:
    
    Matrix       = pd.DataFrame(0, index=np.arange(len(distance.columns)), columns=distance.columns)
    Matrix.index = distance.columns
    
    # Fill the Matrix:
    
    GO_terms = list(Matrix.index)
    
    for cluster in set(clusters.labels_):
        
        cluster_list       = np.where(clusters.labels_ == cluster)[0]
        GO_terms_clustered = [GO_terms[i] for i in cluster_list]
        
        Matrix.loc[GO_terms_clustered, GO_terms_clustered] = 1
    
    return(Matrix)


def Connections(Normal_tissue_list):

    network_path = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    Filter_path  = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    
    for tissue in Normal_tissue_list:
        
        # Load info (here are the p-values that we computed over 100k):
    
        dendogram_cancer  = pd.read_csv(f'/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/Table_Dendogram_{tissue}_Cancer.csv')
        dendogram_control = pd.read_csv(f'/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/Table_Dendogram_{tissue}_Control.csv')
        
        Cancer_FMM  = pd.read_csv(f'{network_path}_GO_Embeddings_GO_PPI_Cancer_{tissue}.csv', index_col = 0) 
        Control_FMM = pd.read_csv(f'{network_path}_GO_Embeddings_GO_PPI_Control_{tissue}.csv', index_col = 0)
        
        common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}.csv')["0"])
        
        Cancer_FMM  = Cancer_FMM[common_set]
        Control_FMM = Control_FMM[common_set]
        
        Cancer_FMM[Cancer_FMM < 0]   = 0
        Control_FMM[Control_FMM < 0] = 0
        
        # Do the hierarchical (just for the clusters, is deterministic):
    
        pv_cancer   = clusterp.PvClust(Cancer_FMM,   method="average", metric="cosine", nboot=1, parallel=True) 
        pv_control  = clusterp.PvClust(Control_FMM,  method="average", metric="cosine", nboot=1, parallel=True)  
        
        AU_list_cancer   = list(dendogram_cancer.AU)
        AU_list_control  = list(dendogram_control.AU)
        
        clusters_cancer  = pv_cancer.clusters
        clusters_control = pv_control.clusters
        
        # Get the group of GO terms in the signf cluster:
        
        clusters_sigf_cancer  = [clusters_cancer[i]  for i in range(len(clusters_cancer))  if AU_list_cancer[i] > 0.95]
        clusters_sigf_control = [clusters_control[i] for i in range(len(clusters_control)) if AU_list_control[i] > 0.95]
        
        del clusters_sigf_cancer[-1]
        del clusters_sigf_control[-1]
           
        GO_cancer  = list(Cancer_FMM.columns)
        GO_control = list(Control_FMM.columns)
        
        final_cancer = []
        
        for clust in clusters_sigf_cancer:
            
            final_cancer.append([GO_cancer[i] for i in clust])
            
        final_control = []
        
        for clust in clusters_sigf_control:
            
            final_control.append([GO_control[i] for i in clust]) 
        
        # Fill a Matrix to do the recount:
                
        Matrix_cluster_cancer       = pd.DataFrame(0, index=np.arange(len(Cancer_FMM.columns)), columns=Cancer_FMM.columns)
        Matrix_cluster_cancer.index = Cancer_FMM.columns
        
        Matrix_cluster_control       = pd.DataFrame(0, index=np.arange(len(Control_FMM.columns)), columns=Control_FMM.columns)
        Matrix_cluster_control.index = Control_FMM.columns
        
        for cluster in final_cancer:
            
            Matrix_cluster_cancer.loc[cluster,cluster] = 1
            
        Matrix_cluster_cancer.values[[np.arange(Matrix_cluster_cancer.shape[0])]*2] = 0
        Matrix_cluster_cancer = list(Matrix_cluster_cancer.values[np.triu_indices(len(Matrix_cluster_cancer),   k = 1)])

        for cluster in final_control:
            
            Matrix_cluster_control.loc[cluster,cluster] = 1
            
        Matrix_cluster_control.values[[np.arange(Matrix_cluster_control.shape[0])]*2] = 0
        Matrix_cluster_control = list(Matrix_cluster_control.values[np.triu_indices(len(Matrix_cluster_control),   k = 1)])
        
        print(f'Functional_Interactions_Cancer_{tissue}  = {sum(Matrix_cluster_cancer)}')
        print(f'Functional_Interactions_Control_{tissue}  = {sum(Matrix_cluster_control)}')    
    
    
def Show_Functional_Interactions_GO_Embeddings(Normal_tissue_list):
    
    # Paths:

    save_cosine = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    path        = "/media/sergio/sershiosdisk/Cancer/Networks/"
    
    # Load GO/gene matrix:
        
    GO_Matrix = pd.read_csv(f'{path}_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv', 
                                    index_col = 0, dtype={0: str}) 
    GO_Matrix.index = GO_Matrix.index.astype(str)
        
    # For each cancer and tissue:
        
    for tissue in Normal_tissue_list:
        
        # Load gene matrrix:
        
        Cancer_genes       = pd.read_table(f'{save_cosine}Cancer_{tissue}_Genes.csv', header = 0, dtype={0: str})
        Cancer_genes_list  = list(Cancer_genes["0"])
        Control_genes      = pd.read_table(f'{save_cosine}Control_{tissue}_Genes.csv', header = 0, dtype={0: str})
        Control_genes_list = list(Control_genes["0"])
        
        # Subset the gene/GO matrix:
        
        GO_Gene_Matrix_Cancer  = GO_Matrix[GO_Matrix.index.isin(Cancer_genes_list)]
        GO_Gene_Matrix_Control = GO_Matrix[GO_Matrix.index.isin(Control_genes_list)]
        
        GO_Gene_Matrix_Cancer  = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
        GO_Gene_Matrix_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]
        
        # Get the informative GO terms:
        
        # Cancer:
    
        GO_Cancer      = GO_Gene_Matrix_Cancer.sum(axis=0)
        GO_Cancer_filt = set(GO_Cancer[GO_Cancer >= 3].index)
        
        # Control:
        
        GO_Control      = GO_Gene_Matrix_Control.sum(axis=0)
        GO_Control_filt = set(GO_Control[GO_Control >= 3].index)
   
        # Load the FMM matrix:
        
        Cancer_Embeddings  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{tissue}.csv',  index_col = 0)
        Control_Embeddings = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}.csv', index_col = 0)
        
        Cancer_Embeddings  = Cancer_Embeddings.loc[GO_Cancer_filt,   GO_Cancer_filt]
        Control_Embeddings = Control_Embeddings.loc[GO_Control_filt, GO_Control_filt]
           
        # K medoids:
                   
        n_clusters_cancer  = round(math.sqrt(len(Cancer_Embeddings.index)/2))
        n_clusters_control = round(math.sqrt(len(Control_Embeddings.index)/2))
        
        grouping_1_cancer  = KMedoids(n_clusters=n_clusters_cancer,  metric='precomputed',).fit(Cancer_Embeddings)
        grouping_2_control = KMedoids(n_clusters=n_clusters_control, metric='precomputed',).fit(Control_Embeddings)
              
        # Matrix of clustering:
        
        Cancer_clust  = GO_Clusters_Matrix(grouping_1_cancer,  Cancer_Embeddings)
        Control_clust = GO_Clusters_Matrix(grouping_2_control, Control_Embeddings)

        # Fit diagonal with 0s:
        
        Cancer_clust.values[[np.arange(Cancer_clust.shape[0])]*2]   = 0
        Control_clust.values[[np.arange(Control_clust.shape[0])]*2] = 0

        # Take the diagonal:        

        diag_Enrch_Cancer  = list(Cancer_clust.values[np.triu_indices(len(Cancer_clust),     k = 1)])
        diag_Enrch_Control = list(Control_clust.values[np.triu_indices(len(Control_clust),   k = 1)])
            
        # Show the number of interactions:
        
        print(f'Functional_Interactions_Cancer_{tissue}  = {sum(diag_Enrch_Cancer)}')
        print(f'Functional_Interactions_Control_{tissue} = {sum(diag_Enrch_Control)}')

def Noel_exp(FMM, enrich):
    
    
    dist_enrich = FMM.values[enrich.values == 1]
    
    # Get the distribution of values of those terms that are not enriched:
    
    dist_no   = FMM.values[enrich.values == 0] 
    
    # Compute the statistical test:
    
    print(mannwhitneyu(list(dist_enrich), list(dist_no), alternative = "less").pvalue)

def Test_Reviwer(Normal_tissue_list):            
            
    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    # Get the information needed:
    
    list_keys = []
    
    for i in Normal_tissue_list:      
        list_keys.append(f'Cancer_{i}')
        list_keys.append(f'Control_{i}')
        
    
    Final_fpr = dict.fromkeys(list_keys)
    Final_tpr = dict.fromkeys(list_keys)
    Final_thr = dict.fromkeys(list_keys)    


    for tissue in Normal_tissue_list:
        
        # Load the enrichments and the FMMs:
        
        Cancer_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Cancer_{tissue}.csv', index_col = 0)
        Cancer_Embeddings  = pd.read_csv(f'{save}Cosine_Cancer_{tissue}.csv', index_col = 0)
        
        Control_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Control_{tissue}.csv', index_col = 0)
        Control_Embeddings  = pd.read_csv(f'{save}Cosine_Control_{tissue}.csv', index_col = 0)
        
        # Keep only those GO terms that are enriched for the comparison:
        
        FMM_Cancer  = Cancer_Embeddings.loc[Cancer_Enrichments.index, Cancer_Enrichments.index]
        FMM_Control = Control_Embeddings.loc[Control_Enrichments.index, Control_Enrichments.index]
        
        # Do the Noel experiment:
        
        Noel_exp(FMM_Cancer, Cancer_Enrichments)
        Noel_exp(FMM_Control, Control_Enrichments)
        
        # Transform the FMM to a similarity:
        
        FMM_Cancer  = 1 - FMM_Cancer
        FMM_Control = 1 - FMM_Control
                

        
        
        
        
        # Take the upper diagonal matrix of each:
        
        # Fit diagonal with 0s:
        
        Cancer_Enrichments.values[[np.arange(Cancer_Enrichments.shape[0])]*2] = 0
        FMM_Cancer.values[[np.arange(FMM_Cancer.shape[0])]*2] = 0
        
        Control_Enrichments.values[[np.arange(Control_Enrichments.shape[0])]*2] = 0
        FMM_Control.values[[np.arange(FMM_Control.shape[0])]*2] = 0
        
        # Take the diagonal:        
        
        diag_FMM_Cancer   = list(FMM_Cancer.values[np.triu_indices(len(FMM_Cancer),   k = 1)])
        diag_Enrch_Cancer = list(Cancer_Enrichments.values[np.triu_indices(len(Cancer_Enrichments),   k = 1)])

        diag_FMM_Control   = list(FMM_Control.values[np.triu_indices(len(FMM_Control),   k = 1)])
        diag_Enrch_Control = list(Control_Enrichments.values[np.triu_indices(len(Control_Enrichments),   k = 1)])
        
        # Calculate the AUROC:
        
        fpr_Cancer,  tpr_Cancer,  thresholds_Cancer  = metrics.roc_curve(diag_Enrch_Cancer,  diag_FMM_Cancer,  pos_label=1)
        fpr_Control, tpr_Control, thresholds_Control = metrics.roc_curve(diag_Enrch_Control, diag_FMM_Control, pos_label=1)
        
        print(metrics.auc(fpr_Cancer, tpr_Cancer))
        print(metrics.auc(fpr_Control, tpr_Control))
        
        Final_fpr[f'Cancer_{tissue}'] = list(fpr_Cancer)
        Final_tpr[f'Cancer_{tissue}'] = list(tpr_Cancer)
        Final_thr[f'Cancer_{tissue}'] = list(thresholds_Cancer)
        
        Final_fpr[f'Control_{tissue}'] = list(fpr_Control)
        Final_tpr[f'Control_{tissue}'] = list(tpr_Control)
        Final_thr[f'Control_{tissue}'] = list(thresholds_Control)   
    
    # Save the results:
    
    json_id = json.dumps(Final_fpr)
    f       = open(f'{save}FPR_Human.json',"w")
    f.write(json_id)
    f.close()    
    
    json_id = json.dumps(Final_tpr)
    f       = open(f'{save}TPR_Human.json',"w")
    f.write(json_id)
    f.close() 

    json_id = json.dumps(Final_thr)
    f       = open(f'{save}THR_Human.json',"w")
    f.write(json_id)
    f.close() 

            
def Plot_AUROC(Normal_tissue_list):            

    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"

    # Load the results previously shown:    

    f           = open (f'{save}FPR_Human.json', "r")
    fpr         = json.loads(f.read())
    f.close()    
    
    f           = open (f'{save}TPR_Human.json', "r")
    tpr         = json.loads(f.read())
    f.close() 

    f           = open (f'{save}THR_Human.json', "r")
    thr         = json.loads(f.read())
    f.close() 
    
    # Do the plots:
    
    fpr_cancer = []
    tpr_cancer = []
    thr_cancer = []
    
    auroc_cancer = []
   
    fpr_control = []
    tpr_control = []
    thr_control = []
    
    auroc_control = []
    
    for tissue in Normal_tissue_list:
        
        fpr_cancer.append(fpr[f'Cancer_{tissue}'])
        tpr_cancer.append(tpr[f'Cancer_{tissue}'])
        thr_cancer.append(thr[f'Cancer_{tissue}'])
        
        fpr_control.append(fpr[f'Cancer_{tissue}'])
        tpr_control.append(tpr[f'Cancer_{tissue}'])
        thr_control.append(thr[f'Cancer_{tissue}'])
        
        auroc_cancer.append(metrics.auc(fpr[f'Cancer_{tissue}'],   tpr[f'Cancer_{tissue}']))
        auroc_control.append(metrics.auc(fpr[f'Control_{tissue}'],  tpr[f'Control_{tissue}']))
        
    # Do the plot:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
        
    fig, axes = plt.subplots(1, 1, figsize=(10,5))
    fig.tight_layout(pad = 0.5)
    
    colors = cycle(["aqua", "darkorange", "cornflowerblue","brown"])

    count = 0
    lw = 2
    for i, color in zip(range(4), colors):
        axes.plot(
            fpr_cancer[i],
            tpr_cancer[i],
            color=color,
            lw=lw,
            label=f'AUROC {Normal_tissue_list[count]} cancer ' + "({1:0.2f})".format(i, auroc_cancer[i]),)  
        count= count+1
        
    axes.plot([0, 1], [0, 1], "k--", lw=lw)
    axes.set_xlim([0.0, 1.0])
    axes.set_ylim([0.0, 1.05])
    axes.set_xlabel("False Positive Rate")
    axes.set_ylabel("True Positive Rate")
    axes.legend(loc="lower right")

    fig.savefig(f'{save}ROC_Cancer.svg', format="svg", dpi=600,
                    bbox_inches='tight') 
           
    # Do the plot:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
        
    fig, axes = plt.subplots(1, 1, figsize=(10,5))
    fig.tight_layout(pad = 0.5)
    
    colors = cycle(["aqua", "darkorange", "cornflowerblue","brown"])

    count = 0
    lw = 2
    for i, color in zip(range(4), colors):
        axes.plot(
            fpr_control[i],
            tpr_control[i],
            color=color,
            lw=lw,
            label=f'AUROC {Normal_tissue_list[count]} control ' + "({1:0.2f})".format(i, auroc_control[i]),)  
        count=count+1
        
    axes.plot([0, 1], [0, 1], "k--", lw=lw)
    axes.set_xlim([0.0, 1.0])
    axes.set_ylim([0.0, 1.05])
    axes.set_xlabel("False Positive Rate")
    axes.set_ylabel("True Positive Rate")
    axes.legend(loc="lower right")

    fig.savefig(f'{save}ROC_Control.svg', format="svg", dpi=600,
                    bbox_inches='tight')         
        
 
def Test_Reviwer_Random(Normal_tissue_list, repetitions = 1000):   
         
    # Save Path:
    
    save = "/media/sergio/sershiosdisk/Cancer/Gene_GO_Test/"
    
    # Get the information needed:
    
    list_keys = []
    
    for i in Normal_tissue_list:      
        list_keys.append(f'Cancer_{i}')
        list_keys.append(f'Control_{i}')
        
    # Result (Real ROC and the final results):
    
    Counter_pd      = pd.DataFrame(pd.DataFrame(0,  columns= ["Counter"], index = list_keys))
    Counter_pd_real = pd.DataFrame(pd.DataFrame([0.90, 0.86, 0.87, 0.91, 0.90, 0.88, 0.93, 0.84], columns= ["Counter"], index = list_keys))
      
    for tissue in Normal_tissue_list:
        
        print(tissue)
        
        # Load the enrichments and the FMMs:
        
        Cancer_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Cancer_{tissue}.csv', index_col = 0)
        Cancer_Embeddings  = pd.read_csv(f'{save}Cosine_Cancer_{tissue}.csv', index_col = 0)
        
        Control_Enrichments = pd.read_csv(f'{save}Matrix_Enrichments_Control_{tissue}.csv', index_col = 0)
        Control_Embeddings  = pd.read_csv(f'{save}Cosine_Control_{tissue}.csv', index_col = 0)
        
        # Keep only those GO terms that are enriched for the comparison:
        
        FMM_Cancer  = Cancer_Embeddings.loc[Cancer_Enrichments.index, Cancer_Enrichments.index]
        FMM_Control = Control_Embeddings.loc[Control_Enrichments.index, Control_Enrichments.index]
        
        # Transform the FMM to a similarity:
        
        FMM_Cancer  = 1 - FMM_Cancer
        FMM_Control = 1 - FMM_Control
        
        # Start the permutations:
        
        for rep in range(repetitions):
            
            if rep % 100 == 0:
                print(rep)
            
            # Randomize the experiment:
            
            Cancer_Enrichments_random  = Cancer_Enrichments.sample(frac=1).reset_index(drop=True)
            Control_Enrichments_random = Control_Enrichments.sample(frac=1).reset_index(drop=True)
            
            Cancer_Enrichments_random.index   = Cancer_Enrichments.index
            Cancer_Enrichments_random.columns = Cancer_Enrichments.columns
            
            Control_Enrichments_random.index   = Control_Enrichments.index
            Control_Enrichments_random.columns = Control_Enrichments.columns
            
            # Compute the ROC of the random:
        
            # Take the upper diagonal matrix of each:
            
            # Fit diagonal with 0s:
            
            Cancer_Enrichments_random.values[[np.arange(Cancer_Enrichments_random.shape[0])]*2] = 0
            FMM_Cancer.values[[np.arange(FMM_Cancer.shape[0])]*2] = 0
            
            Control_Enrichments_random.values[[np.arange(Control_Enrichments_random.shape[0])]*2] = 0
            FMM_Control.values[[np.arange(FMM_Control.shape[0])]*2] = 0
            
            # Take the diagonal:        
            
            diag_FMM_Cancer   = list(FMM_Cancer.values[np.triu_indices(len(FMM_Cancer),   k = 1)])
            diag_Enrch_Cancer = list(Cancer_Enrichments_random.values[np.triu_indices(len(Cancer_Enrichments_random),   k = 1)])
    
            diag_FMM_Control   = list(FMM_Control.values[np.triu_indices(len(FMM_Control),   k = 1)])
            diag_Enrch_Control = list(Control_Enrichments_random.values[np.triu_indices(len(Control_Enrichments_random),   k = 1)])
        
            # Calculate the AUROC:
            
            fpr_Cancer,  tpr_Cancer,  thresholds_Cancer  = metrics.roc_curve(diag_Enrch_Cancer,  diag_FMM_Cancer,  pos_label=1)
            fpr_Control, tpr_Control, thresholds_Control = metrics.roc_curve(diag_Enrch_Control, diag_FMM_Control, pos_label=1)
            
            Random_ROC_Cancer  = metrics.auc(fpr_Cancer, tpr_Cancer)
            Random_ROC_Control = metrics.auc(fpr_Control, tpr_Control)
            
            # Compare with the real:
            
            if  Random_ROC_Cancer >= Counter_pd_real.loc[f'Cancer_{tissue}'].values[0]:
                
                Counter_pd.loc[f'Cancer_{tissue}'] = Counter_pd.loc[f'Cancer_{tissue}'] + 1
            
            if  Random_ROC_Control  >=  Counter_pd_real.loc[f'Control_{tissue}'].values[0]:
                
                Counter_pd.loc[f'Control_{tissue}'] = Counter_pd.loc[f'Control_{tissue}'] + 1
            
                
    # Compute the p-values:
    
    Pvalues = (Counter_pd + 1)/(repetitions + 1)  

    # Save the p_values in a csv file:      
    
    Pvalues.to_csv(f'{save}Pvalues_ROC.csv')
    
    
    

def Venn_Diagram_Cancer_Networks(Cancer_type_list, Normal_Tissue_list, Cell_type_list):
    
    '''
    This fucntions plots the Venn diagram of the tissue specific Networks:
        
        Inputs:
                - Cancer_type   : string list, Name of the Cancer.
                - Normal_Tissue : string list, Name of the normal tissue.
                - Cell_type     : string list, Name of the normal cell in the tissue.
    '''
    
    # Networks path:
    
    network_path = "./Cleaned_version/Data/"
    
    # Plot Propieties:

    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    
    # Set the grid:
    
    fig, axes = plt.subplots(2, 2, figsize=(12,8))
    fig.tight_layout(pad = 0.5)
    
    row_count    = 0
    column_count = 0
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the data:
        
        path_Cancer  = f'{network_path}Cancer_{cancer}_Genes.csv'
        path_Control = f'{network_path}Control_{tissue}_{cell}_Genes.csv'
        
        Cancer_net  = pd.read_csv(path_Cancer)
        Control_net = pd.read_csv(path_Control)
        
        # Get the counts for the Venn Diagram:
            
        cancer_genes   = set(Cancer_net["0"])
        control_genes  = set(Control_net["0"])
        common_genes   =  len(cancer_genes.intersection(control_genes))
        total_genes    =  len(cancer_genes.union(control_genes))
        unique_control = len(control_genes) - common_genes
        unique_cancer  = len(cancer_genes) - common_genes
        
        # Plot the Venn:
        
        venn2_unweighted(subsets = (unique_control, unique_cancer, common_genes), set_labels = ('', ''), alpha = 0.5
                          ,subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total_genes):1.1%}" + ")",
                          ax=axes[row_count,column_count])
        
        axes[row_count,column_count].set_title(f'{tissue}'.capitalize(), fontweight='bold', fontsize = 23, y = 1) 
        
        # Move in the grid:
        
        if column_count != 1:
            column_count = column_count+1
        else:
            column_count = 0
            row_count = row_count +1
    
    # Set the legend for all the plots:
               
    fig.legend(labels= ["Control","Cancer"],
           borderaxespad=0.1,
           bbox_to_anchor=(0.5, 0.54),
           loc="upper center", frameon=True, ncol = 3)
    
    # Save the plot:
    
    fig.savefig(f'{network_path}Venn_Diagrams_Networks.png', format="png", dpi=600,
                    bbox_inches='tight')
    
def edges(A):
    m = A.shape[0]
    r,c = np.triu_indices(m,1)
    A = A[r,c]
    np.count_nonzero(A == 1.0)
    return(np.count_nonzero(A == 1.0))

def Density(Cancer_nodes,Cancer_edges):
    PC = (Cancer_nodes * (Cancer_nodes - 1))/2
    return(Cancer_edges/PC * 100)
    
def Network_Statistics(Cancer_type_list, Normal_Tissue_list, Cell_type_list):
    
    '''
    This function generate a csv with the network statistics (nodes, edges, density)
        
        Inputs:
                - Cancer_type   : string list, Name of the Cancer.
                - Normal_Tissue : string list, Name of the normal tissue.
                - Cell_type     : string list, Name of the normal cell in the tissue.
    '''
    # Networks path:
    
    network_path = "./Cleaned_version/Data/"
    
    # Open the file to write:
    
    with open(network_path + "Network_Statistics.txt", 'a') as the_file:
        
        
        # Write the titles per each column:
        
        the_file.write("# Name"   + "\t")
        the_file.write("# Nodes"   + "\t")
        the_file.write("# Edges"   + "\t")
        the_file.write("% Density" + "\n")
        
        # Start writting the file:
        
        for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
            
            # Load the Adj matrices:
            
            Cancer_adj  = np.load(f'{network_path}Cancer_{cancer}_PPI.npy', allow_pickle= True)
            Control_adj = np.load(f'{network_path}Control_{tissue}_{cell}_PPI.npy', allow_pickle= True)
            
            # Get info of statistics:
            
            Cancer_nodes   = len(Cancer_adj)
            Cancer_edges   = edges(Cancer_adj)
            Cancer_density = round(Density(Cancer_nodes, Cancer_edges),3)
            
            Control_nodes   = len(Control_adj)
            Control_edges   = edges(Control_adj)
            Control_density = round(Density(Control_nodes, Control_edges),3)
            
            # Write the information:
            
            the_file.write(f'{cancer}' + "\t")
            the_file.write(f'{str(Cancer_nodes)}' + "\t")
            the_file.write(f'{str(Cancer_edges)}' + "\t")
            the_file.write(f'{str(Cancer_density)}' + "\n")
        
            the_file.write(f'{tissue}' + "\t")
            the_file.write(f'{str(Control_nodes)}' + "\t")
            the_file.write(f'{str(Control_edges)}' + "\t")
            the_file.write(f'{str(Control_density)}' + "\n")
            
        # Close the file:
        
        the_file.close()
      
            
def Generate_PPMI(Cancer_type_list, Normal_Tissue_list, Cell_type_list, context_window = 10): 
    
    '''
    This function generates the PPMI matrix from the Adj Matrices and save it as a np:
        
        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - context_window  : int, contex window for the PPMI matrix.
    '''    
    
    # Networks path:
    
    network_path = "./Cleaned_version/Data/"
    counter      = 1
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the Adj: 
        
        Cancer_adj  = np.load(f'{network_path}Cancer_{cancer}_PPI.npy', allow_pickle= True)
        Control_adj = np.load(f'{network_path}Control_{tissue}_{cell}_PPI.npy', allow_pickle= True) 
        
        # Generate the PPMI:
        
        PPMI_Cancer  = generate_embedding_spaces_Human.deep_walk_ppmi(Cancer_adj,  context_window)
        PPMI_Control = generate_embedding_spaces_Human.deep_walk_ppmi(Control_adj,  context_window)
        
        # Save the PPMI:
        
        np.save(f'{network_path}Cancer_PPMI_{cancer}', PPMI_Cancer)
        np.save(f'{network_path}Control_PPMI_{tissue}_{cell}', PPMI_Control)
        
        print(f'{counter}/{len(Cancer_type_list)}')
        
        counter = counter + 1


def Generate_Gene_Embeddings(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, max_iter=1000, verbose = 50):
    
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
    
    network_path = "./Cleaned_version/Data/"
    save_path    = "./Cleaned_version/Data/Embeddings/"
    
    counter = 0
    
    # Generate the embeddings:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        Solver = generate_embedding_spaces.SNMTF(max_iter, verbose) 
        
        # Load the information:
        
        Cancer_adj  = np.load(f'{network_path}Cancer_{cancer}_PPI.npy', allow_pickle= True)
        Control_adj = np.load(f'{network_path}Control_{tissue}_{cell}_PPI.npy', allow_pickle= True) 
        
        Cancer_PPMI  = np.load(f'{network_path}Cancer_PPMI_{cancer}.npy', allow_pickle= True)
        Control_PPMI = np.load(f'{network_path}Control_PPMI_{tissue}_{cell}.npy', allow_pickle= True) 
        
        # Do the embeddings:
        
        for dim in dimension_list:
        
            Solver.Solve_MUR(Cancer_PPMI, dim, save_path, network = f'PPI_{cancer}', matrix = "PPMI_Cancer", init="SVD")
            Solver.Solve_MUR(Cancer_adj,  dim, save_path, network = f'PPI_{cancer}', matrix = "Adj_Cancer",  init="SVD")
            
            Solver.Solve_MUR(Control_PPMI, dim, save_path, network = f'PPI_{tissue}_{cell}', matrix = "PPMI_Control", init="SVD")
            Solver.Solve_MUR(Control_adj,  dim, save_path, network = f'PPI_{tissue}_{cell}', matrix = "Adj_Control",  init="SVD")
            
            counter = counter + 4
            
            print(f'{counter}/{len(Cancer_type_list)*4*(len(dimension_list))}')


def Annotate_Gene_Space(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, Annotation = "GO"):
    
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
    
    network_path   = "./Cleaned_version/Data/"
    embedding_path = "./Cleaned_version/Data/Embeddings/"
    save_path      = "./Cleaned_version/Data/Embeddings/"
    
    counter = 0
    
    # Get the annotations:
    
    if Annotation == "GO":
        
        GO_Matrix = pd.read_csv(f'{network_path}_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv', 
                                index_col = 0, dtype={0: str}) 
        GO_Matrix.index = GO_Matrix.index.astype(str)
        
        print("nop")
        
    elif Annotation == "Leaf":
        
        GO_Matrix = pd.read_csv(f'./Cleaned_version/Data/_Matrix_Genes_GO_BP_PPI.csv', 
                                index_col = 0, dtype={0: str}) 
        GO_Matrix.index = GO_Matrix.index.astype(str)  
        
    else:
        print("Error Check the paths to your annotation matrix")
        
    # Start the embeddings:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
            
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
            
            G_Cancer_PPMI  = np.load(f'{embedding_path}_G_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)
            G_Cancer_Adj   = np.load(f'{embedding_path}_G_Matrix_{dim}_PPI_{cancer}_Adj_Cancer.npy', allow_pickle=True)
            
            G_Control_PPMI = np.load(f'{embedding_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True)
            G_Control_Adj  = np.load(f'{embedding_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_Adj_Control.npy', allow_pickle=True)
            
            # Do the embeddings:
            
            GO_Embeddings.GO_Embeddings_Human(Annotation_Cancer,  G_Cancer_PPMI, 
                                              save_path,  GO = Annotation, network = f'PPI_{cancer}', 
                                              matrix = "PPMI_Cancer",  ortogonal = True)
            
            GO_Embeddings.GO_Embeddings_Human(Annotation_Cancer,  G_Cancer_Adj, 
                                              save_path,  GO = Annotation, network = f'PPI_{cancer}', 
                                              matrix = "Adj_Cancer",  ortogonal = True) 
            
            GO_Embeddings.GO_Embeddings_Human(Annotation_Control,  G_Control_PPMI, 
                                              save_path,  GO = Annotation, network = f'PPI_{tissue}_{cell}', 
                                              matrix = "PPMI_Control",  ortogonal = True)
            
            GO_Embeddings.GO_Embeddings_Human(Annotation_Control,  G_Control_Adj, 
                                              save_path,  GO = Annotation, network = f'PPI_{tissue}_{cell}', 
                                              matrix = "Adj_Control",  ortogonal = True) 
            
            counter = counter + 4
            
            print(f'{counter}/{len(Cancer_type_list)*4*(len(dimension_list))}')


def KMedoids_format(matrix, labels):
    
    # Calculate cosine distance:
    
    distance = pairwise_distances(matrix, metric="cosine")
    
    # Number of clusters (thum-rule):
    
    n_clusters        = round(math.sqrt(len(matrix)/2))
    
    # Compute the kmedoids using a precomputed distance matrix:
    
    clustering   = KMedoids(n_clusters=n_clusters, metric='precomputed',random_state=0).fit(distance)
    
    # Get the clusters:
    
    clusters_lab = clustering.labels_
    
    # Adapt it to the correct format:
    
    clusters = [[] for i in range(n_clusters)]
    n = len(clusters_lab)
    
    for i in range(n):
        clusters[clusters_lab[i]].append(labels[i])
    
    # Return:
    
    return clusters  


def Embedding_Structure(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, Gene_Matrix = False, matrix = "PPMI",
                        annotation = "GO"):
    
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
    
    network_path = "./Cleaned_version/Data/Embeddings/"
    save_cosine  = "./Cleaned_version/Data/FMM/"
    
    control = 0
    
    # Get the structure of the spaces:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        for dim in dimension_list:
    
            if Gene_Matrix == False:
                
                if annotation == "GO":
                
                    path_Cancer  = f'{network_path}_GO_Embeddings_GO_PPI_{cancer}_{dim}_{matrix}_Cancer.csv'
                    path_Control = f'{network_path}_GO_Embeddings_GO_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv'
                    
                elif annotation == "Reactome":
                    
                    path_Cancer  = f'{network_path}_GO_Embeddings_Reactome_PPI_{cancer}_{dim}_{matrix}_Cancer.csv'
                    path_Control = f'{network_path}_GO_Embeddings_Reactome_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv'
                
                elif annotation == "Leaf":
                    
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
                
                if annotation == "GO":
                
                    cosine_Cancer.to_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}.csv', header = True, index=True) 
                    cosine_Control.to_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv', header = True, index=True)
                    
                    control = control + 2
                
                elif annotation == "Reactome":
                    
                    cosine_Cancer.to_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Reactome.csv', header = True, index=True) 
                    cosine_Control.to_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Reactome.csv', header = True, index=True)
                    
                    control = control + 2
                
                elif annotation == "Leaf":
                    
                    cosine_Cancer.to_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Leaf.csv', header = True, index=True) 
                    cosine_Control.to_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Leaf.csv', header = True, index=True)
                    
                    control = control + 2
                
                print(f'{control}/{len(Cancer_type_list) * 2 * len(dimension_list)}')
                
            else:
                path_Cancer  = f'{network_path}_G_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy'
                path_Control = f'{network_path}_G_Matrix_{dim}_PPI_{cancer}_Adj_Cancer.npy'
                
                embeddings_Cancer  = np.load(path_Cancer, allow_pickle=True)
                embeddings_Control = np.load(path_Control, allow_pickle=True)
                
                Cancer_cosine  = pairwise_distances(embeddings_Cancer, metric="cosine")
                Control_cosine = pairwise_distances(embeddings_Control, metric="cosine")
                
                Cancer_cosine.to_csv(f'{save_cosine}Gene_Cosine_Cancer_{cancer}_{dim}_{matrix}.csv', header = True, index=True)
                Control_cosine.to_csv(f'{save_cosine}Gene_Cosine_Control_{tissue}_{cell}_{dim}_{matrix}.csv', header = True, index=True)

def Common_GO_Terms(Cancer_type_list, Normal_Tissue_list, Cell_type_list, annotation = "GO"):

    '''
    This gets the common set of GO terms between cancer/control. Then, the annotated space can be filtered for comparisons.
    
        Inputs:
                - Cancer_type_list     : string list, Name of the Cancer.
                - Normal_Tissue_list   : string list, Name of the normal tissue.
                - Cell_type_list       : string list, Name of the normal cell in the tissue.
    ''' 
    
    path = "./Cleaned_version/Data/"
    save = "./Cleaned_version/Data/"
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load GO Matrix:
        
        if annotation == "GO":
        
            GO_Matrix = pd.read_csv(f'{path}_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv', 
                                    index_col = 0, dtype={0: str}) 
            GO_Matrix.index = GO_Matrix.index.astype(str)
        
        else:
            
            GO_Matrix = pd.read_csv(f'/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/_Matrix_Genes_GO_BP_PPI.csv', 
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
        
        # Intersecions:
        
        Common_Annotations = GO_terms_filtered_Cancer.intersection(GO_terms_filtered_Control)
        
        if annotation == "GO":
        
            pd.DataFrame(Common_Annotations).to_csv(f'{save}Common_Set_{tissue}.csv') 
            
        else:
        
            pd.DataFrame(Common_Annotations).to_csv(f'{save}Common_Set_{tissue}_{annotation}.csv') 
        


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


def Compute_Success(FMMs_common_pairs, terms_2):
    
    success = 0
        
    for pair in list(FMMs_common_pairs.index):
        
        if len(set(pair).intersection(set(terms_2))) > 1:
            success = success+1 
    return(success)
        
def Compute_Success_2(FMMs_common_pairs, terms_2):
    
    x = [item for sublist in list(FMMs_common_pairs.index) for item in sublist]
    return(len(set(x).intersection(set(terms_2))))
    

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
    
    save_cosine  = "./Cleaned_version/Data/FMM/"
    
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
            
            if annotation == "GO":
            
                Cancer_Structure_1  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d1}_{matrix}.csv', index_col = 0)
                Cancer_Structure_2  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d2}_{matrix}.csv', index_col = 0)
                
                Control_Structure_1 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d1}_{matrix}.csv', index_col = 0)
                Control_Structure_2 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d2}_{matrix}.csv', index_col = 0)
            
            else:
                
                Cancer_Structure_1  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d1}_{matrix}_{annotation}.csv', index_col = 0)
                Cancer_Structure_2  = pd.read_csv(f'{save_cosine}Cosine_Cancer_{cancer}_{d2}_{matrix}_{annotation}.csv', index_col = 0)
                
                Control_Structure_1 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d1}_{matrix}_{annotation}.csv', index_col = 0)
                Control_Structure_2 = pd.read_csv(f'{save_cosine}Cosine_Control_{tissue}_{cell}_{d2}_{matrix}_{annotation}.csv', index_col = 0)
            
            # If we calculate it not only with the common set:
            
            if filtered == True:
                if annotation == "GO":
                    
                    print("Nop")
                    
                else:
                    print(annotation)
                    BP = json.load(open("./Cleaned_version/Data/gene2go_Human_PPIGO_Specific_BP.json"))
                
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
                     
            
def Plot_Relative_Error(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                   comparison_list, matrix = "PPMI", annotation = "Leaf"):

    # Path:
    
    save_cosine  = "./Cleaned_version/Data/FMM/"
    row_control = 0
    label_count = 0
    
    # Set the plot characteristics:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 15})
    plt.tight_layout()
    fig, axs = plt.subplots(len(Cancer_type_list),2 , figsize=(15,10)) 
    
    # Reference to the Plots:
    
    plot_labels = []
    
    for i in range(len(Cancer_type_list * 2)):
        
        plot_labels.append(ascii_lowercase[i])
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the relatives errors:
        
        if annotation == "GO":
        
            with open(f'{save_cosine}Relative_Cancer_{cancer}_{matrix}.txt', "r") as fp:
                Cancer_Error = json.load(fp)
            with open(f'{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}.txt', "r") as fp:
                Control_Error = json.load(fp)    
        
        else:
            
            with open(f'{save_cosine}Relative_Cancer_{cancer}_{matrix}_{annotation}.txt', "r") as fp:
                Cancer_Error = json.load(fp)
            with open(f'{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}_{annotation}.txt', "r") as fp:
                Control_Error = json.load(fp)   
            
        # Plot them:
        
        if len(Cancer_type_list) == 1:
        
            axs[0].plot(comparison_list, Cancer_Error, marker='o', label=f'Cancer {tissue} {matrix}', linewidth = 5, markersize=8)
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['left'].set_linewidth(4)
            axs[0].spines['bottom'].set_linewidth(4)
            axs[0].spines['left'].set_color("grey")
            axs[0].spines['bottom'].set_color("grey") 
            axs[0].set_ylabel(' ', fontsize=16, fontweight='bold')
            axs[0].set_xlabel(' ', fontsize=16, fontweight='bold')
            axs[0].set_xticks(fontsize=14, rotation=45)
            axs[0].set_yticks(fontsize=14)
        
            axs[1].plot(comparison_list, Control_Error, marker='o', label='Control {Tissue} {matrix}', linewidth = 5, markersize = 8)   
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)
            axs[1].spines['left'].set_linewidth(4)
            axs[1].spines['bottom'].set_linewidth(4)
            axs[1].spines['left'].set_color("grey")
            axs[1].spines['bottom'].set_color("grey") 
            axs[1].set_ylabel(' ', fontsize=16, fontweight='bold')
            axs[1].set_xlabel(' ', fontsize=16, fontweight='bold')
            axs[1].set_xticks(fontsize=14, rotation=45)
            axs[1].set_yticks(fontsize=14)
        
        else:
            
            axs[row_control][0].plot(comparison_list, Cancer_Error, marker='o', label=f'Cancer {tissue} {matrix}', color="#76de61",linewidth = 5, markersize = 8)
            axs[row_control][0].spines['right'].set_visible(False)
            axs[row_control][0].spines['top'].set_visible(False)
            axs[row_control][0].spines['left'].set_linewidth(4)
            axs[row_control][0].spines['bottom'].set_linewidth(4)
            axs[row_control][0].spines['left'].set_color("grey")
            axs[row_control][0].spines['bottom'].set_color("grey") 
            axs[row_control][0].set_ylabel(' ', fontsize=16, fontweight='bold')
            axs[row_control][0].set_xlabel(' ', fontsize=16, fontweight='bold')
            axs[row_control][0].set_title(plot_labels[label_count].capitalize(), fontweight='bold', fontsize = 17, y =1)
            
            label_count = label_count + 1
            
            if row_control == 3:
                axs[row_control][0].set_xticklabels(comparison_list,fontsize=14)
                axs[row_control][0].get_xticklabels()[3].set_color("red")
            else:
                axs[row_control][0].set_xticklabels(' ')
                
            axs[row_control][1].plot(comparison_list, Control_Error, marker='o', label='Control {Tissue} {matrix}', color="#f1948a", linewidth = 5, markersize = 8)   
            axs[row_control][1].spines['right'].set_visible(False)
            axs[row_control][1].spines['top'].set_visible(False)
            axs[row_control][1].spines['left'].set_linewidth(4)
            axs[row_control][1].spines['bottom'].set_linewidth(4)
            axs[row_control][1].spines['left'].set_color("grey")
            axs[row_control][1].spines['bottom'].set_color("grey") 
            axs[row_control][1].set_ylabel(' ', fontsize=16, fontweight='bold')
            axs[row_control][1].set_xlabel(' ', fontsize=16, fontweight='bold')
            axs[row_control][1].set_title(plot_labels[label_count].capitalize(), fontweight='bold', fontsize = 17, y =1)
            
            label_count = label_count + 1

            
            if row_control == 3:
                axs[row_control][1].set_xticklabels(comparison_list,fontsize=14)
                axs[row_control][1].get_xticklabels()[3].set_color("red")
            else:
                axs[row_control][1].set_xticklabels(' ')
        
        # Control the rows:
            
        row_control = row_control + 1
            
    fig.text(0.5, 0.05, 'Dimensions', ha='center', va='center', 
                     rotation='horizontal', fontsize=20, fontweight = "bold")
    
    fig.text(0.07, 0.5, 'RSE', ha='center', va='center', 
                     rotation='vertical', fontsize=20, fontweight = "bold")
    
    fig.legend(labels= ["Control","Cancer"],
           borderaxespad=0.1,
           bbox_to_anchor=(0.5, 0.01),
           loc="upper center", frameon=True, ncol = 3)
                
    # Save the plot:
    
   # fig.savefig(f'{save_cosine}Relative_Error_{annotation}.png', format="png", dpi=600,
    #                bbox_inches='tight') 
    
    
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
    
    G = graph.from_resource("go-basic") # ./Data/go-basic.obo"
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
            b.append(similarity.lin(G, GO1, GO2)) # Computes the Lim's semantic similarity, other similarities can be also computed changing this function.
        except Exception as PGSSLookupError:
            continue   
    
    
    # Relate the distance with the semantic similarity:
    
    #c = [similarity.values[np.triu_indices(len(similarity), k = 1)]]
    
    c = []

    # Give back the results:
    
    if len(b) > 0:
        return(np.mean(a), np.std(a), np.mean(b), np.std(b), np.mean(c), np.std(c), top_100_similar_values,top_100_dissimilar_values)
    else:
        print("error")
        

def Parallel_Functional_Organization(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                    matrix = "PPMI", filtered = True, Common = False, Jaccard = False, number_similar = 500, annotation = "GO"):
    
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
    
    save_cosine     = "./Cleaned_version/Data/FMM/"
    functional_path = "./Cleaned_version/Data/"
    Filter_path     = "./Cleaned_version/Data/"
    
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
        
def Compute_Jaccard_Top_Pairs(Cancer, Control):
       
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
    
def Calculate_limits(list_mean, list_std):
    
    lower_error = []
    upper_error = []
    
    for i in range(len(list_mean)):
        
        lower = list_mean[i] - list_std[i]
        upper = list_mean[i] + list_std[i]
        
        if lower < 0:
            lower = list_mean[i]
        else:
            lower = list_std[i] 
            
        if upper > 1:
            upper = 1 - list_mean[i]
        else:
            upper = list_std[i] 
            
        
        lower_error.append(lower)
        upper_error.append(upper)
    
    result = [lower_error, upper_error]
    
    return(result)
                 
def Plot_Functional_Organization(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list,
                   matrix = "PPMI", Common = False, number_similar = 500,annotation = "GO", note = "Common_Set"):
    
    # Paths:
    
    path_func = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"
    
    if Common == False:
        label = "Individual"
    else:
        label = "Common"
    
    # Controlers:
    
    row_control = 0
    label_count = 0
    
    # Authomatic Labels:
    
    plot_labels = []
    
    for i in range(len(Cancer_type_list * 3)):
        
        plot_labels.append(ascii_lowercase[i])
    
    # Prepare the grid:
        
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    
    fig, axes = plt.subplots(len(Cancer_type_list), 3, figsize=(15,8), gridspec_kw={"width_ratios": [1,1,0.8]})
    
    # Start plotting:

    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Open the file:
        
        if annotation == "GO":
            Organization = open(f'{path_func}Similarity_{tissue}_{note}_Set_{number_similar}.json' ,)
        else:
             Organization = open(f'{path_func}Similarity_{tissue}_{note}_Set_{number_similar}_{annotation}.json' ,)
            
        Organization = json.load(Organization)
        
        # Line Plot Control Similar:
        
        axes[row_control,0].errorbar(dimension_list, Organization["Sim_Control"],  yerr= Calculate_limits(Organization["Sim_Control"], Organization["Sim_Control_std"]),    
            marker = "o", capsize=5, lw = 3, color = "#76de61")
            
        axes[row_control,0].errorbar(dimension_list, Organization["Diss_Control"],  yerr= Calculate_limits(Organization["Diss_Control"], Organization["Diss_Control_std"]),    
            marker = "p",fmt ='--', capsize=5, lw = 3, color = "#76de61")
            
        axes[row_control,0].set_xticks(dimension_list)
        axes[row_control,0].spines['right'].set_visible(False)
        axes[row_control,0].spines['top'].set_visible(False)
        axes[row_control,0].spines['left'].set_linewidth(4)
        axes[row_control,0].spines['bottom'].set_linewidth(4)
        axes[row_control,0].spines['left'].set_color("grey")
        axes[row_control,0].spines['bottom'].set_color("grey")  
        axes[row_control,0].spines['bottom']
        axes[row_control,0].spines['bottom']
        axes[row_control,0].set_ylabel(' ')
        axes[row_control,0].set_xlabel(' ')
        axes[row_control,0].set_title(f'{plot_labels[label_count].capitalize()}', fontweight='bold', fontsize = 17, y = 0.9)
        label_count = label_count + 1
        
        if row_control != len(Cancer_type_list) -1:
            axes[row_control,0].xaxis.set_major_formatter(NullFormatter())
        
        # Line plot with dissimilar:
        
        axes[row_control,1].errorbar(dimension_list, Organization["Sim_Cancer"],  yerr= Calculate_limits(Organization["Sim_Cancer"], Organization["Sim_Cancer_std"]),    
            marker = "o", capsize=5, lw = 3, color = "#f1948a")
        axes[row_control,1].errorbar(dimension_list, Organization["Diss_Cancer"],  yerr= Calculate_limits(Organization["Diss_Cancer"], Organization["Diss_Cancer_std"]),    
            marker = "p", fmt ='--', capsize=5, lw = 3, color = "#f1948a")
            
        axes[row_control,1].set_xticks(dimension_list)
        axes[row_control,1].spines['right'].set_visible(False)
        axes[row_control,1].spines['top'].set_visible(False)
        axes[row_control,1].spines['left'].set_linewidth(4)
        axes[row_control,1].spines['bottom'].set_linewidth(4)
        axes[row_control,1].spines['left'].set_color("grey")
        axes[row_control,1].spines['bottom'].set_color("grey") 
        axes[row_control,1].spines['bottom']
        axes[row_control,1].spines['bottom']
        axes[row_control,1].set_ylabel(' ')
        axes[row_control,1].set_xlabel(' ')
        axes[row_control,1].yaxis.set_major_formatter(NullFormatter())
        axes[row_control,1].set_title(f'{plot_labels[label_count].capitalize()}', fontweight='bold', fontsize = 17, y = 0.9)
        label_count = label_count + 1
        
        if row_control != len(Cancer_type_list) -1:
            axes[row_control,1].xaxis.set_major_formatter(NullFormatter())
            
        # Bar plot for the tope500 Jaccard:
        
        ax = sns.barplot(dimension_list, Organization["Jaccard_Similar"], ax=axes[row_control,2], color = "#55a9f3")  
        ax.set_yticks([0,0.95])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(4)
        ax.spines['bottom'].set_linewidth(4)
        ax.spines['left'].set_color("grey")
        ax.spines['bottom'].set_color("grey")
        ax.yaxis.set_major_formatter(NullFormatter())
        if row_control != len(Cancer_type_list) -1:
            ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_title(f'{plot_labels[label_count].capitalize()}', fontweight='bold', fontsize = 17, y = 0.9)
            
        label_count = label_count + 1
        
        # Next row:
        
        row_control = row_control+1
        
    fig.text(0.55, 0.05, 'Dimensions', ha='center', va='center', 
                     rotation='horizontal', fontsize=20, fontweight = "bold")
    
    fig.text(0.09, 0.5, 'Semantic Similarity', ha='center', va='center', 
                     rotation='vertical', fontsize=20, fontweight = "bold")
    
    fig.text(0.68, 0.5, 'Jaccard Index', ha='center', va='center', 
                     rotation='vertical', fontsize=20, fontweight = "bold")
    
    # Generate a manual legend:
    
    Control_patch    = mpatches.Patch(color='#76de61', label='Control')
    Cancer_patch     = mpatches.Patch(color='#f1948a', label='Cancer')
                                 
    plt.legend(handles=[Control_patch, Cancer_patch,lines.Line2D([0], [0], marker='', ls='-', c='black', label = "Similar"),
                        lines.Line2D([0], [0], marker='', ls='--', c='black', label = "Diimilar")],borderaxespad=0.1,
           bbox_to_anchor=(-0.80, -0.72), loc="upper center", frameon=True, ncol = 2)
    
    if annotation == "GO":
        fig.savefig(f'{path_func}Functional_Organization_{label}_{note}.png', format="png", dpi=600,
                    bbox_inches='tight') 
    else:
        fig.savefig(f'{path_func}Functional_Organization.png_{annotation}_{label}_{note}', format="png", dpi=600,
                    bbox_inches='tight') 
        
        
        
def Movement_Ranking(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI", annotation = "Leaf"):
    
    # Paths:
    
    cosine_path   = "./Cleaned_version/Data/FMM/" 
    Filter_path   = "./Cleaned_version/Data/"   
    movement_path = "./Cleaned_version/Data/"  
    
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
        
        
def Global_moving_in_the_space():
    
    '''
    This function test if the cancer-related annotations are moving statistically more than the rest.
    '''
    
    # Load our set of cancer-related annotations:
    
    with open(f'./Cleaned_version/Data/enriched_in_cancer_gobp_terms_cosmic.txt', 'r') as fp:
        cancer_related_go_terms =json.load(fp)
    
    # Do the test
        
    cancers = ["breast","prostate","lung","colon"]
    
    for i in cancers:

        Ranking = pd.read_csv(f'./Cleaned_version/Data/Rank_movement_{i}_PPMI_Leaf.csv',index_col=0,names=['norm'],header=0)

        # Globally validate the movement is connected with cancer:

        a = Ranking.loc[~Ranking.index.isin(cancer_related_go_terms)]['norm']
        b = Ranking.loc[Ranking.index.isin(cancer_related_go_terms)]['norm']
        print(i,mannwhitneyu(a,b,alternative='less'))


def moving_in_the_space():
    
    '''
    This function performs the enrichment analyses of cancer-related annotations.
    '''
    
    # Load our set of cancer-related annotations:
    
    with open(f'./Cleaned_version/Data/enriched_in_cancer_gobp_terms_cosmic.txt', 'r') as fp:
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
        Ranking = pd.read_csv(f'/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/Rank_movement_{i}_PPMI_Leaf.csv',index_col=0,names=['norm'],header=0)

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

def Reorder(Cancer_structure, Control_structure):
    
    # Do a clustering on the Control to see how it affects the Cancer:
    
    linkage       = hc.linkage(sp.distance.squareform(Control_structure),method='average')
    Control_clust = sns.clustermap(Control_structure,cmap="jet", row_linkage=linkage, col_linkage=linkage)
    order         = Control_clust.dendrogram_row.reordered_ind
        
    order_cancer  = Cancer_structure.index[order]
    cancer_final  = Cancer_structure.reindex(order_cancer)
    cancer_final  = cancer_final[cancer_final.index]
        
    order_control  = Control_structure.index[order]
    control_final  = Control_structure.reindex(order_control)
    control_final  = control_final[control_final.index]
        
    return(control_final, cancer_final)
    
    
def Query_GO_terms(Normal_tissue_list):
    
    '''
    To run this function is mandatory to manually obtain the top 100 moving file is needed. This file can be recomputed manually
    for new runs. The file contains the description of the GO term and its total movement. Since this ranking can variate depending
    on the annotation file that is used (and also on the embeddings) it should be recomputed with the new ranking if needed.
    '''
    
    movement_path = "./Cleaned_version/Data/"
    path_rank     = "./Cleaned_version/Data/"
    
    for tissue in Normal_tissue_list:
        
        file      = pd.read_csv(f'{movement_path}top_100_moving_{tissue}.csv')
        
        rank      = pd.read_csv(f'{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv')
        moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"]
        
   
        file = file.head(len(moving_GO))
        
        cancer_q = f'{tissue} cancer'
        
        counts = []
        
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
    
    
def Compare_Organizations(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI", annotation = "GO"):
    
    # Paths:
    
    cosine_path = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    fucntional  = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"
    Filter_path = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    
    # Set the grid:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    
    fig, axes = plt.subplots(len(Cancer_type_list), 4, figsize=(15,15), gridspec_kw={"width_ratios": [1,1,1,0.05]})
    
    fig.tight_layout()
    
    # Controllers
    
    plot_labels = []
    
    for i in range(len(Cancer_type_list * 3)):
        
        plot_labels.append(ascii_lowercase[i])
        
    control_lables = 0
    control_rows   = 0
    control_it     = 0
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the cosine-cosine matrix:
        
        if annotation == "GO":
        
            Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}.csv',         index_col = 0)
            Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}.csv', index_col = 0)
            
            # Filter the matrices with the common set:
            
            common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}.csv')["0"])
            
            Cancer_structure  = Cancer_structure.loc[common_set, common_set]
            Control_structure = Control_structure.loc[common_set, common_set]
            
        
        elif annotation == "Reactome":
            
            Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_Reactome.csv',         index_col = 0)
            Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_Reactome.csv', index_col = 0)
            
            # Filter the matrices with the common set:
            
            common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Reactome.csv')["0"])
            
            Cancer_structure  = Cancer_structure.loc[common_set, common_set]
            Control_structure = Control_structure.loc[common_set, common_set]
            
        elif annotation == "Leaf":
            
            Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_Leaf.csv',         index_col = 0)
            Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_Leaf.csv', index_col = 0)
            
            # Filter the matrices with the common set:
            
            common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Leaf.csv')["0"])
            
            Cancer_structure  = Cancer_structure.loc[common_set, common_set]
            Control_structure = Control_structure.loc[common_set, common_set]            

        # Reorder the cosine matrices:
        
        control_final, cancer_final = Reorder(Cancer_structure, Control_structure)
        
        diference = Control_structure - Cancer_structure 
        
        # Plot:
        
        im = axes[control_rows,0].matshow(control_final.values, interpolation='nearest', cmap="jet", vmin=0, vmax=1)
        axes[control_rows,0].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
        axes[control_rows,0].axis('off')
        control_lables = control_lables + 1
        
        im = axes[control_rows,1].matshow(cancer_final.values, interpolation='nearest', cmap="jet", vmin=0, vmax=1)
        axes[control_rows,1].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
        axes[control_rows,1].axis('off')
        control_lables = control_lables + 1 
        
        im2 = axes[control_rows,2].matshow(diference.values, interpolation='nearest', cmap="seismic", vmin=-1, vmax=1)
        axes[control_rows,2].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
        axes[control_rows,2].axis('off')
        control_lables = control_lables + 1 
        control_rows = control_rows + 1
        
        control_it = control_it + 1
        
        print(f'{control_it}/{len(Cancer_type_list)}')
        
    # Delete the axes that are not needed:
    
    # Color bars:

    fig.colorbar(im, cax= axes[1,3], 
                         ax=[axes[0,0],axes[0,1],axes[1,0],
                             axes[1,1],axes[2,0],axes[2,1]])

    fig.colorbar(im2, cax= axes[2,3], 
                         ax=[axes[0,2], axes[1,2], axes[2,2]])
    
    # Save the plot:

    fig.savefig(f'{fucntional}Structure_Comparisons_{matrix}_{annotation}.png', format="png", dpi=600,
                    bbox_inches='tight')
    

    
    
def Compare_Organizations_UNI(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI", annotation = "GO"):
    
    # Paths:
    
    cosine_path = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    fucntional  = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"
    Filter_path = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
            
    cancer = "lung cancer"
    tissue = "lung"
    cell   = 'alveolar cells'

    Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_Leaf.csv',         index_col = 0)
    Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_Leaf.csv', index_col = 0)
            
    # Filter the matrices with the common set:
            
    common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Leaf.csv')["0"])
            
    Cancer_structure  = Cancer_structure.loc[common_set, common_set]
    Control_structure = Control_structure.loc[common_set, common_set]            

    # Reorder the cosine matrices:
        
    control_final, cancer_final = Reorder(Cancer_structure, Control_structure)
    
    diference = control_final - cancer_final 

    # Set the grid:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    
    fig, axes = plt.subplots(1, 5, figsize=(15,5), gridspec_kw={"width_ratios": [1,1,1,0.05, 0.05]})    
    fig.tight_layout()

    control_lables = 0
    # Controllers

        
    plot_labels = []
    
    for i in range(len(Cancer_type_list * 3)):
        
        plot_labels.append(ascii_lowercase[i])
        
    # Plot:
        
    im = axes[0].matshow(control_final.values, interpolation='nearest', cmap="jet", vmin=0, vmax=1)
    axes[0].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
    axes[0].axis('off')
    control_lables = control_lables + 1
        
    axes[1].matshow(cancer_final.values, interpolation='nearest', cmap="jet", vmin=0, vmax=1)
    axes[1].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
    axes[1].axis('off')
    control_lables = control_lables + 1 
        
    im2 = axes[2].matshow(diference.values, interpolation='nearest', cmap="seismic", vmin=-1, vmax=1)
    axes[2].set_title(f'{plot_labels[control_lables].capitalize()}', fontweight='bold', fontsize = 17, y =1)
    axes[2].axis('off')
    control_lables = control_lables + 1 
        
        
    # Delete the axes that are not needed:
    
    # Color bars:
    

    fig.colorbar(im, cax= axes[3], ax=[axes[0], axes[1]])

    fig.colorbar(im2, cax= axes[4], ax=[axes[3]])
 
    
    
    # Save the plot:

    fig.savefig(f'{fucntional}UNI_Structure_Comparisons_{matrix}_{annotation}.png', format="png", dpi=600,
                    bbox_inches='tight')

def Movement_Ranking_Genes(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI"):

    # Paths:
    
    movement_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/" 
    network_path  = "/media/sergio/sershiosdisk/Cancer/Networks/"
    
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        
        # Load the gene names:
        
        genes_cancer  = list(pd.read_csv(f'{network_path}Cancer_{cancer}_Genes.csv')["0"])
        genes_control = list(pd.read_csv(f'{network_path}Control_{tissue}_{cell}_Genes.csv')["0"])
        
        common_genes = list(set(genes_cancer).intersection(genes_control))
        
        # Load the information of the embeddings:
            
        Cancer_genes   = np.load(f'{network_path}_P_Matrix_{optimal_dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)
        Cancer_midle   = np.load(f'{network_path}_U_Matrix_{optimal_dim}_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)
                
        Control_genes  = np.load(f'{network_path}_P_Matrix_{optimal_dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True)
        Control_midle  = np.load(f'{network_path}_U_Matrix_{optimal_dim}_PPI_{tissue}_{cell}_PPMI_Control.npy', allow_pickle=True)   
        
        # Get the coordinates:
        
        cancer_coordinates   = Cancer_genes.dot(Cancer_midle)
        control_coordinates  = Control_genes.dot(Control_midle)
        
        # Do the FMM:
    
        cosine_Cancer  = pd.DataFrame(pairwise_distances(cancer_coordinates, metric="cosine"))
        cosine_Control = pd.DataFrame(pairwise_distances(control_coordinates, metric="cosine"))
        
        cosine_Cancer.columns = genes_cancer
        cosine_Cancer.index   = genes_cancer
        
        cosine_Control.columns = genes_control
        cosine_Control.index   = genes_control
        
        # Common genes:
        
        cosine_Cancer  = cosine_Cancer.loc[common_genes,common_genes]
        cosine_Control = cosine_Control.loc[common_genes,common_genes]
        
        # Calculate the movement:
        
        movement = cosine_Control - cosine_Cancer
        movement = movement.apply(np.linalg.norm, axis=1) 
        movement = movement.sort_index(ascending = False)
        
        # Save the ranking:
        
        movement.to_csv(f'{movement_path}Rank_movement_{tissue}_{matrix}_Leaf_Genes.csv') 
        
def Check_Enrichments_Genes_Movement(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI"):
    
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/" 
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        ranking               = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_{matrix}_Leaf_Genes.csv', dtype={0: str}) 
        ranking = ranking.sort_values(by = "0", ascending = False)
               
        mov_threshold = np.percentile(ranking["0"], 25)
        stb_threshold = np.percentile(ranking["0"], 75)
        
        # moving not movint (is correct change the names):
        
        moving_genes = ranking[ranking["0"] >= stb_threshold]
        stb_genes    = ranking[ranking["0"] <= mov_threshold]
        
        # Enrichment:
        
        Cancer_genes = pd.read_csv("/media/sergio/sershiosdisk/Cancer/Cancer_GO/cancer_gene_census_COSMIC.csv",sep = ",", encoding='latin-1')["Entrez GeneId"]
        Cancer_genes = Cancer_genes.dropna()
        Cancer_genes = [str(int(i)) for i in Cancer_genes]
        
        Total_success_mov = sum(moving_genes["Unnamed: 0"].isin(Cancer_genes))
        Total_move        = len(moving_genes)
        
        Total_stb_succ   = sum(stb_genes["Unnamed: 0"].isin(Cancer_genes))
        Total_stb        = len(stb_genes)
        
        Total      = len(ranking)
        Total_succ = sum(ranking["Unnamed: 0"].isin(Cancer_genes))
        
     
        fold_stb, p_value_stb = Fold_enriched(Total, Total_succ, Total_move, Total_success_mov) 
        fold_mv, p_value_mv   = Fold_enriched(Total, Total_succ, Total_stb, Total_stb_succ) 
        
        print(fold_stb, p_value_stb)
        print(fold_mv, p_value_mv)


        
        
        
def Overlapping_Moving_Genes(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI"):
    
    # Paths:
    
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/" 
    prediction_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    
    ensmbl_transformation = pd.read_csv("/media/sergio/sershiosdisk/Cancer/mart_export.txt", sep='\t', dtype={1: str})
    ensmbl_transformation = ensmbl_transformation.dropna()
    ensmbl_transformation = ensmbl_transformation.drop_duplicates()
    
    with open(f'{prediction_path}Movement_Statistics_Leaf_PPMI_Genes_Overlapping.txt', 'a') as the_file:
        
        # Write the titles per each column:
                
        the_file.write("# Cancer"       + "\t")
        the_file.write("# Overlaping"   + "\n")
   
        for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
            
            # Load the Ranking:
            
            ranking      = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_{matrix}_Leaf_Genes.csv', dtype={0: str})        
            moving_genes = ranking[ranking["0"] > (np.mean(ranking["0"]) + 2*np.std(ranking["0"]))]
            
            # Transform to genes:
            
            names        = moving_genes["Unnamed: 0"]
            names_common = list(set(names).intersection(ensmbl_transformation['NCBI gene (formerly Entrezgene) ID']))
            
            ensmbl_transformation = ensmbl_transformation[ensmbl_transformation['NCBI gene (formerly Entrezgene) ID'].isin(names_common)]
            moving_genes          = moving_genes[moving_genes["Unnamed: 0"].isin(names_common)]
            
            ensmbl_transformation = ensmbl_transformation.sort_values( by = "NCBI gene (formerly Entrezgene) ID")
            moving_genes          = moving_genes.sort_values( by = "Unnamed: 0")
            
            moving_genes["Names"] = list(ensmbl_transformation["Gene name"])
            
            # Load the predicted genes based on the GO terms:
            
            predicted_GO = pd.read_csv(f'{prediction_path}Predictions_Rank_{tissue}.csv')
            predicted_GO = list(predicted_GO["name"])
            
            # Evaluate the intersection:
            
            x = len(set(moving_genes["Names"]).intersection(set(predicted_GO)))
            x = round(x/len(moving_genes["Names"])*100,2)
            
            the_file.write(f'{tissue}\t')
            the_file.write(f'{x}\n')
            
    the_file.close()      
        
        
        
def Hallmarks_Cancer_GO():
    
    original = pd.read_excel("/media/sergio/sershiosdisk/Cancer/Cancer_GO/GO_hallmarks.xlsx", sheet_name=None)
    
    hallmark = 'Sustaining Proliferative Signal'
    
    Result = pd.DataFrame()
    
    for hallmark in original.keys():
        
        GO_terms    = original[hallmark]["GO Terms"]
        
        try:
            Description = original[hallmark]["Name"]
        except KeyError:
            Description = original[hallmark]["Names"]
            
        Hallmark    = hallmark 
        
        Result_it = pd.DataFrame({"GO" : GO_terms, "Description" : Description, "Hallmark" : Hallmark})
        Result = pd.concat([Result, Result_it])
        
    Result.to_csv(f'/media/sergio/sershiosdisk/Cancer/Cancer_GO/Hallmarks_Cancer_GO_list.csv') 


def Calculate_Sim_Rank_Hallmarks(Normal_Tissue_list, matrix = "PPMI", annotation = "Leaf"):
    
    # Paths:
    
    movement_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"  

    # Load hallmarks:
    
    hallmarks = pd.read_csv('/media/sergio/sershiosdisk/Cancer/Cancer_GO/Hallmarks_Cancer_GO_list.csv', index_col = 0)
    
    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)
    
    for tissue in Normal_Tissue_list:
        
        # Load ranking:
        
        rank = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_{matrix}_{annotation}.csv', index_col = 0)
        
        # Get the semantic similarity between the rank and the hallmarks:
    
        result  = []
        dec_fin = []
        
        for GO1 in list(rank.index):
            result_it = []
            desc      = []
            for GO2 in list(hallmarks.GO):
                try:
                    semantic = similarity.lin(G, GO1, GO2)
                    result_it.append(semantic)
                    desc.append(hallmarks[hallmarks.GO == GO2].Hallmark[hallmarks[hallmarks.GO == GO2].Hallmark.index[0]])
                except Exception as PGSSLookupError:
                    result_it.append(0)
                    desc.append("NO")
                    continue
            result.append(max(result_it))
            dec_fin.append(desc[result_it.index(max(result_it))])
            
        # Get the tail and the top:
        
        rank["Sim"] =  result 
        rank["Desc"] =  dec_fin  
        
        # Solve the problem with some definitions:
        
        hallmark_list = []
        
        for i in list(rank.index):
            if type(rank.loc[i].Desc) != str:
                hallmark_list.append(list(rank.loc[i].Desc)[0])
            else:
                    hallmark_list.append(rank.loc[i].Desc)
                    
        rank["Desc"] =  hallmark_list 
        
        # Save:
        
        rank.to_csv(f'{movement_path}Rank_Hallmarks_{tissue}.csv')
            

def Plot_Movement_Hallmarks(Normal_Tissue_list, matrix = "PPMI", annotation = "Leaf", number = 500):
    
    # Paths:
    
    movement_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # Plot features:
    
    # Set the grid:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 20})
    plt.rc('font', weight='bold')
    
    fig, axes = plt.subplots(2, 2, figsize=(10,8))
    fig.tight_layout()
    
    # Label control:
    
    plot_labels = []
    
    for i in range(len(Normal_Tissue_list)):
        
        plot_labels.append(ascii_lowercase[i])
        
    # Iteration control:
    
    column_it = 0
    row_it    = 0
    label_it  = 0
     
    for tissue in Normal_Tissue_list:
        
        # Load the ranking with the semantic similarity:
        
        rank= pd.read_csv(f'{movement_path}Rank_Hallmarks_{tissue}.csv', index_col = 0)
        
        top500  = rank["Sim"].head(number)
        tail500 = rank["Sim"].tail(number) 
        
        axes[row_it, column_it].hist(tail500, color = "#37bd70", alpha = 0.7)
        axes[row_it, column_it].hist(top500, color = "#bd37bb", alpha = 0.7)
            
        axes[row_it, column_it].set_xlabel("")
        axes[row_it, column_it].set_ylabel("")
        
        axes[row_it, column_it].spines['right'].set_visible(False)
        axes[row_it, column_it].spines['top'].set_visible(False)
        axes[row_it, column_it].spines['left'].set_linewidth(4)
        axes[row_it, column_it].spines['bottom'].set_linewidth(4)
        axes[row_it, column_it].spines['left'].set_color("grey")
        axes[row_it, column_it].spines['bottom'].set_color("grey")
        axes[row_it,column_it].set_title(f'{plot_labels[label_it].capitalize()}', fontweight='bold', 
            fontsize = 17, y = 1, x = 0)
        
        label_it  +=1
        column_it +=1
        
        if column_it > 1:
            column_it = 0
            row_it +=1
            
        fig.text(0., 0.5, 'Counts', ha='center', va='center',
                  rotation='vertical', fontsize=20, fontweight = "bold")
        
        fig.text(0.5, 0, 'Semantic Similarity', ha='center', va='center',
                  rotation='horizontal', fontsize=20, fontweight = "bold")
        
        fig.legend(labels= ["Stable","Moving"],
        borderaxespad=0.1,
        bbox_to_anchor=(0.5, -0.02),
        loc="upper center", frameon=True, ncol = 3)
        
    fig.savefig(f'{movement_path}Semantic_Dist_Hallamrks_{matrix}_{annotation}_{tissue}.png', format="png", dpi=600,
                    bbox_inches='tight')

def Hallmarks_Movement(Normal_Tissue_list, similarity_thr = 0.5):
    
    # Paths:
    
    movement_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    for tissue in Normal_Tissue_list:
        
        # Load the rank with the semantic information:
        
        rank= pd.read_csv(f'{movement_path}Rank_Hallmarks_{tissue}.csv', index_col = 0)
        
        # Create the hallamark clusters:
        
        hallmark_list = []
        
        for i in list(rank.index):
            if rank.loc[i].Sim <= similarity_thr:
                hallmark_list.append("No Classy")
            else:
                hallmark_list.append(rank.loc[i].Desc)
        
        rank["Desc"] = hallmark_list
        
        # Start the plots:
        


def Statistics_Movement_Cancer(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim, matrix = "PPMI",
                               annotation = "Leaf", number = 500):
    
    '''
    This function produce a preliminary study of the relation between the movemt (norm) of the GO terms and their link to cancer
    related functions.
    '''
    
    # Paths:
    
    path         = "/media/sergio/sershiosdisk/Cancer/Cancer_GO/"
    path_cos     = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    Filter_path  = "/media/sergio/sershiosdisk/Cancer/Common_Set/" 
    save_path    = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # Open the cancer reference:
    
    if (annotation == "GO") | (annotation == "Leaf"):
    
        with open(f'{path}cancer_related_go_terms (1).txt', 'r') as f:
            cancer_related_go_terms = set()
            for l in f:
                l = l.rstrip()
                cancer_related_go_terms.add(l)
                
    elif annotation == "Reactome":
        print("add the reference")
            
    # Test the different samples:
    
    with open(f'{save_path}Movement_Statistics_{annotation}_{matrix}_{number}.txt', 'a') as the_file:
        
        # Write the titles per each column:
                
        the_file.write("# Cancer"   + "\t")
        the_file.write("# Common Set"   + "\t")
        the_file.write("# MannWhit Less"   + "\t")
        the_file.write("# Top500 Cancer"   + "\t")
        the_file.write("# Tail500 Cancer"   + "\n")
    
        for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
            # Load the cosine matrices:
            
            if annotation!= "GO":
                cos_control = pd.read_csv(f'{path_cos}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_{annotation}.csv',index_col=0) 
                cos_cancer  = pd.read_csv(f'{path_cos}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_{annotation}.csv',index_col=0)
                common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_{annotation}.csv')["0"])
            else:
                cos_control = pd.read_csv(f'{path_cos}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}.csv',index_col=0) 
                cos_cancer  = pd.read_csv(f'{path_cos}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}.csv',index_col=0)
                common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}.csv')["0"])                

                # Subset the matrices:
            
                cos_cancer  = cos_cancer.loc[common_set, common_set]
                cos_control = cos_control.loc[common_set, common_set] 
                
                difference = cos_control - cos_cancer
                
                
                # Plottini:
                
                diagonal = difference.mask(np.triu(np.ones(difference.shape)).astype(bool)).stack().sort_values(ascending=True)
                
                # Ranking:
                
                ranking = pd.DataFrame(difference.apply(np.linalg.norm,axis=1),columns=['norm'])
                
                non_cancer     = ranking.loc[~ranking.index.isin(cancer_related_go_terms)]['norm']
                cancer_rel     = ranking.loc[ranking.index.isin(cancer_related_go_terms)]['norm']
                
                # MannWhit test:
                
                p_value = mannwhitneyu(non_cancer,cancer_rel,alternative='less')[1]
                
                # Top/Tail 500 proportions:
                
                top500_moving_norm  = set(ranking.sort_values(by='norm',ascending=False)[0:number].index)
                less500_moving_norm = set(ranking.sort_values(by='norm',ascending=False)[-number:].index)
                
                # Fit the information:
                
                the_file.write(f'{tissue}' + "\t")
                the_file.write(f'{len(common_set)}' + "\t")
                the_file.write(f'{p_value}' + "\t")
                the_file.write(f'{len(top500_moving_norm.intersection(cancer_related_go_terms))/number*100}' + "\t")
                the_file.write(f'{len(less500_moving_norm.intersection(cancer_related_go_terms))/number*100}' + "\n")
            
        the_file.close() 
 

def Plot_movement_PanCancers(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim,
                            matrix = "PPMI", annotation = "Leaf", number = 5):
    
    # Paths:
    
    cosine_path     = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    common_path     = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # Take the pan cancer GO terms:
    
    # Do the dictionary:
    
    ranking_dictionary = {}
    for tissue in Normal_Tissue_list:
        if annotation != "GO":
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI_{annotation}.csv')
        else:
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI.csv')
        
        top = set(Ranking["Unnamed: 0"][:500])
        ranking_dictionary[tissue] = top
    
    # Common PanCancer:
    
    common_pan_cancer = list(ranking_dictionary['breast'].intersection(ranking_dictionary['prostate'],
                                          ranking_dictionary['lung'],ranking_dictionary['colon']))
    
    # Plot characteristics:
    
    row_control    = 0
    column_control = 0
    label_control  = 0
    
    # Prepare the labels for each sub-plot:
    
    plot_labels = []
    
    for i in range((number*4)):
        plot_labels.append(ascii_lowercase[i])
         
    # Start the plot:
            
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 20})
    plt.rc('font', weight='bold')
            
    fig, axes = plt.subplots(number, 4, figsize=(15,15))
    fig.tight_layout()
    
    # Start the plotting
    
    # Itarate by row (GO terms):
    
    for GO1 in common_pan_cancer[0:number]:
        
        # And by cancer type:
        
        for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
            
        # Load the cosine matrices:
        
            Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_{annotation}.csv',         index_col = 0)
            Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_{annotation}.csv', index_col = 0)
            
            # read the common set:
            
            common = list(pd.read_csv(f'{common_path}Common_Set_{tissue}_Leaf.csv')["0"])
            
            Cancer_structure   = Cancer_structure.loc[common][common]
            Control_structure  = Control_structure.loc[common][common]
            
            Cancer_structure_GO1  = Cancer_structure.loc[GO1]
            Control_structure_GO1 = Control_structure.loc[GO1]
            
            # Get the movement of the corresponding pan-cancer GO term
            
            moving_GO1 = Control_structure_GO1 - Cancer_structure_GO1
            
            control = []
            cancer  = []

            for i in n.index:
                row    = i[0]
                column = i[1]
                
                control.append(Control_structure.loc[row, column])
                cancer.append(Cancer_structure.loc[row, column])
                          
            
            # Get the information of this GO term to plot:
            
            colors = []

            n_dots = len(moving_GO1)   
            angs = np.linspace(0, 2*np.pi, n_dots)  
            cx, cy = (0, 0)  
            xs, ys = [], []
            ra = 1   
            
            count = 0
            for ang in angs:
                # radius:
                ra_it = ra - moving_GO1[count]
                
                if ra_it > (1 + 2*np.std(moving_GO1)):
                    colors.append("#2679ea")
                elif ra_it < (1 - 2*np.std(moving_GO1)):
                    colors.append("#bd3751")
                else:
                    colors.append("#50dc1f")
                x = cx + ra_it*np.cos(ang)
                y = cy + ra_it*np.sin(ang)
                xs.append(x)   
                ys.append(y)   
                count+=1
            
            # plot the corresponding cell in the multi-plot:
            
            axes[row_control, column_control].scatter(xs, ys, s=9 ,c=colors)  # plot points 
            circle = plt.Circle([0,0], 0.05, color = "#04fbe1")
            axes[row_control, column_control].add_patch(circle)
            axes[row_control, column_control].spines['right'].set_visible(False)
            axes[row_control, column_control].spines['top'].set_visible(False)
            axes[row_control, column_control].spines['left'].set_linewidth(4)
            axes[row_control, column_control].spines['bottom'].set_linewidth(4) 
            axes[row_control, column_control].spines['left'].set_color("grey")
            axes[row_control, column_control].spines['bottom'].set_color("grey")
            axes[row_control, column_control].set_yticklabels(" ")
            axes[row_control, column_control].set_xticklabels(" ") 
            axes[row_control, column_control].set_xticks([])
            axes[row_control, column_control].set_yticks([])
            axes[row_control, column_control].set_title(f'{plot_labels[label_control].capitalize()}', fontweight='bold', 
            fontsize = 17, y = 1, x = 0.5)
            
            # Add the name of the row:
            
            if column_control == 0:
                axes[row_control, column_control].set_ylabel(f'{GO1}', 
                    fontsize=17, fontweight = "bold")
            
            # Add the name of the column:
            
            if row_control == (number-1):
                axes[row_control, column_control].set_xlabel(f'{tissue} Cancer', 
                    fontsize=17, fontweight = "bold")
                
        # Update the column counters (each column is a cancer):
            
            column_control += 1
            label_control  += 1
        
        # Update the row counter (each row is a GO):
        
        column_control = 0
        row_control   += 1
    
    # Add the legend:
    
    Closer_path     = mpatches.Patch(color='#bd3751', label='Closer')
    Farther_patch   = mpatches.Patch(color='#2679ea', label='Farther')
    Stable_patch    = mpatches.Patch(color='#50dc1f', label='Stable')
                                     
    plt.legend(handles=[Closer_path, Farther_patch,Stable_patch],borderaxespad=0.1,
                            bbox_to_anchor=(-1, -0.8), loc="lower center", frameon=True, ncol = 3)
    
    # Save the plot:
        
    plt.savefig(f'{movement_path}PanCancer_Movement.png', format="png", dpi=600,
                    bbox_inches='tight')

def Plot_top_movement_space(Cancer_type_list, Normal_Tissue_list, Cell_type_list, optimal_dim,
                            matrix = "PPMI", annotation = "Leaf", top = 10):  
    
    # Paths:
    
    cosine_path = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    common_path = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    save_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the cosine matrices:
        
        Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_{annotation}.csv',         index_col = 0)
        Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_{annotation}.csv', index_col = 0)
        
        # read the common set:
        
        common = list(pd.read_csv(f'{common_path}Common_Set_{tissue}_Leaf.csv')["0"])
        
        Cancer_structure   = Cancer_structure.loc[common][common]
        Control_structure  = Control_structure.loc[common][common]
              
        # Labels for the plot:
        
        plot_labels = []
    
        for i in range(top):
        
            plot_labels.append(ascii_lowercase[i])
        
        # Load the top 100 GO terms that move:
        
        moving = pd.read_csv(f'/media/sergio/sershiosdisk/Cancer/Annotations_Movement/top_100_moving_{tissue}.csv', index_col = 0)
        
        # Control Axes, Columns, and titles:
            
        row_control    = 0
        column_control = 0
        label_control  = 0
            
        # Start the plot:
            
        plt.style.use("seaborn-whitegrid")
        plt.rcParams.update({'font.size': 20})
        plt.rc('font', weight='bold')
            
        fig, axes = plt.subplots(4, 3, figsize=(15,15))
        fig.tight_layout()
        
        # For each GO terms in the top 10:
        
        for GO1 in list(moving.index)[0:top]:
            
            print(GO1)
            
            # Get the information for the plot:
            
            Cancer_structure_GO1  = Cancer_structure.loc[GO1]
            Control_structure_GO1 = Control_structure.loc[GO1]
            
            moving_GO1 = Control_structure_GO1 - Cancer_structure_GO1
                                    
            # Colors:
            
            erase = []
            erase2 = []
            
            colors = []
            
            n_dots = len(moving_GO1)   # set number of GO terms
            angs = np.linspace(0, 2*np.pi, n_dots)  # angles of the GO terms to keep the circle.
            cx, cy = (0, 0)  # center of circle
            xs, ys = [], []    # for coordinates of points to plot
            ra = 1         # radius of circle
            
            count = 0
            for ang in angs:
                # radius:
                ra_it = ra - moving_GO1[count]
                
                if ra_it > (1 + 2*np.std(moving_GO1)):
                    colors.append("#2679ea")
                                  
                    
                    erase.append(list(moving_GO1.index)[count])
                                  
                                  
                elif ra_it < (1 - 2*np.std(moving_GO1)):
                    colors.append("#bd3751")
                                  
                    erase2.append(list(moving_GO1.index)[count])
                else:
                    colors.append("#50dc1f")
                x = cx + ra_it*np.cos(ang)
                y = cy + ra_it*np.sin(ang)
                xs.append(x)   # collect x
                ys.append(y)   # collect y
                count+=1
                
            # Shapes:
            
            shapes = []
            size   = []
            
            sorted_GO = moving_GO1.sort_values(ascending = False)
            
            maximum = sorted_GO.head(5)
            minimum = sorted_GO.tail(5)
            
            for i in moving_GO1:
                if (sum(i == maximum) != 0) | (sum(i == minimum) != 0) :
                    shapes.append("^")
                    size.append(50)
                else:
                    shapes.append("o")
                    size.append(9)
            
            for i in range(len(xs)):
                axes[row_control, column_control].scatter(xs[i], ys[i], s=size[i] ,c=colors[i], marker = shapes[i])  # plot points 
            circle = plt.Circle([0,0], 0.05, color = "#04fbe1")
            axes[row_control, column_control].add_patch(circle)
            axes[row_control, column_control].spines['right'].set_visible(False)
            axes[row_control, column_control].spines['top'].set_visible(False)
            axes[row_control, column_control].spines['left'].set_linewidth(4)
            axes[row_control, column_control].spines['bottom'].set_linewidth(4) 
            axes[row_control, column_control].spines['left'].set_color("grey")
            axes[row_control, column_control].spines['bottom'].set_color("grey")
            axes[row_control, column_control].set_yticklabels(" ")
            axes[row_control, column_control].set_xticklabels(" ") 
            axes[row_control, column_control].set_xticks([])
            axes[row_control, column_control].set_yticks([])
            axes[row_control, column_control].set_title(f'{plot_labels[label_control].capitalize()}) {GO1}', fontweight='bold', 
            fontsize = 17, y = 1, x = 0.5)
            
        
            if column_control == 2:
                column_control = 0
                row_control +=1
                label_control +=1
            else:
                column_control +=1
                label_control +=1
            
            # For the empty sets:
            
            if (column_control == 0) & (row_control == 3):
                axes[row_control, column_control].axis('off')
                column_control +=1
                
            if (column_control == 2) & (row_control == 3):
                axes[row_control, column_control].axis('off')
                column_control +=1
                
        Closer_path     = mpatches.Patch(color='#bd3751', label='Closer')
        Farther_patch   = mpatches.Patch(color='#2679ea', label='Farther')
        Stable_patch    = mpatches.Patch(color='#50dc1f', label='Stable')
                                 
        plt.legend(handles=[Closer_path, Farther_patch,Stable_patch],borderaxespad=0.1,
                            bbox_to_anchor=(0.20, 0.8), loc="upper center", frameon=True, ncol = 1)
        
        plt.savefig(f'{save_path}Top10_Movement_{matrix}_{annotation}_{tissue}.png', format="png", dpi=600,
                    bbox_inches='tight')
            
        
def Submit_Query(GO_list, savepath):
    
    # Set the preferences:

    fp = webdriver.FirefoxProfile()
    fp.set_preference("browser.download.folderList", 2)
    fp.set_preference("browser.download.manager.showWhenStarting", False)
    fp.set_preference("browser.download.dir", savepath)
    fp.set_preference('browser.helperApps.neverAsk.saveToDisk', "text/plain, application/vnd.ms-excel, text/csv, text/comma-separated-values, application/octet-stream")

    # Start the browser:

    binary = FirefoxBinary('/usr/bin/firefox')
    driver = webdriver.Firefox(firefox_binary = binary, firefox_profile = fp)
    driver.get("http://revigo.irb.hr/")
    assert "Revigo" in driver.title
    
    # Accept the terms:
    
    press_button = driver.find_element_by_css_selector(
        "html body div.ui-dialog.ui-corner-all.ui-widget.ui-widget-content.ui-front.no-close.ui-dialog-buttons div.ui-dialog-buttonpane.ui-widget-content.ui-helper-clearfix div.ui-dialog-buttonset button.ui-button.ui-corner-all.ui-widget"
         )
    press_button.click()
    
    # Fill the GO terms:
    
    for GO in GO_list:
        
        text_area = driver.find_element_by_css_selector("html body form#aspnetForm div div div textarea#ctl00_MasterContent_txtGOInput")
        text_area.send_keys(str(GO))
        text_area.send_keys("\n")
        
    # Click the small list option:
    
    press_button = driver.find_element_by_css_selector(
        "html body form#aspnetForm div div div p input#ctl00_MasterContent_chkSimilarity0_5"
         )
    press_button.click()
    
    # Submit the job:
    
    press_button_sub = driver.find_element_by_css_selector("html body form#aspnetForm div div div input#ctl00_MasterContent_btnStart.ui-button.ui-corner-all.ui-widget")
    press_button_sub.click()
    
    # Download the file:
    
    time.sleep(30)
        
    link = driver.find_element_by_link_text('Export to CSV1,2')
    link.click()
    
    # Wait to download and quit the browser:
    
    time.sleep(20)
    driver.quit()


def Get_Domains_Cancer(Cancer_type_list, Normal_Tissue_list, Cell_type_list, matrix = "PPMI",
                               annotation = "Leaf"):
    
    # paths:
    
    rank_path       = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    
    numbers = [185, 173, 197, 142]
    x = 0
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load the ranking: 
        
        if annotation != "GO":
        
            Ranking = pd.read_csv(f'{rank_path}Rank_movement_{tissue}_{matrix}_{annotation}.csv')
        
        else:
            
            Ranking = pd.read_csv(f'{rank_path}Rank_movement_{tissue}_{matrix}.csv')
        
        # Top number ranking:
        
        number = numbers[x]
        top_ranking = list(Ranking["Unnamed: 0"][:number])
        x+=1
        
        print(number)
        
        # get the fucntional domains:
        
        save_path = f'{functional_path}Fucntional_Domain_top_{number}_{tissue}_{annotation}.csv'
        
        Submit_Query(top_ranking,save_path)

def Demonstrate_Sets(n_domains = 15, hallmark_thresh = 0.6):
    
    path            = "/media/sergio/sershiosdisk/Cancer/Cancer_GO/Revigo(2).csv"
    hallmarks_path  = "/media/sergio/sershiosdisk/Cancer/Cancer_GO/"
    
    # Open file:
        
    file = pd.read_csv(path, sep = " ")
    term_ids = []
    
    for line in range(len(file)):
        file_line = file.loc[line]
        if file_line["Eliminated"] == True:
            temp_representative = file_line["Representative,"]
            term_ids.append(temp_representative)
        else:
            term_ids.append(file_line["TermID,"])    
    
    # Count the GO terms
        
    word = "GO"
    real_terms = []
    for line in range(len(term_ids)):
        Term = term_ids[line]
        if word in Term:
            real_terms.append(Term)
        else:
            new_term = real_terms[len(real_terms) - 1]
            real_terms.append(new_term)
    
        # Count them:
        
    pd.DataFrame(real_terms).count(0)
        
    Counter_dic = Counter(real_terms)
    Count_pd    = pd.DataFrame.from_dict(Counter_dic, orient='index')
    Count_pd    = Count_pd.sort_values(0, ascending = False)
    
        
    GO_filter     = list(Count_pd.index)
    file.index    = file["TermID,"]
    GO_decription = list(file.loc[GO_filter]["Name,"])
    Counter_list  = Count_pd[0]
        
    domains = pd.DataFrame({"GO": GO_filter, "Description" : GO_decription, 
                                     "Counter" : Counter_list }).reset_index(drop = True)
                 
    # Get Hallmarks info:

    Clusters = Give_clusters(f'/media/sergio/sershiosdisk/Cancer/Cancer_GO/Revigo(2).csv')
    Clusters = Clusters.replace(',','', regex=True)

    hallmarks = pd.read_csv(f'{hallmarks_path}Hallmarks_Cancer_GO_list.csv' , index_col=0)
    
    # Compute the similarity:
    
    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)
    
    result       = pd.DataFrame(0., index=np.arange(len(Clusters.GO_Term)), columns=list(set(hallmarks.Hallmark)))
    result.index = Clusters.GO_Term
    
    for hal in set(hallmarks.Hallmark):
        
        compare_hal = hallmarks[hallmarks.Hallmark == hal].GO
        dist_it     = []
        
        for GO1 in Clusters.GO_Term:
            dist_it = []
            for GO2 in compare_hal:
                try:
                    dist_it.append(similarity.lin(G, GO1, GO2))
                except Exception as PGSSLookupError:
                    dist_it.append(0)
                    continue
            result.loc[GO1][hal] = max(dist_it)
    
    list_hall = []
    value     = []
    
    for GO_hall in result.index:
        
        indice = list(result.loc[GO_hall][result.loc[GO_hall] == max(result.loc[GO_hall])].index)
        
        if len(indice) > 1:
            list_hall.append([indice[0]])
        else:
            list_hall
            list_hall.append(indice)
            
        value.append(max(result.loc[GO_hall]))
    
    # Generate the DataFrame:
        
    list_hall = [val for sublist in list_hall for val in sublist]
    classification = pd.DataFrame({"GO" : result.index, "Hallmark" : list_hall, "value" : value})  

    # Plot:

    domain_pan          = domains.replace(',','', regex=True)
    domain_pan_examples = domains.iloc[:n_domains]
    domain_pan_examples = domain_pan_examples.replace(',','', regex=True)
    Clusters            = Clusters.replace(',','', regex=True)
    
    # Load the Hallmark information:
    
    hallmarks = classification[classification.value >= hallmark_thresh]
    
    # Link hallmarks to the Clusters:
    
    Description_fin = []
    hallmark        = []
    
    for domain in domain_pan_examples.GO:
        
        # Get info:
        
        GO_it          = Clusters[Clusters.Cluster == domain].GO_Term
        hallamark_it   = list(hallmarks[hallmarks.GO.isin(GO_it)].Hallmark)
        Description    = domain_pan_examples[domain_pan_examples.GO == domain].Description[domain_pan_examples[domain_pan_examples.GO == domain].index[0]]
        
        # Adapt info:
        
        if len(hallamark_it) > 0:
        
            Description = [Description]*len(GO_it)
            
            to_add = ["Not Defined"] * int(len(Description) - len(hallamark_it))
            hallamark_it.extend(to_add)
            
            # Upload the info:
            
            Description_fin.extend(Description)
            hallmark.extend(hallamark_it)
        
        else:
            
            Description   = [Description]*len(GO_it)
            hallamark_it  = ["Not Defined"] * len(Description)
        
            # Upload the info:
            
            Description_fin.extend(Description)
            hallmark.extend(hallamark_it)    
    
    
    # Generate the data:
    
    Nested_Pie_Data = pd.DataFrame({"Description" : Description_fin, 
                                    "Hallmark" : hallmark})
            
    Nested_Pie_Data = Nested_Pie_Data.groupby(Nested_Pie_Data.columns.tolist()).size().reset_index().rename(columns={0:'n'})
    
    # Colors for the inner:
    
    cmap        = plt.get_cmap('tab20')
    colors      = [cmap(i) for i in np.linspace(0, 1, n_domains)]
    final_color = colors
    
    # Colors hallmarks:
                 
    colours_hallmarks = {'Deregulating Cellular Energetic': "#8958e5",
                         'Genome Instability and Mutation': '#587ce5',
                         'Not Defined': '#ffffff',
                         'Activating Invasion and Metasta': '#000000',
                         'Enabling Replicative Immortalit' : '#dcdada',
                         'Sustaining Proliferative Signal': '#75bc84',
                         'Inducing Angiogenesis': '#ce5d5d',
                         'Tumor promoting Inflammation': '#d8994c',
                          'Evading Growth Suppressor' : '#935810',
                          'Avoiding Immune Destruction': "pink",
                          'Resist Cell Death': "orange"}
    
    # Start the plot:
    
    hallmarks_out  = Nested_Pie_Data.groupby('Description').sum()
    domains_in     = Nested_Pie_Data.groupby(['Description', 'Hallmark']).sum()
    
    labels = domains_in.index.get_level_values(1)

    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(12,12))
    size = 0.3
    

    ax.pie(domains_in.values.flatten(), radius=1,
           colors = [colours_hallmarks[key] for key in labels],
           wedgeprops=dict(width=size, edgecolor='w'))
        
    text2 = ax.pie(hallmarks_out.values.flatten(), radius=1-size+0.1, colors = final_color ,
           wedgeprops=dict(width=size, edgecolor='w'),autopct='%1.1f%%', pctdistance=0.5)
    
    # Add the numbers:
    
    for perc in range(len(text2[1])):
        text2[1][perc].set_fontsize(18)
        text2[1][perc].set_fontweight("bold")
        
    Total           = sum(domain_pan_examples.Counter) 
    Real_Total      = sum(domain_pan.Counter)
    Real_Percentage = (Total*100)/ Real_Total 
    
    ax.text(-0.4,1.1, f'Coverage : {round(Real_Percentage,2)}%', fontsize = 20, weight="bold")
    
    # Legends:
    
    list_legend = list(domains_in.index.get_level_values(0))
    list_legend = list(dict.fromkeys(list_legend))
    
    Other_Stuff = pd.DataFrame({"Value" : list_legend, "Colors" : final_color})
    
    # Control the length of the descriptions:
        
    jump = "\n"
        
    for index in Other_Stuff["Value"].index:
        if len(Other_Stuff["Value"][index]) > 90:
              it_setence   = Other_Stuff["Value"][index].split()
              mid_pos      = len(it_setence)//2
              new_sentence = it_setence[:mid_pos] + [jump] + it_setence[mid_pos:]
              new_sentence = ' '.join(new_sentence) 
              Other_Stuff["Value"][index] = new_sentence
        
    # Create a legend:
        
    patch_list = []
        
    for patch in range(n_domains):
            
          patch_list.append(mpatches.Patch(color= final_color[patch], label= Other_Stuff["Value"][patch]))
    
    legend1 = plt.legend(handles=patch_list,borderaxespad=0.1,
                   bbox_to_anchor=(0.94, 0.85), loc="upper left", frameon=True, fontsize=14, 
                   title = r"$\bf{Fucntional\ Domain}$")
    
    # Add the second legend:
    
    patch_hallmarks = []
    
    for patch in colours_hallmarks.keys():
        
        if patch != "Not Defined":
            patch_hallmarks.append(mpatches.Patch(color= colours_hallmarks[patch], label= patch)) 
        else:
            continue
    
    legend2 =  plt.legend(handles=patch_hallmarks,borderaxespad=0.1,
                   bbox_to_anchor=(0.94, 0.25), loc="upper left", frameon=True, fontsize=14, 
                   title = r"$\bf{Cancer\ Hallmarks}$")
    
    ax.add_artist(legend1)
    ax.add_artist(legend2)
    
    plt.savefig("/media/sergio/sershiosdisk/Cancer/Cancer_GO/Hallmarks_Cancer_Set.png", format="png", dpi=600,
                    bbox_inches='tight')
    
    
        
def Produce_Domains(Cancer_type_list, Normal_Tissue_list, Cell_type_list, matrix = "PPMI",
                               annotation = "Leaf", number = 500):

    
    # Paths:
    
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"

    b = 0
    numbers = [185, 173, 197, 142]
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        number = numbers[b]
        
        
        domain_path = f'{functional_path}Fucntional_Domain_top_{number}_{tissue}_{annotation}/Revigo.csv' 
       
       # Open file:
       
        file = pd.read_csv(domain_path, sep = " ")
        term_ids = []
    
        for line in range(len(file)):
            file_line = file.loc[line]
            if file_line["Eliminated"] == True:
                temp_representative = file_line["Representative,"]
                term_ids.append(temp_representative)
            else:
                term_ids.append(file_line["TermID,"])
        
        # Count the GO terms
        
        word = "GO"
        real_terms = []
        for line in range(len(term_ids)):
            Term = term_ids[line]
            if word in Term:
                real_terms.append(Term)
            else:
                new_term = real_terms[len(real_terms) - 1]
                real_terms.append(new_term)
    
        # Count them:
        
        pd.DataFrame(real_terms).count(0)
        
        Counter_dic = Counter(real_terms)
        Count_pd    = pd.DataFrame.from_dict(Counter_dic, orient='index')
        Count_pd    = Count_pd.sort_values(0, ascending = False)
        
        # Final result:
        
        GO_filter     = list(Count_pd.index)
        file.index    = file["TermID,"]
        GO_decription = list(file.loc[GO_filter]["Name,"])
        Counter_list  = Count_pd[0]
        
        Final_pandas = pd.DataFrame({"GO": GO_filter, "Description" : GO_decription, 
                                     "Counter" : Counter_list }).reset_index(drop = True)
            
        Final_pandas.to_csv(f'{functional_path}Domains_{tissue}_{matrix}_{annotation}_{number}.csv') 
        
        b = b+1

def Plot_Domains(Cancer_type_list, Normal_Tissue_list, Cell_type_list, matrix = "PPMI",
                               annotation = "Leaf", number = 400, n_domains = 10):
    
    # Paths:
    
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    
    row_control = 0
    label_count = 0

    # Authomatic Labels:
    
    plot_labels = []
    
    for i in range(len(Cancer_type_list * 3)):
        plot_labels.append(ascii_lowercase[i])
    
    # Plot features:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    
    fig, axes = plt.subplots(len(Cancer_type_list), 1, figsize=(18,18))
    fig.tight_layout()

    b = 0
    numbers = [185, 173, 197, 142]
    
    x = domain["Description"]
    x = list(x)
    n = 0
    for i in x:
        if "reg" in i:
            n = n +1
    
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        number = numbers[b]
        
        domain = pd.read_csv(f'{functional_path}Domains_{tissue}_{matrix}_{annotation}_{number}.csv', index_col = 0)
        
        domains_representative = domain.iloc[:n_domains]
        domains_representative = domains_representative.replace(',','', regex=True)
        
        # Information for the PiePlot:
        
        Total           = sum(domains_representative.Counter) 
        percentages     = [i*100/Total for i in domains_representative.Counter]
        Real_Total      = sum(domain.Counter)
        Real_Percentage = (Total*100)/ Real_Total # The coverag
        
        # Pie_Plot:
        
        # Colors:
        
        cmap        = plt.get_cmap('tab20')
        colors      = [cmap(i) for i in np.linspace(0, 1, n_domains)]
        final_color = colors
        
        # Plot:
        
        patches, texts, autotexts = axes[row_control].pie(percentages, autopct='%1.f%%', shadow=False, colors= final_color)
        axes[row_control].set_title(plot_labels[label_count].capitalize(), fontweight='bold', fontsize = 17, y =0.5, x = -0.01)
        circle = plt.Circle(xy=(0,0), radius=0.75, facecolor='white')
        
        axes[row_control].add_patch(circle)
        #plt.gca().add_artist(circle)
        
        for perc in range(len(autotexts)):
            autotexts[perc].set_fontsize(14)
            autotexts[perc].set_fontweight("bold")
        axes[row_control].text(-0.7,1.15,f'Coverage : {round(Real_Percentage,2)}%', fontsize = 15)
        
        # Legend:
        
        list_description = []
        for description in range(n_domains):
            list_description.append(domains_representative.iloc[description]["Description"])
            
        Order_Stuff = pd.DataFrame({"Value" : list_description, "Colors" : final_color})
        
        
        # Control the length of the descriptions:
        
        jump = "\n"
        
        for index in Order_Stuff["Value"].index:
            if len(Order_Stuff["Value"][index]) > 90:
                it_setence   = Order_Stuff["Value"][index].split()
                mid_pos      = len(it_setence)//2
                new_sentence = it_setence[:mid_pos] + [jump] + it_setence[mid_pos:]
                new_sentence = ' '.join(new_sentence)
                
                Order_Stuff["Value"][index] = new_sentence
        
        # Create a legend:
        
        patch_list = []
        
        for patch in range(n_domains):
            
           patch_list.append(mpatches.Patch(color= final_color[patch], label= Order_Stuff["Value"][patch]))
    
        axes[row_control].legend(handles=patch_list,borderaxespad=0.1,
                   bbox_to_anchor=(0.95, 1), loc="upper left", frameon=True)
        
        row_control+=1
        label_count+=1
        b = b +1
        
    fig.savefig(f'{functional_path}Functional_Domains_{matrix}_{number}_{annotation}.png', format="png", dpi=300,
                    bbox_inches='tight') 
        
        
def Analyze_Common_Set(Cancer_type_list, Normal_Tissue_list, Cell_type_list, matrix = "PPMI",
                               annotation = "Leaf", number = 500):
    
    # Paths:
    
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    
    # Do the dictionary:
    
    ranking_dictionary = {}
    for tissue in Normal_Tissue_list:
        if annotation != "GO":
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI_{annotation}.csv')
        else:
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI.csv')
        
        top                        = set(Ranking["Unnamed: 0"][:number])
        ranking_dictionary[tissue] = top
        
    # Calculate the intersections:
    
    intersect = distance.squareform([len(ranking_dictionary[u].intersection(ranking_dictionary[v])) for u,v in combinations(ranking_dictionary.keys(),2)])
    df        = pd.DataFrame(intersect,index=ranking_dictionary.keys(),columns=ranking_dictionary.keys())
    np.fill_diagonal(df.values, number)
    df = df/500
    
    # Plot the Intersection:
    
    plt.figure(figsize=(10,8))
    plot = sns.heatmap(df,annot=True, cmap = "viridis", linewidths=.5)
    plot.set_xticklabels(plot.get_xmajorticklabels(), fontsize = 18)
    plot.set_yticklabels(plot.get_ymajorticklabels(), fontsize = 18)
    cbar = plot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=18)
    
    # Save Plot:
    
    plt.savefig(f'{functional_path}Common_Cancer_{matrix}_{number}_{annotation}.png', format="png", dpi=600,
                    bbox_inches='tight') 

def Domains_Pancer_Cancer(Cancer_type_list, Normal_Tissue_list, Cell_type_list, matrix = "PPMI",
                               annotation = "Leaf", number = 500):
    
    # Paths:
    
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    
    # Do the dictionary:
    
    ranking_dictionary = {}
    for tissue in Normal_Tissue_list:
        if annotation != "GO":
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI_{annotation}.csv')
        else:
            Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI.csv')
        
        top = set(Ranking["Unnamed: 0"][:number])
        ranking_dictionary[tissue] = top
    
    # Common PanCancer:
    
    common_pan_cancer = ranking_dictionary['breast'].intersection(ranking_dictionary['prostate'],
                                          ranking_dictionary['lung'],ranking_dictionary['colon'])
    
    common_pan_cancer = ['GO:0006970', 'GO:0007254','GO:0008543','GO:0009314','GO:0016572','GO:0046833','GO:0051403','GO:1990869']
    
    # Send to create the Domains:
    
    save_path = f'{functional_path}PanCancer_Domains_{matrix}_{number}_{annotation}'
    
    Submit_Query(common_pan_cancer, save_path)
    
    
def Domain_Pan_Cancer_Connection(Cancer_list, Normal_Tissue_list, Control_list, number = 500, matrix = "PPMI", annotation = "Leaf", optimal_dim = 200):
    
    # Paths:
    
    movement_path   = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    common_path     = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    cosine_path     = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"

    # Get the pan-cancer annotations:

    ranking_dictionary = {}
    for tissue in Normal_Tissue_list:
        Ranking = pd.read_csv(f'{movement_path}Rank_movement_{tissue}_PPMI_Leaf.csv')
        top     = set(Ranking["Unnamed: 0"][:number])
        ranking_dictionary[tissue] = top
    
    common_pan_cancer = ranking_dictionary['breast'].intersection(ranking_dictionary['prostate'],
                                          ranking_dictionary['lung'],ranking_dictionary['colon']) 

    # Load the information about all cancers:
    
    moving_list = []
    
    for cancer, tissue, cell in zip(Cancer_list, Normal_Tissue_list,Control_list):
        Cancer_structure  = pd.read_csv(f'{cosine_path}Cosine_Cancer_{cancer}_{optimal_dim}_{matrix}_{annotation}.csv',         index_col = 0)
        Control_structure = pd.read_csv(f'{cosine_path}Cosine_Control_{tissue}_{cell}_{optimal_dim}_{matrix}_{annotation}.csv', index_col = 0)
        common = list(pd.read_csv(f'{common_path}Common_Set_{tissue}_{annotation}.csv')["0"])
        
        Cancer_structure  = Cancer_structure.loc[common][common]
        Control_structure = Control_structure.loc[common][common]
        
        moving = Control_structure  - Cancer_structure
        moving_list.append(moving)
       
    # Now we evaluate the percentage of interactions that are disconnected or connected in all cancers:
    
    # Open the file to write:
    
    with open(movement_path + "Movement_Pancancer_Agreement.txt", 'a') as the_file:
         
        # Write the column names:
         
        the_file.write("# Annotation"  + "\t")
        the_file.write("# Connected"    + "\t")
        the_file.write("# Disconnected" + "\t")
        the_file.write("# Both"         + "\n") 
        
        for Pann_GO in common_pan_cancer:
            
            connected_list_fin    = []
            disconnected_list_fin = []
            both_list_fin         = []
            
            for i in range(len(moving_list)):
                
                connected_list    = []
                disconnected_list = []
                both_list         = []
                
                moving_cancer = moving_list[i]
            
                # Load the movement of the pann cancer annotations
            
                moving_pann = moving_cancer.loc[Pann_GO] 
            
                # Get info (this can be changed but i am keeping the same code to make it faster):
                        
                n_dots = len(moving_pann)
                angs   = np.linspace(0, 2*np.pi, n_dots)
                ra     = 1
            
                # Analyze things that are connnected or disconnected:
                
                count = 0         
                for ang in angs:
                    ra_it = ra - moving_pann[count]
                    
                    # Disconnected:
                    
                    if ra_it > (1 + 2*np.std(moving_pann)):
                        disconnected_list.append(list(moving_pann.index)[count])
                        both_list.append(list(moving_pann.index)[count])
                        count = count + 1
                    
                    # Connected:
                    
                    elif ra_it < (1 - 2*np.std(moving_pann)):
                        connected_list.append(list(moving_pann.index)[count])
                        both_list.append(list(moving_pann.index)[count])   
                        count = count + 1
                    
                    else:
                        count = count + 1
                        
                
                # Update the list per each cancer we will have one list:
                
                connected_list_fin.append(connected_list)
                disconnected_list_fin.append(disconnected_list)
                both_list_fin.append(both_list)  
                
            # Do the calculations:
                
            Total_connected    = len(set.union(*map(set,connected_list_fin)))
            Total_disconnected = len(set.union(*map(set,disconnected_list_fin)))
            Total_both         = len(set.union(*map(set,both_list_fin))) 
            
            Int_connected    = len(set.intersection(*map(set,connected_list_fin)))
            Int_disconnected = len(set.intersection(*map(set,disconnected_list_fin)))
            Int_both         = len(set.intersection(*map(set,both_list_fin))) 
            
            # Write the information to a file:            
                        
            the_file.write(f'{Pann_GO}\t')
            the_file.write(f'{Int_connected} ({Total_connected}, {round((Int_connected/Total_connected)*100,2)}%)\t')
            the_file.write(f'{Int_disconnected} ({Total_disconnected}, {round((Int_disconnected/Total_disconnected)*100,2)}%)\t')
            the_file.write(f'{Int_both} ({Total_both}, {round((Int_both/Total_both)*100,2)}%)\n')
  
    the_file.close()    
    
     
    
def Produce_Domains_PanCancer(matrix = "PPMI", annotation = "Leaf", number = 500):

    # Paths:
    
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    domain_path = f'{functional_path}PanCancer_Domains_{matrix}_{number}_{annotation}/Revigo.csv'
    
    file     = pd.read_csv(domain_path, sep = " ")
    term_ids = []
    
    for line in range(len(file)):
        file_line = file.loc[line]
        if file_line["Eliminated"] == True:
            temp_representative = file_line["Representative,"]
            term_ids.append(temp_representative)
        else:
            term_ids.append(file_line["TermID,"])    
    
        # Count the GO terms
        
    word = "GO"
    real_terms = []
    for line in range(len(term_ids)):
        Term = term_ids[line]
        if word in Term:
            real_terms.append(Term)
        else:
            new_term = real_terms[len(real_terms) - 1]
            real_terms.append(new_term)
    
    # Count them:
        
    pd.DataFrame(real_terms).count(0)
        
    Counter_dic = Counter(real_terms)
    Count_pd    = pd.DataFrame.from_dict(Counter_dic, orient='index')
    Count_pd    = Count_pd.sort_values(0, ascending = False)
        
    # Final result:
    
    GO_filter     = list(Count_pd.index)
    file.index    = file["TermID,"]
    GO_decription = list(file.loc[GO_filter]["Name,"])
    Counter_list  = Count_pd[0]
        
    Final_pandas = pd.DataFrame({"GO": GO_filter, "Description" : GO_decription, 
                                     "Counter" : Counter_list }).reset_index(drop = True)
            
    Final_pandas.to_csv(f'{functional_path}Domains_PannCancer_{matrix}_{number}_{annotation}.csv') 
    
def Give_clusters(path):
    
    # Load the Axes information:
    
    file = pd.read_csv(path, sep = " ")
    
    # Create the clusters:
    
    term_ids = []
    
    for line in range(len(file)):
        file_line = file.loc[line]
        if file_line["Eliminated"] == True:
            temp_representative = file_line["Representative,"]
            term_ids.append(temp_representative)
        else:
            term_ids.append(file_line["TermID,"])
        
    # Clasfy the terms
    
    clusters = []
    word = "GO"
    for line in range(len(term_ids)):
        Term = term_ids[line]
        if word in Term:
            clusters.append(Term)
        else:
            new_term = clusters[len(clusters) - 1]
            clusters.append(new_term)
            
    # Cluster dataframe:
    
    Cluster_pd = pd.DataFrame({"GO_Term" : file["TermID,"], "Cluster" : clusters})
    
    # Return the information
    
    return(Cluster_pd)
    
def Classify_Hallmarks(cancer, matrix = "PPMI", annotation = "Leaf", number = 500, n_domains=10, colon = False):
    
    # Paths:
    
    functional_path      = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    hallmarks_path       = "/media/sergio/sershiosdisk/Cancer/Cancer_GO/"
    
    # Read the information needed:
    
 
    Clusters = Give_clusters(f'{functional_path}PanCancer_Domains_{matrix}_{number}_{annotation}_NoColon/Revigo.csv')
    Clusters = Clusters.replace(',','', regex=True)

    hallmarks = pd.read_csv(f'{hallmarks_path}Hallmarks_Cancer_GO_list.csv' , index_col=0)
    
    # Compute the similarity:
    
    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)
    
    result       = pd.DataFrame(0., index=np.arange(len(Clusters.GO_Term)), columns=list(set(hallmarks.Hallmark)))
    result.index = Clusters.GO_Term
    

    for hal in set(hallmarks.Hallmark):
        
        compare_hal = hallmarks[hallmarks.Hallmark == hal].GO
        dist_it     = []
        
        for GO1 in Clusters.GO_Term:
            dist_it = []
            for GO2 in compare_hal:
                try:
                    dist_it.append(similarity.lin(G, GO1, GO2))
                except Exception as PGSSLookupError:
                    continue
            result.loc[GO1][hal] = max(dist_it)
    
    list_hall = []
    value     = []
    
    for GO_hall in result.index:
        
        indice = list(result.loc[GO_hall][result.loc[GO_hall] == max(result.loc[GO_hall])].index)
        
        if len(indice) > 1:
            list_hall.append([indice[0]])
        else:
            list_hall
            list_hall.append(indice)
            
        value.append(max(result.loc[GO_hall]))
    
    # Generate the DataFrame:
        
    list_hall = [val for sublist in list_hall for val in sublist]
    classification = pd.DataFrame({"GO" : result.index, "Hallmark" : list_hall, "value" : value})
    
    # Save the relation Matrix:
    
    classification.to_csv((f'{functional_path}Hallmarks_PanCancer_Links_{matrix}_{number}_{annotation}.csv'))
    

def Plot_PanCancer_Domains_Hallmarks(matrix = "PPMI", annotation = "Leaf", number = 500, n_domains=10, colon = False,
                                     hallmark_thresh = 0.5):   
    
     
    # Paths:
    
    functional_path = "/media/sergio/sershiosdisk/Cancer/Functional_Domains/"
    
    # Load the basic information:
   
    if colon == False:
        
        domain_pan     = pd.read_csv(f'{functional_path}Domains_PannCancer_{matrix}_{number}_{annotation}.csv', 
                                 index_col = 0)
        Clusters       = Give_clusters(f'{functional_path}PanCancer_Domains_{matrix}_{number}_{annotation}/Revigo(13).csv')
                
        classification = pd.read_csv(f'{functional_path}Hallmarks_PanCancer_Links_{matrix}_{number}_{annotation}.csv', 
                                 index_col = 0)
    else:
        
        domain_pan = pd.read_csv(f'{functional_path}Domains_PannCancer_{matrix}_{number}_{annotation}NoColon.csv', 
                                 index_col = 0)
        Clusters   = Give_clusters(f'{functional_path}PanCancer_Domains_{matrix}_{number}_{annotation}_NoColon/Revigo.csv')
        
        classification = pd.read_csv(f'{functional_path}Hallmarks_PanCancer_Links_{matrix}_{number}_{annotation}.csv', 
                                 index_col = 0)
    
    # Prepare the data:

    domain_pan = domain_pan.replace(',','', regex=True)
    domain_pan_examples = domain_pan.iloc[:n_domains]
    domain_pan_examples = domain_pan_examples.replace(',','', regex=True)
    Clusters            = Clusters.replace(',','', regex=True)
    
    # Load the Hallmark information:
    
    hallmarks = pd.read_csv(f"{functional_path}Hallmarks_PanCancer_Links_PPMI_500_Leaf.csv", index_col=0)
    hallmarks = classification[classification.value >= hallmark_thresh]
    
    # Link hallmarks to the Clusters:
    
    Description_fin = []
    hallmark        = []
    
    for domain in domain_pan_examples.GO:
        
        # Get info:
        
        GO_it          = Clusters[Clusters.Cluster == domain].GO_Term
        hallamark_it   = list(hallmarks[hallmarks.GO.isin(GO_it)].Hallmark)
        Description    = domain_pan_examples[domain_pan_examples.GO == domain].Description[domain_pan_examples[domain_pan_examples.GO == domain].index[0]]
        
        # Adapt info:
        
        if len(hallamark_it) > 0:
        
            Description = [Description]*len(GO_it)
            
            to_add = ["Not Defined"] * int(len(Description) - len(hallamark_it))
            hallamark_it.extend(to_add)
            
            # Upload the info:
            
            Description_fin.extend(Description)
            hallmark.extend(hallamark_it)
        
        else:
            
            Description   = [Description]*len(GO_it)
            hallamark_it  = ["Not Defined"] * len(Description)
        
            # Upload the info:
            
            Description_fin.extend(Description)
            hallmark.extend(hallamark_it)  
            
    # Generate the data:
    
    Nested_Pie_Data = pd.DataFrame({"Description" : Description_fin, 
                                    "Hallmark" : hallmark})
            
    Nested_Pie_Data = Nested_Pie_Data.groupby(Nested_Pie_Data.columns.tolist()).size().reset_index().rename(columns={0:'n'})
    
    # Colors for the inner:
    
    cmap        = plt.get_cmap('tab20')
    colors      = [cmap(i) for i in np.linspace(0, 1, n_domains)]
    final_color = colors
    
    # Colors hallmarks:
                 
    colours_hallmarks = {'Deregulating Cellular Energetic': "#8958e5",
                         'Genome Instability and Mutation': '#587ce5',
                         'Not Defined': '#ffffff',
                         'Activating Invasion and Metasta': '#000000',
                         'Enabling Replicative Immortalit' : '#dcdada',
                         'Sustaining Proliferative Signal': '#75bc84',
                         'Inducing Angiogenesis': '#ce5d5d',
                         'Tumor promoting Inflammation': '#d8994c',
                          'Evading Growth Suppressor' : '#935810'}
    
    # Start the plot:
    
    hallmarks_out  = Nested_Pie_Data.groupby('Description').sum()
    domains_in     = Nested_Pie_Data.groupby(['Description', 'Hallmark']).sum()
    
    labels = domains_in.index.get_level_values(1)

    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(12,12))
    size = 0.3
    

    ax.pie(domains_in.values.flatten(), radius=1,
           colors = [colours_hallmarks[key] for key in labels],
           wedgeprops=dict(width=size, edgecolor='w'))
        
    text2 = ax.pie(hallmarks_out.values.flatten(), radius=1-size+0.1, colors = final_color ,
           wedgeprops=dict(width=size, edgecolor='w'),autopct='%1.1f%%', pctdistance=0.4)
    
    # Add the numbers:
    
    for perc in range(len(text2[1])):
        text2[1][perc].set_fontsize(18)
        text2[1][perc].set_fontweight("bold")
        
    Total           = sum(domain_pan_examples.Counter) 
    Real_Total      = sum(domain_pan.Counter)
    Real_Percentage = (Total*100)/ Real_Total 
    
    ax.text(-0.4,1.1, f'Coverage : {round(Real_Percentage,2)}%', fontsize = 20, weight="bold")
    
    # Legends:
    
    list_legend = list(domains_in.index.get_level_values(0))
    list_legend = list(dict.fromkeys(list_legend))
        
    Other_Stuff = pd.DataFrame({"Value" : list_legend, "Colors" : final_color})
    
    # Control the length of the descriptions:  
        
    jump = "\n"
        
    for index in Other_Stuff["Value"].index:
        if len(Other_Stuff["Value"][index]) > 90:
              it_setence   = Other_Stuff["Value"][index].split()
              mid_pos      = len(it_setence)//2
              new_sentence = it_setence[:mid_pos] + [jump] + it_setence[mid_pos:]
              new_sentence = ' '.join(new_sentence) 
              Other_Stuff["Value"][index] = new_sentence
        
    # Create a legend:
        
    patch_list = []
        
    for patch in range(n_domains):
            
          patch_list.append(mpatches.Patch(color= final_color[patch], label= Other_Stuff["Value"][patch]))
    
    legend1 = plt.legend(handles=patch_list,borderaxespad=0.1,
                   bbox_to_anchor=(0.94, 0.85), loc="upper left", frameon=True, fontsize=14, 
                   title = r"$\bf{Fucntional\ Domain}$")
    
    # Add the second legend:
    
    patch_hallmarks = []
    
    for patch in colours_hallmarks.keys():
        
        if patch != "Not Defined":
            patch_hallmarks.append(mpatches.Patch(color= colours_hallmarks[patch], label= patch)) 
        else:
            continue
    
    legend2 =  plt.legend(handles=patch_hallmarks,borderaxespad=0.1,
                   bbox_to_anchor=(0.94, 0.35), loc="upper left", frameon=True, fontsize=14, 
                   title = r"$\bf{Cancer\ Hallmarks}$")
    
    ax.add_artist(legend1)
    ax.add_artist(legend2)
    
    # Save the plot:
    
    plt.savefig(f'{functional_path}Hallmarks_PanCancer_Domains_{n_domains}_TEST_3.png', format="png", dpi=600,
                    bbox_inches='tight')       
    
def run_process(process):                                                             
    os.system(format(process)) 
    
def Run_Permutations():
    
    print("Running Permutations")
    
    command = "/home/sergio/Project_2/Scripts/python3 Terminal_Permutations_Cancer.py"
    run_process(command)
    
    print("Finished")  

            
def Prepare_Annotations(network_path, GO_Matrix, cancer, tissue, cell, Type = "Cancer"):
    
    if Type == "Cancer":
        genes      = pd.read_table(f'{network_path}Cancer_{cancer}_Genes.csv', header = 0, dtype={0: str})
        genes_list = list(genes["0"]) 
    else:
        genes      = pd.read_table(f'{network_path}Control_{tissue}_{cell}_Genes.csv', header = 0, dtype={0: str})
        genes_list = list(genes["0"])  
    
    GO_Matrix_filt   = GO_Matrix[GO_Matrix.index.isin(genes_list)]
    
    GO               = GO_Matrix_filt.sum(axis=0)
    filter_list      = set(GO[GO > 0].index)
    
    GO_Matrix_final  = GO_Matrix_filt[filter_list]
    GO_Matrix_final  = GO_Matrix_final.loc[genes_list]
    
    return(GO_Matrix_final)
        

def Permutation_pvalues(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dim, matrix, times = 60000, Type = "Cancer"):
    
    # Read the files.
    
    # GO Embeddings from the same k (at this point for back propagation).
    
    print("Cancer and Control Permutations")
    
    # Paths:
    
    save_path      = "/media/sergio/sershiosdisk/Cancer/Axes/"
    network_path   = "/media/sergio/sershiosdisk/Cancer/Networks/"
    GO_matrix_path = '/media/sergio/sershiosdisk/Human/Leaf_Annotation/_Matrix_Genes_GO_BP_PPI.csv'
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
    
        if Type == "Cancer":
            
            gene_path       = f'{network_path}_G_Matrix_{dim}_PPI_{cancer}_{matrix}_Cancer.npy'
            go_path         = f'{network_path}_GO_Embeddings_Leaf_PPI_{cancer}_{dim}_{matrix}_Cancer.csv'
            GO_matrix       = pd.read_csv(GO_matrix_path,index_col = 0, dtype={0: str})
            GO_matrix.index = GO_matrix.index.astype(str)
            GO_Matrix_Fin   = Prepare_Annotations(network_path,GO_matrix, cancer, tissue, cell, Type = "Cancer")
            
        else:
            
            gene_path       = f'{network_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_{matrix}_Control.npy'
            go_path         = f'{network_path}_GO_Embeddings_Leaf_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv'
            GO_matrix       = pd.read_csv(GO_matrix_path,index_col = 0, dtype={0: str})
            GO_matrix.index = GO_matrix.index.astype(str)
            GO_Matrix_Fin   = Prepare_Annotations(network_path,GO_matrix, cancer, tissue, cell, Type = "Control")
            
        # Load the matrices:
        
        Gene_embeddings = np.load(gene_path, allow_pickle=True)
        GO_emebddings   = pd.read_csv(go_path, index_col = 0 ,dtype={0: str})
        
        # Control the order of the GO terms:
        
        GO_embeddings = GO_emebddings[GO_Matrix_Fin.columns]
        
        # Start permutations:
        
        Count_db = pd.DataFrame(0, columns=GO_embeddings.columns, index = GO_embeddings.index)
        
        index_orig   = GO_Matrix_Fin.index
        columns_orig = GO_Matrix_Fin.columns
        
        #    Init the permutation test for N times:
            
        for rep in range(times):
            
            # Control de iteration:
           
            if rep%1000 == 0:
                print(rep)
            
            # Randomly suffle the Gene/GO matrix:
                
            Gene_GO_matrix_it = GO_Matrix_Fin.sample(frac=1).reset_index(drop=True)
            Gene_GO_matrix_it.index = index_orig
            Gene_GO_matrix_it.columns = columns_orig
                
            # Do the Embeddings using the suffled Gene/GO matrix (direct calculation):
            
            gene_embeddings_db_inverse = pd.DataFrame(np.linalg.pinv(Gene_embeddings) , columns = index_orig)
            GO_emebddings_it = gene_embeddings_db_inverse.dot(Gene_GO_matrix_it)
    
            # Compare the new scores vs the original scores:
                
            comparison = GO_emebddings_it >= GO_embeddings
                
            # Update the counts:
                
            Count_db = Count_db + comparison
            
        # Finished:
                
        print("the " + str(rep + 1) +  " iterations are finished")
                 
        # Calculate the p-values:
        
        Count_db_Final = (Count_db + 1)/(times + 1)
        
        # Save the matrix.
        
        Count_db_Final.to_csv(f'{save_path}Counts_Axes_{times}_{tissue}_{Type}_2.csv', header = True, index=True)  
        

def Get_Adjusted_P_values_Filtered(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dim, 
                                   Type = "Cancer", alpha = 0.05, times = 60000, matrix = "PPMI", plot = False):
        
    # Paths:
    
    save_path      = "/media/sergio/sershiosdisk/Cancer/Axes/60K/"
    network_path   = "/media/sergio/sershiosdisk/Cancer/Networks/"
    GO_matrix_path = '/media/sergio/sershiosdisk/Human/Leaf_Annotation/_Matrix_Genes_GO_BP_PPI.csv'
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
    
        if Type == "Cancer":
            
            go_path         = f'{network_path}_GO_Embeddings_Leaf_PPI_{cancer}_{dim}_{matrix}_Cancer.csv'
            GO_matrix       = pd.read_csv(GO_matrix_path,index_col = 0, dtype={0: str})
            GO_matrix.index = GO_matrix.index.astype(str)
            GO_Matrix_Fin   = Prepare_Annotations(network_path,GO_matrix, cancer, tissue, cell, Type = "Cancer")
            count_path      = f'{save_path}Counts_Axes_{times}_{tissue}_{Type}.csv'
            
        else:
            
            go_path         = f'{network_path}_GO_Embeddings_Leaf_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv'
            GO_matrix       = pd.read_csv(GO_matrix_path,index_col = 0, dtype={0: str})
            GO_matrix.index = GO_matrix.index.astype(str)
            GO_Matrix_Fin   = Prepare_Annotations(network_path,GO_matrix, cancer, tissue, cell, Type = "Control")
            count_path      = f'{save_path}Counts_Axes_{times}_{tissue}_{Type}.csv'
    
        # Load data:
        
        GO_emebddings = pd.read_csv(go_path, index_col = 0 ,dtype={0: str})
        GO_embeddings = GO_emebddings[GO_Matrix_Fin.columns]
        count_matrix  = pd.read_csv(count_path,index_col = 0, dtype={0: str})
        
        
        # Total number of GO terms:
        
        NMTF_matrix = GO_embeddings.T
        
        # Prepare an empty matrix:
        
        P_adj_Matrix = pd.DataFrame(0, columns = count_matrix.columns, index = count_matrix.index)
        
        # Information for the control plots:
        
        dimensions_list       = []
        go_scores_sign        = []
        go_scores_list_nosign = []
        
        for GO in count_matrix.columns:
            
            # Get the p-values:
            
            p_values_list = count_matrix[GO]
            
            # Correct the p-values:
    
            p_values_list_corrected = multipletests(p_values_list.values, alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)
            
            # Save the p-values adjusted:
            
            P_adj_Matrix[GO] = p_values_list_corrected[1]
            
        # If the user want the plot:
        
            if plot == True:
            
                # Get the p-values that are significant:
                
                corrected_list = [i for i in p_values_list_corrected[1] if i<= alpha]
                
                # Get the information using these values for lists that at least have one significative value:
        
                if len(corrected_list)>0:
                    
                    # Get the info and the scores:
                    
                    corrected_list_idx = [number for number,i in enumerate(p_values_list_corrected[1]) if i<= alpha]
                    GO_scores          = [NMTF_matrix.loc[GO,i] for i in corrected_list_idx]   
                    GO_scores_nosign   = [NMTF_matrix.loc[GO,i] for i in range(0,len(NMTF_matrix.columns)) if i not in corrected_list_idx] 
                    
                    dimensions_list.append(corrected_list_idx)
                    go_scores_sign.append(GO_scores)
                    go_scores_list_nosign.append(GO_scores_nosign)
                    
                else:
                    
                    GO_scores_nosign   = [NMTF_matrix.loc[GO,i] for i in range(0,len(NMTF_matrix.columns))] 
                    go_scores_list_nosign.append(GO_scores_nosign) 
        
        # Save the p-adjusted values:
        
        P_adj_Matrix.to_csv(f'{save_path}P_values_adjusted_{tissue}_{Type}.csv', header = True, index=True)     
    

def Associate_Annotations_Axes(Normal_Tissue_list, dim, counts = 60000):
    
    # Paths:
    
    Axes_path   = "/media/sergio/sershiosdisk/Cancer/Axes/"
    common_path = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    
    for tissue in Normal_Tissue_list:
        
        # Load the p-values:
    
        Control = pd.read_csv(f'{Axes_path}P_values_adjusted_{tissue}_Control.csv', index_col=0)
        Cancer  = pd.read_csv(f'{Axes_path}P_values_adjusted_{tissue}_Cancer.csv', index_col=0)
        
        # Load the common set of GO terms:
        
        common = list(pd.read_csv(f'{common_path}Common_Set_{tissue}_Leaf.csv')["0"])
        
        Control = Control[common]
        Cancer  = Cancer[common]
        
        # Init the dictionaries:
        
        Control_dic = dict.fromkeys(list(Control.index))
        Cancer_dic  = dict.fromkeys(list(Cancer.index))
    
        for dimension in Control.index:
                
            list_Assoc = list(Control.columns[Control.loc[dimension] <= 0.05])
                
            if not list_Assoc:
                Control_dic[dimension] = []
            else:
                Control_dic[dimension] = list_Assoc

        for dimension in Cancer.index:
                
            list_Assoc = list(Cancer.columns[Cancer.loc[dimension] <= 0.05])
                
            if not list_Assoc:
                Cancer_dic[dimension] = []
            else:
                Cancer_dic[dimension] = list_Assoc    
        
        # Save the dictionaries:
        
        control_file = open(f'{Axes_path}Associations_{tissue}_Control.json', 'w')
        json.dump(Control_dic, control_file)
        
        cancer_file = open(f'{Axes_path}Associations_{tissue}_Cancer.json', 'w')
        json.dump(Cancer_dic, cancer_file)
          

def Get_Axes_Dynamics(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dim, matrix = "PPMI"):
    
    # Paths:
    
    network_path = "/media/sergio/sershiosdisk/Cancer/Networks/"
    axes_path    = "/media/sergio/sershiosdisk/Cancer/Axes/"
    common_path = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    
    # For each cancer:
    
    for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
        
        # Load GO embeddings:
    
        Control = pd.read_csv(f'{network_path}_GO_Embeddings_Leaf_PPI_{tissue}_{cell}_{dim}_{matrix}_Control.csv', index_col=0)
        Cancer  = pd.read_csv(f'{network_path}_GO_Embeddings_Leaf_PPI_{cancer}_{dim}_{matrix}_Cancer.csv', index_col=0)
        
        # Load the common set:
           
        common = list(pd.read_csv(f'{common_path}Common_Set_{tissue}_Leaf.csv')["0"])
        
        Control = Control[common]
        Cancer  = Cancer[common]
        
        # Map the axis from one to the other:

        Control_cano = Control.T
        Cancer_cano  = Cancer.T
                       
        correlat    = np.corrcoef(Control_cano, Cancer_cano, rowvar=False)
        corr_matrix = correlat[dim:, :dim]
        
        corr_matrix = pd.DataFrame(corr_matrix)
        
        # Save the correlation matrix:
        
        corr_matrix.to_csv(f'{axes_path}Correlation_Matrix_{tissue}.csv', header = True, index=True) 
        
        # Do the mapping:
        
        result_list = []
        
        # From Control to Cancer:
        
        for axes in list(Control.index):
            
            correlation     = abs(corr_matrix[axes])
            correlation_idx = correlation.idxmax()
            
            result_list.append([axes, correlation_idx, max(correlation)])
        
        Result_db_Control = pd.DataFrame(result_list)    
        Result_db_Control.columns = ["Control", "Cancer", "Score"]
        
        # Merge both:
        
        Result_db_Control = Result_db_Control.drop_duplicates()
        Result_db_Control = Result_db_Control.reset_index(drop = True)
        
        # Save the infor:
        
        Result_db_Control.to_csv(f'{axes_path}Mapping_Matrix_{tissue}.csv', header = True, index=True)
       
def Plot_Axes_Dynamics(Normal_Tissue_list, dim, head = False):  
    
    # Path:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    
    for tissue in Normal_Tissue_list:
        
        # Load the information
        
        maping = pd.read_csv(f'{axes_path}Mapping_Matrix_{tissue}.csv', index_col=0)
        maping = maping.sort_values("Control", ascending = True)
        maping = maping.reset_index(drop = True)
        
        if head != False:
            maping = maping.head(head)
            
        # Plot:
        
        fig= plt.figure()
        sankey(left=maping["Control"], right=maping["Cancer"], leftWeight=maping["Score"],
               aspect=15, fontsize=0)
        fig = plt.gcf()
        fig.set_size_inches(15, 15)
        fig.set_facecolor("w")
        
        fig.text(0.11, 0.5, 'Control Axes', ha='center', va='center', 
                     rotation='vertical', fontsize=30, fontweight = "bold")
        
        fig.text(0.9, 0.5, 'Cancer Axes', ha='center', va='center', 
                     rotation= 269, fontsize=30, fontweight = "bold")
        
        if head == False:
            fig.savefig(f'{axes_path}Dynamics_Axes_{tissue}_{dim}.png',  format="png", dpi=600,
                        bbox_inches='tight')
        else:
            fig.savefig(f'{axes_path}Dynamics_Axes_{tissue}_{dim}_Head.png',  format="png", dpi=600,
                        bbox_inches='tight')
 

def Jaccar_index(cluster1, cluster2):
    intersection = len(list(set(cluster1).intersection(cluster2)))
    union = (len(cluster1) + len(cluster2)) - intersection
    return float(intersection) / union              

def Compute_Jaccard(Normal_Tissue_list, dim):
    
    # Paths:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    
    for tissue in Normal_Tissue_list:
        
        # Load the information:
        
        control = open(f'{axes_path}Associations_{tissue}_Control.json',)
        control = json.load(control)
        
        cancer = open(f'{axes_path}Associations_{tissue}_Cancer.json',)
        cancer = json.load(cancer)
        
        # Keep only annotated clusters;
        
        control = dict((k,v) for k,v in control.items() if v)
        cancer  = dict((k,v) for k,v in cancer.items() if v)
        
        # DataFrame Results:
        
        Results_Jaccard = pd.DataFrame(np.nan, index= np.arange(len(control.keys())) , 
                                   columns=np.arange(len(cancer.keys()))) 
        
        Results_Jaccard.index   = list(control.keys())
        Results_Jaccard.columns = list(cancer.keys())
        
        # Calculate the JI:

        for i in list(control.keys()):
            list_control = control[i]
            for e in list(cancer.keys()):
                list_cancer = cancer[e]
                jaccard = Jaccar_index(list_control, list_cancer)
                Results_Jaccard.loc[i,e] = jaccard
        
        # Save Jaccard Index:
        
        Results_Jaccard.to_csv(f'{axes_path}Jaccard_Index_{tissue}_{dim}.csv', header = True, index=True)
        
def Plot_Jaccard(Normal_Tissue_list, dim):
    
    # Paths:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    
    # Set plots:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')    
    
    fig, axes = plt.subplots(2, 3, figsize=(10,10), 
                             gridspec_kw={"width_ratios": [1,1,0.05]})
    
    fig.tight_layout()
    
    # Control labels:
    
    plot_labels = []
    for i in range(len(Normal_Tissue_list)):
        plot_labels.append(ascii_lowercase[i])
        
    # Control variables:
    
    control_rows    = 0
    control_columns = 0
    control_labels  = 0
    
    # Start ploting:
    
    for tissue in Normal_Tissue_list:
        
        jaccard = pd.read_csv(f'{axes_path}Jaccard_Index_{tissue}_{dim}.csv', index_col = 0)
        im = axes[control_rows,control_columns].matshow(jaccard.values, interpolation='nearest', cmap="seismic", vmin=0, vmax=1)
        axes[control_rows,control_columns].set_title(f'{plot_labels[control_labels].capitalize()}', fontweight='bold', fontsize = 17, y =1)
        axes[control_rows,control_columns].axis('off')
        
        control_columns +=1
        control_labels  +=1
        
        if control_columns >1:
            control_columns = 0
            control_rows+=1
    
    # Finsh the plot:
    
    axes[0,2].axis("off")
    
    # Color bars:

    fig.colorbar(im, cax= axes[1,2], 
                         ax=[axes[0,0],axes[0,1],axes[1,0],
                             axes[1,1]])
    
    fig.text(0., 0.5, 'Control Axes', ha='center', va='center', 
                     rotation='vertical', fontsize=30, fontweight = "bold")
    
    fig.text(0.5, 0., 'Cancer Axes', ha='center', va='center', 
                     rotation='horizontal', fontsize=30, fontweight = "bold")
    
    fig.savefig(f'{axes_path}Jaccard_Index_{dim}.png',  format="png", dpi=600,
                        bbox_inches='tight')
            
           
def Compute_Semantic_Sim(Normal_Tissue_list, dim):  
    
    # Paths:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    
    for tissue in Normal_Tissue_list:
        
        # Load the information:
        
        control = open(f'{axes_path}Associations_{tissue}_Control.json',)
        control = json.load(control)
        
        cancer = open(f'{axes_path}Associations_{tissue}_Cancer.json',)
        cancer = json.load(cancer)
        
        # Keep only annotated clusters;
        
        control = dict((k,v) for k,v in control.items() if v)
        cancer  = dict((k,v) for k,v in cancer.items() if v)
        
        # DataFrame Results:
        
        Results_ss = pd.DataFrame(0, index= np.arange(len(control.keys())) , 
                                   columns=np.arange(len(cancer.keys()))) 
        
        Results_ss.index   = list(control.keys())
        Results_ss.columns = list(cancer.keys())
        
        # Get info for semantic:
        
        G = graph.from_resource("go-basic")
        similarity.precalc_lower_bounds(G)
        
        for i in list(control.keys()):
            list_control = control[i]
            for j in list(cancer.keys()):
                cluster_mean = []
                list_cancer = cancer[j]
                for GO1 in list_control:
                    for GO2 in list_cancer:
                        try:

                            semantic = similarity.lin(G, GO1, GO2)
                            cluster_mean.append(semantic)
                                                      
                        except Exception as PGSSLookupError:
                            continue
                
                Results_ss.loc[i,j] = np.mean(cluster_mean)
                
        
        # Save Jaccard Index:
        
        Results_ss.to_csv(f'{axes_path}Semantic_Index_{tissue}_{dim}.csv', header = True, index=True)            
            

def Correlation_Mapping_Jaccard(Normal_Tissue_list, dim):

    # Paths:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    
    final = pd.DataFrame()
    
    for tissue in Normal_Tissue_list:
        
        # Load Mapping and Jacard:
        
        mapping = pd.read_csv(f'{axes_path}Mapping_Matrix_{tissue}.csv', index_col=0)
        jaccard = pd.read_csv(f'{axes_path}Jaccard_Index_{tissue}_{dim}.csv', index_col = 0)
        
        # Combine the information:
        
        jc_list = []
        
        for row in range(len(mapping)):
            try:
                j_it = jaccard.loc[mapping.loc[row,"Control"], str(mapping.loc[row,"Cancer"])]
                jc_list.append(j_it)
            except KeyError:
                jc_list.append(None)
                
        mapping["Jaccard"] = jc_list
        
        final = final.append(mapping)
    
    # Calculate correlations:
    
    final = final.sort_values("Score", ascending = False)
    final = final.reset_index(drop = True)
    final = final.dropna()
    
    corr = pearsonr(final["Score"], final["Jaccard"])
    
    xy = np.vstack([list(final["Score"].values),list(final["Jaccard"].values)])
    z  = gaussian_kde(xy)(xy)
    
    # Plot the correlation:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold') 
    
    fig, ax = plt.subplots(1,1, figsize=(6,6), facecolor='w')
    
    ax.scatter(list(final["Score"]), list(final["Jaccard"]), c=z, s=100, edgecolor=None)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    
    ax.spines['left'].set_color("grey")
    ax.spines['bottom'].set_color("grey")
    
    ax.set_xlabel("Mapping Score", fontweight='bold', fontsize= 20)
    ax.set_ylabel("Jaccard index", fontweight='bold', fontsize= 20)   

    ax.text(0.3, 0.9, f'Corr: {round(corr[0],1)}*', ha = "center", size = "medium") 

    # Save the plot:

    plt.savefig(f'{axes_path}Global_Correlation_Jaccard_Pearson.png',  format="png", dpi=600,
                        bbox_inches='tight') 


def Get_GO_lists_Axes(mapping, control_dic, cancer_dic):
    
    # Start the lists:

    Control = []
    Cancer  = []
    
    for row in list(mapping.index):
        
        GO_Control = set(control_dic[str(int(mapping.loc[row]["Control"]))])
        GO_Cancer  = set(cancer_dic[str(int(mapping.loc[row]["Cancer"]))])
        
        Control.append(list(GO_Control))
        Cancer.append(list(GO_Cancer))
        
    # Return the new mapping with the information:
    
    mapping["Control_GO"] = Control
    mapping["Cancer_GO"]  = Cancer
    
    return(mapping)
        

def Get_Complete_Mapping_Info(Normal_Tissue_list, dim):
    
    '''
    This function adds to the mapping produced in Get_Axes_Dynamics, the information of Jaccard Index and the GO terms that are
    associated to each axis.
    
        - Normal_Tissue_list : string, tissue of the sample.
        - dim                : int, number of dimensins
    ''' 
    
    # Paths:
    
    axes_path = "/media/sergio/sershiosdisk/Cancer/Axes/"
    save_path = "/media/sergio/sershiosdisk/Cancer/Axes_Domains/"
    
    # Compute:
    
    for tissue in Normal_Tissue_list:
        
        # Load the mapping between Control and Cancer:
        
        mapping = pd.read_csv(f'{axes_path}Mapping_Matrix_{tissue}.csv', index_col=0)
        
        # Load Associations:
        
        control_dic = open(f'{axes_path}Associations_{tissue}_Control.json',)
        control_dic = json.load(control_dic) 
        
        cancer_dic = open(f'{axes_path}Associations_{tissue}_Cancer.json',)
        cancer_dic = json.load(cancer_dic)
        
        # Load Jaccard Index:
        
        jaccard = pd.read_csv(f'{axes_path}Jaccard_Index_{tissue}_{dim}.csv', index_col = 0)
        
        # Pre-Process the data:

        # Delete non-annotated axes:
        
        control_dic = dict((k,v) for k,v in control_dic.items() if v)
        cancer_dic  = dict((k,v) for k,v in cancer_dic.items() if v)
        
        # Merge information:
            
        jaccard_pairs = []
            
        for row in range(len(mapping)):
            try:
                j_it = jaccard.loc[mapping.loc[row,"Control"], str(mapping.loc[row,"Cancer"])]
                jaccard_pairs.append(j_it)
            except KeyError:
                jaccard_pairs.append(None) 
        
        # Delete non annotated axes from the mapping:
        
        mapping["Jaccard"] = jaccard_pairs 
            
        mapping = mapping.sort_values("Score", ascending = False)
        mapping = mapping.reset_index(drop = True)
        mapping = mapping.dropna()
        
        # Add the GO information to the mapping:
        
        mapping = Get_GO_lists_Axes(mapping, control_dic, cancer_dic)
        
        # Save the complete information:
        
        mapping.to_csv(f'{save_path}Complete_Mapping_{tissue}.csv', header = True, index=True)  

def Send_REVIGO(information, save_path, name, number, tissue):
    
    '''
    This sends the list of GO terms to REVIGO. It produces the functional domains. This function is used in Produce_Axes_Functional_Domains. 
    
        - information : DataFrame, subset of the original mapping using the thresholds.
        - save_path   : string, path to save the output.
        - name        : string, name for the file.
        - number      : int, threshold used to keep the information.
        - tissue      : string.         

    ''' 
    
    # Preocessing
    
    info_cancer  = information["Cancer_GO"].apply(literal_eval)
    info_control = information["Control_GO"].apply(literal_eval)
    
    info_cancer  = list(set([item for sublist in info_cancer for item in sublist]))
    info_control = list(set([item for sublist in info_control for item in sublist]))
    
    # Get the sets:
    
    info_inter      = list(set(info_cancer).intersection(set(info_control)))
    info_cancer_f   = list(set(info_cancer) - set(info_control))
    info_control_f  = list(set(info_control) - set(info_cancer))
    
    # Send Revigo:
    
    save_path_revigo_inter   = f'{save_path}{name}_Domains_Inter_{tissue}_{number}'
    save_path_revigo_cancer  = f'{save_path}{name}_Domains_Cancer_{tissue}_{number}'
    save_path_revigo_control = f'{save_path}{name}_Domains_Control_{tissue}_{number}'
    
    Submit_Query(info_inter,   save_path_revigo_inter)
    Submit_Query(info_cancer_f,  save_path_revigo_cancer)
    Submit_Query(info_control_f, save_path_revigo_control)
    
    
def Produce_Axes_Functional_Domains(Normal_Tissue_list, dim, top = 0.8, tail = 0.5):
    
    '''
    This functions produces the functional domains of the pairs of axes. Per each group of pair axes (control to cancer), it computes the 
    domains (summary) of the GO terms that are only annotated in Control, Cancer, or in both.
    
        - Normal_Tissue_list : string list, names of the tissues (e.g., breast, colon, lung, prostate).
        - dim                : int, dimensionality.
        - top                : int, threshold to consider non-moving axes (0.8 by default).
        - tail               : int, threshold to consider moving axes (0.5 by default).

    ''' 
    # Paths:
    
    save_path = "/media/sergio/sershiosdisk/Cancer/Axes_Domains/" 
    
    for tissue in Normal_Tissue_list:
        
        # Load the mapping:
        
        mapping = pd.read_csv(f'{save_path}Complete_Mapping_{tissue}.csv', index_col = 0)
        
        # Get the corresponding threshold:
        
        mapping_tail = mapping[mapping.Score <= tail]
        mapping_top  = mapping[mapping.Score >= top]
        mapping_mid  = mapping[(mapping.Score > tail) & (mapping.Score < top)]
        
        # Do the query:
        
        Send_REVIGO(mapping_tail, save_path, "TAIL",  str(tail), tissue)
        Send_REVIGO(mapping_top,  save_path, "TOP",   str(top), tissue)
        Send_REVIGO(mapping_mid,  save_path, "MIDLE", f'{top}_to_{tail}', tissue)


def Give_the_Domains(path, distribution, tissue, set_dist, save_path):
    
    
    domain_path = f'{path}Revigo.csv' 
       
    # Open file:
       
    file = pd.read_csv(domain_path, sep = " ")
    term_ids = []
    
    for line in range(len(file)):
        file_line = file.loc[line]
        if file_line["Eliminated"] == True:
            temp_representative = file_line["Representative,"]
            term_ids.append(temp_representative)
        else:
            term_ids.append(file_line["TermID,"])
        
    # Count the GO terms
        
    word = "GO"
    real_terms = []
    for line in range(len(term_ids)):
        Term = term_ids[line]
        if word in Term:
            real_terms.append(Term)
        else:
            new_term = real_terms[len(real_terms) - 1]
            real_terms.append(new_term)
            
    # Count them:
        
    pd.DataFrame(real_terms).count(0)
        
    Counter_dic = Counter(real_terms)
    Count_pd    = pd.DataFrame.from_dict(Counter_dic, orient='index')
    Count_pd    = Count_pd.sort_values(0, ascending = False)
    

    # Final result:
        
    GO_filter     = list(Count_pd.index)
    file.index    = file["TermID,"]
    GO_decription = list(file.loc[GO_filter]["Name,"])
    Counter_list  = Count_pd[0]
        
    Final_pandas = pd.DataFrame({"GO": GO_filter, "Description" : GO_decription, 
                                     "Counter" : Counter_list }).reset_index(drop = True)
            
    Final_pandas.to_csv(f'{save_path}Domains_{tissue}_{distribution}_{set_dist}.csv') 
     
        
def Produce_Mapped_Axes_Domains(Normal_Tissue_list, top = 0.8, tail = 0.5, distribution = "TAIL"):   
    
    # Paths:
    
    domain_path = "/media/sergio/sershiosdisk/Cancer/Axes_Domains/"
    
    for tissue in Normal_Tissue_list:
        
        # Choose the informatio to load:
        
        if distribution == "TAIL":
            
            path_control = f'{domain_path}TAIL_Domains_Control_{tissue}_{tail}/'
            path_cancer  = f'{domain_path}TAIL_Domains_Cancer_{tissue}_{tail}/'
            path_midle   = f'{domain_path}TAIL_Domains_Inter_{tissue}_{tail}/'
            
        elif distribution == "TOP":
            
            path_control = f'{domain_path}TOP_Domains_Control_{tissue}_{top}/'
            path_cancer  = f'{domain_path}TOP_Domains_Cancer_{tissue}_{top}/'
            path_midle   = f'{domain_path}TOP_Domains_Inter_{tissue}_{top}/'
        
        elif distribution == "MIDLE":
            
            path_control = f'{domain_path}MIDLE_Domains_Control_{tissue}_{top}_to_{tail}/'
            path_cancer  = f'{domain_path}MIDLE_Domains_Cancer_{tissue}_{top}_to_{tail}/'
            path_midle   = f'{domain_path}MIDLE_Domains_Inter_{tissue}_{top}_to_{tail}/'
            
        # Compute the domains:
        
        Give_the_Domains(path_control, distribution, tissue, "Control", domain_path)
        Give_the_Domains(path_cancer,  distribution, tissue, "Cancer",  domain_path)
        Give_the_Domains(path_midle,   distribution, tissue, "Inter",   domain_path)
        
def Plot_Mapped_Axes_Domains(Normal_Tissue_list, n_domains, distribution = "TAIL"):
    
    # Paths:
    
    domain_path = "/media/sergio/sershiosdisk/Cancer/Axes_Domains/"    
    
    for tissue in Normal_Tissue_list:
        
        # Choose the informatio to load:
        
        if distribution == "TAIL":
            
            path_control = f'{domain_path}Domains_{tissue}_{distribution}_Control.csv'
            path_cancer  = f'{domain_path}Domains_{tissue}_{distribution}_Cancer.csv'
            path_midle   = f'{domain_path}Domains_{tissue}_{distribution}_Inter.csv'
            
        elif distribution == "TOP":
            
            path_control = f'{domain_path}Domains_{tissue}_{distribution}_Control.csv'
            path_cancer  = f'{domain_path}Domains_{tissue}_{distribution}_Cancer.csv'
            path_midle   = f'{domain_path}Domains_{tissue}_{distribution}_Inter.csv'
        
        elif distribution == "MIDLE":
            
            path_control = f'{domain_path}Domains_{tissue}_{distribution}_Control.csv'
            path_cancer  = f'{domain_path}Domains_{tissue}_{distribution}_Cancer.csv'
            path_midle   = f'{domain_path}Domains_{tissue}_{distribution}_Inter.csv'  
        
        # Start the plots:
    
        # Authomatic Labels:
    
        plot_labels = []
    
        for i in range(3):
            plot_labels.append(ascii_lowercase[i])  
        
        # Plot features:
    
        plt.style.use("seaborn-whitegrid")
        plt.rcParams.update({'font.size': 17})
        plt.rc('font', weight='bold') 
        
        fig, axes = plt.subplots(3, 1, figsize=(18,18))
        fig.tight_layout()
        
        # Read the Domains:
        
        domain_Control = pd.read_csv(path_control, index_col = 0)
        domain_Cancer  = pd.read_csv(path_cancer,  index_col = 0)
        domain_Inter   = pd.read_csv(path_midle,   index_col = 0)
        
        domain_Control_rep  = domain_Control.iloc[:n_domains]
        domain_Control_rep  = domain_Control_rep.replace(',','', regex=True)
    
        domain_Cancer_rep  = domain_Cancer.iloc[:n_domains]
        domain_Cancer_rep  = domain_Cancer_rep.replace(',','', regex=True)   
    
        domain_Inter_rep  = domain_Inter.iloc[:n_domains]
        domain_Inter_rep  = domain_Inter_rep.replace(',','', regex=True)  
        
        domain_tot       = [domain_Control, domain_Cancer, domain_Inter]
        domain_Final_rep = [domain_Control_rep, domain_Cancer_rep, domain_Inter_rep]
        
        # Iterate over the domains:
        
        row_control    = 0
        label_count    = 0
        
        for i in range(len(domain_Final_rep)):
            
            # Information for the PiePlot:
        
            Total           = sum(domain_Final_rep[i].Counter) 
            percentages     = [i*100/Total for i in domain_Final_rep[i].Counter]
            Real_Total      = sum(domain_tot[i].Counter)
            Real_Percentage = (Total*100)/ Real_Total # The coverag
            
            # Colors for the pie plot:
            
            cmap        = plt.get_cmap('tab20')
            colors      = [cmap(i) for i in np.linspace(0, 1, n_domains)]
            final_color = colors            
            
           # Plot:
            
            patches, texts, autotexts = axes[row_control].pie(percentages, autopct='%1.f%%', shadow=False, colors= final_color)
            axes[row_control].set_title(plot_labels[label_count].capitalize(), fontweight='bold', fontsize = 20, y =0.5, x = -0.01)
            circle = plt.Circle(xy=(0,0), radius=0.75, facecolor='white')
            axes[row_control].add_patch(circle)
            
            # Change the text size:

            for perc in range(len(autotexts)):
                autotexts[perc].set_fontsize(17)
                autotexts[perc].set_fontweight("bold")
            axes[row_control].text(-0.5,1.15,f'Coverage : {round(Real_Percentage,1)}%', fontsize = 20)
            
            # Legend:
    
            list_description = []
            for description in range(n_domains):
                list_description.append(domain_Final_rep[i].iloc[description]["Description"])
                
            Order_Stuff = pd.DataFrame({"Value" : list_description, "Colors" : final_color})
        
        
            # Control the length of the descriptions:
            
            jump = "\n"
            
            for index in Order_Stuff["Value"].index:
                if len(Order_Stuff["Value"][index]) > 90:
                    it_setence   = Order_Stuff["Value"][index].split()
                    mid_pos      = len(it_setence)//2
                    new_sentence = it_setence[:mid_pos] + [jump] + it_setence[mid_pos:]
                    new_sentence = ' '.join(new_sentence)
                    
                    Order_Stuff["Value"][index] = new_sentence
        
            # Create a legend:
            
            patch_list = []
            
            for patch in range(n_domains):
                
               patch_list.append(mpatches.Patch(color= final_color[patch], label= Order_Stuff["Value"][patch]))
        
            axes[row_control].legend(handles=patch_list,borderaxespad=0.1,
                       bbox_to_anchor=(0.95, 0.85), loc="upper left", frameon=True)
            
            row_control+=1
            label_count+=1   
        
        # Save the plot:
        
        fig.savefig(f'{domain_path}Functional_Domains_{tissue}_{distribution}_{n_domains}.png', format="png", dpi=600,
                    bbox_inches='tight') 


def Plot_Cancer_Axes_Comparisons(Normal_Tissue_list): 

    # Paths:
        
    save_path = "/media/sergio/sershiosdisk/Cancer/Axes_Domains/"
    
    # Prepare the DataFrame to plot:
    
    mapping = pd.DataFrame()
    
    for tissue in Normal_Tissue_list:
            mapping_it = pd.read_csv(f'{save_path}Complete_Mapping_{tissue}.csv', index_col = 0)
            mapping_it = mapping_it[mapping_it.Score <= 0.6]
            mapping_it.Cancer  = [f'{i}_{tissue}' for i in mapping_it.Cancer]        
            mapping_it = mapping_it[[ "Cancer", "Cancer_GO"]]
            mapping = mapping.append(mapping_it)    
    
    mapping.index = mapping.Cancer
    mapping       = mapping.drop_duplicates()
    mapping["Cancer_GO"] = mapping["Cancer_GO"].apply(literal_eval)    
    
    # Get the JI:
            
    JI_final = pd.DataFrame(np.nan, index= np.arange(len(mapping)) , 
                                       columns=np.arange(len(mapping)))
    
    JI_final.index   = mapping.Cancer
    JI_final.columns = mapping.Cancer    
    
    for i in list(mapping.index):
        axes_1 = mapping.loc[i]["Cancer_GO"]
        for j in list(mapping.index):
            axes_2 = mapping.loc[j]["Cancer_GO"]
            JI_final.loc[i][j] = Jaccar_index(axes_1, axes_2)  
    JI_final = JI_final.fillna(0)  
    
    # Use cluster map to get the order:
    
    linkage       = hc.linkage(sp.distance.squareform(1 - JI_final),method='average')
    Control_clust = sns.clustermap(1- JI_final, row_linkage=linkage, col_linkage=linkage)
    order         = Control_clust.dendrogram_row.reordered_ind
    
    JI_final = JI_final.iloc[order]
    JI_final = JI_final[JI_final.index]
    
    # Now plot using plt, easier to play with:
    
    # plot the main with the zoomed:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 17})
    plt.rc('font', weight='bold')
    fig, ax = plt.subplots(figsize=[10,10])
    fig.tight_layout()
    
    
    im = ax.matshow(JI_final, cmap="seismic")
    axins = zoomed_inset_axes(ax, 5, loc=1, bbox_to_anchor=(950,600))
    axins.matshow(JI_final, cmap="seismic")
    

    # Ticks information:
    
    ax.set_yticks(range(len(JI_final.columns)))
    ax.set_xticks(range(len(JI_final.columns)))  
    axins.set_yticks(range(len(JI_final.columns)))
    axins.set_xticks([], [])
    
    # labels name:
    
    names_ticks = list(JI_final.columns)
    ax.set_yticklabels(names_ticks)
    axins.set_yticklabels(names_ticks)
    axins.yaxis.tick_right()
    
    ax.axis('off')
    
    # Add the ticks for the zoomed:
    
    x1, x2, y1, y2 =   42, 57, 57, 42
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    
    plt.xticks(visible=False)
    plt.yticks(visible=True)
    
    _patch, pp1, pp2 = mark_inset(ax, axins, loc1=4, loc2=4, fc="none", ec="grey", linewidth= 3)
    
    pp1.loc1, pp1.loc2 = 2, 3  
    pp2.loc1, pp2.loc2 = 3, 2
    
    # Add the color bar:
    
    clb = plt.colorbar(im, ax=ax, shrink=0.5, orientation="horizontal")
    clb.set_label('Jaccard Index', labelpad=-60, y=3200, rotation=0, size=20,weight='bold')
    
    plt.draw()
    plt.show()
    

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
    This functions calculates the cosine distance between the gene/GO common space
    '''
    
    # Paths:
    
    network_path = "./Cleaned_version/Data/Embeddings/"
    gene_GO_path = "./Cleaned_version/Data/"
    moving_path  = "./Cleaned_version/Data/"
    FMM_path     = "./Cleaned_version/Data/FMM/"
    
    # For each combination:
    
    for cancer, tissue, cell in zip(Cancer_list,Normal_tissue_list, Cell_type_list):
        
        # Prepare the paths:
        
        # Genes:
        
        genes_Control =  f'{gene_GO_path}Control_{tissue}_{cell}_Genes.csv'
        genes_Cancer  =  f'{gene_GO_path}Cancer_{cancer}_Genes.csv'
        
        # Gene embedded coordinates:
        
        Gene_embedd_Control = f'{network_path}_P_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        Gene_embedd_Cancer  =  f'{network_path}_P_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy'
        
        # Middle factor:
        
        midle_Control = f'{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        midle_Cancer  = f'{network_path}_U_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy'
        
        Gene_midle_Control_db  = pd.DataFrame(np.load(midle_Control, allow_pickle=True))
        Gene_midle_Cancer_db   = pd.DataFrame(np.load(midle_Cancer,  allow_pickle=True))
        
        # GO embedded coordinates:
        
        GO_embedd_Control = f'{FMM_path}_GO_Embeddings_GO_PPI_{tissue}_{cell}_{dim}_PPMI_Control.csv'
        GO_embedd_Cancer  = f'{FMM_path}_GO_Embeddings_GO_PPI_{cancer}_{dim}_PPMI_Cancer.csv'
        
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
        
        Matrix_GO_Gene = pd.read_csv(f'./Cleaned_version/Data/_Matrix_Genes_GO_BP_PPI.csv', 
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
    
    # Paths:

    network_path = "./Cleaned_version/Data/"
    gene_GO_path = "./Cleaned_version/Data/"
    matrix_path  = "./Cleaned_version/Data/"
    
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

def genes_Color(Cancer_Cosine_it, Control_Cosine_it, GO1):

    gene_common = Cancer_Cosine_it.index[Cancer_Cosine_it.index.isin(Control_Cosine_it.index)]
    
    Cancer_Cosine_it_com  = Cancer_Cosine_it.loc[gene_common]
    Control_Cosine_it_com = Control_Cosine_it.loc[gene_common] 

    Cancer_structure_GO1  = Cancer_Cosine_it_com[GO1]
    Control_structure_GO1 = Control_Cosine_it_com[GO1]
    
    move = Control_structure_GO1 - Cancer_structure_GO1
    
    set_close  = []
    set_far    = []
    set_stable = []
    x1= []
    x2= []
    
    for gene in move.index:
        if move[gene] > (np.mean(move) + 2*np.std(move)):
            set_far.append(gene)
            x1.append(move[gene])
        elif move[gene] < (np.mean(move) - 2*np.std(move)):
            set_close.append(gene)
            x2.append(move[gene])
        else:
            set_stable.append(gene)
    
    final = []
    final.extend(set_close)
    final.extend(set_far)
    
    x_final = []
    x_final.extend(x1)
    x_final.extend(x2)
    
    max_value = max(x_final)
    max_index = x_final.index(max_value)
    
    max_value = min(x_final)
    max_index = x_final.index(max_value)
    
    return(final)        

def Plot_Distribution_movement(Normal_tissue_list):
       
    # Paths:
    
    path_rank = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # Iterate by each cancer:
    
    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({'font.size': 20})
    plt.rc('font', weight='bold')
                               
    fig, axes = plt.subplots(4, figsize=(10,10))
    fig.tight_layout()
    
    row_control = 0
            
    for tissue in Normal_tissue_list:
        
        # Load distribution:
        
        rank = pd.read_csv(f'{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv')
        
        threshold = np.mean(rank["0"])
        
        axes[row_control].hist(rank["0"], bins = 20)
        axes[row_control].set_title(f'{tissue}')
        
        row_control+=1
    
    fig.savefig(f'{path_rank}Movement_GO_Distribution.png', format="png", dpi=300,
                    bbox_inches='tight')
        
def PokePlot(Normal_tissue_list):
    
    # Paths:
    
    gene_GO_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    ranking_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # For each tissue:
    
    for tissue in Normal_tissue_list:
        
        # Load the data:
        
         ranking_p    = f'{ranking_path}top_100_moving_{tissue}.csv'
         Ranking_db = pd.read_csv(ranking_p,  index_col = 0)           
         top_10 = Ranking_db.index[0:10]
         
         Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
         Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
            
         Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
         Control_Cosine_it.index  = [str(i) for i in Control_Cosine_it.index] 
         
         # Controllers:
        
         row_control    = 0
         column_control = 0
         label_control  = 0
          
        # Prepare the labels for each sub-plot:
            
         plot_labels = []
         for i in range((10)):
             plot_labels.append(ascii_lowercase[i])
                 
        # Start the plot:
                    
         plt.style.use("seaborn-whitegrid")
         plt.rcParams.update({'font.size': 20})
         plt.rc('font', weight='bold')
                               
         fig, axes = plt.subplots(4, 3, figsize=(18,18))
         fig.tight_layout()
        
         for GO1 in top_10:
            
            Cancer_structure_GO1  = list(Cancer_Cosine_it[GO1])
            Control_structure_GO1 = list(Control_Cosine_it[GO1]) 
            
            gene_Cancer  = list(Cancer_Cosine_it[GO1].index)
            Gene_Control = list(Control_Cosine_it[GO1].index)
            
            gene  = []
            final = []
            
            final.extend(Cancer_structure_GO1)
            final.extend(Control_structure_GO1)
            gene.extend(gene_Cancer)
            gene.extend(Gene_Control)
            
            #  List of top movement genes:
            
            list_movement = genes_Color(Cancer_Cosine_it, Control_Cosine_it, GO1)
            
            # Set the colors:
               
            colors = []

            n_dots = len(final)   # set number of GO terms
            angs = np.linspace(0, 2*np.pi, n_dots)  # angles of the GO terms to keep the circle.
            cx, cy = (0, 0)  # center of circle
            xs, ys = [], []    # for coordinates of points to plot
            ra = 1         # radius of circle
            
            count = 0
            for ang in angs:
                # radius:
                ra_it = ra - (1-final[count])
                
                # Cancer distances:
                
                if count < len(Cancer_Cosine_it):
                    if gene[count] in list_movement:
                        colors.append("#2e2b04")
                    else:
                        colors.append("#2679ea") 
                else:
                    if gene[count] in list_movement:
                        colors.append("#2e2b04")
                    else:
                        colors.append("#bd3751")
                                      
                x = cx + ra_it*np.cos(ang)
                y = cy + ra_it*np.sin(ang)
                xs.append(x)   # collect x
                ys.append(y)   # collect y
                count+=1                                      
                                      
            axes[row_control, column_control].scatter(xs, ys, s=9 ,c=colors, alpha = 0.6)  # plot points 
            circle = plt.Circle([0,0], 0.02, color = "#04fbe1")
            axes[row_control, column_control].add_patch(circle)
            axes[row_control, column_control].spines['right'].set_visible(False)
            axes[row_control, column_control].spines['top'].set_visible(False)
            axes[row_control, column_control].spines['left'].set_linewidth(4)
            axes[row_control, column_control].spines['bottom'].set_linewidth(4) 
            axes[row_control, column_control].spines['left'].set_color("grey")
            axes[row_control, column_control].spines['bottom'].set_color("grey")
            axes[row_control, column_control].set_yticklabels(" ")
            axes[row_control, column_control].set_xticklabels(" ") 
            axes[row_control, column_control].set_xticks([])
            axes[row_control, column_control].set_yticks([])
            axes[row_control, column_control].set_title(f'{GO1}', fontweight='bold', 
            fontsize = 17, y = 1, x = 0.5)
            
            # Add the closest gene in each space:
            
            cancer_min  = Cancer_Cosine_it[GO1].sort_values(ascending = True).index[0]
            control_min = Control_Cosine_it[GO1].sort_values(ascending = True).index[0]
            
            axes[row_control, column_control].text(1.2, 0.5, f'{cancer_min}',  horizontalalignment='center', verticalalignment='center', color = "#2679ea")
            axes[row_control, column_control].text(1.2, 0.3, f'{control_min}', horizontalalignment='center', verticalalignment='center', color = "#bd3751")   
                
            if column_control == 2:
                column_control = 0
                row_control +=1
                label_control +=1
            else:
                column_control +=1
                label_control +=1
                    
            # For the empty sets:
                    
            if (column_control == 0) & (row_control == 3):
                axes[row_control, column_control].axis('off')
                column_control +=1
                        
            if (column_control == 2) & (row_control == 3):
                axes[row_control, column_control].axis('off')
                column_control +=1
            
         Cancer_path     = mpatches.Patch(color='#2679ea', label='Cancer')
         Control_patch   = mpatches.Patch(color='#bd3751', label='Control')
         Move_patch      = mpatches.Patch(color='#2e2b04', label='Top Change')
                                         
         plt.legend(handles=[Cancer_path, Control_patch, Move_patch],borderaxespad=0.1,
                                    bbox_to_anchor=(-0.69, -0.3), loc="lower center", frameon=True, ncol = 3)
            
         # Save the plot:
            
         plt.savefig(f'{gene_GO_path}Common_sace_{tissue}.png', format="png", dpi=600,
                    bbox_inches='tight')
                                
def extract_and(*args): 
    #modificato per kirc %5Bmh%5D (altri tumori mettere %5Btw%5D)!!!
    query = "https://pubmed.ncbi.nlm.nih.gov/?term=%28{}%5Btw%5D%29".format(args[0].replace(" ", "+"))
    for a in args[1:]:
        query += "+AND+%28{}%5Btw%5D%29".format(a.replace(" ", "+"))
    contains = ''
    while contains == '':
        try:
            contains = requests.get(query)
            break
        except:
            time.sleep(5)
            continue
    m = re.search('<meta name="log_resultcount" content="[0-9]+"', str(contains.content))
    if m is not None:
        return int(m.group(0)[38:-1])
    if m is None:
        m1 = re.search('<meta name="log_source_db" content=".+" />',str(contains.content))
        if m1 is not None:
            return 1
        if m1 is None:
            return 0        
        
def Compute_common_genes(Normal_tissue_list):
    
    # Paths:
    
    gene_GO_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    
    # For each tissue:
    
    for tissue in Normal_tissue_list: 
        
        # Load the distances:
        
         Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
         Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
            
         Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
         Control_Cosine_it.index  = [str(i) for i in Control_Cosine_it.index] 
         
         # Take the common set of genes:
        
         gene_common = list(Cancer_Cosine_it.index[Cancer_Cosine_it.index.isin(Control_Cosine_it.index)])
         gene_common = pd.DataFrame(gene_common)
        
         # Save the common set of genes (to translate them later)
        
         gene_common.to_csv(f'{gene_GO_path}Common_Set_Genes_{tissue}.csv',   header = True, index=True)
         
def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='100000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results   
        
def Query_Common_gene_set(Normal_tissue_list):
    
    '''
    To run this function the file Transformed_Common_{tissue}.csv is needed. This file corresponds to a query in
    bioMart. For each gene, we retrieve the description of the gene and different name codes. The file is attached
    but it can be manually recomputed. By using this file, the function performs a literature search for assessing if a
    gene is related to a particular cancer type.
    '''
    
    # Paths:
    
    gene_GO_path = "./Cleaned_version/Data/"
    
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


def Predict_Oncogenes(Normal_tissue_list):
    
    # Paths:
    
    Ranking_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    gene_GO_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    
    # Iterate by tissues:
    
    for tissue in Normal_tissue_list:
            
        # Load Ranking of moving GO:
            
        rank       = pd.read_csv(f'{Ranking_path}Rank_movement_{tissue}_PPMI_Leaf.csv')      
        moving_GO  = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"] 
        
        # Load the bibliography information:
            
        counter_list = pd.read_csv(f'{gene_GO_path}Cancer_Count_{tissue}.csv')
            
        # Load cosines:
            
        Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
        Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
          
        Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
        Control_Cosine_it.index = [str(i) for i in Control_Cosine_it.index]
        
        moving_GO = list(set(moving_GO).intersection(set(Cancer_Cosine_it.columns)))
            
        # List for the predicted genes:
            
        GO_list = []
        
        # Start the experiment:
            
        for GO1 in moving_GO:
                
            # If we need the common for the yellow dots:
                
            gene_common = Cancer_Cosine_it.index[Cancer_Cosine_it.index.isin(Control_Cosine_it.index)]
            
            print(GO1)
                
            Cancer_Cosine_it_com  = Cancer_Cosine_it.loc[gene_common]
            Control_Cosine_it_com = Control_Cosine_it.loc[gene_common] 
            
            Cancer_structure_GO1  = Cancer_Cosine_it_com[GO1]
            Control_structure_GO1 = Control_Cosine_it_com[GO1]
                
            move = Control_structure_GO1 - Cancer_structure_GO1

            set_close  = []
            set_far    = []
            set_stable = []
            x1= []
            x2= []    
                
            for gene in move.index:
                if move[gene] > (np.mean(move) + 2*np.std(move)):
                    set_far.append(gene)
                    x1.append(move[gene])
                elif move[gene] < (np.mean(move) - 2*np.std(move)):
                    set_close.append(gene)
                    x2.append(move[gene])
                else:
                    set_stable.append(gene)
                        
            final = []
            final.extend(set_close)
            final.extend(set_far)
                
            compi = list(counter_list[counter_list["Counts"] >= 0]["initial_alias"])
            compi =  [str(i) for i in compi]
                
            predictions =  set(final).intersection(set(compi))
            
            
            # Going back to gene name
            
            predictions_name = list(counter_list[counter_list["initial_alias"].isin(predictions)].name)
                
            GO_list.append(predictions_name)
        
                
        # Save predictions:
        
        flat_list   = [item for sublist in GO_list for item in sublist]
                     
        Final = pd.DataFrame.from_dict(Counter(flat_list), orient='index').reset_index()
        Final = Final.reset_index(drop = True)
        Final = Final.drop_duplicates("index")
        counter_list = counter_list.drop_duplicates("name")
        Final["Publications"]  = list(counter_list[counter_list["name"].isin(Final["index"])].Counts)  
          
        Final.columns = ["Prediction", "Counts_10","Publications"]
        Final     = Final.sort_values(by = "Counts_10", ascending = False) 
        
        Final.to_csv(f'{gene_GO_path}Complete_Predictions_{tissue}.csv')
        
        # For the gene names:
        
        genes = list(Control_Cosine_it.index)
        genes = pd.DataFrame(counter_list[counter_list["initial_alias"].isin(genes)].name)
        genes = genes.drop_duplicates()
        genes = genes.reset_index(drop = True)
        
        genes.to_csv(f'{gene_GO_path}Genes_All_{tissue}.csv')

def Do_Overlapping_Table(Normal_tissue_list):
    
    # Paths:
    
    gene_GO_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    path_rank    = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    
    # Open the table:
    
    with open(f'{gene_GO_path}Overlapping_Table.txt', 'a') as the_file:
        
        # Write columns:
        
        the_file.write(f'Cancer Type\t')
        
        for i in range(1,68):
            the_file.write(f'{i}\t')
        the_file.write(f'\n')
        
        for tissue in Normal_tissue_list:
            
            the_file.write(f'{tissue}\t')
            
            # Load the information:
            
            info = pd.read_csv(f'{gene_GO_path}Complete_Predictions_{tissue}.csv', index_col = 0)            
            rank      = pd.read_csv(f'{path_rank}Rank_movement_{tissue}_PPMI_Leaf.csv')      
            moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"]
            
            for i in range(1,len(moving_GO)):
                
                info_it = info[info.Counts_10 == i]
                
                if len(info_it) == 0:
                    the_file.write(f'0 (0)\t')
                else:
                    info_perc = round(len(info_it[info_it.Publications > 0])/len(info_it),2)*100
                    the_file.write(f'{round((len(info_it)/len(info))*100,2)}-({info_perc}%)\t')
            
            the_file.write(f'\n') 


def Do_Fold_Enrichment_Complete(Normal_tissue_list):
    
    path_rank    = "./Cleaned_version/Data/"
    gene_GO_path = "./Cleaned_version/Data/"
    
    # Open the file:
    
    with open(f'{gene_GO_path}Fold_Rank_Table_Dos.txt', 'a') as the_file:
        
        # Write the name of the columns:
        
        the_file.write(f'Cancer Type\t')
        the_file.write(f'Moving Genes\t')
        the_file.write(f'Stable Genes\n')
        
        for tissue in Normal_tissue_list:
            
            # Load movement and choose the most moving:
            
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
            
            # A) Moving Vs Not moving (no matter the direction):
            
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
        
            Restuls_abs_top  = Restuls_abs.head(round(len(Restuls_abs)*5/100))
            Restuls_abs_tail = Restuls_abs.tail(round(len(Restuls_abs)*5/100))
              
            # To save the predictions:
            
            gene_names = counter_list[counter_list.initial_alias.isin(Restuls_abs_top.Gene)] 
            gene_names = gene_names.drop(["Unnamed: 0", "initial_alias"], axis=1) 
            gene_names = gene_names.drop_duplicates("name")
            gene_names = gene_names.reset_index(drop = True)
            
            gene_names.to_csv(f'{path_rank}Predictions_Rank_{tissue}.csv', header = True, index=True)
            
            # Get the total set

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
            
        
def Do_Fold_Enrichment(Normal_tissue_list):
    
    # Paths:
    
    Ranking_path = "/media/sergio/sershiosdisk/Cancer/Annotations_Movement/"
    gene_GO_path = "/media/sergio/sershiosdisk/Cancer/Genes_GO/"
    
    # Iterate by tissues:
    
    for tissue in Normal_tissue_list:
        
        with open(f'{gene_GO_path}Fold_Enrichment_{tissue}.txt', 'a') as the_file:
            
            p_values_move = []
            p_values_stb  = []
            
            # Write the columns:
            
            the_file.write(f'GO term\t')
            the_file.write(f'Fold Moving (Over-Represented)\t')
            the_file.write(f'P-value Moving\t')
            the_file.write(f'Fold Stable (Under-Represented)\t')
            the_file.write(f'P-value Stable\t')            
            the_file.write(f'Predicted\t')
            the_file.write(f'Percentage Predicted Biblio\n') 
        
            # Load Ranking of moving GO:
            
            rank      = pd.read_csv(f'{Ranking_path}Rank_movement_{tissue}_PPMI_Leaf.csv')      
            moving_GO = rank[rank["0"] > (np.mean(rank["0"]) + 2*np.std(rank["0"]))]["Unnamed: 0"]            
            
            # Load the bibliography information:
            
            counter_list = pd.read_csv(f'{gene_GO_path}Cancer_Count_{tissue}.csv')
            
            # Load cosines:
            
            Cancer_Cosine_it  = pd.read_csv(f'{gene_GO_path}Cancer_Gene_GO_dist_{tissue}.csv',  index_col = 0, dtype={0: str})
            Control_Cosine_it = pd.read_csv(f'{gene_GO_path}Control_Gene_GO_dist_{tissue}.csv', index_col = 0, dtype={0: str})
          
            Cancer_Cosine_it.index  = [str(i) for i in Cancer_Cosine_it.index] 
            Control_Cosine_it.index = [str(i) for i in Control_Cosine_it.index]

            moving_GO = list(set(moving_GO).intersection(set(Cancer_Cosine_it.columns)))
    
            # Start the experiment:
    
            for GO1 in moving_GO:
                
                # If we need the common for the yellow dots:
                
                gene_common = Cancer_Cosine_it.index[Cancer_Cosine_it.index.isin(Control_Cosine_it.index)]
            
                print(GO1)
                
                Cancer_Cosine_it_com  = Cancer_Cosine_it.loc[gene_common]
                Control_Cosine_it_com = Control_Cosine_it.loc[gene_common] 
            
                Cancer_structure_GO1  = Cancer_Cosine_it_com[GO1]
                Control_structure_GO1 = Control_Cosine_it_com[GO1]
                
                move = Control_structure_GO1 - Cancer_structure_GO1
                
                set_close  = []
                set_far    = []
                set_stable = []
                x1= []
                x2= []    
                
                for gene in move.index:
                    if move[gene] > (np.mean(move) + 2*np.std(move)):
                        set_far.append(gene)
                        x1.append(move[gene])
                    elif move[gene] < (np.mean(move) - 2*np.std(move)):
                        set_close.append(gene)
                        x2.append(move[gene])
                    else:
                        set_stable.append(gene)
                
                final = []
                final.extend(set_close)
                final.extend(set_far)
                
                compi = list(counter_list[counter_list["Counts"] == 0]["initial_alias"])
                compi =  [str(i) for i in compi]
                
                set(final).intersection(set(compi))
                
                x_final = []
                x_final.extend(x1)
                x_final.extend(x2)
                
                # Enrichment analyses:
                
                genes_succes = list(counter_list[counter_list.Counts > 0]["initial_alias"])
                genes_succes = [str(i) for i in genes_succes]
                            
                Total_succes = len(genes_succes)
                Total        = len(counter_list)
                
                # Close:
                
                set_close_succes = len(set(set_close).intersection(set(genes_succes)))
                set_close_total  = len(set_close)
                
                # Far:
                
                set_far_succes = len(set(set_far).intersection(set(genes_succes)))
                set_far_total  = len(set_far)    
                
                # Not move:
                
                set_stable_succes = len(set(set_stable).intersection(set(genes_succes)))
                set_stable_total  = len(set_stable) 
                
                # Both:
                
                both_succes = set_close_succes + set_far_succes
                both_total  = set_close_total + set_far_total
                
                percentage = (both_succes/both_total) * 100
                
                # Test:
                
                fold, p_value = Fold_enriched(Total, Total_succes, both_total, both_succes)
                
                
                p_values_stb.append(p_value)
                               
                the_file.write(f'{GO1}\t')
                the_file.write(f'{fold}\t')
                the_file.write(f'{p_value}\t')
                
                fold, p_value = Fold_enriched(Total, Total_succes, set_stable_total, set_stable_succes) 
                
                p_values_move.append(p_value)
                                
                the_file.write(f'{fold}\t')
                the_file.write(f'{p_value}\t')               
                the_file.write(f'{both_succes}\t')
                the_file.write(f'{percentage}\n')
                            
        # Close file:
        
            p_values_stb_perc = len([i for i in p_values_stb if i < 0.05])/len(p_values_stb)
            p_move_stb_perc = len([i   for i in p_values_move if i < 0.05])/len(p_values_move)
        
            the_file.write(f'Percentage_Stb_Agree: {p_values_stb_perc}\n')
            the_file.write(f'Percentage_Move_Agree: {p_move_stb_perc}\n')

        the_file.close()               


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
            

def Calculate_Clusters_GOGO(Cancer_type_list, Normal_Tissue_list, Cell_type_list, dimension_list, matrix = "PPMI", Common = True): 
    
    # Paths:

    save_cosine    = "/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/FMM/"
    Path_Semantic  = "/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/"
    Filter_path    = "/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/"
    path_organit   = "/media/sergio/sershiosdisk/Scripts/Main_Code/Cleaned_version/Data/FMM/"
       
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
            
            path_Cancer  = f'{save_cosine}Cosine_Cancer_{cancer}_{dim}_{matrix}_Leaf.csv'
            path_Control = f'{save_cosine}Cosine_Control_{tissue}_{cell}_{dim}_{matrix}_Leaf.csv'            
        
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
        f = open(f'{path_organit}Cancer_{tissue}_Cluster_SS_{matrix}.json' ,"w")
        f.write(sa)
        f.close() 

        sa = json.dumps(Result_Control)
        f = open(f'{path_organit}Control_{tissue}_Cluster_SS_{matrix}.json' ,"w")
        f.write(sa)
        f.close()  
            

        
def Fold_Intra_Inter_Cluster(Normal_Tissue_list, dimensions):
       
    # Paths:

    path_organit   = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"

    with open(f'{path_organit}Folds_Intra_Inter.txt', 'a') as the_file: 
        
        # Write the columns
        
        the_file.write(f'Tisse')
             
        for dim in dimensions:
            the_file.write(f'{dim}-D\t')
        the_file.write(f'{dim}-D\n')
        
        # Start the comparisons:
        
        for tissue in Normal_Tissue_list:
               
            file_control  = open(f'{path_organit}Control_{tissue}_Cluster_SS.json', "r") 
            file_cancer   = open(f'{path_organit}Cancer_{tissue}_Cluster_SS.json', "r") 
                
            info_control  = json.load(file_control)
            file_control.close()
            
            info_cancer  = json.load(file_cancer)
            file_cancer.close()
            
            the_file.write(f'Control {tissue}\t') 
            
            for dim in dimensions:
                
                key1 = f'{dim}_Inter'
                key2 = f'{dim}_Intra'
                
                Inter_SS_control = info_control[key1]
                Intra_SS_control = info_control[key2]
                
                the_file.write(f'{np.nanmean(Intra_SS_control)/np.nanmean(Inter_SS_control)}\t') 
            the_file.write(f'\n') 
            
            the_file.write(f'Cancer {tissue}\t')
            
            for dim in dimensions:
                
                key1 = f'{dim}_Inter'
                key2 = f'{dim}_Intra'
                
                Inter_SS_cancer = info_cancer[key1]
                Intra_SS_cancer = info_cancer[key2]
                
                the_file.write(f'{np.nanmean(Intra_SS_cancer)/np.nanmean(Inter_SS_cancer)}\t') 
            the_file.write(f'\n') 
               
    the_file.close()   

        
def Fold_Intra_Inter_Cluster_Paper(Normal_Tissue_list, dimensions, matrix = "PPMI"):
       
    # Paths:

    path_organit   = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"

    with open(f'{path_organit}Paper_Folds_Intra_Inter_{matrix}.txt', 'a') as the_file: 
        
        # Write the columns
        
        the_file.write(f'Tisse\t')
        the_file.write(f'Intra\t')
        the_file.write(f'Inter\t')
        the_file.write(f'Fold\t')
        the_file.write(f'pvalue\n')
     
        # Start the comparisons:
        
        x1 = []
        x2 = []
        
        for tissue in Normal_Tissue_list:
               
            file_control  = open(f'{path_organit}Control_{tissue}_Cluster_SS_{matrix}.json', "r") 
            file_cancer   = open(f'{path_organit}Cancer_{tissue}_Cluster_SS_{matrix}.json', "r") 
            
            file_control_PPMI  = open(f'{path_organit}Control_{tissue}_Cluster_SS.json', "r") 
            file_cancer_PPMI   = open(f'{path_organit}Cancer_{tissue}_Cluster_SS.json', "r")             
                
            info_control  = json.load(file_control)
            file_control.close()
            
            info_cancer  = json.load(file_cancer)
            file_cancer.close()
            
            info_control_PPMI  = json.load(file_control_PPMI)
            file_control_PPMI.close()
            
            info_cancer_PPMI  = json.load(file_cancer_PPMI)
            file_cancer_PPMI.close()  
            
            x1.extend(info_control_PPMI["200_Intra"])
            x1.extend(info_cancer_PPMI["200_Intra"])
            
            x2.extend(info_control["200_Inter"])
            x2.extend(info_cancer["200_Inter"])
            
            X = mannwhitneyu(info_control_PPMI['200_Intra'], info_control['200_Intra'], alternative = "greater").pvalue
            Y= mannwhitneyu(info_control_PPMI['200_Inter'], info_control['200_Inter'], alternative = "greater").pvalue

            Z = mannwhitneyu(info_cancer_PPMI['200_Intra'], info_cancer['200_Intra'], alternative = "greater").pvalue
            K = mannwhitneyu(info_cancer_PPMI['200_Inter'], info_cancer['200_Inter'], alternative = "greater").pvalue 
            
            print(f'Inter_Control {tissue} {X}')
            print(f'Intra_Control {tissue} {Y}')
            print(f'Inter_Cancer {tissue} {Z}')
            print(f'Intra_Cancer {tissue} {K}')
   
            the_file.write(f'Control {tissue}\t')
                       
            for dim in [200]:
                
                key1 = f'{dim}_Inter'
                key2 = f'{dim}_Intra'
                
                Inter_SS_control = info_control[key1]
                Intra_SS_control = info_control[key2]
                
                the_file.write(f'{np.nanmean(Intra_SS_control)}\t')
                the_file.write(f'{np.nanmean(Inter_SS_control)}\t')
                the_file.write(f'{np.nanmean(Intra_SS_control)/np.nanmean(Inter_SS_control)}\t')
                
                x = mannwhitneyu(Intra_SS_control,Inter_SS_control, alternative = "greater").pvalue
                
                the_file.write(f'{x}\t')
                
                
            the_file.write(f'\n') 
            
            the_file.write(f'Cancer {tissue}\t')
            
            for dim in [200]:
                
                key1 = f'{dim}_Inter'
                key2 = f'{dim}_Intra'
                
                Inter_SS_cancer = info_cancer[key1]
                Intra_SS_cancer = info_cancer[key2]
                
                the_file.write(f'{np.nanmean(Intra_SS_cancer)}\t')
                the_file.write(f'{np.nanmean(Inter_SS_cancer)}\t')                
                the_file.write(f'{np.nanmean(Intra_SS_cancer)/np.nanmean(Inter_SS_cancer)}\t')
                
                x = mannwhitneyu(Intra_SS_cancer,Inter_SS_cancer, alternative = "greater").pvalue
                
                the_file.write(f'{x}\t')
                
            the_file.write(f'\n') 
               
    the_file.close()               


def Compute_Similarity(ontology, tissue, sample, Semantic = "Resnik"):
    
    save_path      = "/media/sergio/sershiosdisk/Data_Driven_Ontology"
    
    if Semantic == "Resnik":
        sim, axes = ontology.flatten() 
        np.save(f'{save_path}Resnik_{tissue}_{sample}.npy', sim)
        np.save(f'{save_path}Resnik_Axes_list_{sample}.npy', axes)
                                                      

def Random_Analyses(Cancer_type_list, Normal_Tissue_list, Cell_type_list, times = 1):
    
    # For each cancer type:
    
    network_path  = "/media/sergio/sershiosdisk/Cancer/Networks/"
    Path_Semantic = "/media/sergio/sershiosdisk/Axes_Species/Semantic_Similarity/"
    Filter_path   = "/media/sergio/sershiosdisk/Cancer/Common_Set/"
    path_organit  = "/media/sergio/sershiosdisk/Cancer/Functional_Organization/"
    
    # Load semantic similarity:
    
    Semantic = pd.read_csv(f'{Path_Semantic}Semantic_Human.csv', index_col = 0)
    
    with open(f'{path_organit}Folds_Intra_Inter_Random.txt', 'a') as the_file: 
    
        for cancer, tissue, cell in zip(Cancer_type_list,Normal_Tissue_list,Cell_type_list):
            
            the_file.write(f'Tisse\t')
            the_file.write(f'P-value\t')
            the_file.write(f'Intra\t')
            the_file.write(f'Inter\n')
            
            # For each repetition in the randomization:
            
            for i in range(times):
                
                # Load the Adj matrix and Randomize the adjancencies:
    
                Cancer_adj  = np.load(f'{network_path}Cancer_{cancer}_PPI.npy',         allow_pickle= True)
                Control_adj = np.load(f'{network_path}Control_{tissue}_{cell}_PPI.npy', allow_pickle= True)
                
                Cancer_adj  = np.random.randint(2, size=(len(Cancer_adj), len(Cancer_adj)))
                Control_adj = np.random.randint(2, size=(len(Control_adj), len(Control_adj)))
                
                np.fill_diagonal(Cancer_adj,  0)
                np.fill_diagonal(Control_adj, 0)
                
                # Generate the PPMI for the corresponding random Adj:
                
                PPMI_Cancer  = generate_embedding_spaces_Human.deep_walk_ppmi(Cancer_adj,  10)
                PPMI_Control = generate_embedding_spaces_Human.deep_walk_ppmi(Control_adj, 10)
                
                # Generate the embeddings:
                
                save_path = "/media/sergio/sershiosdisk/Cancer/Random/"
                
                Solver = generate_embedding_spaces.SNMTF(100, 50)   
                Solver.Solve_MUR(PPMI_Cancer,  200, save_path, network = f'Random_PPI_{cancer}', matrix = "PPMI_Cancer",  init="SVD")
                Solver.Solve_MUR(PPMI_Control, 200, save_path, network = f'Random_PPI_{cancer}', matrix = "PPMI_Control",  init="SVD")
                     
                # Generate the GO embeddings:
                
                # Read the annotations:
                
                GO_Matrix = pd.read_csv(f'/media/sergio/sershiosdisk/Human/Leaf_Annotation/_Matrix_Genes_GO_BP_PPI.csv', 
                                        index_col = 0, dtype={0: str})  
                GO_Matrix.index = GO_Matrix.index.astype(str) 
                
                # Read genes:
                
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
                
                # Fix the gene vectorial representation:
            
                G_Cancer_PPMI  = np.load(f'{save_path}_G_Matrix_200_Random_PPI_{cancer}_PPMI_Cancer.npy', allow_pickle=True)
                G_Control_PPMI = np.load(f'{save_path}_G_Matrix_200_Random_PPI_{cancer}_PPMI_Control.npy', allow_pickle=True)
    
                # Do the embeddings:
                
                GO_Embeddings.GO_Embeddings_Human(Annotation_Cancer,  G_Cancer_PPMI, 
                                                  save_path,  GO = "Leaf", network = f'Random_PPI_{cancer}', 
                                                  matrix = "PPMI_Cancer",  ortogonal = True)
                
                GO_Embeddings.GO_Embeddings_Human(Annotation_Control,  G_Control_PPMI, 
                                                  save_path,  GO = "Leaf", network = f'Random_PPI_{tissue}_{cell}', 
                                                  matrix = "PPMI_Control",  ortogonal = True)
                
                # Calculate the FMMs of the random embeddings:
                
                path_Cancer  = f'{save_path}_GO_Embeddings_Leaf_Random_PPI_{cancer}_200_PPMI_Cancer.csv'
                path_Control = f'{save_path}_GO_Embeddings_Leaf_Random_PPI_{tissue}_{cell}_200_PPMI_Control.csv'                    
                
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
                
                # Calculate intra inter cluster distance:
               
                # Filter the cosines:
                
                common_set = list(pd.read_csv(f'{Filter_path}Common_Set_{tissue}_Leaf.csv')["0"])
                
                cosine_Cancer  = cosine_Cancer.loc[common_set, common_set]
                cosine_Control = cosine_Control.loc[common_set, common_set]
                
                # K medoids:
                
                n_clusters_1 = round(math.sqrt(len(cosine_Cancer.index)/2))
                n_clusters_2 = round(math.sqrt(len(cosine_Control.index)/2))
                
                grouping_1 = KMedoids(n_clusters=n_clusters_1, metric='precomputed',).fit(cosine_Cancer)
                grouping_2 = KMedoids(n_clusters=n_clusters_2, metric='precomputed',).fit(cosine_Control)
    
                Cluster_db_1 = pd.DataFrame(cosine_Cancer.index)  
                Cluster_db_2 = pd.DataFrame(cosine_Control.index) 
    
                Cluster_db_1.columns = ["GO"]    
                Cluster_db_2.columns = ["GO"]     
                
                Cluster_db_1["Cluster"] = grouping_1.labels_
                Cluster_db_2["Cluster"] = grouping_2.labels_ 
                
                # For Cancer:
                
                intra_distance_1, inter_distance_1 = Calculater_Intra_Inter_SS(Semantic, Cluster_db_1)
                intra_distance_2, inter_distance_2 = Calculater_Intra_Inter_SS(Semantic, Cluster_db_2)
                                  
                x = mannwhitneyu(intra_distance_1,inter_distance_1, alternative = "greater").pvalue
                y = mannwhitneyu(intra_distance_2,inter_distance_2, alternative = "greater").pvalue
                
                print(f'{cancer}{x}_{i}: Intra: {np.nanmean(intra_distance_1)} / Inter: {np.nanmean(inter_distance_1)}')
                print(f'{tissue}{y}_{i}: Intra: {np.nanmean(intra_distance_2)} / Inter: {np.nanmean(inter_distance_2)}')
                
                the_file.write(f'Cancer_{tissue}\t')
                the_file.write(f'{x}\t')
                the_file.write(f'{intra_distance_1}\t')
                the_file.write(f'{inter_distance_1}\n')
                
                the_file.write(f'Control_{tissue}\t')
                the_file.write(f'{y}\t')
                the_file.write(f'{intra_distance_2}\t')
                the_file.write(f'{inter_distance_2}\n')                
                             
    the_file.close()            
            
            
            
            

 
    
def Calcultate(Cancer_list, Normal_tissue_list, Control_list, dimensions = 150):   
        
    Path_Semantic   = "/media/sergio/sershiosdisk/Cancer/Semantic_Similarity/"
    Path_Cosine     = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
           
    for cancer, tissue, cell in zip(Cancer_list,Normal_tissue_list,Control_list): 
        
        # Load the cosine:
        
        Control_Structure  = pd.read_csv(f'{Path_Cosine}Cosine_Control_{tissue}_{cell}_{dimensions}_PPMI.csv', index_col = 0)
        Cancer_Structure   = pd.read_csv(f'{Path_Cosine}Cosine_Cancer_{cancer}_{dimensions}_PPMI.csv', index_col = 0) 
        
        GO_1 = list(Control_Structure.index)
        GO_1.extend(list(Cancer_Structure.index))
        GO_1 = list(set(GO_1))
        
        # Do an empty semantic matrix:
        
        Semantic_db       = pd.DataFrame(0, index=np.arange(len(GO_1)), columns=GO_1)
        Semantic_db.index = GO_1
        
        # Calculate the Semantic Similarity:
        
        G = graph.from_resource("go-basic")
        similarity.precalc_lower_bounds(G)
        
        print(tissue)
        
        for i in range(len(Semantic_db)):
            for j in range(i+1, len(Semantic_db)):
                
                try:
                    distance = similarity.lin(G, GO_1[i], GO_1[j]) 
                    
                    Semantic_db.loc[GO_1[i],GO_1[j]] = distance 
                    Semantic_db.loc[GO_1[j],GO_1[i]] = distance
                                     
                except Exception as PGSSLookupError:
                        continue  
                
        # Save the corresponding matrix:
        
        Semantic_db.to_csv(f'{Path_Semantic}Semantic_{tissue}.csv')    



        
                       
def Calculate_Similarity_matrix(Cancer_list, Normal_tissue_list, Control_list, dimensions = 150, annotation = "Leaf"): 
    
    # Paths:
    
    Path_Semantic   = "/media/sergio/sershiosdisk/Axes_Species/Semantic_Similarity/"
    Path_Cosine     = "/media/sergio/sershiosdisk/Cancer/Optimal_Dimensionality/"
    Path_Annotation = "/media/sergio/sershiosdisk/Axes_Species/Annotations/"
            
    # Load the BP annotations:
    
    if annotation == "Leaf":   
        BP = json.load(open("/media/sergio/sershiosdisk/Human/Annotation/gene2go_Human_PPIGO_Specific_BP.json"))
    else:    
        BP = json.load(open(f'{Path_Annotation}_BP_Back_Human.json'))
    
    # Load the Semantic Similarity:
    
    Semantic = pd.read_csv(f'{Path_Semantic}Semantic_Human.csv', index_col = 0)
    
    for cancer, tissue, cell in zip(Cancer_list,Normal_tissue_list,Control_list):
        
        # Load the cosine:
        
        if annotation == "Leaf":        
            Control_Structure  = pd.read_csv(f'{Path_Cosine}Cosine_Control_{tissue}_{cell}_{dimensions}_PPMI_Leaf.csv', index_col = 0)
            Cancer_Structure   = pd.read_csv(f'{Path_Cosine}Cosine_Cancer_{cancer}_{dimensions}_PPMI_Leaf.csv', index_col = 0)
        else:        
            Control_Structure  = pd.read_csv(f'{Path_Cosine}Cosine_Control_{tissue}_{cell}_{dimensions}_PPMI.csv', index_col = 0)
            Cancer_Structure   = pd.read_csv(f'{Path_Cosine}Cosine_Cancer_{cancer}_{dimensions}_PPMI.csv', index_col = 0)            
        
        # Filter the cosine:
        
        number_GO = 3
        annotation_list = [name for sublist in BP.values() for name in sublist]
        occurance_of_each_annotation_in_network = Counter(annotation_list)
        terms_filtering = [key for key,values in occurance_of_each_annotation_in_network.items() if values >= number_GO]
        
        terms_filtering_1 = list(set(terms_filtering).intersection(set(Cancer_Structure.index)))
        terms_filtering_2 = list(set(terms_filtering).intersection(set(Control_Structure.index)))
        
        Control_Structure = Control_Structure.loc[terms_filtering_2,terms_filtering_2]
        Cancer_Structure = Cancer_Structure.loc[terms_filtering_1,terms_filtering_1]
        
        # Filter the SS:
        
        Semantic_control = Semantic.loc[terms_filtering_2,terms_filtering_2]
        Semantic_cancer  = Semantic.loc[terms_filtering_1,terms_filtering_1]
        
        # Get the two vectors:
        
        semantic_vec_control =  np.array(Semantic_control)[np.triu_indices(np.array(Semantic_control).shape[0], k = 1)]
        semantic_vec_cancer  =  np.array(Semantic_cancer)[np.triu_indices(np.array(Semantic_cancer).shape[0], k = 1)]
        
        cosine_vec_control = np.array(Control_Structure)[np.triu_indices(np.array(Control_Structure).shape[0], k = 1)]
        cosine_vec_cancer  = np.array(Cancer_Structure)[np.triu_indices(np.array(Cancer_Structure).shape[0], k = 1)]
        
        # Start the comparisons (do a categorical variable):
                
        control_pd = pd.DataFrame()
        cancer_pd  = pd.DataFrame()
        
        control_pd["SS"] = semantic_vec_control
        control_pd["PP"] = cosine_vec_control
        control_pd["SS"] = control_pd["SS"].round(1)
        control_pd       = control_pd.sort_values("SS")
        control_pd["SS"] = control_pd["SS"].astype(str)
        
        cancer_pd["SS"] = semantic_vec_cancer
        cancer_pd["PP"] = cosine_vec_cancer
        cancer_pd["SS"] = cancer_pd["SS"].round(1)
        cancer_pd       = cancer_pd.sort_values("SS")
        cancer_pd["SS"] = cancer_pd["SS"].astype(str)
        
        x_control = list(set(control_pd["SS"]))
        x_cancer  = list(set(cancer_pd["SS"]))
        
        x_control.sort()
        x_cancer.sort()
        
        mean_control = []
        std_control  = []
                    
        for i in x_control:
            count = control_pd[control_pd["SS"] == i]
            mean_control.append(np.mean(count["PP"]))
            std_control.append(np.std(count["PP"]))        
        
        mean_cancer = []
        std_cancer  = []    
        
        for i in x_cancer:
            count = cancer_pd[cancer_pd["SS"] == i]
            mean_cancer.append(np.mean(count["PP"]))
            std_cancer.append(np.std(count["PP"]))  
            
        plt.style.use("seaborn-whitegrid")
        plt.rcParams.update({'font.size': 17})
        plt.rc('font', weight='bold')
        
        fig, axes = plt.subplots(figsize=(7,5))
        fig.tight_layout(pad = 0.5)
            
        axes.errorbar(x_control, mean_control,   yerr= Calculate_limits(mean_control, std_control), marker = "o", capsize=5, lw = 3, color = "#f1948a")
        axes.errorbar(x_cancer, mean_cancer,   yerr= Calculate_limits(mean_cancer, std_cancer), marker = "o", capsize=5, lw = 3, color = "#76de61")
                     
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['left'].set_linewidth(4)
        axes.spines['bottom'].set_linewidth(4)
        axes.spines['left'].set_color("grey")
        axes.spines['bottom'].set_color("grey")  
        
        axes.set_ylabel("Cosine Distance", fontsize=20, fontweight='bold') 
        axes.set_xlabel("SS", fontsize=20, fontweight='bold') 
