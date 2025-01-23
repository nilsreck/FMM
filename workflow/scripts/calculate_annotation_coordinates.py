import sys
import numpy as np
import pandas as pd
from scipy.linalg import svd
import matplotlib.pyplot as plt

# 4 Map Functional annotations (GO terms) to the embedding spaces:

def Score_NMF(A, U, V, norm_A):
    ErR1 = A - (U @ V)
    norm_erA = np.linalg.norm(ErR1, ord="fro")

    rel_erR1 = norm_erA / norm_A

    return norm_erA, rel_erR1


def NMF(A, k1, repetitions, init="SVD", fix_U=True, U=None, fix_V=False, V=None):
    """
    Imputs:
            A     : Data Frame Pandas (gene/GO matrix)
            k1    : Dimensions
            init  : Intialitation (random or singular value decomposition based)
            fix_U : If the matrix U is fixed from a previous decompotition
            U     : Numpy Array Matrix (only if fix_U is True)
            fix_V : Fix the matrix V (GO embeddings)
            V     : Numpy Array Matrix (only if fix_V is True)
    """

    # 1. Get the dimensions values:

    GO_Terms = A.columns  # GO terms
    genes = [str(i) for i in A.index]  # genes
    A = np.array(A)  #
    n, m = np.shape(A)
    norm_A = np.linalg.norm(A, ord="fro")

    # 2. Inizilizate the Matrices U and V:

    if init == "rand":
        print(" * Using random initialization")

        if fix_U == False:
            U = np.random.rand(n, k1) + 1e-5
            V = np.random.rand(m, k1) + 1e-5
            V = V.T
        elif fix_U == True:
            print(" * Fixing U matrix")
            U = U
            V = np.random.rand(m, k1) + 1e-5
            V = V.T
        elif fix_V == True:
            print(" * Fixing V matrix")
            U = np.random.rand(n, k1) + 1e-5
            V = V
            V = V.T

    elif init == "Sam":
        # Initialize matrix factor P, G with
        print(" * Sam Test")
        U = U
        V = np.ones((m, k1))
        V = V * 0.1
        V = V.T
    elif init == "SVD":
        print(" * -- Eig decomposition on R1 to get UxV^T")

        """
        Here we are using the whole annotation GO terms. If the user wants to change 
        to a different annotation, he would have to change the files to the corresponding
        SVD.
        """

        print("Loading the L and J matrices")

        J = np.load(
            "/media/sergio/sershiosdisk/Human/Node2Vec_Exp/_J_Matrix_Matrix_Gene_GO_Human.npy"
        )
        L = np.load(
            "/media/sergio/sershiosdisk/Human/Node2Vec_Exp/_L_Matrix_Matrix_Gene_GO_Human.npy"
        )

        if fix_U == False:
            U = np.zeros((n, k1)) + 1e-5
            V = np.zeros((m, k1)) + 1e-5

            for i in range(n):
                for j in range(k1):
                    if J[i][j] > 0.0:
                        U[i][j] = J[i][j]

            for i in range(m):
                for j in range(k1):
                    if L[i][j] > 0.0:
                        V[i][j] = L[i][j]

            V = V.T

        elif fix_U == True:
            print(" * Fixing U matrix")

            U = U
            V = np.zeros((m, k1)) + 1e-5

            for i in range(m):
                for j in range(k1):
                    if L[i][j] > 0.0:
                        V[i][j] = L[i][j]

            V = V.T

        elif fix_V == True:
            print(" * Fixing V matrix")

            U = np.zeros((n, k1)) + 1e-5
            V = V.T

            for i in range(n):
                for j in range(k1):
                    if J[i][j] > 0.0:
                        U[i][j] = J[i][j]

    # 3. Begining M.U.R.

    error = []
    objective = []

    if fix_U == False and fix_V == False:
        for it in range(1, repetitions + 1):
            VT = np.transpose(V)

            # Update rule for U:

            A_VT_T = A @ VT
            UT_VT_VT_T = (U @ V) @ VT
            Ut1 = np.divide(A_VT_T, UT_VT_VT_T + 1e-14)

            U = np.multiply(U, Ut1) + 1e-14

            # Update rule for V:

            UT_A = U.T @ A
            UT_U_V = U.T @ U @ V
            Vt1 = np.divide(UT_A, UT_U_V + 1e-14)

            # Update:

            V = np.multiply(V, Vt1) + 1e-14

            # 4. Check the error:

            OBJ, REL2 = Score_NMF(A, U, V, norm_A)
            print(" - It %i:\t OBJ:%.4f\t REL2:%.4f" % (it, OBJ, REL2))

            error.append(REL2)
            objective.append(OBJ)

            # Plot the results:

            # 4. Check the error:

            if it == 1 or it % 10 == 0:
                OBJ, REL2 = Score_NMF(A, U, V, norm_A)
                print(" - It %i:\t OBJ:%.4f\t REL2:%.4f" % (it, OBJ, REL2))

            error.append(REL2)
            objective.append(OBJ)

        plt.plot(error)
        plt.plot(objective)

        V = pd.DataFrame(V)
        V.columns = GO_Terms

        return V

    elif fix_U == True and fix_V == False:  # This is the one we use
        for it in range(1, repetitions + 1):
            VT = np.transpose(V)

            # Update rule for V:

            UT_A = U.T @ A
            UT_U_V = U.T @ U @ V
            Vt1 = np.divide(UT_A, UT_U_V + 1e-14)

            # Update:

            V = np.multiply(V, Vt1) + 1e-14

            # 4. Check the error:

            if it == 1 or it % 10 == 0:
                OBJ, REL2 = Score_NMF(A, U, V, norm_A)
                print(" - It %i:\t OBJ:%.4f\t REL2:%.4f" % (it, OBJ, REL2))

            error.append(REL2)
            objective.append(OBJ)

        plt.plot(error)
        plt.plot(objective)

        V = pd.DataFrame(V)
        V.columns = GO_Terms

        return V

    elif fix_V == True and fix_U == False:  # Here we fix V to get U.
        for it in range(1, repetitions + 1):
            VT = np.transpose(V)

            # Update rule for U:

            A_VT_T = A @ VT
            UT_VT_VT_T = (U @ V) @ VT
            Ut1 = np.divide(A_VT_T, UT_VT_VT_T + 1e-14)

            # Update:

            U = np.multiply(U, Ut1) + 1e-14

            # 4. Check the error:

            # Plot the results:

            if it == 1 or it % 10 == 0:
                OBJ, REL2 = Score_NMF(A, U, V, norm_A)
                print(" - It %i:\t OBJ:%.4f\t REL2:%.4f" % (it, OBJ, REL2))

            error.append(REL2)
            objective.append(OBJ)

        plt.plot(error)
        plt.plot(objective)

        U = pd.DataFrame(U)
        U.index = genes

        return U


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
    Inputs:

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
        GO_embeddings = NMF(
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

def Annotate_Gene_Space(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    data_dir,
    Annotation="GO",
):
    """
    This function generates the Annotation of the gene embedding space:

        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - context_window  : int, contex window for the PPMI matrix.
                - Annotation      : string, annotation for the space.
    """

    # Networks path:

    network_path = f"{data_dir}/"
    embedding_path = f"{data_dir}/Embeddings/"
    save_path = f"{data_dir}/Embeddings/"

    counter = 0
    GO_Matrix = None

    # Get the annotations:

    if Annotation == "GO":
        GO_Matrix = pd.read_csv(
            f"{network_path}_Matrix_Genes_GO_BP_Back_Propagation_PPI.csv",
            index_col=0,
            dtype={0: str},
        )
        GO_Matrix.index = GO_Matrix.index.astype(str)

        print("nop")

    elif Annotation == "Leaf":
        GO_Matrix = pd.read_csv(
            f"{data_dir}/_Matrix_Genes_GO_BP_PPI.csv",
            index_col=0,
            dtype={0: str},
        )
        GO_Matrix.index = GO_Matrix.index.astype(str)

    else:
        print("Error Check the paths to your annotation matrix")

    # Start the embeddings:

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        Cancer_genes = pd.read_table(
            f"{network_path}Cancer_{cancer}_Genes.csv", header=0, dtype={0: str}
        )
        Cancer_genes_list = list(Cancer_genes["0"])

        Control_genes = pd.read_table(
            f"{network_path}Control_{tissue}_{cell}_Genes.csv", header=0, dtype={0: str}
        )
        Control_genes_list = list(Control_genes["0"])

        # Subset the annotation to keep only the genes:

        GO_Gene_Matrix_Cancer = GO_Matrix[GO_Matrix.index.isin(Cancer_genes_list)]
        GO_Gene_Matrix_Control = GO_Matrix[GO_Matrix.index.isin(Control_genes_list)]

        # Delete annotations without genes (> 0):

        GO_Cancer = GO_Gene_Matrix_Cancer.sum(axis=0)
        filter_Cancer = set(GO_Cancer[GO_Cancer > 0].index)

        GO_Control = GO_Gene_Matrix_Control.sum(axis=0)
        filter_Control = set(GO_Control[GO_Control > 0].index)

        GO_Gene_Matrix_Cancer = GO_Gene_Matrix_Cancer[filter_Cancer]
        GO_Gene_Matrix_Control = GO_Gene_Matrix_Control[filter_Control]

        # Reorder the Matrices:

        Annotation_Cancer = GO_Gene_Matrix_Cancer.loc[Cancer_genes_list]
        Annotation_Control = GO_Gene_Matrix_Control.loc[Control_genes_list]

        for dim in dimension_list:
            # Fix the gene vectorial representation:

            G_Cancer_PPMI = np.load(
                f"{embedding_path}_G_Matrix_{dim}_PPI_{cancer}_PPMI_Cancer.npy",
                allow_pickle=True,
            )
            G_Cancer_Adj = np.load(
                f"{embedding_path}_G_Matrix_{dim}_PPI_{cancer}_Adj_Cancer.npy",
                allow_pickle=True,
            )

            G_Control_PPMI = np.load(
                f"{embedding_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_PPMI_Control.npy",
                allow_pickle=True,
            )
            G_Control_Adj = np.load(
                f"{embedding_path}_G_Matrix_{dim}_PPI_{tissue}_{cell}_Adj_Control.npy",
                allow_pickle=True,
            )

            # Do the embeddings:

            GO_Embeddings_Human(
                Annotation_Cancer,
                G_Cancer_PPMI,
                save_path,
                GO=Annotation,
                network=f"PPI_{cancer}",
                matrix="PPMI_Cancer",
                ortogonal=True,
            )

            GO_Embeddings_Human(
                Annotation_Cancer,
                G_Cancer_Adj,
                save_path,
                GO=Annotation,
                network=f"PPI_{cancer}",
                matrix="Adj_Cancer",
                ortogonal=True,
            )

            GO_Embeddings_Human(
                Annotation_Control,
                G_Control_PPMI,
                save_path,
                GO=Annotation,
                network=f"PPI_{tissue}_{cell}",
                matrix="PPMI_Control",
                ortogonal=True,
            )

            GO_Embeddings_Human(
                Annotation_Control,
                G_Control_Adj,
                save_path,
                GO=Annotation,
                network=f"PPI_{tissue}_{cell}",
                matrix="Adj_Control",
                ortogonal=True,
            )

            counter = counter + 4

            print(f"{counter}/{len(Cancer_type_list) * 4 * (len(dimension_list))}")


Annotate_Gene_Space(
    snakemake.params.Cancer_list, 
    snakemake.params.Normal_tissue_list, 
    snakemake.params.Control_list, 
    snakemake.params.dimension_list, 
    snakemake.params.data_path,
    Annotation=snakemake.params.annotation
)