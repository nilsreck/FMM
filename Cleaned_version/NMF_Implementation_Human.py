"""
This is a proof of concept of an NMF with the multiplicative update rules.
"""

import numpy as np
import pandas as pd
from scipy.linalg import svd
import matplotlib.pyplot as plt


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
