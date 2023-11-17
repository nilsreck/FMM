""" Non-negative matrix tri-factorization (numpy)"""

import numpy as np
from scipy.linalg import svd
from math import sqrt

#  Input: R, a matrix relating n genes to n genes. Entry R_ij = 1 if gene i interacts with gene j, 0 otherwise

#  Partial NMTF Decomposition

#           m×m                m×k_2
#     ┌───────────────┐        ┌───┐
#     │               │        │   │     k_2×k_2       k_2×m
#     │               │        │   │     ┌───┐   ┌──────────────┐
#     │      R1       │    ≈   │ P │  x  │ U │ x │       G'     │
#     │               │        │   │     └───┘   └──────────────┘
#     │               │        │   │
#     │               │        │   │
#     └───────────────┘        └───┘
#
#
#  With G,U,P >= 0
#  With the added contraints of orthogonality: G'G = I
#


def deep_walk_ppmi(adj_matrix, context_window=10):
    '''
    Input: Adj Matrix as numpy array
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


class SVD_EMBEDDINGS(object):
    def __init__(self,data_matrix):
        self.U, self.S, self.V = svd(data_matrix, full_matrices=True)

    "Compute SVD Based Embeddings"
    def embedding_space_svd(self,k):
        S = [sqrt(i) for i in self.S]
        diag_singular_values = np.diag(S)
        factorized_matrix = np.dot(self.U[:,0:k],diag_singular_values[0:k,0:k])
        return factorized_matrix


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

    def Solve_MUR(self, R1, k1, file_name_array_path, network = "PPI", matrix = "PPMI", init="rand"):
        print("Starting NMTF")
        n,n=np.shape(R1)
        norm_R1 = np.linalg.norm(R1, ord='fro')

        if init == "rand":
            print(" * Using random initialization")
            P = np.random.rand(n,k1) + 1e-5
            U = np.random.rand(k1,k1) + 1e-5
            G = np.random.rand(n,k1) + 1e-5
        elif init == "Sam":
            # Initialize matrix factor P, U, and G with all 0.5
            print(" * Sam Test")
            P = np.ones((n,k1))
            U = np.zeros((k1,k1)) +  1e-5
            np.fill_diagonal(U,0.1)
            G = np.ones((n,k1))
            P =P*0.1
            G = G*0.1
        elif init == 'SVD':
            # Initialize matrix factor P, U, and G using SVD decomposition

            print(" * -- Eig decomposition on R1 to get PxUxG^T")

            J,K,L = svd(R1)
        
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
