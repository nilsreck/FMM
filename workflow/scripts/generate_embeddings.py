import sys
import numpy as np
import pandas as pd
from scipy.linalg import svd
from multiprocessing import Pool

# 3. Generate the Gene embeddings:

def solve(process):

    Solver = SNMTF(*process[0:2])
    dim = process[4]
    save_path = process[5]
    network = process[6]
    init = "SVD"
    ppmi_file = process[2] 
    PPMI = np.load(ppmi_file, allow_pickle=True)
    adj_file = process[3]
    adj = np.load(adj_file, allow_pickle=True)
    ppmi_matrix = f"PPMI_{process[-1]}"
    adj_matrix = f"Adj_{process[-1]}"

    Solver.Solve_MUR(
        PPMI,
        dim,
        save_path,
        network,
        ppmi_matrix,
        init)
    
    Solver.Solve_MUR(
        adj,
        dim,
        save_path,
        network,
        adj_matrix,
        init)

    return 0
              

class SNMTF(object):
    """Compute Partial Orthonormal Non-negative Matrix Tri-Factorization (NMTF) Embeddings"""

    def __init__(self, max_iter=1000, verbose=10):
        super(SNMTF, self).__init__()
        self.max_iter = max_iter
        self.verbose = verbose

    def Score(self, R1, P, U, G, norm_R1):
        GT = np.transpose(G)

        ErR1 = R1 - np.matmul(P, np.matmul(U, GT))
        norm_erR1 = np.linalg.norm(ErR1, ord="fro")

        rel_erR1 = norm_erR1 / norm_R1

        return norm_erR1, rel_erR1

    def Get_Clusters(self, M, nodes):
        n, k = np.shape(M)
        Clusters = [[] for i in range(k)]
        for i in range(n):
            idx = np.argmax(M[i])
            Clusters[idx].append(nodes[i])
        return Clusters

    def Solve_MUR(
        self, R1, k1, file_name_array_path, network="PPI", matrix="PPMI", init="rand"
    ):
        print("Starting NMTF")
        n, n = np.shape(R1)
        norm_R1 = np.linalg.norm(R1, ord="fro")

        if init == "rand":
            print(" * Using random initialization")
            P = np.random.rand(n, k1) + 1e-5
            U = np.random.rand(k1, k1) + 1e-5
            G = np.random.rand(n, k1) + 1e-5
        elif init == "Sam":
            # Initialize matrix factor P, U, and G with all 0.5
            print(" * Sam Test")
            P = np.ones((n, k1))
            U = np.zeros((k1, k1)) + 1e-5
            np.fill_diagonal(U, 0.1)
            G = np.ones((n, k1))
            P = P * 0.1
            G = G * 0.1
        elif init == "SVD":
            # Initialize matrix factor P, U, and G using SVD decomposition

            print(" * -- Eig decomposition on R1 to get PxUxG^T")

            J, K, L = svd(R1)

            P = np.zeros((n, k1)) + 1e-5
            G = np.zeros((n, k1)) + 1e-5
            U = np.zeros((k1, k1)) + 1e-5

            for i in range(n):
                for j in range(k1):
                    if J[i][j] > 0.0:
                        P[i][j] = J[i][j]

            for i in range(k1):
                U[i][i] = abs(K[i])

            for i in range(n):
                for j in range(k1):
                    if L[i][j] > 0.0:
                        G[i][j] = L[i][j]

        else:
            print("Unknown initializer: %s" % (init))
            exit(0)
        OBJ, REL2 = self.Score(R1, P, U, G, norm_R1)
        print(" - Init:\t OBJ:%.4f\t REL2:%.4f" % (OBJ, REL2))

        # Begining M.U.R.
        for it in range(1, self.max_iter + 1):
            R1T = np.transpose(R1)
            UT = np.transpose(U)

            PT = np.transpose(P)
            PT_P = np.matmul(PT, P)

            GT = np.transpose(G)
            GT_G = np.matmul(GT, G)

            # update rule for G

            # R1 matrix

            # update rule for P

            R1_G_UT = np.matmul(np.matmul(R1, G), UT)
            U_GT_G_UT = np.matmul(np.matmul(U, GT_G), UT)

            P_U_GT_G_UT = np.matmul(P, U_GT_G_UT)

            P_mult = np.sqrt(np.divide(R1_G_UT, P_U_GT_G_UT + 1e-14))

            # update rule for U

            PT_R1_G = np.matmul(PT, np.matmul(R1, G))
            PT_P_U_GT_G = np.matmul(PT_P, np.matmul(U, GT_G))

            U_mult = np.sqrt(np.divide(PT_R1_G, PT_P_U_GT_G + 1e-14))

            # update rule for G

            R1T_P_U = np.matmul(np.matmul(R1T, P), U)
            G_GT_R1T_P_U = np.matmul(G, np.matmul(GT, R1T_P_U))

            G_mult = np.sqrt(np.divide(R1T_P_U, G_GT_R1T_P_U + 1e-14))

            # Applying M.U.R.
            P = np.multiply(P, P_mult) + 1e-14
            U = np.multiply(U, U_mult) + 1e-14
            G = np.multiply(G, G_mult) + 1e-14

            if (it % self.verbose == 0) or (it == 1):
                OBJ, REL2 = self.Score(R1, P, U, G, norm_R1)
                print(" - It %i:\t OBJ:%.4f\t REL2:%.4f" % (it, OBJ, REL2))

        # Save the matrices:

        file_name_array_path_1 = (
            file_name_array_path
            + "_P_Matrix_"
            + str(k1)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
        )
        file_name_array_path_2 = (
            file_name_array_path
            + "_U_Matrix_"
            + str(k1)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
        )
        file_name_array_path_3 = (
            file_name_array_path
            + "_G_Matrix_"
            + str(k1)
            + "_"
            + str(network)
            + "_"
            + str(matrix)
        )

        np.save(file_name_array_path_1, P)
        np.save(file_name_array_path_2, U)
        np.save(file_name_array_path_3, G)

        return P, U, G


def Generate_Gene_Embeddings(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    data_dir,
    max_iter=1000,
    verbose=50,
):
    """
    This function generates the Gene Embeddings from the PPMI and Adj:

        Inputs:
                - Cancer_type     : string list, Name of the Cancer.
                - Normal_Tissue   : string list, Name of the normal tissue.
                - Cell_type       : string list, Name of the normal cell in the tissue.
                - dimension_list  : int list, dimensions to do the embeddings
                - context_window  : int, contex window for the PPMI matrix.
                - max_iter        : int, iterations for the NMTF.
                - verbose         : int.
    """
    # Networks path:

    network_path = f"{data_dir}/"
    save_path = f"{data_dir}/Embeddings/"

    counter = 0

    # Generate the embeddings:

    Processes = []

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Solver = SNMTF(max_iter, verbose)

        # Load the information:

        # Cancer_adj = np.load(
        #     f"{network_path}Cancer_{cancer}_PPI.npy", allow_pickle=True
        # )
        # Control_adj = np.load(
        #     f"{network_path}Control_{tissue}_{cell}_PPI.npy", allow_pickle=True
        # )

        # Cancer_PPMI = np.load(
        #     f"{network_path}Cancer_PPMI_{cancer}.npy", allow_pickle=True
        # )
        # Control_PPMI = np.load(
        #     f"{network_path}Control_PPMI_{tissue}_{cell}.npy", allow_pickle=True
        # )

        Cancer_adj = f"{network_path}Cancer_{cancer}_PPI.npy"
        Control_adj = f"{network_path}Control_{tissue}_{cell}_PPI.npy"

        Cancer_PPMI = f"{network_path}Cancer_PPMI_{cancer}.npy"
        Control_PPMI = f"{network_path}Control_PPMI_{tissue}_{cell}.npy"


        # Do the embeddings:

        
        Processes.append((
                max_iter, 
                verbose,
                Cancer_PPMI,
                Cancer_adj,
                dimension_list[0],
                save_path,
                f"PPI_{cancer}",
                "Cancer"))
        Processes.append((
                max_iter, 
                verbose,
                Control_PPMI,
                Control_adj,
                dimension_list[0],
                save_path,
                f"PPI_{tissue}_{cell}",
                "Control"))

        # Cancer_PPMI_process = (max_iter, verbose,
        #         Cancer_PPMI,
        #         dimension_list[0],
        #         save_path,
        #         f"PPI_{cancer}",
        #         "PPMI_Cancer",
        #         "SVD")
        # Cancer_adj_process = (max_iter, verbose,
        #         Cancer_adj,
        #         dimension_list[0],
        #         save_path,
        #         f"PPI_{cancer}",
        #         "Adj_Cancer",
        #         "SVD")
        # Control_PPMI_process = (max_iter, verbose,
        #         Control_PPMI,
        #         dimension_list[0],
        #         save_path,
        #         f"PPI_{tissue}_{cell}",
        #         "PPMI_Control",
        #         "SVD")
        # Control_adj_process = (max_iter, verbose,
        #         Control_adj, 
        #         dimension_list[0],
        #         save_path,
        #         f"PPI_{tissue}_{cell}",
        #         "Adj_Control",
        #         "SVD")

        # with Pool(4) as pool:
            # pool.map(solve, [Cancer_PPMI_process, Cancer_adj_process, Control_PPMI_process, Control_adj_process])          
            # pool.close()
            # pool.join()

        # for dim in dimension_list:
        #     Solver.Solve_MUR(
        #         Cancer_PPMI,
        #         dim,
        #         save_path,
        #         network=f"PPI_{cancer}",
        #         matrix="PPMI_Cancer",
        #         init="SVD",
        #     )
        #     Solver.Solve_MUR(
        #         Cancer_adj,
        #         dim,
        #         save_path,
        #         network=f"PPI_{cancer}",
        #         matrix="Adj_Cancer",
        #         init="SVD",
        #     )

        #     Solver.Solve_MUR(
        #         Control_PPMI,
        #         dim,
        #         save_path,
        #         network=f"PPI_{tissue}_{cell}",
        #         matrix="PPMI_Control",
        #         init="SVD",
        #     )
        #     Solver.Solve_MUR(
        #         Control_adj,
        #         dim,
        #         save_path,
        #         network=f"PPI_{tissue}_{cell}",
        #         matrix="Adj_Control",
        #         init="SVD",
        #     )

        # counter = counter + 4

        # print(f"{counter}/{len(Cancer_type_list) * 4 * (len(dimension_list))}")
    
    with Pool(8) as pool:
        pool.map(solve, Processes)          
        pool.close()
        pool.join()


Generate_Gene_Embeddings(
    snakemake.params.Cancer_list,
    snakemake.params.Normal_tissue_list,
    snakemake.params.Control_list,
    [int(open(snakemake.input.dimension, "w").name.split("_")[-1])],
    snakemake.params.data_path,
    max_iter=snakemake.params.max_iter,
    verbose=50,
)