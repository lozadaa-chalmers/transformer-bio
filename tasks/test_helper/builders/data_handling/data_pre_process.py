import random as rand

import h5py
import numpy as np
import scipy.sparse
from scipy.sparse import csc_array, random


def simulate_csc_h5_data(
        file_path: str = None,
        filename: str = 'simulated_csc_h5_data.h5'
        ) -> (scipy.sparse, tuple[int, int], str):
    """
    ### Function that simulates h5-file and corresponding matrix for comparison and testing of functions.
    ---
    #### Args:
        - file_path (str): Path to folder where h5-file will be stored. Default None
        - filename (str): Name of h5-file. Default 'simulated_csc_h5_data.h5'
    ---
    #### Returns:
        - Transpose of sparse matrix
        - str: Path to saved simulated h5-file
    """

    shape = (rand.randint(50, 1000), rand.randint(50, 1000))
    simulated_matrix = random(m=shape[0],
                              n=shape[1],
                              density=0.1)

    sparse_representation = csc_array(simulated_matrix)
    simulated_indices = sparse_representation.indices
    simulated_indptr = sparse_representation.indptr
    simulated_data = sparse_representation.data
    amount_mt = int(0.05*shape[0])
    gene_names = np.array(['g'] * shape[0], dtype=object)
    cell_names = [f'C{i}' for i in range(1, shape[1] + 1)]
    id_number = [f'Id{i}' for i in range(1, shape[0]+1)]
    feature_type = ['Gene Expression'] * shape[0]
    genome = ['GRCh38'] * shape[0]
    random_indices = rand.sample(range(shape[0]), amount_mt)

    # Change the values at the selected indices to 'MT'
    gene_names[np.array(random_indices, dtype=int)] = 'MT-'
    gene_names = gene_names.tolist()

    with h5py.File(f'{file_path}\\{filename}', 'w') as f:
        matrix = f.create_group('matrix')
        features = matrix.create_group('features')
        features.create_dataset('name', data=gene_names)
        features.create_dataset('id', data=id_number)
        features.create_dataset('feature_type', data=feature_type)
        features.create_dataset('genome', data=genome)
        matrix.create_dataset('barcodes', data=cell_names)
        matrix.create_dataset('data', data=simulated_data)
        matrix.create_dataset('indices', data=simulated_indices)
        matrix.create_dataset('indptr', data=simulated_indptr)
        matrix.create_dataset('shape', data=np.array([shape[0], shape[1]]))

    return simulated_matrix.T, f'{file_path}\\{filename}'
