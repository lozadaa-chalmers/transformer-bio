import random as rand

import h5py
import numpy as np
import scipy.sparse
from scipy.sparse import csr_array, csc_array, random


def simulate_h5_data(mode: str = 'csr',
                     file_path: str = None,
                     filename: str = 'simulated_h5_data.h5') -> (scipy.sparse, tuple[int, int]):

    shape = (rand.randint(50, 1000), rand.randint(50, 1000))
    mode_format = {'csr': 0, 'csc': 1}[mode]
    mode_representation = {0: csr_array, 1: csc_array}[mode_format]
    simulated_matrix = random(m=shape[0],
                              n=shape[1],
                              density=0.1)

    sparse_representation = mode_representation(simulated_matrix)
    simulated_indices = sparse_representation.indices
    simulated_indptr = sparse_representation.indptr
    simulated_data = sparse_representation.data

    with h5py.File(f'{file_path}\\{filename}', 'w') as f:
        group1 = f.create_group('exprs')
        group1.create_dataset('data', data=simulated_data)
        group1.create_dataset('indices', data=simulated_indices)
        group1.create_dataset('indptr', data=simulated_indptr)
        group1.create_dataset('shape', data=np.array([shape[0], shape[1]]))

    return simulated_matrix
