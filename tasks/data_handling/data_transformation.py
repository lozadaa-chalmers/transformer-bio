import h5py
import scipy as scp

'''
read_data(filename, mode)
Function that transforms h5-file to a sparse matrix.

filename: str. Path to h5-file to read. Default None.
mode: str. 'csr' for Compressed Sparse Row or 'csc' for Compressed Sparse Column. Default 'csr'.

Returns: sparse reconstructed matrix, tuples of form ((row, col), value).

Written: ronjah@chalmers.se
'''


def read_data(filename: str = None, mode: str = 'csr') -> scp.sparse.csr_matrix:
    with h5py.File(filename, "r") as f:
        express_handle = f["exprs"]
        shape = express_handle.get('shape')
        indices = express_handle.get('indices')

        if mode == 'csr':
            if not max(indices) <= shape[1]:
                raise ValueError('Not CSR format.')
            reconstructed_matrix = scp.sparse.csr_matrix((express_handle.get('data'),
                                                          express_handle.get('indices'),
                                                          express_handle.get('indptr')),
                                                         shape=express_handle.get('shape'))
        elif mode == 'csc':
            if not max(indices) <= shape[0]:
                raise ValueError('Not CSC format.')
            reconstructed_matrix = scp.sparse.csc_matrix((express_handle.get('data'),
                                                          express_handle.get('indices'),
                                                          express_handle.get('indptr')),
                                                         shape=express_handle.get('shape'))
        else:
            raise ValueError(f'Choose a right value for mode. Current {mode}. Needed: [csr,csc]')

    return reconstructed_matrix
