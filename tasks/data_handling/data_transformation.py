import h5py
import scipy as scp
'''
read_data(filename, csr)
Function that transforms h5-file to a sparse matrix

- str filename: path to h5-file to read.
- bool csr: True if CSR-format. False if CSC-format.

Returns: sparse reconstructed matrix, tuples of form ((row, col), value)
'''

def read_data(filename: str = None, csr: bool = True) -> scp.sparse.csr_matrix:
    with h5py.File(filename, "r") as f:
        express_handle = f["exprs"]

        shape = express_handle.get('shape')
        n_genes = len(f["var_names"])
        n_cells = len(f["obs_names"])

        if csr:
            if n_genes != shape[1] or n_cells != shape[0]:
                raise ValueError('Not CSR format.')
        else:
            if n_genes != shape[0] or n_cells != shape[1]:
                raise ValueError('Not CSC format.')
            
        reconstructed_matrix = scp.sparse.csr_matrix((express_handle.get('data'),
                                                      express_handle.get('indices'),
                                                      express_handle.get('indptr')),
                                                      shape=express_handle.get('shape'))

    return reconstructed_matrix

