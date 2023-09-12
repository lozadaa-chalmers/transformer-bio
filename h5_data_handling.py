import h5py
import numpy as np
import os

def transform_h5_to_matrix(h5_file_directory: str, count_matrix_directory: str, gene_names_directory: str):
    # Directory containing the HDF5 files
    directory = h5_file_directory

    # List all files with the .h5 extension in the directory
    h5_files = [f for f in os.listdir(directory) if f.endswith('.h5')]

    # Extract dataset names by removing the .h5 extension
    dataset_names = [os.path.splitext(filename)[0] for filename in h5_files]

    # List all files with the .h5 extension in the directory
    h5_files = [f for f in os.listdir(directory) if f.endswith('.h5')]

    for h5_filename in h5_files:
        # Construct the full path to the HDF5 file
        h5_filepath = os.path.join(directory, h5_filename)

        # Open the HDF5 file for reading
        with h5py.File(h5_filepath, 'r') as h5_file:
            dataset_name = os.path.splitext(h5_filename)[0]

            # Display information about the HDF5 file
            #print(h5_file.keys())  # Print the keys (datasets/groups)

            # Read data from the HDF5 file
            data = h5_file['/exprs/data'][()]
            obs_names = h5_file['/obs_names'][()]
            var_names = h5_file['/var_names'][()]
            indices = h5_file['/exprs/indices'][()]
            indptr = h5_file['/exprs/indptr'][()]
            shape = h5_file['/exprs/shape'][()]

            M = np.zeros(shape, dtype=np.float32)

            cell_id = 0
            for i in range(len(indices)):
                if i >= indptr[cell_id + 1]:
                    cell_id += 1

                if cell_id >= shape[0]:
                    break

                M[cell_id, indices[i]] = data[i]

            # Save the count matrix and gene names
            count_matrix_filename = os.path.join(count_matrix_directory, f'{dataset_name}_countmatrix.npy')
            gene_names_filename = os.path.join(gene_names_directory, f'{dataset_name}_genenames.npy')

            np.save(count_matrix_filename, M)
            np.save(gene_names_filename, var_names)