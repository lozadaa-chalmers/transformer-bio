import scanpy as sc
import numpy as np
import os
from typing import Tuple, List


def find_common_genes(
        folder_path: str = None,
        save_path: str = None
) -> List:
    """
    ### Function that returns a list of genes available in every h5-files found in folder_path

    ---

    #### Args:
        - folder_path (str): Path to folder with all .h5-files. Default None
        - save_path (str): Where to save the list, if None it is not saved. Default None

    ---
    #### Returns:
        - List of genes that are in all .h5-files from folder_path.
    ---
    Written: ronjah@chalmers.se
    """
    files = os.listdir(folder_path)
    file_names = [file for file in files if file.endswith('.h5')]
    adata = sc.read_h5ad(folder_path + '/' + file_names[0])
    common_genes = adata.var_names.to_list()
    # Checks intersection of genes between .h5-files in folder_path
    for file in file_names:
        adata = sc.read_h5ad(folder_path + '/' + file)
        gene_list_file = adata.var_names.to_list()
        common_genes_set = set(common_genes)
        common_genes = set(gene_list_file) & common_genes_set
    if save_path is not None:
        np.save(save_path, common_genes)
    return common_genes


def get_new_gene2vec(
        gene2vec_path: str = None,
        h5_data_path: str = None,
        common_gene_list: List = None,
        save_path_gene2vec: str = None,
        save_path_finished_gene_list: str = None
) -> Tuple[np.ndarray, List]:
    """
    ### Function that returns a Tuple of the numpy array of the resulting gene2vec and a list of genes in
        corresponding order, alphabetical.

    ---

    #### Args:
        - gene2vec_path (str): Path to the full gene2vec. Default None
        - h5_data_path (str): Path to .h5-file corresponding to gene2vec. Default None
        - common_gene_list (List): Gene list of all genes available in every .h5-file from find_common_genes.
            Default None
        - save_path_gene2vec (str): Save path for the new gene2vec, if None it is not saved. Default None
        - save_path_finished_gene_list (str): Save path for the corresponding gene list, if None it is not saved.
            Default None

    ---
    #### Returns:
        - Numpy array of the new gene2vec
        - List of all corresponding genes in order.
    ---
    Written: ronjah@chalmers.se
    """
    gene2vec = np.load(gene2vec_path)
    adata_h5 = sc.read_h5ad(h5_data_path)
    full_list_genes = adata_h5.var_names.to_list()
    available_genes = set(common_gene_list) & set(full_list_genes)
    # Sort gene_list in alphabetical order
    ordered_common_gene_list = sorted(available_genes)
    new_gene2vec = np.empty([len(ordered_common_gene_list), 200])
    i = 0
    # Create new gene2vec in correct order
    for gene in ordered_common_gene_list:
        index = full_list_genes.index(gene)
        new_gene2vec[i] = gene2vec[index]
        i += 1
    if save_path_gene2vec is not None:
        np.save(save_path_gene2vec, new_gene2vec)
    if save_path_finished_gene_list is not None:
        np.save(save_path_finished_gene_list, ordered_common_gene_list)
    return new_gene2vec, ordered_common_gene_list
