
import json

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch import nn
import h5py
import seaborn as sns

from example.scFormer_assets import utils
from tasks.data_handling import data_pre_processing as dpp

random_embeddings = True

our_data = dpp.create_count_matrix('data/pancreas.h5ad',
                                   make_genes_unique=True)

cell_embeddings = torch.rand(size=(our_data.shape[0], 2)).detach().cpu().numpy()

cell_embeddings = cell_embeddings / (np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
) + 0.000000000000001)
breakpoint()

cell_embeddings_anndata = anndata.AnnData(X=cell_embeddings)
breakpoint()
# our_data.obsm["X_scFormer_alejandro"] = cell_embeddings

sc.pp.neighbors(cell_embeddings_anndata)
breakpoint()
sc.tl.umap(cell_embeddings_anndata, min_dist=0.3)
breakpoint()

'''
batch_key = 'batch'
if not isinstance(our_data.obs[batch_key][0], str):
    cell_embeddings_anndata.obs["str_" + batch_key] = our_data.obs[batch_key].astype(str)
    batch_key = "str_" + batch_key

str_celltype_col = 'louvain'
celltype_col = 'int' + str_celltype_col
cell_embeddings_anndata.obs[str_celltype_col] = utils.category_str2int(
    our_data.obs[str_celltype_col]
)
'''

cell_embeddings_anndata.obs['batch'] = our_data.obs['batch'].values
cell_embeddings_anndata.obs['louvain'] = our_data.obs['louvain'].values
breakpoint()
fig = sc.pl.umap(
    cell_embeddings_anndata,
    color=['batch', 'louvain'],
    ncols=2,
    frameon=False,
    return_fig=True,
)
fig.savefig(
    f"random_cell_embeddings_umap-test-alc-pancreas-5.png",
    bbox_inches="tight",
)
breakpoint()
print('done')
