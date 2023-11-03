"""

We'll make a first run of ScFormer using our own data. Specifically '3013_P_filtered_feature_bc_matrix.h5'

"""
# TODO: look into what is being done with 'sc.pp.neighbors' and 'sc.tl.umap'
# TODO: check what's going on with read_10X 'genome' option.
# TODO: check if random embeddings are as good as the 'pretrained' ones.
# TODO: check if there's any difference between cell embedding strategies.
# TODO: normalize our  data to see if it changes the clustering.

import json

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch import nn

from scFormer_assets import utils
from tasks.data_handling import data_pre_processing as dpp

random_embeddings = True

our_data = dpp.create_count_matrix('data/pancreas.h5ad',
                                   make_genes_unique=True)

our_genes = pd.DataFrame({'genes': our_data.var_names.to_numpy()})

with open('scFormer_assets/vocab.json', "r") as f:
    vocab_dict = json.load(f)

with open('scFormer_assets/args.json', "r") as f:
    model_configs = json.load(f)

df = pd.DataFrame(data={'tokens': vocab_dict.keys(),
                        'indices': vocab_dict.values()})

df_intersection = df.merge(our_genes, how='inner', left_on='tokens', right_on='genes')

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
'''
model = TransformerModel(
    ntoken=len(vocab_dict),
    d_model=model_configs["embsize"],
    nhead=model_configs["nheads"],
    d_hid=model_configs["d_hid"],
    nlayers=model_configs["nlayers"],
    nlayers_cls=model_configs["n_layers_cls"],
    n_cls=1,  # TODO: fix loading this
    vocab=vocab_dict,
    dropout=model_configs["dropout"],
    pad_token=model_configs["pad_token"],
    pad_value=model_configs["pad_value"],
)
model.to(device)

try:
    model.load_state_dict(torch.load('scFormer_assets/best_model.pt', map_location=device))
except:
    params = model.state_dict()
    for key, value in torch.load('scFormer_assets/best_model.pt', map_location=device).items():
        # only load params that are in the current model
        if (
                key in model.state_dict()
                and model.state_dict()[key].shape == value.shape
        ):
            params[key] = value
    model.load_state_dict(params)
    
model.eval()
'''

gene_indices = torch.as_tensor(np.array([vocab_dict.get(key) for key in df_intersection.genes.values]))

# gene_embeddings = model.encoder(gene_indices.to(device))

if random_embeddings:
    gene_embeddings = torch.rand(size=(172, 64))
    # num_samples, embedding_dim = 172, 64  # gene_embeddings.shape
    # min_value = -10  # int(gene_embeddings.min().item())
    # max_value = 10  # int(gene_embeddings.max().item())
    # random_values = torch.randint(min_value, max_value, size=(num_samples, embedding_dim), dtype=nn.Embedding.dtype)
    # print(f"Min value: {min_value}, Max value: {max_value}")
    # gene_embeddings = random_values.to(gene_embeddings.device)

final_data = our_data[:, df_intersection.genes].copy()

cell_embeddings = np.matmul(final_data.X, gene_embeddings.detach().cpu().numpy())
cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)


#cell_embeddings = np.matmul(final_data.X > 0, gene_embeddings.detach().cpu().numpy())

cell_embeddings = cell_embeddings / (np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
) + 0.000000000000001)

final_data.obsm["X_scFormer_alejandro"] = cell_embeddings

sc.pp.neighbors(final_data, use_rep="X_scFormer_alejandro")
sc.tl.umap(final_data, min_dist=0.3)

batch_key = 'batch'
if not isinstance(final_data.obs[batch_key][0], str):
    final_data.obs["str_" + batch_key] = final_data.obs[batch_key].astype(str)
    batch_key = "str_" + batch_key

str_celltype_col = 'louvain'
celltype_col = 'int' + str_celltype_col
final_data.obs[celltype_col] = utils.category_str2int(
    final_data.obs[str_celltype_col]
)

fig = sc.pl.umap(
    final_data,
    color=[batch_key, str_celltype_col],
    ncols=2,
    frameon=False,
    return_fig=True,
)
fig.savefig(
    f"random_embeddings_umap[non-weighted-mean]-test-rkh-pancreas_2.png",
    bbox_inches="tight",
)

print('done')
