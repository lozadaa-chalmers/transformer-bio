"""

We'll make a first run of ScFormer using our own data. Specifically '3013_P_filtered_feature_bc_matrix.h5'

"""
# TODO: check if random embeddings are as good as the 'pretrained' ones.
# TODO: normalize our  data to see if it changes the clustering

import json

import numpy as np
import pandas as pd
import scanpy as sc
import torch

from scFormer_assets.model_scformer import TransformerModel
from tasks.data_handling import data_pre_processing as dpp

our_data = dpp.create_count_matrix('data/3013_P_filtered_feature_bc_matrix.h5',
                                   make_genes_unique=True)

our_genes = pd.Series(our_data.var_names.to_numpy(), name='genes')

with open('scFormer_assets/vocab.json', "r") as f:
    vocab_dict = json.load(f)

with open('scFormer_assets/args.json', "r") as f:
    model_configs = json.load(f)

df = pd.DataFrame(data={'tokens': vocab_dict.keys(),
                        'indices': vocab_dict.values()})

df_intersection = df.merge(our_genes, how='inner', left_on='tokens', right_on='genes')

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

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

gene_indices = torch.as_tensor(np.array([vocab_dict.get(key) for key in df_intersection.genes.values]))

gene_embeddings = model.encoder(gene_indices.to(device))

final_data = our_data[:, df_intersection.genes].copy()

cell_embeddings = np.matmul(final_data.X.A, gene_embeddings.detach().cpu().numpy())
cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)

final_data.obsm["X_scFormer_alejandro"] = cell_embeddings

sc.pp.neighbors(final_data, use_rep="X_scFormer_alejandro")
sc.tl.umap(final_data, min_dist=0.3)

fig = sc.pl.umap(
    final_data,
    frameon=False,
    return_fig=True,
)
fig.savefig(
    f"embeddings_umap[weighted-mean]-test-alc-ourdata.png",
    bbox_inches="tight",
)

print('done')
