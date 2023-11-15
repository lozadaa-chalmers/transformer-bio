"""

We'll make a first run of ScFormer using our own data. Specifically '3013_P_filtered_feature_bc_matrix.h5'

"""
import json

from tasks.data_handling import data_pre_processing as dpp

# TODO: look into what is being done with 'sc.pp.neighbors' and 'sc.tl.umap'
# TODO: check what's going on with read_10X 'genome' option.
# TODO: check if random embeddings are as good as the 'pretrained' ones.
# TODO: check if there's any difference between cell embedding strategies.
# TODO: normalize our  data to see if it changes the clustering.

random_embeddings = True

our_data = dpp.create_count_matrix('data/pancreas.h5ad',
                                   make_genes_unique=True)
# get top expressed genes in vocab
with open('scFormer_assets/vocab.json', "r") as f:
    vocab_dict = json.load(f)
breakpoint()
j = 0
for gene in vocab_dict.keys():
    if any(our_data.var.index == gene):

breakpoint()
