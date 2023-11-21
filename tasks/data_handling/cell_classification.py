import scanpy as sc
import pandas as pd
from typing import Literal


def classify_cells(pre_processed_adata: sc.AnnData = None,
                   class_df: pd.DataFrame = None,
                   classification_key: str = 'suva_class',
                   save_path: str = 'classified_temp.h5'
                   ) -> sc.AnnData:
    """
    ### Function that classifies cells and saves the result as a h5-file.

    ---

    #### Args:
        - pre_processed_adata (AnnData): AnnData object that has gone through pre-processing. Default None
        - class_df (DataFrame): DataFrame with categories 'genes' and 'set'. Default None
        - classification_key (str): Key to stored classification result in pre_processed_adata.obs.
            Default 'suva_class'
        - save_path (str): Save path for h5-file, include file name. Default 'classified_temp.h5'

    #### Returns:
        - AnnData: With class scores and classification

    ---
    Written: ronjah@chalmers.se
    """

    adata = pre_processed_adata.copy()
    categories = class_df['set'].unique().tolist()

    score_categories(pre_processed_adata=adata,
                     unique_categories_list=categories,
                     class_df=class_df)

    annotate_categories(scored_adata=adata,
                        unique_categories_list=categories,
                        classification_key=classification_key)

    # Save the AnnData object to an h5 file
    adata.write_h5ad(save_path)
    return adata


def score_categories(pre_processed_adata: sc.AnnData = None,
                     unique_categories_list: list = [],
                     class_df: pd.DataFrame = ''
                     ):
    for suva_class in unique_categories_list:
        gene_list = class_df.loc[class_df['set'] == suva_class, 'genes'].tolist()
        sc.tl.score_genes(pre_processed_adata, score_name=suva_class, gene_list=gene_list)
        sc.tl.score_genes(pre_processed_adata, score_name=suva_class, gene_list=gene_list)


def annotate_categories(scored_adata: sc.AnnData = None,
                        unique_categories_list: list = [],
                        classification_key: str = 'suva_class'
                        ):
    # Calculate max value in each cell
    max_class_column = scored_adata.obs[unique_categories_list].idxmax(axis=1)

    # Annotate the result in column 'suva_class'
    scored_adata.obs[classification_key] = max_class_column


def get_suva_dict(dictionary: Literal['full', 'abbr'] = 'abbr') -> dict:
    if dictionary == 'abbr':
        suva_dict = {
            'suva_ac_PMID_31327527': 'AC',
            'suva_g1s_PMID_31327527': 'G1S',
            'suva_g2m_PMID_31327527': 'G2M',
            'suva_mes1_PMID_31327527': 'MES1',
            'suva_mes2_PMID_31327527': 'MES2',
            'suva_npc1_PMID_31327527': 'NPC1',
            'suva_npc2_PMID_31327527': 'NPC2',
            'suva_opc_PMID_31327527': 'OPC'
        }
    else:
        suva_dict = {
            'suva_ac_PMID_31327527': 'astrocyte',
            'suva_g1s_PMID_31327527': 'interphase',
            'suva_g2m_PMID_31327527': 'mitotic-phase',
            'suva_mes1_PMID_31327527': 'mesenchymal1',
            'suva_mes2_PMID_31327527': 'mesenchymal2',
            'suva_npc1_PMID_31327527': 'neural-progenitor1',
            'suva_npc2_PMID_31327527': 'neural-progenitor2',
            'suva_opc_PMID_31327527': 'oligodenrocyte-progenitor'
        }
    return suva_dict
