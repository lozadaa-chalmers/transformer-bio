import scanpy as sc
import pandas as pd
from typing import Literal


def classify_cells(pre_processed_adata: sc.AnnData = None,
                   class_df: pd.DataFrame = None,
                   classification_key: str = 'suva_class',
                   save_path: str = 'classified_temp.h5',
                   store_classification_score: bool = True
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
        - store_classification_score (bool): If True the scores of each class is included in the h5-file, if False
            only the classification is included without the scores. Default True

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

    if not store_classification_score:
        adata_temp = pre_processed_adata.copy()
        adata_temp.obs[classification_key] = adata.obs[classification_key]
        adata = adata_temp

    # Save the AnnData object to a h5 file
    adata.write_h5ad(save_path)
    return adata


def score_categories(pre_processed_adata: sc.AnnData = None,
                     unique_categories_list: list = None,
                     class_df: pd.DataFrame = None
                     ):
    """
    ### Function that scores categories in the AnnData object

    ---

    #### Args:
        - pre_processed_adata (AnnData): AnnData object that has gone through pre-processing. Default None
        - unique_categories_list (list): List with the unique categories. Default None
        - class_df (DataFrame): DataFrame with categories 'genes' and 'set'. Default None

    ---
    Written: ronjah@chalmers.se
    """
    for suva_class in unique_categories_list:
        gene_list = class_df.loc[class_df['set'] == suva_class, 'genes'].tolist()
        sc.tl.score_genes(pre_processed_adata, score_name=suva_class, gene_list=gene_list)
        sc.tl.score_genes(pre_processed_adata, score_name=suva_class, gene_list=gene_list)


def annotate_categories(scored_adata: sc.AnnData = None,
                        unique_categories_list: list = None,
                        classification_key: str = 'suva_class'
                        ):
    """
    ### Function that annotates the class with the highest score for each cell

    ---

    #### Args:
        - scored_adata (AnnData): AnnData object that has scored classes. Default None
        - unique_categories_list (list): List with the unique categories. Default None
        - classification_key (str): Key to stored classification result in pre_processed_adata.obs.
            Default 'suva_class'

    ---
    Written: ronjah@chalmers.se
    """
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
