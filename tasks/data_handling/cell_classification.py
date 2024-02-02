import scanpy as sc
import pandas as pd
from typing import Literal


def score_categories(pre_processed_adata: sc.AnnData = None,
                     unique_categories_list: list = None,
                     class_df: pd.DataFrame = None
                     ) -> sc.AnnData:
    """
    ### Function that scores categories in the AnnData object

    ---

    #### Args:
        - pre_processed_adata (AnnData): AnnData object that has gone through pre-processing via
            data_pre_processing.pre_process_data_pipeline(). Default None
        - unique_categories_list (list): List with the unique categories. Default None
        - class_df (DataFrame): DataFrame with categories 'genes' and 'set'. Default None

    #### Return:
        - AnnData: With class scores

    ---
    Written: ronjah@chalmers.se
    """

    temp_adata = pre_processed_adata.copy()

    for suva_class in unique_categories_list:
        gene_list = class_df.loc[class_df['set'] == suva_class, 'genes'].tolist()
        sc.tl.score_genes(temp_adata, score_name=suva_class, gene_list=gene_list)

    return temp_adata


def annotate_categories(scored_adata: sc.AnnData = None,
                        unique_categories_list: list = None,
                        classification_key: str = 'suva_class'
                        ) -> sc.AnnData:
    """
    ### Function that annotates the class with the highest score for each cell

    ---

    #### Args:
        - scored_adata (AnnData): AnnData object that has scored classes. Default None
        - unique_categories_list (list): List with the unique categories. Default None
        - classification_key (str): Key to stored classification result in pre_processed_adata.obs.
            Default 'suva_class'

    #### Return:
        - AnnData: With class scores and classification

    ---
    Written: ronjah@chalmers.se
    """
    temp_adata = scored_adata.copy()

    # Calculate max value in each cell
    max_class_column = temp_adata.obs[unique_categories_list].idxmax(axis=1)

    # Annotate the result in column 'classification_key'
    temp_adata.obs[classification_key] = max_class_column

    return temp_adata


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
            'suva_g1s_PMID_31327527': 'G1/S-phase',
            'suva_g2m_PMID_31327527': 'G2/M-phase',
            'suva_mes1_PMID_31327527': 'mesenchymal1',
            'suva_mes2_PMID_31327527': 'mesenchymal2',
            'suva_npc1_PMID_31327527': 'neural-progenitor1',
            'suva_npc2_PMID_31327527': 'neural-progenitor2',
            'suva_opc_PMID_31327527': 'oligodenrocyte-progenitor'
        }
    return suva_dict


def classify_cells(pre_processed_adata: sc.AnnData = None,
                   class_df: pd.DataFrame = None,
                   classification_key: str = 'suva_class',
                   save_path: str = None,
                   store_classification_score: bool = True,
                   use_dictionary: bool = False,
                   dictionary_type: Literal['full', 'abbr'] = 'abbr'
                   ) -> sc.AnnData:
    """
    ### Function that classifies cells with Scanpy's version of Seurat's AddModuleScore
        and saves the result as a h5-file.

        More information:
            - https://satijalab.org/seurat/reference/addmodulescore
            - https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html

    ---

    #### Args:
        - pre_processed_adata (AnnData): AnnData object that has gone through pre-processing via
            data_pre_processing.pre_process_data_pipeline(). Default None
        - class_df (DataFrame): DataFrame with categories 'genes' and 'set'. Default None
        - classification_key (str): Key to stored classification result in pre_processed_adata.obs.
            Default 'suva_class'
        - save_path (str): Save path for h5-file, include file name. Default None.
        - store_classification_score (bool): If True the scores of each class is included in the h5-file, if False
            only the classification is included without the scores. Default True
        - use_dictionary (bool): Make True if wanting to change the names to a more readable format,
            instead of the raw names from text file. Default False
        - dictionary_type (Literal['full', 'abbr']): Only applied if use_dictionary=True.
            Determines if full names of cell states are used or its abbreviations. Default 'abbr'

    #### Returns:
        - AnnData: With class scores and classification

    ---
    Written: ronjah@chalmers.se
    """

    adata = pre_processed_adata.copy()
    categories = class_df['set'].unique().tolist()

    adata_scored = score_categories(pre_processed_adata=adata,
                                    unique_categories_list=categories,
                                    class_df=class_df)

    adata_classified = annotate_categories(scored_adata=adata_scored,
                                           unique_categories_list=categories,
                                           classification_key=classification_key)

    if use_dictionary:
        class_dict = get_suva_dict(dictionary_type)
        adata_classified.obs[classification_key] = adata_classified.obs[classification_key].replace(class_dict)

    if save_path is not None:
        if store_classification_score:
            # Save the AnnData object to a h5 file
            adata_classified.write_h5ad(save_path)
        else:
            adata_temp = adata_classified.copy()
            adata_temp.obs.drop(columns=categories, inplace=True)
            adata_temp.write_h5ad(save_path)

    return adata_classified
