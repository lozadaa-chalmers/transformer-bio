import unittest

import numpy as np
import pandas as pd
import scanpy as sc

from tasks.data_handling import cell_classification


class TestScoreCategories(unittest.TestCase):
    unique_categories_list = ['c1', 'c2', 'c3']
    d = {'set': ['c1', 'c1', 'c1', 'c2', 'c2', 'c2', 'c3', 'c3', 'c3'],
         'genes': ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9']}
    class_df = pd.DataFrame(data=d)
    X_data = np.random.rand(3, 9)

    adata = sc.AnnData(X=X_data)
    adata.obs_names = ['a', 'b', 'c']
    adata.var_names = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9']
    adata_scored = cell_classification.score_categories(pre_processed_adata=adata,
                                                        unique_categories_list=unique_categories_list,
                                                        class_df=class_df)
    c1_test = np.zeros((3, 1))
    c2_test = np.zeros((3, 1))
    c3_test = np.zeros((3, 1))
    for cell in range(3):
        # breakpoint()
        c1_test[cell] = (sum(X_data[cell, 0:3]) / 3 - sum(X_data[cell, 3:]) / 6)
        c2_test[cell] = (sum(X_data[cell, 3:6]) / 3 - sum(X_data[cell, 0:3] + X_data[cell, 6:]) / 6)
        c3_test[cell] = (sum(X_data[cell, 6:9]) / 3 - sum(X_data[cell, 0:6]) / 6)

    def test_c1_correct(self):
        self.assertTrue(
            (np.round(np.array(self.adata_scored.obs.c1).reshape((3, 1)), 7) == np.round(self.c1_test, 7)).all())

    def test_c2_correct(self):
        self.assertTrue(
            (np.round(np.array(self.adata_scored.obs.c2).reshape((3, 1)), 7) == np.round(self.c2_test, 7)).all())

    def test_c3_correct(self):
        self.assertTrue(
            (np.round(np.array(self.adata_scored.obs.c3).reshape((3, 1)), 7) == np.round(self.c3_test, 7)).all())


class TestAnnotateCategories(unittest.TestCase):
    unique_categories_list = ['c1', 'c2', 'c3']
    d = {'set': ['c1', 'c1', 'c1', 'c2', 'c2', 'c2', 'c3', 'c3', 'c3'],
         'genes': ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9']}
    class_df = pd.DataFrame(data=d)
    X_data = np.random.rand(3, 9)
    adata = sc.AnnData(X=X_data)
    adata.obs_names = ['a', 'b', 'c']
    adata.var_names = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9']
    adata_scored = cell_classification.score_categories(pre_processed_adata=adata,
                                                        unique_categories_list=unique_categories_list,
                                                        class_df=class_df)
    adata_classified = cell_classification.annotate_categories(scored_adata=adata_scored,
                                                               unique_categories_list=unique_categories_list)

    max_class_column = adata_classified.obs[unique_categories_list].idxmax(axis=1)

    def test_classes_correct(self):
        self.assertTrue((self.adata_classified.obs.suva_class == self.max_class_column).all())
