import scanpy as sc
import unittest
import os
import numpy as np
from collections import Counter

import global_settings as gs
from tasks.data_handling import data_pre_processing
from tasks.test_helper.builders.data_handling import data_pre_process as test_data_pre_process
from tasks.test_helper.builders.utils import context_utils as cu


class TestCreateCountMatrix(unittest.TestCase):
    def test_need_file_path(self):
        with self.assertRaises(ValueError) as context:
            result = data_pre_processing.create_count_matrix()

    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                  filename=file_name)
    simulated_adata = data_pre_processing.create_count_matrix(file_path=data_path + '\\' + file_name)

    def test_equal_values_equal_index(self):
        equal = True
        indices = self.simulated_adata.X.nonzero()
        rows = indices[0]
        columns = indices[1]

        for i in range(len(rows)):
            row = rows[i]
            col = columns[i]
            if self.simulated_adata.X[row, col] != self.simulated_matrix.A[row, col]:
                equal = False
        self.assertTrue(equal)

    def test_all_var_names_unique(self):
        name_list = np.array(self.simulated_adata.var_names).tolist()
        names_unique = len(name_list) == len(set(name_list))
        self.assertTrue(names_unique)


class TestQualityControl(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    file_name = 'simulated_csc_h5_data.h5'
    test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                               filename=file_name)
    simulated_adata = data_pre_processing.create_count_matrix(file_path=data_path + '\\' + file_name)
    data_pre_processing.quality_control(simulated_adata)
    n_genes = simulated_adata.n_vars
    n_cells = simulated_adata.n_obs

    def test_n_genes_by_count(self):
        cells = self.simulated_adata.X.nonzero()[0]
        counted_values = Counter(cells)
        equal = True
        for cell, count in counted_values.items():
            if count != self.simulated_adata.obs['n_genes_by_counts'].iloc[cell]:
                equal = False
        self.assertTrue(equal)

    def test_obs_total_counts(self):
        equal = True
        for i in range(self.n_cells):
            if (self.simulated_adata.obs['total_counts'].iloc[i] - sum(self.simulated_adata.X.A[i, :])) > 1e-10:
                equal = False
        self.assertTrue(equal)

    def test_obs_total_counts_mt(self):
        var_mt_indices = np.where(self.simulated_adata.var['mt'].values)[0]
        equal = True
        for cell_index in range(self.n_cells):
            expr = 0
            for mt_index in var_mt_indices:
                expr += self.simulated_adata.X.A[cell_index, mt_index]
            if abs(expr - self.simulated_adata.obs['total_counts_mt'].iloc[cell_index]) > 1e-10:
                equal = False
        self.assertTrue(equal)

    def test_obs_pct_counts_mt(self):
        equal = True
        for cell_index in range(self.n_cells):
            pct = self.simulated_adata.obs['total_counts_mt'].iloc[cell_index] / \
                  self.simulated_adata.obs['total_counts'].iloc[cell_index]
            if abs(100 * pct - self.simulated_adata.obs['pct_counts_mt'].iloc[cell_index]) > 1e-10:
                equal = False
        self.assertTrue(equal)

    def test_var_mt_annotated(self):
        self.assertEqual(self.simulated_adata.var['mt'].all(),
                         self.simulated_adata.var_names.str.startswith('MT').all())

    def test_n_cells_by_count(self):
        genes = self.simulated_adata.X.nonzero()[1]
        counted_values = Counter(genes)
        equal = True
        for gene, count in counted_values.items():
            if count != self.simulated_adata.var['n_cells_by_counts'].iloc[gene]:
                equal = False
        self.assertTrue(equal)

    def test_mean_counts(self):
        equal = True
        for iGene in range(self.simulated_adata.n_vars):
            mean = self.simulated_adata.var['total_counts'].iloc[iGene] / self.n_cells
            if abs(mean - self.simulated_adata.var['mean_counts'].iloc[iGene]) > 1e-10:
                equal = False
        self.assertTrue(equal)

    def test_pct_dropout_by_count(self):
        genes = self.simulated_adata.X.nonzero()[1]
        counted_values = Counter(genes)
        equal = True
        for gene, count in counted_values.items():
            pct = 1 - count / self.n_cells
            if abs(100 * pct - self.simulated_adata.var['pct_dropout_by_counts'].iloc[gene]) > 1e-10:
                equal = False
        self.assertTrue(equal)

    def test_var_total_counts(self):
        equal = True
        for gene_index in range(self.n_genes):
            if abs(self.simulated_adata.var['total_counts'].iloc[gene_index] - sum(
                    self.simulated_adata.X.A[:, gene_index])) > 1e-10:
                equal = False
        self.assertTrue(equal)


class TestRemoveBadCells(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    file_name = 'simulated_csc_h5_data.h5'
    test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                               filename=file_name)
    simulated_adata = data_pre_processing.create_count_matrix(file_path=data_path + '\\' + file_name)
    data_pre_processing.quality_control(simulated_adata)

    def test_too_high_gene_count_removed(self):
        max_n_genes = int(self.simulated_adata.n_vars * 0.1)
        wanted_cells = self.simulated_adata.obs.n_genes_by_counts < max_n_genes
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.simulated_adata,
                                                          max_n_genes=max_n_genes,
                                                          min_n_genes=0,
                                                          mitochondrial_percent=100)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))

    def test_too_low_gene_count_removed(self):
        min_n_genes = int(self.simulated_adata.n_vars * 0.1)
        wanted_cells = self.simulated_adata.obs.n_genes_by_counts > min_n_genes
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.simulated_adata,
                                                          max_n_genes=1000,
                                                          min_n_genes=min_n_genes,
                                                          mitochondrial_percent=100)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))

    def test_too_high_mt_pct_removed(self):
        mt_pct = 5
        wanted_cells = self.simulated_adata.obs.pct_counts_mt < mt_pct
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.simulated_adata,
                                                          max_n_genes=1000,
                                                          min_n_genes=0,
                                                          mitochondrial_percent=mt_pct)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))


# Test normalize data
class TestNormalizeData(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    file_name = 'simulated_csc_h5_data.h5'
    test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                               filename=file_name)
    simulated_adata = data_pre_processing.create_count_matrix(file_path=data_path + '\\' + file_name)
    data_pre_processing.quality_control(simulated_adata)

    max_n_genes = int(simulated_adata.n_vars * 0.2)
    min_n_genes = int(simulated_adata.n_vars * 0.1)
    mt_pct = 5
    target_sum = 1e4

    simulated_adata = data_pre_processing.remove_bad_cells(simulated_adata,
                                                           max_n_genes=max_n_genes,
                                                           min_n_genes=min_n_genes,
                                                           mitochondrial_percent=mt_pct)
    simulated_adata_copy = simulated_adata.copy()
    data_pre_processing.normalize_data(simulated_adata,
                                       target_sum=target_sum,
                                       exclude_highly_expressed=False)

    sc.pp.normalize_total(simulated_adata_copy,
                          target_sum=target_sum,
                          exclude_highly_expressed=False)

    sc.pp.log1p(simulated_adata_copy)

    sc.pp.highly_variable_genes(simulated_adata_copy)

    def test_normalize_data_right_order(self):
        equal = True
        indices = self.simulated_adata.X.nonzero()
        rows = indices[0]
        columns = indices[1]

        for i in range(len(rows)):
            row = rows[i]
            col = columns[i]
            if self.simulated_adata.X[row, col] != self.simulated_adata_copy.X[row, col]:
                equal = False
        self.assertTrue(equal)


class TestPreProcessPipeline(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    file_name = 'simulated_csc_h5_data.h5'
    test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                               filename=file_name)
    simulated_adata = data_pre_processing.create_count_matrix(file_path=data_path + '\\' + file_name)
    data_pre_processing.quality_control(simulated_adata)

    max_n_genes = int(simulated_adata.n_vars * 0.2)
    min_n_genes = int(simulated_adata.n_vars * 0.1)
    mt_pct = 5
    target_sum = 1e4

    simulated_adata = data_pre_processing.remove_bad_cells(simulated_adata,
                                                           max_n_genes=max_n_genes,
                                                           min_n_genes=min_n_genes,
                                                           mitochondrial_percent=mt_pct)
    data_pre_processing.normalize_data(simulated_adata,
                                       target_sum=target_sum,
                                       exclude_highly_expressed=False)
    pipeline_adata = data_pre_processing.pre_process_data_pipeline(file_path=data_path + '\\' + file_name,
                                                                   max_n_genes=max_n_genes,
                                                                   min_n_genes=min_n_genes,
                                                                   mitochondrial_percent=mt_pct,
                                                                   target_sum=target_sum)

    def test_annData_X_equal(self):
        equal = True
        indices = self.simulated_adata.X.nonzero()
        rows = indices[0]
        columns = indices[1]

        for i in range(len(rows)):
            row = rows[i]
            col = columns[i]
            if self.simulated_adata.X[row, col] != self.pipeline_adata.X[row, col]:
                equal = False
        self.assertTrue(equal)


if __name__ == '__main__':
    unittest.main()
