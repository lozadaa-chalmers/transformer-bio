import os
import unittest
from collections import Counter

import numpy as np

import global_settings as gs
from tasks.data_handling import data_pre_processing
from tasks.test_helper.builders.data_handling import data_pre_process as test_data_pre_process
from tasks.test_helper.builders.utils import context_utils as cu


class TestCreateCountMatrix(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix, h5_path = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                           filename=file_name)
    reconstructed_adata_unique = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                         make_genes_unique=True)
    reconstructed_adata_duplicates = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                             make_genes_unique=False)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    def test_need_file_path(self):
        with self.assertRaises(ValueError) as context:
            result = data_pre_processing.create_count_matrix()

    def test_equal_values_equal_index(self):
        self.assertTrue((self.reconstructed_adata_unique.X.toarray() == self.simulated_matrix.A).all())

    def test_all_var_names_unique(self):
        name_list = np.array(self.reconstructed_adata_unique.var_names).tolist()
        self.assertTrue(len(name_list) == len(set(name_list)))

    def test_var_names_duplicates(self):
        name_list = np.array(self.reconstructed_adata_duplicates.var_names).tolist()
        self.assertTrue(len(name_list) != len(set(name_list)))


class TestQualityControl(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix, h5_path = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                           filename=file_name)
    reconstructed_adata = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                  make_genes_unique=True)
    data_pre_processing.quality_control(reconstructed_adata)

    n_genes = reconstructed_adata.n_vars
    n_cells = reconstructed_adata.n_obs

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    def test_n_genes_by_count(self):
        cells = self.reconstructed_adata.X.nonzero()[0]
        counted_values = Counter(cells)
        self.assertEqual(list(self.reconstructed_adata.obs.n_genes_by_counts),
                         list(counted_values.values()))

    def test_obs_total_counts(self):
        self.assertEqual(list(self.reconstructed_adata.obs.total_counts),
                         list(np.sum(self.reconstructed_adata.X, axis=1)))

    def test_obs_total_counts_mt(self):
        true_indices = np.where(self.reconstructed_adata.var.mt)[0]
        diff_array = abs(self.reconstructed_adata.obs.total_counts_mt.values -
                         self.reconstructed_adata.X.A[:, true_indices].sum(axis=1))
        self.assertTrue(all(x < 1e-5 for x in diff_array))

    def test_obs_pct_counts_mt(self):
        pct = 100 * (self.reconstructed_adata.obs.total_counts_mt / self.reconstructed_adata.obs.total_counts)
        self.assertEqual(list(pct),
                         list(self.reconstructed_adata.obs.pct_counts_mt))

    def test_var_mt_annotated(self):
        self.assertEqual(self.reconstructed_adata.var['mt'].all(),
                         self.reconstructed_adata.var_names.str.startswith('MT').all())

    def test_n_cells_by_count(self):
        genes = self.reconstructed_adata.X.nonzero()[1]
        counted_values = Counter(genes)
        counted_list = sorted(counted_values.items())
        n_counts = [tup[1] for tup in counted_list]
        self.assertEqual(self.reconstructed_adata.var.n_cells_by_counts.to_list(),
                         n_counts)

    def test_mean_counts(self):
        mean = self.reconstructed_adata.var.total_counts / self.n_cells
        self.assertEqual(list(self.reconstructed_adata.var.mean_counts),
                         list(mean))

    def test_pct_dropout_by_count(self):
        genes = self.reconstructed_adata.X.nonzero()[1]
        counter = Counter(genes)
        counted_list = sorted(counter.items())
        n_counts = [tup[1] for tup in counted_list]
        dropout_pct = 100 * (1 - np.array(n_counts) / self.n_cells)
        self.assertEqual(list(self.reconstructed_adata.var.pct_dropout_by_counts),
                         list(dropout_pct))

    def test_var_total_counts(self):
        self.assertEqual(self.reconstructed_adata.var.total_counts.to_list(),
                         list(self.reconstructed_adata.X.A.sum(axis=0)))


class TestRemoveBadCells(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix, h5_path = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                           filename=file_name)
    reconstructed_adata = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                  make_genes_unique=True)
    data_pre_processing.quality_control(reconstructed_adata)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    def test_too_high_gene_count_removed(self):
        max_n_genes = int(self.reconstructed_adata.n_vars * 0.1)
        wanted_cells = self.reconstructed_adata.obs.n_genes_by_counts < max_n_genes
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.reconstructed_adata,
                                                          max_n_genes=max_n_genes,
                                                          min_n_genes=0,
                                                          mitochondrial_percent=100)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))

    def test_too_low_gene_count_removed(self):
        min_n_genes = int(self.reconstructed_adata.n_vars * 0.1)
        wanted_cells = self.reconstructed_adata.obs.n_genes_by_counts > min_n_genes
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.reconstructed_adata,
                                                          max_n_genes=1000,
                                                          min_n_genes=min_n_genes,
                                                          mitochondrial_percent=100)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))

    def test_too_high_mt_pct_removed(self):
        mt_pct = 5
        wanted_cells = self.reconstructed_adata.obs.pct_counts_mt < mt_pct
        wanted_cells_indices = np.where(wanted_cells.values)[0]
        temp_adata = data_pre_processing.remove_bad_cells(self.reconstructed_adata,
                                                          max_n_genes=1000,
                                                          min_n_genes=0,
                                                          mitochondrial_percent=mt_pct)
        wanted_cells_array = np.array(wanted_cells.iloc[wanted_cells_indices].keys())
        cells_left_adata_array = np.array(temp_adata.obs_names)
        self.assertTrue(np.array_equal(cells_left_adata_array,
                                       wanted_cells_array))


class TestNormalizeData(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix, h5_path = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                           filename=file_name)
    reconstructed_adata = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                  make_genes_unique=True)
    data_pre_processing.quality_control(reconstructed_adata)

    max_n_genes = int(reconstructed_adata.n_vars * 0.2)
    min_n_genes = int(reconstructed_adata.n_vars * 0.1)
    mt_pct = 5
    target_sum = 1e4

    reconstructed_adata = data_pre_processing.remove_bad_cells(reconstructed_adata,
                                                               max_n_genes=max_n_genes,
                                                               min_n_genes=min_n_genes,
                                                               mitochondrial_percent=mt_pct)

    data_pre_processing.normalize_data(reconstructed_adata,
                                       target_sum=target_sum,
                                       exclude_highly_expressed=False)
    reconstructed_adata_e = np.exp(reconstructed_adata.X.toarray())
    row_sums = np.sum(reconstructed_adata_e, axis=1)

    diff = row_sums - row_sums[0]

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    def test_normalize_data_correct(self):
        self.assertTrue(np.all(abs(self.diff) < 1e-3))


class TestFilterGenes(unittest.TestCase):
    data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
    file_name = 'simulated_csc_h5_data.h5'
    simulated_matrix, h5_path = test_data_pre_process.simulate_csc_h5_data(file_path=data_path,
                                                                           filename=file_name)
    reconstructed_adata = data_pre_processing.create_count_matrix(file_path=h5_path,
                                                                  make_genes_unique=True)
    data_pre_processing.quality_control(reconstructed_adata)

    max_n_genes = int(reconstructed_adata.n_vars * 0.2)
    min_n_genes = int(reconstructed_adata.n_vars * 0.1)
    mt_pct = 5
    target_sum = 1e4

    min_mean = 0.0125
    max_mean = 3
    min_dispersion = 0
    min_n_cells = 5

    reconstructed_adata = data_pre_processing.remove_bad_cells(reconstructed_adata,
                                                               max_n_genes=max_n_genes,
                                                               min_n_genes=min_n_genes,
                                                               mitochondrial_percent=mt_pct)

    data_pre_processing.normalize_data(reconstructed_adata,
                                       target_sum=target_sum,
                                       exclude_highly_expressed=False)

    comparison_adata = reconstructed_adata.copy()
    comparison_adata = comparison_adata[:, (comparison_adata.var.n_cells_by_counts > min_n_cells) &
                                           (comparison_adata.var.mean_counts > min_mean) &
                                           (comparison_adata.var.mean_counts < max_mean)]

    reconstructed_adata_genes_filtered = data_pre_processing.filter_genes(reconstructed_adata,
                                                                          min_mean=min_mean,
                                                                          max_mean=max_mean,
                                                                          min_dispersion=min_dispersion,
                                                                          min_n_cells=min_n_cells)

    def tearDown(self) -> None:
        cu.clean(folder_path=self.data_path)

    def test_correct_genes_removed(self):
        self.assertTrue(
            (self.comparison_adata.X.toarray() == self.reconstructed_adata_genes_filtered.X.toarray()).all())


if __name__ == '__main__':
    unittest.main()
