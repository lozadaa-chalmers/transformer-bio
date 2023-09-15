import os
import unittest

import global_settings as gs
from tasks.data_handling.data_transformation import read_data
from tasks.test_helper.builders.data_handling import data_transformation as test_data_transformation
from tasks.test_helper.builders.utils import context_utils as cu

'''
Test data_handling
'''


class TestReadData(unittest.TestCase):

    def test_matrix_reconstruction(self, mode='csr'):
        data_path = os.path.join('\\'.join(os.getcwd().split('\\')[:-2]), gs.TEST_DATA_PATH)
        cu.clean(folder_path=data_path)
        simulated_matrix = test_data_transformation.simulate_h5_data(mode=mode,
                                                                     file_path=data_path)

        reconstructed_matrix = read_data(filename=os.path.join(data_path, 'simulated_h5_data.h5'),
                                         mode=mode)
        with self.assertRaises(ValueError):
            self.assertEqual(simulated_matrix.toarray(), reconstructed_matrix.toarray())


if __name__ == '__main__':
    unittest.main()
