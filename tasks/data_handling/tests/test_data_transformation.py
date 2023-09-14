import unittest
from data_handling import data_transformation

'''
Test data_handling
'''
class Test_read_data(unittest.TestCase):

    def test_format_check(self):
        # Use assertRaises to check if the function gives an error for the wrong format
        with self.assertRaises(ValueError) as context:
            result = data_transformation.read_data('data_handling/Adam.h5', csr=False)
        
    def test_correct_shape(self):
        # Use assertEqual to check if the function gives the correct shape
        result = data_transformation.read_data('data_handling/Adam.h5', csr=True)
        self.assertEqual([result.shape[0], result.shape[1]], [3660,23797])
    

if __name__ == '__main__':
    unittest.main()
