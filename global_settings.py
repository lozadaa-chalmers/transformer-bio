TEST_DATA_PATH = 'test_helper\\builders\\data'

# Uppsala's data processing parameters
PROCESSING_PARAMETERS = {
    'file_path':  None,
    'plots': False,
    'remove_cells': True,
    'remove_genes': True,
    'normalize': True,
    'max_n_genes': 1e6,
    'min_n_genes': 500,
    'mitochondrial_percent': 20.0,
    'target_sum': 1e4,
    'exclude_highly_expressed': False,
    'max_fraction': 0.05,
    'min_mean': 0.0125,
    'max_mean': 3,
    'min_dispersion': 0.5,
    'n_top_genes': None,
    'min_n_cells': 50,
    'n_top_plot': 20,
    'save_format': 'png'
}
