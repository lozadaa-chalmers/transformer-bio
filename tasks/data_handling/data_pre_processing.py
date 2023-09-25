import scanpy as sc

def create_count_matrix(file_path : str = None) -> sc.AnnData:
    """
    ### Function that transforms h5-file to a count matrix.

    ---
    
    #### Args:
        - file_path (str): Path to a 10x hdf5 file.

    #### Returns: 
        - AnnData: Count matrix from h5-file in file_path with unique variables

    ---

    Annotated data matrix, where observations/cells are named by their barcode and variables/genes by gene name. Stores the following information:

    `~anndata.AnnData.X`
        The data matrix is stored
    `~anndata.AnnData.obs_names`
        Cell names
    `~anndata.AnnData.var_names`
        Gene names
    `~anndata.AnnData.var`\ `['gene_ids']`
        Gene IDs
    `~anndata.AnnData.var`\ `['feature_types']`
        Feature types
    ---

    Written: ronjah@chalmers.se
    """
    if file_path == None:
        raise ValueError('Need a path to h5-file.')
    adata = sc.read_10x_h5(file_path)
    adata.var_names_make_unique()
    return adata


def quality_control(adata: sc.AnnData = None) -> None:
    """
    ### Function that adds QC to AnnData object

    ---
    
    #### Args:
        - adata (AnnData): AnnData object to perform quality control on.

    #### Returns: 
        - None

    ---
    Written: ronjah@chalmers.se
    """
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


def remove_bad_cells(adata: sc.AnnData = None) -> sc.AnnData:
    """
    ### Function that removes bad cells

    ---
    
    #### Args:
        - adata (AnnData): AnnData object that has gone through QC.

    #### Returns: 
        - AnnData: Where bad cells has been removed.

    ---
    Written: ronjah@chalmers.se
    """
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    return adata

# 3: Normalizing the data. LogNormalize?
# Rest is feature selection and dimension reduction

def normalize_data(adata: sc.AnnData) -> None:
    """
    ### Function that normalize the AnnData

    ---
    
    #### Args:
        - adata (AnnData): AnnData object that has gone through QC.

    #### Returns: 
        - None.

    ---
    Written: ronjah@chalmers.se
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def pre_process_data_pipeline(file_path: str=None, plots: bool = False) -> sc.AnnData:
    """
    ### Function that acts as a pipeline for pre-processing the data.

    ---
    
    #### Args:
        - file_path (str): File path to h5-file, default None
        - plots (bool): True if want to show plots, else false. Default False

    #### Returns: 
        - AnnData: Pre-processed

    ---
    ##### Order of operations:
        1. Create count matrix
        2. QC
        3. Remove bad cells
        4. Normalize data

    ---
    Written: ronjah@chalmers.se
    """
    if file_path == None:
        raise ValueError('Need a path to h5-file.')
    
    adata = create_count_matrix(file_path)

    if plots:
        sc.pl.highest_expr_genes(adata, n_top=20)

    quality_control(adata)

    if plots:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
              jitter=0.4, multi_panel=True)
 
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    adata = remove_bad_cells(adata)
    normalize_data(adata)

    if plots:
        sc.pl.highly_variable_genes(adata)

    return adata
