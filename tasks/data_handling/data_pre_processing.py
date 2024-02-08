import warnings

import scanpy as sc


def create_count_matrix(
        file_path: str = None,
        make_genes_unique: bool = False
) -> sc.AnnData:
    """
    ### Function that transforms 10x h5-file to a count matrix.

    ---
    
    #### Args:
        - file_path (str): Path to a 10x hdf5 file. Default None.
        - make_genes_unique (bool): Boolean indicating if genes should be made unique. If True, '-i' suffix is added
          to the variable name, where i is an integer. Default False.

    #### Returns: 
        - AnnData: Count matrix from h5-file in file_path with unique variables.

    ---
    Annotated data matrix, where observations/cells are named by their barcode and variables/genes by gene name.
    Stores the following information:

    `~anndata.AnnData.X`
        The sparse data matrix is stored
    `~anndata.AnnData.obs_names`
        Cell names
    `~anndata.AnnData.var_names`
        Gene names
    `~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    `~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    ---

    Written: ronjah@chalmers.se
    """
    if file_path is None:
        raise ValueError('Need a path to h5-file.')

    if make_genes_unique:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    message="Variable names are not unique. "
                                            "To make them unique, call `.var_names_make_unique`.")
            adata = sc.read_10x_h5(file_path)
            adata.var_names_make_unique()
    else:
        adata = sc.read_10x_h5(file_path)

    return adata


def quality_control(
        adata: sc.AnnData = None
):
    """
    ### Function that adds QC to AnnData object

    ---
    
    #### Args:
        - adata (AnnData): AnnData object to perform quality control on.

    ---
    #### Adds keys:
        - `adata.obs['n_genes_by_counts']` Number of active genes per cell 
        - `adata.obs['total_counts']` Sum of all gene expressions in each cell
        - `adata.obs['total_counts_mt']` Sum of all mitochondrial expressions in each cell
        - `adata.obs['pct_counts_mt']` Percent of total expression in each cell that are from mitochondrial expression
        - `adata.var['mt']` Booleans of which variables are mitochondrial genes
        - `adata.var['n_cells_by_counts']` Number of cells each gene is active
        - `adata.var['total_counts']` Sum of expression of each gene over all cells
        - `adata.var['mean_counts']` Average expression of each gene over all cells
        - `adata.var['pct_dropout_by_counts']` Percentage of cells each gene has no expression in
        
    More information: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html

    ---
    Written: ronjah@chalmers.se
    """
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt'],
                               percent_top=None,
                               log1p=False,
                               inplace=True)


def remove_bad_cells(
        adata: sc.AnnData = None,
        max_n_genes: int = 2500,
        min_n_genes: int = 200,
        mitochondrial_percent: float = 5.0
) -> sc.AnnData:
    """
    ### Function that removes bad cells

    ---
    
    #### Args:
        - adata (AnnData): AnnData object that has gone through quality control (QC). Default None
        - max_n_genes (int): Remove cells with a gene count over max_n_genes. Default 2500
        - min_n_genes (int): Remove cells with a gene count less than min_n_genes. Default 200
        - mitochondrial_percent (float): Remove cells with a mitochondrial expression percent above
            mitochondrial_percent. Default 5.0

    #### Returns: 
        - AnnData: Where bad cells has been removed.

    ---
    Written: ronjah@chalmers.se
    """
    filtered_adata = adata[(adata.obs.n_genes_by_counts < max_n_genes) &
                           (adata.obs.n_genes_by_counts > min_n_genes) &
                           (adata.obs.pct_counts_mt < mitochondrial_percent), :]
    return filtered_adata


def normalize_data(
        adata: sc.AnnData,
        target_sum: float = 1e4,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
):
    """
    ### Function that normalize the AnnData by cell and logarithmize each value by Log(X+1) with the natural logarithm.

    ---
    
    #### Args:
        - adata (AnnData): AnnData object that has gone through quality control (QC).
        - target_sum (float): Float each cell sums up to. Default 1e4
        - exclude_highly_expressed (bool): If true, highly expressed cells are not part of normalization. Default False
        - max_fraction (float): If `exclude_highly_expressed = True`, consider cells as highly expressed that have more
            counts than max_fraction of the original total counts in at least one cell. Default 0.05

    ---
    Written: ronjah@chalmers.se
    """
    sc.pp.normalize_total(adata,
                          target_sum=target_sum,
                          exclude_highly_expressed=exclude_highly_expressed,
                          max_fraction=max_fraction,
                          inplace=True)
    sc.pp.log1p(adata,
                copy=False)


def filter_genes(
        adata: sc.AnnData,
        min_mean: float = 0.0125,
        max_mean: float = 3,
        min_dispersion: float = 0.5,
        n_top_genes: int = None,
        min_n_cells: int = 50
) -> sc.AnnData:
    """
    ### Function that filters out genes from AnnData.

    ---
    
    #### Args:
        - adata (AnnData): AnnData object that has gone through normalization.
        - min_mean (float): Minimum cutoff for means, only used if n_top_genes unequals None. Default 0.0125
        - max_mean (float): Maximum cutoff for means, only used if n_top_genes unequals None. Default 3
        - min_dispersion (float): The threshold for the minimum dispersion a gene must have to be considered highly
            variable, only used if n_top_genes unequals None. Default 0.5
        - n_top_genes (int): Number of highly-variable genes to keep. Default None
        - min_n_cells (int): Number of minimum cells each gene must have an expression in to be kept. Default 50

    ### Returns:
        - AnnData: AnnData object which have filtered out genes and labeled them according to the parameters given.

    ---
    Written: ronjah@chalmers.se
    """

    filtered_adata = adata[:, (adata.var.n_cells_by_counts > min_n_cells)]
    sc.pp.highly_variable_genes(filtered_adata,
                                min_mean=min_mean,
                                max_mean=max_mean,
                                min_disp=min_dispersion,
                                n_top_genes=n_top_genes,
                                inplace=True)
    return filtered_adata


def pre_process_data_pipeline(
        file_path: str = None,
        plots: bool = False,
        max_n_genes: int = 2500,
        min_n_genes: int = 200,
        mitochondrial_percent: float = 5.0,
        target_sum: float = 1e4,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
        min_mean: float = 0.0125,
        max_mean: float = 3,
        min_dispersion: float = 0.5,
        n_top_genes: int = None,
        min_n_cells: int = 50,
        n_top_plot: int = 20,
        save_format: str = 'png'
) -> sc.AnnData:
    """
    ### Function that acts as a pipeline for pre-processing the data.

    ---
    
    #### Args:
        - file_path (str): File path to h5-file, default None
        - plots (bool): True if wanting to save plots, else false. Default False
        - max_n_genes (int): Remove cells with a gene count over max_n_genes. Default 2500
        - min_n_genes (int): Remove cells with a gene count less than min_n_genes. Default 200
        - mitochondrial_percent (float): Remove cells with a mitochondrial percent above mitochondrial_percent.
            Default 5.0
        - target_sum (float): Float each cell sums up to. Default 1e4
        - exclude_highly_expressed (bool): If true, highly expressed cells are not part of normalization. Default False
        - max_fraction (float): If `exclude_highly_expressed = True`, consider cells as highly expressed that have more
            counts than max_fraction of the original total counts in at least one cell. Default 0.05
        - min_mean (float): Minimum cutoff for means, only used if n_top_genes unequals None. Default 0.0125
        - max_mean (float): Maximum cutoff for means, only used if n_top_genes unequals None. Default 3
        - min_dispersion (float): The threshold for the minimum dispersion a gene must have to be considered highly
            variable, only used if n_top_genes unequals None. Default 0.5
        - n_top_genes (int): Number of highly-variable genes to keep. Default None
        - min_n_cells (int): Number of minimum cells each gene must have an expression in to be kept. Default 50
        - n_top_plot (int): The n_top genes with the highest mean fraction over all cells are plotted as boxplot.
            Default 20
        - save_format (str): File format for saved plots. Allowed formats: ['png', 'pdf', 'svg']. Default 'png'

    #### Returns: 
        - AnnData: Pre-processed

    ---
    ##### Order of operations:
        1. Create count matrix
        2. Quality control (QC)
        3. Remove bad cells
        4. Normalize data
        5. Filter out genes and mark if they are highly variable

    ---
    Written: ronjah@chalmers.se
    """
    if file_path is None:
        raise ValueError('Need a path to h5-file.')

    adata = create_count_matrix(file_path)

    if plots:
        allowed_formats = ['png', 'pdf', 'svg']
        if save_format not in allowed_formats:
            raise ValueError(f"Invalid save_format: '{save_format}'. It must be one of {allowed_formats}")

        sc.pl.highest_expr_genes(adata,
                                 n_top=n_top_plot,
                                 save='.' + save_format,
                                 show=False)

    quality_control(adata)

    if plots:
        sc.pl.violin(adata,
                     ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                     jitter=0.4,
                     multi_panel=True,
                     save='.' + save_format,
                     show=False)

        sc.pl.scatter(adata,
                      x='total_counts',
                      y='pct_counts_mt',
                      save='_pct_counts_mt' + '.' + save_format,
                      show=False)

        sc.pl.scatter(adata,
                      x='total_counts',
                      y='n_genes_by_counts',
                      save='_n_genes_by_count' + '.' + save_format,
                      show=False)

    adata_filtered_cells = remove_bad_cells(adata,
                                            max_n_genes=max_n_genes,
                                            min_n_genes=min_n_genes,
                                            mitochondrial_percent=mitochondrial_percent)
    normalize_data(adata_filtered_cells,
                   target_sum=target_sum,
                   exclude_highly_expressed=exclude_highly_expressed,
                   max_fraction=max_fraction)

    adata_filtered_genes = filter_genes(adata=adata_filtered_cells,
                                        min_mean=min_mean,
                                        max_mean=max_mean,
                                        min_dispersion=min_dispersion,
                                        n_top_genes=n_top_genes,
                                        min_n_cells=min_n_cells)

    if plots:
        sc.pl.highly_variable_genes(adata_filtered_genes,
                                    save='.' + save_format,
                                    show=False)

    return adata_filtered_genes
