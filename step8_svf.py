import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import logging

logger = logging.getLogger(__name__)

def run_step8_svf(adata, output_path=None, sample_tag="sample"):
    """
    Step 8: Identify Spatially Variable Features (SVF) using SpatialDE.
    References: notebooks/8_SVF_identification/8_1_batch_processing_SpatialDE_SVF.ipynb
    """
    logger.info("="*60)
    logger.info("STEP 8: SVF Identification (SpatialDE)")
    logger.info("="*60)

    try:
        import SpatialDE
        import NaiveDE
    except ImportError:
        logger.error("SpatialDE or NaiveDE module not found. Please install them.")
        return adata

    # 1. Prepare Data
    # SpatialDE requires raw counts or normalized counts in specific format
    logger.info("  - Preparing data for SpatialDE...")
    
    # Needs a dataframe of counts (Cells x Genes) and Coordinates (Cells x 2)
    if scipy.sparse.issparse(adata.X):
        X_df = adata.X.toarray()
    else:
        X_df = adata.X
        
    # Filter genes with low counts to prevent SVD convergence errors
    gene_counts = np.array(X_df.sum(axis=0)).flatten()
    mask = gene_counts > 10
    X_df = X_df[:, mask]
    adata_subset = adata[:, mask]
    logger.info(f"  - Filtered {sum(~mask)} low-count genes. Remaining: {sum(mask)}")

    counts = pd.DataFrame(X_df, index=adata_subset.obs_names, columns=adata_subset.var_names)
    coords = adata.obs[['x_centroid', 'y_centroid']].astype(float)
    coords.columns = ['x', 'y']
    
    # 2. Normalize (NaiveDE)
    logger.info("  - Normalizing counts (NaiveDE)...")
    # Add total_counts to coords for regression
    coords['total_counts'] = counts.sum(axis=1)
    
    norm_expr = NaiveDE.stabilize(counts.T).T
    
    # Regress out library size effect
    # Note: regress_out expects a string formula for covariates
    resid_expr = NaiveDE.regress_out(coords, norm_expr.T, 'np.log(total_counts)').T
    
    # 3. Run SpatialDE
    logger.info("  - Running SpatialDE test (this may take time)...")
    results = SpatialDE.run(coords[['x', 'y']], resid_expr)
    
    # 4. Identify SVFs
    # Filter by qval < 0.05
    svfs = results[results.qval < 0.05]
    n_svfs = len(svfs)
    logger.info(f"✓ SpatialDE finished: Found {n_svfs} spatially variable genes.")
    
    # Sort by fraction of spatial variance (FSV)
    results = results.sort_values('FSV', ascending=False)
    
    if output_path:
        out_file = output_path / f"{sample_tag}_step8_svf_spatialde.csv"
        results.to_csv(out_file)
        logger.info(f"✓ Saved SVF results to {out_file}")
        
    return results
