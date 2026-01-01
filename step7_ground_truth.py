import scanpy as sc
import pandas as pd
import numpy as np
import cellxgene_census
import logging
from sklearn.metrics import adjusted_rand_score
from pathlib import Path
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def download_reference_data(tissue_type: str, download_dir: str = "data/reference", organ: str = None) -> sc.AnnData:
    """
    Downloads scRNA-seq reference data for a specific tissue from Cellxgene Census.
    
    Args:
        tissue_type (str): Tissue type to filter for (e.g., "brain", "breast").
        download_dir (str): Directory to save/cache the reference data.
        organ (str): Specific organ name if 'tissue_type' is generic.
        
    Returns:
        sc.AnnData: The downloaded and concatenated reference AnnData object.
    """
    download_path = Path(download_dir)
    download_path.mkdir(parents=True, exist_ok=True)
    
    # Check if a cached file exists
    cache_file = download_path / f"reference_{tissue_type}.h5ad"
    if cache_file.exists():
        logger.info(f"Loading cached reference data from {cache_file}...")
        return sc.read_h5ad(cache_file)

    logger.info(f"Downloading reference data for tissue: {tissue_type} from Cellxgene Census...")
    
    try:
        # Open the census
        census = cellxgene_census.open_soma()
        
        # Query for human data
        # Note: This is a simplified query. You might want to filter by assay, disease status, etc.
        # For this implementation, we fetch a relevant slice.
        
        organ_filter = organ if organ else tissue_type
        
        # This is a potentially large download. We verify existence first.
        # Filtering for 'homo_sapiens' and specific tissue.
        # We limit specific datasets to avoid pulling the entire atlas if not needed.
        # For efficiency, let's grab a representative dataset if possible, or filter obs.
        
        if tissue_type == 'brain':
            # Use a specific manageable dataset from Human Brain Cell Atlas v1.0 (Midbrain PAG-DR) ~14k cells
            obs_value_filter = "dataset_id == 'c5cfa2b7-abb1-4a50-908f-707b54ca606b'"
        elif tissue_type == 'breast':
             # Use "Breast cancer 4 patients archival FFPE samples" (~10k cells) - relevant for Xenium FFPE
             obs_value_filter = "dataset_id == 'fbdd8c17-b34a-4cbc-abc4-1aeaa294a538'"
        else:
            obs_value_filter = f"tissue_general == '{organ_filter}'"
        
        # Loading AnnData
        logger.info(f"Querying Census with filter: {obs_value_filter}")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_value_filter=obs_value_filter,
            column_names={"obs": ["cell_type", "tissue", "tissue_general", "assay", "disease"]}
        )
        
        logger.info(f"Downloaded data shape: {adata.shape}")
        
        # Basic preprocessing to ensure it's usable as a reference
        logger.info("Preprocessing reference data...")
        adata.var_names_make_unique()
        
        # Filter for healthy if possible/desired, or just use all
        # adata = adata[adata.obs['disease'] == 'normal']
        
        # Normalize and log transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Identify highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        # PCA for embeddding
        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        
        # Save to cache
        logger.info(f"Saving reference data to {cache_file}...")
        adata.write(cache_file)
        
        return adata
        
    except Exception as e:
        logger.error(f"Failed to download or process reference data: {e}")
        raise

def transfer_labels(adata_xenium: sc.AnnData, adata_reference: sc.AnnData, label_key: str = "cell_type") -> sc.AnnData:
    """
    Transfers cell type labels from the reference dataset to the Xenium dataset using Scanpy's ingest.
    
    Args:
        adata_xenium (sc.AnnData): The target Xenium dataset.
        adata_reference (sc.AnnData): The annotated reference dataset.
        label_key (str): The column in adata_reference.obs containing cell type labels.
        
    Returns:
        sc.AnnData: The annotated Xenium dataset (modified in-place, but returned for convenience).
    """
    logger.info("Starting label transfer using Scanpy Ingest...")
    
    # Ensure reference has neighbors calculated (required for ingest)
    if 'neighbors' not in adata_reference.uns:
        logger.info("Calculating neighbors for reference...")
        sc.pp.pca(adata_reference)
        sc.pp.neighbors(adata_reference)

    # Robust Gene Mapping
    logger.info("checking reference gene names mapping...")
    
    xenium_is_ensembl = False
    if adata_xenium.var_names.shape[0] > 0 and str(adata_xenium.var_names[0]).startswith('ENSG'):
        xenium_is_ensembl = True
        logger.info("Xenium dataset appears to use Ensembl IDs.")
    
    # Strategy 1: Use 'feature_name' or 'gene_symbol' columns if present (Census standard)
    # BUT only if Xenium does NOT use Ensembl, or if Ref is not Ensembl.
    
    # Check if Ref is currently Ensembl
    ref_is_ensembl = False
    if adata_reference.var_names.shape[0] > 0 and str(adata_reference.var_names[0]).startswith('ENSG'):
        ref_is_ensembl = True
        logger.info("Reference dataset uses Ensembl IDs.")
        
    mapped = False
    
    if xenium_is_ensembl and ref_is_ensembl:
        logger.info("Both datasets use Ensembl IDs. Keeping Reference as Ensembl.")
        # Do NOT map to feature_name
    elif xenium_is_ensembl and not ref_is_ensembl:
        # Xenium matches Ensembl, but Reference does not (e.g. integer index or Symbol?).
        # Try to find Ensembl ID column in Reference ('feature_id' usually)
        if 'feature_id' in adata_reference.var.columns:
             logger.info("Xenium uses Ensembl, mapping Reference 'feature_id' to var_names.")
             adata_reference.var_names = adata_reference.var['feature_id']
             adata_reference.var_names_make_unique()
        else:
             logger.warning("Xenium uses Ensembl but Reference has no 'feature_id' column to map to.")
    elif not xenium_is_ensembl:
        # Xenium is likely Symbols. Map Ref to Symbols if possible.
        if 'feature_name' in adata_reference.var.columns:
            logger.info("Using 'feature_name' column for reference gene symbols (matching Xenium format).")
            adata_reference.var_names = adata_reference.var['feature_name']
            adata_reference.var_names_make_unique()
            mapped = True
        elif 'gene_symbol' in adata_reference.var.columns:
            logger.info("Using 'gene_symbol' column for reference gene symbols (matching Xenium format).")
            adata_reference.var_names = adata_reference.var['gene_symbol']
            adata_reference.var_names_make_unique()
            mapped = True
    else:
        # Fallback
        logger.info("Xenium is Ensembl, Ref is Ensembl (logic check passed).")
        
    # Strategy 2: Check intersection
    common_genes = adata_xenium.var_names.intersection(adata_reference.var_names)
    logger.info(f"Initial intersection: {len(common_genes)} genes.")
    
    # Strategy 3: Case insensitive matching if intersection is low
    if len(common_genes) < 10:
        logger.warning("Low intersection. Trying case-normalization...")
        
        # Create map {upper: original} for both
        ref_map = {x.upper(): x for x in adata_reference.var_names}
        xenium_map = {x.upper(): x for x in adata_xenium.var_names}
        
        common_upper = set(ref_map.keys()).intersection(xenium_map.keys())
        logger.info(f"Found {len(common_upper)} common genes via case-insensitive match.")
        
        if len(common_upper) >= 10:
             # Rename reference vars to match Xenium (target) format
             # We want reference var names to match Xenium so intersection works directly
             # This is complex because we need to rename only the matching ones.
             # Easier: Rename Reference to Upper, Xenium to Upper? No, Xenium defines the space.
             # Rename Reference vars to their Upper version?
             # Better: Rename Reference vars to match Xenium vars where possible.
             
             new_names = []
             for name in adata_reference.var_names:
                 upper = name.upper()
                 if upper in xenium_map:
                     new_names.append(xenium_map[upper]) # use Xenium's casing
                 else:
                     new_names.append(name)
             
             adata_reference.var_names = new_names
             adata_reference.var_names_make_unique()
             
             # Re-calc intersection
             common_genes = adata_xenium.var_names.intersection(adata_reference.var_names)
             logger.info(f"New intersection after case normalization: {len(common_genes)}")
        
    # Intersect variables (genes)
    common_genes = adata_xenium.var_names.intersection(adata_reference.var_names)
    if len(common_genes) < 10:
        logger.error(f"Only {len(common_genes)} common genes found. Cannot perform label transfer.")
        return adata_xenium
    
    logger.info(f"Found {len(common_genes)} common genes between Xenium and Reference.")
    
    # Subset Reference to common genes
    adata_ref_sub = adata_reference[:, common_genes].copy()
    
    # Subset Xenium to common genes
    # Re-process Reference on common genes
    logger.info("Re-processing reference on common genes (PCA, Neighbors)...")
    sc.pp.pca(adata_ref_sub, n_comps=min(50, len(common_genes)-1))
    sc.pp.neighbors(adata_ref_sub, n_neighbors=15, n_pcs=min(40, len(common_genes)-1))
    sc.tl.umap(adata_ref_sub) 
    
    # Create a temporary view of Xenium with only common genes for ingest
    adata_xenium_sub = adata_xenium[:, common_genes].copy()
    
    # Run ingest
    logger.info(f"Transferring labels from '{label_key}'...")
    sc.tl.ingest(adata_xenium_sub, adata_ref_sub, obs=label_key)
    
    # Copy transferred labels back to original Xenium object
    if label_key in adata_xenium_sub.obs:
        adata_xenium.obs[label_key] = adata_xenium_sub.obs[label_key]
    else:
        logger.warning(f"Label key '{label_key}' not found after ingest.")
        
    logger.info("Label transfer complete.")
    return adata_xenium

def calculate_ground_truth_ari(adata_xenium: sc.AnnData, cluster_key: str = "leiden_1_4", ground_truth_key: str = "cell_type") -> float:
    """
    Calculates the Adjusted Rand Index (ARI) between the pipeline clustering and the ground truth labels.
    
    Args:
        adata_xenium (sc.AnnData): The Xenium dataset with both clustering and transferred labels.
        cluster_key (str): Column name for Xenium clusters (e.g., 'leiden_1_4').
        ground_truth_key (str): Column name for transferred ground truth labels.
        
    Returns:
        float: The ARI score.
    """
    if cluster_key not in adata_xenium.obs:
        logger.error(f"Cluster key '{cluster_key}' not found in observations.")
        return 0.0
    
    if ground_truth_key not in adata_xenium.obs:
        logger.error(f"Ground truth key '{ground_truth_key}' not found in observations.")
        return 0.0
        
    logger.info(f"Calculating ARI between '{cluster_key}' and '{ground_truth_key}'...")
    
    labels_pred = adata_xenium.obs[cluster_key]
    labels_true = adata_xenium.obs[ground_truth_key]
    
    ari = adjusted_rand_score(labels_true, labels_pred)
    
    logger.info(f"Calculated ARI: {ari:.4f}")
    return ari

if __name__ == "__main__":
    # Test block
    print("This module provides Step 7 functionality.")
