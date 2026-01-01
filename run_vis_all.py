import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import logging
import argparse
import visualize_results as vis

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_all_visualizations(base_dir: str):
    base_path = Path(base_dir)
    # Define datasets to process
    datasets = [d for d in base_path.iterdir() if d.is_dir() and d.name in ['human_brain', 'h_breast_1', 'h_breast_2']]
    
    for dataset_dir in datasets:
        sample_tag = dataset_dir.name
        logger.info(f"=== Processing Visualizations for {sample_tag} ===")
        
        # 1. Load AnnData (Step 1 Output)
        adata_path = dataset_dir / f"{sample_tag}_step1_preprocessed.h5ad"
        if not adata_path.exists():
            logger.warning(f"  [!] Missing H5AD file: {adata_path}. Skipping core visualizations.")
            adata = None
        else:
            logger.info(f"  Loading AnnData: {adata_path}")
            adata = sc.read_h5ad(adata_path)
            
            # --- Generate Step 1 Visualizations (Restoring Missing Plots) ---
            logger.info("  -> Generating Step 1 Plots (QC, PCA, UMAP)...")
            # vis.visualize_step1_qc(adata, dataset_dir, sample_tag)
            # vis.visualize_pca_and_clustering(adata, dataset_dir, sample_tag)
            
            # --- Generate Marker Gene Violin Plots (New Feature) ---
            logger.info("  -> Generating Marker Gene Violin Plots...")
            if hasattr(vis, 'visualize_marker_genes_violin'):
                 vis.visualize_marker_genes_violin(adata, dataset_dir, sample_tag)
            
            # --- Generate Segmentation Overlay (Partial/Full) ---
            logger.info("  -> Generating Segmentation Overlay...")
            # We assume raw data path inference might fail, so we pass adata for fallback
            # Construct theoretical raw path (though mostly missing here)
            raw_path = Path(f"data/unprocessed_adata/{sample_tag}") 
            vis.visualize_segmentation_overlay(raw_path, dataset_dir, sample_tag, adata=adata)

        # 2. Step 2 Visualizations (Distance Plots)
        logger.info("  -> Generating Step 2 Plots...")
        vis.visualize_step2_gene_distances(dataset_dir, sample_tag) # Reads CSV internally

        # 3. Step 4 Visualizations (Optimization Curve)
        logger.info("  -> Generating Step 4 Plots...")
        vis.visualize_step4_optimization(None, dataset_dir, sample_tag)   # Reads CSV internally
        
        # 4. Step 6 Visualizations (ARI Heatmap - Pending Implementation)
        # Check if heatmap function exists or if we need to generic plot from CSV
        # For now, we rely on what is in visualize_results.py
        
        logger.info(f"=== Completed {sample_tag} ===\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Regenerate all visualizations from existing outputs.")
    parser.add_argument("--dir", default="./output", help="Base output directory")
    args = parser.parse_args()
    
    run_all_visualizations(args.dir)
