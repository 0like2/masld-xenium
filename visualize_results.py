
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import scanpy as sc
import argparse
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configure plotting style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def visualize_step2_gene_distances(output_dir, sample_tag):
    """Visualizes gene distance to nucleus distribution."""
    logger.info(f"Visualizing Step 2 for {sample_tag}...")
    
    file_path = output_dir / f"{sample_tag}_step2_gene_distances.csv"
    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        return

    try:
        # Check for raw reads data (Parquet) first for detailed Boxplot
        raw_file_path = output_dir / f"{sample_tag}_step2_reads_sample.parquet"
        
        if raw_file_path.exists():
            logger.info(f"Examples of raw reads found at {raw_file_path}. Generating Box-plots...")
            df_raw = pd.read_parquet(raw_file_path)
            
            # Calculate mean distance to identify top genes
            gene_stats = df_raw.groupby('feature_name')['distance'].mean().sort_values(ascending=False)
            
            # Top 20 genes with highest mean distance (likely cytoplasmic/extracellular)
            top_genes = gene_stats.head(20).index
            plot_data = df_raw[df_raw['feature_name'].isin(top_genes)]
            
            # Sort plot data by mean distance
            plot_data['feature_name'] = pd.Categorical(plot_data['feature_name'], categories=top_genes, ordered=True)
            plot_data = plot_data.sort_values('feature_name')

            plt.figure(figsize=(16, 8))
            sns.boxplot(x='feature_name', y='distance', data=plot_data, fliersize=1)
            plt.xticks(rotation=90, fontsize=10)
            plt.title(f"Distance to Nucleus Distribution (Top 20 Genes) - {sample_tag}")
            plt.ylabel("Distance (µm)")
            plt.xlabel("Gene")
            plt.tight_layout()
            
            save_path = output_dir / f"{sample_tag}_step2_gene_dist_boxplot.png"
            plt.savefig(save_path, dpi=300)
            logger.info(f"✓ Saved boxplot to {save_path}")
            plt.close()
            
        else:
             logger.info("Raw reads not found. Falling back to summary bar plot.")

        # Fallback/Complementary: Summary Bar Plot from CSV
        df = pd.read_csv(file_path, index_col=0)
        df = df.sort_values('mean')
        
        plt.figure(figsize=(15, 6))
        n_display = 50
        if len(df) > 2 * n_display:
            plot_df = pd.concat([df.head(n_display), df.tail(n_display)])
        else:
            plot_df = df
            
        sns.barplot(x=plot_df.index, y=plot_df['mean'])
        plt.xticks(rotation=90, fontsize=8)
        plt.title(f"Average Distance to Nucleus by Gene ({sample_tag})")
        plt.ylabel("Mean Distance (µm)")
        plt.xlabel("Gene")
        plt.tight_layout()
        
        save_path = output_dir / f"{sample_tag}_step2_gene_dist_plot.png"
        plt.savefig(save_path, dpi=300)
        logger.info(f"✓ Saved summary plot to {save_path}")
        plt.close()

    except Exception as e:
        logger.error(f"Error plotting Step 2: {e}")

def visualize_step4_optimal_expansion(output_dir, sample_tag):
    """Visualizes optimal expansion distance results."""
    logger.info(f"Visualizing Step 4 for {sample_tag}...")
    
    file_path = output_dir / f"{sample_tag}_step4_optimal_expansion.csv"
    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        return

    try:
        df = pd.read_csv(file_path, index_col=0)
        
        # Check if it's the new simulation format (has columns Capture_Reads, Global_Purity_Proxy)
        if 'Capture_Reads' in df.columns:
            fig, ax1 = plt.subplots(figsize=(10, 6))
            
            # Plot Capture (Reads) on left y-axis
            color = 'tab:blue'
            ax1.set_xlabel('Expansion Distance (µm)')
            ax1.set_ylabel('Capture (Total Assigned Reads)', color=color)
            ax1.plot(df.index, df['Capture_Reads'], color=color, marker='o', label='Capture')
            ax1.tick_params(axis='y', labelcolor=color)
            ax1.grid(True, alpha=0.3)

            # Plot Purity on right y-axis
            ax2 = ax1.twinx()  
            color = 'tab:red'
            ax2.set_ylabel('Pseudo-Purity (Nuclear Fraction)', color=color)  
            ax2.plot(df.index, df['Global_Purity_Proxy'], color=color, marker='s', linestyle='--', label='Purity')
            ax2.tick_params(axis='y', labelcolor=color)
            ax2.set_ylim(0, 1.0) # Purity is 0-1

            plt.title(f"Optimal Expansion Optimization (Purity vs Capture) - {sample_tag}")
            fig.tight_layout()
            
            save_path = output_dir / f"{sample_tag}_step4_optimization_curve.png"
            plt.savefig(save_path, dpi=300)
            logger.info(f"✓ Saved optimization curve to {save_path}")
            plt.close()
            
        else:
            # Fallback for old format
            metrics = df.index.tolist()
            values = df.iloc[:, 0].tolist()
            
            plt.figure(figsize=(10, 6))
            sns.barplot(x=metrics, y=values)
            plt.title(f"Optimal Cell Expansion Parameters ({sample_tag})")
            plt.ylabel("Distance (µm)")
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            save_path = output_dir / f"{sample_tag}_step4_expansion_plot.png"
            plt.savefig(save_path, dpi=300)
            logger.info(f"✓ Saved plot to {save_path}")
            plt.close()

    except Exception as e:
        logger.error(f"Error plotting Step 4: {e}")

def visualize_pca_and_clustering(adata: sc.AnnData, output_dir: Path, sample_tag: str):
    """Visualizes PCA and clustering results from AnnData."""
    logger.info(f"Visualizing PCA/Clustering for {sample_tag}...")
    
    try:
        # Ensure PCA is computed
        if 'X_pca' not in adata.obsm:
            logger.info("Computing PCA for visualization...")
            sc.pp.pca(adata)
            
        # Plot PCA Variance Ratio
        sc.pl.pca_variance_ratio(adata, show=False)
        plt.title(f"PCA Variance Ratio ({sample_tag})")
        plt.tight_layout()
        plt.savefig(output_dir / f"{sample_tag}_step1_pca_variance.png", dpi=300)
        plt.close()
        
        # Plot PCA 2D
        # Use 'leiden' or 'leiden_1_4' if available, else no color
        color_key = None
        if 'leiden_1_4' in adata.obs:
            color_key = 'leiden_1_4'
        elif 'leiden' in adata.obs:
            color_key = 'leiden'
        
        sc.pl.pca(adata, color=color_key, show=False)
        plt.title(f"PCA 2D Projection ({sample_tag})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(output_dir / f"{sample_tag}_step1_pca_plot.png", dpi=300)
        plt.close()
        
        logger.info("✓ Saved PCA plots.")

    except Exception as e:
        logger.error(f"Error visualizing PCA: {e}")

def visualize_step1_qc(adata: sc.AnnData, output_dir: Path, sample_tag: str):
    """Visualizes QC metrics: Violin plots for counts and genes."""
    logger.info(f"Visualizing Step 1 QC for {sample_tag}...")
    
    try:
        # Calculate QC metrics if not present
        if 'n_genes_by_counts' not in adata.obs:
            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            
        # Violin Plots
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        sns.violinplot(y=adata.obs['n_genes_by_counts'], ax=axs[0], color='skyblue')
        axs[0].set_title(f"Genes per Cell ({sample_tag})")
        axs[0].set_ylabel("Number of Genes")
        
        sns.violinplot(y=adata.obs['total_counts'], ax=axs[1], color='lightgreen')
        axs[1].set_title(f"Counts per Cell ({sample_tag})")
        axs[1].set_ylabel("Total Counts")
        
        plt.tight_layout()
        save_path = output_dir / f"{sample_tag}_step1_qc_violin.png"
        plt.savefig(save_path, dpi=300)
        plt.close()
        logger.info(f"✓ Saved QC plots to {save_path}")

    except Exception as e:
        logger.error(f"Error visualizing QC: {e}")

def visualize_pca_and_clustering(adata: sc.AnnData, output_dir: Path, sample_tag: str):
    """Visualizes PCA, UMAP and clustering results from AnnData."""
    logger.info(f"Visualizing PCA/UMAP/Clustering for {sample_tag}...")
    
    try:
        # Ensure PCA is computed
        if 'X_pca' not in adata.obsm:
            logger.info("Computing PCA for visualization...")
            sc.pp.pca(adata)
            
        # Plot PCA Variance Ratio
        sc.pl.pca_variance_ratio(adata, show=False)
        plt.title(f"PCA Variance Ratio ({sample_tag})")
        plt.tight_layout()
        plt.savefig(output_dir / f"{sample_tag}_step1_pca_variance.png", dpi=300)
        plt.close()
        
        # Plot PCA 2D
        color_key = None
        if 'leiden_1_4' in adata.obs:
            color_key = 'leiden_1_4'
        elif 'leiden' in adata.obs:
            color_key = 'leiden'
        elif 'louvain_1_4' in adata.obs:
            color_key = 'louvain_1_4'
        
        sc.pl.pca(adata, color=color_key, show=False)
        plt.title(f"PCA 2D Projection ({sample_tag})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(output_dir / f"{sample_tag}_step1_pca_plot.png", dpi=300)
        plt.close()
        
        # Plot UMAP
        if 'X_umap' not in adata.obsm:
             logger.info("Computing Neighborhood Graph & UMAP for visualization...")
             sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
             sc.tl.umap(adata)
        
        sc.pl.umap(adata, color=color_key, show=False)
        plt.title(f"UMAP Projection ({sample_tag})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(output_dir / f"{sample_tag}_step1_umap_plot.png", dpi=300)
        plt.close()
        
        logger.info("✓ Saved PCA and UMAP plots.")

    except Exception as e:
        logger.error(f"Error visualizing PCA/UMAP: {e}")

# ... (previous code)

def visualize_marker_genes(adata: sc.AnnData, output_dir: Path, sample_tag: str, n_genes: int = 5):
    """
    Identifies marker genes for each cluster and plots a dotplot.
    High Priority Feature.
    """
    logger.info(f"Visualizing Marker Genes for {sample_tag}...")
    
    try:
        # Check if clustering exists
        cluster_key = None
        if 'leiden_1_4' in adata.obs:
            cluster_key = 'leiden_1_4'
        elif 'leiden' in adata.obs:
            cluster_key = 'leiden'
            
        if not cluster_key:
            logger.warning("No clustering found (leiden/leiden_1_4). Skipping marker genes.")
            return

        # Rank genes
        logger.info(f"Ranking genes for clusters ({cluster_key})...")
        # Use simple t-test or wilcoxon. 't-test' is faster for large data.
        sc.tl.rank_genes_groups(adata, cluster_key, method='t-test')
        
        # Plot Dotplot
        plt.figure(figsize=(12, 6))
        # Swap axes: Genes on X, Clusters on Y is standard scanpy dotplot
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, standard_scale='var', show=False)
        plt.title(f"Top {n_genes} Marker Genes per Cluster - {sample_tag}")
        
        save_path = output_dir / f"{sample_tag}_step1_marker_genes_dotplot.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"✓ Saved Marker Gene Dotplot to {save_path}")

        # Also save the list of markers
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        markers_df = pd.DataFrame(
            {group + '_' + key: result[key][group]
            for group in groups for key in ['names', 'pvals']}
        ).head(n_genes)
        markers_df.to_csv(output_dir / f"{sample_tag}_step1_marker_genes.csv")

    except Exception as e:
        logger.error(f"Error determining marker genes: {e}")

def visualize_specific_genes(adata: sc.AnnData, output_dir: Path, sample_tag: str, genes: list):
    """
    Plots UMAP and Spatial maps for user-defined specific genes.
    Medium Priority Feature.
    """
    if not genes:
        return

    logger.info(f"Visualizing Specific Genes {genes} for {sample_tag}...")
    
    # Resolve gene names (Symbols -> Index)
    # Some datasets (human_brain) have Ensembl IDs as index and Symbols in 'gene_id' column
    genes_to_plot = [] # List of (index_name, display_name)
    
    for g in genes:
        # 1. Check Index
        if g in adata.var_names:
            genes_to_plot.append((g, g))
            continue
            
        # 2. Check columns (gene_id, feature_name)
        found = False
        for col in ['gene_id', 'feature_name', 'Gene', 'Symbol']:
            if col in adata.var.columns:
                # Find rows where col matches g (case insensitive?)
                # Assuming exact match for now to match config
                matches = adata.var.index[adata.var[col] == g].tolist()
                if matches:
                    # Use the first match's index
                    genes_to_plot.append((matches[0], g))
                    found = True
                    break
        if not found:
            logger.warning(f"Gene '{g}' not found in index or gene columns.")

    if not genes_to_plot:
        logger.warning("No valid genes found to plot.")
        return

    gene_dir = output_dir / "gene_plots"
    gene_dir.mkdir(exist_ok=True)

    try:
        # Plot UMAP Expression
        for gene_idx, gene_name in genes_to_plot:
            # Title for plot
            title_suffix = f"- {sample_tag}"
            if gene_idx != gene_name:
                title_suffix = f"({gene_name}) - {sample_tag}"

            # UMAP
            if 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color=gene_idx, show=False, color_map='viridis')
                plt.title(f"{gene_name} Expression (UMAP) {title_suffix}")
                plt.savefig(gene_dir / f"{sample_tag}_umap_{gene_name}.png", dpi=300, bbox_inches='tight')
                plt.close()
            
            # Spatial (if coords exist)
            if 'spatial' in adata.obsm or ('x_centroid' in adata.obs and 'y_centroid' in adata.obs):
                # Construct spatial obsm if missing but coords exist
                if 'spatial' not in adata.obsm and 'x_centroid' in adata.obs:
                     adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

                sc.pl.spatial(adata, color=gene_idx, spot_size=30, show=False, color_map='viridis')
                plt.title(f"{gene_name} Expression (Spatial) {title_suffix}")
                plt.savefig(gene_dir / f"{sample_tag}_spatial_{gene_name}.png", dpi=300, bbox_inches='tight')
                plt.close()
                
        logger.info(f"✓ Saved gene expression plots to {gene_dir}")

    except Exception as e:
        logger.error(f"Error plotting specific genes: {e}")


def visualize_step4_optimization(adata: sc.AnnData, output_dir: Path, sample_tag: str):
    """Visualizes the cell expansion optimization curve (Step 4)."""
    logger.info(f"Visualizing Step 4 Optimization Curve for {sample_tag}...")
    
    try:
        # Note: We expect the simulation data in this specific file
        csv_path = output_dir / f"{sample_tag}_step4_expansion_curve.csv"
        
        if not csv_path.exists():
            # Fallback for legacy runs (or if I didn't re-run yet)
            legacy_path = output_dir / f"{sample_tag}_step4_optimal_expansion.csv"
            if legacy_path.exists():
                 # Check if it has multiple rows
                 df_check = pd.read_csv(legacy_path)
                 if len(df_check) > 2:
                     csv_path = legacy_path
                 else:
                     logger.warning("Step 4 Curve Data not found (only summary exists). Please re-run Step 4.")
                     return
            else:
                 logger.warning(f"Step 4 Data not found: {csv_path}")
                 return

        df = pd.read_csv(csv_path)
        
        # Plotting
        if 'Expansion_Distance' not in df.columns and df.index.name != 'Expansion_Distance':
             if df.columns[0] == 'Unnamed: 0':
                 df.set_index('Unnamed: 0', inplace=True)
                 df.index.name = 'Expansion_Distance'

        fig, ax1 = plt.subplots(figsize=(10, 6))

        color = 'tab:blue'
        ax1.set_xlabel('Expansion Distance (µm)')
        ax1.set_ylabel('Total Transcripts Captured', color=color)
        ax1.plot(df.index, df['Capture_Reads'], color=color, marker='o', label='Transcripts')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, alpha=0.3)

        ax2 = ax1.twinx() 
        color = 'tab:orange'
        ax2.set_ylabel('Purity (Nuclear Fraction)', color=color) 
        ax2.plot(df.index, df['Global_Purity_Proxy'], color=color, linestyle='--', marker='x', label='Purity')
        ax2.tick_params(axis='y', labelcolor=color)

        plt.title(f"Cell Boundary Expansion Optimization ({sample_tag})")
        fig.tight_layout()
        
        output_plot = output_dir / f"{sample_tag}_step4_optimization_curve.png"
        plt.savefig(output_plot, dpi=300)
        plt.close()
        logger.info(f"✓ Saved Optimization Curve: {output_plot}")

    except Exception as e:
        logger.error(f"Error visualizing Step 4: {e}")


def visualize_segmentation_overlay(xenium_input: Path, output_dir: Path, sample_tag: str, roi_size_um: float = 200.0):
    """
    Visualizes actual cell segmentation shapes (boundaries) and transcripts in a specific ROI.
    Paper-Grade Visual (Image 5c style).
    """
    logger.info(f"Visualizing Segmentation Overlay (ROI {roi_size_um}µm) for {sample_tag}...")
    
    try:
        xenium_input = Path(xenium_input)
        
        # 1. Determine ROI center using transcripts summary or centroids if valid
        # We can try reading a subsample of transcripts or just assume center of image?
        # Better: use existing adata centroids if available, else read metadata.
        # Im taking a shortcut: Use transcripts.parquet metadata or just read it (pandas is fast enough for header?)
        # Let's assume we can read transcripts.parquet reasonably fast or use adata if passed.
        # To keep signature simple, let's just peek at transcripts.parquet middle.
        
        transcripts_path = xenium_input / "transcripts.parquet"
        if not transcripts_path.exists():
            transcripts_path = xenium_input / "transcripts.csv" # Fallback
            if not transcripts_path.exists():
                logger.warning(f"Transcripts file not found in {xenium_input}")
                return

        # Load specific columns to find bounds
        # pd.read_parquet supports columns.
        if transcripts_path.suffix == '.parquet':
            df_loc = pd.read_parquet(transcripts_path, columns=['x_location', 'y_location'])
        else:
            df_loc = pd.read_csv(transcripts_path, usecols=['x_location', 'y_location'])
            
        x_min, x_max = df_loc['x_location'].min(), df_loc['x_location'].max()
        y_min, y_max = df_loc['y_location'].min(), df_loc['y_location'].max()
        
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        
        # Define ROI box
        roi_x_min = x_center - (roi_size_um / 2)
        roi_x_max = x_center + (roi_size_um / 2)
        roi_y_min = y_center - (roi_size_um / 2)
        roi_y_max = y_center + (roi_size_um / 2)
        
        logger.info(f"ROI Center: ({x_center:.1f}, {y_center:.1f}), Size: {roi_size_um}µm")
        
        # 2. Load Data in ROI
        # Transcripts
        df_transcripts = df_loc[
            (df_loc['x_location'] >= roi_x_min) & (df_loc['x_location'] <= roi_x_max) &
            (df_loc['y_location'] >= roi_y_min) & (df_loc['y_location'] <= roi_y_max)
        ]
        
        # Boundaries
        # Try 'cell_boundaries.parquet' (newer) or 'nucleus_boundaries.parquet'
        bounds_file = xenium_input / "cell_boundaries.parquet"
        bounds_type = "Cell"
        if not bounds_file.exists():
            bounds_file = xenium_input / "nucleus_boundaries.parquet"
            bounds_type = "Nucleus"
            if not bounds_file.exists():
                logger.warning("No boundary file (cell/nucleus) found.")
                return 

        df_bounds = pd.read_parquet(bounds_file)
        # vertex_x, vertex_y
        df_bounds = df_bounds[
            (df_bounds['vertex_x'] >= roi_x_min) & (df_bounds['vertex_x'] <= roi_x_max) &
            (df_bounds['vertex_y'] >= roi_y_min) & (df_bounds['vertex_y'] <= roi_y_max)
        ]
        
        # 3. Plot
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Plot Transcripts (Dots)
        ax.scatter(df_transcripts['x_location'], df_transcripts['y_location'], 
                   s=1, c='gray', alpha=0.5, label='Transcripts')
        
        # Plot Boundaries (Polygons)
        # Group by cell_id
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        import numpy as np

        patches = []
        # Optimizing: only plot cells that have vertices in ROI
        # (Already filtered df_bounds)
        
        for cell_id, grp in df_bounds.groupby('cell_id'):
            if len(grp) < 3: continue
            poly = Polygon(grp[['vertex_x', 'vertex_y']].values, closed=True)
            patches.append(poly)
            
        p = PatchCollection(patches, alpha=0.4, edgecolor='orange', facecolor='none', linewidth=1.5, label=f'{bounds_type} Boundaries')
        ax.add_collection(p)
        
        # Refine View
        ax.set_xlim(roi_x_min, roi_x_max)
        ax.set_ylim(roi_y_min, roi_y_max)
        ax.set_aspect('equal')
        ax.set_title(f"Segmentation Overlay (ROI Center) - {sample_tag}\n{bounds_type} Boundaries + Transcripts")
        ax.legend(loc='upper right')
        
        output_plot = output_dir / f"{sample_tag}_segmentation_overlay_roi.png"
        plt.savefig(output_plot, dpi=300)
        plt.close()
        logger.info(f"✓ Saved Segmentation Overlay: {output_plot}")

    except Exception as e:
        logger.error(f"Error visualizing segmentation overlay: {e}")

def main():
    parser = argparse.ArgumentParser(description="Visualize Xenium Pipeline Results")
    parser.add_argument("--steps", nargs="+", default=["1", "2", "4", "6"], help="Steps to visualize (1, 2, 4, 6)")
    # Default includes Step 1 now
    
    args = parser.parse_args()
    
    # Define dataset paths
    datasets = [
        # (output_dir, sample_tag)
        (Path("./output/human_brain"), "human_brain"),
        (Path("./output/h_breast_1"), "h_breast_1")
    ]
    
    for output_dir, sample_tag in datasets:
        logger.info(f"--- Visualizing {sample_tag} ---")
        if not output_dir.exists():
            logger.warning(f"Output directory {output_dir} does not exist. Skipping.")
            continue
            
        # 1. Visualize CSV-based results
        if "2" in args.steps:
            visualize_step2_gene_distances(output_dir, sample_tag)
        if "4" in args.steps:
            visualize_step4_optimal_expansion(output_dir, sample_tag)
        if "6" in args.steps:
            visualize_step6_simulation(output_dir, sample_tag)

        # 2. Visualize AnnData-based results (Step 1 QC, PCA, UMAP, Clusters)
        if "1" in args.steps or "qc" in args.steps: # Allow 'qc' as explicit step alias
             adata_path = output_dir / f"{sample_tag}_step1_preprocessed.h5ad"
             
             if adata_path.exists():
                try:
                    logger.info(f"Loading AnnData from {adata_path}...")
                    adata = sc.read_h5ad(adata_path) 
                    
                    if "1" in args.steps or "qc" in args.steps:
                        visualize_step1_qc(adata, output_dir, sample_tag)
                        visualize_pca_and_clustering(adata, output_dir, sample_tag)
                        visualize_spatial_clusters(adata, output_dir, sample_tag)
                        
                        # New Feature: Marker Genes
                        visualize_marker_genes(adata, output_dir, sample_tag)
                        
                        # New Feature: Specific Genes (Example)
                        # In a real pipeline run, these would come from config
                        example_genes = ['Gfap', 'Olig2', 'Slc17a7'] if 'brain' in sample_tag else ['EPCAM', 'CD3D']
                        visualize_specific_genes(adata, output_dir, sample_tag, example_genes)
                    
                except Exception as e:
                    logger.error(f"Failed to load AnnData: {e}")
             else:
                logger.warning(f"AnnData file not found (Step 1 needed): {adata_path}")

def visualize_step6_simulation(output_dir: Path, sample_tag: str):
    """
    Visualize Step 6 Simulation Results: Heatmap of ARI scores.
    """
    print(f"Generating Step 6 Simulation Heatmap for {sample_tag}...")
    results_file = output_dir / "step6_simulation" / "simulation_results.csv"
    if not results_file.exists():
        print(f"Warning: Results file not found at {results_file}")
        return

    try:
        df = pd.read_csv(results_file)
        # Assuming df has columns like 'radius', 'quantile', 'ARI'
        # Pivot for heatmap
        if {'radius', 'quantile', 'ARI'}.issubset(df.columns):
            pivot_df = df.pivot(index='radius', columns='quantile', values='ARI')
            plt.figure(figsize=(10, 8))
            sns.heatmap(pivot_df, annot=True, cmap="viridis", fmt=".3f")
            plt.title(f"Simulation ARI Scores: Radius vs Quantile ({sample_tag})")
            plt.xlabel("Quantile")
            plt.ylabel("Radius")
            plt.tight_layout()
            plt.savefig(output_dir / "step6_simulation" / "simulation_ari_heatmap.png")
            plt.close()
            print("Step 6 Heatmap generated.")
        else:
            print("Warning: Step 6 results dataframe missing expected columns for heatmap.")
    except Exception as e:
        print(f"Error generating Step 6 heatmap: {e}")

def visualize_step8_svf(adata: sc.AnnData, output_dir: Path, sample_tag: str):
    """
    Visualize Step 8 SVF Results: Spatial maps of top spatially variable genes.
    """
    print(f"Generating Step 8 SVF Maps for {sample_tag}...")
    svf_dir = output_dir / "step8_svf"
    svf_dir.mkdir(parents=True, exist_ok=True)
    
    # Load SVF results
    results_file = svf_dir / "spatialde_results.csv"
    if not results_file.exists():
        print(f"Warning: SVF results file not found at {results_file}")
        return

    try:
        svf_df = pd.read_csv(results_file)
        # Get top 9 SVGs by q-value
        if 'qval' in svf_df.columns:
            top_genes = svf_df.sort_values("qval").head(9)['g'].tolist()
            
            # Plot spatial expression
            sc.pl.spatial(adata, color=top_genes, ncols=3, show=False, spot_size=30)
            plt.savefig(svf_dir / "top_svgs_spatial.png")
            plt.close()
            print("Step 8 SVF Maps generated.")
        else:
            print("Warning: SVF results missing 'qval' column.")
    except Exception as e:
        print(f"Error generating Step 8 SVF maps: {e}")
