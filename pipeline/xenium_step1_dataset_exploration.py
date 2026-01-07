# ------------------------------------------------------------
# Xenium Pipeline Step 1: Exploratory Metrics & Analysis
# Calculates statistics, neighborhood patterns (1_7), and plots.
# STANDALONE VERSION
# ------------------------------------------------------------

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns

# Configure local logging if run mainly
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def plotting_settings():
    """Applies visualization settings matching the notebooks"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white', figsize=(10, 10))
    
    # Custom matplotlib/seaborn settings
    # plt.style.use('seaborn-white') # Invalid in newer mpl versions
    sns.set_style("white")
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    plt.rcParams['figure.facecolor'] = 'white'
    # sns.set_style("whitegrid") # Overridden by seaborn-white preference

def calculate_general_stats(adata, output_dir, sample_tag):
    """Calculates general dataset statistics (Ref: 1_1 notebook)"""
    print("\n[Step 1-1] Calculating General Statistics (Ref: 1_1)...")
    
    # Ensure QC metrics are fresh
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    stats = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_genes_per_cell": np.median(adata.obs['n_genes_by_counts']),
        "median_counts_per_cell": np.median(adata.obs['total_counts']),
        "total_counts": np.sum(adata.obs['total_counts'])
    }
    
    # Cell area stats if available
    if 'cell_area' in adata.obs:
         stats["median_cell_area"] = np.median(adata.obs['cell_area'])

    # Save stats to simple text report
    report_file = os.path.join(output_dir, f"{sample_tag}_step1_stats.txt")
    with open(report_file, "w") as f:
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")
    print(f"    - Saved stats to {report_file}")
    
    # Plot Histograms
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    sns.histplot(adata.obs['total_counts'], bins=50, ax=axs[0], color='skyblue')
    axs[0].set_title('Total Counts per Cell')
    sns.histplot(adata.obs['n_genes_by_counts'], bins=50, ax=axs[1], color='lightgreen')
    axs[1].set_title('Genes Detected per Cell')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{sample_tag}_step1_qc_hist.png"))
    plt.close()

    # Plot Highest Expressed Genes (Top 20)
    print("    - Plotting highest expressed genes...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    plt.title("Highest Expressed Genes")
    plt.savefig(os.path.join(output_dir, f"{sample_tag}_step1_highest_expr_genes.png"), bbox_inches='tight')
    plt.close()

def calculate_transcript_dispersion(adata, output_dir, sample_tag):
    """
    Calculates distance of each transcript to its assigned cell centroid (Ref: 1_3).
    Optimized vectorized implementation.
    """
    print("\n[Step 1-X] Deep QC: Transcript Dispersion Analysis (Ref: 1_3)...")
    
    if 'spots' not in adata.uns:
        print("    [WARNING] No 'spots' (transcripts) data found in adata.uns. Skipping.")
        return

    # 1. Prepare Data
    print("    - Preparing transcript data...")
    spots = adata.uns['spots'].copy()
    
    # Needs cell info with centroids
    if 'x_centroid' not in adata.obs.columns or 'y_centroid' not in adata.obs.columns:
        print("    [WARNING] Cell centroids (x_centroid, y_centroid) missing in adata.obs. Skipping.")
        return
        
    # Map cell centroids to spots based on 'cell_id' or index
    # Note: 'spots' usually has 'cell_id' column or index is cell_id? 
    # Let's inspect typical structure. Step0 loads transcripts.csv.
    # Usually 'cell_id' is a column in transcripts.csv if it's the assigned transcripts file.
    
    if 'cell_id' not in spots.columns:
        # Sometimes 'cell_id' is the index if loaded that way, or maybe 'cell_uuid'
        # If standard Xenium format, 'cell_id' should be there for assigned transcripts.
        # If unassigned ones are included (cell_id = UNASSIGNED), we filter them.
        print("    [WARNING] 'cell_id' column missing in spots dataframe. Checking index/content...")
        if spots.index.name == 'cell_id':
            spots = spots.reset_index()
        else:
             print("    [ERROR] Cannot link spots to cells. Skipping.")
             return
             
    # Filter for assigned transcripts only
    n_total = len(spots)
    # Assuming 'cell_id' is string matching adata.obs.index or 'UNASSIGNED'
    # Filter out unassigned or negative IDs
    spots_assigned = pd.DataFrame() # Initialize
    
    # Strategy: Merge on 'cell_id' column if available, else index
    match_col = 'cell_id' if 'cell_id' in adata.obs.columns else None
    
    if match_col:
        print(f"    - Linking spots to cells using 'cell_id' column...")
        # Check overlap
        overlap = spots['cell_id'].isin(adata.obs[match_col]).sum()
        if overlap > 0:
            spots_assigned = spots[spots['cell_id'].isin(adata.obs[match_col])]
        else:
             print("    [WARNING] 'cell_id' column exists but values don't overlap with spots. Trying index...")
             match_col = None
    
    if not match_col:
        # Fallback to index
        if 'cell_id' in spots.columns:
            overlap = spots['cell_id'].isin(adata.obs.index).sum()
            if overlap > 0:
                 print(f"    - Linking spots to cells using adata index...")
                 spots_assigned = spots[spots['cell_id'].isin(adata.obs.index)]
                 
    n_assigned = len(spots_assigned)

    if n_assigned == 0:
        print(f"    [ERROR] Could not link spots to cells (checked 'cell_id' col and index). Deep QC failed.")
        return

    print(f"    - analyzing {n_assigned} assigned transcripts (out of {n_total})...")
    
    # 2. Vectorized Distance Calculation
    print("    - Merging cell coordinates...")
    
    if match_col:
        merged = spots_assigned.merge(
            adata.obs[['x_centroid', 'y_centroid', match_col]], 
            left_on='cell_id', 
            right_on=match_col, 
            how='left'
        )
    else:
        merged = spots_assigned.merge(
            adata.obs[['x_centroid', 'y_centroid']], 
            left_on='cell_id', 
            right_index=True, 
            how='left'
        )
    
    print("    - Computing Euclidean distances...")
    # dist = sqrt((x_spot - x_cell)^2 + (y_spot - y_cell)^2)
    # spots usually have 'x_location', 'y_location'
    dx = merged['x_location'] - merged['x_centroid']
    dy = merged['y_location'] - merged['y_centroid']
    distances = np.sqrt(dx**2 + dy**2)
    
    # 3. Stats & Plotting
    mean_dist = np.mean(distances)
    median_dist = np.median(distances)
    print(f"    - Median distance to centroid: {median_dist:.2f} (unit)")
    
    # Plot Distribution
    plt.figure(figsize=(8, 6))
    sns.histplot(distances, bins=100, kde=True, color='purple')
    plt.title(f"Transcript-Centroid Distance Distribution\nMedian: {median_dist:.2f}")
    plt.xlabel("Distance")
    plt.ylabel("Count")
    
    save_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_dist.png")
    plt.savefig(save_path)
    plt.close()
    print(f"    - Saved Dispersion Plot to {save_path}")
    
    # Save metrics to a file
    with open(os.path.join(output_dir, f"{sample_tag}_step1_dispersion_metrics.txt"), "w") as f:
        f.write(f"Total Transcripts: {n_total}\n")
        f.write(f"Assigned Transcripts: {n_assigned}\n")
        f.write(f"Median Distance: {median_dist}\n")
        f.write(f"Mean Distance: {mean_dist}\n")


def analyze_neighborhoods(adata, output_dir, sample_tag, radius=20.0):
    """Calculates neighborhood metrics (Ref: 1_7 notebook)"""
    print("\n[Step 1-3] Analyzing Neighborhood Architecture (Ref: 1_7)...")
    
    if 'spatial' not in adata.obsm:
        print("    [WARNING] No spatial coordinates found. Skipping neighborhood analysis.")
        return

    # Use Squidpy to compute spatial neighbors
    # Notebook 1_7 uses radius=100 for this specific analysis
    print(f"    - Computing spatial neighbors graph (radius={radius})...")
    sq.gr.spatial_neighbors(adata, coord_type="generic", radius=radius) 
    
    # --- Gap Filled: Custom Neighborhood Metrics (Ref: 1_7) ---
    # Notebook calculates 'neighborhood_diversity' and 'neighborhood_density'
    # Logic: sum of neighbors graph connectivity
    
    # Check if we have the connectivity matrix
    if 'spatial_connectivities' in adata.obsp:
        print("    - Calculating custom neighborhood metrics (diversity, density)...")
        # Density: Sum of connections per cell (degree)
        # Convert sparse matrix to dense sum per row
        adata.obs['neighborhood_density'] = np.array(adata.obsp['spatial_connectivities'].sum(axis=1)).flatten()
        
        # Diversity requires cell types. Check for primary cluster key.
        cluster_key = None
        for key in ['leiden_1.0', 'leiden_0.8', 'leiden', 'graph_clusters']:
            if key in adata.obs:
                cluster_key = key
                break
        
        if cluster_key:
            # Diversity calculation: # of unique clusters in neighborhood
            # This is expensive to compute via pure numpy on large matrices without specific optimization.
            # Notebook 1_7 creates 'adataneigh' by summing one-hot encodings of neighbors
            # Simplified approach: If density exists, we assume density is sufficient for now, 
            # or we can implement a diversity score if strictly required. 
            # Given the request for missing logic:
            
            # Implementation of Nhood Diversity (count of distinct labels in neighborhood)
            # This is effectively: sum( (Adj * OneHot) > 0 )
            
            dummies = pd.get_dummies(adata.obs[cluster_key])
            adj = adata.obsp['spatial_connectivities']
            # matmul: (Cells x Cells) * (Cells x Clusters) = (Cells x Clusters) -> weighted count of each cluster in nhood
            # We want binary presence:
            # Note: adj is typically 1s and 0s. 
            
            # Optimized sparse multiplication
            from scipy import sparse
            if not sparse.issparse(adj):
                adj = sparse.csr_matrix(adj)
            
            cluster_counts = adj.dot(sparse.csr_matrix(dummies.values))
            # Diversity = count of non-zero entries per row
            adata.obs['neighborhood_diversity'] = np.array((cluster_counts > 0).sum(axis=1)).flatten()
            
            print(f"      - Diversity calculated using '{cluster_key}'")
    # ----------------------------------------------------------
    
    # Calculate Neighborhood Enrichment (which cell types sit next to each other?)
    # Use machine-provided clusters (graph_clusters) if available
    cluster_key = None
    possible_keys = ['graph_clusters', 'kmeans10_clusters', 'leiden']
    
    for key in possible_keys:
        if key in adata.obs.keys():
            cluster_key = key
            print(f"    - Using cluster key '{cluster_key}' for enrichment analysis.")
            break
            
    if cluster_key:
        print(f"    - Computing neighborhood enrichment for '{cluster_key}'...")
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
        
        # Plot Enrichment
        plt.figure(figsize=(8, 8))
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_key)
        plt.title(f"Neighborhood Enrichment ({cluster_key})")
        save_path = os.path.join(output_dir, f"{sample_tag}_step1_grad_enrichment.png")
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"    - Saved enrichment plot to {save_path}")
        
        # Calculate Centrality scores (how 'connected' is each cell type?)
        print(f"    - Computing centrality scores...")
        sq.gr.centrality_scores(adata, cluster_key=cluster_key)
        
        # Plot Centrality
        plt.figure(figsize=(10, 5))
        sq.pl.centrality_scores(adata, cluster_key=cluster_key)
        plt.title("Centrality Scores")
        save_path_cent = os.path.join(output_dir, f"{sample_tag}_step1_centrality.png")
        plt.savefig(save_path_cent, bbox_inches='tight')
        plt.close() 
    else:
        print("    [WARNING] No cluster key found. Skipping enrichment/centrality.")


def perform_clustering_and_annotation(adata, output_dir, sample_tag, config):
    """
    Performs Clustering (Leiden) and Marker Ranking (Ref: 1_2 notebook)
    Merged from xenium_step3_nuclei_annotation.py
    """
    print("\n[Step 1-2] Running Clustering & Annotation (Ref: 1_2)...")

    # Ensure PCA and Neighbors are present (sanity check for clustering)
    if 'X_pca' not in adata.obsm:
        print("    - PCA not found. Running PCA...")
        sc.pp.pca(adata)
    if 'neighbors' not in adata.uns:
        print("    - Neighbors not found. Computing Neighbors...")
        sc.pp.neighbors(adata)

    # --- Gap Filled: HVG Calculation & Plotting (User Request) ---
    print("    - Calculating Highly Variable Genes (ref request)...")
    try:
        # Parameters from user request: min_mean=0.3, max_mean=7, min_disp=-0.5
        # Note: These values might need adjustment if log1p was done or not. 
        # Typically run on logged counts.
        if 'log1p' not in adata.uns:
             sc.pp.log1p(adata)
             
        sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
        
        # Plot HVG
        sc.pl.highly_variable_genes(adata, show=False)
        plt.title("Highly Variable Genes")
        hvg_file = os.path.join(output_dir, f"{sample_tag}_step1_hvg.png")
        plt.savefig(hvg_file)
        plt.close()
        print(f"    - Saved HVG Plot to: {hvg_file}")
    except Exception as e:
        print(f"    [WARNING] HVG calculation/plotting failed: {e}")
    # -------------------------------------------------------------

    # Run Leiden at multiple resolutions
    # Note: Using 'exploration' section as source now, fallback to 'annotation' if missing
    resolutions = config.get("exploration", {}).get("resolutions", 
                   config.get("annotation", {}).get("resolutions", [0.5, 0.8, 1.0]))
    
    for res in resolutions:
        key = f"leiden_{res}"
        print(f"    - Running Leiden clustering (resolution={res})...")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        print(f"      - Found {len(adata.obs[key].unique())} clusters.")

    # Select primary resolution
    primary_res = config.get("exploration", {}).get("primary_resolution", 
                    config.get("annotation", {}).get("primary_resolution", 1.0))
    primary_key = f"leiden_{primary_res}"
    if primary_key not in adata.obs:
         print(f"    - Primary resolution {primary_res} missing. Computing...")
         sc.tl.leiden(adata, resolution=primary_res, key_added=primary_key)
    
    print(f"    - Using '{primary_key}' as primary clustering.")

    # Rank genes
    print("    - Ranking marker genes...")
    sc.tl.rank_genes_groups(adata, groupby=primary_key, method='wilcoxon')
    
    # Extract top markers
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    markers_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']}
    ).head(5)
    
    # Save markers
    markers_file = os.path.join(output_dir, f"{sample_tag}_step1_markers_res{primary_res}.csv")
    markers_df.to_csv(markers_file)
    print(f"    - Saved markers to: {markers_file}")

    # Visualizations
    
    # 1. Dotplot of Top Markers
    # We need to construct a dict of markers for the dotplot
    try:
        # Get top 3 markers per cluster
        markers_dict = {}
        for group in groups:
            markers_dict[group] = result['names'][group][:3].tolist()
            
        sc.pl.dotplot(adata, markers_dict, groupby=primary_key, dendrogram=False, standard_scale='var', show=False)
        plt.title(f"Top Markers ({primary_key})")
        dotplot_file = os.path.join(output_dir, f"{sample_tag}_step1_markers_dotplot.png")
        plt.savefig(dotplot_file, bbox_inches='tight')
        plt.close()
        print(f"    - Saved Marker Dotplot to: {dotplot_file}")
    except Exception as e:
        print(f"    [WARNING] Failed to plot dotplot: {e}")

    # 2. UMAP
    # Custom Style: legend_loc='on data', etc.
    print("    - Plotting UMAP with custom styles...")
    try:
        sc.pl.umap(adata, color=[primary_key], 
                   size=1, 
                   legend_loc='on data', 
                   legend_fontsize=8, 
                   legend_fontoutline=2, 
                   show=False)
        plt.title(f"UMAP (Leiden {primary_res})")
        umap_file = os.path.join(output_dir, f"{sample_tag}_step1_umap_res{primary_res}.png")
        plt.savefig(umap_file, bbox_inches='tight')
        plt.close()
        print(f"    - Saved UMAP to: {umap_file}")
    except Exception as e:
        print(f"    [WARNING] Failed to plot UMAP: {e}")

    # Spatial Map
    if 'spatial' in adata.obsm:
        spatial_df = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'], index=adata.obs.index)
        spatial_df['cluster'] = adata.obs[primary_key]
        
        plt.figure(figsize=(10, 10))
        sns.scatterplot(data=spatial_df, x='x', y='y', hue='cluster', s=2, linewidth=0, palette='tab20')
        plt.title(f"Spatial Map (Leiden {primary_res})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., markerscale=5)
        plt.axis('equal')
        
        spatial_file = os.path.join(output_dir, f"{sample_tag}_step1_spatial_res{primary_res}.png")
        plt.savefig(spatial_file, bbox_inches='tight')
        plt.close()
        print(f"    - Saved Spatial Map to: {spatial_file}")
        
    return primary_key


def run_step1(config):
    """Main execution for Step 1 (Exploration)"""
    input_dir = config["output_dir"]
    sample_tag = config["sample_tag"]
    # Get exploration parameters with defaults
    expl_params = config.get("exploration", {})
    # Default changed to 100 to match notebook 1_7 usually, or keep 20 if preferred. 
    # User requested gap filling, notebook uses 100.
    neighbor_radius = expl_params.get("neighbor_radius", 100.0)
    run_dispersion = expl_params.get("run_transcript_dispersion", False)
    
    plotting_settings()
    
    # Should use the Output from Step 0 (Formatted) DIRECTLY
    # Check for explicit previous step path
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        input_file_step0 = config['previous_step_adata_path']
        print(f"\n[Step 1-0] Loading Step 0 Data from explicit path: {input_file_step0}")
    else:
        # Fallback to current dir or parent structure
        input_file_step0 = os.path.join(input_dir, f"{sample_tag}.h5ad")
        if not os.path.exists(input_file_step0):
             # Try parent directory structure
             parent_dir = os.path.dirname(input_dir)
             step0_path = os.path.join(parent_dir, "step0_formatting", f"{sample_tag}.h5ad")
             if os.path.exists(step0_path):
                 input_file_step0 = step0_path
    
    if os.path.exists(input_file_step0):
        print(f"  > Reading: {input_file_step0}")
        adata = sc.read_h5ad(input_file_step0)
    else:
        logger.error(f"No Step 0 input found at {input_file_step0}. Stop.")
        return None

    # 1. General Stats (1_1)
    calculate_general_stats(adata, input_dir, sample_tag)
    
    # 1-X. Deep QC (Optional)
    if run_dispersion:
        calculate_transcript_dispersion(adata, input_dir, sample_tag)
    else:
        print("\n[Step 1-X] Deep QC skipped (enable 'run_transcript_dispersion' in config to run)")

    # 2. Clustering & Annotation (1_2) 
    # This generates 'leiden' clusters which can be used for Neighborhood Analysis
    primary_cluster_key = perform_clustering_and_annotation(adata, input_dir, sample_tag, config)
    
    # 3. Neighborhood Analysis (1_7)
    # We pass the newly generated primary cluster key to prioritize it over machine clusters
    # But analyze_neighborhoods logic searches for keys. We can rely on it finding 'leiden_X' or 'graph_clusters'.
    analyze_neighborhoods(adata, input_dir, sample_tag, radius=neighbor_radius)
    
    # 4. Save Final Annotated Object
    # We won't overwrite Step 2 file to avoid confusion, but we could save a 'step3' file if needed.
    # For now, we trust the plots/text reports are the main output.
    
    output_file = os.path.join(input_dir, f"{sample_tag}_step1_exploration.h5ad")
    print(f"\n[Step 1-4] Saving Annotated Data to {output_file}")
    adata.write(output_file)
    
    print("\n=== Step 1 Analysis Complete ===")
    return adata

if __name__ == '__main__':
    import yaml
    print("Xenium Pipeline Step 3: Exploration (Standalone)")
    
    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = yaml.safe_load(f)
        run_step1(config)
    else:
        print(" [ERROR] config.yaml not found.")
