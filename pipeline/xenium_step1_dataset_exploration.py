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

def _decode_bytes(df):
    """Decodes byte columns to utf-8 strings."""
    for col in df.columns:
        if df[col].dtype == 'object':
            # Check first non-null element to see if it is bytes
            first_valid = df[col].dropna().iloc[0] if not df[col].dropna().empty else None
            if isinstance(first_valid, bytes):
                print(f"    - Decoding bytes column: {col}")
                # Decode bytes to string
                df[col] = df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
    return df

def load_transcripts_sidecar(adata_path, sample_tag):
    """Loads transcripts from sidecar parquet/csv if not in adata.uns"""
    # Try Step 0 directory first (standard path)
    if adata_path:
        step0_dir = os.path.dirname(adata_path)
    else:
        # Fallback to current dir ? 
        step0_dir = "."

    df = None
    # Try Parquet first (new standard)
    parquet_path = os.path.join(step0_dir, f"{sample_tag}_transcripts.parquet")
    if os.path.exists(parquet_path):
        print(f"    - Loading transcripts from sidecar: {parquet_path}")
        df = pd.read_parquet(parquet_path)
    
    # Check for legacy standard name
    elif os.path.exists(os.path.join(step0_dir, "transcripts.parquet")):
        parquet_path_simple = os.path.join(step0_dir, "transcripts.parquet")
        df = pd.read_parquet(parquet_path_simple)
        
    # Check for legacy CSV
    elif os.path.exists(os.path.join(step0_dir, "transcripts.csv")):
        csv_path = os.path.join(step0_dir, "transcripts.csv")
        print(f"    - Loading transcripts from sidecar CSV: {csv_path}")
        df = pd.read_csv(csv_path, low_memory=False)
        
    if df is not None:
        # Decode bytes if present
        df = _decode_bytes(df)
        return df
        
    return None

def calculate_general_stats(adata, output_dir, sample_tag, spots=None):
    """Calculates general dataset statistics (Ref: 1_1)..."""
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
    
    # [NEW] Add Spot-level QC metrics from metadata if available (calculated in Step 0)
    # Step 0 computes total_counts_control, total_counts_raw
    if 'pct_counts_control' in adata.obs:
        stats['median_pct_control_reads'] = np.median(adata.obs['pct_counts_control'])
        stats['mean_pct_control_reads'] = np.mean(adata.obs['pct_counts_control'])
    
    # [NEW] Add Spot-level QC metrics from raw spots table if available
    if spots is not None:
        n_total_reads = len(spots)
        stats['total_reads'] = n_total_reads
        
        # Quality Value > 20
        if 'qv' in spots.columns:
            prop_qv20 = (spots['qv'] > 20).mean()
            stats['reads_prop_qv>20'] = prop_qv20
        
        # Reads in Panel (already filtered controls in Step 0)
        # Check against valid genes in filtered adata
        valid_genes = set(adata.var_names)
        if 'feature_name' in spots.columns:
            # Count how many of the *raw* spots match the *filtered* gene list
            n_in_panel = spots['feature_name'].isin(valid_genes).sum()
            stats['prop_reads_in_panel'] = n_in_panel / n_total_reads
            
        # Assigned to cells
        if 'cell_id' in spots.columns:
             # Check if adata has 'cell_id' column (preferred) or use index
             if 'cell_id' in adata.obs.columns:
                 valid_cells = set(adata.obs['cell_id'])
             else:
                 valid_cells = set(adata.obs_names)
                 
             n_assigned = spots['cell_id'].isin(valid_cells).sum()
             stats['prop_reads_assigned_to_cells'] = n_assigned / n_total_reads
        
        # Proportion cells > 10 reads        
        # Proportion cells > 10 reads
        n_cells_gt_10 = (adata.obs['total_counts'] > 10).sum()
        stats['proportion_cells>10reads'] = n_cells_gt_10 / adata.n_obs
    else:
        print("    [WARNING] No transcripts data available for Spot-level QC stats.")

    # Save statistics
    print("    - Stats:", stats)
    with open(os.path.join(output_dir, f"{sample_tag}_step1_stats.txt"), "w") as f:
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")


def calculate_transcript_dispersion(adata, output_dir, sample_tag, spots=None):
    """
    Calculates dispersion metrics for transcripts (Ref: 1_3).
    Prioritizes 'nucleus_distance' if available, otherwise calculates Euclidean distance to cell centroid.
    """
    print("\n[Step 1-X] Deep QC: Transcript Dispersion Analysis (Ref: 1_3)...")

    if spots is None:
        print("    [WARNING] No transcripts data provided. Skipping dispersion analysis.")
        return

    # 1. Prepare Data
    print("    - Preparing transcript data...")
    # Don't copy if not needed, but safe
    spots = spots.copy()
    
    # Needs cell info
    if 'x_centroid' not in adata.obs.columns or 'y_centroid' not in adata.obs.columns:
        print("    [WARNING] Cell centroids (x_centroid, y_centroid) missing in adata.obs. Skipping Euclidean calculation.")
    
    # Usually 'cell_id' is a column in transcripts.csv
    if 'cell_id' not in spots.columns:
        if spots.index.name == 'cell_id':
            spots = spots.reset_index()
        else:
             print("    [ERROR] Cannot link spots to cells (missing 'cell_id'). Skipping.")
             return
             
    # Filter for assigned transcripts only
    # Match cell_id to adata.obs.index OR adata.obs['cell_id']
    if 'cell_id' in adata.obs.columns:
        valid_cells = set(adata.obs['cell_id'])
    else:
        valid_cells = set(adata.obs.index)

    spots_assigned = spots[spots['cell_id'].isin(valid_cells)]
    
    n_total = len(spots)
    n_assigned = len(spots_assigned)

    if n_assigned == 0:
        print(f"    [WARNING] No assigned transcripts found matching filtered cells. Skipping dispersion.")
        return

    print(f"    - Analyzing {n_assigned} assigned transcripts (out of {n_total})...")
    
    # 2. Logic Branch: Metric Selection
    distances = None
    metric_name = "distance_to_centroid"
    
    # Option A: Pre-calculated Nucleus Distance
    if 'nucleus_distance' in spots_assigned.columns:
        print("    - Found 'nucleus_distance' column. Using pre-calculated values.")
        distances = spots_assigned['nucleus_distance']
        metric_name = "distance_to_nucleus"
        
    # Option B: Calculate Euclidean Distance to Centroid
    elif 'x_centroid' in adata.obs.columns:
        print("    - 'nucleus_distance' not found. Calculating Euclidean distance to Cell Centroid...")
        
        # Merge cell coordinates
        if 'cell_id' in adata.obs.columns:
            right_on_key = 'cell_id'
            use_index = False
        else:
             right_on_key = None
             use_index = True
             
        merged = spots_assigned.merge(
            adata.obs[['x_centroid', 'y_centroid'] + ([right_on_key] if right_on_key else [])], 
            left_on='cell_id', 
            right_on=right_on_key,
            right_index=use_index, 
            how='left'
        )
    
        dx = merged['x_location'] - merged['x_centroid']
        dy = merged['y_location'] - merged['y_centroid']
        distances = np.sqrt(dx**2 + dy**2)
    
    else:
        print("    [ERROR] Missing required columns for dispersion calculation.")
        return
        
    # 3. Stats & Plotting
    if distances is None or len(distances) == 0:
        try:
            # Fallback if distances is None but somehow we got here
            pass
        except:
             return
        print("    [WARNING] No valid distances computed.")
        return
        
    # Handle NaNs
    distances = distances.dropna()
    
    mean_dist = np.mean(distances)
    median_dist = np.median(distances)
    print(f"    - Median {metric_name}: {median_dist:.2f}")
    
    # Plot Distribution
    plt.figure(figsize=(8, 6))
    sns.histplot(distances, bins=100, kde=True, color='purple')
    plt.title(f"Transcript Dispersion Distribution\nMetric: {metric_name} | Median: {median_dist:.2f}")
    plt.xlabel(f"{metric_name} (pixels/microns)")
    plt.ylabel("Count")
    
    # Optional: Clip outliers for better visualization
    q99 = np.percentile(distances, 99)
    plt.xlim(0, q99)
    
    save_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_dist.png")
    plt.savefig(save_path)
    plt.close()
    print(f"    - Saved Dispersion Plot to {save_path}")
    
    # Save metrics to a file
    with open(os.path.join(output_dir, f"{sample_tag}_step1_dispersion_metrics.txt"), "w") as f:
        f.write(f"Metric Used: {metric_name}\n")
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
    
    if input_file_step0 and os.path.exists(input_file_step0):
        print(f"  > Reading: {input_file_step0}")
        adata = sc.read_h5ad(input_file_step0)
    else:
        logger.error(f"No Step 0 input found at {input_file_step0} or path check failed. Stop.")
        return None

    # Load Transcripts Sidecar
    spots = None
    if 'spots' in adata.uns:
        # Backward compatibility or if user forced it back
        print("  > Found 'spots' in adata.uns (Legacy format).")
        spots = adata.uns['spots']
    else:
        spots = load_transcripts_sidecar(input_file_step0, sample_tag)

    # 1. General Stats (1_1)
    # Pass spots explicitely
    calculate_general_stats(adata, input_dir, sample_tag, spots=spots)
    
    # 1-X. Deep QC (Optional)
    if run_dispersion:
        calculate_transcript_dispersion(adata, input_dir, sample_tag, spots=spots)
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
