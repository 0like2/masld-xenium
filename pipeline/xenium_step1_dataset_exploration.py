# Step 1: Dataset Exploration (Ref: notebooks/1_datasets_exploration/)
# Calculates dataset statistics, transcript dispersion, clustering, and
# neighborhood architecture metrics.
#
# Flow:
#   1-1. General statistics & summary
#   1-2. Transcript dispersion analysis
#     1-2a. Distance histogram + KDE
#     1-2b. ECDF (overall + per-gene top 10)
#     1-2c. Violin plot (top 20 genes)
#   1-3. KS tests (gene-pair distance distributions)
#   1-4. Clustering & marker annotation
#     1-4a. PCA + neighbors
#     1-4b. HVG selection
#     1-4c. Leiden clustering (multi-resolution)
#     1-4d. Marker genes (Wilcoxon) → dotplot + heatmap
#     1-4e. UMAP + spatial scatter
#   1-5. Neighborhood analysis
#     1-5a. Spatial neighbors graph
#     1-5b. Neighborhood diversity
#     1-5c. Enrichment analysis
#     1-5d. Centrality scores

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def plotting_settings():
    """Applies visualization settings matching the notebooks."""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white', figsize=(10, 10))

    sns.set_style("white")
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    plt.rcParams['figure.facecolor'] = 'white'

def _decode_bytes(df):
    """Decodes byte columns to utf-8 strings."""
    for col in df.columns:
        if df[col].dtype == 'object':
            first_valid = df[col].dropna().iloc[0] if not df[col].dropna().empty else None
            if isinstance(first_valid, bytes):
                print(f"    - Decoding bytes column: {col}")
                df[col] = df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
    return df

def load_transcripts_sidecar(adata_path, sample_tag):
    """Loads transcripts from sidecar parquet/csv if not in adata.uns."""
    if adata_path:
        step0_dir = os.path.dirname(adata_path)
    else:
        step0_dir = "."

    df = None
    # Try parquet first, then legacy formats
    parquet_path = os.path.join(step0_dir, f"{sample_tag}_transcripts.parquet")
    if os.path.exists(parquet_path):
        print(f"    - Loading transcripts from sidecar: {parquet_path}")
        df = pd.read_parquet(parquet_path)

    elif os.path.exists(os.path.join(step0_dir, "transcripts.parquet")):
        parquet_path_simple = os.path.join(step0_dir, "transcripts.parquet")
        df = pd.read_parquet(parquet_path_simple)

    elif os.path.exists(os.path.join(step0_dir, "transcripts.csv")):
        csv_path = os.path.join(step0_dir, "transcripts.csv")
        print(f"    - Loading transcripts from sidecar CSV: {csv_path}")
        df = pd.read_csv(csv_path, low_memory=False)

    if df is not None:
        df = _decode_bytes(df)
        return df

    return None


# --- 1-1. General statistics & summary ---

def calculate_general_stats(adata, output_dir, sample_tag, spots=None):
    """Calculates general dataset statistics (Ref: 1_1)."""
    print("\n[Step 1-1] Calculating General Statistics (Ref: 1_1)...")

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    stats = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_genes_per_cell": np.median(adata.obs['n_genes_by_counts']),
        "median_counts_per_cell": np.median(adata.obs['total_counts']),
        "total_counts": np.sum(adata.obs['total_counts'])
    }

    # Control-probe QC (computed in Step 0)
    if 'pct_counts_control' in adata.obs:
        stats['median_pct_control_reads'] = np.median(adata.obs['pct_counts_control'])
        stats['mean_pct_control_reads'] = np.mean(adata.obs['pct_counts_control'])

    # Spot-level QC from raw transcripts table
    if spots is not None:
        n_total_reads = len(spots)
        stats['total_reads'] = n_total_reads

        if 'qv' in spots.columns:
            prop_qv20 = (spots['qv'] > 20).mean()
            stats['reads_prop_qv>20'] = prop_qv20

        valid_genes = set(adata.var_names)
        if 'feature_name' in spots.columns:
            n_in_panel = spots['feature_name'].isin(valid_genes).sum()
            stats['prop_reads_in_panel'] = n_in_panel / n_total_reads

        if 'cell_id' in spots.columns:
             if 'cell_id' in adata.obs.columns:
                 valid_cells = set(adata.obs['cell_id'])
             else:
                 valid_cells = set(adata.obs_names)

             n_assigned = spots['cell_id'].isin(valid_cells).sum()
             stats['prop_reads_assigned_to_cells'] = n_assigned / n_total_reads

        n_cells_gt_10 = (adata.obs['total_counts'] > 10).sum()
        stats['proportion_cells>10reads'] = n_cells_gt_10 / adata.n_obs
    else:
        print("    [WARNING] No transcripts data available for Spot-level QC stats.")

    print("    - Stats:", stats)
    with open(os.path.join(output_dir, f"{sample_tag}_step1_stats.txt"), "w") as f:
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")

    stats_df = pd.DataFrame([stats])
    stats_csv_path = os.path.join(output_dir, f"{sample_tag}_step1_stats.csv")
    stats_df.to_csv(stats_csv_path, index=False)
    print(f"    - Saved stats CSV to: {stats_csv_path}")

    # Normalized heatmap of stats
    try:
        numeric_stats = stats_df.select_dtypes(include=[np.number])
        if numeric_stats.shape[1] > 1:
            plt.figure(figsize=(max(12, numeric_stats.shape[1] * 1.2), 4))
            norm_stats = numeric_stats.copy()
            for col in norm_stats.columns:
                col_min = norm_stats[col].min()
                col_max = norm_stats[col].max()
                if col_max != col_min:
                    norm_stats[col] = (norm_stats[col] - col_min) / (col_max - col_min)
                else:
                    norm_stats[col] = 1.0
            sns.heatmap(norm_stats, annot=numeric_stats.values, fmt='.4g',
                        cmap='YlOrRd', xticklabels=numeric_stats.columns,
                        yticklabels=[sample_tag], cbar_kws={'label': 'Normalized Value'})
            plt.title(f"Dataset Statistics Summary - {sample_tag}")
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            heatmap_path = os.path.join(output_dir, f"{sample_tag}_step1_stats_heatmap.png")
            plt.savefig(heatmap_path, bbox_inches='tight')
            plt.close()
            print(f"    - Saved stats heatmap to: {heatmap_path}")
    except Exception as e:
        print(f"    [WARNING] Failed to generate stats heatmap: {e}")


# --- 1-2. Transcript dispersion analysis ---

def _compute_distances(adata, spots):
    """
    Computes transcript-to-centroid distances from spots dataframe.
    Returns (spots_assigned, distances, metric_name) or (None, None, None) on failure.
    """
    spots = spots.copy()

    if 'cell_id' not in spots.columns:
        if spots.index.name == 'cell_id':
            spots = spots.reset_index()
        else:
            print("    [ERROR] Cannot link spots to cells (missing 'cell_id'). Skipping.")
            return None, None, None

    if 'cell_id' in adata.obs.columns:
        valid_cells = set(adata.obs['cell_id'])
    else:
        valid_cells = set(adata.obs.index)

    spots_assigned = spots[spots['cell_id'].isin(valid_cells)].copy()

    if len(spots_assigned) == 0:
        print("    [WARNING] No assigned transcripts found matching filtered cells.")
        return None, None, None

    distances = None
    metric_name = "distance_to_centroid"

    # Option A: pre-calculated nucleus distance
    if 'nucleus_distance' in spots_assigned.columns:
        print("    - Found 'nucleus_distance' column. Using pre-calculated values.")
        distances = spots_assigned['nucleus_distance']
        metric_name = "distance_to_nucleus"

    # Option B: Euclidean distance to cell centroid
    elif 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        print("    - Calculating Euclidean distance to Cell Centroid...")
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
        return None, None, None

    if distances is None or len(distances) == 0:
        return None, None, None

    distances = distances.dropna()
    spots_assigned = spots_assigned.loc[distances.index]
    spots_assigned['distance'] = distances.values

    return spots_assigned, distances, metric_name


def calculate_transcript_dispersion(adata, output_dir, sample_tag, spots=None):
    """
    Calculates dispersion metrics for transcripts (Ref: 1_3).
    Prioritizes 'nucleus_distance' if available, otherwise calculates
    Euclidean distance to cell centroid. Generates histogram, ECDF, violin plots.
    """
    print("\n[Step 1-2] Transcript Dispersion Analysis (Ref: 1_3)...")

    if spots is None:
        print("    [WARNING] No transcripts data provided. Skipping dispersion analysis.")
        return

    print("    - Preparing transcript data...")
    spots_assigned, distances, metric_name = _compute_distances(adata, spots)

    if spots_assigned is None or distances is None:
        print("    [WARNING] Could not compute distances. Skipping dispersion analysis.")
        return

    n_total = len(spots)
    n_assigned = len(spots_assigned)
    print(f"    - Analyzing {n_assigned} assigned transcripts (out of {n_total})...")

    mean_dist = np.mean(distances)
    median_dist = np.median(distances)
    print(f"    - Median {metric_name}: {median_dist:.2f}")

    # --- 1-2a. Distance histogram + KDE ---
    plt.figure(figsize=(8, 6))
    sns.histplot(distances, bins=100, kde=True, color='purple')
    plt.title(f"Transcript Dispersion Distribution\nMetric: {metric_name} | Median: {median_dist:.2f}")
    plt.xlabel(f"{metric_name} (pixels/microns)")
    plt.ylabel("Count")

    q99 = np.percentile(distances, 99)
    plt.xlim(0, q99)

    save_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_dist.png")
    plt.savefig(save_path)
    plt.close()
    print(f"    - Saved Dispersion Histogram to {save_path}")

    # --- 1-2b. ECDF (overall + per-gene top 10) ---
    try:
        print("    - Generating ECDF plots...")

        plt.figure(figsize=(10, 6))
        if 'source' in spots_assigned.columns and spots_assigned['source'].nunique() > 1:
            sns.ecdfplot(data=spots_assigned, x='distance', hue='source', complementary=False)
            plt.title("ECDF of Transcript Distance to Centroid (by Source)")
        else:
            sns.ecdfplot(data=spots_assigned, x='distance', complementary=False, color='steelblue')
            plt.title("ECDF of Transcript Distance to Centroid")
        plt.xlabel("Distance (um)")
        plt.ylabel("Cumulative Proportion")
        plt.xlim(0, q99)
        ecdf_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_ecdf.png")
        plt.savefig(ecdf_path, bbox_inches='tight')
        plt.close()
        print(f"    - Saved ECDF plot to {ecdf_path}")

        if 'feature_name' in spots_assigned.columns:
            top10_genes = spots_assigned['feature_name'].value_counts().head(10).index
            spots_top10 = spots_assigned[spots_assigned['feature_name'].isin(top10_genes)]
            if len(spots_top10) > 0:
                plt.figure(figsize=(10, 6))
                sns.ecdfplot(data=spots_top10, x='distance', hue='feature_name', complementary=False)
                plt.title("ECDF of Transcript Distance to Centroid (Top 10 Genes)")
                plt.xlabel("Distance (um)")
                plt.ylabel("Cumulative Proportion")
                plt.xlim(0, q99)
                plt.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
                ecdf_gene_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_ecdf_genes.png")
                plt.savefig(ecdf_gene_path, bbox_inches='tight')
                plt.close()
                print(f"    - Saved gene-level ECDF plot to {ecdf_gene_path}")
    except Exception as e:
        print(f"    [WARNING] ECDF plotting failed: {e}")

    # --- 1-2c. Violin plot (top 20 genes) ---
    try:
        if 'feature_name' in spots_assigned.columns:
            print("    - Generating violin plot for top 20 genes...")
            top20_genes = spots_assigned['feature_name'].value_counts().head(20).index
            spots_top20 = spots_assigned[spots_assigned['feature_name'].isin(top20_genes)]
            if len(spots_top20) > 0:
                plt.figure(figsize=(14, 6))
                sns.violinplot(data=spots_top20, x='feature_name', y='distance',
                               cut=0, scale='width', order=top20_genes)
                plt.xticks(rotation=90)
                plt.title("Distance Distribution per Gene (Top 20)")
                plt.xlabel("Gene")
                plt.ylabel(f"{metric_name} (um)")
                plt.tight_layout()
                violin_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_violin.png")
                plt.savefig(violin_path, bbox_inches='tight')
                plt.close()
                print(f"    - Saved violin plot to {violin_path}")
    except Exception as e:
        print(f"    [WARNING] Violin plot generation failed: {e}")

    # Save dispersion metrics summary
    with open(os.path.join(output_dir, f"{sample_tag}_step1_dispersion_metrics.txt"), "w") as f:
        f.write(f"Metric Used: {metric_name}\n")
        f.write(f"Total Transcripts: {n_total}\n")
        f.write(f"Assigned Transcripts: {n_assigned}\n")
        f.write(f"Median Distance: {median_dist}\n")
        f.write(f"Mean Distance: {mean_dist}\n")


def run_ks_tests(spots_assigned, output_dir, sample_tag, n_top_genes=20):
    """
    Runs pairwise Kolmogorov-Smirnov tests on transcript distance distributions
    between top genes. Saves p-value matrix as CSV and heatmap of -log10(p-values).

    Parameters
    ----------
    spots_assigned : pd.DataFrame
        Assigned transcripts with 'feature_name' and 'distance' columns.
    output_dir : str
        Directory to save outputs.
    sample_tag : str
        Sample identifier for file naming.
    n_top_genes : int
        Number of top genes (by transcript count) to include.
    """
    print(f"\n[Step 1-KS] Running Pairwise KS Tests (top {n_top_genes} genes)...")

    if spots_assigned is None or 'feature_name' not in spots_assigned.columns or 'distance' not in spots_assigned.columns:
        print("    [WARNING] spots_assigned missing required columns ('feature_name', 'distance'). Skipping KS tests.")
        return

    top_genes = spots_assigned['feature_name'].value_counts().head(n_top_genes).index.tolist()
    n_genes = len(top_genes)

    if n_genes < 2:
        print("    [WARNING] Fewer than 2 genes available for KS test. Skipping.")
        return

    print(f"    - Testing {n_genes} genes pairwise ({n_genes * (n_genes - 1) // 2} pairs)...")

    gene_distances = {}
    for gene in top_genes:
        dists = spots_assigned.loc[spots_assigned['feature_name'] == gene, 'distance'].dropna().values
        if len(dists) > 0:
            gene_distances[gene] = dists

    valid_genes = [g for g in top_genes if g in gene_distances]
    n_valid = len(valid_genes)

    if n_valid < 2:
        print("    [WARNING] Fewer than 2 genes with valid distances. Skipping KS tests.")
        return

    pval_matrix = pd.DataFrame(np.ones((n_valid, n_valid)), index=valid_genes, columns=valid_genes)
    stat_matrix = pd.DataFrame(np.zeros((n_valid, n_valid)), index=valid_genes, columns=valid_genes)

    for i in range(n_valid):
        for j in range(i + 1, n_valid):
            gene_a = valid_genes[i]
            gene_b = valid_genes[j]
            ks_stat, p_val = ks_2samp(gene_distances[gene_a], gene_distances[gene_b])
            pval_matrix.loc[gene_a, gene_b] = p_val
            pval_matrix.loc[gene_b, gene_a] = p_val
            stat_matrix.loc[gene_a, gene_b] = ks_stat
            stat_matrix.loc[gene_b, gene_a] = ks_stat

    pval_csv_path = os.path.join(output_dir, f"{sample_tag}_step1_ks_pvalues.csv")
    pval_matrix.to_csv(pval_csv_path)
    print(f"    - Saved KS p-value matrix to: {pval_csv_path}")

    stat_csv_path = os.path.join(output_dir, f"{sample_tag}_step1_ks_statistics.csv")
    stat_matrix.to_csv(stat_csv_path)
    print(f"    - Saved KS statistic matrix to: {stat_csv_path}")

    try:
        pval_clipped = pval_matrix.clip(lower=1e-300)
        neglog10_pval = -np.log10(pval_clipped)
        np.fill_diagonal(neglog10_pval.values, 0)

        plt.figure(figsize=(max(10, n_valid * 0.6), max(8, n_valid * 0.5)))
        sns.heatmap(neglog10_pval, cmap='YlOrRd', square=True,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={'label': '-log10(p-value)'})
        plt.title(f"Pairwise KS Test: -log10(p-value)\n(Top {n_valid} Genes by Transcript Count)")
        plt.xticks(rotation=90, fontsize=8)
        plt.yticks(fontsize=8)
        plt.tight_layout()

        heatmap_path = os.path.join(output_dir, f"{sample_tag}_step1_ks_heatmap.png")
        plt.savefig(heatmap_path, bbox_inches='tight')
        plt.close()
        print(f"    - Saved KS heatmap to: {heatmap_path}")
    except Exception as e:
        print(f"    [WARNING] KS heatmap generation failed: {e}")


# --- 1-4. Clustering & marker annotation ---

def perform_clustering_and_annotation(adata, output_dir, sample_tag, config):
    """
    Performs Leiden clustering, HVG selection, and marker gene ranking (Ref: 1_2).
    """
    print("\n[Step 1-3] Running Clustering & Annotation (Ref: 1_2)...")

    # --- 1-4a. PCA + neighbors ---
    if 'X_pca' not in adata.obsm:
        print("    - PCA not found. Running PCA...")
        sc.pp.pca(adata)
    if 'neighbors' not in adata.uns:
        print("    - Neighbors not found. Computing Neighbors...")
        sc.pp.neighbors(adata)

    # --- 1-4b. HVG selection ---
    print("    - Calculating Highly Variable Genes...")
    try:
        if 'log1p' not in adata.uns:
             sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)

        sc.pl.highly_variable_genes(adata, show=False)
        plt.title("Highly Variable Genes")
        hvg_file = os.path.join(output_dir, f"{sample_tag}_step1_hvg.png")
        plt.savefig(hvg_file)
        plt.close()
        print(f"    - Saved HVG Plot to: {hvg_file}")
    except Exception as e:
        print(f"    [WARNING] HVG calculation/plotting failed: {e}")

    # --- 1-4c. Leiden clustering (multi-resolution) ---
    resolutions = config.get("exploration", {}).get("resolutions",
                   config.get("annotation", {}).get("resolutions", [0.5, 0.8, 1.0]))

    for res in resolutions:
        key = f"leiden_{res}"
        print(f"    - Running Leiden clustering (resolution={res})...")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        print(f"      - Found {len(adata.obs[key].unique())} clusters.")

    # Select primary resolution for downstream analyses
    primary_res = config.get("exploration", {}).get("primary_resolution",
                    config.get("annotation", {}).get("primary_resolution", 1.0))
    primary_key = f"leiden_{primary_res}"
    if primary_key not in adata.obs:
         print(f"    - Primary resolution {primary_res} missing. Computing...")
         sc.tl.leiden(adata, resolution=primary_res, key_added=primary_key)

    print(f"    - Using '{primary_key}' as primary clustering.")

    # --- 1-4d. Marker genes (Wilcoxon) → dotplot + heatmap ---
    print("    - Ranking marker genes...")
    sc.tl.rank_genes_groups(adata, groupby=primary_key, method='wilcoxon')

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    markers_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']}
    ).head(5)

    markers_file = os.path.join(output_dir, f"{sample_tag}_step1_markers_res{primary_res}.csv")
    markers_df.to_csv(markers_file)
    print(f"    - Saved markers to: {markers_file}")

    # Dotplot of top 3 markers per cluster
    try:
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

    # Heatmap of top marker genes
    try:
        print("    - Generating marker gene ranking heatmap...")
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby=primary_key, show=False)
        heatmap_file = os.path.join(output_dir, f"{sample_tag}_step1_markers_heatmap.png")
        plt.savefig(heatmap_file, bbox_inches='tight')
        plt.close()
        print(f"    - Saved Marker Heatmap to: {heatmap_file}")
    except Exception as e:
        print(f"    [WARNING] Failed to plot marker heatmap: {e}")

    # --- 1-4e. UMAP + spatial scatter ---
    print("    - Plotting UMAP...")
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

    # Spatial scatter map
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


# --- 1-5. Neighborhood analysis ---

def analyze_neighborhoods(adata, output_dir, sample_tag, radius=20.0):
    """Calculates neighborhood metrics: diversity, density, enrichment, centrality (Ref: 1_7)."""
    print("\n[Step 1-4] Analyzing Neighborhood Architecture (Ref: 1_7)...")

    if 'spatial' not in adata.obsm:
        print("    [WARNING] No spatial coordinates found. Skipping neighborhood analysis.")
        return

    # --- 1-5a. Spatial neighbors graph ---
    print(f"    - Computing spatial neighbors graph (radius={radius})...")
    sq.gr.spatial_neighbors(adata, coord_type="generic", radius=radius)

    # --- 1-5b. Neighborhood diversity ---
    if 'spatial_connectivities' in adata.obsp:
        print("    - Calculating custom neighborhood metrics (diversity, density)...")
        # Density = degree (number of neighbors per cell)
        adata.obs['neighborhood_density'] = np.array(adata.obsp['spatial_connectivities'].sum(axis=1)).flatten()

        # Diversity requires cluster labels
        cluster_key = None
        for key in ['leiden_1.0', 'leiden_0.8', 'leiden', 'graph_clusters']:
            if key in adata.obs:
                cluster_key = key
                break

        if cluster_key:
            # Diversity = count of distinct cluster labels in each cell's neighborhood
            # Computed via: (Adjacency * OneHot) > 0, then sum per row
            dummies = pd.get_dummies(adata.obs[cluster_key])
            adj = adata.obsp['spatial_connectivities']

            from scipy import sparse
            if not sparse.issparse(adj):
                adj = sparse.csr_matrix(adj)

            cluster_counts = adj.dot(sparse.csr_matrix(dummies.values))
            adata.obs['neighborhood_diversity'] = np.array((cluster_counts > 0).sum(axis=1)).flatten()

            print(f"      - Diversity calculated using '{cluster_key}'")

    # --- 1-5c. Enrichment analysis ---
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

        plt.figure(figsize=(8, 8))
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_key)
        plt.title(f"Neighborhood Enrichment ({cluster_key})")
        save_path = os.path.join(output_dir, f"{sample_tag}_step1_grad_enrichment.png")
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"    - Saved enrichment plot to {save_path}")

        # --- 1-5d. Centrality scores ---
        print(f"    - Computing centrality scores...")
        sq.gr.centrality_scores(adata, cluster_key=cluster_key)

        plt.figure(figsize=(10, 5))
        sq.pl.centrality_scores(adata, cluster_key=cluster_key)
        plt.title("Centrality Scores")
        save_path_cent = os.path.join(output_dir, f"{sample_tag}_step1_centrality.png")
        plt.savefig(save_path_cent, bbox_inches='tight')
        plt.close()
    else:
        print("    [WARNING] No cluster key found. Skipping enrichment/centrality.")


# --- Main Entry Point ---

def run_step1(config):
    """Main execution for Step 1 (Exploration)."""
    input_dir = config["output_dir"]
    sample_tag = config["sample_tag"]
    expl_params = config.get("exploration", {})
    neighbor_radius = expl_params.get("neighbor_radius", 100.0)
    run_dispersion = expl_params.get("run_transcript_dispersion", True)

    plotting_settings()

    # Load Step 0 output
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        input_file_step0 = config['previous_step_adata_path']
        print(f"\n[Step 1-0] Loading Step 0 Data from explicit path: {input_file_step0}")
    else:
        input_file_step0 = os.path.join(input_dir, f"{sample_tag}.h5ad")
        if not os.path.exists(input_file_step0):
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

    # Load transcripts sidecar
    spots = None
    if 'spots' in adata.uns:
        print("  > Found 'spots' in adata.uns (Legacy format).")
        spots = adata.uns['spots']
    else:
        spots = load_transcripts_sidecar(input_file_step0, sample_tag)

    # --- 1-1. General statistics ---
    calculate_general_stats(adata, input_dir, sample_tag, spots=spots)

    # --- 1-2. Transcript dispersion ---
    if run_dispersion:
        calculate_transcript_dispersion(adata, input_dir, sample_tag, spots=spots)

        # --- 1-3. KS tests ---
        if spots is not None:
            spots_assigned, distances, metric_name = _compute_distances(adata, spots)
            if spots_assigned is not None and 'distance' in spots_assigned.columns:
                run_ks_tests(spots_assigned, input_dir, sample_tag, n_top_genes=20)
            else:
                print("    [INFO] Skipping KS tests - could not compute distances.")
    else:
        print("\n[Step 1-2] Dispersion analysis skipped (enable 'run_transcript_dispersion' in config)")

    # --- 1-4. Clustering & annotation ---
    primary_cluster_key = perform_clustering_and_annotation(adata, input_dir, sample_tag, config)

    # --- 1-5. Neighborhood analysis ---
    analyze_neighborhoods(adata, input_dir, sample_tag, radius=neighbor_radius)

    # --- 5. Save ---
    output_file = os.path.join(input_dir, f"{sample_tag}_step1_exploration.h5ad")
    print(f"\n[Step 1-5] Saving Annotated Data to {output_file}")
    adata.write(output_file)

    print("\n=== Step 1 Analysis Complete ===")
    return adata

if __name__ == '__main__':
    import yaml
    print("Xenium Pipeline Step 1: Exploration (Standalone)")

    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = yaml.safe_load(f)
        run_step1(config)
    else:
        print(" [ERROR] config.yaml not found.")
