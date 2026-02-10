# Step 2: Segmentation-Free Analysis (Ref: notebooks/2_segmentation_free_analysis/)
# Analyzes transcript spatial patterns without cell boundaries.
#
# Flow:
#   2-1. Points2Regions clustering
#   2-2. Overlaps analysis (ovrlpy, optional)
#   2-3. Distance metrics
#     2-3a. Distance to centroid
#     2-3b. Distance to boundary (EDT)
#   2-4. Centroid distance plots
#     2-4a. Extreme genes stripplot + boxplot
#     2-4b. Color-coded boxplot
#     2-4c. ECDF
#     2-4d. Welch's t-test heatmap
#   2-5. Summary distance histograms
#   2-6. Save final AnnData

import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy import stats
from skimage.draw import polygon
from skimage.filters.rank import entropy as img_entropy
from skimage.morphology import disk

# Points2Regions — pip package (pip install points2regions)
from points2regions import Points2Regions

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# --- 2-3b. Distance to boundary (EDT) — helpers ---

def load_boundaries(input_dir):
    """
    Loads nucleus boundaries from parquet or csv files.
    Returns a dataframe with columns ['cell_id', 'vertex_x', 'vertex_y'].
    """
    candidates = [
        "nucleus_boundaries.parquet",
        "nucleus_boundaries.csv.gz",
        "nucleus_boundaries.csv"
    ]

    boundary_df = None

    for fname in candidates:
        fpath = os.path.join(input_dir, fname)
        if os.path.exists(fpath):
            logger.info(f"Loading boundaries from {fpath}...")
            if fname.endswith('.parquet'):
                boundary_df = pd.read_parquet(fpath)
            else:
                boundary_df = pd.read_csv(fpath)
            break

    if boundary_df is not None:
        return boundary_df
    else:
        logger.warning(f"No boundary files found in {input_dir}. Skipping boundary-based analysis.")
        return None

def calculate_distance_to_boundary(adata, input_dir, pixel_size=0.5, padding=50):
    """
    Calculates distance from each transcript to the nearest nucleus boundary
    using rasterization + Euclidean Distance Transform (EDT).
    """
    logger.info("Starting Distance to Nucleus Analysis (Rasterization Method)...")

    boundaries = load_boundaries(input_dir)
    if boundaries is None:
        return adata

    # Determine grid dimensions from transcript + boundary extents
    spots = adata.uns['spots']
    min_x, max_x = spots['x_location'].min(), spots['x_location'].max()
    min_y, max_y = spots['y_location'].min(), spots['y_location'].max()

    b_min_x, b_max_x = boundaries['vertex_x'].min(), boundaries['vertex_x'].max()
    b_min_y, b_max_y = boundaries['vertex_y'].min(), boundaries['vertex_y'].max()

    min_x = min(min_x, b_min_x) - padding
    min_y = min(min_y, b_min_y) - padding
    max_x = max(max_x, b_max_x) + padding
    max_y = max(max_y, b_max_y) + padding

    width = int((max_x - min_x) / pixel_size) + 1
    height = int((max_y - min_y) / pixel_size) + 1

    logger.info(f"Rasterizing boundaries onto {height}x{width} grid (Pixel Size: {pixel_size} um)...")

    # Rasterize nucleus polygons onto binary mask
    mask = np.zeros((height, width), dtype=bool)

    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    groups = boundaries.groupby('cell_id')

    for cell_id, grp in groups:
        vx = (grp['vertex_x'].values - min_x) / pixel_size
        vy = (grp['vertex_y'].values - min_y) / pixel_size
        rr, cc = polygon(vy, vx, shape=mask.shape)
        mask[rr, cc] = True

    # EDT on inverted mask: background pixels get distance to nearest nucleus
    logger.info("Computing Euclidean Distance Transform...")
    dist_map = distance_transform_edt(~mask)
    dist_map = dist_map * pixel_size

    # Map each transcript location to the distance grid
    logger.info("Mapping transcripts to distance map...")
    px = (spots['x_location'].values - min_x) / pixel_size
    py = (spots['y_location'].values - min_y) / pixel_size

    valid = (px >= 0) & (px < width) & (py >= 0) & (py < height)

    dists = np.full(len(spots), np.nan)
    dists[valid] = map_coordinates(dist_map, [py[valid], px[valid]], order=1)

    adata.uns['spots']['dist_to_nucleus'] = dists

    logger.info("Distance to boundary analysis complete.")
    return adata


# --- 2-3a. Distance to centroid — helper ---

def calculate_distance_to_centroid(adata, input_dir):
    """
    Calculates Euclidean distance from each transcript to its assigned cell's centroid.
    Implements the approach from xb.calculating.dispersion().

    Requires:
        - adata.obs: 'cell_id', 'x_centroid', 'y_centroid'
        - adata.uns['spots']: 'cell_id', 'x_location', 'y_location', 'feature_name'

    Stores:
        - adata.uns['spots']['dist_to_centroid']
        - adata.uns['gene_distance_stats'] (per-gene mean, median, std)
    """
    logger.info("Starting Distance to Centroid Analysis (Notebook 2_1 approach)...")

    spots = adata.uns['spots'].copy()
    cells_metadata = adata.obs

    if 'cell_id' not in spots.columns:
        logger.error("'cell_id' column not found in spots dataframe. Cannot compute centroid distance.")
        return adata

    if 'x_centroid' not in cells_metadata.columns or 'y_centroid' not in cells_metadata.columns:
        logger.error("'x_centroid' and/or 'y_centroid' not found in adata.obs. Cannot compute centroid distance.")
        return adata

    cells_metadata_filt = cells_metadata.loc[
        cells_metadata['cell_id'].isin(spots['cell_id']), :
    ]
    spots = spots[spots['cell_id'].isin(cells_metadata_filt['cell_id'])].copy()

    logger.info(f"Computing centroid distance for {len(spots)} transcripts "
                f"across {len(cells_metadata_filt)} cells...")

    dict_x = dict(zip(cells_metadata_filt['cell_id'], cells_metadata_filt['x_centroid']))
    dict_y = dict(zip(cells_metadata_filt['cell_id'], cells_metadata_filt['y_centroid']))

    spots['x_cell'] = spots['cell_id'].map(dict_x).astype(float)
    spots['y_cell'] = spots['cell_id'].map(dict_y).astype(float)

    spots['dist_to_centroid'] = np.sqrt(
        (spots['x_location'] - spots['x_cell']) ** 2 +
        (spots['y_location'] - spots['y_cell']) ** 2
    )

    # Write back; unassigned transcripts get NaN
    adata.uns['spots']['dist_to_centroid'] = np.nan
    adata.uns['spots'].loc[spots.index, 'dist_to_centroid'] = spots['dist_to_centroid'].values

    # Per-gene distance statistics (excluding control probes)
    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    spots_filtered = spots[gene_mask]

    gene_stats = spots_filtered.groupby('feature_name')['dist_to_centroid'].agg(
        ['mean', 'median', 'std', 'count']
    ).rename(columns={
        'mean': 'mean_distance',
        'median': 'median_distance',
        'std': 'std_distance',
        'count': 'transcript_count'
    })
    gene_stats = gene_stats.sort_values('mean_distance')

    adata.uns['gene_distance_stats'] = gene_stats

    logger.info(f"Distance to centroid analysis complete. "
                f"Gene stats computed for {len(gene_stats)} genes.")
    return adata


# --- 2-4. Centroid distance plots ---

def plot_centroid_distance_analysis(adata, output_dir, sample_tag, n_top_genes=20):
    """
    Creates per-gene centroid distance visualizations (Ref: notebook 2_1):
      1. Stripplot + boxplot for extreme genes (closest/farthest from centroid)
      2. ECDF per gene
      3. Welch's t-test p-value heatmap
      4. Gene-level distance stats CSV
    """
    logger.info("Generating per-gene centroid distance visualizations...")

    spots = adata.uns['spots']

    if 'dist_to_centroid' not in spots.columns:
        logger.warning("dist_to_centroid not found in spots. Skipping centroid distance plots.")
        return

    spots_valid = spots.dropna(subset=['dist_to_centroid']).copy()

    # Exclude control probes
    gene_mask = ~spots_valid['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    spots_valid = spots_valid[gene_mask]

    if spots_valid.empty:
        logger.warning("No valid transcripts with centroid distance. Skipping plots.")
        return

    # 10% subsample for efficiency (as in notebook)
    np.random.seed(0)
    n_subsample = max(1, len(spots_valid) // 10)
    random_indices = np.random.permutation(len(spots_valid))[:n_subsample]
    rp = spots_valid.iloc[random_indices].copy()

    # Select extreme genes: top 5 closest + top 5 farthest among high-count genes
    dist_by_gene = rp.groupby('feature_name')['dist_to_centroid'].mean().sort_values()
    gene_counts = rp.groupby('feature_name')['dist_to_centroid'].count().sort_values(ascending=False)
    high_exp = gene_counts.index[:100]
    dist_by_gene = dist_by_gene[dist_by_gene.index.isin(high_exp)]

    n_extreme = min(5, len(dist_by_gene) // 2)
    if n_extreme < 1:
        logger.warning("Not enough genes for extreme gene analysis. Skipping plots.")
        return

    extreme_genes = list(dist_by_gene.index[:n_extreme]) + list(dist_by_gene.index[-n_extreme:])
    top_genes = list(gene_counts.index[:n_top_genes])

    # Mean cell border distance (from nucleus-overlapping transcripts, if available)
    if 'overlaps_nucleus' in rp.columns:
        rp_nuc = rp[rp['overlaps_nucleus'] == 1]
        if not rp_nuc.empty:
            dist_max = rp_nuc.groupby('cell_id')['dist_to_centroid'].max()
            mean_cellborder = np.mean(dist_max)
        else:
            mean_cellborder = None
    else:
        mean_cellborder = None

    # --- 2-4a. Extreme genes stripplot + boxplot ---
    d_extreme = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_extreme.empty:
        # Balance group sizes
        min_group_size = d_extreme['feature_name'].value_counts().min()
        d_extreme = d_extreme.groupby('feature_name').head(min_group_size).reset_index(drop=True)

        plt.figure(figsize=(12, 6))
        sns.stripplot(
            x=d_extreme['dist_to_centroid'], y=d_extreme['feature_name'],
            s=0.15, order=extreme_genes, jitter=0.4
        )
        sns.boxplot(
            x=d_extreme['dist_to_centroid'], y=d_extreme['feature_name'],
            order=extreme_genes, saturation=0, width=0.3,
            fliersize=0, whis=1, linewidth=0.9, color='white'
        )
        if mean_cellborder is not None:
            plt.axvline(x=mean_cellborder, color='grey', linestyle='--', label='Mean cell border')
        plt.title(f'Distance to Centroid - Extreme Genes ({sample_tag})')
        plt.xlabel('Distance to Centroid (um)')
        plt.ylabel('Gene')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_stripplot.png"), dpi=150)
        plt.close()
        logger.info("Saved centroid distance stripplot.")

        # --- 2-4b. Color-coded boxplot ---
        color_dict = {}
        for g in dist_by_gene.index[:n_extreme]:
            color_dict[g] = 'red'
        for g in dist_by_gene.index[-n_extreme:]:
            color_dict[g] = 'blue'

        plt.figure(figsize=(10, 5))
        sns.boxplot(
            x=d_extreme['dist_to_centroid'], y=d_extreme['feature_name'],
            order=extreme_genes, width=0.6, fliersize=0, whis=1,
            saturation=0.3, linewidth=0.9, palette=color_dict
        )
        if mean_cellborder is not None:
            plt.axvline(x=mean_cellborder, color='grey', linestyle='--', label='Mean cell border')
        plt.title(f'Distance to Centroid - Boxplot ({sample_tag})')
        plt.xlabel('Distance to Centroid (um)')
        plt.ylabel('Gene')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_boxplot.png"), dpi=150)
        plt.close()
        logger.info("Saved centroid distance boxplot.")

    # --- 2-4c. ECDF ---
    d_ecdf = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_ecdf.empty:
        fig = sns.displot(
            x=d_ecdf['dist_to_centroid'], kind='ecdf', hue=d_ecdf['feature_name']
        )
        if mean_cellborder is not None:
            plt.axvline(x=mean_cellborder, color='grey', linestyle='--', label='Mean cell border')
        plt.title(f'ECDF of Distance to Centroid ({sample_tag})')
        plt.xlabel('Distance to Centroid (um)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_ecdf.png"), dpi=150)
        plt.close()
        logger.info("Saved centroid distance ECDF plot.")

    # --- 2-4d. Welch's t-test heatmap ---
    d_ttest = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_ttest.empty and len(extreme_genes) > 1:
        pval_matrix = pd.DataFrame(columns=extreme_genes, index=extreme_genes, dtype=float)

        for gene1 in extreme_genes:
            for gene2 in extreme_genes:
                g1_dists = d_ttest.loc[d_ttest['feature_name'] == gene1, 'dist_to_centroid']
                g2_dists = d_ttest.loc[d_ttest['feature_name'] == gene2, 'dist_to_centroid']
                if len(g1_dists) > 1 and len(g2_dists) > 1:
                    pval_matrix.loc[gene1, gene2] = stats.ttest_ind(
                        g1_dists, g2_dists, equal_var=False
                    ).pvalue
                else:
                    pval_matrix.loc[gene1, gene2] = np.nan

        pval_path = os.path.join(output_dir, f"{sample_tag}_step2_centroid_welchs_pvals.csv")
        pval_matrix.to_csv(pval_path)
        logger.info(f"Saved Welch's t-test p-value matrix to {pval_path}")

        plt.figure(figsize=(10, 8))
        pval_float = pval_matrix.astype(float)
        with np.errstate(divide='ignore'):
            log_pvals = -np.log10(pval_float.values.astype(float))
        log_pvals_df = pd.DataFrame(log_pvals, index=extreme_genes, columns=extreme_genes)
        sns.heatmap(log_pvals_df, annot=True, fmt='.1f', cmap='YlOrRd',
                    square=True, cbar_kws={'label': '-log10(p-value)'})
        plt.title(f"Welch's t-test: -log10(p-value) ({sample_tag})")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_welchs_heatmap.png"), dpi=150)
        plt.close()
        logger.info("Saved Welch's t-test heatmap.")

    # Save gene-level distance statistics
    if 'gene_distance_stats' in adata.uns:
        stats_path = os.path.join(output_dir, f"{sample_tag}_step2_gene_distance_stats.csv")
        adata.uns['gene_distance_stats'].to_csv(stats_path)
        logger.info(f"Saved gene distance statistics to {stats_path}")

    # Save plotting subset for extreme genes
    plotting_subset = spots[
        spots['feature_name'].isin(extreme_genes)
    ][['x_location', 'y_location', 'feature_name', 'dist_to_centroid']].copy()
    subset_path = os.path.join(output_dir, f"{sample_tag}_step2_centroid_plotting_subset.csv")
    plotting_subset.to_csv(subset_path, index=False)
    logger.info(f"Saved plotting subset to {subset_path}")


# --- Entry Point ---

def run_step2(config):
    logger.info("--- Step 2: Segmentation-Free Analysis ---")

    sample_tag = config['sample_tag']
    prev_adata_path = config['previous_step_adata_path']
    input_dir = config['input_path']
    output_dir = config['output_dir']

    distance_metric = config['segmentation_free'].get('distance_metric', 'centroid')

    logger.info(f"Loading data from {prev_adata_path}...")
    adata = sc.read_h5ad(prev_adata_path)

    # Ensure spots are loaded
    if 'spots' not in adata.uns:
        logger.warning("'spots' dataframe not found in adata.uns. Attempting to reload from raw CSV...")
        try:
            transcripts_path = os.path.join(input_dir, "transcripts.csv")
            df = pd.read_csv(transcripts_path)
            if 'feature_name' not in df.columns:
                 if 'gene' in df.columns: df.rename(columns={'gene': 'feature_name'}, inplace=True)

            adata.uns['spots'] = df
        except Exception as e:
            logger.error(f"Failed to load transcripts: {e}")
            raise

    # --- 2-1. Points2Regions ---
    if config['segmentation_free'].get('run_points2regions', True):
        logger.info("Running Points2Regions...")
        spots = adata.uns['spots']
        xy = spots[['x_location', 'y_location']].values
        genes = spots['feature_name'].values

        p2r_params = config['segmentation_free']['points2regions']

        unique_genes = np.unique(genes)
        gene_map = {g: i for i, g in enumerate(unique_genes)}
        gene_labels = np.array([gene_map[g] for g in genes])

        sigma = p2r_params.get('sigma', 2.0)
        p2r = Points2Regions(
            xy, gene_labels,
            pixel_width=sigma / 3.0,
            pixel_smoothing=sigma,
            min_num_pts_per_pixel=p2r_params.get('min_genes_per_bin', 5),
        )
        p2r_adata = p2r.fit_predict(
            num_clusters=p2r_params.get('n_clusters', 10),
            output='anndata',
            seed=42,
            adata_cluster_key='points2regions',
        )

        cluster_col = p2r_adata.uns['reads']['points2regions']
        adata.uns['spots']['points2regions'] = cluster_col.values

        p2r_out_path = os.path.join(output_dir, f"{sample_tag}_step2_points2regions_bins.h5ad")
        p2r_adata.write_h5ad(p2r_out_path)
        logger.info(f"Saved Points2Regions bin data to {p2r_out_path}")

    # --- 2-2. Overlaps Analysis (ovrlpy, optional) ---
    if config['segmentation_free'].get('run_overlaps', False):
        logger.info("Starting Signal Overlaps (Incoherence) Analysis using ovrlpy...")
        try:
            import ovrlpy
            logger.info("ovrlpy imported successfully. Proceeding with analysis...")
            # Placeholder: ovrlpy API call would go here once API is finalized
            pass

        except ImportError as e:
            logger.critical("'ovrlpy' library is required for Overlaps analysis but is not installed.")
            logger.critical("Please install it or disable 'run_overlaps' in config.yaml.")
            raise e

    # --- 2-3. Distance Metrics ---
    logger.info(f"Distance metric mode: '{distance_metric}'")

    if distance_metric in ('centroid', 'both'):
        adata = calculate_distance_to_centroid(adata, input_dir)
        plot_centroid_distance_analysis(adata, output_dir, sample_tag)

    if distance_metric in ('boundary', 'both'):
        adata = calculate_distance_to_boundary(adata, input_dir)

    # --- 2-5. Summary distance histograms ---

    # Boundary-based distance histogram
    if 'dist_to_nucleus' in adata.uns['spots'].columns:
        plt.figure(figsize=(10, 6))
        dists = adata.uns['spots']['dist_to_nucleus']
        dists = dists[~np.isnan(dists)]
        plt.hist(dists, bins=100, log=True)
        plt.title('Distribution of Transcript Distances to Nuclei (Boundary)')
        plt.xlabel('Distance (um)')
        plt.ylabel('Count (Log Scale)')
        plt.axvline(0, color='r', linestyle='--', label='Nucleus Boundary')
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_distance_boundary_dist.png"))
        plt.close()

    # Centroid-based distance histogram
    if 'dist_to_centroid' in adata.uns['spots'].columns:
        plt.figure(figsize=(10, 6))
        dists = adata.uns['spots']['dist_to_centroid']
        dists = dists[~np.isnan(dists)]
        if len(dists) > 0:
            plt.hist(dists, bins=100, log=True)
            plt.title('Distribution of Transcript Distances to Centroid')
            plt.xlabel('Distance to Centroid (um)')
            plt.ylabel('Count (Log Scale)')
            plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_distance_centroid_dist.png"))
            plt.close()

    # --- 2-6. Save final AnnData ---
    output_file = os.path.join(output_dir, f"{sample_tag}_step2_points2regions.h5ad")
    adata.write_h5ad(output_file)
    logger.info(f"Step 2 Completed. Output saved to {output_file}")

    return adata

if __name__ == "__main__":
    pass
