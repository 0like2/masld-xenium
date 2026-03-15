# Step 2: Segmentation-Free Analysis (Ref: notebooks/2_segmentation_free_analysis/)
# Analyzes transcript spatial patterns without cell boundaries.
#
# Flow:
#   2-1.  Points2Regions clustering
#   2-1b. P2R subcellular classification (nuclear/cyto)          [NEW]
#   2-1b2.P2R spatial map (Fig 1i / Ext Data Fig 3A)            [NEW]
#   2-1c. P2R cluster x cell type heatmap                        [NEW]
#   2-1d. P2R top genes per subcellular cluster (line + bar)     [NEW]
#   2-1e. P2R DE marker genes + violin + ECDF + proportion      [NEW]
#   2-2.  Overlaps analysis (ovrlpy) + ROI vertical cuts        [NEW]
#   2-3.  Distance metrics
#     2-3a. Distance to centroid
#     2-3b. Distance to boundary (EDT, SIGNED)                   [UPDATED]
#   2-4.  Centroid distance plots
#     2-4a. Extreme genes stripplot + boxplot
#     2-4b. Color-coded boxplot
#     2-4c. ECDF
#     2-4d. Welch's t-test heatmap
#   2-5.  P2R boundary distance box plots per cell type          [NEW]
#   2-6.  Summary distance histograms
#   2-7.  SSAM analysis (optional)                               [NEW]
#   2-8.  Save final AnnData

import os
import logging
import warnings
import colorsys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy import stats
from scipy.cluster import hierarchy
from skimage.draw import polygon

# Points2Regions -- pip package (pip install points2regions)
from points2regions import Points2Regions

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Publication-quality figure settings
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# ============================================================
# Boundary helpers
# ============================================================

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

    for fname in candidates:
        fpath = os.path.join(input_dir, fname)
        if os.path.exists(fpath):
            logger.info(f"Loading boundaries from {fpath}...")
            if fname.endswith('.parquet'):
                return pd.read_parquet(fpath)
            else:
                return pd.read_csv(fpath)

    logger.warning(f"No boundary files found in {input_dir}. Skipping boundary-based analysis.")
    return None


# ============================================================
# 2-1b. P2R Subcellular Classification
# Ref: DEPRECATED_2_2, points2regions/compute_colors.ipynb
# ============================================================

def classify_p2r_subcellular(adata, output_dir, sample_tag):
    """
    Classifies Points2Regions clusters as nuclear or cytoplasmic
    based on the fraction of transcripts overlapping nuclei.

    Creates:
        - adata.uns['spots']['p2r_name']:  "{cluster}_{celltype}_{nuclei|cyto}"
        - adata.uns['spots']['p2r_compartment']:  "nuclei" or "cyto"
        - adata.uns['spots']['p2r_celltype']:  majority cell type per cluster
        - adata.uns['p2r_classification']:  DataFrame with cluster metadata
    """
    logger.info("Classifying P2R clusters as nuclear/cytoplasmic...")

    spots = adata.uns['spots']

    if 'points2regions' not in spots.columns:
        logger.warning("P2R clusters not found in spots. Skipping subcellular classification.")
        return adata

    if 'overlaps_nucleus' not in spots.columns:
        logger.warning("'overlaps_nucleus' not in spots. Skipping subcellular classification.")
        return adata

    # Resolve cell type column (from spots or adata.obs)
    celltype_col = _resolve_celltype_column(spots, adata)

    # Per-cluster statistics
    p2r_groups = spots.groupby('points2regions')

    cluster_stats = []
    for cluster_id, group in p2r_groups:
        nuc_frac = group['overlaps_nucleus'].mean()
        compartment = 'nuclei' if nuc_frac > 0.5 else 'cyto'

        # Majority cell type (prefer non-background)
        ct_counts = group[celltype_col].value_counts()
        valid_cts = ct_counts.drop(
            ['Background', 'Unknown', 'unassigned'], errors='ignore'
        )
        majority_ct = (valid_cts.index[0] if len(valid_cts) > 0
                       else (ct_counts.index[0] if len(ct_counts) > 0 else 'Unknown'))

        cluster_stats.append({
            'cluster_id': cluster_id,
            'majority_celltype': majority_ct,
            'compartment': compartment,
            'nuc_fraction': nuc_frac,
            'n_transcripts': len(group),
            'p2r_name': f"{cluster_id}_{majority_ct}_{compartment}"
        })

    classification_df = pd.DataFrame(cluster_stats)
    adata.uns['p2r_classification'] = classification_df

    # Map back to spots
    name_map = dict(zip(classification_df['cluster_id'], classification_df['p2r_name']))
    comp_map = dict(zip(classification_df['cluster_id'], classification_df['compartment']))
    ct_map = dict(zip(classification_df['cluster_id'], classification_df['majority_celltype']))

    spots['p2r_name'] = spots['points2regions'].map(name_map)
    spots['p2r_compartment'] = spots['points2regions'].map(comp_map)
    spots['p2r_celltype'] = spots['points2regions'].map(ct_map)
    adata.uns['spots'] = spots

    # Save
    class_path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_classification.csv")
    classification_df.to_csv(class_path, index=False)

    n_nuc = (classification_df['compartment'] == 'nuclei').sum()
    n_cyto = (classification_df['compartment'] == 'cyto').sum()
    logger.info(f"P2R subcellular: {n_nuc} nuclear, {n_cyto} cytoplasmic clusters")

    return adata


def _normalize_cell_ids(series):
    """Normalize cell IDs to consistent int-strings (e.g. '1.0' → '1', '2' → '2')."""
    s = series.astype(str)
    # Strip trailing '.0' from float-string representations
    return s.str.replace(r'\.0$', '', regex=True)


def _resolve_celltype_column(spots, adata):
    """Find or create a cell type column in spots.

    Handles type mismatches between spots['cell_id'] (often int/float) and
    adata.obs index (often str) by normalizing both sides to int-strings.
    """
    # Build candidate column list — include any leiden_* variants
    base_cols = ['Class', 'celltype', 'cell_type', 'leiden']
    leiden_variants = sorted([c for c in list(spots.columns) + list(adata.obs.columns)
                              if c.startswith('leiden_')], key=lambda x: x)
    candidate_cols = base_cols + leiden_variants

    for col in candidate_cols:
        if col in spots.columns:
            n_valid = spots[col].notna().sum()
            n_unknown = spots[col].isin(['Unknown', 'Background', 'unassigned']).sum()
            if (n_valid - n_unknown) > 0:
                return col

    # Try mapping from adata.obs via cell_id
    if 'cell_id' in spots.columns:
        spots_cid = _normalize_cell_ids(spots['cell_id'])
        for col in candidate_cols:
            if col in adata.obs.columns:
                # Build mapping — normalize keys to int-strings to match spots_cid
                if 'cell_id' in adata.obs.columns:
                    mapping = dict(zip(_normalize_cell_ids(adata.obs['cell_id']), adata.obs[col]))
                else:
                    mapping = dict(zip(_normalize_cell_ids(adata.obs.index.to_series()), adata.obs[col]))
                spots['_celltype_mapped'] = spots_cid.map(mapping).fillna('Background')
                n_mapped = (spots['_celltype_mapped'] != 'Background').sum()
                logger.info(f"Mapped {n_mapped}/{len(spots)} transcripts to cell types via '{col}' "
                            f"(mapping keys: {len(mapping)}, spots cell_ids: {spots_cid.nunique()})")
                if n_mapped > 0:
                    return '_celltype_mapped'

    # Fallback: spatial nearest-neighbor mapping
    # If cell_id mapping fails (e.g. unassigned transcripts), use spatial proximity
    if 'x_location' in spots.columns and 'y_location' in spots.columns:
        for col in candidate_cols:
            if col in adata.obs.columns and 'spatial' in adata.obsm:
                try:
                    from scipy.spatial import cKDTree
                    cell_coords = adata.obsm['spatial']
                    tree = cKDTree(cell_coords)
                    spot_coords = spots[['x_location', 'y_location']].values
                    dists, idxs = tree.query(spot_coords, k=1)
                    labels = adata.obs[col].values[idxs]
                    # Only assign if within reasonable distance (e.g., 50µm)
                    max_dist = 50.0
                    labels[dists > max_dist] = 'Background'
                    spots['_celltype_mapped'] = labels
                    n_mapped = (spots['_celltype_mapped'] != 'Background').sum()
                    logger.info(f"Spatial nearest-cell mapping: {n_mapped}/{len(spots)} "
                                f"transcripts within {max_dist}µm via '{col}'")
                    if n_mapped > 0:
                        return '_celltype_mapped'
                except Exception as e:
                    logger.warning(f"Spatial nearest-cell mapping failed: {e}")

    logger.warning("Could not resolve cell type column — all transcripts labelled 'Unknown'")
    spots['_celltype_mapped'] = 'Unknown'
    return '_celltype_mapped'


def assign_p2r_colors(adata):
    """
    Assigns HLS-based colors to P2R clusters, varying lightness
    by nuclear proximity within each cell type group.
    Ref: points2regions/compute_colors.ipynb
    """
    if 'p2r_classification' not in adata.uns:
        return {}

    classification_df = adata.uns['p2r_classification']

    # Base palette per cell type
    unique_cts = classification_df['majority_celltype'].unique()
    base_palette = sns.color_palette('husl', n_colors=max(len(unique_cts), 1))
    ct_colors = dict(zip(unique_cts, [matplotlib.colors.rgb2hex(c) for c in base_palette]))

    color_palette = {}
    for ct in unique_cts:
        ct_clusters = classification_df[classification_df['majority_celltype'] == ct].copy()
        ct_clusters = ct_clusters.sort_values('nuc_fraction', ascending=False)
        n = len(ct_clusters)

        base_hex = ct_colors[ct].lstrip('#')
        base_rgb = tuple(int(base_hex[i:i+2], 16) / 255.0 for i in (0, 2, 4))
        base_hls = colorsys.rgb_to_hls(*base_rgb)

        for idx, (_, row) in enumerate(ct_clusters.iterrows()):
            if n == 1:
                lightness = 0.5
                hue_shift = 0
            else:
                lightness = 0.4 + 0.4 * (idx / (n - 1))
                hue_shift = 0.2 * (idx / (n - 1)) - 0.1

            hue = (base_hls[0] + hue_shift) % 1.0
            rgb = colorsys.hls_to_rgb(hue, lightness, base_hls[2])
            hex_color = '#%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            color_palette[row['p2r_name']] = hex_color

    adata.uns['p2r_colors'] = color_palette
    return color_palette


# ============================================================
# 2-1b2. P2R Spatial Map
# Ref: Fig 1i (celltype-specific P2R clusters on tissue),
#      Extended Data Fig 3A (full tissue P2R cluster map + ROI)
# ============================================================

def plot_p2r_spatial_map(adata, output_dir, sample_tag):
    """Spatial scatter of transcripts coloured by P2R cluster/celltype.

    Generates:
    1. Full tissue view coloured by P2R celltype (Fig 1i style)
    2. Zoomed 500µm × 500µm ROI coloured by P2R cluster (Ext Data Fig 3A style)
    """
    logger.info("Generating P2R spatial map (Fig 1i / Ext Data Fig 3A)...")

    spots = adata.uns.get('spots')
    if spots is None or 'points2regions' not in spots.columns:
        logger.warning("No P2R cluster data in spots. Skipping P2R spatial map.")
        return

    if 'x_location' not in spots.columns or 'y_location' not in spots.columns:
        logger.warning("No spatial coordinates in spots. Skipping P2R spatial map.")
        return

    color_palette = adata.uns.get('p2r_colors', {})

    # --- Figure 1: Full tissue coloured by P2R celltype ---
    try:
        fig, axes = plt.subplots(1, 2, figsize=(22, 10))

        # Subsample for performance
        n_tx = min(len(spots), 500_000)
        sub = spots.sample(n=n_tx, random_state=42) if len(spots) > n_tx else spots

        # Panel 1: coloured by P2R celltype (Fig 1i style)
        ct_col = 'p2r_celltype' if 'p2r_celltype' in sub.columns else None
        if ct_col:
            unique_cts = sub[ct_col].dropna().unique()
            ct_cmap = plt.cm.get_cmap('tab20', max(len(unique_cts), 1))
            ct_color_map = {ct: ct_cmap(i) for i, ct in enumerate(sorted(unique_cts))}

            for ct in sorted(unique_cts):
                mask = sub[ct_col] == ct
                ct_sub = sub[mask]
                axes[0].scatter(ct_sub['x_location'], ct_sub['y_location'],
                                s=0.05, alpha=0.3, c=[ct_color_map[ct]],
                                label=ct, rasterized=True)

            # Legend with larger markers
            handles, labels = axes[0].get_legend_handles_labels()
            legend = axes[0].legend(handles, labels, loc='upper right',
                                    fontsize=6, markerscale=20, ncol=2,
                                    framealpha=0.8, title='Cell Type',
                                    title_fontsize=8)
            axes[0].set_title(f'P2R cell types ({n_tx:,} transcripts)', fontsize=12)
        else:
            axes[0].scatter(sub['x_location'], sub['y_location'],
                            s=0.05, alpha=0.2, c='steelblue', rasterized=True)
            axes[0].set_title(f'Transcripts ({n_tx:,})', fontsize=12)

        axes[0].invert_yaxis()
        axes[0].set_aspect('equal')
        axes[0].axis('off')

        # Panel 2: coloured by P2R cluster with p2r_name colors
        if 'p2r_name' in sub.columns and color_palette:
            for name in sorted(sub['p2r_name'].dropna().unique()):
                mask = sub['p2r_name'] == name
                c = color_palette.get(name, '#999999')
                axes[1].scatter(sub.loc[mask, 'x_location'],
                                sub.loc[mask, 'y_location'],
                                s=0.05, alpha=0.3, c=c, rasterized=True)
            axes[1].set_title(f'P2R clusters ({sub["p2r_name"].nunique()} clusters)',
                              fontsize=12)
        else:
            # Fallback: colour by cluster int
            clusters = sub['points2regions']
            unique_k = sorted(clusters.dropna().unique())
            k_cmap = plt.cm.get_cmap('tab20', max(len(unique_k), 1))
            for i, k in enumerate(unique_k):
                mask = clusters == k
                axes[1].scatter(sub.loc[mask, 'x_location'],
                                sub.loc[mask, 'y_location'],
                                s=0.05, alpha=0.3, c=[k_cmap(i % 20)],
                                rasterized=True)
            axes[1].set_title(f'P2R clusters ({len(unique_k)} clusters)', fontsize=12)

        axes[1].invert_yaxis()
        axes[1].set_aspect('equal')
        axes[1].axis('off')

        fig.suptitle('P2R Spatial Map — Full Tissue', fontsize=14, fontweight='bold')
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_spatial_map.png")
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved P2R spatial map: {path}")
    except Exception as e:
        logger.warning(f"P2R spatial map (full tissue) failed: {e}")
        plt.close('all')

    # --- Figure 2: Zoomed ROI (500 µm × 500 µm) ---
    try:
        x_all = spots['x_location'].values
        y_all = spots['y_location'].values
        x_med, y_med = np.nanmedian(x_all), np.nanmedian(y_all)

        roi_half = 250  # 500 µm × 500 µm
        x_min, x_max = x_med - roi_half, x_med + roi_half
        y_min, y_max = y_med - roi_half, y_med + roi_half

        roi_mask = ((x_all >= x_min) & (x_all <= x_max) &
                    (y_all >= y_min) & (y_all <= y_max))
        roi_spots = spots[roi_mask]

        if len(roi_spots) < 10:
            logger.warning("Too few spots in ROI for P2R spatial zoom.")
            return

        fig, axes = plt.subplots(1, 2, figsize=(18, 8))

        # Panel 1: ROI coloured by P2R celltype
        ct_col = 'p2r_celltype' if 'p2r_celltype' in roi_spots.columns else None
        if ct_col:
            unique_cts = roi_spots[ct_col].dropna().unique()
            ct_cmap = plt.cm.get_cmap('tab20', max(len(unique_cts), 1))
            ct_color_map = {ct: ct_cmap(i) for i, ct in enumerate(sorted(unique_cts))}

            for ct in sorted(unique_cts):
                mask = roi_spots[ct_col] == ct
                ct_sub = roi_spots[mask]
                axes[0].scatter(ct_sub['x_location'], ct_sub['y_location'],
                                s=1.5, alpha=0.5, c=[ct_color_map[ct]],
                                label=ct, rasterized=True)

            axes[0].legend(loc='upper right', fontsize=6, markerscale=8,
                           ncol=2, framealpha=0.8, title='Cell Type',
                           title_fontsize=8)
        else:
            axes[0].scatter(roi_spots['x_location'], roi_spots['y_location'],
                            s=1.5, alpha=0.4, c='steelblue', rasterized=True)

        axes[0].set_title(f'P2R cell types — ROI ({len(roi_spots):,} transcripts)',
                          fontsize=11)
        axes[0].set_xlim(x_min, x_max)
        axes[0].set_ylim(y_max, y_min)
        axes[0].set_aspect('equal')
        axes[0].axis('off')

        # Panel 2: ROI coloured by P2R cluster
        if 'p2r_name' in roi_spots.columns and color_palette:
            for name in sorted(roi_spots['p2r_name'].dropna().unique()):
                mask = roi_spots['p2r_name'] == name
                c = color_palette.get(name, '#999999')
                axes[1].scatter(roi_spots.loc[mask, 'x_location'],
                                roi_spots.loc[mask, 'y_location'],
                                s=1.5, alpha=0.5, c=c, rasterized=True)
        else:
            clusters = roi_spots['points2regions']
            unique_k = sorted(clusters.dropna().unique())
            k_cmap = plt.cm.get_cmap('tab20', max(len(unique_k), 1))
            for i, k in enumerate(unique_k):
                mask = clusters == k
                axes[1].scatter(roi_spots.loc[mask, 'x_location'],
                                roi_spots.loc[mask, 'y_location'],
                                s=1.5, alpha=0.5, c=[k_cmap(i % 20)],
                                rasterized=True)

        axes[1].set_title(f'P2R clusters — ROI ({roi_spots["points2regions"].nunique()} clusters)',
                          fontsize=11)
        axes[1].set_xlim(x_min, x_max)
        axes[1].set_ylim(y_max, y_min)
        axes[1].set_aspect('equal')
        axes[1].axis('off')

        fig.suptitle(f'P2R Spatial Map — ROI (500 µm × 500 µm)', fontsize=14,
                     fontweight='bold')
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_spatial_zoomed.png")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved P2R spatial zoomed ROI: {path}")
    except Exception as e:
        logger.warning(f"P2R spatial map (zoomed ROI) failed: {e}")
        plt.close('all')


# ============================================================
# 2-1c. P2R Cluster x Cell Type Heatmap
# Ref: points2regions/figures.ipynb, Extended Data Fig. 3b
# ============================================================

def plot_p2r_heatmap(adata, output_dir, sample_tag):
    """Heatmap of P2R cluster distribution per cell type (normalized by row)."""
    logger.info("Generating P2R x cell type distribution heatmap...")

    spots = adata.uns['spots']
    if 'p2r_name' not in spots.columns:
        logger.warning("P2R classification not found. Skipping heatmap.")
        return

    celltype_col = _resolve_celltype_column(spots, adata)
    pivot = pd.crosstab(spots['p2r_name'], spots[celltype_col])
    if pivot.empty:
        return

    normalized = pivot.div(pivot.sum(axis=1), axis=0)

    # Hierarchical clustering for rows
    if len(normalized) > 2:
        linkage = hierarchy.linkage(normalized.values, method='average', metric='euclidean')
        sorted_indices = hierarchy.leaves_list(linkage)
        sorted_table = normalized.iloc[sorted_indices]
    else:
        sorted_table = normalized

    # Sort columns by position of their max in row order
    val = sorted_table.values - 0.1
    max_row_per_col = np.argmax(val, axis=0)
    max_row_per_col[val.max(axis=0) < 0] = len(sorted_table)
    col_order = np.argsort(max_row_per_col)
    sorted_table = sorted_table.iloc[:, col_order]

    # Row colors
    color_palette = adata.uns.get('p2r_colors', {})
    row_colors = [color_palette.get(n, '#888888') for n in sorted_table.index]

    fig_h = max(8, len(sorted_table) * 0.22)
    fig_w = max(10, len(sorted_table.columns) * 0.45)

    try:
        g = sns.clustermap(
            sorted_table, row_cluster=True, col_cluster=False,
            method='average', metric='euclidean',
            cmap='viridis', figsize=(fig_w, fig_h),
            row_colors=row_colors, yticklabels=1,
            cbar_pos=(0.02, 0.8, 0.03, 0.15)
        )

        # Mark nuclear clusters with circles
        classification_df = adata.uns.get('p2r_classification', pd.DataFrame())
        nuclear_names = set(
            classification_df[classification_df['compartment'] == 'nuclei']['p2r_name']
        )
        for i, row_label in enumerate(sorted_table.index):
            if row_label in nuclear_names:
                max_col = sorted_table.iloc[i].argmax()
                g.ax_heatmap.scatter(
                    [max_col + 0.5], [i + 0.5],
                    marker='o', s=30, linewidth=1, facecolors='none', edgecolors='k'
                )

        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=7)
        g.ax_heatmap.set_title(f'P2R Cluster Distribution per Cell Type ({sample_tag})')

        path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_celltype_heatmap.png")
        g.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved P2R heatmap to {path}")
    except Exception as e:
        logger.warning(f"P2R heatmap failed: {e}")


# ============================================================
# 2-1d. P2R Top Genes per Subcellular Cluster
# Ref: points2regions/figures.ipynb, Extended Data Fig. 3d
# ============================================================

def plot_p2r_top_genes(adata, output_dir, sample_tag, n_top_genes=7):
    """Top differentially expressed genes per P2R subcellular cluster."""
    logger.info("Generating P2R top genes per subcellular cluster...")

    spots = adata.uns['spots']
    if 'p2r_name' not in spots.columns:
        logger.warning("P2R classification not found. Skipping top genes.")
        return

    classification_df = adata.uns.get('p2r_classification', pd.DataFrame())
    if classification_df.empty:
        return

    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    spots_filt = spots[gene_mask].copy()

    # Cell types with >= 3 P2R clusters (interesting for subcellular comparison)
    ct_counts = classification_df.groupby('majority_celltype').size()
    interesting_cts = ct_counts[ct_counts >= 3].index.tolist()
    if not interesting_cts:
        interesting_cts = ct_counts.nlargest(3).index.tolist()

    for selected_ct in interesting_cts[:3]:
        ct_clusters = classification_df[classification_df['majority_celltype'] == selected_ct]
        ct_names = ct_clusters['p2r_name'].tolist()
        subset = spots_filt[spots_filt['p2r_name'].isin(ct_names)]
        if subset.empty:
            continue

        expression = pd.crosstab(subset['feature_name'], subset['p2r_name'])
        min_expr = max(200, expression.values.max() * 0.01)
        expr_filt = expression[expression.max(axis=1) > min_expr]
        if expr_filt.empty or len(expr_filt) < 3:
            continue

        # Normalize per-cluster then per-gene
        expr_norm = expr_filt.div(expr_filt.sum(axis=0), axis=1)
        expr_row_norm = expr_norm.div(expr_norm.sum(axis=1), axis=0)

        # Sort columns: nuclei first (high nuc_fraction) then cyto
        col_order = sorted(
            expr_row_norm.columns,
            key=lambda x: ct_clusters.loc[
                ct_clusters['p2r_name'] == x, 'nuc_fraction'
            ].values[0] if len(ct_clusters[ct_clusters['p2r_name'] == x]) > 0 else 0,
            reverse=True
        )
        expr_row_norm = expr_row_norm[col_order]

        # --- Line plot: top genes per cluster (grid layout for readability) ---
        n_cols_data = len(expr_row_norm.columns)
        max_per_row = 10  # max subplots per row to prevent overly wide figures
        n_grid_cols = min(n_cols_data, max_per_row)
        n_grid_rows = int(np.ceil(n_cols_data / n_grid_cols))
        fig, axes = plt.subplots(n_grid_rows, n_grid_cols,
                                 figsize=(1.8 * n_grid_cols, 4 * n_grid_rows), sharey=True)
        axes_flat = np.atleast_1d(axes).flatten()

        for i, cluster_name in enumerate(expr_row_norm.columns):
            top = expr_row_norm[cluster_name].nlargest(n_top_genes)
            axes_flat[i].plot(range(len(top)), top.values)
            short_name = cluster_name.replace(selected_ct, selected_ct[:8])
            axes_flat[i].set_title(short_name, fontsize=7)
            axes_flat[i].set_xticks(range(len(top)))
            axes_flat[i].set_xticklabels(top.index, rotation=90, fontsize=6)
            axes_flat[i].set_xlim([-0.5, n_top_genes - 0.5])
        # Hide unused subplots
        for j in range(n_cols_data, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle(f'Top Genes per P2R Cluster - {selected_ct} ({sample_tag})', fontsize=10)
        fig.tight_layout()

        safe_ct = selected_ct.replace(' ', '_').replace('/', '-')
        fig.savefig(
            os.path.join(output_dir, f"{sample_tag}_step2_p2r_topgenes_{safe_ct}.png"),
            dpi=150, bbox_inches='tight'
        )
        plt.close(fig)

        # --- Bar chart: top genes per cluster (paper Fig 1k style) ---
        try:
            n_bar_clusters = min(8, n_cols_data)  # limit panels for readability
            fig_bar, axes_bar = plt.subplots(1, n_bar_clusters,
                                             figsize=(2.5 * n_bar_clusters, 4))
            axes_bar = np.atleast_1d(axes_bar)
            color_palette = adata.uns.get('p2r_colors', {})
            for i, cluster_name in enumerate(expr_row_norm.columns[:n_bar_clusters]):
                top = expr_row_norm[cluster_name].nlargest(n_top_genes)
                bar_color = color_palette.get(cluster_name, '#4477AA')
                axes_bar[i].barh(range(len(top)), top.values[::-1], color=bar_color,
                                 edgecolor='white', linewidth=0.3)
                axes_bar[i].set_yticks(range(len(top)))
                axes_bar[i].set_yticklabels(top.index[::-1], fontsize=6)
                short_name = cluster_name.replace(selected_ct, selected_ct[:8])
                axes_bar[i].set_title(short_name, fontsize=7)
                axes_bar[i].set_xlim(0, 1)
            fig_bar.suptitle(f'Top Genes (bar) — {selected_ct} ({sample_tag})', fontsize=10)
            fig_bar.tight_layout()
            fig_bar.savefig(
                os.path.join(output_dir,
                             f"{sample_tag}_step2_p2r_topgenes_bar_{safe_ct}.png"),
                dpi=150, bbox_inches='tight'
            )
            plt.close(fig_bar)
        except Exception as e:
            logger.warning(f"P2R top genes bar chart for {selected_ct} failed: {e}")

        # --- Gene expression heatmap ---
        max_col = expr_row_norm.idxmax(axis=1)
        expr_sorted = expr_row_norm.iloc[max_col.argsort()]
        fig_h = max(6, len(expr_sorted) * 0.2)

        try:
            heatmap_w = min(30, max(6, n_cols_data * 1.5))  # cap width at 30 inches
            g = sns.clustermap(
                expr_sorted, row_cluster=False, col_cluster=False,
                figsize=(heatmap_w, fig_h), cmap='viridis'
            )
            plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
            g.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_p2r_gene_heatmap_{safe_ct}.png"),
                dpi=150, bbox_inches='tight'
            )
            plt.close()
        except Exception as e:
            logger.warning(f"P2R gene heatmap for {selected_ct} failed: {e}")

        logger.info(f"Saved P2R top genes for {selected_ct}")


# ============================================================
# 2-1e. P2R DE Marker Genes + Compartment Violin
# Ref: Paper Fig. 1f, Extended Data Fig. 3
# ============================================================

def extract_p2r_de_markers(adata, output_dir, sample_tag, n_top=15,
                           pval_thresh=0.05, lfc_thresh=0.5):
    """
    Wilcoxon rank-sum DE on pseudo-bulk P2R clusters (nuclear vs cytoplasmic).

    Produces:
      - TSV of marker genes per compartment
      - Grouped violin plot of marker gene distance-to-centroid
      - ECDF of distance-to-centroid aggregated by compartment class
    """
    logger.info("Extracting P2R DE marker genes (Wilcoxon rank-sum)...")

    spots = adata.uns['spots']
    classification_df = adata.uns.get('p2r_classification', pd.DataFrame())
    if classification_df.empty:
        logger.warning("No P2R classification found. Skipping DE markers.")
        return

    if 'p2r_name' not in spots.columns or 'p2r_compartment' not in spots.columns:
        logger.warning("P2R compartment labels not in spots. Skipping DE markers.")
        return

    # Need at least 2 clusters per compartment for DE
    comp_counts = classification_df['compartment'].value_counts()
    if comp_counts.get('nuclei', 0) < 2 or comp_counts.get('cyto', 0) < 2:
        logger.warning(
            f"Need >=2 clusters per compartment for DE (nuclei={comp_counts.get('nuclei', 0)}, "
            f"cyto={comp_counts.get('cyto', 0)}). Skipping."
        )
        return

    # --- 1. Build pseudo-bulk AnnData (gene x p2r_cluster) ---
    gene_mask = ~spots['feature_name'].str.contains(
        'BLANK|NegControl', case=False, na=False
    )
    spots_filt = spots[gene_mask]

    crosstab = pd.crosstab(spots_filt['feature_name'], spots_filt['p2r_name'])
    # Each column is a P2R cluster; transpose so clusters are observations
    X = crosstab.values.T.astype(np.float32)
    cluster_names = crosstab.columns.tolist()
    gene_names = crosstab.index.tolist()

    # Map cluster name -> compartment
    name_to_comp = dict(zip(classification_df['p2r_name'], classification_df['compartment']))
    obs_df = pd.DataFrame({
        'p2r_name': cluster_names,
        'compartment': [name_to_comp.get(n, 'unknown') for n in cluster_names]
    }, index=cluster_names)

    # Drop clusters with unknown compartment
    valid_mask = obs_df['compartment'].isin(['nuclei', 'cyto'])
    obs_df = obs_df[valid_mask]
    X = X[valid_mask.values]

    p2r_adata = sc.AnnData(
        X=X,
        obs=obs_df,
        var=pd.DataFrame(index=gene_names)
    )

    # Normalize (library-size + log1p) for DE
    sc.pp.normalize_total(p2r_adata, target_sum=1e4)
    sc.pp.log1p(p2r_adata)

    # --- 2. Wilcoxon rank-sum DE ---
    sc.tl.rank_genes_groups(p2r_adata, groupby='compartment', method='wilcoxon')

    marker_rows = []
    for group in ['nuclei', 'cyto']:
        result = sc.get.rank_genes_groups_df(p2r_adata, group=group)
        # Strict filter first
        sig = result[
            (result['pvals_adj'] < pval_thresh) &
            (result['logfoldchanges'] > lfc_thresh)
        ].head(n_top)
        if sig.empty:
            # Fallback: score-ranked with lfc filter only (small n → padj underpowered)
            sig = result[result['logfoldchanges'] > lfc_thresh].head(n_top)
        for _, row in sig.iterrows():
            marker_rows.append({
                'gene': row['names'],
                'class': 'nuclear' if group == 'nuclei' else 'cytoplasmic',
                'logfoldchange': round(row['logfoldchanges'], 4),
                'pval_adj': row['pvals_adj'],
                'score': round(row['scores'], 4),
            })

    if not marker_rows:
        logger.warning("No DE markers found (lfc>0.5). Skipping plots.")
        return

    # Note if padj filtering was bypassed
    padj_sig = sum(1 for r in marker_rows if r['pval_adj'] < pval_thresh)
    if padj_sig == 0:
        logger.info(
            f"No genes passed padj<{pval_thresh} (small n={len(p2r_adata)}). "
            f"Using top {len(marker_rows)} by score+lfc instead."
        )

    markers_df = pd.DataFrame(marker_rows)
    tsv_path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_de_markers.tsv")
    markers_df.to_csv(tsv_path, sep='\t', index=False)
    logger.info(f"Saved {len(markers_df)} DE markers to {tsv_path}")

    marker_genes = markers_df['gene'].tolist()
    gene_to_class = dict(zip(markers_df['gene'], markers_df['class']))

    # --- 3. Resolve distance column ---
    dist_col = None
    for candidate in ['dist_to_centroid', 'dist_to_nucleus']:
        if candidate in spots.columns:
            dist_col = candidate
            break
    if dist_col is None:
        logger.warning("No distance column found. Skipping violin/ECDF plots.")
        return

    # --- 4. Compartment-grouped violin plot ---
    plot_spots = spots_filt[
        spots_filt['feature_name'].isin(marker_genes) &
        spots_filt['p2r_compartment'].isin(['nuclei', 'cyto'])
    ].copy()
    plot_spots['distance'] = plot_spots[dist_col].abs()

    if len(plot_spots) < 10:
        logger.warning("Too few transcripts for violin plot. Skipping.")
        return

    # Sort genes by median distance within their compartment class
    plot_spots['marker_class'] = plot_spots['feature_name'].map(gene_to_class)
    gene_order = (
        plot_spots.groupby('feature_name')['distance']
        .median()
        .sort_values()
        .index.tolist()
    )

    n_genes = len(gene_order)
    fig_w = max(8, n_genes * 0.9)
    fig, ax = plt.subplots(figsize=(fig_w, 5))

    palette = {'nuclei': '#4477AA', 'cyto': '#EE7733'}
    sns.violinplot(
        data=plot_spots, x='feature_name', y='distance',
        hue='p2r_compartment', order=gene_order,
        palette=palette, split=False, inner='quartile',
        density_norm='width', cut=0, linewidth=0.5, ax=ax
    )

    # Threshold guidelines (paper Fig. 1f)
    ax.axhline(5, color='grey', linestyle='--', linewidth=0.8, alpha=0.6, label='5 µm')
    ax.axhline(10, color='grey', linestyle=':', linewidth=0.8, alpha=0.6, label='10 µm')

    ax.set_xlabel('')
    ax.set_ylabel(f'Distance to centroid (µm)')
    ax.set_title(f'P2R DE Marker Genes — Nuclear vs Cytoplasmic ({sample_tag})')
    ax.tick_params(axis='x', rotation=45)
    plt.setp(ax.get_xticklabels(), ha='right', fontsize=7)
    ax.legend(title='Compartment', fontsize=8, title_fontsize=8)
    fig.tight_layout()
    violin_path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_de_marker_violin.png")
    fig.savefig(violin_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DE marker violin to {violin_path}")

    # --- 5. Summary ECDF by compartment ---
    fig, ax = plt.subplots(figsize=(6, 4))

    nuc_genes = [g for g, c in gene_to_class.items() if c == 'nuclear']
    cyto_genes = [g for g, c in gene_to_class.items() if c == 'cytoplasmic']

    for genes, label, color in [
        (nuc_genes, 'Nuclear markers', '#4477AA'),
        (cyto_genes, 'Cytoplasmic markers', '#EE7733'),
    ]:
        dists = spots_filt.loc[
            spots_filt['feature_name'].isin(genes), dist_col
        ].abs().dropna()
        if len(dists) == 0:
            continue
        sorted_d = np.sort(dists.values)
        ecdf_y = np.arange(1, len(sorted_d) + 1) / len(sorted_d)
        ax.plot(sorted_d, ecdf_y, label=f'{label} (n={len(genes)})', color=color, linewidth=1.5)

    ax.axvline(5, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(10, color='grey', linestyle=':', linewidth=0.8, alpha=0.6)
    ax.set_xlabel(f'Distance to centroid (µm)')
    ax.set_ylabel('ECDF')
    ax.set_title(f'DE Marker Distance ECDF — Nuclear vs Cytoplasmic ({sample_tag})')
    ax.legend(fontsize=8)
    ax.set_xlim(left=0)
    fig.tight_layout()
    ecdf_path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_de_marker_ecdf.png")
    fig.savefig(ecdf_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DE marker ECDF to {ecdf_path}")

    # --- 6. Stacked area chart: gene proportion per P2R cluster (Ext Data Fig 3D) ---
    try:
        _plot_p2r_de_proportion_area(spots_filt, classification_df, markers_df,
                                     output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R DE proportion area chart failed: {e}")


def _plot_p2r_de_proportion_area(spots_filt, classification_df, markers_df,
                                  output_dir, sample_tag):
    """Stacked area chart of DE gene proportions per P2R cluster.

    Matches paper Extended Data Fig 3D: each cluster on x-axis, stacked
    proportions of top DE genes on y-axis, coloured by gene.
    """
    marker_genes = markers_df['gene'].unique().tolist()
    if len(marker_genes) < 2:
        return

    # Build cluster × gene proportion matrix
    cluster_gene = pd.crosstab(
        spots_filt[spots_filt['feature_name'].isin(marker_genes)]['p2r_name'],
        spots_filt[spots_filt['feature_name'].isin(marker_genes)]['feature_name']
    )
    if cluster_gene.empty:
        return

    # Normalize rows to proportions
    proportions = cluster_gene.div(cluster_gene.sum(axis=1), axis=0)

    # Sort clusters by compartment (nuclei first) then nuc_fraction
    name_to_comp = dict(zip(classification_df['p2r_name'], classification_df['compartment']))
    name_to_nuc = dict(zip(classification_df['p2r_name'], classification_df['nuc_fraction']))
    valid_idx = [i for i in proportions.index if i in name_to_comp]
    proportions = proportions.loc[valid_idx]
    sort_order = sorted(range(len(proportions)),
                        key=lambda i: (-1 if name_to_comp.get(proportions.index[i]) == 'nuclei' else 1,
                                       -name_to_nuc.get(proportions.index[i], 0)))
    proportions = proportions.iloc[sort_order]

    # Limit to top 15 genes by total count for readability
    top_genes = proportions.sum().nlargest(15).index.tolist()
    proportions = proportions[top_genes]

    fig, ax = plt.subplots(figsize=(max(8, len(proportions) * 0.5), 5))
    x = np.arange(len(proportions))
    colors = plt.cm.get_cmap('tab20', len(top_genes))
    bottom = np.zeros(len(proportions))

    for i, gene in enumerate(top_genes):
        vals = proportions[gene].values
        ax.bar(x, vals, bottom=bottom, width=0.8, label=gene,
               color=colors(i), edgecolor='white', linewidth=0.3)
        bottom += vals

    # Compartment separator line
    n_nuclei = sum(1 for n in proportions.index if name_to_comp.get(n) == 'nuclei')
    if 0 < n_nuclei < len(proportions):
        ax.axvline(n_nuclei - 0.5, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
        ax.text(n_nuclei / 2, 1.02, 'Nuclear', ha='center', fontsize=8,
                transform=ax.get_xaxis_transform())
        ax.text((n_nuclei + len(proportions)) / 2, 1.02, 'Cytoplasmic',
                ha='center', fontsize=8, transform=ax.get_xaxis_transform())

    ax.set_xticks(x)
    ax.set_xticklabels([n.split('_')[0] for n in proportions.index],
                       rotation=90, fontsize=7)
    ax.set_xlabel('P2R Cluster')
    ax.set_ylabel('Gene Proportion')
    ax.set_title(f'DE Gene Proportions per P2R Cluster ({sample_tag})\n'
                 f'(Ext Data Fig 3D style)', fontsize=11)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=6, ncol=3, loc='upper right', bbox_to_anchor=(1.3, 1.0),
              title='Gene', title_fontsize=7)
    fig.tight_layout()
    path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_de_proportion.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved P2R DE proportion area chart: {path}")


# ============================================================
# 2-2. Overlaps Analysis (ovrlpy) -- z-axis signal coherence
# Ref: notebooks/2_3_brain_cell_overlaps.ipynb
# Paper: Fig 1e, Extended Data Fig 2b-e
# ============================================================

def run_ovrlpy_analysis(adata, config, output_dir, sample_tag):
    """
    Runs ovrlpy v1.0+ signal coherence (VSI) analysis to detect overlapping
    cells using z-axis transcript information.

    The analysis:
      1. Builds KDE vector fields split by z-coordinate (top/bottom)
      2. Computes vertical signal integrity (VSI = cosine similarity)
      3. Detects doublet / overlap candidates
      4. Generates signal integrity map and pseudocell plots
    """
    logger.info("Starting ovrlpy signal coherence analysis...")

    try:
        import ovrlpy
    except ImportError:
        logger.warning(
            "ovrlpy not installed. Install with: pip install ovrlpy\n"
            "Skipping overlaps analysis."
        )
        return adata

    spots = adata.uns['spots']

    if 'z_location' not in spots.columns:
        logger.warning("z_location not found in spots. Skipping ovrlpy analysis.")
        return adata

    ovrlpy_cfg = config['segmentation_free'].get('overlaps', {})
    um_per_pixel = ovrlpy_cfg.get('um_per_pixel', 2.0)
    bw = ovrlpy_cfg.get('bw', 1)
    radius = ovrlpy_cfg.get('radius', 50.0)

    # Prepare transcript dataframe (ovrlpy expects x, y, z, gene columns)
    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    coords = spots[gene_mask][['x_location', 'y_location', 'z_location', 'feature_name']].copy()
    coords.columns = ['x', 'y', 'z', 'gene']
    coords['x'] = coords['x'] / um_per_pixel
    coords['y'] = coords['y'] / um_per_pixel
    coords['gene'] = coords['gene'].astype('category')
    coords = coords.reset_index(drop=True)

    logger.info(f"Running ovrlpy on {len(coords)} transcripts...")

    try:
        # --- Initialize Ovrlp object ---
        ovrlp_obj = ovrlpy.Ovrlp(
            coords,
            KDE_bandwidth=bw,
            min_distance=8,
            n_components=30,
            gene_key='gene',
            coordinate_keys=('x', 'y', 'z'),
        )

        # --- Run full analysis: coordinate processing + KDE + VSI ---
        logger.info("Running ovrlpy.analyse (coordinate processing + KDE + VSI)...")
        ovrlp_obj.analyse(gridsize=1, min_transcripts=10, fit_umap=True)

        # --- Signal integrity (coherence) map ---
        path_integrity = os.path.join(output_dir, f"{sample_tag}_step2_coherence_map.png")
        fig = ovrlpy.plot_signal_integrity(ovrlp_obj, signal_threshold=2)
        fig.suptitle(f'Vertical Signal Integrity ({sample_tag})', y=1.02)
        fig.savefig(path_integrity, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved signal integrity map to {path_integrity}")

        # --- Pseudocell plot ---
        try:
            path_pseudo = os.path.join(output_dir, f"{sample_tag}_step2_ovrlpy_pseudocells.png")
            fig_ps = ovrlpy.plot_pseudocells(ovrlp_obj)
            fig_ps.savefig(path_pseudo, dpi=150, bbox_inches='tight')
            plt.close(fig_ps)
            logger.info(f"Saved pseudocell plot to {path_pseudo}")
        except Exception as e:
            logger.warning(f"Pseudocell plot failed: {e}")

        # --- Tissue overview ---
        try:
            path_tissue = os.path.join(output_dir, f"{sample_tag}_step2_ovrlpy_tissue.png")
            fig_t, ax_t = plt.subplots(figsize=(12, 10))
            ovrlpy.plot_tissue(ovrlp_obj, ax=ax_t)
            fig_t.savefig(path_tissue, dpi=150, bbox_inches='tight')
            plt.close(fig_t)
            logger.info(f"Saved tissue plot to {path_tissue}")
        except Exception as e:
            logger.warning(f"Tissue plot failed: {e}")

        # --- ROI vertical cuts: z-axis coherence cross-sections (Fig 1e top) ---
        try:
            _plot_ovrlpy_roi_vertical_cuts(
                ovrlp_obj, coords, um_per_pixel, output_dir, sample_tag
            )
        except Exception as e:
            logger.warning(f"ovrlpy ROI vertical cuts failed: {e}")

        # --- Detect doublets (overlap candidates) ---
        try:
            doublets = ovrlp_obj.detect_doublets(
                min_distance=10, min_integrity=0.7, min_signal=3
            )
            if doublets is not None and len(doublets) > 0:
                doublets_pd = doublets.to_pandas() if hasattr(doublets, 'to_pandas') else doublets
                doublets_pd.to_csv(
                    os.path.join(output_dir, f"{sample_tag}_step2_ovrlpy_doublets.csv"),
                    index=False
                )
                logger.info(f"Detected {len(doublets)} potential doublet/overlap regions")
            else:
                logger.info("No doublets detected")
                doublets = None
        except Exception as e:
            logger.warning(f"Doublet detection failed: {e}")
            doublets = None

        # --- Per-cell incoherence ---
        if hasattr(ovrlp_obj, 'integrity_map') and hasattr(ovrlp_obj, 'signal_map'):
            _compute_and_save_cell_incoherence(
                adata, config, ovrlp_obj.integrity_map, ovrlp_obj.signal_map,
                um_per_pixel, output_dir, sample_tag
            )

        # --- Summary stats ---
        n_doublets = len(doublets) if doublets is not None else 0
        integrity = ovrlp_obj.integrity_map if hasattr(ovrlp_obj, 'integrity_map') else None
        signal = ovrlp_obj.signal_map if hasattr(ovrlp_obj, 'signal_map') else None

        mean_integrity = np.nan
        if integrity is not None and signal is not None:
            high_signal = signal > 2
            if high_signal.any():
                mean_integrity = float(np.nanmean(integrity[high_signal]))

        adata.uns['ovrlpy_coherence'] = {
            'n_doublets': n_doublets,
            'mean_integrity': mean_integrity,
        }

    except Exception as e:
        logger.error(f"ovrlpy analysis failed: {e}")
        import traceback
        traceback.print_exc()

    return adata


def _compute_and_save_cell_incoherence(
    adata, config, integrity_map, signal_map, um_per_pixel, output_dir, sample_tag
):
    """Per-cell incoherence metrics using nucleus boundaries and ovrlpy integrity map."""
    try:
        input_dir = config['input_path']
        boundaries = load_boundaries(input_dir)
        if boundaries is None:
            return

        if 'cell_id' not in adata.obs.columns:
            return

        logger.info("Computing per-cell incoherence metrics...")

        cell_meta = adata.obs
        cell_ids = cell_meta['cell_id'].unique()

        # Limit for performance (full dataset can be very slow)
        max_cells = min(len(cell_ids), 2000)
        cell_ids = cell_ids[:max_cells]

        # Incoherence = 1 - integrity
        incoherence_map = 1.0 - integrity_map

        results = []
        for cid in cell_ids:
            cell_bounds = boundaries[boundaries['cell_id'] == cid]
            if len(cell_bounds) < 3:
                continue

            vx = cell_bounds['vertex_x'].values / um_per_pixel
            vy = cell_bounds['vertex_y'].values / um_per_pixel

            # Bounding box
            min_x, max_x = int(np.floor(vx.min())), int(np.ceil(vx.max()))
            min_y, max_y = int(np.floor(vy.min())), int(np.ceil(vy.max()))

            # Clip to map bounds (numpy: array[row=y, col=x])
            min_x = max(0, min(min_x, incoherence_map.shape[1] - 1))
            max_x = max(0, min(max_x, incoherence_map.shape[1]))
            min_y = max(0, min(min_y, incoherence_map.shape[0] - 1))
            max_y = max(0, min(max_y, incoherence_map.shape[0]))

            if min_x >= max_x or min_y >= max_y:
                continue

            region = incoherence_map[min_y:max_y, min_x:max_x]
            if region.size == 0:
                continue

            results.append({
                'cell_id': cid,
                'incoherence_mean': float(np.nanmean(region)),
                'incoherence_max': float(np.nanmax(region)),
                'incoherence_median': float(np.nanmedian(region)),
            })

        if not results:
            return

        inc_df = pd.DataFrame(results).set_index('cell_id')
        adata.uns['cell_incoherence'] = inc_df

        inc_path = os.path.join(output_dir, f"{sample_tag}_step2_cell_incoherence.csv")
        inc_df.to_csv(inc_path)

        # Histogram
        vals = inc_df['incoherence_mean'].dropna()
        if len(vals) > 0:
            plt.figure(figsize=(8, 5))
            plt.hist(vals, bins=50, range=(0, 1), edgecolor='k', alpha=0.7)
            plt.axvline(0.2, color='grey', linestyle='--', label='Low coherence threshold')
            low_pct = (vals < 0.2).mean() * 100
            plt.title(
                f'Per-Cell Incoherence ({sample_tag})\n'
                f'{low_pct:.1f}% cells with low coherence (<0.2)'
            )
            plt.xlabel('Mean incoherence (1 - VSI)')
            plt.ylabel('Count')
            plt.legend()
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_incoherence_hist.png"),
                dpi=150
            )
            plt.close()
            logger.info(f"Saved incoherence histogram ({low_pct:.1f}% low coherence)")

    except Exception as e:
        logger.warning(f"Per-cell incoherence failed: {e}")


def _plot_ovrlpy_roi_vertical_cuts(ovrlp_obj, coords, um_per_pixel,
                                    output_dir, sample_tag):
    """ROI vertical cross-sections showing z-axis signal coherence.

    Paper Fig 1e (top): x-cut and y-cut through a 500µm ROI showing
    transcript z-positions coloured by coherence/integrity.
    Also plots ROI celltype map (top/bottom z-halves).
    """
    integrity = getattr(ovrlp_obj, 'integrity_map', None)
    signal = getattr(ovrlp_obj, 'signal_map', None)
    if integrity is None or signal is None:
        logger.info("No integrity/signal maps for ROI vertical cuts.")
        return

    # Pick ROI: center of tissue, 250µm radius (500µm box)
    x_all = coords['x'].values * um_per_pixel  # back to µm
    y_all = coords['y'].values * um_per_pixel
    x_med, y_med = np.nanmedian(x_all), np.nanmedian(y_all)
    roi_half = 250
    x_min, x_max = x_med - roi_half, x_med + roi_half
    y_min, y_max = y_med - roi_half, y_med + roi_half

    roi_mask = ((x_all >= x_min) & (x_all <= x_max) &
                (y_all >= y_min) & (y_all <= y_max))
    roi = coords[roi_mask].copy()
    roi['x_um'] = roi['x'].values * um_per_pixel
    roi['y_um'] = roi['y'].values * um_per_pixel

    if len(roi) < 100:
        logger.warning("Too few transcripts in ROI for vertical cuts.")
        return

    z_vals = roi['z'].values
    z_med = np.median(z_vals)

    # Lookup integrity for each transcript (pixel coords)
    roi_px_x = np.clip(roi['x'].values.astype(int), 0, integrity.shape[1] - 1)
    roi_px_y = np.clip(roi['y'].values.astype(int), 0, integrity.shape[0] - 1)
    roi['integrity'] = integrity[roi_px_y, roi_px_x]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel (0,0): x-cut — y vs z, coloured by integrity
    y_band_center = y_med
    y_band_half = 25  # 50µm strip
    x_strip = roi[(roi['y_um'] >= y_band_center - y_band_half) &
                   (roi['y_um'] <= y_band_center + y_band_half)]
    if len(x_strip) > 0:
        sc0 = axes[0, 0].scatter(x_strip['x_um'], x_strip['z'],
                                  s=0.5, c=x_strip['integrity'],
                                  cmap='RdYlGn', vmin=0, vmax=1,
                                  alpha=0.5, rasterized=True)
        axes[0, 0].set_xlabel('x (µm)')
        axes[0, 0].set_ylabel('z (µm)')
        axes[0, 0].set_title(f'X-cut (y={y_band_center:.0f}±{y_band_half}µm)')
        plt.colorbar(sc0, ax=axes[0, 0], label='Integrity', shrink=0.8)
    else:
        axes[0, 0].set_title('X-cut (no transcripts in strip)')
        axes[0, 0].text(0.5, 0.5, 'No data', ha='center', va='center',
                         transform=axes[0, 0].transAxes, fontsize=12, color='gray')

    # Panel (0,1): y-cut — x vs z, coloured by integrity
    x_band_center = x_med
    x_band_half = 25
    y_strip = roi[(roi['x_um'] >= x_band_center - x_band_half) &
                   (roi['x_um'] <= x_band_center + x_band_half)]
    if len(y_strip) > 0:
        sc1 = axes[0, 1].scatter(y_strip['y_um'], y_strip['z'],
                                  s=0.5, c=y_strip['integrity'],
                                  cmap='RdYlGn', vmin=0, vmax=1,
                                  alpha=0.5, rasterized=True)
        axes[0, 1].set_xlabel('y (µm)')
        axes[0, 1].set_ylabel('z (µm)')
        axes[0, 1].set_title(f'Y-cut (x={x_band_center:.0f}±{x_band_half}µm)')
        plt.colorbar(sc1, ax=axes[0, 1], label='Integrity', shrink=0.8)
    else:
        axes[0, 1].set_title('Y-cut (no transcripts in strip)')
        axes[0, 1].text(0.5, 0.5, 'No data', ha='center', va='center',
                         transform=axes[0, 1].transAxes, fontsize=12, color='gray')

    # Panel (1,0): ROI top half (z > median) — spatial coloured by integrity
    top_half = roi[roi['z'] >= z_med]
    if len(top_half) > 0:
        sc2 = axes[1, 0].scatter(top_half['x_um'], top_half['y_um'],
                                  s=0.3, c=top_half['integrity'],
                                  cmap='RdYlGn', vmin=0, vmax=1,
                                  alpha=0.4, rasterized=True)
        axes[1, 0].set_title(f'Top z-half (z≥{z_med:.1f}µm, n={len(top_half):,})')
        axes[1, 0].invert_yaxis()
        axes[1, 0].set_aspect('equal')
        plt.colorbar(sc2, ax=axes[1, 0], label='Integrity', shrink=0.8)

    # Panel (1,1): ROI bottom half (z < median)
    bot_half = roi[roi['z'] < z_med]
    if len(bot_half) > 0:
        sc3 = axes[1, 1].scatter(bot_half['x_um'], bot_half['y_um'],
                                  s=0.3, c=bot_half['integrity'],
                                  cmap='RdYlGn', vmin=0, vmax=1,
                                  alpha=0.4, rasterized=True)
        axes[1, 1].set_title(f'Bottom z-half (z<{z_med:.1f}µm, n={len(bot_half):,})')
        axes[1, 1].invert_yaxis()
        axes[1, 1].set_aspect('equal')
        plt.colorbar(sc3, ax=axes[1, 1], label='Integrity', shrink=0.8)

    fig.suptitle(f'ROI Vertical Cuts — Z-axis Signal Coherence ({sample_tag})\n'
                 f'(Paper Fig 1e style, 500µm × 500µm ROI)',
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    path = os.path.join(output_dir, f"{sample_tag}_step2_ovrlpy_roi_vertical_cuts.png")
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved ovrlpy ROI vertical cuts: {path}")


# ============================================================
# 2-3a. Distance to centroid
# Ref: xb.calculating.dispersion(), notebook 2_1
# ============================================================

def calculate_distance_to_centroid(adata, input_dir):
    """
    Calculates Euclidean distance from each transcript to its assigned
    cell's centroid.  Implements the approach from xb.calculating.dispersion().
    """
    logger.info("Starting Distance to Centroid Analysis...")

    spots = adata.uns['spots'].copy()
    cells_metadata = adata.obs

    if 'cell_id' not in spots.columns:
        logger.error("'cell_id' not found in spots. Cannot compute centroid distance.")
        return adata

    if 'x_centroid' not in cells_metadata.columns or 'y_centroid' not in cells_metadata.columns:
        logger.error("Centroid columns not found in adata.obs. Cannot compute centroid distance.")
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

    logger.info(f"Distance to centroid complete. Stats for {len(gene_stats)} genes.")
    return adata


# ============================================================
# 2-3b. Distance to boundary (EDT) -- SIGNED
# Ref: points2regions/compute_distance.ipynb
# Convention: negative = inside nucleus, positive = outside
# ============================================================

def calculate_distance_to_boundary(adata, input_dir, pixel_size=0.5, padding=50):
    """
    Calculates SIGNED distance from each transcript to the nearest
    nucleus boundary using rasterization + EDT.

    Sign convention (matching notebook compute_distance.ipynb):
        overlaps_nucleus == 1  -->  negative distance (inside nucleus)
        overlaps_nucleus == 0  -->  positive distance (outside nucleus)
    """
    logger.info("Starting Signed Distance to Nucleus Boundary (EDT)...")

    boundaries = load_boundaries(input_dir)
    if boundaries is None:
        return adata

    spots = adata.uns['spots']
    min_x = min(spots['x_location'].min(), boundaries['vertex_x'].min()) - padding
    min_y = min(spots['y_location'].min(), boundaries['vertex_y'].min()) - padding
    max_x = max(spots['x_location'].max(), boundaries['vertex_x'].max()) + padding
    max_y = max(spots['y_location'].max(), boundaries['vertex_y'].max()) + padding

    width = int((max_x - min_x) / pixel_size) + 1
    height = int((max_y - min_y) / pixel_size) + 1

    logger.info(f"Rasterizing onto {height}x{width} grid (pixel_size={pixel_size} um)...")

    mask = np.zeros((height, width), dtype=bool)
    warnings.filterwarnings("ignore", category=FutureWarning)

    for cell_id, grp in boundaries.groupby('cell_id'):
        vx = (grp['vertex_x'].values - min_x) / pixel_size
        vy = (grp['vertex_y'].values - min_y) / pixel_size
        rr, cc = polygon(vy, vx, shape=mask.shape)
        mask[rr, cc] = True

    # Signed EDT
    logger.info("Computing signed EDT...")
    dist_outside = distance_transform_edt(~mask) * pixel_size   # positive outside
    dist_inside = distance_transform_edt(mask) * pixel_size      # positive inside

    signed_map = np.where(mask, -dist_inside, dist_outside)

    # Map transcripts
    logger.info("Mapping transcripts to signed distance map...")
    px = (spots['x_location'].values - min_x) / pixel_size
    py = (spots['y_location'].values - min_y) / pixel_size
    valid = (px >= 0) & (px < width) & (py >= 0) & (py < height)

    dists = np.full(len(spots), np.nan)
    dists[valid] = map_coordinates(signed_map, [py[valid], px[valid]], order=1)

    adata.uns['spots']['dist_to_nucleus'] = dists

    logger.info("Signed distance to boundary complete.")
    return adata


# ============================================================
# 2-4. Centroid distance plots
# Ref: notebook 2_1
# ============================================================

def plot_centroid_distance_analysis(adata, output_dir, sample_tag, n_top_genes=20):
    """
    Per-gene centroid distance visualizations:
      1. Stripplot + boxplot for extreme genes
      2. Color-coded boxplot
      3. ECDF
      4. Welch's t-test p-value heatmap
    """
    logger.info("Generating per-gene centroid distance visualizations...")

    spots = adata.uns['spots']
    if 'dist_to_centroid' not in spots.columns:
        logger.warning("dist_to_centroid not found. Skipping centroid distance plots.")
        return

    spots_valid = spots.dropna(subset=['dist_to_centroid']).copy()
    gene_mask = ~spots_valid['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    spots_valid = spots_valid[gene_mask]

    if spots_valid.empty:
        logger.warning("No valid transcripts with centroid distance.")
        return

    # 10% subsample (as in notebook)
    np.random.seed(0)
    n_sub = max(1, len(spots_valid) // 10)
    rp = spots_valid.iloc[np.random.permutation(len(spots_valid))[:n_sub]].copy()

    # Extreme genes: top-5 closest + top-5 farthest (among top-100 expressed)
    dist_by_gene = rp.groupby('feature_name')['dist_to_centroid'].mean().sort_values()
    gene_counts = rp.groupby('feature_name')['dist_to_centroid'].count().sort_values(ascending=False)
    high_exp = gene_counts.index[:100]
    dist_by_gene = dist_by_gene[dist_by_gene.index.isin(high_exp)]

    n_extreme = min(5, len(dist_by_gene) // 2)
    if n_extreme < 1:
        logger.warning("Not enough genes for extreme gene analysis.")
        return

    extreme_genes = list(dist_by_gene.index[:n_extreme]) + list(dist_by_gene.index[-n_extreme:])

    # Mean cell border distance (from overlaps_nucleus)
    mean_cellborder = None
    if 'overlaps_nucleus' in rp.columns:
        rp_nuc = rp[rp['overlaps_nucleus'] == 1]
        if not rp_nuc.empty:
            dist_max = rp_nuc.groupby('cell_id')['dist_to_centroid'].max()
            mean_cellborder = np.mean(dist_max)

    # --- 2-4a. Stripplot + boxplot ---
    d_ext = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_ext.empty:
        min_grp = d_ext['feature_name'].value_counts().min()
        d_ext = d_ext.groupby('feature_name').head(min_grp).reset_index(drop=True)

        plt.figure(figsize=(12, 6))
        sns.stripplot(
            x=d_ext['dist_to_centroid'], y=d_ext['feature_name'],
            s=0.15, order=extreme_genes, jitter=0.4
        )
        sns.boxplot(
            x=d_ext['dist_to_centroid'], y=d_ext['feature_name'],
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
        logger.info("Saved centroid stripplot.")

        # --- 2-4b. Color-coded boxplot ---
        color_dict = {}
        for g in dist_by_gene.index[:n_extreme]:
            color_dict[g] = 'red'
        for g in dist_by_gene.index[-n_extreme:]:
            color_dict[g] = 'blue'

        plt.figure(figsize=(10, 5))
        sns.boxplot(
            x=d_ext['dist_to_centroid'], y=d_ext['feature_name'],
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
        logger.info("Saved centroid boxplot.")

    # --- 2-4c. ECDF ---
    d_ecdf = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_ecdf.empty:
        sns.displot(x=d_ecdf['dist_to_centroid'], kind='ecdf', hue=d_ecdf['feature_name'])
        if mean_cellborder is not None:
            plt.axvline(x=mean_cellborder, color='grey', linestyle='--', label='Mean cell border')
        plt.title(f'ECDF of Distance to Centroid ({sample_tag})')
        plt.xlabel('Distance to Centroid (um)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_ecdf.png"), dpi=150)
        plt.close()
        logger.info("Saved centroid ECDF.")

    # --- 2-4d. Welch's t-test heatmap ---
    d_tt = rp[rp['feature_name'].isin(extreme_genes)].copy()
    if not d_tt.empty and len(extreme_genes) > 1:
        pval_matrix = pd.DataFrame(columns=extreme_genes, index=extreme_genes, dtype=float)

        for g1 in extreme_genes:
            for g2 in extreme_genes:
                v1 = d_tt.loc[d_tt['feature_name'] == g1, 'dist_to_centroid']
                v2 = d_tt.loc[d_tt['feature_name'] == g2, 'dist_to_centroid']
                if len(v1) > 1 and len(v2) > 1:
                    pval_matrix.loc[g1, g2] = stats.ttest_ind(v1, v2, equal_var=False).pvalue
                else:
                    pval_matrix.loc[g1, g2] = np.nan

        pval_path = os.path.join(output_dir, f"{sample_tag}_step2_centroid_welchs_pvals.csv")
        pval_matrix.to_csv(pval_path)

        plt.figure(figsize=(10, 8))
        pval_float = pval_matrix.astype(float)
        with np.errstate(divide='ignore'):
            log_pvals = -np.log10(pval_float.values.astype(float))
        log_df = pd.DataFrame(log_pvals, index=extreme_genes, columns=extreme_genes)
        sns.heatmap(log_df, annot=True, fmt='.1f', cmap='YlOrRd',
                    square=True, cbar_kws={'label': '-log10(p-value)'})
        plt.title(f"Welch's t-test: -log10(p-value) ({sample_tag})")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_tag}_step2_centroid_welchs_heatmap.png"), dpi=150)
        plt.close()
        logger.info("Saved Welch's t-test heatmap.")

    # Save gene stats + plotting subset
    if 'gene_distance_stats' in adata.uns:
        adata.uns['gene_distance_stats'].to_csv(
            os.path.join(output_dir, f"{sample_tag}_step2_gene_distance_stats.csv")
        )
    plotting_subset = spots[
        spots['feature_name'].isin(extreme_genes)
    ][['x_location', 'y_location', 'feature_name', 'dist_to_centroid']].copy()
    plotting_subset.to_csv(
        os.path.join(output_dir, f"{sample_tag}_step2_centroid_plotting_subset.csv"),
        index=False
    )


# ============================================================
# 2-5. P2R Boundary Distance Box Plots
# Ref: points2regions/figures.ipynb, Extended Data Fig. 3c
# ============================================================

def plot_p2r_boundary_distance(adata, output_dir, sample_tag):
    """
    Box plots of signed distance to nuclear edge per P2R cluster,
    grouped by cell type.  Negative = inside nucleus.
    """
    logger.info("Generating P2R boundary distance box plots...")

    spots = adata.uns['spots']
    if 'dist_to_nucleus' not in spots.columns:
        logger.warning("dist_to_nucleus not found. Skipping P2R boundary plots.")
        return
    if 'p2r_name' not in spots.columns:
        logger.warning("P2R classification not found. Skipping boundary plots.")
        return

    classification_df = adata.uns.get('p2r_classification', pd.DataFrame())
    if classification_df.empty:
        return

    color_palette = adata.uns.get('p2r_colors', {})

    ct_counts = classification_df.groupby('majority_celltype').size()
    interesting_cts = ct_counts[ct_counts >= 3].index.tolist()
    if not interesting_cts:
        interesting_cts = ct_counts.nlargest(2).index.tolist()

    for selected_ct in interesting_cts[:3]:
        ct_clusters = classification_df[classification_df['majority_celltype'] == selected_ct]
        ct_names = ct_clusters.sort_values('nuc_fraction', ascending=False)['p2r_name'].tolist()

        subset = spots[spots['p2r_name'].isin(ct_names)].dropna(subset=['dist_to_nucleus'])
        if subset.empty:
            continue

        median_dists = subset.groupby('p2r_name')['dist_to_nucleus'].median()
        order = median_dists.sort_values().index.tolist()

        plt.figure(figsize=(8, max(4, len(order) * 0.6)))
        palette = {n: color_palette.get(n, '#888888') for n in order}

        # Clip extreme distances for readable box plots (99th percentile)
        dist_q01 = subset['dist_to_nucleus'].quantile(0.01)
        dist_q99 = subset['dist_to_nucleus'].quantile(0.99)
        subset_clipped = subset[
            (subset['dist_to_nucleus'] >= dist_q01) & (subset['dist_to_nucleus'] <= dist_q99)
        ]

        sns.boxplot(
            data=subset_clipped, y='p2r_name', x='dist_to_nucleus',
            order=order, showfliers=False, palette=palette,
            width=0.6, linewidth=0.9
        )
        plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Nuclear edge')
        plt.xlim(dist_q01 * 1.1, dist_q99 * 1.1)
        plt.xlabel('Distance to nuclear edge (um)')
        plt.ylabel('P2R cluster')
        plt.title(f'Distance to Nuclear Edge - {selected_ct} ({sample_tag})')
        plt.legend(loc='lower right')
        plt.tight_layout()

        safe_ct = selected_ct.replace(' ', '_').replace('/', '-')
        plt.savefig(
            os.path.join(output_dir, f"{sample_tag}_step2_p2r_boundary_{safe_ct}.png"),
            dpi=150, bbox_inches='tight'
        )
        plt.close()
        logger.info(f"Saved P2R boundary distance plot for {selected_ct}")


# ============================================================
# 2-5b. P2R Mean Distance Scatter Plot
# Ref: points2regions/compute_colors.ipynb
# ============================================================

def plot_p2r_mean_distance_scatter(adata, output_dir, sample_tag):
    """
    Scatter plot of mean signed distance to nucleus per P2R cluster,
    colored by HLS palette.  X-axis: cluster index, Y-axis: mean distance.
    Ref: compute_colors.ipynb visualization.
    """
    logger.info("Generating P2R mean distance scatter plot...")

    classification_df = adata.uns.get('p2r_classification', pd.DataFrame())
    if classification_df.empty:
        logger.warning("P2R classification not found. Skipping mean distance scatter.")
        return

    spots = adata.uns['spots']
    if 'dist_to_nucleus' not in spots.columns:
        logger.warning("dist_to_nucleus not found. Skipping mean distance scatter.")
        return
    if 'points2regions' not in spots.columns:
        return

    # Compute mean distance per cluster
    mean_dist = spots.dropna(subset=['dist_to_nucleus']).groupby(
        'points2regions'
    )['dist_to_nucleus'].mean()

    classification_df = classification_df.copy()
    classification_df['mean_distance'] = classification_df['cluster_id'].map(mean_dist)
    classification_df = classification_df.dropna(subset=['mean_distance'])

    if classification_df.empty:
        return

    # Sort by mean distance
    classification_df = classification_df.sort_values('mean_distance')

    color_palette = adata.uns.get('p2r_colors', {})
    colors = [color_palette.get(n, '#888888') for n in classification_df['p2r_name']]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(
        range(len(classification_df)),
        classification_df['mean_distance'].values,
        c=colors, s=40, edgecolors='k', linewidth=0.3, zorder=3
    )
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Nuclear edge')

    # Label extreme clusters
    n_label = min(3, len(classification_df))
    for idx in list(range(n_label)) + list(range(-n_label, 0)):
        row = classification_df.iloc[idx]
        ax.annotate(
            row['p2r_name'], (idx if idx >= 0 else len(classification_df) + idx, row['mean_distance']),
            fontsize=5, rotation=30, ha='left', va='bottom'
        )

    ax.set_xlabel('Cluster (sorted by mean distance)')
    ax.set_ylabel('Mean signed distance to nucleus (um)')
    ax.set_title(f'P2R Cluster Mean Distance to Nuclear Edge ({sample_tag})')
    ax.legend(loc='upper left')
    fig.tight_layout()

    path = os.path.join(output_dir, f"{sample_tag}_step2_p2r_mean_distance_scatter.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved P2R mean distance scatter to {path}")


# ============================================================
# 2-7. SSAM Analysis (optional)
# Ref: notebooks/2_4_brain_ssam.ipynb, Extended Data Fig. 2
# ============================================================

def _patch_zarr_v3_compat():
    """Patch zarr v3 Group.array to accept zarr v2-style calls (data= without shape/dtype).
    SSAM v1.1.3 uses zarr v2 API: zarr_group.array(name=..., data=...).
    zarr v3 enforces shape/dtype via typing_extensions wrapper on Group.array.
    Redirect to Group.create_array which accepts data-only calls.
    Must be called BEFORE importing ssam.
    """
    import zarr
    if int(zarr.__version__.split('.')[0]) < 3:
        return  # No patch needed for zarr v2

    if hasattr(zarr.Group, '_ssam_patched'):
        return  # Already patched

    def _compat_array(self, name, shape=None, dtype=None, *, data=None, **kwargs):
        if data is not None:
            data = np.asarray(data)
            return self.create_array(name=name, data=data, overwrite=True)
        return self.create_array(name=name, shape=shape, dtype=dtype, overwrite=True, **kwargs)

    zarr.Group.array = _compat_array

    # Patch parse_shapelike to accept numpy int64 values in shape tuples.
    # SSAM passes (np.int64(W), np.int64(H), np.int64(D), N) which zarr v3 rejects.
    try:
        from zarr.core.common import parse_shapelike as _orig_parse_shapelike

        def _patched_parse_shapelike(data):
            if isinstance(data, (list, tuple)):
                data = tuple(int(x) for x in data)
            return _orig_parse_shapelike(data)

        # Patch ALL zarr modules that import parse_shapelike
        for _mod_name in [
            'zarr.core.common', 'zarr.core.array', 'zarr.core.array_spec',
            'zarr.core.chunk_grids', 'zarr.core.group',
            'zarr.core.metadata.v2', 'zarr.core.metadata.v3',
            'zarr.codecs.sharding',
        ]:
            try:
                _mod = __import__(_mod_name, fromlist=['parse_shapelike'])
                if hasattr(_mod, 'parse_shapelike'):
                    setattr(_mod, 'parse_shapelike', _patched_parse_shapelike)
            except (ImportError, AttributeError):
                pass
    except ImportError:
        pass  # zarr v3 structure changed, skip

    zarr.Group._ssam_patched = True


def run_ssam_analysis(adata, config, output_dir, sample_tag):
    """
    SSAM v1.1+ segmentation-free cell type mapping
    (Ref: notebooks/2_4_brain_ssam.ipynb, Extended Data Fig. 2):
      1. Build KDE vector field from transcript coordinates (run_kde)
      2. Set expression/norm thresholds, detect local maxima
      3. Normalize vectors, cluster via Leiden
      4. Map clusters to cell types via reference signatures
      5. Generate spatial map and UMAP
    """
    logger.info("Starting SSAM analysis...")

    # Patch zarr v3 before importing ssam
    _patch_zarr_v3_compat()

    try:
        import ssam
    except ImportError:
        logger.warning("ssam not installed (pip install ssam). Skipping SSAM analysis.")
        return adata

    ssam_cfg = config['segmentation_free'].get('ssam', {})
    um_per_px = ssam_cfg.get('ssam_vf_um_per_px', 2.0)
    kde_bw = ssam_cfg.get('kde_bandwidth_um', 2.5)
    norm_thres = ssam_cfg.get('norm_thres', 5)
    exp_thres = ssam_cfg.get('exp_thres', 0.2)
    min_dist = ssam_cfg.get('min_dist', 7)

    spots = adata.uns['spots']
    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    coords = spots[gene_mask][['x_location', 'y_location', 'feature_name']].copy()
    coords.columns = ['x', 'y', 'gene']

    # Shift to origin and convert to pixel coordinates
    coords[['x', 'y']] = coords[['x', 'y']] - coords[['x', 'y']].min()
    coords['x'] = coords['x'] / um_per_px
    coords['y'] = coords['y'] / um_per_px
    coords['gene'] = coords['gene'].astype(str)

    genes = sorted(coords['gene'].unique().tolist())
    logger.info(f"Preparing SSAM: {len(genes)} genes, {len(coords)} transcripts")

    # Image dimensions in pixels
    width_px = int(np.ceil(coords['x'].max())) + 1
    height_px = int(np.ceil(coords['y'].max())) + 1

    try:
        # SSAM 1.0.2 API: SSAMDataset needs per-gene location arrays
        # Build per-gene locations list: list of N_i × 2 arrays
        locations_list = []
        for gene in genes:
            gene_coords = coords[coords['gene'] == gene][['x', 'y']].values
            locations_list.append(gene_coords)

        ds = ssam.SSAMDataset(genes, locations_list, width_px, height_px, depth=1)

        # Persistent save_dir enables SSAM's built-in per-gene KDE npy/pkl caching
        ssam_kde_cache_dir = os.path.join(output_dir, 'ssam_kde_cache')
        os.makedirs(ssam_kde_cache_dir, exist_ok=True)
        analysis = ssam.SSAMAnalysis(ds, ncores=4, verbose=True,
                                     save_dir=ssam_kde_cache_dir)

        # KDE (bandwidth in pixel units)
        # Check if ssam_adata already exists (from a previous run) → skip KDE entirely
        kde_bw_px = kde_bw / um_per_px
        ssam_h5ad_path = os.path.join(output_dir, f"{sample_tag}_step2_ssam.h5ad")
        if os.path.exists(ssam_h5ad_path):
            logger.info(f"SSAM h5ad cache found ({ssam_h5ad_path}), loading cached results...")
            ssam_adata = sc.read_h5ad(ssam_h5ad_path)
            logger.info(f"  Loaded cached: {ssam_adata.n_obs} pseudo-cells, {ssam_adata.n_vars} genes")

            # Still need to run KDE + downstream for ds.plot_* calls
            # SSAM caches per-gene KDE in save_dir as npy files — subsequent runs load from cache
            logger.info(f"Running SSAM KDE for visualization (bandwidth={kde_bw_px:.2f} px)...")
            analysis.run_kde(bandwidth=kde_bw_px, sampling_distance=1.0)

            analysis.set_thresholds(
                expression_threshold=exp_thres,
                norm_threshold=norm_thres,
            )
            search_size = max(3, int(np.round(min_dist / um_per_px)))
            if search_size % 2 == 0:
                search_size += 1
            analysis.find_localmax(search_size=search_size)
            analysis.normalize_vectors(normalize_vector=True)
            analysis.scale_vectors()

            try:
                analysis.cluster_vectors()
                ds.run_umap()
            except Exception as e:
                logger.warning(f"SSAM cluster/umap failed: {e}")

            # Rebuild pixel-level map from cached leiden signatures (needed for ds.plot_celltypes_map)
            mapping_result = {'has_reference': 'leiden_assignment' in ssam_adata.obs.columns}
            if mapping_result['has_reference']:
                unique_types = sorted(ssam_adata.obs['leiden_assignment'].cat.categories)
                n_colors = max(len(unique_types), 1)
                palette = {ct: plt.cm.tab20(i / n_colors) for i, ct in enumerate(unique_types)}
                mapping_result['palette'] = palette
                # assignment_names must be per-leiden-cluster (one per centroid),
                # not per-unique-celltype, so ds.plot_celltypes_map colors align.
                if 'ssam_celltype_assignment_names' in ssam_adata.uns:
                    mapping_result['assignment_names'] = list(ssam_adata.uns['ssam_celltype_assignment_names'])
                else:
                    # Reconstruct from leiden → leiden_assignment mapping
                    leiden_cats = ssam_adata.obs['leiden'].cat.categories
                    cluster_map = ssam_adata.obs.groupby('leiden')['leiden_assignment'].first()
                    mapping_result['assignment_names'] = [cluster_map[cl] for cl in leiden_cats]
                if 'ssam_leiden_celltype_correlations' in ssam_adata.uns:
                    mapping_result['corr_matrix'] = ssam_adata.uns['ssam_leiden_celltype_correlations']

                # Re-run pixel-level map_celltypes so ds.plot_celltypes_map works
                try:
                    common_genes = [g for g in genes if g in ssam_adata.var_names]
                    ssam_sub = ssam_adata[:, common_genes]
                    leiden_cats = ssam_adata.obs['leiden'].cat.categories
                    lei_sig = np.zeros((len(common_genes), len(leiden_cats)))
                    for j, cl in enumerate(leiden_cats):
                        mask = ssam_adata.obs['leiden'] == cl
                        X = ssam_sub[mask].X
                        X = X.toarray() if hasattr(X, 'toarray') else X
                        lei_sig[:, j] = X.mean(axis=0)
                    analysis.map_celltypes(centroids=lei_sig.T)
                    analysis.filter_celltypemaps(min_norm=3, min_r=0.2)
                    logger.info("Re-built pixel-level cell type map from cached signatures.")
                except Exception as e:
                    logger.warning(f"Pixel-level map rebuild failed (cache): {e}")

            _ssam_plots(ssam_adata, ds, mapping_result, output_dir, sample_tag)

            ssam_summary = {
                'n_localmax': ssam_adata.n_obs,
                'n_leiden_clusters': len(ssam_adata.obs['leiden'].cat.categories),
                'celltype_mapping': 'cached',
            }
            adata.uns['ssam_summary'] = ssam_summary
            logger.info("SSAM plots regenerated from cache.")
            return adata

        logger.info(f"Running SSAM KDE (bandwidth={kde_bw_px:.2f} px, grid {width_px}x{height_px})...")
        analysis.run_kde(bandwidth=kde_bw_px, sampling_distance=1.0)

        # Set thresholds before finding local maxima
        analysis.set_thresholds(
            expression_threshold=exp_thres,
            norm_threshold=norm_thres,
        )

        # Local maxima (search_size in pixels, must be odd integer)
        search_size = max(3, int(np.round(min_dist / um_per_px)))
        if search_size % 2 == 0:
            search_size += 1
        logger.info(f"Finding SSAM local maxima (search_size={search_size})...")
        analysis.find_localmax(search_size=search_size)

        # Normalize and scale vectors (both steps required per SSAM tutorial)
        analysis.normalize_vectors(normalize_vector=True)
        analysis.scale_vectors()

        # SSAM internal clustering + UMAP (needed for ds.plot_diagnostic_plot)
        try:
            analysis.cluster_vectors()
            ds.run_umap()
            logger.info("SSAM internal clustering + UMAP complete.")
        except Exception as e:
            logger.warning(f"SSAM internal cluster_vectors/run_umap failed: {e}")

        # Build AnnData from SSAM vectors at local maxima
        normalized = ds.normalized_vectors
        local_maxs = ds.local_maxs
        if normalized is None or local_maxs is None or len(local_maxs[0]) == 0:
            logger.warning("SSAM found no local maxima. Skipping downstream analysis.")
            return adata

        # Get gene names from zarr store (zarr v3 returns numpy arrays)
        if 'genes' in ds.zarr_group:
            ssam_genes = [str(g) for g in ds.zarr_group['genes'][:]]
        else:
            ssam_genes = genes

        ssam_adata = sc.AnnData(
            np.array(normalized),
            var=pd.DataFrame(index=ssam_genes),
            obs=pd.DataFrame({'x': local_maxs[0], 'y': local_maxs[1]})
        )

        logger.info(f"SSAM: {ssam_adata.n_obs} local maxima (pseudo-cells)")

        # Scanpy pipeline for clustering (tutorial: normalize → log1p → PCA → ...)
        sc.pp.normalize_total(ssam_adata, target_sum=1e4)
        sc.pp.log1p(ssam_adata)
        sc.tl.pca(ssam_adata, svd_solver='arpack')
        n_pcs = min(40, ssam_adata.n_vars - 1, ssam_adata.n_obs - 1)
        sc.pp.neighbors(ssam_adata, n_neighbors=15, n_pcs=n_pcs)
        sc.tl.umap(ssam_adata, min_dist=0.02, random_state=42)
        sc.tl.leiden(ssam_adata, resolution=2.0, random_state=42)

        n_clusters = len(ssam_adata.obs['leiden'].cat.categories)
        logger.info(f"SSAM Leiden: {n_clusters} clusters")

        # Cell type mapping via reference signatures
        try:
            mapping_result = _ssam_map_celltypes(analysis, ssam_adata, config, ssam_genes)
        except Exception as e:
            logger.warning(f"SSAM cell type mapping failed: {e}. Proceeding with de novo clusters only.")
            mapping_result = {'has_reference': False}

        # Save
        ssam_path = os.path.join(output_dir, f"{sample_tag}_step2_ssam.h5ad")
        ssam_adata.write_h5ad(ssam_path)
        logger.info(f"Saved SSAM AnnData to {ssam_path}")

        # Plots
        _ssam_plots(ssam_adata, ds, mapping_result, output_dir, sample_tag)

        ssam_summary = {
            'n_localmax': ssam_adata.n_obs,
            'n_leiden_clusters': n_clusters,
        }
        if mapping_result.get('has_reference', False):
            ssam_summary['celltype_mapping'] = 'reference-correlated'
            ssam_summary['n_common_genes'] = mapping_result['n_common_genes']
            ssam_summary['cluster_assignments'] = mapping_result['cluster_assignments']
        else:
            ssam_summary['celltype_mapping'] = 'de_novo_only'
        adata.uns['ssam_summary'] = ssam_summary

    except Exception as e:
        logger.error(f"SSAM analysis failed: {e}")
        import traceback
        traceback.print_exc()

    return adata


def _ssam_map_celltypes(analysis, ssam_adata, config, ssam_genes):
    """Map SSAM leiden clusters to cell types via scRNA-seq reference correlation.

    Workflow (matching notebook cells 30-35):
      1. Load scRNA-seq reference
      2. Build reference signatures per cell type (log1p normalized)
      3. Build leiden cluster mean expression signatures (already log1p via scanpy)
      4. Compute Pearson correlation matrix (reference celltypes × leiden clusters)
      5. Assign best-matching reference type per leiden cluster via argmax
      6. Write leiden_assignment to ssam_adata.obs
      7. Call analysis.map_celltypes with leiden sigs for pixel-level map
      8. Store correlation matrix in ssam_adata.uns

    Returns dict with mapping results, or {'has_reference': False} on failure.
    """
    import glob as glob_module
    from scipy.stats import pearsonr

    sc_cfg = config.get('sc_reference', {})
    ref_dir = sc_cfg.get('dest_dir', None)
    ct_key = sc_cfg.get('celltype_key', 'subclass_label')

    # --- Load reference ---
    ref_adata = None
    if ref_dir:
        h5ad_files = glob_module.glob(os.path.join(ref_dir, '*.h5ad'))
        if h5ad_files:
            try:
                ref_adata = sc.read_h5ad(h5ad_files[0])
            except Exception:
                pass

    if ref_adata is None or ct_key not in ref_adata.obs.columns:
        logger.info("No scRNA-seq reference for SSAM cell typing. Using de novo clusters only.")
        return {'has_reference': False}

    common_genes = [g for g in ssam_genes if g in ref_adata.var_names]
    if len(common_genes) < 10:
        logger.warning(f"Only {len(common_genes)} common genes with reference. Skipping mapping.")
        return {'has_reference': False}

    logger.info(f"Mapping SSAM clusters via {len(common_genes)} common genes...")

    # --- Build reference signatures (n_genes × n_celltypes), log1p ---
    ref_sub = ref_adata[:, common_genes]
    # Drop NaN cell types (unannotated cells in reference)
    ref_sub = ref_sub[ref_sub.obs[ct_key].notna()]
    celltypes = sorted(ref_sub.obs[ct_key].dropna().unique())
    ref_signatures = pd.DataFrame(index=common_genes, columns=celltypes, dtype=float)
    for ct in celltypes:
        ct_data = ref_sub[ref_sub.obs[ct_key] == ct]
        X = ct_data.X.toarray() if hasattr(ct_data.X, 'toarray') else ct_data.X
        ref_signatures[ct] = X.mean(axis=0)
    ref_signatures = np.log1p(ref_signatures)

    # --- Build leiden cluster signatures (n_genes × n_clusters), log1p ---
    ssam_sub = ssam_adata[:, [g for g in common_genes if g in ssam_adata.var_names]]
    leiden_cats = ssam_adata.obs['leiden'].cat.categories
    leiden_signatures = pd.DataFrame(
        index=[g for g in common_genes if g in ssam_adata.var_names],
        columns=leiden_cats, dtype=float
    )
    for cl in leiden_cats:
        mask = ssam_adata.obs['leiden'] == cl
        X = ssam_sub[mask].X
        X = X.toarray() if hasattr(X, 'toarray') else X
        leiden_signatures[cl] = X.mean(axis=0)
    # ssam_adata.X is already log1p-transformed (scanpy pipeline), so skip double log
    # leiden_signatures = np.log1p(leiden_signatures)  # removed: already log-space

    # Align gene sets
    shared_genes = leiden_signatures.index.intersection(ref_signatures.index)
    ref_sig = ref_signatures.loc[shared_genes]
    lei_sig = leiden_signatures.loc[shared_genes]

    # --- Pearson correlation matrix (celltypes × clusters) ---
    corr_matrix = pd.DataFrame(
        index=celltypes, columns=leiden_cats, dtype=float
    )
    for ct in celltypes:
        for cl in leiden_cats:
            r, _ = pearsonr(ref_sig[ct].values, lei_sig[cl].values)
            corr_matrix.loc[ct, cl] = r

    # --- Assign best-matching reference type per cluster ---
    cluster_assignments = {}
    for cl in leiden_cats:
        best_ct = corr_matrix[cl].astype(float).idxmax()
        cluster_assignments[cl] = best_ct

    assignment_names = [cluster_assignments[cl] for cl in leiden_cats]
    logger.info(f"Cluster → celltype assignments: {cluster_assignments}")

    # Write to ssam_adata
    ssam_adata.obs['leiden_assignment'] = (
        ssam_adata.obs['leiden'].map(cluster_assignments).astype('category')
    )
    ssam_adata.uns['ssam_leiden_celltype_correlations'] = corr_matrix
    ssam_adata.uns['ssam_celltype_assignment_names'] = assignment_names

    # --- Pixel-level map via SSAM using leiden cluster signatures ---
    try:
        # map_celltypes expects centroids: shape (n_clusters, n_genes)
        # Use leiden signatures so pixel map reflects our clusters
        analysis.map_celltypes(centroids=lei_sig.astype(float).values.T)
        analysis.filter_celltypemaps(min_norm=3, min_r=0.2)
        logger.info("SSAM pixel-level cell type mapping complete.")
    except Exception as e:
        logger.warning(f"SSAM pixel-level mapping failed: {e}")

    # --- Color palette ---
    unique_types = sorted(ssam_adata.obs['leiden_assignment'].cat.categories)
    n_colors = max(len(unique_types), 1)
    palette = {ct: plt.cm.tab20(i / n_colors) for i, ct in enumerate(unique_types)}

    return {
        'has_reference': True,
        'corr_matrix': corr_matrix,
        'cluster_assignments': cluster_assignments,
        'assignment_names': assignment_names,
        'n_common_genes': len(common_genes),
        'palette': palette,
    }


def _ssam_plots(ssam_adata, ds, mapping_result, output_dir, sample_tag):
    """Generate SSAM UMAP, spatial map, highest-expressed genes, and reference-mapped plots."""
    has_ref = mapping_result.get('has_reference', False)

    try:
        # UMAP colored by leiden (Step 1 style: on-data labels, small dots, clean look)
        sc.pl.umap(
            ssam_adata, color=['leiden'],
            size=8,
            legend_loc='on data',
            legend_fontsize=9,
            legend_fontoutline=2,
            title=f'SSAM Leiden Clusters ({sample_tag})',
            frameon=False, show=False,
        )
        plt.savefig(
            os.path.join(output_dir, f"{sample_tag}_step2_ssam_umap.png"),
            dpi=200, bbox_inches='tight'
        )
        plt.close()

        # Spatial scatter — color by leiden_assignment if available, else leiden
        color_col = 'leiden_assignment' if has_ref else 'leiden'
        fig, ax = plt.subplots(figsize=(14, 12))
        ax.set_facecolor('white')
        categories = ssam_adata.obs[color_col].cat.categories
        # Use a nicer palette: tab20 for many categories
        n_cats = len(categories)
        if n_cats <= 10:
            cmap = plt.cm.tab10
        elif n_cats <= 20:
            cmap = plt.cm.tab20
        else:
            cmap = plt.cm.gist_ncar
        cat_colors = {cat: cmap(i / max(n_cats - 1, 1)) for i, cat in enumerate(categories)}
        # Override with reference palette if available
        ref_palette = mapping_result.get('palette', None)
        if ref_palette:
            cat_colors.update({k: v for k, v in ref_palette.items() if k in cat_colors})

        for cat in categories:
            mask = ssam_adata.obs[color_col] == cat
            ax.scatter(
                ssam_adata.obs.loc[mask, 'x'], ssam_adata.obs.loc[mask, 'y'],
                c=[cat_colors[cat]], s=6, alpha=0.8, label=str(cat), rasterized=True,
                edgecolors='none',
            )
        ax.invert_yaxis()
        ax.set_title(f'SSAM Spatial Cell Type Map ({sample_tag})', fontsize=14, fontweight='bold')
        ax.set_xlabel('x (pixels)', fontsize=11)
        ax.set_ylabel('y (pixels)', fontsize=11)
        ax.set_aspect('equal')
        ax.set_facecolor('#f8f8f8')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        leg = ax.legend(markerscale=5, fontsize=9, loc='upper right', ncol=2,
                        frameon=True, fancybox=True, framealpha=0.9)
        leg.get_frame().set_edgecolor('lightgray')
        plt.tight_layout()
        plt.savefig(
            os.path.join(output_dir, f"{sample_tag}_step2_ssam_spatial.png"),
            dpi=200, bbox_inches='tight', facecolor='white'
        )
        plt.close()

        # Highest expressed genes
        sc.pl.highest_expr_genes(ssam_adata, n_top=20, show=False)
        plt.savefig(
            os.path.join(output_dir, f"{sample_tag}_step2_ssam_highest_genes.png"),
            dpi=150, bbox_inches='tight'
        )
        plt.close()

        # L1 Norm + Local Maxima overlay (SSAM tutorial: ds.plot_l1norm + ds.plot_localmax)
        try:
            fig, ax = plt.subplots(figsize=(14, 12))
            try:
                ds.plot_l1norm(ax=ax)
                ds.plot_localmax(s=0.3, c='red', ax=ax)
            except TypeError:
                plt.close()
                fig, ax = plt.subplots(figsize=(14, 12))
                ds.plot_l1norm()
                ds.plot_localmax(s=0.3, c='red')
                ax = plt.gca()
            ax.set_title(f'SSAM L1 Norm + Local Maxima ({sample_tag})')
            plt.tight_layout()
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_ssam_l1norm_localmax.png"),
                dpi=150, bbox_inches='tight'
            )
            plt.close()
            logger.info("SSAM L1 norm + local maxima plot saved.")
        except Exception as e:
            logger.warning(f"SSAM L1 norm plot failed: {e}")

        # SSAM-native UMAP (ds.run_umap result, different from scanpy UMAP)
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            try:
                ds.plot_umap(s=1, ax=ax)
            except TypeError:
                plt.close()
                fig, ax = plt.subplots(figsize=(10, 8))
                ds.plot_umap(s=1)
                ax = plt.gca()
            ax.set_title(f'SSAM Native UMAP ({sample_tag})')
            plt.tight_layout()
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_ssam_native_umap.png"),
                dpi=150, bbox_inches='tight'
            )
            plt.close()
            logger.info("SSAM native UMAP plot saved.")
        except Exception as e:
            logger.warning(f"SSAM native UMAP plot failed: {e}")

        # PCA variance plot
        try:
            sc.pl.pca(ssam_adata, show=False)
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_ssam_pca.png"),
                dpi=150, bbox_inches='tight'
            )
            plt.close()
            logger.info("SSAM PCA plot saved.")
        except Exception as e:
            logger.warning(f"SSAM PCA plot failed: {e}")

        # Per-cluster diagnostic plots (top 15 clusters max, SSAM tutorial pattern)
        try:
            n_diag = min(15, len(ssam_adata.obs['leiden'].cat.categories))
            for i in range(n_diag):
                try:
                    fig = plt.figure(figsize=(30, 5))
                    ds.plot_diagnostic_plot(i, use_embedding='umap')
                    plt.savefig(
                        os.path.join(output_dir, f"{sample_tag}_step2_ssam_diagnostic_{i}.png"),
                        dpi=100, bbox_inches='tight'
                    )
                    plt.close()
                except Exception:
                    plt.close()
            logger.info(f"SSAM diagnostic plots saved ({n_diag} clusters)")
        except Exception as e:
            logger.warning(f"SSAM diagnostic plots failed: {e}")

        # 3D KDE surface plot (PDF slide 61: wireframe of KDE density field)
        try:
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
            vf_norm = ds.vf_norm
            if vf_norm is not None:
                # Compute or access the L1/L2 norm; subsample for 3D rendering
                import dask.array as da
                norm_arr = vf_norm.compute() if isinstance(vf_norm, da.Array) else np.array(vf_norm)
                norm_2d = norm_arr.squeeze()
                # Subsample to ~200x200 for performance
                step_x = max(1, norm_2d.shape[0] // 200)
                step_y = max(1, norm_2d.shape[1] // 200)
                sub = norm_2d[::step_x, ::step_y]
                X, Y = np.meshgrid(np.arange(sub.shape[1]), np.arange(sub.shape[0]))
                fig = plt.figure(figsize=(14, 10))
                ax3d = fig.add_subplot(111, projection='3d')
                ax3d.plot_wireframe(X, Y, sub, rstride=1, cstride=1,
                                    linewidth=0.3, color='navy', alpha=0.7)
                ax3d.set_xlabel('X (subsampled)')
                ax3d.set_ylabel('Y (subsampled)')
                ax3d.set_zlabel('KDE Density (L2 norm)')
                ax3d.set_title(f'SSAM 3D KDE Surface ({sample_tag})')
                ax3d.view_init(elev=35, azim=225)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_kde_3d.png"),
                    dpi=150, bbox_inches='tight'
                )
                plt.close()
                logger.info("SSAM 3D KDE surface plot saved.")
        except Exception as e:
            logger.warning(f"SSAM 3D KDE surface plot failed: {e}")

        # Downsampling comparison: L1 norm heatmap vs L1 maxima scatter (PDF slide 63)
        try:
            import dask.array as da
            vf_norm = ds.vf_norm
            if vf_norm is not None:
                norm_arr = vf_norm.compute() if isinstance(vf_norm, da.Array) else np.array(vf_norm)
                norm_2d = norm_arr.squeeze()
                local_maxs = ds.local_maxs

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
                # Left: L1 norm heatmap (full vector field)
                ax1.imshow(norm_2d.T, origin='lower', cmap='magma', aspect='auto')
                h, w = norm_2d.shape
                ax1.set_title(f'Full Vector Field ({h}x{w} = {h*w:,} vectors)')
                ax1.set_xlabel('X (pixels)')
                ax1.set_ylabel('Y (pixels)')

                # Right: L1 maxima scatter (downsampled)
                ax2.scatter(local_maxs[0], local_maxs[1], s=0.3, c='steelblue', alpha=0.5)
                ax2.set_xlim(0, h)
                ax2.set_ylim(0, w)
                ax2.invert_yaxis()
                ax2.set_title(f'L1 Maxima ({len(local_maxs[0]):,} pseudo-cells)')
                ax2.set_xlabel('X (pixels)')
                ax2.set_ylabel('Y (pixels)')
                ax2.set_aspect('equal')

                plt.suptitle(f'SSAM Downsampling: Vector Field → Local Maxima ({sample_tag})',
                             fontsize=14, y=1.02)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_downsampling.png"),
                    dpi=150, bbox_inches='tight'
                )
                plt.close()
                logger.info("SSAM downsampling comparison plot saved.")
        except Exception as e:
            logger.warning(f"SSAM downsampling plot failed: {e}")

        # Per-gene KDE heatmap (top 4 expressed genes, PDF slide 62 concept)
        try:
            import dask.array as da
            vf = ds.vf  # shape: (n_genes, width, height) or similar
            if vf is not None:
                # Find top 4 genes by total expression
                gene_names = ssam_adata.var_names.tolist()
                gene_totals = np.array(ssam_adata.X.sum(axis=0)).flatten()
                top_idx = np.argsort(gene_totals)[-4:][::-1]

                fig, axes = plt.subplots(1, 4, figsize=(24, 6))
                for idx, gene_i in enumerate(top_idx):
                    gene_name = gene_names[gene_i]
                    try:
                        gene_field = vf[gene_i]
                        gene_arr = gene_field.compute() if isinstance(gene_field, da.Array) else np.array(gene_field)
                        gene_2d = gene_arr.squeeze()
                        axes[idx].imshow(gene_2d.T, origin='lower', cmap='hot', aspect='auto')
                    except Exception:
                        axes[idx].text(0.5, 0.5, 'N/A', ha='center', va='center',
                                       transform=axes[idx].transAxes)
                    axes[idx].set_title(f'{gene_name}', fontsize=11)
                    axes[idx].set_xlabel('X')
                    axes[idx].set_ylabel('Y')
                plt.suptitle(f'SSAM Per-Gene KDE Fields (top 4 genes, {sample_tag})', fontsize=13)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_gene_kde.png"),
                    dpi=150, bbox_inches='tight'
                )
                plt.close()
                logger.info("SSAM per-gene KDE heatmap saved.")
        except Exception as e:
            logger.warning(f"SSAM per-gene KDE heatmap failed: {e}")

        # Vector field gradient visualization (composite multi-gene KDE + gradient arrows)
        try:
            import dask.array as da
            vf_norm = ds.vf_norm
            if vf_norm is not None:
                norm_arr = vf_norm.compute() if isinstance(vf_norm, da.Array) else np.array(vf_norm)
                norm_2d = norm_arr.squeeze()

                fig, axes = plt.subplots(1, 2, figsize=(22, 10))

                # Left: composite KDE density with gradient overlay
                ax = axes[0]
                im = ax.imshow(norm_2d.T, origin='lower', cmap='inferno', aspect='equal')
                plt.colorbar(im, ax=ax, label='KDE density (L2 norm)', shrink=0.8)

                # Compute gradient and overlay as quiver arrows (subsampled)
                step = max(norm_2d.shape[0] // 30, 1)
                gy, gx = np.gradient(norm_2d)
                Y_grid, X_grid = np.mgrid[0:norm_2d.shape[0]:step, 0:norm_2d.shape[1]:step]
                U = gx[::step, ::step]
                V = gy[::step, ::step]
                mag = np.sqrt(U**2 + V**2)
                mag_safe = np.where(mag > 0, mag, 1)
                ax.quiver(X_grid, Y_grid, V / mag_safe, U / mag_safe,
                          mag, cmap='cool', alpha=0.7, scale=40, width=0.003)
                ax.set_title('KDE Density + Gradient Field', fontsize=13, fontweight='bold')
                ax.set_xlabel('X (px)')
                ax.set_ylabel('Y (px)')

                # Right: top 6 genes as RGB composite
                vf = ds.vf
                if vf is not None:
                    gene_names = ssam_adata.var_names.tolist()
                    gene_totals = np.array(ssam_adata.X.sum(axis=0)).flatten()
                    top6 = np.argsort(gene_totals)[-6:][::-1]
                    # Make 2-row × 3-col panels for top 6 genes
                    # Replace right panel with inset grid
                    axes[1].remove()
                    gs = fig.add_gridspec(2, 3, left=0.55, right=0.98,
                                          top=0.92, bottom=0.08,
                                          wspace=0.15, hspace=0.25)
                    for panel_i, gene_i in enumerate(top6):
                        row, col = divmod(panel_i, 3)
                        ax_g = fig.add_subplot(gs[row, col])
                        gene_name = gene_names[gene_i]
                        try:
                            gene_field = vf[gene_i]
                            gene_arr = gene_field.compute() if isinstance(gene_field, da.Array) else np.array(gene_field)
                            gene_2d = gene_arr.squeeze()
                            ax_g.imshow(gene_2d.T, origin='lower', cmap='magma', aspect='equal')
                        except Exception:
                            ax_g.text(0.5, 0.5, 'N/A', ha='center', va='center',
                                      transform=ax_g.transAxes)
                        ax_g.set_title(gene_name, fontsize=9, fontweight='bold')
                        ax_g.set_xticks([])
                        ax_g.set_yticks([])

                plt.suptitle(f'SSAM Vector Field & Per-Gene KDE ({sample_tag})',
                             fontsize=14, fontweight='bold')
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_vector_field.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close()
                logger.info("SSAM vector field + gradient plot saved.")
        except Exception as e:
            logger.warning(f"SSAM vector field plot failed: {e}")

        # --- Reference-mapped plots ---
        if has_ref:
            # 1. Pixel-level SSAM cell type map
            try:
                palette = mapping_result['palette']
                assignment_names = mapping_result['assignment_names']
                colors = [palette.get(n, (0.5, 0.5, 0.5, 1.0)) for n in assignment_names]
                plt.figure(figsize=(14, 12))
                ds.plot_celltypes_map(colors=colors)
                plt.title(f'SSAM Pixel-Level Cell Type Map ({sample_tag})')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_celltype_map.png"),
                    dpi=150, bbox_inches='tight'
                )
                plt.close()
            except Exception as e:
                logger.warning(f"SSAM pixel-level map plot failed: {e}")

            # 1b. Zoomed ROI cell-type map (dense region crop)
            try:
                import dask.array as da
                vf_norm = ds.vf_norm
                if vf_norm is not None:
                    norm_arr = vf_norm.compute() if isinstance(vf_norm, da.Array) else np.array(vf_norm)
                    norm_2d = norm_arr.squeeze()
                    # Find densest 500x500 region
                    from scipy.ndimage import uniform_filter
                    smoothed = uniform_filter(norm_2d.astype(float), size=min(500, min(norm_2d.shape) // 2))
                    peak = np.unravel_index(np.argmax(smoothed), smoothed.shape)
                    roi_size = min(500, min(norm_2d.shape) // 3)
                    half = roi_size // 2
                    x0 = max(0, min(peak[0] - half, norm_2d.shape[0] - roi_size))
                    y0 = max(0, min(peak[1] - half, norm_2d.shape[1] - roi_size))

                    # Render full map to get image data, then compose 2-panel manually
                    plt.figure(figsize=(14, 12))
                    ds.plot_celltypes_map(colors=colors)
                    # Get rendered image from current figure
                    fig_tmp = plt.gcf()
                    fig_tmp.canvas.draw()
                    buf = np.frombuffer(fig_tmp.canvas.buffer_rgba(), dtype=np.uint8)
                    buf = buf.reshape(fig_tmp.canvas.get_width_height()[::-1] + (4,))
                    plt.close()

                    # Compose 2-panel: full + zoomed
                    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(20, 10))
                    ax_full.imshow(buf)
                    ax_full.set_title(f'Full Cell Type Map ({sample_tag})')
                    ax_full.axis('off')

                    # For zoom: re-render and crop
                    plt.figure(figsize=(14, 12))
                    ds.plot_celltypes_map(colors=colors)
                    ax_tmp = plt.gca()
                    ax_tmp.set_xlim(y0, y0 + roi_size)
                    ax_tmp.set_ylim(x0 + roi_size, x0)
                    fig_tmp2 = plt.gcf()
                    fig_tmp2.canvas.draw()
                    buf2 = np.frombuffer(fig_tmp2.canvas.buffer_rgba(), dtype=np.uint8)
                    buf2 = buf2.reshape(fig_tmp2.canvas.get_width_height()[::-1] + (4,))
                    plt.close()

                    ax_zoom.imshow(buf2)
                    ax_zoom.set_title(f'Zoomed ROI ({roi_size}×{roi_size} px)')
                    ax_zoom.axis('off')

                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(output_dir, f"{sample_tag}_step2_ssam_celltype_map_zoomed.png"),
                        dpi=150, bbox_inches='tight'
                    )
                    plt.close()
                    logger.info("SSAM zoomed cell-type map saved.")
            except Exception as e:
                logger.warning(f"SSAM zoomed cell-type map failed: {e}")

            # 2. UMAP colored by leiden_assignment (Step 1 style)
            try:
                sc.pl.umap(
                    ssam_adata, color=['leiden_assignment'],
                    size=8,
                    legend_loc='on data',
                    legend_fontsize=8,
                    legend_fontoutline=2,
                    title=f'SSAM Cell Types ({sample_tag})',
                    frameon=False, show=False,
                )
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_umap_celltypes.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close()
            except Exception as e:
                logger.warning(f"SSAM celltype UMAP plot failed: {e}")

            # 3. Correlation heatmap (leiden clusters × reference celltypes) — paper style
            try:
                corr_matrix = mapping_result['corr_matrix']
                corr_float = corr_matrix.astype(float)

                # Highlight best match per cluster (column)
                best_mask = pd.DataFrame(False, index=corr_float.index, columns=corr_float.columns)
                for col in corr_float.columns:
                    best_row = corr_float[col].idxmax()
                    best_mask.loc[best_row, col] = True

                n_clusters = len(corr_float.columns)
                n_types = len(corr_float.index)
                fig, ax = plt.subplots(figsize=(max(12, n_clusters * 0.7 + 2),
                                                max(6, n_types * 0.55 + 2)))

                # Annotation: bold best match per cluster
                annot_arr = corr_float.values.copy()
                annot_strs = np.empty_like(annot_arr, dtype=object)
                for i in range(annot_arr.shape[0]):
                    for j in range(annot_arr.shape[1]):
                        val = annot_arr[i, j]
                        if best_mask.iloc[i, j]:
                            annot_strs[i, j] = f'{val:.2f}*'
                        else:
                            annot_strs[i, j] = f'{val:.2f}'

                sns.heatmap(
                    corr_float, annot=annot_strs, fmt='',
                    cmap='RdBu_r', center=0, vmin=-0.5, vmax=1.0,
                    linewidths=0.5, linecolor='white',
                    xticklabels=True, yticklabels=True, ax=ax,
                    cbar_kws={'label': 'Pearson r', 'shrink': 0.8},
                )
                ax.set_xlabel('Leiden Cluster', fontsize=12, fontweight='bold')
                ax.set_ylabel('Reference Cell Type', fontsize=12, fontweight='bold')
                ax.set_title(f'SSAM Leiden–Reference Correlation ({sample_tag})',
                             fontsize=14, fontweight='bold', pad=12)
                ax.tick_params(axis='x', rotation=0, labelsize=10)
                ax.tick_params(axis='y', rotation=0, labelsize=10)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(output_dir, f"{sample_tag}_step2_ssam_correlation.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close()
                logger.info("SSAM correlation heatmap saved (paper style).")
            except Exception as e:
                logger.warning(f"SSAM correlation heatmap failed: {e}")

        logger.info("SSAM plots saved.")
    except Exception as e:
        logger.warning(f"SSAM plotting failed: {e}")


# ============================================================
# 2-9. Boundary overlay on spatial maps (Ref: 2_3_brain_cell_overlaps)
# ============================================================

def plot_boundary_spatial_overlay(adata, config, output_dir, sample_tag):
    """Spatial map with nucleus boundary polygons overlaid.

    Matching notebook 2_3: plots boundary contours + centroids on a spatial
    scatter background coloured by cell type or P2R cluster.

    Generates two figures:
    1. Transcript scatter + boundary polygons (full tissue, downsampled)
    2. Zoomed ROI with boundary detail + centroids
    """
    bd_path = config.get('nucleus_boundaries_path')
    if not bd_path:
        # Try loading from input_path
        input_dir = config.get('input_path', '')
        for fname in ['nucleus_boundaries.parquet', 'nucleus_boundaries.csv.gz']:
            candidate = os.path.join(input_dir, fname)
            if os.path.exists(candidate):
                bd_path = candidate
                break
    if not bd_path or not os.path.exists(bd_path):
        logger.info("No boundary file available; skipping boundary overlay.")
        return

    try:
        if bd_path.endswith('.parquet'):
            boundaries = pd.read_parquet(bd_path)
        else:
            boundaries = pd.read_csv(bd_path)
    except Exception as e:
        logger.warning(f"Failed to load boundaries for overlay: {e}")
        return

    if 'vertex_x' not in boundaries.columns or 'vertex_y' not in boundaries.columns:
        logger.warning("Boundary file missing vertex_x/vertex_y columns.")
        return

    spots = adata.uns.get('spots')
    if spots is None or 'x_location' not in spots.columns:
        logger.warning("No spots with x_location for boundary overlay.")
        return

    # --- Figure 1: Full tissue boundary overview ---
    try:
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))

        # Panel 1: Transcripts coloured by cell type
        ct_col = None
        for col in ['_celltype_mapped', 'celltype', 'Class', 'leiden']:
            if col in spots.columns:
                ct_col = col
                break

        # Subsample transcripts for speed
        n_tx = min(len(spots), 100000)
        sub = spots.sample(n=n_tx, random_state=42)

        if ct_col and ct_col in sub.columns:
            cats = sub[ct_col].astype(str)
            unique_cats = cats.unique()
            cmap = plt.cm.get_cmap('tab20', len(unique_cats))
            color_map = {c: cmap(i) for i, c in enumerate(unique_cats)}
            colors = [color_map.get(c, (0.5, 0.5, 0.5, 0.5)) for c in cats]
            axes[0].scatter(sub['x_location'], sub['y_location'],
                            s=0.1, alpha=0.3, c=colors, rasterized=True)
            axes[0].set_title(f'Transcripts by {ct_col} ({n_tx:,} shown)')
        else:
            axes[0].scatter(sub['x_location'], sub['y_location'],
                            s=0.1, alpha=0.2, c='steelblue', rasterized=True)
            axes[0].set_title(f'Transcripts ({n_tx:,} shown)')
        axes[0].invert_yaxis()
        axes[0].set_aspect('equal')
        axes[0].axis('off')

        # Panel 2: Transcripts + boundary polygons (subsample boundaries)
        axes[1].scatter(sub['x_location'], sub['y_location'],
                        s=0.05, alpha=0.1, c='gray', rasterized=True)

        unique_cells = boundaries['cell_id'].unique()
        n_bd = min(len(unique_cells), 1000)
        np.random.seed(42)
        sample_cells = np.random.choice(unique_cells, n_bd, replace=False)

        for cid in sample_cells:
            cell_bd = boundaries[boundaries['cell_id'] == cid]
            axes[1].plot(cell_bd['vertex_x'], cell_bd['vertex_y'],
                         c='cyan', linewidth=0.3, alpha=0.5)

        axes[1].set_title(f'Nucleus boundaries ({n_bd:,} / {len(unique_cells):,})')
        axes[1].invert_yaxis()
        axes[1].set_aspect('equal')
        axes[1].axis('off')

        fig.suptitle('Spatial transcript map with nucleus boundaries', fontsize=14)
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step2_boundary_overview.png")
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved boundary overview: {path}")
    except Exception as e:
        logger.warning(f"Boundary overview plot failed: {e}")
        plt.close('all')

    # --- Figure 2: Zoomed ROI with boundary detail + centroids ---
    try:
        x_all = spots['x_location'].values
        y_all = spots['y_location'].values
        x_med, y_med = np.nanmedian(x_all), np.nanmedian(y_all)

        # 500 µm x 500 µm ROI
        roi_half = 250
        x_min, x_max = x_med - roi_half, x_med + roi_half
        y_min, y_max = y_med - roi_half, y_med + roi_half

        # Filter spots in ROI
        roi_mask = ((x_all >= x_min) & (x_all <= x_max) &
                    (y_all >= y_min) & (y_all <= y_max))
        roi_spots = spots[roi_mask]

        # Filter boundaries in ROI
        bd_in_roi = boundaries[
            (boundaries['vertex_x'] >= x_min) & (boundaries['vertex_x'] <= x_max) &
            (boundaries['vertex_y'] >= y_min) & (boundaries['vertex_y'] <= y_max)
        ]
        roi_cells = bd_in_roi['cell_id'].unique()

        fig, axes = plt.subplots(1, 3, figsize=(21, 7))

        # Panel 1: Transcripts in ROI
        axes[0].scatter(roi_spots['x_location'], roi_spots['y_location'],
                        s=1, alpha=0.3, c='steelblue', rasterized=True)
        axes[0].set_title(f'Transcripts ({len(roi_spots):,})')
        axes[0].set_xlim(x_min, x_max)
        axes[0].set_ylim(y_max, y_min)
        axes[0].set_aspect('equal')
        axes[0].axis('off')

        # Panel 2: Boundaries + transcripts
        axes[1].scatter(roi_spots['x_location'], roi_spots['y_location'],
                        s=0.5, alpha=0.15, c='gray', rasterized=True)
        for cid in roi_cells:
            cell_bd = boundaries[boundaries['cell_id'] == cid]
            axes[1].plot(cell_bd['vertex_x'], cell_bd['vertex_y'],
                         c='cyan', linewidth=0.8, alpha=0.7)
        axes[1].set_title(f'Boundaries ({len(roi_cells):,} nuclei)')
        axes[1].set_xlim(x_min, x_max)
        axes[1].set_ylim(y_max, y_min)
        axes[1].set_aspect('equal')
        axes[1].axis('off')

        # Panel 3: Boundaries + centroids
        for cid in roi_cells:
            cell_bd = boundaries[boundaries['cell_id'] == cid]
            axes[2].plot(cell_bd['vertex_x'], cell_bd['vertex_y'],
                         c='cyan', linewidth=0.5, alpha=0.5)

        # Plot centroids if available
        if 'x_centroid' in adata.obs.columns:
            centroids = adata.obs[['x_centroid', 'y_centroid']].astype(float)
            c_roi = centroids[
                (centroids['x_centroid'] >= x_min) & (centroids['x_centroid'] <= x_max) &
                (centroids['y_centroid'] >= y_min) & (centroids['y_centroid'] <= y_max)
            ]
            axes[2].scatter(c_roi['x_centroid'], c_roi['y_centroid'],
                            s=8, c='red', marker='x', linewidths=0.5, zorder=5)
            axes[2].set_title(f'Boundaries + centroids ({len(c_roi):,})')
        else:
            axes[2].set_title(f'Boundaries (no centroids)')
        axes[2].set_xlim(x_min, x_max)
        axes[2].set_ylim(y_max, y_min)
        axes[2].set_aspect('equal')
        axes[2].axis('off')

        fig.suptitle(f'Zoomed ROI (500 µm × 500 µm) — Nucleus boundaries', fontsize=14)
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step2_boundary_zoomed_roi.png")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved boundary zoomed ROI: {path}")
    except Exception as e:
        logger.warning(f"Boundary zoomed ROI failed: {e}")
        plt.close('all')


# ============================================================
# Entry Point
# ============================================================

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
        logger.warning("'spots' not found in adata.uns. Reloading from raw CSV...")
        try:
            transcripts_path = os.path.join(input_dir, "transcripts.csv")
            df = pd.read_csv(transcripts_path)
            if 'feature_name' not in df.columns:
                if 'gene' in df.columns:
                    df.rename(columns={'gene': 'feature_name'}, inplace=True)
            adata.uns['spots'] = df
        except Exception as e:
            logger.error(f"Failed to load transcripts: {e}")
            raise

    # Pre-populate celltype on spots via cell_id → adata.obs mapping (matches notebook pattern)
    spots = adata.uns['spots']
    if 'cell_id' in spots.columns:
        # Find best leiden/celltype column in adata.obs
        ct_col = None
        for c in ['Class', 'celltype', 'cell_type', 'leiden'] + sorted(
                [c for c in adata.obs.columns if c.startswith('leiden_')]):
            if c in adata.obs.columns:
                ct_col = c
                break
        if ct_col is not None:
            spots_cid = _normalize_cell_ids(spots['cell_id'])
            if 'cell_id' in adata.obs.columns:
                obs_cid = _normalize_cell_ids(adata.obs['cell_id'])
            else:
                obs_cid = _normalize_cell_ids(adata.obs.index.to_series())
            cid_to_ct = dict(zip(obs_cid, adata.obs[ct_col]))
            spots['celltype'] = spots_cid.map(cid_to_ct).fillna('Background')
            n_mapped = (spots['celltype'] != 'Background').sum()
            logger.info(f"Pre-mapped {n_mapped}/{len(spots)} spots to cell types via '{ct_col}'")
            adata.uns['spots'] = spots

    # ==========================================================
    # 2-1. Points2Regions
    # ==========================================================
    if config['segmentation_free'].get('run_points2regions', True):
        logger.info("Running Points2Regions...")
        spots = adata.uns['spots']
        xy = spots[['x_location', 'y_location']].values
        genes = spots['feature_name'].values

        p2r_params = config['segmentation_free']['points2regions']

        unique_genes = np.unique(genes)
        gene_map = {g: i for i, g in enumerate(unique_genes)}
        gene_labels = np.array([gene_map[g] for g in genes])

        sigma = p2r_params.get('sigma', 3.0)

        n_clusters_cfg = p2r_params.get('n_clusters', 100)
        n_clusters_list = n_clusters_cfg if isinstance(n_clusters_cfg, list) else [n_clusters_cfg]

        primary_p2r_adata = None
        for nc in n_clusters_list:
            logger.info(f"  Points2Regions: fitting k={nc}")
            p2r = Points2Regions(
                xy, gene_labels,
                pixel_width=sigma / 3.0,
                pixel_smoothing=sigma,
                min_num_pts_per_pixel=p2r_params.get('min_genes_per_bin', 15),
            )
            p2r_adata = p2r.fit_predict(
                num_clusters=int(nc),
                output='anndata',
                seed=42,
                adata_cluster_key=f'points2regions_{nc}',
            )

            p2r_path = os.path.join(output_dir, f"{sample_tag}_step2_points2regions_k{nc}_bins.h5ad")
            p2r_adata.write_h5ad(p2r_path)
            logger.info(f"  Saved P2R (k={nc}) to {p2r_path}")

            if primary_p2r_adata is None:
                primary_p2r_adata = p2r_adata

        # Store primary clusters
        if primary_p2r_adata is not None:
            primary_key = f'points2regions_{n_clusters_list[0]}'
            cluster_col = primary_p2r_adata.uns['reads'][primary_key]
            adata.uns['spots']['points2regions'] = cluster_col.values
            logger.info(f"Primary P2R (k={n_clusters_list[0]}) stored in spots")

        # --- 2-1b. Subcellular classification ---
        adata = classify_p2r_subcellular(adata, output_dir, sample_tag)
        assign_p2r_colors(adata)

    # ==========================================================
    # 2-2. Overlaps Analysis (ovrlpy)
    # ==========================================================
    if config['segmentation_free'].get('run_overlaps', False):
        try:
            adata = run_ovrlpy_analysis(adata, config, output_dir, sample_tag)
        except Exception as e:
            logger.warning(f"ovrlpy analysis failed: {e}")

    # ==========================================================
    # 2-3. Distance Metrics
    # ==========================================================
    logger.info(f"Distance metric mode: '{distance_metric}'")

    if distance_metric in ('centroid', 'both'):
        adata = calculate_distance_to_centroid(adata, input_dir)

    if distance_metric in ('boundary', 'both'):
        adata = calculate_distance_to_boundary(adata, input_dir)

    # --- Checkpoint ---
    output_file = os.path.join(output_dir, f"{sample_tag}_step2_points2regions.h5ad")
    adata.write_h5ad(output_file)
    logger.info(f"  Checkpoint: adata saved to {output_file}")

    # ==========================================================
    # Visualization (failures won't lose data)
    # ==========================================================

    # --- 2-4. Centroid distance plots ---
    try:
        if distance_metric in ('centroid', 'both'):
            plot_centroid_distance_analysis(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"Centroid distance plots failed: {e}")

    # --- 2-1c. P2R heatmap ---
    try:
        plot_p2r_heatmap(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R heatmap failed: {e}")

    # --- 2-1b2. P2R spatial map (Fig 1i / Ext Data Fig 3A) ---
    try:
        plot_p2r_spatial_map(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R spatial map failed: {e}")

    # --- 2-1d. P2R top genes ---
    try:
        plot_p2r_top_genes(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R top genes failed: {e}")

    # --- 2-1e. P2R DE marker extraction & violin ---
    try:
        extract_p2r_de_markers(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R DE marker extraction failed: {e}")

    # --- 2-5. P2R boundary distance ---
    try:
        if distance_metric in ('boundary', 'both'):
            plot_p2r_boundary_distance(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R boundary distance plots failed: {e}")

    # --- 2-5b. P2R mean distance scatter ---
    try:
        if distance_metric in ('boundary', 'both'):
            plot_p2r_mean_distance_scatter(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"P2R mean distance scatter failed: {e}")

    # --- 2-6. Summary distance histograms ---
    try:
        _plot_distance_histograms(adata, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"Distance histograms failed: {e}")

    # --- 2-9. Boundary spatial overlay (Ref: 2_3_brain_cell_overlaps) ---
    try:
        plot_boundary_spatial_overlay(adata, config, output_dir, sample_tag)
    except Exception as e:
        logger.warning(f"Boundary spatial overlay failed: {e}")

    # ==========================================================
    # 2-7. SSAM Analysis (optional)
    # ==========================================================
    if config['segmentation_free'].get('run_ssam', False):
        try:
            adata = run_ssam_analysis(adata, config, output_dir, sample_tag)
        except Exception as e:
            logger.warning(f"SSAM analysis failed: {e}")

    # ==========================================================
    # 2-8. Final save
    # ==========================================================
    adata.write_h5ad(output_file)
    logger.info(f"Step 2 Completed. Final output: {output_file}")

    return adata


def _plot_distance_histograms(adata, output_dir, sample_tag):
    """Summary histograms for boundary and centroid distances."""
    spots = adata.uns['spots']

    if 'dist_to_nucleus' in spots.columns:
        dists = spots['dist_to_nucleus'].dropna()
        if len(dists) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # Signed
            axes[0].hist(dists, bins=100, log=True, edgecolor='k', alpha=0.7)
            axes[0].axvline(0, color='r', linestyle='--', label='Nuclear edge')
            axes[0].set_title('Signed Distance to Nucleus Boundary')
            axes[0].set_xlabel('Distance (um, negative=inside)')
            axes[0].set_ylabel('Count (log)')
            axes[0].legend()

            # Unsigned
            axes[1].hist(np.abs(dists), bins=100, log=True, edgecolor='k', alpha=0.7)
            axes[1].set_title('Absolute Distance to Nucleus Boundary')
            axes[1].set_xlabel('Distance (um)')
            axes[1].set_ylabel('Count (log)')

            plt.tight_layout()
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_distance_boundary_dist.png"),
                dpi=150
            )
            plt.close()

    if 'dist_to_centroid' in spots.columns:
        dists = spots['dist_to_centroid'].dropna()
        if len(dists) > 0:
            plt.figure(figsize=(10, 6))
            plt.hist(dists, bins=100, log=True, edgecolor='k', alpha=0.7)
            plt.title('Distribution of Transcript Distances to Centroid')
            plt.xlabel('Distance to Centroid (um)')
            plt.ylabel('Count (log)')
            plt.savefig(
                os.path.join(output_dir, f"{sample_tag}_step2_distance_centroid_dist.png"),
                dpi=150
            )
            plt.close()


if __name__ == "__main__":
    pass
