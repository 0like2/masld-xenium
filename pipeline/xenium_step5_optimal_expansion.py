# Step 5: Optimal Expansion (Ref: notebooks/4_optimal_expansion/4_1)
# Assigns unassigned transcripts to nearest annotated cell domain via KDTree,
# then computes correlation-based turnover to find optimal expansion radius.
#
# Flow:
#   5-1. Load original transcripts (Step 0)
#   5-2. Load annotated cells (Step 1/2/4)
#   5-3. Map domain assignments to reads
#   5-4. KDTree nearest-domain expansion
#   5-5. Distance threshold filter (optional)
#   5-6. Spatial map + distance histogram
#   5-7. Save expanded transcripts
#   5-8. Turnover analysis
#     5-8a. Nuclear vs background expression profiles
#     5-8b. Distance-bin correlation curves
#     5-8c. Turnover distance detection
#     5-8d. Optimal expansion = turnover − nuclei_size
#     5-8e. Summary barplot + CSVs

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, ConvexHull
import random as rd

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def dist_nuc(reads_ctdsub):
    """Median distance to nuclei edges via ConvexHull vertices."""
    allds = []
    for g, n in reads_ctdsub.groupby('cell_id'):
        try:
            if len(n) < 3:
                continue
            hull = ConvexHull(np.array(n.loc[:, ['x_location', 'y_location']]))
            if 'distance' not in n.columns:
                logger.warning("dist_nuc(): 'distance' column missing for cell_id=%s", g)
                continue
            allds.append(np.mean(n.iloc[hull.vertices]['distance']))
        except Exception:
            pass
    if len(allds) > 0:
        median_dist = np.median(allds)
    else:
        median_dist = np.nan
    return median_dist


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def _load_p2r_spatial_domains(adata, output_dir, sample_tag, n_clusters=50, max_dist_um=500.0):
    """Map cells to nearest P2R bin cluster as spatial domain proxy.

    Loads P2R bin-level h5ad from Step 2, maps each cell to its nearest
    bin via KDTree, assigns the bin's cluster label as the cell's domain.
    Cells farther than max_dist_um from any bin get NaN.

    If n_clusters is smaller than the smallest available P2R file (e.g., 8
    requested but only k=50 exists), the bins are spatially meta-clustered:
    KMeans on bin spatial coordinates groups nearby bins into n_clusters
    contiguous spatial regions.  This produces anatomically meaningful
    domains where each region spans multiple cell types, yielding a clear
    gap between nuclear and background expression profiles for turnover
    detection.
    """
    parent_dir = os.path.dirname(output_dir)
    p2r_dir = os.path.join(parent_dir, "step2_segmentation_free")

    # --- Find the best available P2R file ---
    p2r_path = os.path.join(p2r_dir,
                            f"{sample_tag}_step2_points2regions_k{n_clusters}_bins.h5ad")
    source_k = n_clusters  # cluster count in the file we load
    needs_metaclustering = False

    if not os.path.exists(p2r_path):
        # Search for available P2R files and pick the smallest k >= n_clusters,
        # or the smallest k available if all are larger than n_clusters.
        import glob as _glob
        candidates = sorted(_glob.glob(
            os.path.join(p2r_dir, f"{sample_tag}_step2_points2regions_k*_bins.h5ad")))
        available_ks = []
        for c in candidates:
            try:
                k_str = os.path.basename(c).split('_k')[1].split('_bins')[0]
                available_ks.append((int(k_str), c))
            except (IndexError, ValueError):
                continue

        if not available_ks:
            print(f"    - No P2R bins files found in {p2r_dir}")
            return None

        available_ks.sort(key=lambda x: x[0])
        # Prefer the smallest available k (finest resolution for meta-clustering)
        source_k, p2r_path = available_ks[0]
        needs_metaclustering = (n_clusters < source_k)
        print(f"    - Requested k={n_clusters} not found. Using k={source_k} "
              f"{'(will meta-cluster → ' + str(n_clusters) + ' domains)' if needs_metaclustering else ''}")

    print(f"    - Loading P2R bins from: {p2r_path}")
    p2r_bins = sc.read_h5ad(p2r_path)
    bin_coords = p2r_bins.obsm['spatial']  # (n_bins, 2)

    cluster_key = f'points2regions_{source_k}'
    if cluster_key not in p2r_bins.obs.columns:
        # Try first matching column
        candidates_cols = [c for c in p2r_bins.obs.columns if 'points2regions' in c]
        if not candidates_cols:
            print(f"    - No points2regions column found in P2R bins. Available: {list(p2r_bins.obs.columns)}")
            return None
        cluster_key = candidates_cols[0]

    bin_clusters = p2r_bins.obs[cluster_key].values.copy()

    # --- Spatial meta-clustering: merge bins → n_clusters by position ---
    if needs_metaclustering and n_clusters < len(np.unique(bin_clusters)):
        from sklearn.cluster import KMeans

        n_orig = len(np.unique(bin_clusters))
        print(f"    - Spatial meta-clustering {n_orig} P2R clusters → {n_clusters} domains "
              f"(KMeans on bin spatial coordinates)")

        # Cluster bins by their spatial position, not expression.
        # This produces contiguous spatial regions where each domain
        # contains a diverse mix of cell types — critical for turnover.
        km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        spatial_labels = km.fit_predict(bin_coords)

        bin_clusters = np.array([f"region_{s}" for s in spatial_labels])

        # Report cluster sizes
        meta_sizes = pd.Series(spatial_labels).value_counts().sort_index()
        print(f"    - Spatial domain sizes (# bins each): "
              f"{dict(zip([f'region_{i}' for i in meta_sizes.index], meta_sizes.values))}")

    # Get cell coordinates: obsm['spatial'] or obs x/y_centroid columns
    if 'spatial' in adata.obsm:
        cell_coords = adata.obsm['spatial']  # (n_cells, 2) — µm
    elif 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        cell_coords = np.column_stack([
            adata.obs['x_centroid'].values.astype(float),
            adata.obs['y_centroid'].values.astype(float)
        ])
    else:
        print("    - No spatial coordinates found in adata (obsm['spatial'] or obs x/y_centroid)")
        return None

    # Auto-detect pixel vs µm mismatch: if cell range >> bin range, convert px→µm
    cell_range = cell_coords[:, 0].max() - cell_coords[:, 0].min()
    bin_range = bin_coords[:, 0].max() - bin_coords[:, 0].min()
    if cell_range > bin_range * 2:
        scale = cell_range / bin_range
        print(f"    - Coordinate mismatch detected (cell_range={cell_range:.0f}, "
              f"bin_range={bin_range:.0f}, ratio={scale:.2f}x). Converting cell coords px→µm.")
        cell_coords = cell_coords / scale

    tree = cKDTree(bin_coords)
    dists, idxs = tree.query(cell_coords, k=1)

    domains = pd.Series(bin_clusters[idxs], index=adata.obs.index, dtype=str)
    domains[dists > max_dist_um] = np.nan  # Too far from any P2R bin

    n_assigned = domains.notna().sum()
    n_unique = len(np.unique(bin_clusters))
    print(f"    - P2R domain mapping: {n_assigned}/{len(domains)} cells assigned "
          f"({n_assigned/len(domains):.1%}), {n_unique} domains")

    return domains


# --- Correlation-based Turnover / Crossover Analysis ---

def calculate_turnover(reads_assigned, reads_not_assigned, output_dir, sample_tag,
                       min_reads_per_domain=5000, diff_threshold=0.1,
                       min_reads_per_bin=1, celltype_colors=None,
                       threshold_mode="absolute"):
    """
    Correlation-based turnover/crossover for each cell-type x domain.

    Algorithm per (cell-type, domain) pair with enough reads:
      1. Build distance x gene crosstab from assigned reads.
      2. Get nuclear-expression profile (overlaps_nucleus == 1).
      3. Get background profile from unassigned reads.
      4. Correlate expression at each distance bin with nuclear and background.
      5. Turnover distance detection (mode-dependent, see threshold_mode).
      6. Compute nuclei_size and cell_size via dist_nuc.

    Parameters
    ----------
    celltype_colors : dict or None
        Mapping from cell type name to color (hex). If available from
        adata.uns['{celltype_key}_colors'], passed here for consistent plots.
    threshold_mode : str
        Turnover detection mode:
        - "absolute": diff < diff_threshold (notebook method, default)
        - "crossover": diff < 0 (paper method, exact intersection)
        - "relative": diff < diff_threshold * peak_diff (half-life method)

    Returns (turnover_summary, per_celltype_df, optimal_expansion_value).
    """
    print("\n    [Turnover] Starting correlation-based turnover analysis ...")

    # Ensure overlaps_nucleus exists; use distance-based proxy if missing
    has_overlap_col = 'overlaps_nucleus' in reads_assigned.columns
    used_proxy = False
    proxy_threshold = None
    if not has_overlap_col:
        print("    [Turnover] 'overlaps_nucleus' column missing -- using distance-based proxy.")
        med_dist = reads_assigned['distance'].median()
        reads_assigned = reads_assigned.copy()
        reads_assigned['overlaps_nucleus'] = (reads_assigned['distance'] < med_dist).astype(int)
        print(f"    [Turnover] Proxy threshold (median distance): {med_dist:.2f}")
        used_proxy = True
        proxy_threshold = med_dist

    # Background expression from unassigned reads
    reads_not_assigned_clean = reads_not_assigned.dropna(subset=['domain', 'feature_name'])
    if len(reads_not_assigned_clean) == 0:
        print("    [Turnover] No usable background reads -- aborting turnover calculation.")
        return None, None, np.nan

    # --- 5-8a. Nuclear vs background expression profiles ---
    background_express = pd.crosstab(
        reads_not_assigned_clean['domain'],
        reads_not_assigned_clean['feature_name']
    )

    unique_celltypes = reads_assigned['initial_annotation'].dropna().unique()
    unique_domains = reads_assigned['domain'].dropna().unique()

    turnover_summ = pd.DataFrame(
        columns=unique_celltypes,
        index=unique_domains,
        dtype=float
    )

    meand_celltype = []
    nuclimall = []
    cellmall = []
    ct_labels = []

    figures_dir = os.path.join(output_dir, "crossover_plots")
    os.makedirs(figures_dir, exist_ok=True)

    # --- Main loop: cell type x domain ---
    for celltype in unique_celltypes:
        if pd.isna(celltype):
            continue
        print(f"    [Turnover] Processing cell type: {celltype}")

        reads_ct = reads_assigned[reads_assigned['initial_annotation'] == celltype]

        domain_counts = reads_ct.groupby('domain').size()
        good_domains = domain_counts[domain_counts > min_reads_per_domain].index

        if len(good_domains) == 0:
            print(f"      - No domains with > {min_reads_per_domain} reads. Skipping.")
            continue

        tdistanceall = []
        nuclim = []
        cellm = []

        for dom in good_domains:
            reads_ctd = reads_ct[reads_ct['domain'] == dom]

            # 1. Distance x gene crosstab
            expression_distances = pd.crosstab(
                reads_ctd['distance'], reads_ctd['feature_name']
            )

            # 2. Nuclear expression profile
            reads_ctd_nucl = pd.crosstab(
                reads_ctd['overlaps_nucleus'], reads_ctd['feature_name']
            )
            if 1 not in reads_ctd_nucl.index:
                continue

            # 3. Background expression for this domain
            if dom not in background_express.index:
                continue
            bck_sub = background_express

            # Align columns across all three matrices
            common_genes = (
                set(reads_ctd_nucl.columns) &
                set(bck_sub.columns) &
                set(expression_distances.columns)
            )
            if len(common_genes) < 5:
                continue

            common_genes = sorted(common_genes)

            # Remove dominant housekeeping genes that compress correlation range.
            # Genes like MTRNR2L12/MTRNR2L8 can account for >50% of reads,
            # making all profiles converge to ~0.97 correlation regardless of
            # cell-type composition.  Excluding them widens the dynamic range
            # from ~0.03 to ~0.2, enabling reliable turnover detection.
            bck_total = bck_sub.loc[dom, common_genes].values.astype(float)
            bck_sum = bck_total.sum()
            if bck_sum > 0:
                bck_frac = bck_total / bck_sum
                dominant_mask = bck_frac > 0.15  # genes > 15% of background reads
                if dominant_mask.sum() > 0:
                    keep_genes = [g for g, m in zip(common_genes, dominant_mask) if not m]
                    if len(keep_genes) >= 20:
                        common_genes = keep_genes

            reads_ctd_nucl = reads_ctd_nucl[common_genes]
            bck_sub = bck_sub[common_genes]
            expression_distances = expression_distances[common_genes]

            nuc_vec = reads_ctd_nucl.loc[1, :].values.astype(float)
            bck_vec = bck_sub.loc[dom, :].values.astype(float)

            # --- 5-8b. Distance-bin correlation curves ---
            corr_nuc_list = []
            corr_bck_list = []
            for dist in expression_distances.index:
                dist_vec = expression_distances.loc[dist, :].values.astype(float)
                n_reads_at_dist = dist_vec.sum()
                # Skip bins with too few reads — correlation is unreliable
                if n_reads_at_dist < min_reads_per_bin:
                    corr_nuc_list.append(np.nan)
                    corr_bck_list.append(np.nan)
                    continue
                if np.std(dist_vec) == 0 or np.std(nuc_vec) == 0:
                    corr_nuc_list.append(np.nan)
                else:
                    corr_nuc_list.append(np.corrcoef(dist_vec, nuc_vec)[0, 1])
                if np.std(dist_vec) == 0 or np.std(bck_vec) == 0:
                    corr_bck_list.append(np.nan)
                else:
                    corr_bck_list.append(np.corrcoef(dist_vec, bck_vec)[0, 1])

            summary = pd.DataFrame({
                'corr_nuc': corr_nuc_list,
                'corr_back': corr_bck_list
            }, index=expression_distances.index)
            summary.index.name = 'distance'
            summary['corr_nuc'] = summary['corr_nuc'].astype(float)
            summary['corr_back'] = summary['corr_back'].astype(float)
            summary['diff'] = summary['corr_nuc'] - summary['corr_back']

            # Per-domain crossover plot
            try:
                fig, ax = plt.subplots(figsize=(6, 4))
                summary_plot = summary.reset_index()
                sns.lineplot(data=summary_plot, x='distance', y='corr_nuc', label='corr_nuc', ax=ax)
                sns.lineplot(data=summary_plot, x='distance', y='corr_back', label='corr_back', ax=ax)
                dom_str = str(dom).replace('/', '_')
                ct_str = str(celltype).replace('/', '_')
                ax.set_title(f"Domain {dom}, celltype {celltype}")
                ax.set_ylabel("Correlation")
                ax.set_xlabel("Distance (rounded)")
                fig_path = os.path.join(figures_dir, f"{dom_str}_{ct_str}.png")
                fig.savefig(fig_path, dpi=150, bbox_inches='tight')
                plt.close(fig)
            except Exception as e:
                logger.debug(f"      Plot failed for {dom}/{celltype}: {e}")
                plt.close('all')

            # --- 5-8c. Turnover distance detection ---
            summary_valid = summary.dropna(subset=['diff'])

            n_valid_bins = len(summary_valid)
            n_total_bins = len(summary)

            if n_valid_bins < 3:
                tdistance = np.nan
                print(f"      - Domain {dom}: {n_valid_bins}/{n_total_bins} valid bins — too few, skipping")
            else:
                # Light smoothing (window=3) for noise reduction
                diff_smooth = summary_valid['diff'].rolling(
                    window=3, min_periods=1, center=True).mean()

                peak_diff = diff_smooth.max()
                if peak_diff <= 0.01:
                    # No detectable nuclear-vs-background signal
                    tdistance = np.nan
                    print(f"      - Domain {dom}: {n_valid_bins}/{n_total_bins} valid bins, "
                          f"peak_diff={peak_diff:.3f} — no signal, skipping")
                else:
                    peak_dist = diff_smooth.idxmax()
                    # Only look AFTER the peak for the drop
                    after_peak = diff_smooth.loc[summary_valid.index >= peak_dist]

                    if threshold_mode == "crossover":
                        # Paper method: exact intersection where corr_nuc < corr_back
                        below_after = after_peak < 0
                        effective_thresh = 0
                    elif threshold_mode == "relative":
                        # Half-life method: diff < diff_threshold * peak_diff
                        effective_thresh = diff_threshold * peak_diff
                        below_after = after_peak < effective_thresh
                    else:
                        # "absolute" (default, notebook method): diff < diff_threshold
                        effective_thresh = diff_threshold
                        below_after = after_peak < effective_thresh

                    if below_after.any():
                        tdistance = below_after.index[below_after.values].min()
                    else:
                        tdistance = summary_valid.index.max()

                    print(f"      - Domain {dom}: {n_valid_bins}/{n_total_bins} valid bins, "
                          f"peak_diff={peak_diff:.3f}@d={peak_dist}, "
                          f"mode={threshold_mode}, thresh={effective_thresh:.4f}, turnover={tdistance}")

            turnover_summ.loc[dom, celltype] = tdistance
            tdistanceall.append(tdistance)

            # 6. nuclei_size and cell_size
            reads_ctdsub_nuc = reads_ctd[reads_ctd['overlaps_nucleus'] == 1]
            cellm.append(dist_nuc(reads_ctd))
            nuclim.append(dist_nuc(reads_ctdsub_nuc))

        ct_labels.append(celltype)
        meand_celltype.append(np.nanmean(tdistanceall) if len(tdistanceall) > 0 else np.nan)
        nuclimall.append(np.nanmean(nuclim) if len(nuclim) > 0 else np.nan)
        cellmall.append(np.nanmean(cellm) if len(cellm) > 0 else np.nan)

    # --- Per-celltype summary ---
    per_celltype = pd.DataFrame({
        'celltype': ct_labels,
        'turnover': meand_celltype,
        'nuclei_size': nuclimall,
        'cell_size': cellmall
    })

    # --- 5-8d. Optimal expansion = turnover − nuclei_size ---
    # Use nanmean (matching notebook cell 32: np.nanmean)
    valid_turnover = [t for t in meand_celltype if not np.isnan(t)]
    valid_nuc = [n for n in nuclimall if not np.isnan(n)]
    n_valid_ct = len(valid_turnover)
    n_total_ct = len(meand_celltype)

    if n_valid_ct == 0:
        logger.warning(
            f"No cell types had a valid crossover out of {n_total_ct} total. "
            f"Correlation range likely too compressed for this gene panel. "
            f"Consider lowering diff_threshold (currently {diff_threshold}) "
            f"or changing threshold_mode (currently '{threshold_mode}')."
        )
        mean_turnover = np.nan
        mean_nuc_size = np.nanmean(nuclimall) if valid_nuc else np.nan
        optimal_expansion = 0.0
    else:
        mean_turnover = np.nanmean(meand_celltype)
        mean_nuc_size = np.nanmean(nuclimall)
        optimal_expansion = mean_turnover - mean_nuc_size

    print(f"\n    [Turnover] Valid cell types       : {n_valid_ct}/{n_total_ct}")
    print(f"    [Turnover] Mean turnover dist     : {mean_turnover}")
    print(f"    [Turnover] Mean nuclei size       : {mean_nuc_size}")
    print(f"    [Turnover] >>> Optimal expansion   : {optimal_expansion}")

    if optimal_expansion < 0:
        logger.warning(
            f"Optimal expansion is negative ({optimal_expansion:.3f}). "
            f"Turnover ({mean_turnover:.3f}) < nuclei_size ({mean_nuc_size:.3f}). "
            f"Clamping to 0."
        )
        optimal_expansion = 0.0

    # --- 5-8e. Summary barplot + CSVs (matching notebook cells 33-39) ---
    try:
        per_ct_sorted = per_celltype.dropna(subset=['turnover']).sort_values('turnover')
        if len(per_ct_sorted) > 0:
            # Build per-domain score table for error-bar barplot (notebook tfmerge style)
            tball_frames = []
            for ct in turnover_summ.columns:
                tb = pd.DataFrame(turnover_summ[ct]).copy()
                tb.columns = ['score']
                tb['cluster'] = ct
                tb = tb.reset_index(drop=True)
                tball_frames.append(tb)
            if tball_frames:
                tball = pd.concat(tball_frames, axis=0, ignore_index=True)
            else:
                tball = pd.DataFrame(columns=['score', 'cluster'])

            # Merge per-domain scores with per-celltype summary
            tf = per_ct_sorted.copy()
            tf['cluster'] = tf['celltype']
            tfmerge = tball.merge(tf, on='cluster', how='inner')
            tfmerge = tfmerge.sort_values(by='turnover')

            # Determine palette: use custom colors if available, else tab20
            if celltype_colors and len(celltype_colors) > 0:
                palette = [celltype_colors.get(ct, '#999999') for ct in tf.sort_values('turnover')['celltype']]
            else:
                palette = 'tab20'

            # --- Main barplot (all cell types) ---
            fig, ax = plt.subplots(figsize=(10, max(4, len(per_ct_sorted) * 0.5)))
            sns.scatterplot(
                data=tfmerge, y='cluster', x='cell_size',
                edgecolor='gray', color='black', s=80, label='cell_size', ax=ax, zorder=3
            )
            sns.barplot(
                data=tfmerge, y='cluster', x='score',
                palette=palette, alpha=0.9, errorbar='sd', ax=ax
            )
            sns.scatterplot(
                data=tfmerge, y='cluster', x='nuclei_size',
                color='#D83066', edgecolor=None, alpha=0.7, s=80, label='nuclei_size', ax=ax, zorder=3
            )
            ax.set_xlabel("Distance")
            ax.set_ylabel("Cell type")
            ax.set_title("Turnover per cell type (per-domain scores)")
            ax.legend(loc='lower right')
            summary_plot_path = os.path.join(output_dir, f"{sample_tag}_step5_turnover_barplot.png")
            fig.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"    [Turnover] Saved barplot: {summary_plot_path}")

            # --- Filtered barplot: cell types with >5 domains (notebook cells 38-39) ---
            domain_counts_per_ct = tball.dropna(subset=['score']).groupby('cluster').size()
            ct_with_many_domains = domain_counts_per_ct[domain_counts_per_ct > 5].index
            if len(ct_with_many_domains) > 1:
                tfmerge_sub = tfmerge[tfmerge['cluster'].isin(ct_with_many_domains)]
                tf_sub = tf[tf['celltype'].isin(ct_with_many_domains)]
                if celltype_colors and len(celltype_colors) > 0:
                    palette_sub = [celltype_colors.get(ct, '#999999')
                                   for ct in tf_sub.sort_values('turnover')['celltype']]
                else:
                    palette_sub = 'tab20'

                fig2, ax2 = plt.subplots(figsize=(10, max(4, len(tf_sub) * 0.5)))
                sns.scatterplot(
                    data=tfmerge_sub, y='cluster', x='cell_size',
                    edgecolor='gray', color='black', s=80, label='cell_size', ax=ax2, zorder=3
                )
                sns.barplot(
                    data=tfmerge_sub, y='cluster', x='score',
                    palette=palette_sub, alpha=0.9, errorbar='sd', ax=ax2
                )
                sns.scatterplot(
                    data=tfmerge_sub, y='cluster', x='nuclei_size',
                    color='#D83066', edgecolor=None, alpha=0.7, s=80, label='nuclei_size', ax=ax2, zorder=3
                )
                ax2.set_xlabel("Distance")
                ax2.set_ylabel("Cell type")
                ax2.set_title("Turnover (cell types with >5 domains)")
                ax2.legend(loc='lower right')
                filt_plot_path = os.path.join(output_dir, f"{sample_tag}_step5_turnover_barplot_filtered.png")
                fig2.savefig(filt_plot_path, dpi=300, bbox_inches='tight')
                plt.close(fig2)
                print(f"    [Turnover] Saved filtered barplot: {filt_plot_path}")
    except Exception as e:
        logger.warning(f"    [Turnover] Summary barplot failed: {e}")
        plt.close('all')

    # --- Save CSVs ---
    csv_summ = os.path.join(output_dir, f"{sample_tag}_step5_turnover_summary.csv")
    turnover_summ.to_csv(csv_summ)
    print(f"    [Turnover] Saved turnover summary matrix: {csv_summ}")

    csv_ct = os.path.join(output_dir, f"{sample_tag}_step5_turnover_per_celltype.csv")
    per_celltype.to_csv(csv_ct, index=False)
    print(f"    [Turnover] Saved per-celltype table:      {csv_ct}")

    txt_opt = os.path.join(output_dir, f"{sample_tag}_step5_optimal_expansion.txt")
    with open(txt_opt, 'w') as fh:
        fh.write(f"optimal_expansion\t{optimal_expansion}\n")
        fh.write(f"median_turnover\t{mean_turnover}\n")
        fh.write(f"median_nuclei_size\t{mean_nuc_size}\n")
        fh.write(f"valid_celltypes\t{n_valid_ct}/{n_total_ct}\n")
        fh.write(f"overlaps_nucleus_proxy\t{used_proxy}\n")
        if used_proxy:
            fh.write(f"proxy_threshold_median_distance\t{proxy_threshold}\n")
    print(f"    [Turnover] Saved optimal expansion value:  {txt_opt}")

    return turnover_summ, per_celltype, optimal_expansion


# --- Label Transfer from scRNA-seq Reference ---

def _transfer_celltype_labels(adata, sc_ref_path, ref_celltype_key='subclass_label',
                               n_neighbors=15):
    """Transfer cell type labels from scRNA-seq reference via kNN in shared gene PCA space.

    Adds 'celltype' and 'celltype_confidence' columns to adata.obs.
    Returns the column name on success, None on failure.
    """
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import issparse

    print(f"\n    [Label Transfer] Loading scRNA reference: {sc_ref_path}")
    adata_ref = sc.read_h5ad(sc_ref_path)

    if ref_celltype_key not in adata_ref.obs.columns:
        logger.warning("    Reference missing '%s' column. Available: %s",
                        ref_celltype_key, list(adata_ref.obs.columns[:10]))
        del adata_ref
        return None

    # Shared genes
    shared_genes = sorted(adata.var_names.intersection(adata_ref.var_names))
    print(f"    [Label Transfer] {len(shared_genes)} shared genes "
          f"(Xenium {adata.n_vars}, reference {adata_ref.n_vars})")

    if len(shared_genes) < 20:
        logger.warning("    Too few shared genes (%d) for reliable transfer.", len(shared_genes))
        del adata_ref
        return None

    # Prepare reference subset
    ref_sub = adata_ref[:, shared_genes].copy()
    valid_mask = ref_sub.obs[ref_celltype_key].notna()
    if hasattr(ref_sub.obs[ref_celltype_key], 'str'):
        valid_mask = valid_mask & (ref_sub.obs[ref_celltype_key].astype(str) != 'nan')
    ref_sub = ref_sub[valid_mask].copy()
    n_types = ref_sub.obs[ref_celltype_key].nunique()
    print(f"    [Label Transfer] Reference: {ref_sub.n_obs} cells, {n_types} cell types")

    if 'raw' in ref_sub.layers:
        ref_sub.X = ref_sub.layers['raw'].copy()
    if issparse(ref_sub.X):
        ref_sub.X = np.asarray(ref_sub.X.todense())
    sc.pp.normalize_total(ref_sub, target_sum=1e4)
    sc.pp.log1p(ref_sub)
    sc.pp.scale(ref_sub, max_value=10)
    sc.pp.pca(ref_sub)

    # Project target into reference PCA space
    tgt_sub = adata[:, shared_genes].copy()
    if 'raw' in tgt_sub.layers:
        tgt_sub.X = tgt_sub.layers['raw'].copy()
    if issparse(tgt_sub.X):
        tgt_sub.X = np.asarray(tgt_sub.X.todense())
    sc.pp.normalize_total(tgt_sub, target_sum=1e4)
    sc.pp.log1p(tgt_sub)

    tgt_X = np.array(tgt_sub.X)
    ref_mean = ref_sub.var['mean'].values
    ref_std = ref_sub.var['std'].values.copy()
    ref_std[ref_std == 0] = 1.0
    tgt_X = np.clip((tgt_X - ref_mean) / ref_std, -10, 10)
    tgt_pca = tgt_X @ ref_sub.varm['PCs']

    # kNN in aligned PCA space
    print(f"    [Label Transfer] Running kNN (k={n_neighbors}) in PCA space...")
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(ref_sub.obsm['X_pca'])
    distances, indices = nn.kneighbors(tgt_pca)

    ref_labels = ref_sub.obs[ref_celltype_key].values
    cell_labels = []
    confidence_scores = []
    for idx_row in indices:
        neighbour_labels = ref_labels[idx_row]
        values, counts = np.unique(neighbour_labels, return_counts=True)
        winner_idx = np.argmax(counts)
        cell_labels.append(values[winner_idx])
        confidence_scores.append(counts[winner_idx] / len(neighbour_labels))

    adata.obs['celltype'] = cell_labels
    adata.obs['celltype_confidence'] = confidence_scores

    assigned_types = len(np.unique(cell_labels))
    mean_conf = np.mean(confidence_scores)
    print(f"    [Label Transfer] Done: {assigned_types} cell types, mean confidence {mean_conf:.2f}")
    print(f"    [Label Transfer] Cell types: {sorted(pd.Series(cell_labels).unique())}")

    del adata_ref, ref_sub, tgt_sub
    return 'celltype'


# --- Main entry point ---

def run_step5(config):
    """Step 5 entry point: domain expansion + turnover analysis."""
    print("\n" + "=" * 60)
    print("[Step 5] Optimal Expansion (Ref: 4_1 Notebook)")
    print("=" * 60)

    output_dir = config["output_dir"]
    sample_tag = config["sample_tag"]

    exp_config = config.get("optimal_expansion", {})
    run_exp = exp_config.get("run_expansion", True)
    subsample_frac = exp_config.get("subsample_fraction", 0.01)
    dist_threshold = exp_config.get("distance_threshold", None)
    min_reads_per_domain = exp_config.get("min_reads_per_domain", 5000)
    diff_threshold = exp_config.get("diff_threshold", 0.1)
    threshold_mode = exp_config.get("threshold_mode", "absolute")

    if not run_exp:
        print("    - [Info] 'run_expansion' is False in config. Skipping Step 5.")
        return

    # --- 5-1. Load original transcripts (Step 0) ---
    print("\n[Step 5-1] Loading Data...")

    step0_file = os.path.join(output_dir, f"{sample_tag}.h5ad")
    if not os.path.exists(step0_file):
        parent_dir = os.path.dirname(output_dir)
        step0_file_alt = os.path.join(parent_dir, "step0_formatting", f"{sample_tag}.h5ad")
        if os.path.exists(step0_file_alt):
            step0_file = step0_file_alt

    if not os.path.exists(step0_file):
        logging.error(f"    - Step 0 output not found at {step0_file} or in sibling dir. Cannot run expansion.")
        return

    print(f"    - Loading Step 0 Data from: {step0_file}")
    adata_step0 = sc.read_h5ad(step0_file)

    reads_original = None
    # 1) Legacy: uns['spots'] DataFrame (older Step 0 outputs)
    if 'spots' in adata_step0.uns:
        reads_original = adata_step0.uns['spots'].copy()
    else:
        # 2) uns['spots_path'] — Step 0 stores parquet path here after del uns['spots']
        spots_path = adata_step0.uns.get('spots_path')
        if spots_path and os.path.exists(str(spots_path)):
            print(f"    - Loading transcripts from sidecar: {spots_path}")
            reads_original = pd.read_parquet(str(spots_path))
        else:
            # 3) Convention-based sidecar search
            step0_dir = os.path.dirname(step0_file)
            for candidate in [
                os.path.join(step0_dir, f"{sample_tag}_transcripts.parquet"),
                os.path.join(step0_dir, "transcripts.parquet"),
                os.path.join(step0_dir, "transcripts.csv"),
            ]:
                if os.path.exists(candidate):
                    print(f"    - Loading transcripts from sidecar: {candidate}")
                    if candidate.endswith('.parquet'):
                        reads_original = pd.read_parquet(candidate)
                    else:
                        reads_original = pd.read_csv(candidate, low_memory=False)
                    break

    # Decode bytes columns from parquet (Xenium parquet stores strings as bytes)
    if reads_original is not None:
        for col in reads_original.columns:
            if (reads_original[col].dtype == object
                    and len(reads_original) > 0
                    and isinstance(reads_original[col].iloc[0], bytes)):
                reads_original[col] = reads_original[col].str.decode('utf-8')

    if reads_original is None:
        logging.error(
            "    - Transcripts not found: uns['spots'], uns['spots_path'], "
            "or sidecar files all missing."
        )
        return

    print(f"    - Loaded {len(reads_original)} original reads.")

    # --- 5-2. Load annotated cells (Step 1/2/4) ---
    # Load annotated cells (from Step 2, 1, or 4)
    input_adata_file = None
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        input_adata_file = config['previous_step_adata_path']

    if not input_adata_file or not os.path.exists(input_adata_file):
        print("    - [Info] 'previous_step_adata_path' missing or invalid. Searching...")
        parent_dir = os.path.dirname(output_dir)
        possible_paths = [
            os.path.join(parent_dir, "step3_resegmentation", f"{sample_tag}_step3_resegmented.h5ad"),
            os.path.join(parent_dir, "step1_exploration", f"{sample_tag}_step1_exploration.h5ad"),
            os.path.join(parent_dir, "step2_segmentation_free", f"{sample_tag}_step2_points2regions.h5ad")
        ]
        for p in possible_paths:
            if os.path.exists(p):
                input_adata_file = p
                break

    print(f"    - Loading Annotated Data from: {input_adata_file}")
    if not input_adata_file or not os.path.exists(input_adata_file):
        logging.error("    - Annotated adata not found. Run Step 1, 2, or 4 first.")
        return

    adata_annotated = sc.read_h5ad(input_adata_file)
    print(f"    - Loaded annotated cells: {adata_annotated.shape}")

    # --- Label Transfer: ensure cell-type annotation exists ---
    celltype_check = ['Class', 'celltype', 'cell_type', 'celltype_majority', 'initial_annotation']
    has_celltype = any(k in adata_annotated.obs.columns for k in celltype_check
                       if k not in ['leiden', 'cluster', 'graph_clusters'])
    if not has_celltype:
        sc_ref_path = (config.get("comparison", {}).get("sc_reference_path")
                       or config.get("benchmark", {}).get("reference_adata")
                       or config.get("benchmark", {}).get("sc_reference_path"))
        ref_ct_key = config.get("sc_reference", {}).get("celltype_key", "subclass_label")

        if sc_ref_path and os.path.exists(str(sc_ref_path)):
            ct_col = _transfer_celltype_labels(
                adata_annotated, str(sc_ref_path), ref_celltype_key=ref_ct_key)
            if ct_col:
                # Save annotated adata back so downstream steps can reuse it
                print(f"    - Saving annotated adata back to: {input_adata_file}")
                adata_annotated.write_h5ad(input_adata_file)
        else:
            print("    - [WARNING] No cell-type annotation and no scRNA reference available.")
            print("      Turnover analysis will use domain as cell-type proxy (may give trivial results).")

    # --- 5-3. Map domain assignments to reads ---
    print("\n[Step 5-2] Identifying Domains & Unassigned Reads...")

    # Domain key (spatial regions)
    domain_key = None

    # Priority 1: True spatial annotations
    for key in ['spatial_annotation', 'region_annotation']:
        if key in adata_annotated.obs.columns:
            domain_key = key
            break

    # Priority 2: P2R spatial domains from Step 2
    if domain_key is None:
        p2r_domains = _load_p2r_spatial_domains(
            adata_annotated, output_dir, sample_tag,
            n_clusters=exp_config.get('p2r_n_clusters', 50),
            max_dist_um=exp_config.get('p2r_max_dist_um', 100.0))
        if p2r_domains is not None:
            adata_annotated.obs['p2r_domain'] = p2r_domains
            domain_key = 'p2r_domain'

    # Priority 3: Expression-based fallbacks (less ideal but functional)
    if domain_key is None:
        for key in ['leiden', 'cluster', 'graph_clusters']:
            if key in adata_annotated.obs.columns:
                domain_key = key
                break

    # Priority 4: Run leiden as last resort
    if not domain_key:
        print("    - No domain/cluster key found. Running leiden clustering as fallback...")
        try:
            adata_tmp = adata_annotated.copy()
            sc.pp.normalize_total(adata_tmp, target_sum=100)
            sc.pp.log1p(adata_tmp)
            sc.pp.pca(adata_tmp)
            sc.pp.neighbors(adata_tmp, n_neighbors=15)
            sc.tl.leiden(adata_tmp, resolution=1.0, key_added='leiden')
            adata_annotated.obs['leiden'] = adata_tmp.obs['leiden']
            domain_key = 'leiden'
            print(f"    - Fallback leiden clustering complete: {adata_annotated.obs['leiden'].nunique()} clusters")
            del adata_tmp
        except Exception as e:
            logging.error(f"    - Fallback leiden clustering failed: {e}")
            return
    print(f"    - Using '{domain_key}' as domain source.")

    # Cell-type key (distinct from domain — used for per-celltype turnover)
    celltype_key = None
    celltype_priority = ['Class', 'celltype', 'cell_type', 'celltype_majority', 'initial_annotation', 'ct_majority']
    for key in celltype_priority:
        if key in adata_annotated.obs.columns:
            celltype_key = key
            break
    if celltype_key is None:
        celltype_key = domain_key  # fallback: use domain as cell-type proxy
        logger.warning("    - No cell-type column found; using domain key '%s' as cell-type proxy.", domain_key)
    print(f"    - Using '{celltype_key}' as cell-type annotation source.")

    if 'cell_id' in adata_annotated.obs.columns:
        annotated_ids = adata_annotated.obs['cell_id']
    else:
        annotated_ids = adata_annotated.obs.index

    domain_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
    ct_map = dict(zip(annotated_ids, adata_annotated.obs[celltype_key]))

    print("    - Mapping existing domains to reads...")
    if 'cell_id' not in reads_original.columns:
        if reads_original.index.name == 'cell_id':
            reads_original = reads_original.reset_index()
        else:
            logging.error("    - 'cell_id' missing in reads. Cannot link to cells.")
            return

    # --- Primary: cell_id-based mapping ---
    reads_original['domain'] = reads_original['cell_id'].map(domain_map)
    reads_original['initial_annotation'] = reads_original['cell_id'].map(ct_map)

    # --- Fallback: bridge via Step 3 reseg transcripts if cell_id spaces don't match
    #     (e.g. Cellpose mask labels vs original Xenium barcodes) ---
    match_rate = reads_original['domain'].notna().mean()
    if match_rate < 0.01:
        print(f"    - [WARNING] Cell ID match rate {match_rate:.1%} — ID spaces differ.")
        step_dir = os.path.dirname(input_adata_file) if input_adata_file else None
        if step_dir:
            step3_tx = os.path.join(
                step_dir, f"{sample_tag}_step3_transcripts_resegmented.csv")
            if os.path.exists(step3_tx):
                print(f"    - Bridging via reseg transcripts: {step3_tx}")
                df_bridge = pd.read_csv(step3_tx, low_memory=False)
                reseg_col = next(
                    (c for c in ['closest_cell', 'cell_id_reseg', 'in_cell']
                     if c in df_bridge.columns), None)
                if reseg_col and 'x_location' in df_bridge.columns:
                    df_bridge['_dom'] = df_bridge[reseg_col].astype(str).map(domain_map)
                    df_bridge['_ann'] = df_bridge[reseg_col].astype(str).map(ct_map)
                    # Coordinate-based lookup (string key for fast matching)
                    bkey = (df_bridge['x_location'].round(2).astype(str) + '_'
                            + df_bridge['y_location'].round(2).astype(str))
                    dom_lookup = dict(zip(bkey, df_bridge['_dom']))
                    ann_lookup = dict(zip(bkey, df_bridge['_ann']))
                    rkey = (reads_original['x_location'].round(2).astype(str) + '_'
                            + reads_original['y_location'].round(2).astype(str))
                    reads_original['domain'] = rkey.map(dom_lookup)
                    reads_original['initial_annotation'] = rkey.map(ann_lookup)
                    new_rate = reads_original['domain'].notna().mean()
                    print(f"    - After bridging: {new_rate:.1%} reads matched to domains")
                    del df_bridge

        if reads_original['domain'].notna().mean() < 0.01:
            logging.error(
                "    - Domain mapping failed. Annotated cell IDs don't match "
                "original transcript cell IDs, and no Step 3 transcripts CSV found. "
                "Consider using Step 1 data instead of Step 3 for Step 5."
            )
            return

    nancells = reads_original[reads_original['domain'].isna()]
    annotatedcells = reads_original[~reads_original['domain'].isna()]

    n_assigned = len(annotatedcells)
    n_unassigned = len(nancells)
    print(f"    - Assigned Reads (with domain): {n_assigned}")
    print(f"    - Unassigned Reads (to be expanded): {n_unassigned}")

    if n_unassigned == 0:
        print("    - No unassigned reads found. Expansion not needed.")
        return

    # --- 5-4. KDTree nearest-domain expansion ---
    print(f"\n[Step 5-3] Building cKDTree (Subsample Fraction: {subsample_frac})...")

    n_sample = int(np.round(n_assigned * subsample_frac))
    if n_sample < 100:
        n_sample = min(n_assigned, 1000)

    print(f"    - Subsampling {n_sample} reads as anchors...")
    annotated_sub = annotatedcells.sample(n=n_sample, random_state=42)

    coords1 = annotated_sub[['x_location', 'y_location']].values

    print("    - Constructing KDTree...")
    tree = cKDTree(coords1)

    # Query unassigned reads against KDTree
    print(f"\n[Step 5-4] Assigning {n_unassigned} reads to nearest domains...")

    coords2 = nancells[['x_location', 'y_location']].values

    chunk_size = 1000000
    n_chunks = int(np.ceil(len(coords2) / chunk_size))

    all_indices = []
    all_distances = []

    print(f"    - Processing in {n_chunks} chunks...")
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, len(coords2))
        chunk = coords2[start:end]
        dists, idxs = tree.query(chunk, k=1)
        all_indices.append(idxs)
        all_distances.append(dists)

    concat_indices = np.concatenate(all_indices)
    concat_distances = np.concatenate(all_distances)

    # --- 5-5. Distance threshold filter (optional) ---
    print("\n[Step 5-5] Finalizing Assignment...")

    anchor_domains = annotated_sub['domain'].values
    assigned_domains = anchor_domains[concat_indices]

    if dist_threshold:
        print(f"    - Applying distance threshold: {dist_threshold}")
        mask_too_far = concat_distances > dist_threshold
        params_assigned = assigned_domains.copy()
        params_assigned[mask_too_far] = np.nan
        assigned_domains = params_assigned
        n_filtered = np.sum(mask_too_far)
        print(f"    - {n_filtered} reads exceeded distance threshold.")

    anchor_annotations = annotated_sub['initial_annotation'].values
    assigned_annotations = anchor_annotations[concat_indices]

    nancells = nancells.copy()
    nancells.loc[:, 'domain'] = assigned_domains
    nancells.loc[:, 'initial_annotation'] = assigned_annotations
    reads_original.loc[nancells.index, 'domain'] = assigned_domains
    reads_original.loc[nancells.index, 'initial_annotation'] = assigned_annotations

    print("    - Merging results complete.")

    # --- 5-6. Spatial map + distance histogram ---
    print("\n[Step 5-6] Generating Visualization...")
    plot_df = reads_original.sample(n=min(len(reads_original), 100000), random_state=42)

    fig, ax = plt.subplots(figsize=(12, 10))
    plot_df_clean = plot_df.dropna(subset=['domain'])
    n_domains = plot_df_clean['domain'].nunique()
    sns.scatterplot(
        data=plot_df_clean,
        x='x_location',
        y='y_location',
        hue='domain',
        s=1,
        linewidth=0,
        palette='tab20',
        legend='brief' if n_domains <= 25 else False,
        ax=ax
    )
    ax.set_title(f"Optimal Expansion Result (Subsample {subsample_frac})")
    ax.set_aspect('equal')
    if n_domains <= 25:
        ax.legend(markerscale=5, fontsize=6, loc='center left', bbox_to_anchor=(1, 0.5),
                  title='Domain', title_fontsize=7)

    plot_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_map.png")
    fig.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    - Saved Expansion Map to: {plot_file}")

    # Xenium x_location/y_location are already in µm, so KDTree distances
    # are natively in µm — no pixel conversion needed.
    distances_um = concat_distances

    fig_dist, ax_dist = plt.subplots(figsize=(6, 4))
    counts, bin_edges, patches = ax_dist.hist(distances_um, bins=50, color='orange', alpha=0.7, edgecolor='white')
    sns.kdeplot(distances_um, color='darkorange', linewidth=1.5, ax=ax_dist)

    # Annotate peak bin
    peak_idx = int(np.argmax(counts))
    peak_count = int(counts[peak_idx])
    peak_dist = (bin_edges[peak_idx] + bin_edges[peak_idx + 1]) / 2
    ax_dist.annotate(
        f'Peak: {peak_dist:.1f} µm, n={peak_count:,}',
        xy=(peak_dist, peak_count),
        xytext=(peak_dist + (distances_um.max() - distances_um.min()) * 0.15, peak_count * 0.9),
        arrowprops=dict(arrowstyle='->', color='black', lw=1.2),
        fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9)
    )

    ax_dist.set_title("Distance to Nearest Domain-Anchor")
    ax_dist.set_xlabel("Distance (µm)")
    ax_dist.set_ylabel("Count")
    dist_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_distances.png")
    fig_dist.savefig(dist_file, dpi=300, bbox_inches='tight')
    plt.close(fig_dist)
    print(f"    - Saved Distance Histogram to: {dist_file}")

    # --- 5-7. Save expanded transcripts ---
    output_file = os.path.join(output_dir, f"{sample_tag}_step5_expanded_transcripts.csv")
    print(f"    - Saving Expanded Transcripts to: {output_file}")
    reads_original.to_csv(output_file, index=False)

    # --- 5-8. Turnover analysis ---
    print("\n[Step 5-8] Running Turnover/Crossover Analysis...")

    # Separate reads by cell assignment:
    #   reads_not_assigned = cell_id == -1 or 'UNASSIGNED' (background)
    #   reads_assigned = valid cell_id (segmented cells)
    print("    - Preparing reads for turnover analysis...")

    if reads_original['cell_id'].dtype in [np.int64, np.int32, np.float64, np.float32, int, float]:
        mask_not_assigned = reads_original['cell_id'] == -1
    else:
        mask_not_assigned = (
            reads_original['cell_id'].astype(str).str.upper().isin(['-1', 'UNASSIGNED', 'NAN', ''])
        )

    reads_not_assigned = reads_original[mask_not_assigned].copy()
    reads_assigned_all = reads_original[~mask_not_assigned].copy()

    # Filter out negative control probes (NegControlProbe, NegControlCodeword, BLANK, antisense)
    # These have no cell-type specificity and add noise to correlation-based turnover analysis.
    if 'feature_name' in reads_not_assigned.columns:
        _ctrl = reads_not_assigned['feature_name'].str.contains(
            'BLANK|NegControl|antisense', case=False, na=False)
        n_ctrl_bg = _ctrl.sum()
        if n_ctrl_bg > 0:
            reads_not_assigned = reads_not_assigned[~_ctrl]
            print(f"    - Filtered {n_ctrl_bg} negative control transcripts from background reads")
    if 'feature_name' in reads_assigned_all.columns:
        _ctrl = reads_assigned_all['feature_name'].str.contains(
            'BLANK|NegControl|antisense', case=False, na=False)
        n_ctrl_fg = _ctrl.sum()
        if n_ctrl_fg > 0:
            reads_assigned_all = reads_assigned_all[~_ctrl]
            print(f"    - Filtered {n_ctrl_fg} negative control transcripts from assigned reads")

    print(f"    - Reads assigned to cells:   {len(reads_assigned_all)}")
    print(f"    - Reads not assigned (background): {len(reads_not_assigned)}")

    # Compute distance from each assigned read to its cell centroid
    # Try each centroid source and verify cell_id match rate before committing
    cx_map, cy_map = None, None

    # 1) Try annotated adata (e.g. Step 3) centroids
    if 'x_centroid' in adata_annotated.obs.columns and 'y_centroid' in adata_annotated.obs.columns:
        _cx = dict(zip(annotated_ids, adata_annotated.obs['x_centroid']))
        _cy = dict(zip(annotated_ids, adata_annotated.obs['y_centroid']))
        _match = reads_assigned_all['cell_id'].isin(_cx.keys()).mean()
        if _match > 0.01:
            cx_map, cy_map = _cx, _cy
            print(f"    - Using annotated adata centroids (match rate={_match:.1%})")

    # 2) Fall back to Step 0 centroids if annotated IDs don't match reads
    if cx_map is None and 'x_centroid' in adata_step0.obs.columns and 'y_centroid' in adata_step0.obs.columns:
        if 'cell_id' in adata_step0.obs.columns:
            _cx = dict(zip(adata_step0.obs['cell_id'], adata_step0.obs['x_centroid']))
            _cy = dict(zip(adata_step0.obs['cell_id'], adata_step0.obs['y_centroid']))
        else:
            _cx = dict(zip(adata_step0.obs.index, adata_step0.obs['x_centroid']))
            _cy = dict(zip(adata_step0.obs.index, adata_step0.obs['y_centroid']))
        _match = reads_assigned_all['cell_id'].isin(_cx.keys()).mean()
        if _match > 0.01:
            cx_map, cy_map = _cx, _cy
            print(f"    - Using Step 0 centroids (match rate={_match:.1%})")

    if cx_map is None:
        print("    - [Warning] No centroids match read cell_ids. Skipping Turnover Analysis.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    reads_assigned_val = reads_assigned_all[reads_assigned_all['cell_id'].isin(cx_map.keys())].copy()
    print(f"    - Reads with centroid info: {len(reads_assigned_val)}")

    if len(reads_assigned_val) == 0:
        print("    - [Warning] No reads matched to cells with centroids. Skipping Turnover.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    reads_assigned_val['x_cell'] = reads_assigned_val['cell_id'].map(cx_map).astype(float)
    reads_assigned_val['y_cell'] = reads_assigned_val['cell_id'].map(cy_map).astype(float)
    reads_assigned_val['distance'] = np.sqrt(
        (reads_assigned_val['x_location'] - reads_assigned_val['x_cell']) ** 2 +
        (reads_assigned_val['y_location'] - reads_assigned_val['y_cell']) ** 2
    )
    reads_assigned_val['distance'] = reads_assigned_val['distance'].round(0)

    # --- QC scatter: reads + centroids overlay (notebook cell 22) ---
    try:
        fig_qc, ax_qc = plt.subplots(figsize=(10, 10))
        sub_reads = reads_assigned_val.sample(n=min(len(reads_assigned_val), 100000), random_state=42)
        ax_qc.scatter(sub_reads['x_location'], sub_reads['y_location'], s=1, alpha=0.3, label='reads')
        centroid_df = sub_reads[['x_cell', 'y_cell']].drop_duplicates()
        ax_qc.scatter(centroid_df['x_cell'], centroid_df['y_cell'],
                       s=0.5, color='red', alpha=0.5, label='centroids')
        ax_qc.set_title("Reads vs Cell Centroids")
        ax_qc.legend(markerscale=5)
        ax_qc.axis('equal')
        qc_path = os.path.join(output_dir, f"{sample_tag}_step5_reads_vs_centroids.png")
        fig_qc.savefig(qc_path, dpi=150, bbox_inches='tight')
        plt.close(fig_qc)
        print(f"    - Saved reads vs centroids QC plot: {qc_path}")
    except Exception as e:
        logger.warning(f"    - Reads vs centroids QC plot failed: {e}")
        plt.close('all')

    if 'initial_annotation' not in reads_assigned_val.columns or reads_assigned_val['initial_annotation'].isna().all():
        reads_assigned_val['initial_annotation'] = reads_assigned_val['cell_id'].map(
            dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
        )

    if 'feature_name' not in reads_assigned_val.columns:
        print("    - [Warning] 'feature_name' missing in reads. Skipping turnover.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    # Extract custom colors from adata if available (notebook uses Class_colors)
    celltype_colors = None
    if celltype_key:
        color_key = f'{celltype_key}_colors'
        if color_key in adata_annotated.uns:
            cats = adata_annotated.obs[celltype_key].cat.categories if hasattr(
                adata_annotated.obs[celltype_key], 'cat') else adata_annotated.obs[celltype_key].unique()
            colors = adata_annotated.uns[color_key]
            if len(colors) >= len(cats):
                celltype_colors = dict(zip(cats, colors[:len(cats)]))
                print(f"    - Using custom '{color_key}' palette ({len(celltype_colors)} colors)")

    turnover_csv = os.path.join(output_dir, f"{sample_tag}_step5_turnover_per_celltype.csv")
    turnover_summ_csv = os.path.join(output_dir, f"{sample_tag}_step5_turnover_summary.csv")

    if os.path.exists(turnover_csv) and os.path.exists(turnover_summ_csv):
        print(f"\n    [CACHE HIT] Turnover results found. Skipping recomputation.")
        print(f"      - {turnover_csv}")
        print(f"      - {turnover_summ_csv}")
        per_celltype = pd.read_csv(turnover_csv)
    else:
        turnover_summ, per_celltype, optimal_expansion = calculate_turnover(
            reads_assigned=reads_assigned_val,
            reads_not_assigned=reads_not_assigned,
            output_dir=output_dir,
            sample_tag=sample_tag,
            min_reads_per_domain=min_reads_per_domain,
            diff_threshold=diff_threshold,
            celltype_colors=celltype_colors,
            threshold_mode=threshold_mode
        )

    if per_celltype is not None:
        print(f"\n    Turnover per-celltype summary:")
        print(per_celltype.to_string(index=False))

    print("\n=== Step 5 Optimal Expansion Complete ===")


if __name__ == "__main__":
    pass
