# Step 4: Techniques Comparison & Validation (Ref: notebooks/3_techniques_comparison/3_3-3_7)
# Compares original vs resegmented data using efficiency, specificity,
# positivity, and diffusion metrics.
#
# Flow:
#   4-1. Load resegmented + original data
#   4-2. Efficiency analysis
#     4-2a. Transcripts/genes per cell histograms
#     4-2b. Expression ratio (ST vs scRNAseq, CPM)
#     4-2c. Region-based efficiency breakdown
#   4-3. Specificity analysis
#     4-3a. Negative marker purity (coexpression)
#     4-3b. Gene-gene correlation heatmap
#   4-4. Positivity analysis
#     4-4a. Positivity histogram
#     4-4b. Preprocessing + Leiden clustering
#     4-4c. Violin plots per cluster
#     4-4d. UMAP per top gene
#   4-5. Diffusion analysis
#     4-5a. Transcript-to-centroid distance (px→µm)
#     4-5b. Complementary CDF plot
#     4-5c. Per-gene ECDF subplots
#     4-5d. Gene×Method distance heatmap

import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm


# --- Inlined from xb/calculating.py (pipeline self-contained) ---

def _coexpression_calculation(exp, min_exp=0):
    """Gene-gene co-expression matrix: fraction of positive cells per gene pair."""
    coexpression = pd.DataFrame(index=exp.columns, columns=exp.columns)
    for col in tqdm(exp.columns):
        sel = exp.loc[:, col] > min_exp
        positive_cells = exp.loc[sel, :]
        coexpression.loc[:, col] = np.sum(positive_cells > min_exp) / positive_cells.shape[0]
    coexpression = coexpression.fillna(1)
    return coexpression


def _negative_marker_purity_coexpression(adata_sp, adata_sc, key='celltype', pipeline_output=True, minexp=0.0):
    """Negative marker purity via co-expression (from xb.calculating)."""
    min_number_cells = 10
    minimum_exp = 0.05

    adata_sp = adata_sp[:, adata_sp.var_names.isin(adata_sc.var_names)]
    adata_sc = adata_sc[:, adata_sp.var_names]

    try:
        exp_sc = pd.DataFrame(adata_sc.X.todense(), columns=adata_sc.var_names)
    except Exception:
        exp_sc = pd.DataFrame(adata_sc.X, columns=adata_sc.var_names)
    try:
        exp_sp = pd.DataFrame(adata_sp.X.todense(), columns=adata_sp.var_names)
    except Exception:
        exp_sp = pd.DataFrame(adata_sp.X, columns=adata_sp.var_names)

    mean_celltype_sc = _coexpression_calculation(exp_sc, min_exp=minexp)
    mean_celltype_sp = _coexpression_calculation(exp_sp, min_exp=minexp)

    mean_ct_sc_rel = mean_celltype_sc
    mean_ct_sp_rel = mean_celltype_sp
    mean_ct_sc_norm = mean_celltype_sc
    mean_ct_sp_norm = mean_celltype_sp

    neg_marker_mask = np.array(mean_ct_sc_rel < minimum_exp)

    if np.sum(neg_marker_mask) < 1:
        print("No negative markers were found in the sc data reference.")
        negative_marker_purity = 'nan'
        if pipeline_output:
            return negative_marker_purity
        else:
            return negative_marker_purity, None, None

    lowvals_sc = mean_ct_sc_norm.values[neg_marker_mask]
    lowvals_sp = mean_ct_sp_norm.values[neg_marker_mask]

    lowvals_diff = (lowvals_sp - lowvals_sc)
    lowvals_diff[lowvals_diff < 0] = 0
    negative_marker_purity = 1 - np.mean(lowvals_diff)

    if pipeline_output:
        return negative_marker_purity
    else:
        purities = (mean_ct_sp_norm - mean_ct_sc_norm)
        purities[~neg_marker_mask] = np.nan
        purities = purities.loc[~(purities.isnull().all(axis=1)), ~(purities.isnull().all(axis=0))]
        purity_per_gene = purities.mean(axis=0, skipna=True)
        purity_per_celltype = purities.mean(axis=1, skipna=True)
        return negative_marker_purity, purity_per_gene, purity_per_celltype

# Pixel-to-um conversion factors by technology
PIXEL_TO_UM = {
    'xenium': 4.70588,
    'cosmx': 8.3333,
    'vizgen': 9.20586,
    'merfish': 9.28,
    'hybriss': 3.11,
    'resolvedbio': 7.24,
}


def run_step4(config, adata_path, transcripts_path, output_dir, original_adata_path=None, original_transcripts_path=None):
    """Step 4 entry point. Loads resegmented + original data, runs enabled analyses."""
    print("------------------------------------------------")
    print("Starting Step 4: Techniques Comparison & Validation")
    print("------------------------------------------------")

    # --- 4-1. Load resegmented + original data ---
    print(f"Loading Resegmented AnnData from {adata_path}...")
    try:
        adata_reseg = sc.read_h5ad(adata_path)
    except Exception as e:
        print(f"Error loading AnnData: {e}")
        return

    adata_orig = None
    if original_adata_path and os.path.exists(original_adata_path):
        print(f"Loading Original AnnData (for comparison) from {original_adata_path}...")
        try:
            adata_orig = sc.read_h5ad(original_adata_path)
        except Exception as e:
            print(f"Error loading Original AnnData: {e}")

    print(f"Loading Resegmented Transcripts from {transcripts_path}...")
    try:
        df_assigned = pd.read_csv(transcripts_path)
    except Exception as e:
        print(f"Error loading transcripts (needed for diffusion): {e}")
        df_assigned = None

    df_assigned_orig = None
    if original_transcripts_path and os.path.exists(original_transcripts_path):
        print(f"Loading Original Transcripts (for comparison) from {original_transcripts_path}...")
        try:
            df_assigned_orig = pd.read_csv(original_transcripts_path)
        except Exception as e:
            print(f"Error loading Original Transcripts: {e}")

    # --- 2. Config ---
    comp_config = config.get('comparison', {})
    run_efficiency = comp_config.get('run_efficiency', True)
    run_specificity = comp_config.get('run_specificity', True)
    run_positivity = comp_config.get('run_positivity', True)
    run_diffusion = comp_config.get('run_diffusion', True)

    figs_dir = os.path.join(output_dir, 'figures', '4_techniques_comparison')
    os.makedirs(figs_dir, exist_ok=True)

    sample_tag = config.get('sample_tag', 'sample')
    technology = comp_config.get('technology', 'xenium').lower()

    # --- 4-2. Efficiency analysis ---
    if run_efficiency:
        print("Running Efficiency Analysis...")
        try:
           analyze_efficiency(adata_reseg, output_dir, adata_orig, comp_config)
        except Exception as e:
            print(f"Error in Efficiency Analysis: {e}")

    # --- 4-3. Specificity analysis ---
    if run_specificity:
        print("Running Specificity Analysis (NMP Proxy)...")
        try:
            analyze_specificity(adata_reseg, comp_config, figs_dir, adata_orig)
        except Exception as e:
            print(f"Error in Specificity Analysis: {e}")

    # --- 4-4. Positivity analysis ---
    if run_positivity:
        print("Running Positivity Analysis (Validation)...")
        try:
            analyze_positivity(adata_reseg, figs_dir, sample_tag, adata_orig)
        except Exception as e:
            print(f"Error in Positivity Analysis: {e}")

    # --- 4-5. Diffusion analysis ---
    if run_diffusion and df_assigned is not None:
        print("Running Diffusion Analysis (Validation)...")
        try:
            analyze_diffusion(df_assigned, adata_reseg, figs_dir, df_assigned_orig, adata_orig, technology=technology)
        except Exception as e:
            print(f"Error in Diffusion Analysis: {e}")

    print("Step 4: Techniques Comparison & Validation Completed.")


def analyze_efficiency(adata_reseg, output_dir, adata_orig=None, comp_config=None):
    """Transcripts/genes per cell distributions + optional expression ratio (ST/scRNAseq)."""
    if comp_config is None:
        comp_config = {}

    def prep_qc(ad):
        if 'total_counts' not in ad.obs.columns or 'n_genes_by_counts' not in ad.obs.columns:
            sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)

    prep_qc(adata_reseg)

    data_list = []
    df_reseg = pd.DataFrame({
        'Transcripts per Cell': adata_reseg.obs['total_counts'],
        'Genes per Cell': adata_reseg.obs['n_genes_by_counts'],
        'Dataset': 'Resegmented (Step 3)'
    })
    data_list.append(df_reseg)

    if adata_orig is not None:
        prep_qc(adata_orig)
        df_orig = pd.DataFrame({
            'Transcripts per Cell': adata_orig.obs['total_counts'],
            'Genes per Cell': adata_orig.obs['n_genes_by_counts'],
            'Dataset': 'Original (Step 0/1)'
        })
        data_list.append(df_orig)

    df_plot = pd.concat(data_list)

    # --- 4-2a. Transcripts/genes per cell histograms ---
    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Transcripts per Cell', hue='Dataset', kde=True, bins=50, element="step")
    plt.title('Comparison: Transcripts per Cell')
    plt.xlabel('Transcripts per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_transcripts_per_cell_comparison.png'))
    plt.close()

    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Genes per Cell', hue='Dataset', kde=True, bins=50, element="step")
    plt.title('Comparison: Genes Detected per Cell')
    plt.xlabel('Genes per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_genes_per_cell_comparison.png'))
    plt.close()

    stats_list = []
    for label, ad in [('Resegmented', adata_reseg), ('Original', adata_orig)]:
        if ad is None: continue
        stats_list.append({
            'Dataset': label,
            'median_transcripts_per_cell': np.median(ad.obs['total_counts']),
            'mean_transcripts_per_cell': np.mean(ad.obs['total_counts']),
            'median_genes_per_cell': np.median(ad.obs['n_genes_by_counts']),
            'mean_genes_per_cell': np.mean(ad.obs['n_genes_by_counts'])
        })

    pd.DataFrame(stats_list).to_csv(os.path.join(output_dir, 'efficiency_metrics_comparison.csv'), index=False)
    print(f"  > Efficiency comparison saved to {output_dir}")

    # --- 4-2b. Expression ratio (ST vs scRNAseq, CPM) ---
    # ratio = ST_expression / scRNAseq_expression (both CPM-normalized)
    sc_ref_path = comp_config.get('sc_reference_path')
    if sc_ref_path and os.path.exists(sc_ref_path):
        print("  > scRNAseq reference found. Computing per-gene expression ratio (ST / scRNAseq)...")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)
            ad_st = adata_reseg.copy()
            sc.pp.normalize_total(ad_st, target_sum=1e6)
            ad_ref = adata_sc.copy()
            sc.pp.normalize_total(ad_ref, target_sum=1e6)

            common_genes = list(set(ad_st.var_names) & set(ad_ref.var_names))
            if len(common_genes) > 0:
                st_mean = np.array(ad_st[:, common_genes].X.mean(axis=0)).flatten()
                ref_mean = np.array(ad_ref[:, common_genes].X.mean(axis=0)).flatten()

                ref_mean_safe = np.where(ref_mean > 0, ref_mean, np.nan)
                ratio = st_mean / ref_mean_safe

                df_ratio = pd.DataFrame({
                    'gene': common_genes,
                    'st_mean_cpm': st_mean,
                    'sc_mean_cpm': ref_mean,
                    'expression_ratio': ratio,
                }).dropna().sort_values('expression_ratio', ascending=False)

                df_ratio.to_csv(os.path.join(output_dir, 'efficiency_expression_ratio.csv'), index=False)

                plt.figure(figsize=(8, 6))
                plt.hist(df_ratio['expression_ratio'].clip(upper=5), bins=50, edgecolor='black')
                plt.axvline(x=1.0, color='red', linestyle='--', label='Ratio = 1')
                plt.title('Per-Gene Expression Ratio (ST / scRNAseq)')
                plt.xlabel('Expression Ratio (clipped at 5)')
                plt.ylabel('Number of Genes')
                plt.legend()
                plt.savefig(os.path.join(output_dir, 'efficiency_expression_ratio.png'))
                plt.close()
                print(f"  > Expression ratio analysis saved ({len(df_ratio)} common genes).")
            else:
                print("  > No common genes found between ST and scRNAseq reference.")
        except Exception as e:
            print(f"  > Error computing expression ratio: {e}")
    else:
        print("  > NOTE: scRNAseq reference not available. Skipping expression ratio analysis.")
        print("    Set comparison.sc_reference_path in config to enable per-gene expression ratio.")

    # --- 4-2c. Region-based efficiency breakdown ---
    region_col = None
    for candidate in ['spatial_annotation', 'region', 'tissue_region', 'domain']:
        if candidate in adata_reseg.obs.columns:
            region_col = candidate
            break

    if region_col is not None:
        print(f"  > Region column '{region_col}' found. Computing region-based efficiency...")
        try:
            region_stats = []
            for region_name, sub_obs in adata_reseg.obs.groupby(region_col):
                region_stats.append({
                    'region': region_name,
                    'n_cells': len(sub_obs),
                    'median_transcripts': np.median(sub_obs['total_counts']),
                    'mean_transcripts': np.mean(sub_obs['total_counts']),
                    'median_genes': np.median(sub_obs['n_genes_by_counts']),
                    'mean_genes': np.mean(sub_obs['n_genes_by_counts']),
                })
            df_region = pd.DataFrame(region_stats)
            df_region.to_csv(os.path.join(output_dir, 'efficiency_region_breakdown.csv'), index=False)

            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sns.boxplot(data=adata_reseg.obs, x=region_col, y='total_counts', ax=axes[0])
            axes[0].set_title('Transcripts per Cell by Region')
            axes[0].tick_params(axis='x', rotation=45)
            sns.boxplot(data=adata_reseg.obs, x=region_col, y='n_genes_by_counts', ax=axes[1])
            axes[1].set_title('Genes per Cell by Region')
            axes[1].tick_params(axis='x', rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'efficiency_region_breakdown.png'))
            plt.close()
            print(f"  > Region-based efficiency saved.")
        except Exception as e:
            print(f"  > Error in region-based efficiency: {e}")
    else:
        print("  > No region/spatial_annotation column found. Skipping region-based analysis.")


def analyze_specificity(adata_reseg, config, output_dir, adata_orig=None):
    """NMP score if scRNAseq reference available, otherwise gene-gene correlation proxy."""
    datasets = [('Resegmented', adata_reseg)]
    if adata_orig is not None:
        datasets.append(('Original', adata_orig))

    # --- 4-3a. Negative marker purity (coexpression) ---
    sc_ref_path = config.get('sc_reference_path')
    if sc_ref_path and os.path.exists(sc_ref_path):
        print(f"  > Reference scRNAseq found at {sc_ref_path}. Calculating Negative Marker Purity (NMP)...")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)

            with open(os.path.join(output_dir, 'specificity_nmp_score.txt'), 'w') as f:
                f.write(f"Reference: {sc_ref_path}\n")

                for label, ad in datasets:
                    nmp_score = _negative_marker_purity_coexpression(ad, adata_sc, pipeline_output=True)
                    print(f"  > [{label}] NMP Score: {nmp_score}")
                    f.write(f"[{label}] Negative Marker Purity (NMP) Score: {nmp_score}\n")

        except Exception as e:
            print(f"  > Error calculating NMP: {e}. Proceeding with correlation proxy.")
    else:
        print("  > Using Gene Correlation proxy (Reference scRNAseq not available for NMP).")

    # --- 4-3b. Gene-gene correlation heatmap ---
    print("  > Calculating Co-expression (Gene-Gene correlation)...")

    # Top 50 expressed genes from resegmented data for consistent comparison
    if 'total_counts' not in adata_reseg.var.columns:
         adata_reseg.var['total_counts'] = adata_reseg.X.sum(axis=0).A1 if hasattr(adata_reseg.X, 'toarray') else adata_reseg.X.sum(axis=0)

    top_genes = adata_reseg.var['total_counts'].sort_values(ascending=False).head(50).index

    for label, ad in datasets:
        adata_subset = ad[:, top_genes]
        if isinstance(adata_subset.X, np.ndarray):
            X = adata_subset.X
        else:
            try:
                X = adata_subset.X.toarray()
            except:
                X = adata_subset.X

        corr_matrix = np.corrcoef(X, rowvar=False)

        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, xticklabels=top_genes, yticklabels=top_genes, cmap='coolwarm', center=0)
        plt.title(f'Gene-Gene Correlation ({label})')
        safe_label = label.lower().replace(' ', '_')
        plt.savefig(os.path.join(output_dir, f'specificity_gene_correlation_{safe_label}.png'))
        plt.close()

    print("  > Gene correlation heatmaps saved.")


def analyze_positivity(adata_reseg, output_dir, sample_tag, adata_orig=None):
    """Gene detection rates (fraction of positive cells) + cluster-level violin plots."""
    datasets = [('Resegmented', adata_reseg)]
    if adata_orig is not None:
        datasets.append(('Original', adata_orig))

    plot_data = []

    for label, ad in datasets:
        if 'n_cells_by_counts' not in ad.var.columns:
            sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)

        ad.var['positivity'] = ad.var['n_cells_by_counts'] / ad.n_obs

        top_positive = ad.var.sort_values('positivity', ascending=False).head(50)
        safe_label = label.lower().replace(' ', '_')
        top_positive[['n_cells_by_counts', 'positivity']].to_csv(os.path.join(output_dir, f"{sample_tag}_{safe_label}_gene_positivity.csv"))

        df = pd.DataFrame({'Positivity': ad.var['positivity'], 'Dataset': label})
        plot_data.append(df)

    df_plot = pd.concat(plot_data)

    # --- 4-4a. Positivity histogram ---
    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Positivity', hue='Dataset', bins=50, kde=True, element="step")
    plt.title("Gene Positivity Distribution (Comparison)")
    plt.xlabel("Fraction of Positive Cells")
    plt.ylabel("Number of Genes")
    plt.savefig(os.path.join(output_dir, f"{sample_tag}_positivity_dist_comparison.png"))
    plt.close()

    # --- 4-4b. Preprocessing + Leiden clustering ---
    print("  > Running preprocessing for cluster-level positivity analysis...")
    try:
        ad_proc = adata_reseg.copy()
        sc.pp.normalize_total(ad_proc)
        sc.pp.log1p(ad_proc)
        sc.pp.pca(ad_proc)
        sc.pp.neighbors(ad_proc)
        sc.tl.leiden(ad_proc, resolution=1.0, key_added='leiden')
        sc.tl.umap(ad_proc)

        top_pos_genes = adata_reseg.var.sort_values('positivity', ascending=False).head(10).index.tolist()

        # --- 4-4c. Violin plots per cluster ---
        if len(top_pos_genes) > 0:
            fig, axes = plt.subplots(
                len(top_pos_genes), 1,
                figsize=(12, 3 * len(top_pos_genes)),
                squeeze=False,
            )
            for i, gene in enumerate(top_pos_genes):
                sc.pl.violin(ad_proc, keys=gene, groupby='leiden', ax=axes[i, 0], show=False)
                axes[i, 0].set_title(f'{gene} (positivity={adata_reseg.var.loc[gene, "positivity"]:.3f})')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{sample_tag}_positivity_violin_clusters.png"), dpi=150)
            plt.close()
            print(f"  > Violin plot for top {len(top_pos_genes)} positive genes saved.")

        # --- 4-4d. UMAP per top gene ---
        n_umap_genes = min(4, len(top_pos_genes))
        if n_umap_genes > 0:
            fig, axes = plt.subplots(1, n_umap_genes, figsize=(5 * n_umap_genes, 4))
            if n_umap_genes == 1:
                axes = [axes]
            for i, gene in enumerate(top_pos_genes[:n_umap_genes]):
                sc.pl.umap(ad_proc, color=gene, ax=axes[i], show=False, title=gene)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{sample_tag}_positivity_umap_top_genes.png"), dpi=150)
            plt.close()
            print("  > UMAP plots for top positive genes saved.")

        # Optimal cluster per gene (highest mean expression)
        cluster_means = pd.DataFrame(index=ad_proc.var_names)
        for cluster_id in ad_proc.obs['leiden'].unique():
            mask = ad_proc.obs['leiden'] == cluster_id
            subset_X = ad_proc[mask].X
            if hasattr(subset_X, 'toarray'):
                subset_X = subset_X.toarray()
            cluster_means[cluster_id] = np.mean(subset_X, axis=0)

        optimal_cluster = cluster_means.idxmax(axis=1)
        optimal_expr = cluster_means.max(axis=1)
        df_optimal = pd.DataFrame({
            'gene': optimal_cluster.index,
            'optimal_cluster': optimal_cluster.values,
            'max_cluster_mean_expr': optimal_expr.values,
        })
        df_optimal.to_csv(os.path.join(output_dir, f"{sample_tag}_positivity_optimal_cluster.csv"), index=False)
        print("  > Optimal cluster per gene table saved.")

    except Exception as e:
        print(f"  > Error in cluster-level positivity analysis: {e}")

    print(f"  > Positivity analysis done.")


def analyze_diffusion(df_reseg, adata_reseg, output_dir, df_orig=None, adata_orig=None, technology='xenium'):
    """Transcript-to-centroid distance distributions with complementary CDF, per-gene ECDF, and heatmap."""

    conversion_factor = PIXEL_TO_UM.get(technology.lower(), PIXEL_TO_UM['xenium'])
    print(f"  > Technology: {technology}, pixel-to-um factor: {conversion_factor}")

    def calculate_distances(df, adata, label):
        """Per-transcript distance to assigned cell centroid, converted to micrometers."""
        if 'x_global_px' in df.columns:
            x_col, y_col = 'x_global_px', 'y_global_px'
        elif 'global_x' in df.columns:
            x_col, y_col = 'global_x', 'global_y'
        elif 'x_location' in df.columns:
             x_col, y_col = 'x_location', 'y_location'
        else:
             print(f"  > [{label}] Skipping diffusion: cannot find coordinate columns in transcripts.")
             return None

        if 'x_centroid' not in adata.obs.columns:
            print(f"  > [{label}] Skipping diffusion: x_centroid not in adata.obs")
            return None

        centroid_map_x = adata.obs['x_centroid'].to_dict()
        centroid_map_y = adata.obs['y_centroid'].to_dict()

        cell_col = None
        if 'cell_id_reseg' in df.columns:
            cell_col = 'cell_id_reseg'
        elif 'cell_id' in df.columns:
            cell_col = 'cell_id'
        else:
             print(f"  > [{label}] Skipping: cell id column not found.")
             return None

        # Align index types (int vs str) between transcripts and adata
        first_idx = next(iter(centroid_map_x))
        df_proc = df.copy()

        map_col = cell_col
        if isinstance(first_idx, str) and pd.api.types.is_numeric_dtype(df_proc[cell_col]):
             df_proc['cell_temp_idx'] = df_proc[cell_col].astype(str)
             map_col = 'cell_temp_idx'

        df_proc['cell_x'] = df_proc[map_col].map(centroid_map_x)
        df_proc['cell_y'] = df_proc[map_col].map(centroid_map_y)

        df_proc = df_proc.dropna(subset=['cell_x', 'cell_y'])

        y_vals = df_proc[y_col].values
        x_vals = df_proc[x_col].values

        dy = y_vals - df_proc['cell_y'].values
        dx = x_vals - df_proc['cell_x'].values

        dists = np.sqrt(dy**2 + dx**2)
        dists = dists / conversion_factor

        gene_col = None
        for candidate in ['feature_name', 'gene', 'gene_name', 'target']:
            if candidate in df_proc.columns:
                gene_col = candidate
                break

        gene_names = df_proc[gene_col].values if gene_col else None

        return pd.DataFrame({
            'Distance_um': dists,
            'Dataset': label,
            'Gene': gene_names if gene_names is not None else 'unknown',
        })

    # --- 4-5a. Transcript-to-centroid distance (px→µm) ---
    data_list = []

    df_dists_reseg = calculate_distances(df_reseg, adata_reseg, 'Resegmented')
    if df_dists_reseg is not None:
        data_list.append(df_dists_reseg)
        df_dists_reseg['Distance_um'].describe().to_csv(os.path.join(output_dir, 'reseg_diffusion_stats.csv'))

    if df_orig is not None and adata_orig is not None:
        df_dists_orig = calculate_distances(df_orig, adata_orig, 'Original')
        if df_dists_orig is not None:
            data_list.append(df_dists_orig)
            df_dists_orig['Distance_um'].describe().to_csv(os.path.join(output_dir, 'original_diffusion_stats.csv'))

    if not data_list:
        return

    df_all = pd.concat(data_list, ignore_index=True)

    # --- 4-5b. Complementary CDF plot ---
    df_plot = df_all.copy()
    if len(df_plot) > 50000:
        df_plot = df_plot.sample(50000, random_state=42)

    plt.figure(figsize=(8, 6))
    sns.ecdfplot(data=df_plot, x='Distance_um', hue='Dataset', complementary=True)
    plt.title("Transcript-to-Centroid Distance: Complementary CDF")
    plt.xlabel("Distance (um)")
    plt.ylabel("1 - CDF")
    plt.savefig(os.path.join(output_dir, 'diffusion_complementary_cdf_comparison.png'), dpi=150)
    plt.close()
    print("  > Complementary CDF plot saved.")

    # --- 4-5c. Per-gene ECDF subplots ---
    if 'Gene' in df_all.columns and df_all['Gene'].nunique() > 1:
        gene_counts = df_all['Gene'].value_counts()
        top_genes = gene_counts.head(9).index.tolist()

        if len(top_genes) > 0:
            n_genes = len(top_genes)
            ncols = min(3, n_genes)
            nrows = int(np.ceil(n_genes / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

            for idx, gene in enumerate(top_genes):
                row, col = divmod(idx, ncols)
                ax = axes[row][col]
                df_gene = df_plot[df_plot['Gene'] == gene]
                if len(df_gene) > 0:
                    sns.ecdfplot(data=df_gene, x='Distance_um', hue='Dataset', complementary=True, ax=ax)
                ax.set_title(f'{gene} (n={gene_counts[gene]})')
                ax.set_xlabel('Distance (um)')
                ax.set_ylabel('1 - CDF')

            for idx in range(n_genes, nrows * ncols):
                row, col = divmod(idx, ncols)
                axes[row][col].set_visible(False)

            plt.suptitle('Per-Gene Complementary CDF of Transcript-to-Centroid Distance', y=1.02)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'diffusion_per_gene_ecdf.png'), dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  > Per-gene ECDF subplots saved ({len(top_genes)} genes).")

        # --- 4-5d. Gene×Method distance heatmap ---
        datasets_in_data = df_all['Dataset'].unique().tolist()
        top_heatmap_genes = gene_counts.head(30).index.tolist()

        if len(top_heatmap_genes) > 0 and len(datasets_in_data) > 0:
            df_subset = df_all[df_all['Gene'].isin(top_heatmap_genes)]
            pivot = df_subset.groupby(['Gene', 'Dataset'])['Distance_um'].median().reset_index()
            pivot_wide = pivot.pivot(index='Gene', columns='Dataset', values='Distance_um')
            pivot_wide = pivot_wide.loc[pivot_wide.mean(axis=1).sort_values().index]

            plt.figure(figsize=(max(6, len(datasets_in_data) * 2), max(8, len(top_heatmap_genes) * 0.35)))
            sns.heatmap(
                pivot_wide,
                annot=True,
                fmt='.1f',
                cmap='YlOrRd',
                xticklabels=True,
                yticklabels=True,
            )
            plt.title('Median Transcript-to-Centroid Distance (um): Gene x Method')
            plt.xlabel('Method')
            plt.ylabel('Gene')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'diffusion_gene_method_heatmap.png'), dpi=150, bbox_inches='tight')
            plt.close()

            pivot_wide.to_csv(os.path.join(output_dir, 'diffusion_gene_method_median_distances.csv'))
            print(f"  > Gene x Method distance heatmap saved ({len(top_heatmap_genes)} genes).")
    else:
        print("  > Gene column not found or only one gene; skipping per-gene diffusion plots.")

    print("  > Diffusion analysis completed.")
