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


# --- Self-contained NMP/co-expression functions (originally from xb/calculating.py) ---

def _coexpression_calculation(exp, min_exp=0, min_cells=10):
    """Gene-gene co-expression matrix: fraction of positive cells per gene pair."""
    coexpression = pd.DataFrame(index=exp.columns, columns=exp.columns)
    for col in tqdm(exp.columns):
        sel = exp.loc[:, col] > min_exp
        n_expressing = sel.sum()
        if n_expressing < min_cells:
            coexpression.loc[:, col] = np.nan
            continue
        positive_cells = exp.loc[sel, :]
        coexpression.loc[:, col] = np.sum(positive_cells > min_exp) / positive_cells.shape[0]
    coexpression = coexpression.fillna(1)
    return coexpression


def _negative_marker_purity_coexpression(adata_sp, adata_sc, key='celltype', pipeline_output=True, minexp=0.0):
    """Negative marker purity via co-expression.

    Returns
    -------
    If pipeline_output=True: scalar NMP score
    If pipeline_output=False: (nmp, purity_per_gene, purity_per_celltype,
                                lowvals_sc, lowvals_sp, commongenes)
    """
    minimum_exp = 0.05

    # Normalize gene names to lowercase (notebook 3_4, cell-20 pattern)
    adata_sp = adata_sp.copy()
    adata_sc = adata_sc.copy()
    adata_sp.var_names = [g.lower() for g in adata_sp.var_names]
    adata_sc.var_names = [g.lower() for g in adata_sc.var_names]
    # Remove duplicates after lowering
    adata_sp = adata_sp[:, ~adata_sp.var_names.duplicated()]
    adata_sc = adata_sc[:, ~adata_sc.var_names.duplicated()]

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

    commongenes = mean_ct_sc_rel.index
    neg_marker_mask = np.array(mean_ct_sc_rel < minimum_exp)

    if np.sum(neg_marker_mask) < 1:
        print("No negative markers were found in the sc data reference.")
        negative_marker_purity = 'nan'
        if pipeline_output:
            return negative_marker_purity
        else:
            return negative_marker_purity, None, None, None, None, commongenes

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
        return negative_marker_purity, purity_per_gene, purity_per_celltype, lowvals_sc, lowvals_sp, commongenes

from pipeline.utils.spatial_utils import PIXEL_TO_UM_FACTORS as PIXEL_TO_UM


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

    # --- C2-a. Spatial ROI map (cells colored by region) ---
    try:
        region_col = None
        for candidate in ['region_annotation', 'spatial_annotation', 'region', 'tissue_region', 'domain']:
            if candidate in adata_reseg.obs.columns:
                region_col = candidate
                break
        if region_col and 'x_centroid' in adata_reseg.obs.columns and 'y_centroid' in adata_reseg.obs.columns:
            fig, ax = plt.subplots(figsize=(10, 10))
            regions = adata_reseg.obs[region_col].unique()
            cmap = plt.cm.get_cmap('tab20', len(regions))
            for i, r in enumerate(sorted(regions)):
                mask = adata_reseg.obs[region_col] == r
                ax.scatter(adata_reseg.obs.loc[mask, 'x_centroid'],
                           adata_reseg.obs.loc[mask, 'y_centroid'],
                           s=0.5, alpha=0.5, color=cmap(i), label=str(r), rasterized=True)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_title('Spatial ROI Map (cells by region)')
            ax.legend(markerscale=10, fontsize=7, loc='center left', bbox_to_anchor=(1, 0.5))
            fig.tight_layout()
            fig.savefig(os.path.join(figs_dir, 'spatial_roi_map.png'), dpi=150, bbox_inches='tight')
            plt.close(fig)
            print("  > Spatial ROI map saved.")
    except Exception as e:
        print(f"  > Spatial ROI map failed: {e}")
        plt.close('all')

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
            analyze_positivity(adata_reseg, figs_dir, sample_tag, adata_orig, comp_config)
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

    # --- 4-2b. Expression ratio (ST_mean / scRNA_median, raw counts, notebook method) ---
    sc_ref_path = comp_config.get('sc_reference_path')
    minreads = comp_config.get('efficiency_minreads', 1)
    if sc_ref_path and os.path.exists(sc_ref_path):
        print(f"  > scRNAseq reference found. Computing per-gene expression ratio (raw, minreads={minreads})...")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)

            # Use raw counts (no CPM) — notebook method
            ad_st_raw = adata_reseg.copy()
            if 'raw' in ad_st_raw.layers:
                ad_st_raw.X = ad_st_raw.layers['raw'].copy()
            ad_ref_raw = adata_sc.copy()
            if 'raw' in ad_ref_raw.layers:
                ad_ref_raw.X = ad_ref_raw.layers['raw'].copy()

            common_genes = sorted(set(ad_st_raw.var_names) & set(ad_ref_raw.var_names))
            if len(common_genes) > 0:
                X_st = ad_st_raw[:, common_genes].X
                X_sc = ad_ref_raw[:, common_genes].X
                if hasattr(X_st, 'toarray'):
                    X_st = X_st.toarray()
                if hasattr(X_sc, 'toarray'):
                    X_sc = X_sc.toarray()

                # Per-gene: ST mean over expressing cells, scRNA median over expressing cells
                st_means = []
                sc_medians = []
                for gi in range(len(common_genes)):
                    st_col = X_st[:, gi]
                    sc_col = X_sc[:, gi]
                    st_expr = st_col[st_col >= minreads]
                    sc_expr = sc_col[sc_col >= minreads]
                    st_means.append(np.mean(st_expr) if len(st_expr) > 0 else 0.0)
                    sc_medians.append(np.median(sc_expr) if len(sc_expr) > 0 else 0.0)

                st_means = np.array(st_means)
                sc_medians = np.array(sc_medians)
                sc_medians_safe = np.where(sc_medians > 0, sc_medians, np.nan)
                ratio = st_means / sc_medians_safe

                df_ratio = pd.DataFrame({
                    'gene': common_genes,
                    'st_mean_raw': st_means,
                    'sc_median_raw': sc_medians,
                    'expression_ratio': ratio,
                }).dropna().sort_values('expression_ratio', ascending=False)

                df_ratio.to_csv(os.path.join(output_dir, 'efficiency_expression_ratio.csv'), index=False)

                plt.figure(figsize=(8, 6))
                plt.hist(df_ratio['expression_ratio'].clip(upper=5), bins=50, edgecolor='black')
                plt.axvline(x=1.0, color='red', linestyle='--', label='Ratio = 1')
                plt.title('Per-Gene Expression Ratio (ST mean / scRNA median, raw)')
                plt.xlabel('Expression Ratio (clipped at 5)')
                plt.ylabel('Number of Genes')
                plt.legend()
                plt.savefig(os.path.join(output_dir, 'efficiency_expression_ratio.png'))
                plt.close()

                # C2-b: ST vs scRNAseq expression scatter (log-log + identity line)
                try:
                    fig, ax = plt.subplots(figsize=(8, 8))
                    mask_pos = (df_ratio['st_mean_raw'] > 0) & (df_ratio['sc_median_raw'] > 0)
                    df_pos = df_ratio[mask_pos]
                    ax.scatter(np.log10(df_pos['sc_median_raw']), np.log10(df_pos['st_mean_raw']),
                               s=10, alpha=0.5, edgecolors='none')
                    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                            max(ax.get_xlim()[1], ax.get_ylim()[1])]
                    ax.plot(lims, lims, 'r--', alpha=0.7, label='identity')
                    ax.set_xlabel('log10(scRNA median, raw)')
                    ax.set_ylabel('log10(ST mean, raw)')
                    ax.set_title('ST vs scRNAseq Expression (log-log)')
                    ax.legend()
                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, 'efficiency_st_vs_sc_scatter.png'), dpi=150)
                    plt.close(fig)
                    print("  > ST vs scRNAseq scatter saved.")
                except Exception as e:
                    print(f"  > ST vs scRNAseq scatter failed: {e}")
                    plt.close('all')

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
    for candidate in ['region_annotation', 'spatial_annotation', 'region', 'tissue_region', 'domain']:
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

        # C2-c: Region-specific expression ratio boxplot (log2 scale)
        ratio_path = os.path.join(output_dir, 'efficiency_expression_ratio.csv')
        if os.path.exists(ratio_path) and region_col is not None:
            try:
                df_ratio_loaded = pd.read_csv(ratio_path)
                df_ratio_loaded = df_ratio_loaded[df_ratio_loaded['expression_ratio'] > 0]
                df_ratio_loaded['log2_ratio'] = np.log2(df_ratio_loaded['expression_ratio'])

                fig, ax = plt.subplots(figsize=(10, 6))
                sns.boxplot(data=df_ratio_loaded, y='log2_ratio', ax=ax,
                            boxprops=dict(alpha=0.3))
                sns.stripplot(data=df_ratio_loaded, y='log2_ratio', ax=ax,
                              size=2, alpha=0.4, jitter=0.3)
                ax.axhline(y=0, color='red', linestyle='--', alpha=0.7)
                ax.set_ylabel('log2(ST mean / scRNA median)')
                ax.set_title('Expression Ratio Distribution (log2 scale)')
                fig.tight_layout()
                fig.savefig(os.path.join(output_dir, 'efficiency_ratio_boxplot_log2.png'), dpi=150)
                plt.close(fig)
                print("  > Region ratio boxplot (log2) saved.")
            except Exception as e:
                print(f"  > Region ratio boxplot failed: {e}")
                plt.close('all')
    else:
        print("  > No region/spatial_annotation column found. Skipping region-based analysis.")

    # C2-d: Reseg vs Original boxplot (genes/cell, counts/cell)
    if adata_orig is not None:
        try:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            box_data = pd.concat([
                pd.DataFrame({'value': adata_reseg.obs['total_counts'], 'Dataset': 'Resegmented'}),
                pd.DataFrame({'value': adata_orig.obs['total_counts'], 'Dataset': 'Original'}),
            ])
            sns.boxplot(data=box_data, x='Dataset', y='value', ax=axes[0])
            axes[0].set_ylabel('Counts per Cell')
            axes[0].set_title('Total Counts per Cell')

            box_data_g = pd.concat([
                pd.DataFrame({'value': adata_reseg.obs['n_genes_by_counts'], 'Dataset': 'Resegmented'}),
                pd.DataFrame({'value': adata_orig.obs['n_genes_by_counts'], 'Dataset': 'Original'}),
            ])
            sns.boxplot(data=box_data_g, x='Dataset', y='value', ax=axes[1])
            axes[1].set_ylabel('Genes per Cell')
            axes[1].set_title('Genes Detected per Cell')
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, 'efficiency_reseg_vs_original_boxplot.png'), dpi=150)
            plt.close(fig)
            print("  > Reseg vs Original boxplot saved.")
        except Exception as e:
            print(f"  > Reseg vs Original boxplot failed: {e}")
            plt.close('all')


def analyze_specificity(adata_reseg, config, output_dir, adata_orig=None):
    """NMP score if scRNAseq reference available, otherwise gene-gene correlation proxy."""
    datasets = [('Resegmented', adata_reseg)]
    if adata_orig is not None:
        datasets.append(('Original', adata_orig))

    # --- 4-3a. Negative marker purity (coexpression) ---
    sc_ref_path = config.get('sc_reference_path')
    purity_dfs = []  # Collect per-gene purity for Eff vs Spec scatter
    if sc_ref_path and os.path.exists(sc_ref_path):
        print(f"  > Reference scRNAseq found at {sc_ref_path}. Calculating Negative Marker Purity (NMP)...")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)

            # Load efficiency ratio CSV for gene filtering (ratio < 10) per notebook
            ratio_path = os.path.join(output_dir, 'efficiency_expression_ratio.csv')
            ratio_filter_genes = None
            if os.path.exists(ratio_path):
                df_ratio = pd.read_csv(ratio_path)
                ratio_filter_genes = set(df_ratio[df_ratio['expression_ratio'] < 10]['gene'].values)
                print(f"  > Filtering to {len(ratio_filter_genes)} genes with efficiency ratio < 10")

            with open(os.path.join(output_dir, 'specificity_nmp_score.txt'), 'w') as f:
                f.write(f"Reference: {sc_ref_path}\n")

                for label, ad in datasets:
                    # Apply gene filtering if available
                    ad_filtered = ad
                    if ratio_filter_genes is not None:
                        common = [g for g in ad.var_names if g in ratio_filter_genes]
                        if len(common) > 10:
                            ad_filtered = ad[:, common]
                            print(f"  > [{label}] Using {len(common)} filtered genes for NMP")

                    # Full return: nmp, purity_per_gene, purity_per_celltype, lowvals_sc, lowvals_sp, commongenes
                    result = _negative_marker_purity_coexpression(
                        ad_filtered, adata_sc, pipeline_output=False)
                    nmp_score = result[0]
                    purity_per_gene = result[1]
                    purity_per_celltype = result[2]

                    print(f"  > [{label}] NMP Score: {nmp_score}")
                    f.write(f"[{label}] Negative Marker Purity (NMP) Score: {nmp_score}\n")

                    # Save per-gene purity CSV
                    if purity_per_celltype is not None:
                        safe_label = label.lower().replace(' ', '_')
                        df_purity = pd.DataFrame({
                            'gene': purity_per_celltype.index,
                            'purity': 1 - purity_per_celltype.values,
                            'method': label,
                        })
                        purity_csv = os.path.join(output_dir, f'specificity_nmp_per_gene_{safe_label}.csv')
                        df_purity.to_csv(purity_csv, index=False)
                        print(f"  > [{label}] Per-gene purity saved to {purity_csv}")
                        purity_dfs.append(df_purity)

            # --- 4-3a-2. Per-gene NMP boxplot (Reseg vs Original) ---
            if purity_dfs:
                df_purity_all = pd.concat(purity_dfs, ignore_index=True)
                df_purity_all.to_csv(os.path.join(output_dir, 'specificity_nmp_per_gene_all.csv'), index=False)

                plt.figure(figsize=(10, 5))
                sns.boxplot(data=df_purity_all, x='method', y='purity', boxprops=dict(alpha=0.3))
                sns.stripplot(data=df_purity_all, x='method', y='purity',
                              edgecolor='black', linewidth=0.1, s=3, jitter=0.2)
                plt.ylim([max(0, df_purity_all['purity'].min() - 0.05), 1.05])
                plt.title('Negative Marker Purity (NMP) per Gene')
                plt.ylabel('Purity Score')
                plt.xlabel('Dataset')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'specificity_nmp_per_gene_boxplot.png'), dpi=150)
                plt.close()
                print("  > NMP per-gene boxplot saved.")

            # --- 4-3a-3. Efficiency vs Specificity scatter ---
            if purity_dfs and os.path.exists(ratio_path):
                try:
                    df_ratio = pd.read_csv(ratio_path)
                    for df_pur in purity_dfs:
                        label = df_pur['method'].iloc[0]
                        merged = df_ratio.merge(df_pur, left_on='gene', right_on='gene', how='inner')
                        if len(merged) > 0:
                            plt.figure(figsize=(8, 8))
                            plt.scatter(merged['expression_ratio'], merged['purity'],
                                        s=15, alpha=0.7, edgecolors='black', linewidth=0.3)
                            plt.xlabel('Efficiency Ratio (ST / scRNAseq)')
                            plt.ylabel('Specificity (Purity)')
                            plt.title(f'Efficiency vs Specificity [{label}]')
                            plt.axhline(y=1.0, color='grey', linestyle='--', alpha=0.5)
                            plt.axvline(x=1.0, color='red', linestyle='--', alpha=0.5)
                            safe_label = label.lower().replace(' ', '_')
                            plt.tight_layout()
                            plt.savefig(os.path.join(output_dir,
                                        f'specificity_vs_efficiency_scatter_{safe_label}.png'), dpi=150)
                            plt.close()
                            print(f"  > [{label}] Efficiency vs Specificity scatter saved ({len(merged)} genes).")
                except Exception as e:
                    print(f"  > Error creating Efficiency vs Specificity scatter: {e}")

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


def analyze_positivity(adata_reseg, output_dir, sample_tag, adata_orig=None, comp_config=None):
    """Gene detection rates (fraction of positive cells) + cluster-level violin plots."""
    if comp_config is None:
        comp_config = {}
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
    # Configurable params matching notebook defaults: n_neighbors=8, n_pcs=0, resolution=2.2, min_dist=0.1
    pos_config = comp_config.get('positivity', {})
    n_neighbors = pos_config.get('n_neighbors', 8)
    n_pcs = pos_config.get('n_pcs', 0)
    leiden_resolution = pos_config.get('leiden_resolution', 2.2)
    umap_min_dist = pos_config.get('umap_min_dist', 0.1)
    min_counts = pos_config.get('min_counts', 10)
    min_genes = pos_config.get('min_genes', 3)

    print(f"  > Running preprocessing (n_neighbors={n_neighbors}, n_pcs={n_pcs}, "
          f"resolution={leiden_resolution}, min_dist={umap_min_dist})...")
    try:
        ad_proc = adata_reseg.copy()
        # Cell filtering (matching notebook)
        sc.pp.filter_cells(ad_proc, min_counts=min_counts)
        sc.pp.filter_cells(ad_proc, min_genes=min_genes)
        print(f"  > After filtering: {ad_proc.n_obs} cells (min_counts={min_counts}, min_genes={min_genes})")

        ad_proc.layers['raw'] = ad_proc.X.copy()
        sc.pp.normalize_total(ad_proc, target_sum=None)
        sc.pp.log1p(ad_proc)
        if n_pcs > 0:
            sc.pp.pca(ad_proc)
            sc.pp.neighbors(ad_proc, n_neighbors=n_neighbors, n_pcs=n_pcs)
        else:
            sc.pp.neighbors(ad_proc, n_neighbors=n_neighbors, n_pcs=0)
        sc.tl.leiden(ad_proc, resolution=leiden_resolution, key_added='leiden')
        sc.tl.umap(ad_proc, min_dist=umap_min_dist)

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

    # --- 4-5a+. Per-gene diffusion summary CSV ---
    if 'Gene' in df_all.columns and df_all['Gene'].nunique() > 1:
        gene_distance_summary = df_all.groupby(['Gene', 'Dataset'])['Distance_um'].agg(
            mean_distance_um='mean',
            median_distance_um='median',
            std_distance_um='std',
            n_transcripts='count'
        ).reset_index()
        gene_distance_summary.columns = ['feature_name', 'method', 'mean_distance_um',
                                          'median_distance_um', 'std_distance_um', 'n_transcripts']
        summary_path = os.path.join(output_dir, "diffusion_per_gene_summary.csv")
        gene_distance_summary.to_csv(summary_path, index=False)
        print(f"  > Per-gene diffusion summary saved to {summary_path} ({len(gene_distance_summary)} rows).")

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
            pivot = df_subset.groupby(['Gene', 'Dataset'])['Distance_um'].mean().reset_index()
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
            plt.title('Mean Transcript-to-Centroid Distance (um): Gene x Method')
            plt.xlabel('Method')
            plt.ylabel('Gene')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'diffusion_gene_method_heatmap.png'), dpi=150, bbox_inches='tight')
            plt.close()

            pivot_wide.to_csv(os.path.join(output_dir, 'diffusion_gene_method_mean_distances.csv'))
            print(f"  > Gene x Method distance heatmap saved ({len(top_heatmap_genes)} genes).")
    else:
        print("  > Gene column not found or only one gene; skipping per-gene diffusion plots.")

    # --- 4-5e. Assigned reads stacked barplot ---
    try:
        assignment_data = []
        for label, df in [('Resegmented', df_reseg), ('Original', df_orig)]:
            if df is None:
                continue
            # Use in_cell column if available, otherwise cell_id_reseg/cell_id > 0
            if 'in_cell' in df.columns:
                n_assigned = (df['in_cell'] > 0).sum()
                n_total = len(df)
            elif 'cell_id_reseg' in df.columns:
                n_assigned = (df['cell_id_reseg'] > 0).sum()
                n_total = len(df)
            elif 'cell_id' in df.columns:
                n_assigned = (df['cell_id'] > 0).sum()
                n_total = len(df)
            else:
                continue
            assignment_data.append({
                'Dataset': label,
                'In Cell': n_assigned / n_total,
                'Unassigned': 1 - (n_assigned / n_total),
            })

        if assignment_data:
            df_assign = pd.DataFrame(assignment_data).set_index('Dataset')
            df_assign.plot(kind='bar', stacked=True, figsize=(6, 5),
                           color=['#4DA1A9', '#e8e8e8'], edgecolor='black')
            plt.title('Proportion of Reads Assigned to Cells')
            plt.ylabel('Fraction')
            plt.ylim([0, 1.05])
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'diffusion_assigned_reads_barplot.png'), dpi=150)
            plt.close()
            print("  > Assigned reads stacked barplot saved.")
    except Exception as e:
        print(f"  > Error creating assigned reads barplot: {e}")

    print("  > Diffusion analysis completed.")
