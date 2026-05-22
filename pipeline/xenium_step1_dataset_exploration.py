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
#   1-3.5. MTRNR gene proportion visualization
#   1-4. Clustering & marker annotation
#     1-4a. Normalization (normalize_total + log1p)
#     1-4b. HVG selection
#     1-4c. PCA + neighbors (on normalized data)
#     1-4d. Leiden clustering (multi-resolution)
#     1-4e. Marker genes (Wilcoxon) → dotplot + heatmap
#     1-4f. UMAP + spatial scatter
#   1-5. Neighborhood analysis
#     1-5a. Spatial neighbors graph
#     1-5b. Neighborhood diversity
#     1-5c. Enrichment analysis
#     1-5d. Centrality scores (barplot)
#     1-5e. Spatial graph visualization (cropped region, node-edge)
#     1-5f. Cluster network diagram (nodes=clusters, size=centrality)
#     1-5g. Spatial centrality heatmap (tissue coords colored by centrality)

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

def _add_ecdf_guidelines(ax, distances, thresholds=(0.5, 0.9),
                         colors=('red', 'blue')):
    """Adds horizontal + vertical guide lines at ECDF quantile thresholds."""
    sorted_d = np.sort(distances)
    for thresh, color in zip(thresholds, colors):
        q_val = np.percentile(sorted_d, thresh * 100)
        # Horizontal line
        ax.axhline(thresh, color=color, linestyle='--', alpha=0.6, linewidth=1)
        # Vertical dotted line from threshold down to x-axis
        ax.plot([q_val, q_val], [0, thresh], color=color, linestyle=':',
                alpha=0.7, linewidth=1.2)
        # Dot at intersection
        ax.plot(q_val, thresh, 'o', color=color, markersize=5, zorder=5)
        # Distance label below x-axis
        ax.annotate(f'{q_val:.1f} µm',
                    xy=(q_val, 0), xytext=(q_val, -0.04),
                    fontsize=8, color=color, fontweight='bold',
                    ha='center', va='top', annotation_clip=False)
        # Threshold label on the right
        ax.annotate(f'{thresh:.0%}',
                    xy=(ax.get_xlim()[1] * 0.98, thresh),
                    fontsize=8, color=color, fontweight='bold',
                    ha='right', va='bottom')


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


def _load_marker_localization_table(marker_file, allowed_classes=None):
    """Loads marker localization table with required columns: gene, class."""
    if not marker_file:
        print("    [INFO] marker_qc.marker_file is empty. Skipping marker-based QC.")
        return None

    if not os.path.exists(marker_file):
        print(f"    [WARNING] Marker file not found: {marker_file}")
        return None

    try:
        marker_df = pd.read_csv(marker_file, sep=None, engine='python')
    except Exception as e:
        print(f"    [WARNING] Failed to read marker file '{marker_file}': {e}")
        return None

    marker_df.columns = [str(c).strip().lower() for c in marker_df.columns]
    required_cols = {'gene', 'class'}
    if not required_cols.issubset(set(marker_df.columns)):
        print(f"    [WARNING] Marker file missing required columns {required_cols}.")
        return None

    marker_df['gene'] = marker_df['gene'].astype(str).str.strip()
    marker_df['class'] = marker_df['class'].astype(str).str.strip().str.lower()
    marker_df = marker_df[marker_df['gene'] != ""].drop_duplicates(subset=['gene'], keep='first')

    if allowed_classes:
        allowed = {c.lower() for c in allowed_classes}
        marker_df = marker_df[marker_df['class'].isin(allowed)]
        if marker_df.empty:
            print("    [WARNING] No markers left after filtering by allowed_classes.")
            return None

    print(f"    - Loaded marker table: {len(marker_df)} genes from {marker_file}")
    return marker_df


def _prepare_marker_qc_table(spots_assigned, marker_df, min_transcripts_per_gene=30):
    """Builds transcript-level table joined with marker class and gene counts."""
    if marker_df is None:
        return None

    if 'feature_name' not in spots_assigned.columns:
        print("    [WARNING] 'feature_name' not available in transcripts table. Skipping marker-based QC.")
        return None

    df = spots_assigned[['feature_name', 'distance']].copy()
    df['gene'] = df['feature_name'].astype(str)

    gene_counts = df['gene'].value_counts().rename('n_transcripts_gene')
    df = df.join(gene_counts, on='gene')
    df = df[df['n_transcripts_gene'] >= int(min_transcripts_per_gene)].copy()
    if df.empty:
        print(f"    [WARNING] No transcripts left after min_transcripts_per_gene={min_transcripts_per_gene}.")
        return None

    marker_map = marker_df.rename(columns={'class': 'marker_class'})
    df = df.merge(marker_map[['gene', 'marker_class']], on='gene', how='inner')
    if df.empty:
        print("    [WARNING] No overlap between marker genes and transcript genes.")
        return None

    return df


def _compute_marker_qc_metrics(df, thresholds, pass_rules, qc_classes):
    """Computes class-level localization proportions and PASS/WARN status."""
    nuclear_max = float(thresholds.get('nuclear_max', 5.0))
    cyto_max = float(thresholds.get('cyto_max', 10.0))
    qc_classes = {c.lower() for c in qc_classes}

    rows = []
    for marker_class, sub in df.groupby('marker_class', sort=True):
        n_transcripts = len(sub)
        n_genes = sub['gene'].nunique()
        p_le_5 = (sub['distance'] <= nuclear_max).mean()
        p_5_10 = ((sub['distance'] > nuclear_max) & (sub['distance'] <= cyto_max)).mean()
        p_gt_10 = (sub['distance'] > cyto_max).mean()

        status = 'INFO'
        messages = []
        if marker_class in qc_classes:
            tail_max = float(pass_rules.get('tail_gt10um_max', 0.20))
            if marker_class == 'nuclear':
                target = float(pass_rules.get('nuclear_in_5um_min', 0.70))
                if p_le_5 < target:
                    messages.append(f"P<=5µm {p_le_5:.2f} < {target:.2f}")
            elif marker_class == 'cytoplasmic':
                target = float(pass_rules.get('cyto_in_5_10um_min', 0.45))
                if p_5_10 < target:
                    messages.append(f"P5-10µm {p_5_10:.2f} < {target:.2f}")
            if p_gt_10 > tail_max:
                messages.append(f"P>10µm {p_gt_10:.2f} > {tail_max:.2f}")
            status = 'PASS' if not messages else 'WARN'
        else:
            messages.append("Reference class (no PASS/WARN rule)")

        rows.append({
            'class': marker_class,
            'n_genes': int(n_genes),
            'n_transcripts': int(n_transcripts),
            'p_le_5': p_le_5,
            'p_5_10': p_5_10,
            'p_gt_10': p_gt_10,
            'status': status,
            'rule_message': '; '.join(messages),
        })

    metrics_df = pd.DataFrame(rows).sort_values(by='class')
    return metrics_df


def _plot_marker_group_violin(df, output_dir, sample_tag, thresholds):
    """Plots per-gene violin grouped by marker class."""
    if df.empty:
        return None

    gene_counts = df['gene'].value_counts()
    gene_order = (
        df.groupby(['marker_class', 'gene'])['distance']
          .median()
          .reset_index()
          .sort_values(['marker_class', 'distance'])
    )
    order = gene_order['gene'].tolist()
    class_palette = {
        'nuclear': '#1f77b4',
        'cytoplasmic': '#ff7f0e',
        'mitochondrial': '#2ca02c',
        'secreted': '#d62728',
    }

    plot_df = df.copy()
    plot_df['gene_label'] = plot_df['gene'].map(lambda g: f"{g} (n={gene_counts.get(g, 0):,})")
    label_order = [f"{g} (n={gene_counts.get(g, 0):,})" for g in order]

    fig, ax = plt.subplots(figsize=(max(14, len(label_order) * 0.5), 7))
    sns.violinplot(
        data=plot_df,
        x='gene_label',
        y='distance',
        hue='marker_class',
        order=label_order,
        cut=0,
        scale='width',
        dodge=False,
        palette=class_palette,
        ax=ax,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_xlabel("Marker Gene (with transcript count)")
    ax.set_ylabel("distance_to_nucleus (µm)")
    ax.set_title("Marker-Based Localization QC: Distance Distribution by Gene")

    nuclear_max = float(thresholds.get('nuclear_max', 5.0))
    cyto_max = float(thresholds.get('cyto_max', 10.0))
    for y, label in [(nuclear_max, f"{nuclear_max:g} µm"), (cyto_max, f"{cyto_max:g} µm")]:
        ax.axhline(y=y, color='gray', linestyle='--', linewidth=1, alpha=0.6)
        ax.text(ax.get_xlim()[1] + 0.2, y, label, va='center', fontsize=8, color='gray', clip_on=False)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, title='Marker class', bbox_to_anchor=(1.02, 1), loc='upper left')

    fig.subplots_adjust(right=0.85, bottom=0.3)
    save_path = os.path.join(output_dir, f"{sample_tag}_step1_marker_localization_violin.png")
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"    - Saved marker localization violin to: {save_path}")
    return save_path


def _plot_marker_group_ecdf(df, output_dir, sample_tag, thresholds):
    """Plots class-level ECDF with threshold lines."""
    if df.empty:
        return None

    class_palette = {
        'nuclear': '#1f77b4',
        'cytoplasmic': '#ff7f0e',
        'mitochondrial': '#2ca02c',
        'secreted': '#d62728',
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.ecdfplot(data=df, x='distance', hue='marker_class', ax=ax, palette=class_palette, linewidth=2)
    nuclear_max = float(thresholds.get('nuclear_max', 5.0))
    cyto_max = float(thresholds.get('cyto_max', 10.0))
    for x in [nuclear_max, cyto_max]:
        ax.axvline(x=x, color='gray', linestyle='--', linewidth=1, alpha=0.6)

    ax.set_title("Marker-Based Localization QC: ECDF by Marker Class")
    ax.set_xlabel("distance_to_nucleus (µm)")
    ax.set_ylabel("Cumulative proportion")
    ax.set_xlim(left=0)
    ax.legend(title='Marker class', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.subplots_adjust(right=0.82)

    save_path = os.path.join(output_dir, f"{sample_tag}_step1_marker_localization_ecdf.png")
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"    - Saved marker localization ECDF to: {save_path}")
    return save_path


def _plot_marker_threshold_stacked_bar(metrics_df, output_dir, sample_tag):
    """Plots class-level stacked proportions for <=5, 5-10, >10 µm zones."""
    if metrics_df.empty:
        return None

    bar_df = metrics_df[['class', 'p_le_5', 'p_5_10', 'p_gt_10']].set_index('class')
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ['#4c78a8', '#f58518', '#e45756']
    bar_df.plot(kind='bar', stacked=True, color=colors, ax=ax, width=0.75)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("Marker class")
    ax.set_title("Marker Localization QC: Threshold Zone Proportions")
    ax.legend(['<=5 µm', '5-10 µm', '>10 µm'], title='Distance zone', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.subplots_adjust(right=0.8)

    save_path = os.path.join(output_dir, f"{sample_tag}_step1_marker_localization_thresholds.png")
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"    - Saved marker localization threshold bar plot to: {save_path}")
    return save_path


def _save_marker_qc_reports(metrics_df, output_dir, sample_tag, thresholds, pass_rules):
    """Saves marker QC metrics as csv + txt."""
    csv_path = os.path.join(output_dir, f"{sample_tag}_step1_marker_localization_qc.csv")
    metrics_df.to_csv(csv_path, index=False)

    txt_path = os.path.join(output_dir, f"{sample_tag}_step1_marker_localization_qc.txt")
    with open(txt_path, 'w') as f:
        f.write("Marker-based localization QC summary\n")
        f.write(f"Thresholds (um): {thresholds}\n")
        f.write(f"Pass rules: {pass_rules}\n\n")
        for _, row in metrics_df.iterrows():
            f.write(
                f"[{row['class']}] status={row['status']} | "
                f"genes={int(row['n_genes'])}, transcripts={int(row['n_transcripts'])}, "
                f"P<=5={row['p_le_5']:.3f}, P5-10={row['p_5_10']:.3f}, P>10={row['p_gt_10']:.3f} | "
                f"{row['rule_message']}\n"
            )
    print(f"    - Saved marker QC reports to: {csv_path}, {txt_path}")


def _run_marker_localization_qc(spots_assigned, output_dir, sample_tag, marker_qc_cfg):
    """Runs marker-based localization QC workflow."""
    enabled = bool(marker_qc_cfg.get('enabled', False))
    if not enabled:
        print("    - Marker-based localization QC is disabled.")
        return

    marker_file = marker_qc_cfg.get('marker_file')
    allowed_classes = marker_qc_cfg.get('allowed_classes', ['nuclear', 'cytoplasmic', 'mitochondrial', 'secreted'])
    qc_classes = marker_qc_cfg.get('qc_classes', ['nuclear', 'cytoplasmic'])
    min_genes_per_class = int(marker_qc_cfg.get('min_genes_per_class', 5))
    min_transcripts_per_gene = int(marker_qc_cfg.get('min_transcripts_per_gene', 30))
    thresholds = marker_qc_cfg.get('thresholds_um', {'nuclear_max': 5.0, 'cyto_max': 10.0})
    pass_rules = marker_qc_cfg.get(
        'pass_rules',
        {'nuclear_in_5um_min': 0.70, 'cyto_in_5_10um_min': 0.45, 'tail_gt10um_max': 0.20}
    )

    print("    - Running marker-based localization QC...")
    marker_df = _load_marker_localization_table(marker_file, allowed_classes=allowed_classes)
    if marker_df is None:
        return

    marker_qc_df = _prepare_marker_qc_table(
        spots_assigned,
        marker_df=marker_df,
        min_transcripts_per_gene=min_transcripts_per_gene,
    )
    if marker_qc_df is None:
        return

    class_counts = marker_qc_df.groupby('marker_class')['gene'].nunique()
    low_classes = class_counts[class_counts < min_genes_per_class]
    if len(low_classes) > 0:
        print(f"    [WARNING] Classes with <{min_genes_per_class} genes: {dict(low_classes)}")

    metrics_df = _compute_marker_qc_metrics(
        marker_qc_df,
        thresholds=thresholds,
        pass_rules=pass_rules,
        qc_classes=qc_classes,
    )
    if metrics_df.empty:
        print("    [WARNING] No marker class metrics computed.")
        return

    _plot_marker_group_violin(marker_qc_df, output_dir, sample_tag, thresholds)
    _plot_marker_group_ecdf(marker_qc_df, output_dir, sample_tag, thresholds)
    _plot_marker_threshold_stacked_bar(metrics_df, output_dir, sample_tag)
    _save_marker_qc_reports(metrics_df, output_dir, sample_tag, thresholds, pass_rules)


def calculate_transcript_dispersion(adata, output_dir, sample_tag, spots=None, marker_qc_cfg=None):
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
    q99 = np.percentile(distances, 99)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(distances, bins=100, kde=True, color='purple', ax=ax)
    ax.set_xlim(0, q99)

    # Find KDE peak and annotate
    kde_line = [c for c in ax.get_children()
                if isinstance(c, plt.matplotlib.lines.Line2D)]
    peak_info = ""
    if kde_line:
        kde_x, kde_y = kde_line[0].get_xdata(), kde_line[0].get_ydata()
        peak_idx = np.argmax(kde_y)
        peak_dist = kde_x[peak_idx]
        bin_counts, bin_edges = np.histogram(distances, bins=100, range=(0, q99))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        nearest_bin = np.argmin(np.abs(bin_centers - peak_dist))
        peak_hist_count = bin_counts[nearest_bin]

        ax.axvline(peak_dist, color='red', linestyle='--', alpha=0.7, linewidth=1)
        ax.annotate(f'Peak: {peak_dist:.1f} µm\n(bin count ≈ {peak_hist_count:,})',
                    xy=(peak_dist, kde_y[peak_idx]),
                    xytext=(peak_dist + q99 * 0.08, kde_y[peak_idx] * 0.9),
                    fontsize=9, color='red', fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.2))
        peak_info = f" | Peak: {peak_dist:.1f} µm"

    ax.set_title(f"Transcript Dispersion Distribution\n"
                 f"Metric: {metric_name} | Median: {median_dist:.2f} µm{peak_info} | "
                 f"N={n_assigned:,} transcripts")
    ax.set_xlabel(f"{metric_name} (µm)")
    ax.set_ylabel("Count")

    save_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_dist.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
    print(f"    - Saved Dispersion Histogram to {save_path}")

    # --- 1-2b. ECDF (overall + per-gene top 10) ---
    try:
        print("    - Generating ECDF plots...")

        fig_ecdf, ax_ecdf = plt.subplots(figsize=(10, 6))
        if 'source' in spots_assigned.columns and spots_assigned['source'].nunique() > 1:
            sns.ecdfplot(data=spots_assigned, x='distance', hue='source',
                         complementary=False, ax=ax_ecdf)
            ax_ecdf.set_title("ECDF of Transcript Distance to Centroid (by Source)")
        else:
            sns.ecdfplot(data=spots_assigned, x='distance', complementary=False,
                         color='steelblue', ax=ax_ecdf)
            ax_ecdf.set_title("ECDF of Transcript Distance to Centroid")
        ax_ecdf.set_xlabel("Distance (µm)")
        ax_ecdf.set_ylabel("Cumulative Proportion")
        ax_ecdf.set_xlim(0, q99)
        _add_ecdf_guidelines(ax_ecdf, distances.values)
        ecdf_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_ecdf.png")
        plt.savefig(ecdf_path, bbox_inches='tight')
        plt.close()
        print(f"    - Saved ECDF plot to {ecdf_path}")

        if 'feature_name' in spots_assigned.columns:
            top10_genes = spots_assigned['feature_name'].value_counts().head(10).index
            spots_top10 = spots_assigned[spots_assigned['feature_name'].isin(top10_genes)]
            if len(spots_top10) > 0:
                fig_eg, ax_eg = plt.subplots(figsize=(12, 6))
                sns.ecdfplot(data=spots_top10, x='distance', hue='feature_name',
                             complementary=False, ax=ax_eg)
                ax_eg.set_title("ECDF of Transcript Distance to Centroid (Top 10 Genes)")
                ax_eg.set_xlabel("Distance (µm)")
                ax_eg.set_ylabel("Cumulative Proportion")
                ax_eg.set_xlim(0, q99)
                _add_ecdf_guidelines(ax_eg, distances.values)
                ax_eg.legend(title='Gene', bbox_to_anchor=(1.02, 1),
                             loc='upper left', fontsize=8, title_fontsize=9,
                             framealpha=0.9)
                ecdf_gene_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_ecdf_genes.png")
                plt.savefig(ecdf_gene_path, dpi=200, bbox_inches='tight')
                plt.close()
                print(f"    - Saved gene-level ECDF plot to {ecdf_gene_path}")
    except Exception as e:
        print(f"    [WARNING] ECDF plotting failed: {e}")

    # --- 1-2c. Violin plot (top 20 genes by transcript count) ---
    try:
        if 'feature_name' in spots_assigned.columns:
            print("    - Generating violin plot for top 20 genes (by transcript count)...")
            gene_counts = spots_assigned['feature_name'].value_counts()
            top20_genes = gene_counts.head(20).index
            max_count = gene_counts.iloc[0]
            min_count = gene_counts.iloc[min(19, len(gene_counts) - 1)]
            spots_top20 = spots_assigned[spots_assigned['feature_name'].isin(top20_genes)]
            if len(spots_top20) > 0:
                fig, ax = plt.subplots(figsize=(14, 6))
                sns.violinplot(data=spots_top20, x='feature_name', y='distance',
                               cut=0, scale='width', order=top20_genes, ax=ax)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

                # Guidelines at 5 µm and 10 µm
                for guideline_y, label in [(5, '5 µm'), (10, '10 µm')]:
                    ax.axhline(y=guideline_y, color='gray', linestyle='--',
                               alpha=0.6, linewidth=1)
                    ax.text(ax.get_xlim()[1] + 0.2, guideline_y, label,
                            va='center', fontsize=8, color='gray',
                            fontweight='bold', clip_on=False)

                ax.set_title("Distance to Centroid per Gene (Top 20 by Transcript Count)")
                ax.text(0.5, 1.02,
                        f"Transcript counts: {max_count:,} – {min_count:,}",
                        transform=ax.transAxes, ha='center', fontsize=9, color='gray')
                ax.set_xlabel("Gene")
                ax.set_ylabel(f"{metric_name} (µm)")
                fig.subplots_adjust(right=0.92)
                violin_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_violin.png")
                plt.savefig(violin_path, bbox_inches='tight')
                plt.close()
                print(f"    - Saved violin plot to {violin_path}")

                # Localization-grouped violin (nuclear vs cytoplasmic)
                _plot_dispersion_by_localization(
                    spots_top20, output_dir, sample_tag,
                    distance_col='distance', threshold_um=5.0, n_top_genes=20)
    except Exception as e:
        print(f"    [WARNING] Violin plot generation failed: {e}")

    # Save dispersion metrics summary
    with open(os.path.join(output_dir, f"{sample_tag}_step1_dispersion_metrics.txt"), "w") as f:
        f.write(f"Metric Used: {metric_name}\n")
        f.write(f"Total Transcripts: {n_total}\n")
        f.write(f"Assigned Transcripts: {n_assigned}\n")
        f.write(f"Median Distance: {median_dist}\n")
        f.write(f"Mean Distance: {mean_dist}\n")

    # --- 1-2d. Marker-based localization QC ---
    try:
        _run_marker_localization_qc(
            spots_assigned=spots_assigned,
            output_dir=output_dir,
            sample_tag=sample_tag,
            marker_qc_cfg=(marker_qc_cfg or {}),
        )
    except Exception as e:
        print(f"    [WARNING] Marker-based localization QC failed: {e}")


def _plot_dispersion_by_localization(spots_assigned, output_dir, sample_tag,
                                     distance_col='distance', threshold_um=5.0,
                                     n_top_genes=20):
    """
    Violin plot of transcript distance, genes grouped by nuclear vs cytoplasmic
    localization (median distance < or >= threshold).
    Nuclear genes on left (blue), cytoplasmic on right (orange).
    """
    print("    - Generating localization-grouped dispersion violin plot...")

    gene_counts = spots_assigned['feature_name'].value_counts()
    top_genes = gene_counts.head(n_top_genes).index.tolist()
    spots_top = spots_assigned[spots_assigned['feature_name'].isin(top_genes)]

    if len(spots_top) == 0:
        print("    [WARNING] No data for localization violin plot.")
        return

    # classify by median distance
    medians = spots_top.groupby('feature_name')[distance_col].median()
    nuclear = sorted([g for g in top_genes if medians[g] < threshold_um],
                     key=lambda g: medians[g])
    cytoplasmic = sorted([g for g in top_genes if medians[g] >= threshold_um],
                         key=lambda g: medians[g])

    gene_order = cytoplasmic + nuclear   # cytoplasmic left, nuclear right
    n_nuc = len(nuclear)
    n_cyto = len(cytoplasmic)

    # Build label with transcript count
    gene_counts = spots_top['feature_name'].value_counts()
    label_map = {g: f"{g} (n={gene_counts.get(g, 0):,})" for g in gene_order}
    spots_plot = spots_top.copy()
    spots_plot['gene_label'] = spots_plot['feature_name'].map(label_map)
    spots_plot['marker_class'] = spots_plot['feature_name'].apply(
        lambda g: 'cytoplasmic' if g in cytoplasmic else 'nuclear')
    label_order = [label_map[g] for g in gene_order]

    class_palette = {'cytoplasmic': '#ff7f0e', 'nuclear': '#1f77b4'}

    fig, ax = plt.subplots(figsize=(max(14, len(gene_order) * 0.5), 7))
    sns.violinplot(data=spots_plot, x='gene_label', y=distance_col,
                   hue='marker_class', order=label_order,
                   cut=0, density_norm='width', dodge=False,
                   palette=class_palette, inner='box', ax=ax)

    # Guidelines at 5 µm and 10 µm
    for guideline_y, label in [(threshold_um, f'{threshold_um:.0f} µm'),
                                (10, '10 µm')]:
        ax.axhline(y=guideline_y, color='gray', linestyle='--',
                   alpha=0.6, linewidth=1)
        ax.text(ax.get_xlim()[1] + 0.3, guideline_y, label,
                va='center', fontsize=8, color='gray',
                fontweight='bold', clip_on=False)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_title(f"Transcript Dispersion by Localization\n"
                 f"(Top {n_top_genes} Genes, threshold={threshold_um} µm)")
    ax.set_xlabel("Marker Gene (with transcript count)")
    ax.set_ylabel("distance_to_nucleus (µm)")

    # Legend outside plot area
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, title='Marker class',
                  bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.subplots_adjust(right=0.85, bottom=0.3)

    out_path = os.path.join(output_dir, f"{sample_tag}_step1_dispersion_localization.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"    - Saved localization violin plot to {out_path}")
    print(f"      Cytoplasmic (>={threshold_um} µm): {', '.join(cytoplasmic)}")
    print(f"      Nuclear (<{threshold_um} µm): {', '.join(nuclear)}")


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


# --- 1-3.5. MTRNR gene proportion visualization ---

def _plot_mtrnr_proportion(adata, output_dir, sample_tag):
    """
    Visualizes the proportion of MTRNR2L pseudogene expression that can
    dominate clustering/UMAP.  Produces a 2-panel figure:
      Left  – Top-20 genes by total counts (MTRNR highlighted in red)
      Right – Per-cell violin of pct_counts for each MTRNR gene
    """
    print("\n[Step 1-3.5] MTRNR Gene Proportion Visualization...")

    mtrnr_genes = [g for g in ["MTRNR2L12", "MTRNR2L8"] if g in adata.var_names]
    if not mtrnr_genes:
        print("    [INFO] No MTRNR2L12/MTRNR2L8 found in var_names – skipping.")
        return

    # --- compute total counts per gene ---
    from scipy.sparse import issparse
    X = adata.X
    if issparse(X):
        gene_totals = np.array(X.sum(axis=0)).ravel()
    else:
        gene_totals = np.array(X.sum(axis=0)).ravel()
    total_all = gene_totals.sum()

    gene_total_series = pd.Series(gene_totals, index=adata.var_names)
    top20 = gene_total_series.nlargest(20)

    # global proportions
    mtrnr_props = {}
    for g in mtrnr_genes:
        mtrnr_props[g] = gene_total_series[g] / total_all * 100
    combined_pct = sum(mtrnr_props.values())

    # per-cell pct_counts
    cell_totals = np.array(X.sum(axis=1)).ravel()
    cell_totals[cell_totals == 0] = 1  # avoid div-by-zero
    pct_records = []
    for g in mtrnr_genes:
        idx = list(adata.var_names).index(g)
        if issparse(X):
            gene_vals = np.array(X[:, idx].todense()).ravel()
        else:
            gene_vals = X[:, idx].ravel()
        pcts = gene_vals / cell_totals * 100
        for v in pcts:
            pct_records.append({"gene": g, "pct_counts": v})
    pct_df = pd.DataFrame(pct_records)

    # --- figure ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: top-20 barplot
    ax = axes[0]
    colors = ["#d62728" if g in mtrnr_genes else "#999999" for g in top20.index]
    ax.barh(range(len(top20)), top20.values, color=colors)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20.index, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Total counts")
    ax.set_title("Top 20 Genes by Total Expression")

    textstr = "\n".join(
        [f"{g}: {mtrnr_props[g]:.1f}% of total reads" for g in mtrnr_genes]
        + [f"Combined: {combined_pct:.1f}% of total reads"]
    )
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))

    # Right panel: violin plot
    ax = axes[1]
    sns.violinplot(data=pct_df, x="gene", y="pct_counts", inner="quartile",
                   palette=["#d62728"] * len(mtrnr_genes), ax=ax)
    for g in mtrnr_genes:
        subset = pct_df.loc[pct_df["gene"] == g, "pct_counts"]
        med = subset.median()
        ax.axhline(med, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_ylabel("% of cell total counts")
    ax.set_xlabel("")
    ax.set_title("Per-cell MTRNR Expression (%)")

    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{sample_tag}_step1_mtrnr_proportion.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    - Saved MTRNR proportion plot to: {out_path}")
    for g in mtrnr_genes:
        print(f"    - {g}: {mtrnr_props[g]:.2f}% of total reads")
    print(f"    - Combined: {combined_pct:.2f}% of total reads")


# --- 1-4. Clustering & marker annotation ---

def perform_clustering_and_annotation(adata, output_dir, sample_tag, config):
    """
    Performs Leiden clustering, HVG selection, and marker gene ranking (Ref: 1_2).
    """
    print("\n[Step 1-3] Running Clustering & Annotation (Ref: 1_2)...")

    # --- 1-4a. Normalization (must precede PCA) ---
    if 'raw' not in adata.layers:
        print("    - Saving raw counts backup to adata.layers['raw']...")
        adata.layers['raw'] = adata.X.copy()

    if 'log1p' not in adata.uns:
        print("    - Normalizing per-cell (target_sum=1e4) + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("    - log1p already applied, skipping normalization.")

    # Save normalized data as .raw (used by rank_genes_groups, plotting)
    if adata.raw is None:
        adata.raw = adata.copy()

    # --- 1-4b. HVG selection ---
    print("    - Calculating Highly Variable Genes...")
    try:
        sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)

        # --- Custom HVG scatter plot ---
        hvg_mask = adata.var['highly_variable']
        means = adata.var['means']
        dispersions_norm = adata.var['dispersions_norm']

        n_hvg = int(hvg_mask.sum())
        n_total = len(hvg_mask)
        ratio = n_hvg / n_total if n_total > 0 else 0.0
        n_high_mean_excluded = int((means > 7).sum())

        # HVG plot deferred to after marker gene computation (see below)
        print(f"    - HVG selected: {n_hvg}/{n_total} ({ratio:.1%})")
        print(f"    - High-mean excluded (>7): {n_high_mean_excluded}")
    except Exception as e:
        print(f"    [WARNING] HVG calculation/plotting failed: {e}")

    # --- 1-4c. PCA + neighbors (on normalized data) ---
    if 'X_pca' not in adata.obsm:
        print("    - Running PCA on normalized data...")
        sc.pp.pca(adata)
    if 'neighbors' not in adata.uns:
        print("    - Computing Neighbors...")
        sc.pp.neighbors(adata)

    # --- 1-4d. Leiden clustering (multi-resolution) ---
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

        sc.pl.dotplot(adata, markers_dict, groupby=primary_key, dendrogram=False,
                      standard_scale='var', show=False)
        # Rotate gene labels 45 degrees for readability
        for ax_item in plt.gcf().get_axes():
            for label in ax_item.get_xticklabels():
                if label.get_text():
                    label.set_rotation(45)
                    label.set_ha('right')
        plt.suptitle(f"Top Markers ({primary_key})\n"
                     f"Rows: clusters | Columns: top 3 marker genes per cluster",
                     fontsize=10, y=1.02)
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

    # --- 1-4b (deferred). Enhanced HVG scatter with marker gene overlay ---
    if 'highly_variable' in adata.var.columns:
        try:
            hvg_mask = adata.var['highly_variable']
            means = adata.var['means']
            dispersions_norm = adata.var['dispersions_norm']
            n_hvg = int(hvg_mask.sum())
            n_total = len(hvg_mask)
            ratio = n_hvg / n_total if n_total > 0 else 0.0
            n_high_mean_excluded = int((means > 7).sum())

            # Collect unique marker genes from rank_genes_groups
            marker_genes = set()
            if 'rank_genes_groups' in adata.uns:
                rgg = adata.uns['rank_genes_groups']
                for grp in rgg['names'].dtype.names:
                    marker_genes.update(rgg['names'][grp][:5].tolist())

            hvg_names = set(adata.var_names[hvg_mask])
            markers_in_hvg = marker_genes & hvg_names
            markers_missing = marker_genes - hvg_names
            # Keep only markers that are actually in the gene list
            markers_missing = {g for g in markers_missing if g in adata.var_names}

            fig, ax = plt.subplots(figsize=(10, 7))
            ax.scatter(
                means[~hvg_mask], dispersions_norm[~hvg_mask],
                c='gray', s=20, alpha=0.5, label='Other genes', rasterized=True,
            )
            ax.scatter(
                means[hvg_mask], dispersions_norm[hvg_mask],
                c='red', s=40, alpha=0.8, label='Highly variable', rasterized=True,
            )
            # Overlay marker genes missing from HVG as blue triangles
            if markers_missing:
                miss_idx = [adata.var_names.get_loc(g) for g in markers_missing]
                ax.scatter(
                    means.iloc[miss_idx], dispersions_norm.iloc[miss_idx],
                    c='blue', s=60, alpha=0.9, marker='^',
                    label=f'Markers not in HVG ({len(markers_missing)})',
                    zorder=5,
                )

            ax.set_xlabel('Mean expression')
            ax.set_ylabel('Normalized dispersion')
            ax.set_title(f"Highly Variable Genes: {n_hvg}/{n_total} selected ({ratio:.1%})")
            ax.legend(loc='upper right', fontsize=8)

            # Info text box with marker gene coverage
            lines = [f"High mean excluded (>7): {n_high_mean_excluded}"]
            if marker_genes:
                n_mk = len(marker_genes)
                n_in = len(markers_in_hvg)
                lines.append(f"Marker genes in HVG: {n_in}/{n_mk} ({n_in/n_mk:.0%})")
            textstr = '\n'.join(lines)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=9,
                    verticalalignment='top', bbox=props)

            hvg_file = os.path.join(output_dir, f"{sample_tag}_step1_hvg.png")
            fig.savefig(hvg_file, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"    - Saved HVG Plot to: {hvg_file}")
        except Exception as e:
            print(f"    [WARNING] HVG plot failed: {e}")

    # --- 1-4d-2. Cluster composition summary ---
    try:
        print("    - Generating cluster composition summary...")
        cluster_counts = adata.obs[primary_key].value_counts().sort_index()
        
        comp_rows = []
        for cluster_id in cluster_counts.index:
            mask = adata.obs[primary_key] == cluster_id
            n_cells = mask.sum()
            pct = n_cells / adata.n_obs * 100
            
            # Top expressed genes in this cluster (mean expression)
            cluster_expr = adata[mask].X
            if hasattr(cluster_expr, 'toarray'):
                cluster_expr = cluster_expr.toarray()
            gene_means = np.mean(cluster_expr, axis=0)
            top5_idx = np.argsort(gene_means)[::-1][:5]
            top5_genes = adata.var_names[top5_idx].tolist()
            top5_expr = gene_means[top5_idx].tolist()
            
            # Top marker genes from rank_genes_groups (if available)
            top_markers = []
            if 'rank_genes_groups' in adata.uns:
                try:
                    top_markers = adata.uns['rank_genes_groups']['names'][cluster_id][:3].tolist()
                except (KeyError, IndexError):
                    pass
            
            comp_rows.append({
                'cluster': cluster_id,
                'n_cells': n_cells,
                'pct_cells': round(pct, 1),
                'median_counts': np.median(adata.obs.loc[mask, 'total_counts']),
                'median_genes': np.median(adata.obs.loc[mask, 'n_genes_by_counts']),
                'top5_expressed': ', '.join(top5_genes),
                'top5_expr_values': ', '.join([f'{v:.1f}' for v in top5_expr]),
                'top3_markers': ', '.join(top_markers) if top_markers else 'N/A'
            })
        
        comp_df = pd.DataFrame(comp_rows)
        comp_csv = os.path.join(output_dir, f"{sample_tag}_step1_cluster_composition.csv")
        comp_df.to_csv(comp_csv, index=False)
        print(f"    - Saved cluster composition to: {comp_csv}")
        
        # Print summary to console
        for _, row in comp_df.iterrows():
            print(f"      Cluster {row['cluster']}: {row['n_cells']} cells ({row['pct_cells']}%) | "
                  f"markers: {row['top3_markers']}")
    except Exception as e:
        print(f"    [WARNING] Cluster composition summary failed: {e}")

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
    possible_keys = ['leiden_1.0', 'leiden_0.8', 'leiden', 'graph_clusters']

    for key in possible_keys:
        if key in adata.obs.keys():
            cluster_key = key
            print(f"    - Using cluster key '{cluster_key}' for enrichment analysis.")
            break

    if cluster_key:
        print(f"    - Computing neighborhood enrichment for '{cluster_key}'...")
        try:
            import seaborn as sns
            sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, seed=42)

            zscore_matrix = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']
            count_matrix = adata.uns[f'{cluster_key}_nhood_enrichment']['count']
            cluster_names = adata.obs[cluster_key].cat.categories.tolist()
            n_clusters = len(cluster_names)
            z_abs_max = max(abs(np.nanmin(zscore_matrix)), abs(np.nanmax(zscore_matrix)))

            # --- Plot 1: Clean z-score heatmap (seaborn, no annotation clutter) ---
            fig, ax = plt.subplots(figsize=(8, 7))
            sns.heatmap(
                zscore_matrix,
                xticklabels=cluster_names,
                yticklabels=cluster_names,
                cmap='coolwarm', center=0,
                vmin=-z_abs_max, vmax=z_abs_max,
                linewidths=0.8, linecolor='white',
                square=True,
                cbar_kws={'label': 'Z-score', 'shrink': 0.8},
                ax=ax,
            )
            ax.set_xlabel(cluster_key, fontsize=11)
            ax.set_ylabel(cluster_key, fontsize=11)
            ax.set_title(f"Neighborhood Enrichment Z-score ({cluster_key})", fontsize=12, pad=10)
            ax.tick_params(axis='both', labelsize=10)
            save_path = os.path.join(output_dir, f"{sample_tag}_step1_nhood_enrichment.png")
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    - Saved enrichment z-score plot to {save_path}")

            # --- Plot 2: Annotated z-score (supplementary, with values) ---
            fig, ax = plt.subplots(figsize=(9, 8))
            sns.heatmap(
                zscore_matrix,
                xticklabels=cluster_names,
                yticklabels=cluster_names,
                annot=True, fmt='.0f', annot_kws={'size': 9},
                cmap='coolwarm', center=0,
                vmin=-z_abs_max, vmax=z_abs_max,
                linewidths=0.8, linecolor='white',
                square=True,
                cbar_kws={'label': 'Z-score', 'shrink': 0.8},
                ax=ax,
            )
            ax.set_xlabel(cluster_key, fontsize=11)
            ax.set_ylabel(cluster_key, fontsize=11)
            ax.set_title(
                f"Neighborhood Enrichment Z-score ({cluster_key})\n"
                f"Diagonal = within-cluster | Off-diagonal = between-cluster",
                fontsize=11, pad=10,
            )
            ax.tick_params(axis='both', labelsize=10)
            save_path_annot = os.path.join(output_dir, f"{sample_tag}_step1_nhood_enrichment_annotated.png")
            plt.savefig(save_path_annot, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    - Saved annotated enrichment plot to {save_path_annot}")

            # --- Plot 3: Count heatmap (clean, no annotations — just colors) ---
            fig, ax = plt.subplots(figsize=(8, 7))
            sns.heatmap(
                count_matrix,
                xticklabels=cluster_names,
                yticklabels=cluster_names,
                cmap='YlGnBu',
                linewidths=0.8, linecolor='white',
                square=True,
                cbar_kws={'label': 'Neighbor count', 'shrink': 0.8},
                ax=ax,
            )
            ax.set_xlabel(cluster_key, fontsize=11)
            ax.set_ylabel(cluster_key, fontsize=11)
            ax.set_title(f"Neighborhood Enrichment Counts ({cluster_key})", fontsize=12, pad=10)
            ax.tick_params(axis='both', labelsize=10)
            save_path_count = os.path.join(output_dir, f"{sample_tag}_step1_nhood_enrichment_count.png")
            plt.savefig(save_path_count, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    - Saved enrichment count plot to {save_path_count}")

            colocalized_pairs = []
            for i in range(n_clusters):
                for j in range(i + 1, n_clusters):
                    z = zscore_matrix[i, j]
                    if z > 2.0:
                        colocalized_pairs.append((cluster_names[i], cluster_names[j], z))

            colocalized_pairs.sort(key=lambda x: x[2], reverse=True)

            if colocalized_pairs:
                print(f"    - Co-localized cluster pairs (z > 2.0):")
                for c1, c2, z in colocalized_pairs[:10]:  # top 10
                    print(f"      {c1} \u2194 {c2}: z={z:.1f}")
                if len(colocalized_pairs) > 10:
                    print(f"      ... and {len(colocalized_pairs) - 10} more pairs")
            else:
                print(f"    - No co-localized cluster pairs found (z > 2.0).")

            # Print diagonal z-scores (within-cluster enrichment)
            print(f"    - Diagonal z-scores (within-cluster spatial co-localization):")
            for i, name in enumerate(cluster_names):
                diag_z = zscore_matrix[i, i]
                print(f"      Cluster {name}: z={diag_z:.1f}")

            # Save co-localized pairs to text file
            coloc_path = os.path.join(output_dir, f"{sample_tag}_step1_colocalized_clusters.txt")
            with open(coloc_path, 'w') as f_coloc:
                f_coloc.write(f"# Neighborhood Enrichment Analysis ({cluster_key})\n")
                f_coloc.write(f"# Spatial neighbor graph: radius={radius}um\n")
                f_coloc.write(f"# Permutation test: 1000 permutations, seed=42\n\n")
                f_coloc.write("## Within-cluster enrichment (diagonal z-scores):\n")
                for i, name in enumerate(cluster_names):
                    f_coloc.write(f"  Cluster {name}: z={zscore_matrix[i, i]:.1f}\n")
                f_coloc.write(f"\n## Co-localized pairs (z > 2.0, off-diagonal):\n")
                for c1, c2, z in colocalized_pairs:
                    f_coloc.write(f"  Cluster {c1} \u2194 Cluster {c2}: z={z:.1f}\n")
            print(f"    - Saved enrichment report to {coloc_path}")

        except Exception as e:
            print(f"    [WARNING] Neighborhood enrichment failed: {e}")

        # --- 1-5d. Centrality scores ---
        print(f"    - Computing centrality scores...")
        sq.gr.centrality_scores(adata, cluster_key=cluster_key)

        plt.figure(figsize=(10, 5))
        sq.pl.centrality_scores(adata, cluster_key=cluster_key)
        plt.title("Centrality Scores")
        save_path_cent = os.path.join(output_dir, f"{sample_tag}_step1_centrality.png")
        plt.savefig(save_path_cent, bbox_inches='tight')
        plt.close()
        print(f"    - Saved centrality barplot to {save_path_cent}")

        # --- 1-5e. Spatial graph visualization (cropped region) ---
        try:
            print("    - Generating spatial graph visualization (cropped region)...")
            from scipy import sparse
            coords = adata.obsm['spatial']
            adj = adata.obsp['spatial_connectivities']
            if not sparse.issparse(adj):
                adj = sparse.csr_matrix(adj)

            # Pick a dense region: find the cell with the most neighbors, crop around it
            degree = np.array(adj.sum(axis=1)).flatten()
            center_idx = np.argmax(degree)
            cx, cy = coords[center_idx]
            crop_radius = 250  # 500x500 µm region

            # Select cells within the crop window
            in_crop = (
                (coords[:, 0] >= cx - crop_radius) & (coords[:, 0] <= cx + crop_radius) &
                (coords[:, 1] >= cy - crop_radius) & (coords[:, 1] <= cy + crop_radius)
            )
            crop_indices = np.where(in_crop)[0]

            if len(crop_indices) > 50:
                crop_coords = coords[crop_indices]
                crop_adj = adj[crop_indices][:, crop_indices]

                # Get cluster labels and colors
                cluster_labels = adata.obs[cluster_key].values[crop_indices]
                categories = adata.obs[cluster_key].cat.categories
                color_key = f'{cluster_key}_colors'
                if color_key in adata.uns:
                    palette = {cat: col for cat, col in zip(categories, adata.uns[color_key])}
                else:
                    cmap = plt.cm.get_cmap('tab20', len(categories))
                    palette = {cat: cmap(i) for i, cat in enumerate(categories)}
                node_colors = [palette.get(lbl, '#999999') for lbl in cluster_labels]

                fig, ax = plt.subplots(figsize=(10, 10))
                # Draw edges
                cx_arr, cy_arr = crop_adj.nonzero()
                for i, j in zip(cx_arr, cy_arr):
                    if i < j:  # avoid drawing each edge twice
                        ax.plot(
                            [crop_coords[i, 0], crop_coords[j, 0]],
                            [crop_coords[i, 1], crop_coords[j, 1]],
                            color='#cccccc', linewidth=0.3, alpha=0.4, zorder=1
                        )
                # Draw nodes
                ax.scatter(
                    crop_coords[:, 0], crop_coords[:, 1],
                    c=node_colors, s=15, edgecolors='k', linewidths=0.2, zorder=2
                )
                ax.set_title(
                    f"Spatial Neighbor Graph (r={radius}µm)\n"
                    f"Cropped {2*crop_radius}×{2*crop_radius}µm region | "
                    f"{len(crop_indices)} cells | Colored by {cluster_key}",
                    fontsize=11
                )
                ax.set_xlabel("x (µm)")
                ax.set_ylabel("y (µm)")
                ax.set_aspect('equal')
                ax.invert_yaxis()
                # Legend (max 20 categories)
                shown_cats = sorted(set(cluster_labels))[:20]
                handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=palette.get(c, '#999999'), markersize=6, label=c)
                           for c in shown_cats]
                ax.legend(handles=handles, title=cluster_key, loc='upper right',
                          fontsize=7, title_fontsize=8, ncol=max(1, len(shown_cats)//10))
                save_path_graph = os.path.join(output_dir, f"{sample_tag}_step1_spatial_graph.png")
                plt.savefig(save_path_graph, dpi=200, bbox_inches='tight')
                plt.close()
                print(f"    - Saved spatial graph plot to {save_path_graph}")
            else:
                print("    [WARNING] Crop region too sparse for graph visualization.")
        except Exception as e:
            print(f"    [WARNING] Spatial graph visualization failed: {e}")

        # --- 1-5f. Cluster network diagram (centrality as node size) ---
        try:
            print("    - Generating cluster network diagram...")
            import networkx as nx
            adj = adata.obsp['spatial_connectivities']
            labels = adata.obs[cluster_key].values
            categories = adata.obs[cluster_key].cat.categories.tolist()
            n_clust = len(categories)
            cat_to_idx = {c: i for i, c in enumerate(categories)}

            # Build cluster-level adjacency matrix (count inter-cluster edges)
            cluster_adj = np.zeros((n_clust, n_clust))
            rows, cols = adj.nonzero()
            for r, c in zip(rows, cols):
                ci, cj = cat_to_idx.get(labels[r]), cat_to_idx.get(labels[c])
                if ci is not None and cj is not None:
                    cluster_adj[ci, cj] += 1
            # Symmetrize
            cluster_adj = (cluster_adj + cluster_adj.T) / 2

            # Get centrality scores
            cent_key = f'{cluster_key}_centrality_scores'
            if cent_key in adata.uns:
                cent_df = adata.uns[cent_key]
                closeness_cent = cent_df['closeness_centrality'].values if 'closeness_centrality' in cent_df.columns else np.ones(n_clust)
            else:
                closeness_cent = np.ones(n_clust)

            # Get enrichment z-scores for edge coloring
            enrich_key = f'{cluster_key}_nhood_enrichment'
            has_enrich = enrich_key in adata.uns and 'zscore' in adata.uns[enrich_key]

            # Cell counts per cluster for node sizing
            cell_counts = adata.obs[cluster_key].value_counts()

            # Build networkx graph
            G = nx.Graph()
            for i, cat in enumerate(categories):
                G.add_node(cat, cell_count=cell_counts.get(cat, 0))
            # Add edges with z-score coloring
            max_weight = cluster_adj.max() if cluster_adj.max() > 0 else 1
            for i in range(n_clust):
                for j in range(i + 1, n_clust):
                    w = cluster_adj[i, j]
                    if w > max_weight * 0.02:  # edges > 2% of max
                        zscore = float(adata.uns[enrich_key]['zscore'][i, j]) if has_enrich else 0.0
                        G.add_edge(categories[i], categories[j], weight=w, zscore=zscore)

            # Layout — kamada_kawai gives better separation than spring
            pos = nx.kamada_kawai_layout(G)

            # Node sizes proportional to cell count
            counts_arr = np.array([G.nodes[n]['cell_count'] for n in G.nodes()])
            min_size, max_size = 400, 3000
            if counts_arr.max() > counts_arr.min():
                norm_counts = (counts_arr - counts_arr.min()) / (counts_arr.max() - counts_arr.min())
            else:
                norm_counts = np.ones(len(counts_arr)) * 0.5
            node_sizes = min_size + norm_counts * (max_size - min_size)

            # Node colors
            color_key = f'{cluster_key}_colors'
            if color_key in adata.uns:
                node_colors = list(adata.uns[color_key][:n_clust])
            else:
                cmap = plt.cm.get_cmap('tab20', n_clust)
                node_colors = [cmap(i) for i in range(n_clust)]

            # Edge widths and colors (red = co-localized, blue = segregated)
            edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
            max_ew = max(edge_weights) if edge_weights else 1
            edge_widths = [0.5 + 4.0 * (w / max_ew) for w in edge_weights]

            if has_enrich:
                edge_zscores = [G[u][v]['zscore'] for u, v in G.edges()]
                z_abs = max(abs(min(edge_zscores)), abs(max(edge_zscores))) if edge_zscores else 1
                edge_cmap = plt.cm.coolwarm
                edge_colors = [edge_cmap(0.5 + 0.5 * z / z_abs) if z_abs > 0 else '#888888' for z in edge_zscores]
            else:
                edge_colors = ['#888888'] * len(G.edges())

            fig, ax = plt.subplots(figsize=(12, 10))
            nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.6, edge_color=edge_colors)
            nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                                   node_color=node_colors, edgecolors='black', linewidths=1.2)
            # Labels with centrality annotation
            node_labels = {cat: f"{cat}\n({closeness_cent[i]:.2f})" for i, cat in enumerate(categories)}
            nx.draw_networkx_labels(G, pos, labels=node_labels, ax=ax, font_size=8, font_weight='bold')

            ax.set_title(
                f"Cluster Network Diagram ({cluster_key})\n"
                f"Node size = cell count | Label = closeness centrality\n"
                f"Edge color: red = co-localized, blue = segregated (enrichment z-score)\n"
                f"Graph built from spatial neighbors (r={radius}µm)",
                fontsize=11
            )
            ax.axis('off')

            # Add centrality summary box
            cent_text = "Closeness centrality:\n" + "\n".join(
                [f"  Cluster {cat}: {closeness_cent[i]:.3f}" for i, cat in enumerate(categories)]
            )
            ax.text(0.02, 0.02, cent_text, transform=ax.transAxes, fontsize=7,
                    verticalalignment='bottom', fontfamily='monospace',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

            save_path_net = os.path.join(output_dir, f"{sample_tag}_step1_cluster_network.png")
            plt.savefig(save_path_net, dpi=200, bbox_inches='tight')
            plt.close()
            print(f"    - Saved cluster network diagram to {save_path_net}")
        except Exception as e:
            print(f"    [WARNING] Cluster network diagram failed: {e}")

        # --- 1-5g. Spatial centrality heatmap (3-panel) ---
        try:
            print("    - Generating spatial centrality heatmap (3-panel)...")
            cent_key = f'{cluster_key}_centrality_scores'
            enrich_key = f'{cluster_key}_nhood_enrichment'
            if cent_key in adata.uns:
                cent_df = adata.uns[cent_key]
                metric_col = 'closeness_centrality' if 'closeness_centrality' in cent_df.columns else cent_df.columns[0]
                cent_map = dict(zip(cent_df.index.astype(str), cent_df[metric_col].values))

                # Map closeness centrality to each cell (fix: astype(str) to avoid Categorical error)
                cent_scores = adata.obs[cluster_key].astype(str).map(cent_map).fillna(0.0).values

                # Map within-cluster enrichment (diagonal z-scores) to each cell
                has_enrich = enrich_key in adata.uns and 'zscore' in adata.uns[enrich_key]
                if has_enrich:
                    zscore_data = adata.uns[enrich_key]['zscore']
                    diag_zscores = np.diag(zscore_data)
                    categories_list = adata.obs[cluster_key].cat.categories.tolist()
                    diag_map = dict(zip([str(c) for c in categories_list], diag_zscores))
                    enrich_scores = adata.obs[cluster_key].astype(str).map(diag_map).fillna(0.0).values
                else:
                    enrich_scores = np.zeros(len(adata))

                coords = adata.obsm['spatial']
                fig, axes = plt.subplots(1, 3, figsize=(28, 8))

                # Panel 1: Cluster identity on spatial coords
                categories = adata.obs[cluster_key].cat.categories
                color_key_name = f'{cluster_key}_colors'
                if color_key_name in adata.uns:
                    palette = {cat: col for cat, col in zip(categories, adata.uns[color_key_name])}
                else:
                    cmap_cat = plt.cm.get_cmap('tab20', len(categories))
                    palette = {cat: cmap_cat(i) for i, cat in enumerate(categories)}
                cell_colors = [palette.get(lbl, '#999999') for lbl in adata.obs[cluster_key].values]
                axes[0].scatter(coords[:, 0], coords[:, 1], c=cell_colors, s=0.3, alpha=0.5)
                axes[0].set_title(f"Spatial Map — {cluster_key}", fontsize=11)
                axes[0].set_xlabel("x (µm)")
                axes[0].set_ylabel("y (µm)")
                axes[0].set_aspect('equal')
                axes[0].invert_yaxis()

                # Panel 2: Within-cluster enrichment z-score
                if has_enrich:
                    z_abs = max(abs(np.nanmin(enrich_scores)), abs(np.nanmax(enrich_scores)))
                    im2 = axes[1].scatter(coords[:, 0], coords[:, 1], c=enrich_scores,
                                          cmap='coolwarm', s=0.3, alpha=0.5,
                                          vmin=-z_abs if z_abs > 0 else -1, vmax=z_abs if z_abs > 0 else 1)
                    axes[1].set_title("Within-cluster Enrichment Z-score\n(diagonal of nhood_enrichment)", fontsize=11)
                    plt.colorbar(im2, ax=axes[1], label='z-score', shrink=0.7)
                else:
                    axes[1].text(0.5, 0.5, "No enrichment data", transform=axes[1].transAxes,
                                 ha='center', va='center', fontsize=12)
                    axes[1].set_title("Within-cluster Enrichment Z-score\n(not available)", fontsize=11)
                axes[1].set_xlabel("x (µm)")
                axes[1].set_ylabel("y (µm)")
                axes[1].set_aspect('equal')
                axes[1].invert_yaxis()

                # Panel 3: Rescaled closeness centrality (percentile-based for compressed ranges)
                p5, p95 = np.percentile(cent_scores[cent_scores > 0], [5, 95]) if (cent_scores > 0).any() else (0, 1)
                if p95 > p5:
                    rescaled = np.clip((cent_scores - p5) / (p95 - p5), 0, 1)
                else:
                    rescaled = cent_scores
                im3 = axes[2].scatter(coords[:, 0], coords[:, 1], c=rescaled,
                                      cmap='YlOrRd', s=0.3, alpha=0.5, vmin=0, vmax=1)
                axes[2].set_title(
                    f"Closeness Centrality (rescaled)\n"
                    f"Raw range: [{cent_scores.min():.3f}, {cent_scores.max():.3f}] → percentile-normalized",
                    fontsize=11
                )
                axes[2].set_xlabel("x (µm)")
                axes[2].set_ylabel("y (µm)")
                axes[2].set_aspect('equal')
                axes[2].invert_yaxis()
                plt.colorbar(im3, ax=axes[2], label='rescaled centrality', shrink=0.7)

                plt.suptitle(
                    f"Cluster identity vs. Spatial centrality & Enrichment\n"
                    f"Centrality computed on spatial neighbor graph (r={radius}µm), NOT on UMAP",
                    fontsize=12, fontweight='bold', y=1.02
                )
                save_path_spatial = os.path.join(output_dir, f"{sample_tag}_step1_spatial_centrality.png")
                plt.savefig(save_path_spatial, dpi=150, bbox_inches='tight')
                plt.close()
                print(f"    - Saved spatial centrality heatmap (3-panel) to {save_path_spatial}")
            else:
                print("    [WARNING] No centrality scores found for spatial heatmap.")
        except Exception as e:
            print(f"    [WARNING] Spatial centrality heatmap failed: {e}")

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
        calculate_transcript_dispersion(
            adata,
            input_dir,
            sample_tag,
            spots=spots,
            marker_qc_cfg=expl_params.get("marker_qc", {}),
        )

        # --- 1-3. KS tests ---
        if spots is not None:
            spots_assigned, distances, metric_name = _compute_distances(adata, spots)
            if spots_assigned is not None and 'distance' in spots_assigned.columns:
                run_ks_tests(spots_assigned, input_dir, sample_tag, n_top_genes=20)
            else:
                print("    [INFO] Skipping KS tests - could not compute distances.")
    else:
        print("\n[Step 1-2] Dispersion analysis skipped (enable 'run_transcript_dispersion' in config)")

    # --- 1-3.5. MTRNR gene proportion visualization ---
    _plot_mtrnr_proportion(adata, input_dir, sample_tag)

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
