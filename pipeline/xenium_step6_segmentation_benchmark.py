# Step 6: Segmentation Benchmark (Ref: notebooks/5_segmentation_benchmark/)
# Benchmarks multiple segmentation methods (nuclei, Cellpose, expansion, Baysor)
# by concatenating results, preprocessing, and computing comparison metrics.
#
# Flow:
#   6-1. Run Baysor (optional)
#   6-2. Load segmentation results (nuclei/cellpose/expansion/baysor)
#   6-3. Concatenate & preprocess (normalize/PCA/Leiden)
#   6-4. Annotation transfer (majority voting from reference)
#   6-5. Benchmark metrics (n_cells, median_reads, assigned_prop)
#   6-6. Visualizations
#     6-6a. UMAP per method
#     6-6b. Spatial scatter
#     6-6c. Cell type barplot
#     6-6d. Counts violin

import os
import logging
import shutil
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from scipy.sparse import issparse
from pipeline.benchmark_utils import metrics

logger = logging.getLogger("Step6_Benchmark")
logging.basicConfig(level=logging.INFO)


# --- 6-1. Baysor data preparation (helper) ---

def prep_xenium_data_for_baysor(xenium_dir, out_dir, crop=False, coords=None):
    """Format Xenium transcripts into Baysor-compatible CSV (gene, x, y)."""
    out_dir = Path(out_dir)
    xenium_dir = Path(xenium_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Preparing Xenium data for Baysor in {out_dir}...")

    transcripts_path = xenium_dir / "transcripts.csv"
    if not transcripts_path.exists():
        transcripts_path = xenium_dir / "transcripts.parquet"
        if transcripts_path.exists():
            spots = pd.read_parquet(transcripts_path)
        else:
            logger.error(f"Transcripts file not found in {xenium_dir}")
            return False
    else:
        spots = pd.read_csv(transcripts_path)

    logger.info(f"Loaded {len(spots)} transcripts.")

    col_map = {
        'feature_name': 'gene',
        'x_location': 'x',
        'y_location': 'y',
        'z_location': 'z'
    }

    spots_baysor = spots.rename(columns=col_map)

    req_cols = ['gene', 'x', 'y']
    if not all(col in spots_baysor.columns for col in req_cols):
        logger.error(f"Missing required columns for Baysor. Found: {spots_baysor.columns}")
        col_map_alt = {'feature_name': 'gene', 'x': 'x', 'y': 'y'}
        spots_baysor = spots.rename(columns=col_map_alt)
        if not all(col in spots_baysor.columns for col in req_cols):
            return False

    spots_file = out_dir / "spots.csv"
    spots_baysor[req_cols].to_csv(spots_file, index=False)
    logger.info(f"Saved spots to {spots_file}")

    return True


# --- 6-1. Baysor execution ---

def run_baysor(config, xenium_input_dir):
    """Execute Baysor segmentation. Controlled by config["baysor"]:
    enabled, executable_path, dry_run, params, prior_segmentation_tif."""
    logger.info("--- Starting Baysor Segmentation Phase ---")

    baysor_conf = config.get("baysor", {})
    if not baysor_conf.get("enabled", False):
        logger.info("Baysor disabled in config.")
        return None

    executable = baysor_conf.get("executable_path", "baysor")
    params = baysor_conf.get("params", {})
    dry_run = baysor_conf.get("dry_run", False)

    step6_out = config.get("output_dir", "xenium-output/step6_benchmark")
    baysor_out_dir = os.path.join(step6_out, "baysor_run")
    os.makedirs(baysor_out_dir, exist_ok=True)

    prep_dir = os.path.join(baysor_out_dir, "prep")
    if not prep_xenium_data_for_baysor(xenium_input_dir, prep_dir):
        logger.error("Failed to prepare data for Baysor.")
        return None

    cmd = [executable, "run"]

    if "scale" in params:
        cmd.extend(["-s", str(params["scale"])])
    if "min_molecules_per_cell" in params:
        cmd.extend(["-m", str(params["min_molecules_per_cell"])])
    if "prior_segmentation_confidence" in params:
        cmd.extend(["--prior-segmentation-confidence", str(params["prior_segmentation_confidence"])])

    segmentation_csv = os.path.join(baysor_out_dir, "segmentation.csv")
    cmd.extend(["-o", segmentation_csv])

    spots_file = os.path.join(prep_dir, "spots.csv")
    cmd.append(spots_file)

    # Prior segmentation TIF (from Step 3 or legacy input_advanced)
    prior_seg_path = baysor_conf.get("prior_segmentation_tif")
    if not prior_seg_path:
        val = config.get("benchmark", {}).get("input_advanced")
        if str(val).endswith('.tif'):
            prior_seg_path = val

    if prior_seg_path and str(prior_seg_path).endswith('.tif') and os.path.exists(prior_seg_path):
        logger.info(f"Using prior segmentation TIF: {prior_seg_path}")
        cmd.append(prior_seg_path)
    else:
        logger.warning(
            f"Prior segmentation TIF not valid or missing (Path: {prior_seg_path}). "
            "Running Baysor without prior."
        )

    logger.info(f"Running Baysor command: {' '.join(cmd)}")

    # --- Dry-run gate ---
    if dry_run:
        logger.warning(
            "dry_run=True: Baysor execution is SIMULATED. "
            "A minimal empty CSV will be written instead of real output."
        )
        if not os.path.exists(segmentation_csv):
            with open(segmentation_csv, "w") as f:
                f.write("transcript_id,cell\n")
        return segmentation_csv

    # --- Real execution ---
    if not shutil.which(executable):
        logger.error(
            f"Baysor binary not found at '{executable}'. "
            "Install Baysor or set baysor.executable_path in config. "
            "To skip execution, set baysor.dry_run=True."
        )
        return None

    try:
        subprocess.run(cmd, check=True)
        logger.info("Baysor executed successfully.")
        return segmentation_csv
    except subprocess.CalledProcessError as e:
        logger.error(f"Baysor execution failed with return code {e.returncode}: {e}")
        return None
    except FileNotFoundError:
        logger.error(
            f"Baysor binary '{executable}' could not be launched. "
            "Ensure it is installed and on PATH."
        )
        return None
    except Exception as e:
        logger.error(f"Unexpected error running Baysor: {e}")
        return None


# --- 6-2. Transcript aggregation (helper) ---

def load_transcripts_as_adata(csv_path, cell_col='cell', gene_col='gene',
                              current_method_name='unknown', ref_adata_path=None):
    """Load a transcripts CSV and aggregate into AnnData (cells x genes count matrix)."""
    logger.info(f"Loading transcripts from {csv_path} for method '{current_method_name}'...")
    try:
        df = pd.read_csv(csv_path)

        # Resolve column names (multiple naming conventions)
        cols = df.columns
        actual_cell_col = None
        if cell_col in cols:
            actual_cell_col = cell_col
        elif 'cell_id' in cols:
            actual_cell_col = 'cell_id'
        elif 'cell' in cols:
            actual_cell_col = 'cell'

        actual_gene_col = None
        if gene_col in cols:
            actual_gene_col = gene_col
        elif 'feature_name' in cols:
            actual_gene_col = 'feature_name'
        elif 'gene' in cols:
            actual_gene_col = 'gene'

        if not actual_cell_col or not actual_gene_col:
            logger.error(f"Missing required columns (cell/gene) in {csv_path}. Found: {list(cols)}")
            return None

        # Filter unassigned transcripts
        if pd.api.types.is_numeric_dtype(df[actual_cell_col]):
            df_assigned = df[df[actual_cell_col] > 0].copy()
        else:
            mask = (
                (df[actual_cell_col] != '0')
                & (df[actual_cell_col].astype(str).str.lower() != 'unassigned')
                & (df[actual_cell_col].astype(str).str.lower() != 'noise')
            )
            df_assigned = df[mask].copy()

        # Aggregate via crosstab (faster than loop + value_counts)
        counts = pd.crosstab(df_assigned[actual_cell_col], df_assigned[actual_gene_col])
        adata = ad.AnnData(counts)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        return adata

    except Exception as e:
        logger.error(f"Error aggregating csv {csv_path}: {e}")
        return None


# --- 6-4. Annotation transfer ---

def annotate_by_majority_voting(adata_target, adata_reference,
                                cluster_key='leiden', ref_key='celltype',
                                n_neighbors=15):
    """Transfer cell type labels from reference to target via kNN majority voting.
    Adds 'celltype_majority' (per-cell) and 'celltype_cluster' (per-cluster consensus)."""
    from sklearn.neighbors import NearestNeighbors

    if 'X_pca' not in adata_target.obsm or 'X_pca' not in adata_reference.obsm:
        logger.error(
            "PCA embeddings (X_pca) must be present in both target and "
            "reference AnnData objects before annotation transfer."
        )
        return adata_target

    if ref_key not in adata_reference.obs.columns:
        logger.error(f"Reference AnnData is missing the '{ref_key}' column.")
        return adata_target

    shared_genes = adata_target.var_names.intersection(adata_reference.var_names)
    if len(shared_genes) == 0:
        logger.error("No shared genes between target and reference.")
        return adata_target

    logger.info(
        f"Annotation transfer: {len(shared_genes)} shared genes, "
        f"using {n_neighbors} neighbours."
    )

    # Fit kNN on reference PCA, query with target PCA
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(adata_reference.obsm['X_pca'])
    distances, indices = nn.kneighbors(adata_target.obsm['X_pca'])

    # Per-cell majority vote
    ref_labels = adata_reference.obs[ref_key].values
    cell_labels = []
    for idx_row in indices:
        neighbour_labels = ref_labels[idx_row]
        values, counts = np.unique(neighbour_labels, return_counts=True)
        cell_labels.append(values[np.argmax(counts)])

    adata_target.obs['celltype_majority'] = cell_labels

    # Per-cluster consensus
    cluster_labels = []
    for cluster_id in adata_target.obs[cluster_key].unique():
        mask = adata_target.obs[cluster_key] == cluster_id
        cluster_celltypes = adata_target.obs.loc[mask, 'celltype_majority']
        values, counts = np.unique(cluster_celltypes, return_counts=True)
        winner = values[np.argmax(counts)]
        cluster_labels.append((cluster_id, winner))

    cluster_map = dict(cluster_labels)
    adata_target.obs['celltype_cluster'] = (
        adata_target.obs[cluster_key].map(cluster_map).astype(str)
    )

    logger.info("Annotation transfer complete. Added 'celltype_majority' and 'celltype_cluster'.")
    return adata_target


# --- 6-3. Preprocessing ---

def preprocess_benchmark(adata, config):
    """Preprocess concatenated benchmark AnnData (normalize, PCA, Leiden).
    Parameters from config["benchmark"]["preprocessing"], defaults match xb/preprocessing.py."""
    pp = config.get("benchmark", {}).get("preprocessing", {})

    target_sum = pp.get("target_sum", 100)
    min_counts = pp.get("min_counts", 40)
    min_genes = pp.get("min_genes", 15)
    n_neighbors = pp.get("n_neighbors", 15)
    n_pcs = pp.get("n_pcs", 0)          # 0 = use ALL PCs
    umap_min_dist = pp.get("umap_min_dist", 0.1)
    resolution = pp.get("resolution", 1.0)
    scale = pp.get("scale", False)

    logger.info(
        f"Preprocessing params: target_sum={target_sum}, min_counts={min_counts}, "
        f"min_genes={min_genes}, n_neighbors={n_neighbors}, n_pcs={n_pcs}, "
        f"umap_min_dist={umap_min_dist}, resolution={resolution}, scale={scale}"
    )

    adata.layers['raw'] = adata.X.copy()

    # QC filtering
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata.raw = adata

    # Re-store raw layer after filtering (cell count may have changed)
    adata.layers['raw'] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    if scale:
        sc.pp.scale(adata)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=umap_min_dist)

    # Notebook labels this Louvain but actually uses Leiden
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

    return adata


# --- 6-6. Visualizations ---

# --- 6-6a. UMAP per method ---
def _save_umap(adata, output_dir):
    """UMAP coloured by segmentation method and Leiden cluster."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    sc.pl.umap(adata, color='segmentation', ax=axes[0], show=False,
               title='Segmentation method', frameon=False)
    sc.pl.umap(adata, color='leiden', ax=axes[1], show=False,
               title='Leiden clusters', frameon=False)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "umap_benchmark.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved UMAP plot to {out_path}")


# --- 6-6b. Spatial scatter ---
def _save_spatial_map(adata, output_dir):
    """Spatial scatter plot coloured by segmentation method.
    Requires x/y coordinates in adata.obs; skipped if not available."""
    x_col, y_col = None, None
    for xc, yc in [('x_centroid', 'y_centroid'),
                    ('x_location', 'y_location'),
                    ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        logger.warning(
            "Spatial coordinates not found in adata.obs. "
            "Skipping spatial segmentation map."
        )
        return

    methods = adata.obs['segmentation'].unique()
    n_methods = len(methods)
    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 6))
    if n_methods == 1:
        axes = [axes]

    for ax, method in zip(axes, methods):
        sub = adata[adata.obs['segmentation'] == method]
        ax.scatter(
            sub.obs[x_col].values, sub.obs[y_col].values,
            s=0.3, alpha=0.5, rasterized=True
        )
        ax.set_title(method)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_aspect('equal')
        ax.invert_yaxis()

    fig.suptitle("Spatial map by segmentation method", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "spatial_segmentation_map.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved spatial segmentation map to {out_path}")


# --- 6-6c. Cell type barplot ---
def _save_celltype_barplot(adata, output_dir, ct_key='celltype_majority'):
    """Bar plot of cell type proportions across segmentation methods."""
    if ct_key not in adata.obs.columns:
        logger.warning(
            f"'{ct_key}' not in adata.obs -- skipping cell type frequency barplot."
        )
        return

    ct_counts = (
        adata.obs
        .groupby(['segmentation', ct_key])
        .size()
        .reset_index(name='count')
    )
    totals = ct_counts.groupby('segmentation')['count'].transform('sum')
    ct_counts['proportion'] = ct_counts['count'] / totals

    methods = ct_counts['segmentation'].unique()
    celltypes = ct_counts[ct_key].unique()
    n_ct = len(celltypes)

    x = np.arange(n_ct)
    width = 0.8 / len(methods)

    fig, ax = plt.subplots(figsize=(max(10, n_ct * 0.8), 6))
    for i, method in enumerate(methods):
        sub = ct_counts[ct_counts['segmentation'] == method].set_index(ct_key)
        heights = sub.reindex(celltypes, fill_value=0)['proportion']
        ax.bar(x + i * width, heights.values, width, label=method)

    ax.set_xticks(x + width * (len(methods) - 1) / 2)
    ax.set_xticklabels(celltypes, rotation=45, ha='right')
    ax.set_ylabel("Proportion")
    ax.set_title("Cell type frequency by segmentation method")
    ax.legend(title="Method")
    fig.tight_layout()
    out_path = os.path.join(output_dir, "celltype_frequency_barplot.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved cell type frequency barplot to {out_path}")


# --- 6-6d. Counts violin ---
def _save_counts_violin(adata, output_dir):
    """Violin plot of total counts per cell grouped by segmentation method."""
    if 'raw' not in adata.layers:
        logger.warning("adata.layers['raw'] not found -- skipping counts violin plot.")
        return

    raw_X = adata.layers['raw']
    if issparse(raw_X):
        total_counts = np.asarray(raw_X.sum(axis=1)).flatten()
    else:
        total_counts = np.sum(raw_X, axis=1)

    plot_df = pd.DataFrame({
        'total_counts': total_counts,
        'segmentation': adata.obs['segmentation'].values
    })

    methods = plot_df['segmentation'].unique()
    fig, ax = plt.subplots(figsize=(max(6, len(methods) * 2), 6))

    parts = ax.violinplot(
        [plot_df.loc[plot_df['segmentation'] == m, 'total_counts'].values
         for m in methods],
        positions=range(len(methods)),
        showmeans=True, showmedians=True
    )

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods)
    ax.set_ylabel("Total counts per cell")
    ax.set_title("Counts distribution by segmentation method")
    fig.tight_layout()
    out_path = os.path.join(output_dir, "counts_violin_by_method.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved counts violin plot to {out_path}")


# --- Entry Point ---

def run_step6(config):
    """Step 6 orchestrator: Baysor -> load methods -> concat -> preprocess -> metrics -> plots."""
    logger.info("Starting Step 6: Segmentation Benchmark")

    output_dir = config.get("output_dir", "xenium-output/benchmark")
    os.makedirs(output_dir, exist_ok=True)

    # --- 6-1. Run Baysor (optional) ---
    baysor_result_csv = None
    if config.get("baysor", {}).get("enabled", False):
        xenium_input_dir = config.get("input_dir")
        if xenium_input_dir:
            baysor_result_csv = run_baysor(config, xenium_input_dir)

    # --- 6-2. Load segmentation results ---
    logger.info("--- Starting Benchmarking Phase ---")

    input_nuclei = config.get("benchmark", {}).get("input_nuclei")
    input_advanced = config.get("benchmark", {}).get("input_advanced")
    input_expansion = config.get("benchmark", {}).get("input_expansion")

    adata_list = []

    if input_nuclei and os.path.exists(input_nuclei):
        logger.info(f"Loading Nuclei Data: {input_nuclei}")
        nuclei_adata = sc.read(input_nuclei)
        nuclei_adata.obs['segmentation'] = 'nuclei'
        adata_list.append(nuclei_adata)

    if input_advanced and os.path.exists(input_advanced) and str(input_advanced).endswith('.h5ad'):
        logger.info(f"Loading Cellpose Data: {input_advanced}")
        adv_adata = sc.read(input_advanced)
        adv_adata.obs['segmentation'] = 'cellpose'
        adata_list.append(adv_adata)

    if input_expansion and os.path.exists(input_expansion):
        logger.info(f"Loading Expansion Data (Step 5): {input_expansion}")
        exp_adata = load_transcripts_as_adata(
            input_expansion, current_method_name='expansion'
        )
        if exp_adata is not None:
            exp_adata.obs['segmentation'] = 'expansion'
            adata_list.append(exp_adata)

    if baysor_result_csv:
        baysor_adata = load_transcripts_as_adata(
            baysor_result_csv, cell_col='cell', gene_col='gene',
            current_method_name='baysor'
        )
        if baysor_adata is not None:
            baysor_adata.obs['segmentation'] = 'baysor'
            adata_list.append(baysor_adata)

    if not adata_list:
        logger.error("No valid input datasets found for benchmarking.")
        return

    # --- 6-3. Concatenate & preprocess ---
    logger.info(f"Concatenating {len(adata_list)} datasets for comparison...")
    adata = sc.concat(adata_list)

    logger.info("Preprocessing and Clustering...")
    adata = preprocess_benchmark(adata, config)

    # --- 6-4. Annotation transfer ---
    ref_adata_path = config.get("benchmark", {}).get("reference_adata")
    if ref_adata_path and os.path.exists(ref_adata_path):
        logger.info(f"Loading reference AnnData for annotation transfer: {ref_adata_path}")
        adata_ref = sc.read(ref_adata_path)

        if 'X_pca' not in adata_ref.obsm:
            logger.info("Computing PCA on reference AnnData for annotation transfer...")
            sc.pp.pca(adata_ref)

        ref_key = config.get("benchmark", {}).get("ref_celltype_key", "celltype")
        pp = config.get("benchmark", {}).get("preprocessing", {})
        n_neighbors = pp.get("n_neighbors", 15)

        adata = annotate_by_majority_voting(
            adata, adata_ref,
            cluster_key='leiden', ref_key=ref_key,
            n_neighbors=n_neighbors
        )
    else:
        logger.info(
            "No reference AnnData provided (benchmark.reference_adata). "
            "Skipping annotation transfer."
        )

    combined_output = os.path.join(output_dir, "benchmark_combined.h5ad")
    adata.write_h5ad(combined_output)
    logger.info(f"Saved combined benchmark AnnData to {combined_output}")

    # --- 6-5. Benchmark metrics ---
    logger.info("Calculating Benchmark Metrics...")
    results = {}

    for seg_method in adata.obs['segmentation'].unique():
        subset = adata[adata.obs['segmentation'] == seg_method]
        n_cells = subset.shape[0]
        results[f'{seg_method}_n_cells'] = n_cells

        if 'raw' in subset.layers:
            results[f'{seg_method}_median_reads'] = metrics.median_reads_cells(subset)
            results[f'{seg_method}_median_genes'] = metrics.median_genes_cells(subset)
        else:
            logger.warning(
                f"layers['raw'] missing for '{seg_method}' -- "
                "skipping median_reads / median_genes metrics."
            )

        if 'spots' in subset.uns:
            results[f'{seg_method}_assigned_prop'] = metrics.proportion_of_assigned_reads(subset)

    metrics_df = pd.DataFrame([results])
    metrics_file = os.path.join(output_dir, "benchmark_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    logger.info(f"Saved metrics to {metrics_file}")

    # --- 6-6. Visualizations ---
    logger.info("Generating visualizations...")
    _save_umap(adata, output_dir)
    _save_spatial_map(adata, output_dir)
    _save_celltype_barplot(adata, output_dir, ct_key='celltype_majority')
    _save_counts_violin(adata, output_dir)

    logger.info("Step 6 Completed Successfully.")


if __name__ == "__main__":
    pass
