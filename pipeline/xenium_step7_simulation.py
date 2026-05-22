# Step 7: Simulation & Benchmarking (Ref: notebooks/6_simulating_preprocessing/)
# Simulates Xenium data from scRNAseq reference, then benchmarks preprocessing
# parameter combinations using clustering quality metrics (ARI, NMI, FMI, VI).
#
# Flow:
#   7-1. Reference acquisition (CellxGene Census)
#   7-2. Simulation
#     7-2a. Subsample + HVG filter
#     7-2b. Rank marker genes (Wilcoxon)
#     7-2c. Standard simulation (base noise/misseg)
#     7-2d. High-noise simulation (2× noise/misseg)
#   7-3. Benchmarking
#     7-3a. Build preprocessing grid
#     7-3b. Grid search (normalize/log1p/hvg/scale/PCA/Leiden)
#     7-3c. Clustering metrics (NMI, ARI, FMI, VI)
#     7-3d. Benchmark plots (bar, heatmap, importance, boxplot)
#   7-4. Perturbation analysis (single-param sensitivity)

import os
import json
import logging
import itertools
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import random
import math
from tqdm import tqdm
from scipy.spatial import ConvexHull
from sklearn.metrics import (
    mutual_info_score,
    silhouette_score,
    fowlkes_mallows_score,
    normalized_mutual_info_score,
    adjusted_rand_score,
)
import alphashape
from shapely.geometry import Point, Polygon
import warnings

try:
    import cellxgene_census
except ImportError:
    cellxgene_census = None

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
warnings.filterwarnings('ignore')

MAX_GRID_COMBINATIONS = 50

# Paper/notebook default preprocessing params (confirmed as best by simulation)
DEFAULT_PARAMS = {
    'normalize': True, 'target_sum': 100, 'log1p': True,
    'scale': True, 'hvg': False, 'n_neighbors': 16,
    'n_pcs': 0, 'resolution': 1.0,
}

# Paper's simulation confirmed these as optimal
PAPER_BEST_PARAMS = dict(DEFAULT_PARAMS)


# --- 7-2. Simulation helpers (inlined from xb/simulating.py) ---

def missegmentation_simulation(adata_sc_sub, missegmentation_percentage=0.1):
    """Simulate missegmentation by mixing expression profiles between random cell pairs."""
    exp = adata_sc_sub.to_df()
    if missegmentation_percentage > 0:
        cells_affected = int(exp.shape[0] * (missegmentation_percentage / 100))
        indices = list(exp.index)
        sources = random.sample(indices, cells_affected)
        targets = random.sample(indices, cells_affected)

        for s, t in zip(sources, targets):
             factor = random.random()
             # xb uses discrete importance levels: 0.0, 0.1, ... 0.9
             importance = random.randint(0, 9) * 0.1
             if importance > 0:
                exp.loc[t] += exp.loc[s] * importance

    adata_sc_sub.X = np.array(exp.values)
    return adata_sc_sub

def noise_adder(adata_sc, percentage_of_noise=0.1):
    """Add random +/-1 noise to count matrix at the given percentage of total counts."""
    total_counts = np.sum(adata_sc.X)
    noise_events = int(total_counts * (percentage_of_noise / 100))

    if noise_events > 0:
        rows = np.random.randint(0, adata_sc.X.shape[0], noise_events)
        cols = np.random.randint(0, adata_sc.X.shape[1], noise_events)
        changes = np.random.choice([-1, 1], noise_events)

        if hasattr(adata_sc.X, "tocoo"):
             adata_sc.X = adata_sc.X.toarray()

        # np.add.at handles duplicate indices correctly (unlike +=)
        np.add.at(adata_sc.X, (rows, cols), changes)
        adata_sc.X[adata_sc.X < 0] = 0

    return adata_sc

def subset_of_single_cell(adata_sc_sub, markers, random_markers_percentage=0,
                          reads_x_cell=None, number_of_markers=200,
                          percentage_of_noise=0.1, ms_percentage=0.1):
    """Transform scRNAseq data into simulated spatial data with marker selection,
    missegmentation, noise, and optional downsampling."""

    # Collect unique markers from ranked genes DataFrame
    mk = []
    unique_markers = np.unique(markers.values.flatten())
    mk = list(unique_markers)

    # Add random non-marker genes
    all_genes = list(adata_sc_sub.var.index)
    if random_markers_percentage > 0:
         n_random = int(len(mk) * (random_markers_percentage / 100))
         random_genes = random.sample(all_genes, n_random)
         mk.extend(random_genes)

    if len(mk) > number_of_markers:
        mk = random.sample(mk, number_of_markers)

    adata_sc = adata_sc_sub[:, adata_sc_sub.var.index.isin(mk)].copy()

    if adata_sc.n_vars == 0:
        print("    [Sim] Warning: No marker genes found in adata. Check naming convention.")
        return adata_sc_sub

    adata_sc = missegmentation_simulation(adata_sc, missegmentation_percentage=ms_percentage)
    adata_sc = noise_adder(adata_sc, percentage_of_noise=percentage_of_noise)

    if reads_x_cell is not None:
         if np.median(np.sum(adata_sc.X, axis=1)) > reads_x_cell:
             sc.pp.downsample_counts(adata_sc, counts_per_cell=reads_x_cell)

    return adata_sc


# --- 7-3c. Clustering metrics helpers ---

def entropy(clustering):
    """Shannon entropy of a clustering assignment."""
    _, counts = np.unique(clustering, return_counts=True)
    proportions = counts / len(clustering)
    return -np.sum(proportions * np.log(proportions))

def compute_vi(ground_truth, predicted):
    """Variation of Information."""
    mi = mutual_info_score(ground_truth, predicted)
    h_gt = entropy(ground_truth)
    h_pred = entropy(predicted)
    return h_gt + h_pred - 2 * mi

def compute_nmi(ground_truth, predicted):
    """Normalized Mutual Information."""
    return normalized_mutual_info_score(ground_truth, predicted)

def compute_fmi(ground_truth, predicted):
    """Fowlkes-Mallows Index."""
    return fowlkes_mallows_score(ground_truth, predicted)


# --- 7-3a. Build preprocessing grid ---

def build_preprocessing_grid(sim_config):
    """Build parameter grid for preprocessing benchmark. Sub-samples to
    MAX_GRID_COMBINATIONS using stratified sampling if the full grid is too large."""

    grid_spec = sim_config.get("preprocessing_grid", {})

    grid = {
        'normalize':    grid_spec.get('normalize',   [True, False]),
        'target_sum':   grid_spec.get('target_sum',  [100, 1000, 10000]),
        'log1p':        grid_spec.get('log1p',       [True]),
        'scale':        grid_spec.get('scale',       [False, True]),
        'hvg':          grid_spec.get('hvg',         [False]),
        'n_neighbors':  grid_spec.get('n_neighbors', [5, 15, 30]),
        'n_pcs':        grid_spec.get('n_pcs',       [0, 10, 30]),
        'resolution':   grid_spec.get('resolution',  [0.5, 1.0, 1.5]),
    }

    keys = list(grid.keys())
    values = [grid[k] for k in keys]
    all_combos = [dict(zip(keys, combo)) for combo in itertools.product(*values)]

    logger.info(f"Full grid has {len(all_combos)} combinations (max allowed: {MAX_GRID_COMBINATIONS}).")

    if len(all_combos) <= MAX_GRID_COMBINATIONS:
        return all_combos

    # Seed subset with one sample per unique value per dimension for coverage
    rng = np.random.default_rng(42)
    selected_indices = set()

    for key_idx, key in enumerate(keys):
        for val in grid[key]:
            matching = [i for i, c in enumerate(all_combos) if c[key] == val]
            if matching:
                selected_indices.add(rng.choice(matching))

    # Fill remainder randomly
    remaining = list(set(range(len(all_combos))) - selected_indices)
    n_extra = MAX_GRID_COMBINATIONS - len(selected_indices)
    if n_extra > 0 and remaining:
        extra = rng.choice(remaining, size=min(n_extra, len(remaining)), replace=False)
        selected_indices.update(extra.tolist())

    selected = [all_combos[i] for i in sorted(selected_indices)]
    logger.info(f"Sub-sampled to {len(selected)} representative combinations.")
    return selected


def _params_to_run_name(params):
    """Create a short unique label for a parameter combination."""
    norm_tag = "Norm" if params['normalize'] else "Raw"
    ts_tag = f"TS{params['target_sum']}" if params['normalize'] else "TS-"
    scale_tag = "Sc" if params['scale'] else "NoSc"
    hvg_tag = "HVG" if params['hvg'] else "NoHVG"
    return (
        f"{norm_tag}_{ts_tag}_{scale_tag}_{hvg_tag}"
        f"_N{params['n_neighbors']}_P{params['n_pcs']}_R{params['resolution']}"
    )


def _apply_preprocessing(adata, params, compute_umap=False):
    """Apply one preprocessing workflow + Leiden clustering. Returns modified copy.

    When n_pcs=0, PCA is still computed (matching notebook behaviour) and all
    principal components are passed to ``sc.pp.neighbors`` via ``n_pcs=None``.
    """
    ad = adata.copy()

    if hasattr(ad.X, "toarray"):
        ad.X = ad.X.toarray()

    if params['normalize']:
        sc.pp.normalize_total(ad, target_sum=params['target_sum'])

    if params['log1p']:
        sc.pp.log1p(ad)

    if params['hvg']:
        if ad.n_vars > 2000:
            sc.pp.highly_variable_genes(ad, n_top_genes=min(2000, ad.n_vars))
            ad = ad[:, ad.var['highly_variable']].copy()

    if params['scale']:
        sc.pp.scale(ad, max_value=10)

    n_pcs = params['n_pcs']

    # Always run PCA when there are enough features
    n_comps = min(max(n_pcs, 50), ad.n_vars - 1, ad.n_obs - 1)
    if n_comps >= 2:
        sc.pp.pca(ad, n_comps=n_comps)
        # n_pcs=0 means "use all PCs" (pass None to neighbors)
        neighbor_pcs = n_pcs if n_pcs > 0 else None
        sc.pp.neighbors(ad, n_neighbors=params['n_neighbors'], n_pcs=neighbor_pcs)
    else:
        sc.pp.neighbors(ad, n_neighbors=params['n_neighbors'], use_rep='X')

    sc.tl.leiden(ad, resolution=params['resolution'], key_added='leiden')

    if compute_umap:
        sc.tl.umap(ad)

    return ad


# --- 7-1. Data Acquisition ---

def run_step7_1_acquisition(config):
    """Download scRNAseq reference from CellxGene Census (or load existing)."""
    print("\n[Step 7-1] Acquiring Reference Data (Ref: Notebook 6_1)...")
    output_dir = config["output_dir"]
    sim_config = config.get("simulation", {})
    tissue = sim_config.get("census_tissue", "brain")
    organism = sim_config.get("census_organism", "mus_musculus")

    sc_file = os.path.join(output_dir, f"census_{tissue}_reference.h5ad")

    if os.path.exists(sc_file):
        print(f"    - Found existing reference: {sc_file}")
        return sc_file

    if cellxgene_census is None:
        print("    - Error: 'cellxgene_census' not installed. Cannot download.")
        print("    - Please install or provide 'census_brain_reference.h5ad' manually.")
        return None

    print("    - Downloading from CellxGene Census (this may take time)...")
    try:
        census = cellxgene_census.open_soma()
        adata_ref = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            obs_value_filter=f"tissue_general == '{tissue}'",
            column_names={"obs": ["cell_type", "tissue", "tissue_general", "dataset_id"]}
        )
        adata_ref.write_h5ad(sc_file)
        print(f"    - Downloaded and saved {adata_ref.shape} cells.")
        return sc_file
    except Exception as e:
        print(f"    - Census download failed: {e}")
        return None


# --- 7-2. Simulation ---

def run_step7_2_simulation(config, ref_file):
    """Generate synthetic Xenium datasets with missegmentation and noise."""
    print("\n[Step 7-2] Generating Simulated Datasets (Ref: Notebook 6_2)...")
    if not ref_file or not os.path.exists(ref_file):
        print("    - Error: valid reference file required.")
        return []

    output_dir = config["output_dir"]
    sample_tag = config.get("sample_tag", "sample")
    sim_config = config.get("simulation", {})

    n_markers = sim_config.get("n_markers", 200)
    noise_pct = sim_config.get("noise_percentage", 10.0)
    ms_pct = sim_config.get("missegmentation_percentage", 10.0)
    ref_key = sim_config.get("reference_key", "cell_type")
    max_cells = sim_config.get("max_cells", 20000)
    min_celltypes = sim_config.get("min_celltypes", 2)
    max_celltypes = sim_config.get("max_celltypes", 30)

    adata_ref = sc.read_h5ad(ref_file)
    print(f"    - Loaded Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    # --- 7-2a. Subsample + HVG filter ---
    if adata_ref.n_obs > max_cells:
        print(f"    - Sub-sampling from {adata_ref.n_obs} to {max_cells} cells...")
        sc.pp.subsample(adata_ref, n_obs=max_cells)
        print(f"    - After sub-sampling: {adata_ref.n_obs} cells")

    # Filter to top HVGs for efficiency
    if adata_ref.n_vars > 5000:
        print("    - Filtering to top 5000 HVGs for efficiency...")
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
        sc.pp.highly_variable_genes(adata_ref, n_top_genes=5000)
        adata_ref = adata_ref[:, adata_ref.var['highly_variable']].copy()

    # Resolve reference key for cell type grouping
    if ref_key not in adata_ref.obs.columns:
        print(f"    - Warning: Reference key '{ref_key}' not found. Using fallback.")
        for k in ['cell_type', 'cell_type_ontology_term_id', 'Cluster']:
            if k in adata_ref.obs.columns:
                ref_key = k
                break

    # Filter small clusters and validate cell type range
    if ref_key in adata_ref.obs.columns:
        counts = adata_ref.obs[ref_key].value_counts()
        valid_cts = counts[counts > 50].index
        adata_ref = adata_ref[adata_ref.obs[ref_key].isin(valid_cts)].copy()

        n_celltypes = adata_ref.obs[ref_key].nunique()
        if n_celltypes < min_celltypes or n_celltypes > max_celltypes:
            logger.warning(
                f"Dataset has {n_celltypes} cell types (expected {min_celltypes}-{max_celltypes}). "
                "Proceeding but results may be unreliable."
            )
        else:
            print(f"    - Cell types after filtering: {n_celltypes}")

        # Ensure normalization before rank_genes_groups
        if np.max(adata_ref.X) > 50:
             sc.pp.normalize_total(adata_ref, target_sum=1e4)
             sc.pp.log1p(adata_ref)

        # --- 7-2b. Rank marker genes (Wilcoxon) ---
        print(f"    - Identifying markers for simulation (groupby='{ref_key}')...")
        sc.tl.rank_genes_groups(adata_ref, groupby=ref_key, method='t-test')
        markers_df = pd.DataFrame(adata_ref.uns['rank_genes_groups']['names']).head(50)
    else:
        print("    - Error: No valid grouping key found. Cannot identify markers.")
        return []

    generated_files = []

    # --- 7-2c. Standard simulation (base noise/misseg) ---
    print(f"    - Creating Standard Simulation (Noise={noise_pct}%, MS={ms_pct}%)...")
    adata_sim1 = subset_of_single_cell(
        adata_ref.copy(),
        markers=markers_df,
        number_of_markers=n_markers,
        percentage_of_noise=noise_pct,
        ms_percentage=ms_pct
    )
    sim1_file = os.path.join(output_dir, f"{sample_tag}_simulated_standard.h5ad")
    adata_sim1.write_h5ad(sim1_file)
    generated_files.append(sim1_file)
    print(f"      > Saved: {sim1_file}")

    # --- 7-2d. High-noise simulation (2× noise/misseg) ---
    print(f"    - Creating High Noise Simulation (Noise={noise_pct*2}%)...")
    adata_sim2 = subset_of_single_cell(
        adata_ref.copy(),
        markers=markers_df,
        number_of_markers=n_markers,
        percentage_of_noise=noise_pct*2,
        ms_percentage=ms_pct*2
    )
    sim2_file = os.path.join(output_dir, f"{sample_tag}_simulated_noisy.h5ad")
    adata_sim2.write_h5ad(sim2_file)
    generated_files.append(sim2_file)

    return generated_files


# --- 7-3. Benchmarking ---

def run_step7_3_benchmarking(config, sim_file):
    """Grid-search preprocessing parameters and evaluate with ARI, NMI, FMI, VI."""
    print("\n[Step 7-3] Benchmarking Preprocessing (Ref: Notebook 6_3/6_4)...")
    if not sim_file or not os.path.exists(sim_file):
        print("    - Error: Simulation file required for benchmarking.")
        return None

    output_dir = config["output_dir"]
    sample_tag = config.get("sample_tag", "sample")
    sim_config = config.get("simulation", {})
    ref_key = sim_config.get("reference_key", "cell_type")

    adata = sc.read_h5ad(sim_file)
    if ref_key not in adata.obs.columns:
        print(f"    - Warning: Ground truth '{ref_key}' missing. Using first categorical obs.")
        cats = [c for c in adata.obs.columns if adata.obs[c].dtype.name == 'category']
        if cats:
            ref_key = cats[0]
        else:
             print("    - Error: No ground truth available for assessment.")
             return None

    gt_labels = adata.obs[ref_key].astype(str)

    # --- 7-3a. Build preprocessing grid ---
    grid_combos = build_preprocessing_grid(sim_config)
    n_combos = len(grid_combos)

    print(f"    - Running Grid Search ({n_combos} combinations)...")

    results = []

    # --- 7-3b. Grid search ---
    for idx, params in enumerate(tqdm(grid_combos, desc="    Benchmarking", unit="combo")):
        run_name = _params_to_run_name(params)

        try:
            ad_sub = _apply_preprocessing(adata, params)
            pred_labels = ad_sub.obs['leiden'].astype(str)

            nmi = compute_nmi(gt_labels, pred_labels)
            ari = adjusted_rand_score(gt_labels, pred_labels)
            fmi = compute_fmi(gt_labels, pred_labels)
            vi = compute_vi(gt_labels, pred_labels)
        except Exception as e:
            logger.warning(f"Combination {run_name} failed: {e}")
            nmi = ari = fmi = vi = np.nan

        results.append({
            'Run': run_name,
            'normalize': params['normalize'],
            'target_sum': params['target_sum'],
            'log1p': params['log1p'],
            'scale': params['scale'],
            'hvg': params['hvg'],
            'n_neighbors': params['n_neighbors'],
            'n_pcs': params['n_pcs'],
            'resolution': params['resolution'],
            'NMI': nmi,
            'ARI': ari,
            'FMI': fmi,
            'VI': vi,
        })

    bench_df = pd.DataFrame(results)
    bench_file = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_results.csv")
    bench_df.to_csv(bench_file, index=False)
    print(f"    - Saved Results: {bench_file}")

    # Best workflow by ARI
    best_row = bench_df.loc[bench_df['ARI'].idxmax()]
    best_params = {k: best_row[k] for k in [
        'normalize', 'target_sum', 'log1p', 'scale', 'hvg',
        'n_neighbors', 'n_pcs', 'resolution',
    ]}
    # Convert numpy types to native Python for JSON serialization
    best_params = {k: v.item() if hasattr(v, 'item') else v for k, v in best_params.items()}
    print(f"    - Best ARI = {best_row['ARI']:.4f} | {best_row['Run']}")

    # Save best params
    _save_best_params(bench_df, output_dir, sample_tag)

    # --- 7-3d. Benchmark plots ---
    _generate_benchmark_plots(bench_df, output_dir, sample_tag)

    run_perturbation_analysis(config, sim_file, best_params, gt_labels, adata, output_dir, sample_tag)

    return bench_df, best_params


# --- 7-3d. Benchmark visualization helpers ---

def _generate_benchmark_plots(bench_df, output_dir, sample_tag):
    """Generate all benchmark plots: bar charts, heatmap, importance, boxplot."""

    # --- ARI bar chart ---
    plt.figure(figsize=(max(10, len(bench_df) * 0.35), 6))
    order = bench_df.sort_values('ARI', ascending=False)['Run']
    sns.barplot(data=bench_df, x='Run', y='ARI', order=order)
    plt.xticks(rotation=90, fontsize=7)
    plt.title("ARI Score - Preprocessing Benchmarking")
    plt.tight_layout()
    plot_file_ari = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_ari.png")
    plt.savefig(plot_file_ari, dpi=150)
    plt.close()

    # --- NMI bar chart ---
    plt.figure(figsize=(max(10, len(bench_df) * 0.35), 6))
    sns.barplot(data=bench_df, x='Run', y='NMI', order=order)
    plt.xticks(rotation=90, fontsize=7)
    plt.title("NMI Score - Preprocessing Benchmarking")
    plt.tight_layout()
    plot_file_nmi = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_nmi.png")
    plt.savefig(plot_file_nmi, dpi=150)
    plt.close()

    # --- FMI bar chart ---
    plt.figure(figsize=(max(10, len(bench_df) * 0.35), 6))
    sns.barplot(data=bench_df, x='Run', y='FMI', order=order)
    plt.xticks(rotation=90, fontsize=7)
    plt.title("FMI Score - Preprocessing Benchmarking")
    plt.tight_layout()
    plot_file_fmi = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_fmi.png")
    plt.savefig(plot_file_fmi, dpi=150)
    plt.close()

    # --- ARI heatmap (n_neighbors x n_pcs, averaged over other params) ---
    try:
        pivot = bench_df.groupby(['n_neighbors', 'n_pcs'])['ARI'].mean().reset_index()
        heatmap_data = pivot.pivot(index='n_neighbors', columns='n_pcs', values='ARI')
        plt.figure(figsize=(8, 6))
        sns.heatmap(heatmap_data, annot=True, fmt=".3f", cmap="YlOrRd")
        plt.title("Mean ARI (n_neighbors x n_pcs, averaged over other params)")
        plt.tight_layout()
        plot_file_heat = os.path.join(output_dir, f"{sample_tag}_step7_ari_heatmap.png")
        plt.savefig(plot_file_heat, dpi=150)
        plt.close()
    except Exception as e:
        logger.warning(f"Could not generate ARI heatmap: {e}")

    # --- Parameter importance (ARI range per parameter) ---
    try:
        importance = {}
        param_cols = ['normalize', 'target_sum', 'scale', 'hvg', 'n_neighbors', 'n_pcs', 'resolution']
        for col in param_cols:
            unique_vals = bench_df[col].unique()
            if len(unique_vals) > 1:
                group_means = bench_df.groupby(col)['ARI'].mean()
                importance[col] = group_means.max() - group_means.min()
            else:
                importance[col] = 0.0

        imp_df = pd.DataFrame({
            'Parameter': list(importance.keys()),
            'ARI_range': list(importance.values()),
        }).sort_values('ARI_range', ascending=True)

        plt.figure(figsize=(8, 5))
        plt.barh(imp_df['Parameter'], imp_df['ARI_range'], color='steelblue')
        plt.xlabel("ARI Range (max group mean - min group mean)")
        plt.title("Parameter Importance for ARI")
        plt.tight_layout()
        plot_file_imp = os.path.join(output_dir, f"{sample_tag}_step7_param_importance.png")
        plt.savefig(plot_file_imp, dpi=150)
        plt.close()
    except Exception as e:
        logger.warning(f"Could not generate parameter importance plot: {e}")

    # --- Metric distributions boxplot ---
    try:
        metric_cols = ['ARI', 'NMI', 'FMI', 'VI']
        present_cols = [c for c in metric_cols if c in bench_df.columns]
        melted = bench_df[present_cols].melt(var_name='Metric', value_name='Score')

        plt.figure(figsize=(8, 5))
        sns.boxplot(data=melted, x='Metric', y='Score')
        plt.title("Metric Distributions Across Preprocessing Grid")
        plt.tight_layout()
        plot_file_box = os.path.join(output_dir, f"{sample_tag}_step7_metric_boxplot.png")
        plt.savefig(plot_file_box, dpi=150)
        plt.close()
    except Exception as e:
        logger.warning(f"Could not generate metric box plot: {e}")

    print(f"    - Saved benchmark plots to {output_dir}")


# --- 7-4. Perturbation Analysis ---

def run_perturbation_analysis(config, sim_file, best_params, gt_labels, adata, output_dir, sample_tag):
    """Vary each parameter one-at-a-time from the best workflow to assess sensitivity."""
    print("\n    [Perturbation] Analysing sensitivity of best workflow...")

    sim_config = config.get("simulation", {})
    grid_spec = sim_config.get("preprocessing_grid", {})

    param_ranges = {
        'normalize':   grid_spec.get('normalize',   [True, False]),
        'target_sum':  grid_spec.get('target_sum',  [100, 1000, 10000]),
        'log1p':       grid_spec.get('log1p',       [True]),
        'scale':       grid_spec.get('scale',       [False, True]),
        'hvg':         grid_spec.get('hvg',         [False]),
        'n_neighbors': grid_spec.get('n_neighbors', [5, 15, 30]),
        'n_pcs':       grid_spec.get('n_pcs',       [0, 10, 30]),
        'resolution':  grid_spec.get('resolution',  [0.5, 1.0, 1.5]),
    }

    perturb_results = []

    for param_name, values in param_ranges.items():
        for val in values:
            test_params = dict(best_params)
            test_params[param_name] = val

            run_label = f"{param_name}={val}"

            try:
                ad_sub = _apply_preprocessing(adata, test_params)
                pred_labels = ad_sub.obs['leiden'].astype(str)

                ari = adjusted_rand_score(gt_labels, pred_labels)
                nmi = compute_nmi(gt_labels, pred_labels)
                fmi = compute_fmi(gt_labels, pred_labels)
            except Exception as e:
                logger.warning(f"Perturbation {run_label} failed: {e}")
                ari = nmi = fmi = np.nan

            is_best = (val == best_params.get(param_name))
            perturb_results.append({
                'parameter': param_name,
                'value': val,
                'is_best': is_best,
                'ARI': ari,
                'NMI': nmi,
                'FMI': fmi,
            })

    perturb_df = pd.DataFrame(perturb_results)
    perturb_file = os.path.join(output_dir, f"{sample_tag}_step7_perturbation.csv")
    perturb_df.to_csv(perturb_file, index=False)
    print(f"    [Perturbation] Saved results: {perturb_file}")

    # --- Perturbation grouped bar chart ---
    try:
        perturb_df['label'] = perturb_df['parameter'] + '=' + perturb_df['value'].astype(str)
        fig, ax = plt.subplots(figsize=(max(10, len(perturb_df) * 0.4), 6))
        x = np.arange(len(perturb_df))
        width = 0.25
        ax.bar(x - width, perturb_df['ARI'], width, label='ARI', color='steelblue')
        ax.bar(x,         perturb_df['NMI'], width, label='NMI', color='coral')
        ax.bar(x + width, perturb_df['FMI'], width, label='FMI', color='seagreen')

        # Highlight best-parameter bars
        for i, row in perturb_df.iterrows():
            if row['is_best']:
                ax.axvspan(i - 0.4, i + 0.4, alpha=0.12, color='gold')

        ax.set_xticks(x)
        ax.set_xticklabels(perturb_df['label'], rotation=90, fontsize=7)
        ax.set_ylabel("Score")
        ax.set_title("Perturbation Analysis (yellow = best workflow value)")
        ax.legend()
        plt.tight_layout()
        plot_file = os.path.join(output_dir, f"{sample_tag}_step7_perturbation.png")
        plt.savefig(plot_file, dpi=150)
        plt.close()
        print(f"    [Perturbation] Saved plot: {plot_file}")
    except Exception as e:
        logger.warning(f"Could not generate perturbation plot: {e}")


# --- 7-5. Best params persistence & real data validation ---

def _save_best_params(bench_df, output_dir, sample_tag):
    """Extract best-ARI row from benchmark results and save as JSON."""
    best_row = bench_df.loc[bench_df['ARI'].idxmax()]
    param_keys = ['normalize', 'target_sum', 'log1p', 'scale', 'hvg',
                  'n_neighbors', 'n_pcs', 'resolution']
    best = {}
    for k in param_keys:
        v = best_row[k]
        best[k] = v.item() if hasattr(v, 'item') else v
    best['best_ari'] = float(best_row['ARI'])

    out_path = os.path.join(output_dir, f"{sample_tag}_step7_best_params.json")
    with open(out_path, 'w') as f:
        json.dump(best, f, indent=2)
    print(f"    - Saved best params: {out_path}")
    return best


def _load_best_params(output_dir, sample_tag):
    """Load best params JSON; fall back to PAPER_BEST_PARAMS."""
    path = os.path.join(output_dir, f"{sample_tag}_step7_best_params.json")
    if os.path.exists(path):
        with open(path, 'r') as f:
            params = json.load(f)
        # Remove metadata keys
        params.pop('best_ari', None)
        logger.info(f"Loaded best params from {path}")
        return params
    logger.info("No best params JSON found — using PAPER_BEST_PARAMS")
    return dict(PAPER_BEST_PARAMS)


def _run_single_param_variants(adata, default_labels):
    """Vary one param at a time from DEFAULT_PARAMS (notebook 1_6 nv==1 logic).

    Returns DataFrame with columns: variant, param, value, ARI, NMI, FMI, VI.
    """
    variant_grid = {
        'n_neighbors': [6, 12, 20],
        'n_pcs':       [15, 25],
        'target_sum':  [10, 1000, None],
        'scale':       [False],
        'hvg':         [True],
        'normalize':   [False],
        'log1p':       [False],
    }

    results = []
    for param_name, values in variant_grid.items():
        for val in values:
            # Skip if same as default
            if val == DEFAULT_PARAMS.get(param_name):
                continue
            test_params = dict(DEFAULT_PARAMS)
            test_params[param_name] = val
            # normalize=False ⇒ target_sum irrelevant
            if not test_params['normalize']:
                test_params['target_sum'] = 100  # placeholder

            label = f"{param_name}={val}"
            try:
                ad = _apply_preprocessing(adata, test_params)
                pred = ad.obs['leiden'].astype(str)
                ari = adjusted_rand_score(default_labels, pred)
                nmi = compute_nmi(default_labels, pred)
                fmi = compute_fmi(default_labels, pred)
                vi = compute_vi(default_labels, pred)
            except Exception as e:
                logger.warning(f"Variant {label} failed: {e}")
                ari = nmi = fmi = vi = np.nan

            results.append({
                'variant': label,
                'param': param_name,
                'value': val if val is not None else 'None',
                'ARI': ari, 'NMI': nmi, 'FMI': fmi, 'VI': vi,
            })

    return pd.DataFrame(results)


def _generate_real_validation_plots(ad_default, ad_optimized, metrics,
                                     variant_df, best_params,
                                     output_dir, sample_tag):
    """Generate 5 real-validation plots."""

    # 1. UMAP comparison: default vs optimized
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(ad_default, color='leiden', ax=axes[0], show=False, title='Default params')
    sc.pl.umap(ad_optimized, color='leiden', ax=axes[1], show=False, title='Optimized params')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{sample_tag}_step7_real_umap_comparison.png"), dpi=150)
    plt.close(fig)

    # 2. Grouped bar: ARI/NMI/FMI comparison
    metric_names = ['ARI', 'NMI', 'FMI']
    metric_vals = [metrics.get(m, 0) for m in metric_names]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(metric_names, metric_vals, color=['steelblue', 'coral', 'seagreen'])
    ax.set_ylabel('Score')
    ax.set_title('Default vs Optimized Clustering Agreement')
    ax.set_ylim(0, 1.05)
    for i, v in enumerate(metric_vals):
        ax.text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{sample_tag}_step7_real_default_vs_optimized.png"), dpi=150)
    plt.close(fig)

    # 3. Single-param variant ARI barplot (Fig 4e equivalent)
    if variant_df is not None and len(variant_df) > 0:
        fig, ax = plt.subplots(figsize=(8, max(4, len(variant_df) * 0.4)))
        variant_sorted = variant_df.sort_values('ARI', ascending=True)
        ax.barh(variant_sorted['variant'], variant_sorted['ARI'], color='steelblue')
        ax.set_xlabel('ARI vs Default')
        ax.set_title('Single-Parameter Variant Impact (Fig 4e)')
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{sample_tag}_step7_real_variant_ari.png"), dpi=150)
        plt.close(fig)

        # 4. Variant metrics heatmap
        try:
            heat_data = variant_sorted.set_index('variant')[['ARI', 'NMI', 'FMI', 'VI']]
            fig, ax = plt.subplots(figsize=(6, max(4, len(variant_df) * 0.35)))
            sns.heatmap(heat_data, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax)
            ax.set_title('Variant Metrics Heatmap')
            plt.tight_layout()
            fig.savefig(os.path.join(output_dir, f"{sample_tag}_step7_real_variant_metrics_heatmap.png"), dpi=150)
            plt.close(fig)
        except Exception as e:
            logger.warning(f"Could not generate variant heatmap: {e}")

    # 5. Param diff table
    try:
        rows = []
        for k in DEFAULT_PARAMS:
            rows.append({
                'Parameter': k,
                'Default': str(DEFAULT_PARAMS[k]),
                'Best': str(best_params.get(k, DEFAULT_PARAMS[k])),
            })
        diff_df = pd.DataFrame(rows)
        diff_df['Changed'] = diff_df['Default'] != diff_df['Best']

        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis('off')
        table = ax.table(cellText=diff_df.values, colLabels=diff_df.columns,
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.4)
        # Highlight changed rows
        for i, changed in enumerate(diff_df['Changed']):
            if changed:
                for j in range(len(diff_df.columns)):
                    table[i + 1, j].set_facecolor('#FFFFCC')
        ax.set_title('Default vs Best Parameters', fontsize=12, pad=20)
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{sample_tag}_step7_real_param_diff.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        logger.warning(f"Could not generate param diff table: {e}")

    print(f"    - Saved real validation plots to {output_dir}")


def run_step7_5_real_validation(config, real_adata_path, best_params):
    """Apply best preprocessing params to real Xenium data (Paper Fig 4e).

    1. Load real adata (Step 0 raw counts)
    2. Apply DEFAULT_PARAMS → ad_default (with UMAP)
    3. Apply best_params → ad_optimized (with UMAP)
    4. Compute ARI/NMI/FMI/VI between default vs optimized clustering
    5. Run single-param variants (if config enabled)
    6. Generate plots
    7. Save optimized adata with layers['raw'] preserved
    8. Save variant metrics CSV
    9. Return optimized adata path
    """
    print("\n[Step 7-5] Real Data Validation (Ref: Notebook 1_6)...")

    output_dir = config["output_dir"]
    sample_tag = config.get("sample_tag", "sample")
    sim_config = config.get("simulation", {})
    rv_config = sim_config.get("real_validation", {})

    if not rv_config.get("enabled", True):
        print("    - Skipping 7-5 (real_validation.enabled=false)")
        return None

    if not real_adata_path or not os.path.exists(real_adata_path):
        print(f"    - Error: Real adata not found at {real_adata_path}")
        return None

    # Load real data
    print(f"    - Loading real data: {real_adata_path}")
    adata = sc.read_h5ad(real_adata_path)
    print(f"    - Real data: {adata.n_obs} cells, {adata.n_vars} genes")

    # Preserve raw counts
    if 'raw' not in adata.layers:
        if hasattr(adata.X, 'toarray'):
            adata.layers['raw'] = adata.X.toarray().copy()
        else:
            adata.layers['raw'] = adata.X.copy()

    # Apply DEFAULT preprocessing
    print("    - Applying DEFAULT preprocessing...")
    ad_default = _apply_preprocessing(adata, DEFAULT_PARAMS, compute_umap=True)
    default_labels = ad_default.obs['leiden'].astype(str)

    # Apply BEST preprocessing
    print(f"    - Applying BEST preprocessing: {best_params}")
    ad_optimized = _apply_preprocessing(adata, best_params, compute_umap=True)
    optimized_labels = ad_optimized.obs['leiden'].astype(str)

    # Compute agreement metrics (default vs optimized)
    metrics = {
        'ARI': adjusted_rand_score(default_labels, optimized_labels),
        'NMI': compute_nmi(default_labels, optimized_labels),
        'FMI': compute_fmi(default_labels, optimized_labels),
        'VI':  compute_vi(default_labels, optimized_labels),
    }
    print(f"    - Default vs Optimized: ARI={metrics['ARI']:.4f}, NMI={metrics['NMI']:.4f}")

    # Single-param variants
    variant_df = None
    if rv_config.get("run_single_param_variants", True):
        print("    - Running single-param variants...")
        variant_df = _run_single_param_variants(adata, default_labels)
        variant_csv = os.path.join(output_dir, f"{sample_tag}_step7_real_variants.csv")
        variant_df.to_csv(variant_csv, index=False)
        print(f"    - Saved variant metrics: {variant_csv}")

    # Generate plots
    _generate_real_validation_plots(
        ad_default, ad_optimized, metrics, variant_df,
        best_params, output_dir, sample_tag,
    )

    # Save optimized adata (preserve raw layer)
    ad_optimized.layers['raw'] = adata.layers['raw']
    opt_path = os.path.join(output_dir, f"{sample_tag}_step7_optimized.h5ad")
    ad_optimized.write_h5ad(opt_path)
    print(f"    - Saved optimized adata: {opt_path}")

    # Save metrics summary
    metrics_file = os.path.join(output_dir, f"{sample_tag}_step7_real_metrics.json")
    with open(metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)

    return opt_path


def _cleanup_simulation_files(config, sim_files, ref_file):
    """Delete simulated h5ad files after benchmarking to save disk space."""
    sim_config = config.get("simulation", {})
    rv_config = sim_config.get("real_validation", {})

    if rv_config.get("cleanup_intermediates", True):
        for f in sim_files:
            if os.path.exists(f):
                os.remove(f)
                print(f"    - Cleaned up: {f}")

    if rv_config.get("cleanup_reference", False):
        if ref_file and os.path.exists(ref_file):
            os.remove(ref_file)
            print(f"    - Cleaned up reference: {ref_file}")


# --- Main Orchestrator ---

def run_step7(config, real_adata_path=None):
    print("\n" + "="*60)
    print("[Step 7] Simulation & Benchmarking (Notebooks 6_1 - 6_4)")
    print("="*60)

    output_dir = config["output_dir"]
    sample_tag = config.get("sample_tag", "sample")
    sim_config = config.get("simulation", {})

    best_params = None
    sim_files = []
    ref_file = None

    if sim_config.get("run_simulation", False):
        ref_file = run_step7_1_acquisition(config)
        if not ref_file:
            print("    - Step 7-1 Failed. Stopping simulation.")
        else:
            sim_files = run_step7_2_simulation(config, ref_file)
            if sim_files:
                # Benchmark the standard simulation (first file)
                result = run_step7_3_benchmarking(config, sim_files[0])
                if result is not None:
                    _bench_df, best_params = result
            else:
                print("    - Step 7-2 Failed. Stopping simulation.")
    else:
        print("    - Skipping simulation (run_simulation=False)")

    # Fallback to paper best params if simulation didn't run or failed
    if best_params is None:
        best_params = _load_best_params(output_dir, sample_tag)

    # 7-5: Real data validation
    optimized_adata_path = None
    if real_adata_path:
        optimized_adata_path = run_step7_5_real_validation(
            config, real_adata_path, best_params,
        )

    # Cleanup simulation intermediates
    if sim_files:
        _cleanup_simulation_files(config, sim_files, ref_file)

    print("\n=== Step 7 Complete ===")
    return optimized_adata_path

if __name__ == "__main__":
    import yaml
