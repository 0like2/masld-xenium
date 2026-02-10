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


def _apply_preprocessing(adata, params):
    """Apply one preprocessing workflow + Leiden clustering. Returns modified copy."""
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

    if n_pcs > 0:
        n_comps = min(n_pcs, ad.n_vars - 1, ad.n_obs - 1)
        if n_comps < 2:
            n_pcs = 0
        else:
            sc.pp.pca(ad, n_comps=n_comps)

    if n_pcs > 0:
        sc.pp.neighbors(ad, n_neighbors=params['n_neighbors'], n_pcs=n_pcs)
    else:
        sc.pp.neighbors(ad, n_neighbors=params['n_neighbors'], use_rep='X')

    sc.tl.leiden(ad, resolution=params['resolution'], key_added='leiden')
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
    print(f"    - Best ARI = {best_row['ARI']:.4f} | {best_row['Run']}")

    # --- 7-3d. Benchmark plots ---
    _generate_benchmark_plots(bench_df, output_dir, sample_tag)

    run_perturbation_analysis(config, sim_file, best_params, gt_labels, adata, output_dir, sample_tag)

    return bench_df


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


# --- Main Orchestrator ---

def run_step7(config):
    print("\n" + "="*60)
    print("[Step 7] Simulation & Benchmarking (Notebooks 6_1 - 6_4)")
    print("="*60)

    sim_config = config.get("simulation", {})
    if not sim_config.get("run_simulation", False):
        print("    - Skipping Step 7 (run_simulation=False)")
        return

    ref_file = run_step7_1_acquisition(config)
    if not ref_file:
         print("    - Step 7-1 Failed. Stopping.")
         return

    sim_files = run_step7_2_simulation(config, ref_file)
    if not sim_files:
         print("    - Step 7-2 Failed. Stopping.")
         return

    # Benchmark the standard simulation (first file)
    standard_sim = sim_files[0]
    run_step7_3_benchmarking(config, standard_sim)

    print("\n=== Step 7 Complete ===")

if __name__ == "__main__":
    import yaml
