# pipeline/xenium_step7_simulation.py
# ------------------------------------------------------------
# Xenium Pipeline Step 7: Simulation & Benchmarking
# Simulates Xenium data from scRNAseq reference and benchmarks preprocessing.
# Ref: Notebooks 6_1 to 6_4
# ------------------------------------------------------------

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import random
import math
from tqdm import tqdm
from scipy.spatial import ConvexHull
from sklearn.metrics import mutual_info_score, silhouette_score, fowlkes_mallows_score, normalized_mutual_info_score, adjusted_rand_score
import alphashape
from shapely.geometry import Point, Polygon
import warnings

# Try importing cellxgene_census (for 6_1)
try:
    import cellxgene_census
except ImportError:
    cellxgene_census = None

# Configure local logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
warnings.filterwarnings('ignore')

# ------------------------------------------------------------
# 1. XB Library Functions (Inlined for Stability)
# Ref: xb/simulating.py, xb/calculating.py
# ------------------------------------------------------------

def missegmentation_simulation(adata_sc_sub, missegmentation_percentage=0.1):
    """ Simulate missegmentation mixing profiles. """
    exp = adata_sc_sub.to_df() # Dense dataframe
    if missegmentation_percentage > 0:
        cells_affected = int(exp.shape[0] * (missegmentation_percentage / 100))
        # Use simple sampling
        indices = list(exp.index)
        # Vectorized mixing might be hard, loop is explicit in source
        # To speed up, we sample all at once
        sources = random.sample(indices, cells_affected)
        targets = random.sample(indices, cells_affected)
        
        # Mixing factor
        # Source uses random 0-1 * 0.1? Source says: random.sample(range(0,10),1)[0]*0.1
        # Let's simplify
        
        # Iterating might be slow for large N.
        # But simulation usually is on subset.
        for s, t in zip(sources, targets):
             factor = random.random() # 0-1
             # exp.loc[t] += exp.loc[s] * factor (simplified logic from xb)
             # xb logic: exp.loc[target,:]=exp.loc[target,:]+(exp.loc[source,:]*missegmentation_importance)
             # missegmentation_importance is 0.0, 0.1, ... 0.9?
             importance = random.randint(0, 9) * 0.1
             if importance > 0:
                exp.loc[t] += exp.loc[s] * importance
                
    adata_sc_sub.X = np.array(exp.values) # Keep simple
    return adata_sc_sub

def noise_adder(adata_sc, percentage_of_noise=0.1):
    """ Add random noise to counts. """
    # Total counts * percentage
    total_counts = np.sum(adata_sc.X)
    noise_events = int(total_counts * (percentage_of_noise / 100))
    
    if noise_events > 0:
        # Randomly choose coordinates to perturb
        # If matrix is sparse, this is slow. Expect dense for simulation.
        rows = np.random.randint(0, adata_sc.X.shape[0], noise_events)
        cols = np.random.randint(0, adata_sc.X.shape[1], noise_events)
        changes = np.random.choice([-1, 1], noise_events)
        
        # Apply changes (naive loop or at once?)
        # Numpy advanced indexing allow duplicates?
        # Better iteratively or just accept advanced indexing overwrites (statistically consistent)
        # xb loops. We verify matrix type.
        if hasattr(adata_sc.X, "tocoo"):
             adata_sc.X = adata_sc.X.toarray()
             
        # Add noise
        # Note: X[rows, cols] += changes works but if same index picked twice, only one adds.
        # Given sparsity, collision unlikely.
        # To be safe and fast:
        np.add.at(adata_sc.X, (rows, cols), changes)
        
        # Clip negative
        adata_sc.X[adata_sc.X < 0] = 0
        
    return adata_sc

def subset_of_single_cell(adata_sc_sub, markers, random_markers_percentage=0,
                          reads_x_cell=None, number_of_markers=200,
                          percentage_of_noise=0.1, ms_percentage=0.1):
    """ Transform sc data to simulated spatial data. """
    
    # Select Genes (Markers + Random)
    # 1. Flatten markers
    mk = []
    # Logic from xb: iterate columns, pick top N
    # Simplified: Get all unique markers
    # markers is DataFrame where columns are clusters, index is rank.
    unique_markers = np.unique(markers.values.flatten())
    mk = list(unique_markers)
    
    # 2. Add Random
    all_genes = list(adata_sc_sub.var.index)
    if random_markers_percentage > 0:
         n_random = int(len(mk) * (random_markers_percentage / 100))
         random_genes = random.sample(all_genes, n_random)
         mk.extend(random_genes)
         
    # 3. Cap number
    if len(mk) > number_of_markers:
        mk = random.sample(mk, number_of_markers)
        
    # Filter Genes
    adata_sc = adata_sc_sub[:, adata_sc_sub.var.index.isin(mk)].copy()
    
    # Verify we have data
    if adata_sc.n_vars == 0:
        print("    [Sim] Warning: No marker genes found in adata. Check naming convention.")
        return adata_sc_sub # Fallback
        
    # Simulation: Missegmentation + Noise
    adata_sc = missegmentation_simulation(adata_sc, missegmentation_percentage=ms_percentage)
    adata_sc = noise_adder(adata_sc, percentage_of_noise=percentage_of_noise)
    
    # Downsample
    if reads_x_cell is not None:
         # Check if we have enough counts
         if np.median(np.sum(adata_sc.X, axis=1)) > reads_x_cell:
             sc.pp.downsample_counts(adata_sc, counts_per_cell=reads_x_cell)
             
    return adata_sc

def entropy(clustering):
    """ Compute entropy """
    _, counts = np.unique(clustering, return_counts=True)
    proportions = counts / len(clustering)
    return -np.sum(proportions * np.log(proportions))

def compute_vi(ground_truth, predicted):
    """ Variation of Information """
    mi = mutual_info_score(ground_truth, predicted)
    h_gt = entropy(ground_truth)
    h_pred = entropy(predicted)
    return h_gt + h_pred - 2 * mi

def compute_nmi(ground_truth, predicted):
    """ Normalized Mutual Information """
    return normalized_mutual_info_score(ground_truth, predicted)

def compute_fmi(ground_truth, predicted):
    """ Fowlkes-Mallows Index """
    return fowlkes_mallows_score(ground_truth, predicted)


# ------------------------------------------------------------
# 2. Sub-step Functions (Modularized by Notebook)
# ------------------------------------------------------------

def run_step7_1_acquisition(config):
    """
    [Step 7-1] Data Acquisition (Ref: Notebook 6_1)
    Goal: Obtain scRNA-seq reference data (Census or Manual).
    """
    print("\n[Step 7-1] Acquiring Reference Data (Ref: Notebook 6_1)...")
    output_dir = config["output_dir"]
    sim_config = config.get("simulation", {})
    tissue = sim_config.get("census_tissue", "brain")
    organism = sim_config.get("census_organism", "mus_musculus")
    
    sc_file = os.path.join(output_dir, f"census_{tissue}_reference.h5ad")
    
    if os.path.exists(sc_file):
        print(f"    - Found existing reference: {sc_file}")
        return sc_file
        
    # Download logic
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
        # Save raw
        adata_ref.write_h5ad(sc_file)
        print(f"    - Downloaded and saved {adata_ref.shape} cells.")
        return sc_file
    except Exception as e:
        print(f"    - Census download failed: {e}")
        return None

def run_step7_2_simulation(config, ref_file):
    """
    [Step 7-2] Simulation Generation (Ref: Notebook 6_2)
    Goal: Generate synthetic Xenium datasets using 'xb' logic (Missegmentation, Noise).
    """
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
    
    adata_ref = sc.read_h5ad(ref_file)
    print(f"    - Loaded Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    # QC / Filtering of Reference (Internal Step 1 logic)
    if adata_ref.n_vars > 5000:
        print("    - Filtering to top 5000 HVGs for efficiency...")
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
        sc.pp.highly_variable_genes(adata_ref, n_top_genes=5000)
        adata_ref = adata_ref[:, adata_ref.var['highly_variable']].copy()
        
    # Find markers (Requirement for Simulation)
    if ref_key not in adata_ref.obs.columns:
        print(f"    - Warning: Reference key '{ref_key}' not found. Using fallback.")
        for k in ['cell_type', 'cell_type_ontology_term_id', 'Cluster']:
            if k in adata_ref.obs.columns:
                ref_key = k
                break
    
    # Filter small clusters
    if ref_key in adata_ref.obs.columns:
        counts = adata_ref.obs[ref_key].value_counts()
        valid_cts = counts[counts > 50].index
        adata_ref = adata_ref[adata_ref.obs[ref_key].isin(valid_cts)].copy()
        
        # Ensure Normalization before ranking
        if np.max(adata_ref.X) > 50:
             sc.pp.normalize_total(adata_ref, target_sum=1e4)
             sc.pp.log1p(adata_ref)
             
        print(f"    - Identifying markers for simulation (groupby='{ref_key}')...")
        sc.tl.rank_genes_groups(adata_ref, groupby=ref_key, method='t-test')
        markers_df = pd.DataFrame(adata_ref.uns['rank_genes_groups']['names']).head(50)
    else:
        print("    - Error: No valid grouping key found. Cannot identify markers.")
        return []

    # Run Simulation
    generated_files = []
    
    # 1. Standard Simulation
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
    
    # 2. High Noise Simulation (Optional comparison)
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

def run_step7_3_benchmarking(config, sim_file):
    """
    [Step 7-3] Benchmarking & Assessment (Ref: Notebook 6_3/6_4)
    Goal: Test preprocessing parameter combinations and assess using ARI, NMI.
    """
    print("\n[Step 7-3] Benchmarking Preprocessing (Ref: Notebook 6_3/6_4)...")
    if not sim_file or not os.path.exists(sim_file):
        print("    - Error: Simulation file required for benchmarking.")
        return

    output_dir = config["output_dir"]
    sample_tag = config.get("sample_tag", "sample")
    sim_config = config.get("simulation", {})
    ref_key = sim_config.get("reference_key", "cell_type")

    # Load Simulation
    adata = sc.read_h5ad(sim_file)
    if ref_key not in adata.obs.columns:
        print(f"    - Warning: Ground truth '{ref_key}' missing. Using first categorical obs.")
        # Fallback
        cats = [c for c in adata.obs.columns if adata.obs[c].dtype.name == 'category']
        if cats: ref_key = cats[0]
        else:
             print("    - Error: No ground truth available for assessment.")
             return

    gt_labels = adata.obs[ref_key].astype(str)
    
    # Define Grid (Notebook 6_3)
    grid_neigh = [5, 15, 30]
    grid_pca = [10, 30]
    grid_res = [0.5, 1.0, 1.5]
    
    print(f"    - Running Grid Search ({len(grid_neigh)*len(grid_pca)*len(grid_res)} combinations)...")
    
    results = []
    idx = 0
    
    for n_neigh in grid_neigh:
        for n_pca in grid_pca:
            for res in grid_res:
                idx += 1
                run_name = f"N{n_neigh}_P{n_pca}_R{res}"
                
                # Preprocess & Cluster (Notebook 6_3 Scanpy Logic)
                ad_sub = adata.copy()
                sc.pp.normalize_total(ad_sub, target_sum=1e4)
                sc.pp.log1p(ad_sub)
                sc.pp.scale(ad_sub, max_value=10)
                sc.pp.pca(ad_sub, n_comps=n_pca)
                sc.pp.neighbors(ad_sub, n_neighbors=n_neigh, n_pcs=n_pca)
                sc.tl.leiden(ad_sub, resolution=res, key_added='leiden')
                
                # Assess (Notebook 6_4 Logic)
                pred_labels = ad_sub.obs['leiden'].astype(str)
                
                nmi = compute_nmi(gt_labels, pred_labels)
                ari = adjusted_rand_score(gt_labels, pred_labels)
                fmi = compute_fmi(gt_labels, pred_labels)
                vi = compute_vi(gt_labels, pred_labels)
                
                results.append({
                    'Run': run_name,
                    'Neighbors': n_neigh,
                    'PCA': n_pca,
                    'Resolution': res,
                    'NMI': nmi,
                    'ARI': ari,
                    'FMI': fmi,
                    'VI': vi
                })
                
    # Save Results
    bench_df = pd.DataFrame(results)
    bench_file = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_results.csv")
    bench_df.to_csv(bench_file, index=False)
    print(f"    - Saved Results: {bench_file}")
    
    # Visualization (Notebook 6_4)
    # 1. ARI Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=bench_df, x='Run', y='ARI')
    plt.xticks(rotation=45, ha='right')
    plt.title("ARI Score - Preprocessing Benchmarking")
    plt.tight_layout()
    plot_file_ari = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_ari.png")
    plt.savefig(plot_file_ari)
    plt.close()
    
    # 2. NMI Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=bench_df, x='Run', y='NMI')
    plt.xticks(rotation=45, ha='right')
    plt.title("NMI Score - Preprocessing Benchmarking")
    plt.tight_layout()
    plot_file_nmi = os.path.join(output_dir, f"{sample_tag}_step7_benchmark_nmi.png")
    plt.savefig(plot_file_nmi)
    plt.close()
    
    print(f"    - Saved Plots: {plot_file_ari}, {plot_file_nmi}")


# ------------------------------------------------------------
# 3. Main Step 7 Orchestrator
# ------------------------------------------------------------

def run_step7(config):
    print("\n" + "="*60)
    print("[Step 7] Simulation & Benchmarking (Notebooks 6_1 - 6_4)")
    print("="*60)

    sim_config = config.get("simulation", {})
    if not sim_config.get("run_simulation", False):
        print("    - Skipping Step 7 (run_simulation=False)")
        return

    # [Step 7-1] Acquisition (Notebook 6_1)
    ref_file = run_step7_1_acquisition(config)
    if not ref_file:
         print("    - Step 7-1 Failed. Stopping.")
         return

    # [Step 7-2] Simulation (Notebook 6_2)
    sim_files = run_step7_2_simulation(config, ref_file)
    if not sim_files:
         print("    - Step 7-2 Failed. Stopping.")
         return
    
    # [Step 7-3] Benchmarking (Notebook 6_3/6_4)
    # We benchmark distinct simulations. usually focused on the standard one.
    standard_sim = sim_files[0]
    run_step7_3_benchmarking(config, standard_sim)
    
    print("\n=== Step 7 Complete ===")

if __name__ == "__main__":
    import yaml
    # Test run
    pass
    pass
