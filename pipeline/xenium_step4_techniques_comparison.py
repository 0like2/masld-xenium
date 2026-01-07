import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import xb.comparing as comparing
import xb.calculating as calculating
import xb.plotting as plotting

def run_step4(config, adata_path, transcripts_path, output_dir):
    """
    Step 4: Techniques Comparison & Validation
    Calculates efficiency, specificity, and validates resegmentation (Positivity, Diffusion).
    Ref: Notebooks 3.3, 3.4, 3.5, 3.6
    """
    print("------------------------------------------------")
    print("Starting Step 4: Techniques Comparison & Validation")
    print("------------------------------------------------")

    # 1. Load Data (Resegmented Data from Step 3)
    # adata_path passed from pipeline_main should point to Step 3 output
    print(f"Loading Resegmented AnnData from {adata_path}...")
    try:
        adata = sc.read_h5ad(adata_path)
    except Exception as e:
        print(f"Error loading AnnData: {e}")
        return

    # Load Resegmented Transcripts (Needed for Diffusion Analysis)
    # We expect 'transcripts_path' to point to the resegmented csv from Step 3
    print(f"Loading Resegmented Transcripts from {transcripts_path}...")
    try:
        df_assigned = pd.read_csv(transcripts_path)
    except Exception as e:
        print(f"Error loading transcripts (needed for diffusion): {e}")
        # We might proceed without diffusion if this fails, but better to return or skip diffusion
        df_assigned = None
    
    # 2. Config
    comp_config = config.get('comparison', {})
    # We can use 'comparison' config for all metrics, or check specific flags
    # Default to running all for validation
    run_efficiency = comp_config.get('run_efficiency', True)
    run_specificity = comp_config.get('run_specificity', True)
    run_positivity = comp_config.get('run_positivity', True) # New flag option
    run_diffusion = comp_config.get('run_diffusion', True)   # New flag option
    
    figs_dir = os.path.join(output_dir, 'figures', '4_techniques_comparison')
    os.makedirs(figs_dir, exist_ok=True)
    
    sample_tag = config.get('sample_tag', 'sample')

    # 3. Efficiency Analysis (3.3)
    if run_efficiency:
        print("Running Efficiency Analysis...")
        try:
           analyze_efficiency(adata, figs_dir)
        except Exception as e:
            print(f"Error in Efficiency Analysis: {e}")

    # 4. Specificity Analysis (3.4)
    if run_specificity:
        print("Running Specificity Analysis (NMP Proxy)...")
        try:
            analyze_specificity(adata, comp_config, figs_dir)
        except Exception as e:
            print(f"Error in Specificity Analysis: {e}")
            
    # 5. Positivity Analysis (3.5)
    if run_positivity:
        print("Running Positivity Analysis (Validation)...")
        try:
            analyze_positivity(adata, figs_dir, sample_tag)
        except Exception as e:
            print(f"Error in Positivity Analysis: {e}")

    # 6. Diffusion Analysis (3.6)
    if run_diffusion and df_assigned is not None:
        print("Running Diffusion Analysis (Validation)...")
        try:
            analyze_diffusion(df_assigned, adata, figs_dir)
        except Exception as e:
            print(f"Error in Diffusion Analysis: {e}")

    print("Step 4: Techniques Comparison & Validation Completed.")


def analyze_efficiency(adata, output_dir):
    """
    Calculates and plots efficiency metrics (Transcripts per cell, Genes per cell).
    Ref: Notebook 3.3
    """
    if 'total_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
    plt.figure(figsize=(6, 4))
    sns.histplot(adata.obs['total_counts'], kde=True, bins=50)
    plt.title('Efficiency: Transcripts per Cell')
    plt.xlabel('Transcripts per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_transcripts_per_cell.png'))
    plt.close()

    plt.figure(figsize=(6, 4))
    sns.histplot(adata.obs['n_genes_by_counts'], kde=True, bins=50)
    plt.title('Efficiency: Genes Detected per Cell')
    plt.xlabel('Genes per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_genes_per_cell.png'))
    plt.close()
    
    stats = {
        'median_transcripts_per_cell': np.median(adata.obs['total_counts']),
        'mean_transcripts_per_cell': np.mean(adata.obs['total_counts']),
        'median_genes_per_cell': np.median(adata.obs['n_genes_by_counts']),
        'mean_genes_per_cell': np.mean(adata.obs['n_genes_by_counts'])
    }
    pd.DataFrame([stats]).to_csv(os.path.join(output_dir, 'efficiency_metrics.csv'), index=False)
    print(f"  > Efficiency metrics saved to {output_dir}")


def analyze_specificity(adata, config, output_dir):
    """
    Calculates Negative Marker Purity (NMP) proxy via Gene Correlation.
    Ref: Notebook 3.4
    """
    print("  > Calculating Co-expression (Gene-Gene correlation)...")
    
    # Use top 50 expressed genes for correlation
    if 'total_counts' not in adata.var.columns:
         # summing columns if not present
         adata.var['total_counts'] = adata.X.sum(axis=0).A1 if hasattr(adata.X, 'toarray') else adata.X.sum(axis=0)

    top_genes = adata.var['total_counts'].sort_values(ascending=False).head(50).index
    adata_subset = adata[:, top_genes]
    
    if isinstance(adata_subset.X, np.ndarray):
        X = adata_subset.X
    else:
        X = adata_subset.X.toarray()
        
    corr_matrix = np.corrcoef(X, rowvar=False)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, xticklabels=top_genes, yticklabels=top_genes, cmap='coolwarm', center=0)
    plt.title('Gene-Gene Correlation (Top 50 Expressed)')
    plt.savefig(os.path.join(output_dir, 'specificity_gene_correlation.png'))
    plt.close()
    
    print("  > Gene correlation heatmap saved.")


def analyze_positivity(adata, output_dir, sample_tag):
    """
    Calculates responsiveness/positivity of cells to genes.
    Ref: Notebook 3.5
    """
    if 'n_cells_by_counts' not in adata.var.columns:
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
    adata.var['positivity'] = adata.var['n_cells_by_counts'] / adata.n_obs
    
    top_positive = adata.var.sort_values('positivity', ascending=False).head(50)
    top_positive[['n_cells_by_counts', 'positivity']].to_csv(os.path.join(output_dir, f"{sample_tag}_reseg_gene_positivity.csv"))
    
    plt.figure(figsize=(8, 6))
    sns.histplot(adata.var['positivity'], bins=50, kde=True, color='orange')
    plt.title("Gene Positivity Distribution (Resegmented)")
    plt.xlabel("Fraction of Positive Cells")
    plt.ylabel("Number of Genes")
    plt.savefig(os.path.join(output_dir, f"{sample_tag}_reseg_positivity_dist.png"))
    plt.close()
    
    print(f"  > Positivity analysis done.")


def analyze_diffusion(df_assigned, adata, output_dir):
    """
    Calculates distance of transcripts to cell centroids.
    Ref: Notebook 3.6
    """
    # Helper: Find x/y cols
    if 'x_global_px' in df_assigned.columns:
        x_col, y_col = 'x_global_px', 'y_global_px'
    elif 'global_x' in df_assigned.columns:
        x_col, y_col = 'global_x', 'global_y'
    elif 'x_location' in df_assigned.columns:
         x_col, y_col = 'x_location', 'y_location'
    else:
         print("  > Skipping diffusion: cannot find coordinate columns in transcripts.")
         return

    # Check for centroids in adata
    if 'x_centroid' not in adata.obs.columns:
        print("  > Skipping diffusion: x_centroid not in adata.obs")
        return
        
    centroid_map_x = adata.obs['x_centroid'].to_dict()
    centroid_map_y = adata.obs['y_centroid'].to_dict()
    
    # We assume df_assigned has 'cell_id_reseg' from Step 3
    if 'cell_id_reseg' not in df_assigned.columns:
         print("  > Warning: 'cell_id_reseg' col not found. Using 'cell_id' if available or skipping.")
         if 'cell_id' in df_assigned.columns:
             cell_col = 'cell_id'
         else:
             return
    else:
        cell_col = 'cell_id_reseg'

    # Map centroids
    # Ensure indices match types (int vs str)
    # adata index usually string, map keys based on that.
    # Check if map keys need conversion
    first_idx = next(iter(centroid_map_x))
    if isinstance(first_idx, str) and pd.api.types.is_numeric_dtype(df_assigned[cell_col]):
         # Convert df col to string
         df_assigned['cell_temp_idx'] = df_assigned[cell_col].astype(str)
         map_col = 'cell_temp_idx'
    else:
         map_col = cell_col

    df_assigned['cell_x'] = df_assigned[map_col].map(centroid_map_x)
    df_assigned['cell_y'] = df_assigned[map_col].map(centroid_map_y)
    
    # Drop rows where cell_x is NaN (unassigned or not in adata)
    df_assigned = df_assigned.dropna(subset=['cell_x', 'cell_y'])
    
    y_vals = df_assigned[y_col].values
    x_vals = df_assigned[x_col].values
    
    dy = y_vals - df_assigned['cell_y'].values
    dx = x_vals - df_assigned['cell_x'].values
    
    df_assigned['dist_to_centroid'] = np.sqrt(dy**2 + dx**2)
    
    dist_stats = df_assigned['dist_to_centroid'].describe()
    dist_stats.to_csv(os.path.join(output_dir, 'reseg_diffusion_stats.csv'))
    
    if len(df_assigned) > 10000:
         plot_data = df_assigned['dist_to_centroid'].sample(10000)
    else:
         plot_data = df_assigned['dist_to_centroid']
         
    plt.figure()
    sns.ecdfplot(data=plot_data)
    plt.title("Distance to Centroid CDF")
    plt.xlabel("Distance")
    plt.savefig(os.path.join(output_dir, 'reseg_diffusion_cdf.png'))
    plt.close()
    
    print("  > Diffusion analysis completed.")

