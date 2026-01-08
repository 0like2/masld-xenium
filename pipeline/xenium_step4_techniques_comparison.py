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

def run_step4(config, adata_path, transcripts_path, output_dir, original_adata_path=None, original_transcripts_path=None):
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
        adata_reseg = sc.read_h5ad(adata_path)
    except Exception as e:
        print(f"Error loading AnnData: {e}")
        return

    # Load Original Data (for Comparison)
    adata_orig = None
    if original_adata_path and os.path.exists(original_adata_path):
        print(f"Loading Original AnnData (for comparison) from {original_adata_path}...")
        try:
            adata_orig = sc.read_h5ad(original_adata_path)
        except Exception as e:
            print(f"Error loading Original AnnData: {e}")

    # Load Resegmented Transcripts (Needed for Diffusion Analysis)
    # We expect 'transcripts_path' to point to the resegmented csv from Step 3
    print(f"Loading Resegmented Transcripts from {transcripts_path}...")
    try:
        df_assigned = pd.read_csv(transcripts_path)
    except Exception as e:
        print(f"Error loading transcripts (needed for diffusion): {e}")
        # We might proceed without diffusion if this fails, but better to return or skip diffusion
        df_assigned = None

    # Load Original Transcripts (for Diffusion Comparison)
    df_assigned_orig = None
    if original_transcripts_path and os.path.exists(original_transcripts_path):
        print(f"Loading Original Transcripts (for comparison) from {original_transcripts_path}...")
        try:
            df_assigned_orig = pd.read_csv(original_transcripts_path)
        except Exception as e:
            print(f"Error loading Original Transcripts: {e}")
    
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
           analyze_efficiency(adata_reseg, output_dir, adata_orig)
        except Exception as e:
            print(f"Error in Efficiency Analysis: {e}")

    # 4. Specificity Analysis (3.4)
    if run_specificity:
        print("Running Specificity Analysis (NMP Proxy)...")
        try:
            analyze_specificity(adata_reseg, comp_config, figs_dir, adata_orig)
        except Exception as e:
            print(f"Error in Specificity Analysis: {e}")
            
    # 5. Positivity Analysis (3.5)
    if run_positivity:
        print("Running Positivity Analysis (Validation)...")
        try:
            analyze_positivity(adata_reseg, figs_dir, sample_tag, adata_orig)
        except Exception as e:
            print(f"Error in Positivity Analysis: {e}")

    # 6. Diffusion Analysis (3.6)
    if run_diffusion and df_assigned is not None:
        print("Running Diffusion Analysis (Validation)...")
        try:
            analyze_diffusion(df_assigned, adata_reseg, figs_dir, df_assigned_orig, adata_orig)
        except Exception as e:
            print(f"Error in Diffusion Analysis: {e}")

    print("Step 4: Techniques Comparison & Validation Completed.")


def analyze_efficiency(adata_reseg, output_dir, adata_orig=None):
    """
    Calculates and plots efficiency metrics (Transcripts per cell, Genes per cell).
    Ref: Notebook 3.3
    If adata_orig is provided, plots a comparison.
    """
    # Helper to prep adata
    def prep_qc(ad):
        if 'total_counts' not in ad.obs.columns or 'n_genes_by_counts' not in ad.obs.columns:
            sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)
            
    prep_qc(adata_reseg)
    
    data_list = []
    # Resegmented Data
    df_reseg = pd.DataFrame({
        'Transcripts per Cell': adata_reseg.obs['total_counts'],
        'Genes per Cell': adata_reseg.obs['n_genes_by_counts'],
        'Dataset': 'Resegmented (Step 3)'
    })
    data_list.append(df_reseg)
    
    # Original Data
    if adata_orig is not None:
        prep_qc(adata_orig)
        df_orig = pd.DataFrame({
            'Transcripts per Cell': adata_orig.obs['total_counts'],
            'Genes per Cell': adata_orig.obs['n_genes_by_counts'],
            'Dataset': 'Original (Step 0/1)'
        })
        data_list.append(df_orig)
        
    df_plot = pd.concat(data_list)
        
    # Plot 1: Transcripts per Cell
    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Transcripts per Cell', hue='Dataset', kde=True, bins=50, element="step")
    plt.title('Comparison: Transcripts per Cell')
    plt.xlabel('Transcripts per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_transcripts_per_cell_comparison.png'))
    plt.close()

    # Plot 2: Genes per Cell
    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Genes per Cell', hue='Dataset', kde=True, bins=50, element="step")
    plt.title('Comparison: Genes Detected per Cell')
    plt.xlabel('Genes per Cell')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_dir, 'efficiency_genes_per_cell_comparison.png'))
    plt.close()
    
    # Stats
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


def analyze_specificity(adata_reseg, config, output_dir, adata_orig=None):
    """
    Calculates Negative Marker Purity (NMP) if reference is provided, otherwise proxy via Gene Correlation.
    Ref: Notebook 3.4
    Supports comparison if adata_orig is provided.
    """
    datasets = [('Resegmented', adata_reseg)]
    if adata_orig is not None:
        datasets.append(('Original', adata_orig))

    # 1. Authentic NMP Analysis (if reference available)
    sc_ref_path = config.get('sc_reference_path')
    if sc_ref_path and os.path.exists(sc_ref_path):
        print(f"  > Reference scRNAseq found at {sc_ref_path}. Calculating Negative Marker Purity (NMP)...")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)
            
            with open(os.path.join(output_dir, 'specificity_nmp_score.txt'), 'w') as f:
                f.write(f"Reference: {sc_ref_path}\n")
                
                for label, ad in datasets:
                    nmp_score = calculating.negative_marker_purity_coexpression(ad, adata_sc, pipeline_output=True)
                    print(f"  > [{label}] NMP Score: {nmp_score}")
                    f.write(f"[{label}] Negative Marker Purity (NMP) Score: {nmp_score}\n")
                
        except Exception as e:
            print(f"  > Error calculating NMP: {e}. Proceeding with correlation proxy.")
    else:
        print("  > Using Gene Correlation proxy (Reference scRNAseq not available for NMP).")

    # 2. Proxy/Visualization: Gene-Gene Correlation
    print("  > Calculating Co-expression (Gene-Gene correlation)...")
    
    # Use top 50 expressed genes from Resegmented data for consistent comparison
    if 'total_counts' not in adata_reseg.var.columns:
         adata_reseg.var['total_counts'] = adata_reseg.X.sum(axis=0).A1 if hasattr(adata_reseg.X, 'toarray') else adata_reseg.X.sum(axis=0)

    top_genes = adata_reseg.var['total_counts'].sort_values(ascending=False).head(50).index

    for label, ad in datasets:
        # subset
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
    """
    Calculates responsiveness/positivity of cells to genes.
    Ref: Notebook 3.5
    Supports comparison.
    """
    datasets = [('Resegmented', adata_reseg)]
    if adata_orig is not None:
        datasets.append(('Original', adata_orig))
        
    plot_data = []

    for label, ad in datasets:
        if 'n_cells_by_counts' not in ad.var.columns:
            sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)
            
        ad.var['positivity'] = ad.var['n_cells_by_counts'] / ad.n_obs
        
        # Save CSV for top positive
        top_positive = ad.var.sort_values('positivity', ascending=False).head(50)
        safe_label = label.lower().replace(' ', '_')
        top_positive[['n_cells_by_counts', 'positivity']].to_csv(os.path.join(output_dir, f"{sample_tag}_{safe_label}_gene_positivity.csv"))
        
        # Prepare for plot
        df = pd.DataFrame({'Positivity': ad.var['positivity'], 'Dataset': label})
        plot_data.append(df)
        
    df_plot = pd.concat(plot_data)

    plt.figure(figsize=(8, 6))
    sns.histplot(data=df_plot, x='Positivity', hue='Dataset', bins=50, kde=True, element="step")
    plt.title("Gene Positivity Distribution (Comparison)")
    plt.xlabel("Fraction of Positive Cells")
    plt.ylabel("Number of Genes")
    plt.savefig(os.path.join(output_dir, f"{sample_tag}_positivity_dist_comparison.png"))
    plt.close()
    
    print(f"  > Positivity analysis done.")


def analyze_diffusion(df_reseg, adata_reseg, output_dir, df_orig=None, adata_orig=None):
    """
    Calculates distance of transcripts to cell centroids.
    Ref: Notebook 3.6
    Supports comparison.
    """
    
    def calculate_distances(df, adata, label):
        # Helper: Find x/y cols
        if 'x_global_px' in df.columns:
            x_col, y_col = 'x_global_px', 'y_global_px'
        elif 'global_x' in df.columns:
            x_col, y_col = 'global_x', 'global_y'
        elif 'x_location' in df.columns:
             x_col, y_col = 'x_location', 'y_location'
        else:
             print(f"  > [{label}] Skipping diffusion: cannot find coordinate columns in transcripts.")
             return None

        # Check for centroids in adata
        if 'x_centroid' not in adata.obs.columns:
            print(f"  > [{label}] Skipping diffusion: x_centroid not in adata.obs")
            return None
            
        centroid_map_x = adata.obs['x_centroid'].to_dict()
        centroid_map_y = adata.obs['y_centroid'].to_dict()
        
        # Determine cell ID column
        cell_col = None
        if 'cell_id_reseg' in df.columns:
            cell_col = 'cell_id_reseg'
        elif 'cell_id' in df.columns:
            cell_col = 'cell_id'
        else:
             print(f"  > [{label}] Skipping: cell id column not found.")
             return None

        # Map centroids
        # Ensure indices match types (int vs str)
        first_idx = next(iter(centroid_map_x))
        # Create a copy to avoid SettingWithCopy
        df_proc = df.copy()
        
        map_col = cell_col
        if isinstance(first_idx, str) and pd.api.types.is_numeric_dtype(df_proc[cell_col]):
             df_proc['cell_temp_idx'] = df_proc[cell_col].astype(str)
             map_col = 'cell_temp_idx'

        df_proc['cell_x'] = df_proc[map_col].map(centroid_map_x)
        df_proc['cell_y'] = df_proc[map_col].map(centroid_map_y)
        
        # Drop rows where cell_x is NaN (unassigned or not in adata)
        df_proc = df_proc.dropna(subset=['cell_x', 'cell_y'])
        
        y_vals = df_proc[y_col].values
        x_vals = df_proc[x_col].values
        
        dy = y_vals - df_proc['cell_y'].values
        dx = x_vals - df_proc['cell_x'].values
        
        dists = np.sqrt(dy**2 + dx**2)
        return pd.DataFrame({'Distance': dists, 'Dataset': label})

    # Execute
    data_list = []
    
    # Resegmented
    df_dists_reseg = calculate_distances(df_reseg, adata_reseg, 'Resegmented')
    if df_dists_reseg is not None:
        data_list.append(df_dists_reseg)
        # Stats
        df_dists_reseg['Distance'].describe().to_csv(os.path.join(output_dir, 'reseg_diffusion_stats.csv'))
        
    # Original
    if df_orig is not None and adata_orig is not None:
        df_dists_orig = calculate_distances(df_orig, adata_orig, 'Original')
        if df_dists_orig is not None:
            data_list.append(df_dists_orig)
            df_dists_orig['Distance'].describe().to_csv(os.path.join(output_dir, 'original_diffusion_stats.csv'))
            
    if not data_list:
        return

    # Compare Plot (CDF)
    df_plot = pd.concat(data_list)
    
    # Subsample for plotting efficiency if too large
    if len(df_plot) > 50000:
        df_plot = df_plot.sample(50000)

    plt.figure(figsize=(8, 6))
    sns.ecdfplot(data=df_plot, x='Distance', hue='Dataset')
    plt.title("Distance to Centroid CDF (Comparison)")
    plt.xlabel("Distance")
    plt.savefig(os.path.join(output_dir, 'diffusion_cdf_comparison.png'))
    plt.close()
    
    print("  > Diffusion analysis completed.")

