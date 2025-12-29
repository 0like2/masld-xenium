
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from xb.calculating import dispersion
import scanpy as sc

# Configure plotting style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def visualize_step2_gene_distances(output_dir, sample_tag):
    """Visualizes gene distance to nucleus distribution."""
    print(f"Visualizing Step 2 for {sample_tag}...")
    
    file_path = os.path.join(output_dir, f"{sample_tag}_step2_gene_distances.csv")
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    df = pd.read_csv(file_path, index_col=0)
    
    # Sort by mean distance
    df = df.sort_values('mean')
    
    # Create plot
    plt.figure(figsize=(15, 6))
    
    # Select top/bottom genes for cleaner visualization if too many
    n_display = 50
    if len(df) > 2 * n_display:
        plot_df = pd.concat([df.head(n_display), df.tail(n_display)])
    else:
        plot_df = df
        
    sns.barplot(x=plot_df.index, y=plot_df['mean'])
    plt.xticks(rotation=90, fontsize=8)
    plt.title(f"Average Distance to Nucleus by Gene ({sample_tag})")
    plt.ylabel("Mean Distance (µm)")
    plt.xlabel("Gene")
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"{sample_tag}_step2_gene_dist_plot.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot to {save_path}")
    plt.close()

def visualize_step4_optimal_expansion(output_dir, sample_tag):
    """Visualizes optimal expansion distance results."""
    print(f"Visualizing Step 4 for {sample_tag}...")
    
    file_path = os.path.join(output_dir, f"{sample_tag}_step4_optimal_expansion.csv")
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    df = pd.read_csv(file_path, index_col=0)
    
    # Data for plotting
    try:
        # Check available columns based on previous `view_file` output of `step4`
        # It had 'optimal_from_center', 'optimal_from_nucleus_border', 'mean_nucleus_distance', 'xenium_default'
        # The CSV is transposed, so index is the metric name, column is 0
        
        metrics = df.index.tolist()
        values = df.iloc[:, 0].tolist() # Assuming single row transposed to single col
        
        plt.figure(figsize=(10, 6))
        sns.barplot(x=metrics, y=values)
        plt.title(f"Optimal Cell Expansion Parameters ({sample_tag})")
        plt.ylabel("Distance (µm)")
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        save_path = os.path.join(output_dir, f"{sample_tag}_step4_expansion_plot.png")
        plt.savefig(save_path, dpi=300)
        print(f"Saved plot to {save_path}")
        plt.close()
    except Exception as e:
        print(f"Error plotting Step 4: {e}")

def visualize_step6_simulation(output_dir, sample_tag):
    """Visualizes simulation results: ARI vs Silhouette."""
    print(f"Visualizing Step 6 for {sample_tag}...")
    
    # Load summary if exists, else merge individual files
    summary_path = os.path.join(output_dir, f"{sample_tag}_step6_summary.csv")
    ari_path = os.path.join(output_dir, f"{sample_tag}_step6_ari_scores.csv")
    sil_path = os.path.join(output_dir, f"{sample_tag}_step6_silhouette_scores.csv")
    
    if os.path.exists(summary_path):
        df = pd.read_csv(summary_path, index_col=0)
    elif os.path.exists(ari_path) and os.path.exists(sil_path):
        ari_df = pd.read_csv(ari_path, index_col=0)
        sil_df = pd.read_csv(sil_path, index_col=0)
        df = pd.merge(ari_df, sil_df, left_index=True, right_index=True)
    else:
        print("Step 6 results not found.")
        return

    if df.empty:
        print("Step 6 dataframe is empty.")
        return

    # Add 'Recommended' flag logic similar to xenium_pipeline_main.py logger
    # Or just scatter plot ARI vs Silhouette
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x='ARI', y='Silhouette', alpha=0.6)
    
    # Highlight top candidates
    # Top 5 ARI
    top_ari = df.nlargest(5, 'ARI')
    # Top 5 Silhouette
    top_sil = df.nlargest(5, 'Silhouette')
    
    sns.scatterplot(data=top_ari, x='ARI', y='Silhouette', color='red', label='Top ARI', s=100)
    sns.scatterplot(data=top_sil, x='ARI', y='Silhouette', color='green', label='Top Silhouette', s=100)
    
    # Label top points? Might be too crowded. Just label the very best.
    best_ari = top_ari.iloc[0]
    best_sil = top_sil.iloc[0]
    
    plt.annotate(f"Best ARI\n{best_ari.name}", 
                 (best_ari['ARI'], best_ari['Silhouette']), 
                 xytext=(10, -20), textcoords='offset points', 
                 arrowprops=dict(arrowstyle="->", color='red'))
                 
    plt.annotate(f"Best Silhouette\n{best_sil.name}", 
                 (best_sil['ARI'], best_sil['Silhouette']), 
                 xytext=(-10, 20), textcoords='offset points', 
                 arrowprops=dict(arrowstyle="->", color='green'))

    plt.title(f"Preprocessing Simulation: Stability (ARI) vs Quality (Silhouette) - {sample_tag}")
    plt.xlabel("ARI (Stability vs Default)")
    plt.ylabel("Silhouette Score (Cluster Separation)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"{sample_tag}_step6_simulation_plot.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot to {save_path}")
    plt.close()

if __name__ == "__main__":
    datasets = [
        ("human_brain", "./output/human_brain"),
        ("h_breast_1", "./output/h_breast_1")
    ]
    
    for tag, path in datasets:
        if os.path.exists(path):
            visualize_step2_gene_distances(path, tag)
            visualize_step4_optimal_expansion(path, tag)
            visualize_step6_simulation(path, tag)
