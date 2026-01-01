
import scanpy as sc
from pathlib import Path
import visualize_results as vis
import logging

# Setup basic logging
logging.basicConfig(level=logging.INFO)

# Config
datasets = ["human_brain", "h_breast_1", "h_breast_2"]

def main():
    for sample_tag in datasets:
        print(f"\n{'='*40}")
        print(f"Processing {sample_tag}...")
        print(f"{'='*40}")
        
        output_dir = Path(f"./output/{sample_tag}")
        adata_path = output_dir / f"{sample_tag}_step1_preprocessed.h5ad"
        
        if not adata_path.exists():
            print(f"File not found: {adata_path}")
            continue

        print(f"Loading {adata_path}...")
        try:
            adata = sc.read_h5ad(adata_path)
        except Exception as e:
            print(f"Failed to load: {e}")
            continue

        print("Checking for Ground Truth labels...")
        if 'ground_truth_celltype' in adata.obs:
            print("✓ Labels found.")
            print(adata.obs['ground_truth_celltype'].value_counts().head())
            
            print("Generating Comparison Plot...")
            vis.visualize_gt_comparison(adata, output_dir, sample_tag)
        else:
            print("✗ No ground truth labels found. Skipping plot.")

    print("\nAll datasets processed.")

if __name__ == "__main__":
    main()
