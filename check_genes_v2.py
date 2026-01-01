
import scanpy as sc
import os

files = [
    "./output/human_brain/human_brain_step1_preprocessed.h5ad",
    "./output/h_breast_1/h_breast_1_step1_preprocessed.h5ad"
]

for f in files:
    if os.path.exists(f):
        print(f"Checking {f}...")
        try:
            adata = sc.read_h5ad(f)
            print(f"  Var Names (Index): {adata.var_names[:5].tolist()}")
            print(f"  Var Columns: {adata.var.columns.tolist()}")
            
            # Check if symbols are in a column
            if 'feature_name' in adata.var.columns:
                print(f"  'feature_name' example: {adata.var['feature_name'][:5].tolist()}")
                
            if 'gene_ids' in adata.var.columns:
                print(f"  'gene_ids' example: {adata.var['gene_ids'][:5].tolist()}")

        except Exception as e:
            print(f"  Error reading {f}: {e}")
