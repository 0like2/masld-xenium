
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
            genes = adata.var_names[:20].tolist()
            print(f"  First 20 genes: {genes}")
            
            # Check for specific markers
            check_list = ['GFAP', 'Gfap', 'Gapdh', 'GAPDH', 'ESR1', 'Esr1']
            found = [g for g in check_list if g in adata.var_names]
            print(f"  Found markers: {found}")
        except Exception as e:
            print(f"  Error reading {f}: {e}")
