import scanpy as sc
from pathlib import Path

adata_path = Path("output/human_brain/human_brain_step1_preprocessed.h5ad")
if adata_path.exists():
    adata = sc.read_h5ad(adata_path)
    print(f"Obs columns: {adata.obs.columns.tolist()}")
    if 'leiden' in adata.obs:
        print("Leiden found!")
    else:
        print("Leiden NOT found.")
else:
    print("File not found.")
