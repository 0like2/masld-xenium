
import scanpy as sc
adata = sc.read_h5ad("./output/human_brain/human_brain_step1_preprocessed.h5ad")
print(f"Index: {adata.var_names[:3]}")
if 'gene_id' in adata.var.columns:
    print(f"gene_id col: {adata.var['gene_id'][:3].tolist()}")
if 'feature_name' in adata.var.columns:
    print(f"feature_name col: {adata.var['feature_name'][:3].tolist()}")
