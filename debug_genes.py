
import scanpy as sc
import warnings

warnings.filterwarnings("ignore")

def check_genes():
    print("Loading Xenium data...")
    try:
        adata_x = sc.read_h5ad("./data/unprocessed_adata/human_brain.h5ad")
        print(f"Xenium var names (first 10): {list(adata_x.var_names[:10])}")
    except:
        print("Xenium data load failed")

    print("\nLoading Reference data...")
    try:
        adata_ref = sc.read_h5ad("data/reference/reference_brain.h5ad")
        print(f"Reference var names (first 10): {list(adata_ref.var_names[:10])}")
        print(f"Reference var columns: {list(adata_ref.var.columns)}")
        
        # Check intersection
        common = adata_x.var_names.intersection(adata_ref.var_names)
        print(f"\nCommon genes: {len(common)}")
        
        if len(common) == 0:
            print("Trying case insensitive match...")
            x_upper = [x.upper() for x in adata_x.var_names]
            ref_upper = [x.upper() for x in adata_ref.var_names]
            common_upper = set(x_upper).intersection(ref_upper)
            print(f"Common genes (case insensitive): {len(common_upper)}")
            
            # Check for Ensembl IDs (ENSG...)
            if any(x.startswith("ENSG") for x in adata_ref.var_names):
                print("Reference uses Ensembl IDs.")
            
            # Check for feature_name
            if 'feature_name' in adata_ref.var.columns:
                 print(f"feature_name first 10: {list(adata_ref.var['feature_name'][:10])}")
            if 'gene_symbol' in adata_ref.var.columns:
                 print(f"gene_symbol first 10: {list(adata_ref.var['gene_symbol'][:10])}")

    except Exception as e:
        print(f"Reference data load failed: {e}")

if __name__ == "__main__":
    check_genes()
