import scanpy as sc
import visualize_results as vis
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    datasets = [
        ("human_brain", Path("output/human_brain")),
        ("h_breast_1", Path("output/h_breast_1"))
    ]
    
    for tag, output_dir in datasets:
        logger.info(f"Processing {tag}...")
        adata_path = output_dir / f"{tag}_step1_preprocessed.h5ad"
        
        if not adata_path.exists():
            logger.warning(f"File not found: {adata_path}")
            continue
            
        logger.info(f"Loading {adata_path}...")
        adata = sc.read_h5ad(adata_path)
        
        logger.info("Calling visualize_marker_genes_violin directly...")
        try:
            vis.visualize_marker_genes_violin(adata, output_dir, tag)
            logger.info("Success!")
        except AttributeError:
            logger.error("Function visualize_marker_genes_violin NOT FOUND in module!")
        except Exception as e:
            logger.error(f"Error: {e}")

if __name__ == "__main__":
    main()
