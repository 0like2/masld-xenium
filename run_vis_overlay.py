

import visualize_results as vis
from pathlib import Path
import logging
import scanpy as sc

logging.basicConfig(level=logging.INFO)

# Config
# xenium_input = Path("/data1/project/20rak/masld_xenium/human_brain") # Raw input path (Not found)
output_dir = Path("./output/human_brain")
sample_tag = "human_brain"
adata_path = Path("./output/human_brain/human_brain_step1_preprocessed.h5ad")

print(f"Loading AnnData from {adata_path} for overlay fallback...")
try:
    adata = sc.read_h5ad(adata_path)
    # Call with adata for fallback
    vis.visualize_segmentation_overlay(None, output_dir, sample_tag, adata=adata)
except Exception as e:
    print(f"Error: {e}")

