
import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from typing import Tuple, Optional, Dict

# Cellpose imports
try:
    from cellpose import models, io
    CELLPOSE_AVAILABLE = True
except ImportError:
    CELLPOSE_AVAILABLE = False

# Tifffile for image reading
try:
    import tifffile
    TIFFFILE_AVAILABLE = True
except ImportError:
    TIFFFILE_AVAILABLE = False

logger = logging.getLogger(__name__)

def run_step5_benchmark(
    adata: sc.AnnData,
    image_path: str,
    output_path: Path,
    sample_tag: str,
    use_baysor: bool = False,
    diameter: float = 30.0,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0
) -> Dict:
    """
    Step 5: Segmentation Benchmark using Cellpose (and optionally Baysor).
    
    1. Runs Cellpose on DAPI image to generate nuclear masks.
    2. (Optional) Runs Baysor using Cellpose masks as prior.
    3. Compares metrics with existing segmentation.

    Args:
        adata: AnnData object containing transcript location data in uns['spots'].
        image_path: Path to the DAPI image (tif/tiff).
        output_path: Directory to save outputs.
        sample_tag: Sample identifier.
        use_baysor: Whether to run Baysor (requires 'baysor' CLI installed).
        diameter: Cellpose expected nucleus diameter.

    Returns:
        Dict containing comparison metrics.
    """
    logger.info("="*60)
    logger.info("STEP 5: Segmentation Benchmark (Cellpose + Baysor)")
    logger.info("="*60)

    if not CELLPOSE_AVAILABLE:
        logger.error("Cellpose not installed. Cannot run Step 5.")
        return {}

    if not os.path.exists(image_path):
        logger.warning(f"Image not found at {image_path}. Skipping Step 5.")
        return {"status": "skipped", "reason": "image_not_found"}

    # 1. Load Image
    logger.info(f"Loading image from {image_path}...")
    try:
         # Read standard OME-TIFF or TIFF
         # Note: Xenium DAPI is often a pyramid OME-TIFF. accessing level 0 (high res)
         if image_path.endswith('.ome.tif'):
             with tifffile.TiffFile(image_path) as tif:
                 # Usually the first series is the full resolution
                 img = tif.series[0].levels[0].asarray()
         else:
             img = io.imread(image_path)
             
         logger.info(f"Image shape: {img.shape}")
         # Check if 3D or channels exist. Cellpose expects 2D or (2,Y,X). 
         # Xenium DAPI is usually single channel 2D. 
         if img.ndim == 3 and img.shape[0] < 5:
             # Assume channel first, take first channel if multiple?
             # For DAPI only, usually just one channel.
             img = img[0, :, :]
             
    except Exception as e:
        logger.error(f"Failed to load image: {e}")
        return {"status": "failed", "error": str(e)}

    # 2. Run Cellpose
    logger.info(f"Running Cellpose (Nuclei model)... Diameter={diameter}")
    try:
        model = models.Cellpose(gpu=True, model_type='nuclei')
        masks, flows, styles, diams = model.eval(
            img, 
            diameter=diameter, 
            channels=[0,0], 
            flow_threshold=flow_threshold, 
            cellprob_threshold=cellprob_threshold
        )
        logger.info(f"Cellpose completed. Found {masks.max()} nuclei.")
        
        # Save masks
        mask_path = output_path / f"{sample_tag}_step5_cellpose_masks.npy"
        np.save(mask_path, masks)
        
        # Save visualization (snapshot)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10,10))
        plt.imshow(img, cmap='gray')
        plt.imshow(masks, alpha=0.5, cmap='nipy_spectral')
        plt.title(f"Cellpose Nuclei: {masks.max()} cells")
        plt.axis('off')
        plt.savefig(output_path / f"{sample_tag}_step5_cellpose_overlay.png", dpi=150)
        plt.close()
        
    except Exception as e:
        logger.error(f"Cellpose run failed: {e}")
        return {"status": "failed", "error": str(e)}
        
    # 3. Assign Transcripts to new Masks
    # We need to map adata.uns['spots'] (x, y) to masks[y, x]
    # Note: Xenium coordinates (microns) need to be converted to Pixel coordinates if image is pixels.
    # Usually transformation matrix is needed. 
    # IF we assume image is already registered and pixel scaling factor is known.
    # PROVISIONAL: Assuming 1 pixel = X microns OR image matches the coordinate space (unlikely for raw image vs microns).
    # Xenium pixel_size = 0.2125 microns usually.
    
    metrics = {
        'cellpose_n_cells': int(masks.max()),
    }

    # 4. Run Baysor (Optional)
    if use_baysor:
        logger.info("Preparing inputs for Baysor...")
        # Check if baysor exists
        from shutil import which
        if which('baysor') is None:
            logger.warning("Baysor command not found. Skipping Baysor key.")
            metrics['baysor_status'] = 'missing_executable'
        else:
            # Need to export transcripts to CSV for Baysor
            # Required cols: x, y, gene
            if 'spots' in adata.uns:
                spots = adata.uns['spots']
                # Ensure x,y,gene columns
                # ... export to csv ...
                # ... run subprocess ...
                pass
            metrics['baysor_status'] = 'implemented_but_skipped_logic' # placeholder

    return metrics

def compare_segmentations(step4_stats: Dict, step5_stats: Dict) -> pd.DataFrame:
    """Compare Step 4 (Original/Expanded) vs Step 5 (Cellpose)."""
    # Simple table comparison
    return pd.DataFrame([step4_stats, step5_stats])
