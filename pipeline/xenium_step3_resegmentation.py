# Step 3: Resegmentation (Ref: notebooks/3_techniques_comparison/3_1, 3_2)
# Re-segments nuclei using Cellpose and re-assigns transcripts to new cells.
#
# Flow:
#   3-1. Device detection (CUDA/MPS/CPU)
#   3-2. Load DAPI image
#   3-3. Cellpose nuclei segmentation
#   3-4. Save masks → TIF
#   3-5. Map transcripts to cell masks
#   3-6. Build cell×gene matrix → AnnData
#   3-7. Extract cell centroids (regionprops)
#   3-8. Save AnnData + resegmented transcripts
#   3-9. Domain assignment (optional, JSON polygons)

import os
import numpy as np
import pandas as pd
import tifffile as tf
import scanpy as sc
import torch
from cellpose import models
from skimage.measure import label, regionprops
from skimage.segmentation import expand_labels
from tqdm import tqdm
import gc
import json
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt


def _auto_tile_size(device_index, min_tile=4096, max_tile=16384):
    """Compute safe tile size for Cellpose based on available GPU memory.

    Flow dynamics allocates (1, 2, H, W) float32 on GPU ≈ 8 bytes/pixel.
    We budget 40 % of free VRAM for this tensor + overhead.
    """
    if not torch.cuda.is_available():
        return max_tile
    free, _ = torch.cuda.mem_get_info(device_index if device_index is not None else 0)
    tile = int(np.sqrt(free * 0.4 / 8))
    tile = (tile // 256) * 256
    return max(min_tile, min(tile, max_tile))


def _run_cellpose_tiled(model, image, tile_size=None, overlap=512,
                        device_index=None, **eval_kwargs):
    """Run Cellpose with image tiling to avoid GPU OOM during flow dynamics.

    The Cellpose post-processing (``dynamics.steps_interp``) allocates a
    full-image-sized tensor on GPU.  For very large images (e.g. 40 k × 37 k)
    this easily exceeds VRAM.  This helper splits the image into overlapping
    tiles, runs Cellpose on each tile independently, and stitches the masks
    using only the *core* (non-overlap) region of each tile.
    """
    H, W = image.shape[:2]

    if tile_size is None:
        tile_size = _auto_tile_size(device_index)

    # Small enough – run directly
    if H <= tile_size and W <= tile_size:
        print(f"  Image ({H}x{W}) fits in one tile ({tile_size}x{tile_size}), "
              "running directly.")
        return model.eval(image, **eval_kwargs)

    stride = tile_size - overlap
    margin = overlap // 2

    all_masks = np.zeros((H, W), dtype=np.int32)
    label_offset = 0

    y_starts = list(range(0, H, stride))
    x_starts = list(range(0, W, stride))
    n_tiles = len(y_starts) * len(x_starts)
    print(f"  Tiling: {H}x{W} image -> {n_tiles} tiles of <={tile_size}x"
          f"{tile_size} (stride={stride}, overlap={overlap})")

    for idx, (y0, x0) in enumerate(
            [(y, x) for y in y_starts for x in x_starts], 1):
        y1 = min(y0 + tile_size, H)
        x1 = min(x0 + tile_size, W)
        tile = image[y0:y1, x0:x1]
        th, tw = tile.shape[:2]
        print(f"  Tile {idx}/{n_tiles}: [{y0}:{y1}, {x0}:{x1}] ({th}x{tw})")

        tile_masks, _, _ = model.eval(tile, **eval_kwargs)

        # Core region (exclude overlap margins; keep full extent at edges)
        cy0 = 0 if y0 == 0 else margin
        cx0 = 0 if x0 == 0 else margin
        cy1 = th if y1 >= H else th - margin
        cx1 = tw if x1 >= W else tw - margin
        if cy1 <= cy0 or cx1 <= cx0:
            del tile, tile_masks
            continue

        core = tile_masks[cy0:cy1, cx0:cx1]
        ulabels = np.unique(core)
        ulabels = ulabels[ulabels > 0]
        n_cells = len(ulabels)

        if n_cells > 0:
            lut = np.zeros(int(ulabels.max()) + 1, dtype=np.int32)
            lut[ulabels] = np.arange(label_offset + 1,
                                     label_offset + n_cells + 1,
                                     dtype=np.int32)
            all_masks[y0 + cy0:y0 + cy1, x0 + cx0:x0 + cx1] = lut[core]
            label_offset += n_cells

        del tile, tile_masks, core
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(f"    -> {n_cells} cells (running total: {label_offset})")

    print(f"  Tiling complete: {label_offset} cells from {n_tiles} tiles")
    return all_masks, None, None


def _detect_device():
    """Auto-detect best available device: pick CUDA GPU with most free memory > MPS > CPU."""
    if torch.cuda.is_available():
        n_gpus = torch.cuda.device_count()
        best_gpu = 0
        best_free = 0
        for i in range(n_gpus):
            free, total = torch.cuda.mem_get_info(i)
            name = torch.cuda.get_device_name(i)
            print(f"  GPU {i}: {name} — {free/1e9:.1f}GB free / {total/1e9:.1f}GB total")
            if free > best_free:
                best_free = free
                best_gpu = i
        torch.cuda.set_device(best_gpu)
        name = torch.cuda.get_device_name(best_gpu)
        print(f"  Selected GPU {best_gpu}: {name} ({best_free/1e9:.1f}GB free)")
        return torch.device(f"cuda:{best_gpu}"), True
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        print("  Device: MPS (Apple Silicon)")
        return torch.device("mps"), False
    print("  Device: CPU")
    return torch.device("cpu"), False

# --- DAPI Visualization Helpers ---

def _load_dapi_downsampled(dapi_image_path, max_size=2000):
    """Load a DAPI image and downsample for visualization.

    Returns (dapi_ds, step) where step is the downsample factor,
    or (None, 1) if loading fails.
    """
    try:
        dapi = tf.imread(dapi_image_path)
        if dapi.ndim == 3:
            dapi = dapi[0]  # take first channel/z-plane
        step = max(1, max(dapi.shape) // max_size)
        dapi_ds = dapi[::step, ::step]
        del dapi
        gc.collect()
        return dapi_ds, step
    except Exception as e:
        print(f"  > Could not load DAPI for visualization: {e}")
        return None, 1


def _save_dapi_overview(dapi_ds, output_dir, sample_tag):
    """Save DAPI overview image (matching notebook: plt.imshow(dapi, vmax=3000))."""
    try:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(dapi_ds, cmap='gray', vmax=np.percentile(dapi_ds, 99.5))
        ax.set_title('DAPI morphology overview')
        ax.axis('off')
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step3_dapi_overview.png")
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  > DAPI overview saved: {path}")
    except Exception as e:
        print(f"  > Warning: DAPI overview failed: {e}")
        plt.close('all')


def _save_dapi_mask_overlay(dapi_ds, masks, nuclei_masks, step, output_dir, sample_tag):
    """DAPI background + colored segmentation mask contours overlay.

    Shows 2-3 panels: DAPI alone | DAPI + expanded masks | DAPI + nuclei masks (if available).
    Matching notebook pattern: imshow(dapi) + colored contour overlay.
    """
    from skimage.color import label2rgb
    from skimage.segmentation import find_boundaries

    try:
        masks_ds = masks[::step, ::step]
        n_panels = 3 if nuclei_masks is not None else 2
        fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 7))

        # Panel 1: Raw DAPI
        vmax = np.percentile(dapi_ds, 99.5)
        axes[0].imshow(dapi_ds, cmap='gray', vmax=vmax)
        axes[0].set_title('DAPI')
        axes[0].axis('off')

        # Panel 2: DAPI + expanded mask overlay
        # Create colored boundaries on DAPI
        boundaries = find_boundaries(masks_ds, mode='thick')
        dapi_rgb = np.stack([dapi_ds / vmax] * 3, axis=-1).clip(0, 1)
        # Color boundaries in cyan
        dapi_overlay = dapi_rgb.copy()
        dapi_overlay[boundaries, 0] = 0
        dapi_overlay[boundaries, 1] = 1
        dapi_overlay[boundaries, 2] = 1
        axes[1].imshow(dapi_overlay)
        axes[1].set_title(f'DAPI + expanded masks ({masks.max():,} cells)')
        axes[1].axis('off')

        # Panel 3: DAPI + nuclei masks overlay
        if nuclei_masks is not None:
            nuc_ds = nuclei_masks[::step, ::step]
            nuc_boundaries = find_boundaries(nuc_ds, mode='thick')
            nuc_overlay = dapi_rgb.copy()
            nuc_overlay[nuc_boundaries, 0] = 1
            nuc_overlay[nuc_boundaries, 1] = 0.3
            nuc_overlay[nuc_boundaries, 2] = 0
            axes[2].imshow(nuc_overlay)
            axes[2].set_title(f'DAPI + nuclei ({nuclei_masks.max():,} nuclei)')
            axes[2].axis('off')

        fig.suptitle('Segmentation overlay on DAPI', fontsize=14)
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step3_dapi_mask_overlay.png")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  > DAPI + mask overlay saved: {path}")
    except Exception as e:
        print(f"  > Warning: DAPI + mask overlay failed: {e}")
        plt.close('all')


def _save_dapi_zoomed_roi(dapi_ds, masks, nuclei_masks, step, output_dir, sample_tag,
                          roi_fraction=0.15):
    """Zoomed ROI showing DAPI + mask detail at full resolution.

    Picks a region near the median nucleus position and shows:
    - DAPI + nuclei boundaries (orange)
    - DAPI + expanded boundaries (cyan)
    - DAPI + label2rgb mask (semi-transparent)
    """
    from skimage.color import label2rgb
    from skimage.segmentation import find_boundaries

    try:
        masks_ds = masks[::step, ::step]
        h, w = masks_ds.shape

        # Find center of a dense region
        nuc_ds = nuclei_masks[::step, ::step] if nuclei_masks is not None else masks_ds
        nuc_present = np.argwhere(nuc_ds > 0)
        if len(nuc_present) == 0:
            print("  > No nuclei found for ROI selection; skipping zoomed view.")
            return
        cy, cx = nuc_present[len(nuc_present) // 2]

        # Crop
        rh, rw = int(h * roi_fraction), int(w * roi_fraction)
        r0, r1 = max(0, cy - rh // 2), min(h, cy + rh // 2)
        c0, c1 = max(0, cx - rw // 2), min(w, cx + rw // 2)

        dapi_crop = dapi_ds[r0:r1, c0:c1]
        mask_crop = masks_ds[r0:r1, c0:c1]

        vmax = np.percentile(dapi_ds, 99.5)
        dapi_rgb = np.stack([dapi_crop / vmax] * 3, axis=-1).clip(0, 1)

        fig, axes = plt.subplots(1, 3, figsize=(21, 7))

        # Panel 1: DAPI zoomed
        axes[0].imshow(dapi_crop, cmap='gray', vmax=vmax)
        axes[0].set_title('DAPI (zoomed ROI)')
        axes[0].axis('off')

        # Panel 2: DAPI + mask boundaries
        boundaries = find_boundaries(mask_crop, mode='thick')
        overlay = dapi_rgb.copy()
        overlay[boundaries, 0] = 0
        overlay[boundaries, 1] = 1
        overlay[boundaries, 2] = 1
        if nuclei_masks is not None:
            nuc_crop = nuc_ds[r0:r1, c0:c1]
            nuc_bd = find_boundaries(nuc_crop, mode='thick')
            overlay[nuc_bd, 0] = 1
            overlay[nuc_bd, 1] = 0.5
            overlay[nuc_bd, 2] = 0
        axes[1].imshow(overlay)
        axes[1].set_title('DAPI + boundaries (cyan=expanded, orange=nuclei)')
        axes[1].axis('off')

        # Panel 3: DAPI + colored mask (semi-transparent)
        mask_rgb = label2rgb(mask_crop, bg_label=0)
        blended = 0.5 * dapi_rgb + 0.5 * mask_rgb
        axes[2].imshow(blended.clip(0, 1))
        axes[2].set_title('DAPI + colored masks (blended)')
        axes[2].axis('off')

        fig.suptitle('Zoomed ROI — Segmentation on DAPI', fontsize=14)
        fig.tight_layout()
        path = os.path.join(output_dir, f"{sample_tag}_step3_dapi_zoomed_roi.png")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  > DAPI zoomed ROI saved: {path}")
    except Exception as e:
        print(f"  > Warning: DAPI zoomed ROI failed: {e}")
        plt.close('all')


def run_step3(config, dapi_image_path, transcripts_path, output_dir):
    """
    Step 3: Resegmentation (Cellpose)
    Resegments nuclei using Cellpose and re-assigns transcripts.
    """
    print("------------------------------------------------")
    print("Starting Step 3: Resegmentation (Cellpose)")
    print("------------------------------------------------")

    reseg_config = config.get('resegmentation', {})
    if not reseg_config.get('run_resegmentation', False):
        print("Resegmentation is disabled in config. Skipping.")
        return

    # --- 3-2 / 3-3 / 3-4. Cellpose segmentation (with mask caching) ---

    sample_tag = config.get('sample_tag', 'sample')
    mask_out_path = os.path.join(output_dir, f"{sample_tag}_step3_resegmented_masks.tif")
    nuclei_mask_path = os.path.join(output_dir, f"{sample_tag}_step3_nuclei_masks.tif")
    expansion_distance = reseg_config.get('expansion_distance', 400)

    # Check for cached masks — skip Cellpose + expansion if already on disk
    if os.path.exists(mask_out_path):
        print(f"  [CACHE HIT] Loading cached masks from {mask_out_path}")
        masks = tf.imread(mask_out_path)
        nuclei_masks = None
        if expansion_distance > 0 and os.path.exists(nuclei_mask_path):
            nuclei_masks = tf.imread(nuclei_mask_path)
            print(f"  [CACHE HIT] Loaded nuclei masks from {nuclei_mask_path}")
        print(f"  > Expanded masks: {masks.max()} labels, shape {masks.shape}")
        if nuclei_masks is not None:
            print(f"  > Nuclei masks: {nuclei_masks.max()} labels")
    else:
        # --- 3-2. Load DAPI image ---
        print(f"Loading DAPI image from {dapi_image_path}...")
        if not os.path.exists(dapi_image_path):
            print(f"Error: DAPI image not found at {dapi_image_path}")
            return

        try:
            dapi_image = tf.imread(dapi_image_path)
        except Exception as e:
            print(f"Error reading DAPI image: {e}")
            return

        print(f"DAPI Image Shape: {dapi_image.shape}")

        # --- 3-3. Cellpose nuclei segmentation ---
        cp_params = reseg_config.get('cellpose', {})
        diameter = cp_params.get('diameter', None)

        device, use_gpu = _detect_device()

        batch_size = cp_params.get('batch_size', 8)
        if use_gpu and torch.cuda.is_available():
            free_mem, _ = torch.cuda.mem_get_info(device.index or 0)
            free_gb = free_mem / 1e9
            if free_gb >= 30:
                batch_size = 32
            elif free_gb >= 16:
                batch_size = 16
            else:
                batch_size = 8
            print(f"  GPU memory available: {free_gb:.1f}GB → batch_size={batch_size}")

        print(f"Initializing Cellpose (device={device}, model='nuclei')...")
        model = models.CellposeModel(gpu=use_gpu, pretrained_model='nuclei', device=device)

        print(f"Running Cellpose eval (batch_size={batch_size})...")
        tile_size = cp_params.get('tile_size', None)
        masks, flows, styles = _run_cellpose_tiled(
            model, dapi_image,
            tile_size=tile_size,
            overlap=cp_params.get('tile_overlap', 512),
            device_index=device.index if hasattr(device, 'index') else None,
            diameter=diameter, batch_size=batch_size,
        )

        del dapi_image, flows, styles, model
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # --- 3-3b. Label expansion ---
        masks = label(masks)
        print(f"  > Labeled {masks.max()} nuclei objects")

        nuclei_masks = masks  # Keep original nuclei masks for dual assignment
        if expansion_distance > 0:
            print(f"Expanding nuclei masks by {expansion_distance} pixels...")
            masks = expand_labels(nuclei_masks, distance=expansion_distance)
            n_before = len(np.unique(nuclei_masks)) - 1
            n_after = len(np.unique(masks)) - 1
            print(f"  > Nuclei labels: {n_before}, Expanded labels: {n_after}")
        else:
            nuclei_masks = None

        # --- 3-4. Save masks → TIF (checkpoint for future reruns) ---
        print(f"Saving Cellpose masks to TIF: {mask_out_path}")
        tf.imwrite(mask_out_path, masks)

    if nuclei_masks is not None and not os.path.exists(nuclei_mask_path):
        nuclei_mask_path = os.path.join(output_dir, f"{sample_tag}_step3_nuclei_masks.tif")
        print(f"Saving nuclei-only masks to TIF: {nuclei_mask_path}")
        tf.imwrite(nuclei_mask_path, nuclei_masks)

    # --- 3-4b. QC visualization of masks ---
    try:
        from skimage.color import label2rgb
        sample_tag = config.get('sample_tag', 'sample')
        qc_path = os.path.join(output_dir, f"{sample_tag}_step3_mask_qc.png")
        print(f"Generating mask QC plot: {qc_path}")

        fig, axes = plt.subplots(1, 2 if nuclei_masks is not None else 1,
                                 figsize=(12 if nuclei_masks is not None else 6, 6))
        if nuclei_masks is None:
            axes = [axes]

        # Nuclei / primary masks (downsampled for speed)
        step = max(1, masks.shape[0] // 2000)
        axes[0].imshow(label2rgb(masks[::step, ::step], bg_label=0), interpolation='nearest')
        axes[0].set_title(f'Expanded masks ({masks.max()} labels)' if nuclei_masks is not None
                          else f'Nuclei masks ({masks.max()} labels)')
        axes[0].axis('off')

        if nuclei_masks is not None:
            # For nuclei masks: use a zoomed-in crop for visibility (nuclei are
            # too small to survive aggressive downsampling of the full image)
            nuc_ds = nuclei_masks[::step, ::step]
            h_ds, w_ds = nuc_ds.shape
            # Find a region with nuclei and crop to ~25% of the image
            nuc_present = np.argwhere(nuc_ds > 0)
            if len(nuc_present) > 0:
                cy, cx = nuc_present[len(nuc_present) // 2]  # median nucleus position
                crop_h, crop_w = h_ds // 4, w_ds // 4
                r0, r1 = max(0, cy - crop_h // 2), min(h_ds, cy + crop_h // 2)
                c0, c1 = max(0, cx - crop_w // 2), min(w_ds, cx + crop_w // 2)
                nuc_crop = nuc_ds[r0:r1, c0:c1]
            else:
                nuc_crop = nuc_ds
            axes[1].imshow(label2rgb(nuc_crop, bg_label=0), interpolation='nearest')
            axes[1].set_title(f'Nuclei masks — zoomed ({nuclei_masks.max()} labels)')
            axes[1].axis('off')

        plt.tight_layout()
        plt.savefig(qc_path, dpi=150)
        plt.close()
        print("  > Mask QC plot saved.")
    except Exception as e:
        print(f"  > Warning: Could not generate mask QC plot: {e}")

    # --- 3-4c. DAPI + segmentation overlay visualizations ---
    print("Generating DAPI overlay visualizations...")
    dapi_ds, dapi_step = _load_dapi_downsampled(dapi_image_path)
    if dapi_ds is not None:
        sample_tag_dapi = config.get('sample_tag', 'sample')
        _save_dapi_overview(dapi_ds, output_dir, sample_tag_dapi)
        _save_dapi_mask_overlay(dapi_ds, masks, nuclei_masks, dapi_step, output_dir, sample_tag_dapi)
        _save_dapi_zoomed_roi(dapi_ds, masks, nuclei_masks, dapi_step, output_dir, sample_tag_dapi)
        del dapi_ds
        gc.collect()
    else:
        print("  > Skipping DAPI overlays (image not available).")

    # --- 3-5. Map transcripts to cell masks ---

    print("Loading Transcripts...")
    try:
        df_transcripts = pd.read_csv(transcripts_path)
    except Exception as e:
        print(f"Error loading transcripts: {e}")
        return

    if 'x_global_px' in df_transcripts.columns:
        x_col, y_col = 'x_global_px', 'y_global_px'
    elif 'global_x' in df_transcripts.columns:
        x_col, y_col = 'global_x', 'global_y'
    else:
        x_col, y_col = 'x_location', 'y_location'

    # CosMx: x and y coordinates are swapped compared to column names
    technology = config.get('comparison', {}).get('technology', 'xenium').lower()
    if technology == 'cosmx':
        print(f"  > CosMx detected: swapping x/y coordinate columns")
        x_col, y_col = y_col, x_col

    print(f"Mapping transcripts to cells using columns: {x_col}, {y_col}...")

    max_h, max_w = masks.shape

    # Convert µm → pixel when coordinates are in microns (x_location/y_location)
    # Pixel columns (x_global_px etc.) are already in pixel space.
    um_per_pixel_inv = reseg_config.get('um_per_pixel_inv', 4.70588)  # px/µm
    needs_um_to_px = x_col in ('x_location', 'y_location')
    if needs_um_to_px:
        print(f"  > Converting µm → px (factor={um_per_pixel_inv:.4f} px/µm)")

    try:
        y_vals = df_transcripts[y_col].values
        x_vals = df_transcripts[x_col].values
        if needs_um_to_px:
            y_coords = (y_vals * um_per_pixel_inv).astype(int)
            x_coords = (x_vals * um_per_pixel_inv).astype(int)
        else:
            y_coords = y_vals.astype(int)
            x_coords = x_vals.astype(int)
    except KeyError:
        print(f"Error: Columns {x_col}, {y_col} not found in transcripts.")
        return

    valid_mask = (y_coords >= 0) & (y_coords < max_h) & (x_coords >= 0) & (x_coords < max_w)

    if not np.any(valid_mask):
        print("Error: No transcripts fall within image bounds.")
        return

    df_valid = df_transcripts[valid_mask].copy()
    y_coords = y_coords[valid_mask]
    x_coords = x_coords[valid_mask]

    cell_labels = masks[y_coords, x_coords]
    df_valid['cell_id_reseg'] = cell_labels

    # Dual assignment: nuclei (in_cell) vs expanded (closest_cell) per notebook convention
    if nuclei_masks is not None:
        nuclei_labels = nuclei_masks[y_coords, x_coords]
        df_valid['in_cell'] = nuclei_labels
        df_valid['closest_cell'] = cell_labels
    else:
        df_valid['in_cell'] = cell_labels
        df_valid['closest_cell'] = cell_labels

    df_assigned = df_valid[df_valid['cell_id_reseg'] > 0].copy()

    print(f"  > Assigned {len(df_assigned)} transcripts to {df_assigned['cell_id_reseg'].nunique()} cells.")
    if nuclei_masks is not None:
        n_in_nuclei = (df_assigned['in_cell'] > 0).sum()
        n_in_expanded = (df_assigned['closest_cell'] > 0).sum()
        print(f"  > In nuclei: {n_in_nuclei}, In expanded cells: {n_in_expanded}")

    # C1-a: Mask + transcript scatter overlay (binary mask bg + DAPI bg versions)
    try:
        from skimage.color import label2rgb
        sample_tag_vis = config.get('sample_tag', 'sample')
        sub = df_assigned.sample(n=min(len(df_assigned), 50000), random_state=42)
        step_ds = max(1, masks.shape[0] // 2000)
        if needs_um_to_px:
            sx = sub[x_col].values * um_per_pixel_inv / step_ds
            sy = sub[y_col].values * um_per_pixel_inv / step_ds
        else:
            sx = sub[x_col].values / step_ds
            sy = sub[y_col].values / step_ds

        # Version 1: binary mask background
        fig, ax = plt.subplots(figsize=(10, 10))
        binary_mask = (masks[::step_ds, ::step_ds] > 0).astype(np.uint8)
        ax.imshow(binary_mask, cmap='gray', vmin=0, vmax=2, interpolation='nearest', alpha=0.5)
        ax.scatter(sx, sy, s=0.5, alpha=0.6, c='red', rasterized=True)
        ax.set_title(f'Mask + transcript overlay ({len(sub)} transcripts)')
        ax.axis('off')
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{sample_tag_vis}_step3_mask_transcript_overlay.png"), dpi=150)
        plt.close(fig)
        print("  > Mask + transcript overlay saved.")

        # Version 2: DAPI background (notebook-style)
        dapi_ds_tx, _ = _load_dapi_downsampled(dapi_image_path, max_size=2000)
        if dapi_ds_tx is not None:
            fig, ax = plt.subplots(figsize=(10, 10))
            vmax = np.percentile(dapi_ds_tx, 99.5)
            ax.imshow(dapi_ds_tx, cmap='gray', vmax=vmax)
            ax.scatter(sx, sy, s=0.3, alpha=0.5, c='red', rasterized=True)
            ax.set_title(f'DAPI + transcripts ({len(sub)} transcripts)')
            ax.axis('off')
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, f"{sample_tag_vis}_step3_dapi_transcript_overlay.png"), dpi=200)
            plt.close(fig)
            del dapi_ds_tx
            gc.collect()
            print("  > DAPI + transcript overlay saved.")
    except Exception as e:
        print(f"  > Warning: mask + transcript overlay failed: {e}")
        plt.close('all')

    # --- 3-6a. Extract cell centroids + morphology (regionprops) ---
    # Compute regionprops BEFORE AnnData creation so centroids are available for distance_to_centroid
    print("Calculating centroids and cell morphology...")
    try:
        # Centroids from nuclei masks (notebook convention); area/perimeter from expanded masks
        if nuclei_masks is not None:
            nuc_props = regionprops(nuclei_masks)
            centroid_dict = {p.label: p.centroid for p in nuc_props}
            exp_props = regionprops(masks)
            area_dict = {p.label: p.area for p in exp_props}
            perimeter_dict = {p.label: p.perimeter for p in exp_props}
        else:
            props = regionprops(masks)
            centroid_dict = {p.label: p.centroid for p in props}
            area_dict = {p.label: p.area for p in props}
            perimeter_dict = {p.label: p.perimeter for p in props}
    except Exception as e:
        print(f"Error computing regionprops: {e}")
        return

    # --- 3-6b. Compute per-transcript distance_to_centroid ---
    print("Computing per-transcript distance to centroid...")
    df_assigned['closest_cell_x'] = df_assigned['closest_cell'].map(
        lambda c: centroid_dict.get(c, (np.nan, np.nan))[1])
    df_assigned['closest_cell_y'] = df_assigned['closest_cell'].map(
        lambda c: centroid_dict.get(c, (np.nan, np.nan))[0])
    df_assigned['distance_to_centroid'] = np.sqrt(
        (df_assigned[x_col] - df_assigned['closest_cell_x'])**2 +
        (df_assigned[y_col] - df_assigned['closest_cell_y'])**2)
    n_with_dist = df_assigned['distance_to_centroid'].notna().sum()
    print(f"  > Computed distance_to_centroid for {n_with_dist} transcripts")

    # --- 3-6c. Build cell×gene matrix → AnnData ---
    print("Generating Cell x Gene Matrix...")
    try:
        gene_col = 'feature_name' if 'feature_name' in df_assigned.columns else 'gene'

        # Filter out control probes (NegControlProbe, NegControlCodeword, BLANK, antisense)
        # These inflate gene counts and are not biologically meaningful
        ctrl_mask = df_assigned[gene_col].str.contains(
            'NegControl|BLANK|antisense', case=False, na=False)
        n_ctrl = ctrl_mask.sum()
        if n_ctrl > 0:
            print(f"  > Filtering {n_ctrl} control probe transcripts "
                  f"({n_ctrl/len(df_assigned)*100:.1f}%)")
            df_assigned = df_assigned[~ctrl_mask]

        # Notebook uses nuclear labels (in_cell) for the AnnData count matrix when expansion is active
        if nuclei_masks is not None:
            df_nuclear = df_assigned[df_assigned['in_cell'] > 0]
            cell_gene_matrix = pd.crosstab(df_nuclear['in_cell'], df_nuclear[gene_col])
        else:
            df_for_matrix = df_assigned[df_assigned['cell_id_reseg'] > 0]
            cell_gene_matrix = pd.crosstab(df_for_matrix['cell_id_reseg'], df_for_matrix[gene_col])

        adata = sc.AnnData(cell_gene_matrix)
        print(f"  > Cell×Gene matrix: {adata.n_obs} cells, {adata.n_vars} genes")

        # --- 3-7. Map centroids + morphology to AnnData ---
        cell_indices = adata.obs.index.astype(int)
        centroids = np.array([centroid_dict.get(idx, (np.nan, np.nan)) for idx in cell_indices])
        adata.obs['y_centroid'] = centroids[:, 0]
        adata.obs['x_centroid'] = centroids[:, 1]

        # Cell area and perimeter (pixel units + µm conversion)
        adata.obs['cell_area_px'] = [area_dict.get(idx, np.nan) for idx in cell_indices]
        adata.obs['cell_perimeter_px'] = [perimeter_dict.get(idx, np.nan) for idx in cell_indices]

        pixel_to_um = 1.0 / reseg_config.get('um_per_pixel_inv', 4.70588)
        adata.obs['cell_area_um2'] = adata.obs['cell_area_px'] * (pixel_to_um ** 2)
        adata.obs['cell_perimeter_um'] = adata.obs['cell_perimeter_px'] * pixel_to_um

        print(f"  > Median cell area: {adata.obs['cell_area_um2'].median():.1f} µm²")

        # --- 3-7b. Flag edge artifact cells (BEFORE QC so indices match) ---
        # Cells whose expanded masks touch image borders are likely artifacts
        if nuclei_masks is not None and expansion_distance > 0:
            h, w = masks.shape
            edge_labels = set()
            edge_labels.update(np.unique(masks[0, :]))
            edge_labels.update(np.unique(masks[-1, :]))
            edge_labels.update(np.unique(masks[:, 0]))
            edge_labels.update(np.unique(masks[:, -1]))
            edge_labels.discard(0)
            adata.obs['is_edge_cell'] = [idx in edge_labels for idx in cell_indices]
            n_edge = adata.obs['is_edge_cell'].sum()
            print(f"  > Edge artifact cells (touching image border): {n_edge}")

        # --- 3-7c. QC filtering (match Step 0 thresholds) ---
        fmt_config = config.get('formatting', {})
        min_counts = fmt_config.get('mincounts', 10)
        min_genes = fmt_config.get('mingenes', 3)
        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        n_after = adata.n_obs
        print(f"  > QC filter (min_counts={min_counts}, min_genes={min_genes}): "
              f"{n_before} → {n_after} cells ({n_before - n_after} removed)")

        # --- 3-8. Save AnnData + resegmented transcripts ---
        adata_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_resegmented.h5ad")
        transcripts_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_transcripts_resegmented.csv")

        print(f"Saving AnnData to {adata_out_path}...")
        print(f"  > Final: {adata.n_obs} cells, {adata.n_vars} genes, "
              f"median counts/cell={int(np.median(adata.X.sum(axis=1)))}")
        adata.write(adata_out_path)
        print(f"Saving re-assigned transcripts to {transcripts_out_path}...")
        df_assigned.to_csv(transcripts_out_path)

    except Exception as e:
        print(f"Error creating AnnData: {e}")
        return

    # --- 3-9. Domain assignment (optional, JSON polygons) ---

    print("Running Domain Assignment (Polygon Overlap)...")

    domain_map_path = reseg_config.get('domain_assignment', {}).get('domain_map_path', None)

    if domain_map_path:
        if not os.path.exists(domain_map_path):
             pass

    if domain_map_path and os.path.exists(domain_map_path):
        print(f"  > Loading Domain Map from: {domain_map_path}")
        try:
            with open(domain_map_path, 'r') as f:
                domain_data = json.load(f)

            adata.obs['region_annotation'] = 'None'

            count_assigned = 0

            print(f"  > Processing {len(domain_data)} domains...")

            for region in tqdm(domain_data, desc="Assigning Domains"):
                region_name = region.get('name', 'Unknown')

                if 'coordinates' not in region or not region['coordinates']:
                    continue

                # coordinates[0] is the main polygon ring; coords are [y, x] per notebook convention
                raw_coords = region['coordinates'][0]

                if raw_coords[0] != raw_coords[-1]:
                    raw_coords.append(raw_coords[0])

                poly_coords = [(pt[0], pt[1]) for pt in raw_coords]
                poly = Polygon(poly_coords)

                # Bounding-box pre-filter before expensive point-in-polygon checks
                min_x, min_y, max_x, max_y = poly.bounds

                # Poly coords are (y, x) so bounds map to (y_min, x_min, y_max, x_max)
                candidates = adata.obs[
                    (adata.obs['y_centroid'] >= min_x) & (adata.obs['y_centroid'] <= max_x) &
                    (adata.obs['x_centroid'] >= min_y) & (adata.obs['x_centroid'] <= max_y)
                ].index

                for cell_id in candidates:
                    y = adata.obs.loc[cell_id, 'y_centroid']
                    x = adata.obs.loc[cell_id, 'x_centroid']
                    pnt = Point(y, x)

                    if pnt.within(poly):
                        adata.obs.loc[cell_id, 'region_annotation'] = region_name
                        count_assigned += 1

            print(f"  > Domain Assignment Complete. Assigned {count_assigned} cells.")

            # C1-b: Domain polygon + centroid overlay
            try:
                fig, ax = plt.subplots(figsize=(10, 10))
                regions_assigned = adata.obs['region_annotation'].unique()
                cmap_dom = plt.cm.get_cmap('tab20', len(regions_assigned))
                color_map = {r: cmap_dom(i) for i, r in enumerate(sorted(regions_assigned))}
                for r in sorted(regions_assigned):
                    mask_r = adata.obs['region_annotation'] == r
                    ax.scatter(adata.obs.loc[mask_r, 'x_centroid'],
                               adata.obs.loc[mask_r, 'y_centroid'],
                               s=0.5, alpha=0.5, color=color_map[r], label=str(r), rasterized=True)
                # Overlay domain polygon outlines
                for region in domain_data:
                    rname = region.get('name', 'Unknown')
                    if 'coordinates' not in region or not region['coordinates']:
                        continue
                    raw_c = region['coordinates'][0]
                    # Polygon coords are (y, x) — plot as (x, y)
                    xs = [pt[1] for pt in raw_c] + [raw_c[0][1]]
                    ys = [pt[0] for pt in raw_c] + [raw_c[0][0]]
                    ax.plot(xs, ys, linewidth=1, alpha=0.8,
                            color=color_map.get(rname, 'gray'))
                ax.set_aspect('equal')
                ax.invert_yaxis()
                ax.set_title('Domain Polygons + Cell Centroids')
                ax.legend(markerscale=10, fontsize=7, loc='center left', bbox_to_anchor=(1, 0.5))
                fig.tight_layout()
                sample_tag_dom = config.get('sample_tag', 'sample')
                fig.savefig(os.path.join(output_dir, f"{sample_tag_dom}_step3_domain_polygon_overlay.png"),
                            dpi=150, bbox_inches='tight')
                plt.close(fig)
                print("  > Domain polygon + centroid overlay saved.")
            except Exception as e:
                print(f"  > Warning: domain polygon overlay failed: {e}")
                plt.close('all')

            print(f"Saving AnnData with Domains to {adata_out_path}...")
            adata.write(adata_out_path)

        except Exception as e:
            print(f"Error in Domain Assignment: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Domain map not found at {domain_map_path}. Skipping Domain Assignment.")

    print("Step 3: Resegmentation & Domain Assignment Completed.")
