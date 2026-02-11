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


def _detect_device():
    """Auto-detect best available device: CUDA > MPS > CPU."""
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        print(f"  Device: CUDA ({name})")
        return torch.device("cuda"), True
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        print("  Device: MPS (Apple Silicon)")
        return torch.device("mps"), False
    print("  Device: CPU")
    return torch.device("cpu"), False

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

    # --- 3-1. Device detection (CUDA/MPS/CPU) ---
    device, use_gpu = _detect_device()
    print(f"Initializing Cellpose (device={device}, model='nuclei')...")
    model = models.CellposeModel(gpu=use_gpu, model_type='nuclei', device=device)

    tile_size = 2000
    if dapi_image.shape[0] > tile_size or dapi_image.shape[1] > tile_size:
        print("Image larger than tile size. Running standard Cellpose eval (handles tiling internally)...")
        masks, flows, styles, diams = model.eval(dapi_image, diameter=diameter, channels=[0,0], tile=True)
    else:
        masks, flows, styles, diams = model.eval(dapi_image, diameter=diameter, channels=[0,0])

    del dapi_image, flows, styles
    gc.collect()

    # --- 3-3b. Label expansion (nuclei → cytoplasm approximation) ---
    # Ensure contiguous labels via connected-component analysis (matches notebook label() call)
    masks = label(masks)
    print(f"  > Labeled {masks.max()} nuclei objects")

    nuclei_masks = masks  # Keep original nuclei masks for dual assignment
    expansion_distance = reseg_config.get('expansion_distance', 400)
    if expansion_distance > 0:
        print(f"Expanding nuclei masks by {expansion_distance} pixels...")
        masks = expand_labels(nuclei_masks, distance=expansion_distance)
        n_before = len(np.unique(nuclei_masks)) - 1
        n_after = len(np.unique(masks)) - 1
        print(f"  > Nuclei labels: {n_before}, Expanded labels: {n_after}")
    else:
        nuclei_masks = None  # No dual assignment needed

    # --- 3-4. Save masks → TIF ---
    mask_out_path = os.path.join(output_dir, f"{config.get('sample_tag', 'sample')}_step3_resegmented_masks.tif")
    print(f"Saving Cellpose masks to TIF: {mask_out_path}")
    tf.imwrite(mask_out_path, masks)

    if nuclei_masks is not None:
        nuclei_mask_path = os.path.join(output_dir, f"{config.get('sample_tag', 'sample')}_step3_nuclei_masks.tif")
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
            axes[1].imshow(label2rgb(nuclei_masks[::step, ::step], bg_label=0), interpolation='nearest')
            axes[1].set_title(f'Nuclei masks ({nuclei_masks.max()} labels)')
            axes[1].axis('off')

        plt.tight_layout()
        plt.savefig(qc_path, dpi=150)
        plt.close()
        print("  > Mask QC plot saved.")
    except Exception as e:
        print(f"  > Warning: Could not generate mask QC plot: {e}")

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
    try:
        y_coords = df_transcripts[y_col].values.astype(int)
        x_coords = df_transcripts[x_col].values.astype(int)
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

        # Notebook uses nuclear labels (in_cell) for the AnnData count matrix when expansion is active
        if nuclei_masks is not None:
            df_nuclear = df_assigned[df_assigned['in_cell'] > 0]
            cell_gene_matrix = pd.crosstab(df_nuclear['in_cell'], df_nuclear[gene_col])
        else:
            cell_gene_matrix = pd.crosstab(df_assigned['cell_id_reseg'], df_assigned[gene_col])

        adata = sc.AnnData(cell_gene_matrix)

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

        # --- 3-8. Save AnnData + resegmented transcripts ---
        adata_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_resegmented.h5ad")
        transcripts_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_transcripts_resegmented.csv")

        print(f"Saving AnnData to {adata_out_path}...")
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

            print(f"Saving AnnData with Domains to {adata_out_path}...")
            adata.write(adata_out_path)

        except Exception as e:
            print(f"Error in Domain Assignment: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Domain map not found at {domain_map_path}. Skipping Domain Assignment.")

    print("Step 3: Resegmentation & Domain Assignment Completed.")
