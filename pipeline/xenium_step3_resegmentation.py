import os
import numpy as np
import pandas as pd
import tifffile as tf
import scanpy as sc
from cellpose import models
from skimage.measure import label, regionprops
from skimage.segmentation import expand_labels
from tqdm import tqdm
import gc

def run_step3(config, dapi_image_path, transcripts_path, output_dir):
    """
    Step 3: Resegmentation (Cellpose)
    Resegments nuclei using Cellpose and re-assigns transcripts.
    Ref: Notebook 3.1 & 3.2
    """
    print("------------------------------------------------")
    print("Starting Step 3: Resegmentation (Cellpose)")
    print("------------------------------------------------")

    reseg_config = config.get('resegmentation', {})
    if not reseg_config.get('run_resegmentation', False):
        print("Resegmentation is disabled in config. Skipping.")
        return

    # 1. Load Data
    print(f"Loading DAPI image from {dapi_image_path}...")
    if not os.path.exists(dapi_image_path):
        print(f"Error: DAPI image not found at {dapi_image_path}")
        return
        
    try:
        with tf.TiffFile(dapi_image_path) as tif:
            dapi_image = tf.imread(dapi_image_path)
    except Exception as e:
        print(f"Error reading DAPI image: {e}")
        return

    print(f"DAPI Image Shape: {dapi_image.shape}")
    
    # 2. Cellpose Segmentation
    cp_params = reseg_config.get('cellpose', {})
    use_gpu = cp_params.get('use_gpu', False)
    diameter = cp_params.get('diameter', None) 
    
    print(f"Initializing Cellpose (GPU={use_gpu}, model='nuclei')...")
    model = models.CellposeModel(gpu=use_gpu, model_type='nuclei')
    
    tile_size = 2000
    if dapi_image.shape[0] > tile_size or dapi_image.shape[1] > tile_size:
        print("Image larger than tile size. Running standard Cellpose eval (handles tiling internaly)...")
        masks, flows, styles, diams = model.eval(dapi_image, diameter=diameter, channels=[0,0], tile=True)
    else:
        masks, flows, styles, diams = model.eval(dapi_image, diameter=diameter, channels=[0,0])
        
    del dapi_image, flows, styles
    gc.collect()

    # Save Masks as TIF (Required for Baysor in Step 6)
    mask_out_path = os.path.join(output_dir, f"{config.get('sample_tag', 'sample')}_step3_resegmented_masks.tif")
    print(f"Saving Cellpose masks to TIF: {mask_out_path}")
    tf.imwrite(mask_out_path, masks)
    
    # 3. Transcript Assignment
    print("Loading Transcripts...")
    try:
        df_transcripts = pd.read_csv(transcripts_path)
    except Exception as e:
        print(f"Error loading transcripts: {e}")
        return
        
    um_per_pixel = reseg_config.get('um_per_pixel', 1.0)
    
    if 'x_global_px' in df_transcripts.columns:
        x_col, y_col = 'x_global_px', 'y_global_px'
    elif 'global_x' in df_transcripts.columns:
        x_col, y_col = 'global_x', 'global_y'
    else:
        x_col, y_col = 'x_location', 'y_location'

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
    df_assigned = df_valid[df_valid['cell_id_reseg'] > 0].copy()
    
    print(f"  > Assigned {len(df_assigned)} transcripts to {df_assigned['cell_id_reseg'].nunique()} cells.")
    
    # 4. Generate AnnData
    print("Generating Cell x Gene Matrix...")
    try:
        gene_col = 'feature_name' if 'feature_name' in df_assigned.columns else 'gene'
        cell_gene_matrix = pd.crosstab(df_assigned['cell_id_reseg'], df_assigned[gene_col])
        adata = sc.AnnData(cell_gene_matrix)
        
        print("Calculating centroids...")
        props = regionprops(masks)
        centroid_dict = {p.label: p.centroid for p in props}
        
        centroids = [centroid_dict.get(idx, (np.nan, np.nan)) for idx in adata.obs.index.astype(int)]
        centroids = np.array(centroids)
        adata.obs['y_centroid'] = centroids[:, 0]
        adata.obs['x_centroid'] = centroids[:, 1]
        
        # Save Outputs for Step 4
        adata_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_resegmented.h5ad")
        transcripts_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_transcripts_resegmented.csv")
        
        print(f"Saving AnnData to {adata_out_path}...")
        adata.write(adata_out_path)
        print(f"Saving re-assigned transcripts to {transcripts_out_path}...")
        df_assigned.to_csv(transcripts_out_path)
        
    except Exception as e:
        print(f"Error creating AnnData: {e}")
        return
    
    # 5. Domain Assignment
    print("Running Domain Assignment (Spatial Clustering)...")
    if adata.n_obs > 50:
        try:
            domain_config = config.get('resegmentation', {}).get('domain_assignment', {})
            resolution = domain_config.get('resolution', 0.5)
            
            print(f"  > Computing spatial neighbors...")
            adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
            sc.pp.neighbors(adata, use_rep='spatial', n_neighbors=30, key_added='spatial')
            
            print(f"  > Running Leiden clustering for domains (res={resolution})...")
            sc.tl.leiden(adata, resolution=resolution, key_added='spatial_domain', neighbors_key='spatial')
            print(f"  > Domains assigned. Found {adata.obs['spatial_domain'].nunique()} domains.")
            
            # Save again with domains
            adata.write(adata_out_path)
            
        except Exception as e:
            print(f"Error in Domain Assignment: {e}")
    else:
        print("Not enough cells for Domain Assignment.")
        
    print("Step 3: Resegmentation & Domain Assignment Completed.")

