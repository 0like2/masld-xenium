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
import json
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt

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
    
    # 5. Domain Assignment (Polygon Overlap)
    print("Running Domain Assignment (Polygon Overlap)...")
    
    # Get Domain Map Path
    domain_map_path = reseg_config.get('domain_assignment', {}).get('domain_map_path', None)
    
    # Adjust path if relative (Config usually relative to root)
    if domain_map_path:
        # Check absolute or relative to root
        if not os.path.exists(domain_map_path):
             # Try prepending project root (../..) assumption or just rely on user provided path
             pass 
    
    if domain_map_path and os.path.exists(domain_map_path):
        print(f"  > Loading Domain Map from: {domain_map_path}")
        try:
            with open(domain_map_path, 'r') as f:
                domain_data = json.load(f)
                
            # Initialize region annotation
            adata.obs['region_annotation'] = 'None'
            
            # Iterate through regions in JSON
            # Structure from Notebook: List of objects with 'coordinates' and 'name'
            count_assigned = 0
            
            print(f"  > Processing {len(domain_data)} domains...")
            
            for region in tqdm(domain_data, desc="Assigning Domains"):
                region_name = region.get('name', 'Unknown')
                
                # Check for coordinates
                if 'coordinates' not in region or not region['coordinates']:
                    continue
                    
                # Notebook logic implies coordinates[0] is the main polygon ring
                # Coords are [[y, x], ...] based on Notebook dataframe construction
                # Notebook: output.loc[nu,:] = [ob['coordinates'][0][num][0], ob['coordinates'][0][num][1]...] -> y, x
                
                raw_coords = region['coordinates'][0] 
                
                # Ensure closed polygon
                if raw_coords[0] != raw_coords[-1]:
                    raw_coords.append(raw_coords[0])
                    
                # Create Shapely Polygon
                # Note: Shapely uses (x, y) usually, but here we construct Point(y, x) later to match.
                # Consistent coordinate system is key.
                poly_coords = [(pt[0], pt[1]) for pt in raw_coords]
                poly = Polygon(poly_coords)
                
                # Assign cells
                # Vectorized point checking is hard with pure Shapely, check bounding box first?
                # Optimization: Check bounding box of polygon vs cells
                min_x, min_y, max_x, max_y = poly.bounds
                
                # Select candidate cells (Approximate filter)
                # Note: We store y_centroid, x_centroid
                # Poly coords are (y, x) based on notebook logic
                
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
            
            # Save again with domains
            print(f"Saving AnnData with Domains to {adata_out_path}...")
            adata.write(adata_out_path)
            
        except Exception as e:
            print(f"Error in Domain Assignment: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Domain map not found at {domain_map_path}. Skipping Domain Assignment.")
        
    print("Step 3: Resegmentation & Domain Assignment Completed.")

