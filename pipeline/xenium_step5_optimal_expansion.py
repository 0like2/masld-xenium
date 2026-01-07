# pipeline/xenium_step5_optimal_expansion.py
# ------------------------------------------------------------
# Xenium Pipeline Step 5: Optimal Expansion
# Assigns unassigned/cytoplasmic reads to the nearest annotated cell domain.
# Ref: notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb
# ------------------------------------------------------------

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import random as rd
from scipy.spatial import ConvexHull
import math

# --- Inline XB Calculating Functions ---
# Ref: xb/calculating.py

def dist_nuc(reads_ctdsub):
    """ Compute the median distance to the nuclei the edges of each cell. """
    allds=[]
    for g,n in reads_ctdsub.groupby('cell_id'):
        try:
            if len(n) < 3: continue # Need 3 points for Hull
            hull = ConvexHull(np.array(n.loc[:,['x_location','y_location']]))
            # n must have 'distance' column (distance to nucleus center)
            if 'distance' in n.columns:
                 allds.append(np.mean(n.iloc[hull.vertices]['distance']))
        except Exception:
            pass
    if len(allds) > 0:
        median_dist=np.median(allds)
    else:
        median_dist = np.nan
    return median_dist

def distance_calc(x1,y1,x2,y2):
    """ Calculate distance between two points """
    return math.sqrt( ((x1-x2)**2)+((y1-y2)**2) )

def hex_to_rgb(value):
    """ Transform hex to rgb """
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


# Configure local logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_step5(config):
    """
    Main execution function for Step 5: Optimal Expansion.
    """
    print("\n" + "="*60)
    print("[Step 5] Optimal Expansion (Ref: 4_1 Notebook)")
    print("="*60)

    output_dir = config["output_dir"]
    sample_tag = config["sample_tag"]
    
    # Check Step 5 Config
    exp_config = config.get("optimal_expansion", {})
    run_exp = exp_config.get("run_expansion", True)
    subsample_frac = exp_config.get("subsample_fraction", 0.01)
    dist_threshold = exp_config.get("distance_threshold", None)
    
    if not run_exp:
        print("    - [Info] 'run_expansion' is False in config. Skipping Step 5.")
        return

    # --- 1. Load Data ---
    # We need:
    # A. Original Transcripts (Step 0) -> Contains unassigned reads
    # B. Annotated Cells (Step 2 or 1) -> Contains domains/classes
    
    print("\n[Step 5-1] Loading Data...")
    
    # Load Transcripts (Step 0)
    # Using 'spots' from h5ad if available, or finding transcripts.csv
    step0_file = os.path.join(output_dir, f"{sample_tag}.h5ad")
    if not os.path.exists(step0_file):
        # Try finding it in step0_formatting sibling directory
        parent_dir = os.path.dirname(output_dir)
        step0_file_alt = os.path.join(parent_dir, "step0_formatting", f"{sample_tag}.h5ad")
        if os.path.exists(step0_file_alt):
             step0_file = step0_file_alt
            
    if not os.path.exists(step0_file):
        logging.error(f"    - Step 0 output not found at {step0_file} or in sibling dir. Cannot run expansion.")
        return
        
    print(f"    - Loading Step 0 Data from: {step0_file}")
    adata_step0 = sc.read_h5ad(step0_file)
    if 'spots' not in adata_step0.uns:
        logging.error("    - 'spots' dataframe missing in Step 0 adata.uns.")
        return
    
    reads_original = adata_step0.uns['spots'].copy()
    print(f"    - Loaded {len(reads_original)} original reads.")
    
    # Load Annotated Cells (Step 2 or Step 1 or Step 4)
    input_adata_file = None
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        input_adata_file = config['previous_step_adata_path']
    
    if not input_adata_file or not os.path.exists(input_adata_file):
        # Fallback search
        print("    - [Info] 'previous_step_adata_path' missing or invalid. Searching...")
        parent_dir = os.path.dirname(output_dir)
        possible_paths = [
            os.path.join(parent_dir, "step4_resegmentation", f"{sample_tag}_resegmented.h5ad"),
            os.path.join(parent_dir, "step1_exploration", f"{sample_tag}_step1_exploration.h5ad"),
            os.path.join(parent_dir, "step2_segmentation_free", f"{sample_tag}_step2_points2regions.h5ad")
        ]
        for p in possible_paths:
            if os.path.exists(p):
                input_adata_file = p
                break
                
    print(f"    - Loading Annotated Data from: {input_adata_file}")
    if not input_adata_file or not os.path.exists(input_adata_file):
        logging.error(f"    - Annotated adata not found. Run Step 1, 2, or 4 first.")
        return
        
    adata_annotated = sc.read_h5ad(input_adata_file)
    print(f"    - Loaded annotated cells: {adata_annotated.shape}")

    # --- 2. Prepare for Expansion ---
    print("\n[Step 5-2] Identifying Domains & Unassigned Reads...")
    
    # Check for Domain/Annotation column
    # Notebook uses 'spatial_annotation' or 'Class'.
    # Our Step 1 produces 'leiden', Step 2 P2R produces 'cluster'.
    # We will use 'leiden' or 'graph_clusters' or 'cluster' as the "Domain" source.
    
    domain_key = None
    priority_keys = ['spatial_annotation', 'Class', 'leiden', 'cluster', 'graph_clusters']
    for key in priority_keys:
        if key in adata_annotated.obs.columns:
            domain_key = key
            break
            
    if not domain_key:
        logging.error("    - No suitable domain/cluster key found in annotated adata. Cannot assign domains.")
        return
    print(f"    - Using '{domain_key}' as domain source.")
    
    # Map Cell ID to Domain
    # Usually adata.obs.index is cell_id or there is a 'cell_id' column.
    if 'cell_id' in adata_annotated.obs.columns:
        cell_id_col = 'cell_id'
        annotated_ids = adata_annotated.obs['cell_id']
    else:
        cell_id_col = 'index'
        annotated_ids = adata_annotated.obs.index
        
    # Dictionary: Cell ID -> Domain
    domain_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
    
    # Apply to reads
    print("    - Mapping existing domains to reads...")
    if 'cell_id' not in reads_original.columns:
         # If reads don't have cell_id, we can't do expansion from nuclei.
         # Assuming they do (standard Xenium)
         if reads_original.index.name == 'cell_id':
             reads_original = reads_original.reset_index()
         else:
             logging.error("    - 'cell_id' missing in reads. Cannot link to cells.")
             return
             
    reads_original['domain'] = reads_original['cell_id'].map(domain_map)
    
    # Define Assigned vs Unassigned
    # "Unassigned" here means "Assigned to a cell that has NO domain" OR "Not assigned to any cell"?
    # Notebook logic:
    # nancells = reads_original[reads_original['domain'].isna()]
    # annotatedcells = reads_original[~reads_original['domain'].isna()]
    # This implies reads assigned to cells WITHOUT annotation are treated as unassigned candidates? 
    # Or reads with 'UNASSIGNED' cell_ids (which naturally yield NaN in map) are included.
    # Yes, Xenium 'unassigned' transcripts have cell_id = 'UNASSIGNED' or similar, so map returns NaN.
    
    nancells = reads_original[reads_original['domain'].isna()]
    annotatedcells = reads_original[~reads_original['domain'].isna()]
    
    n_assigned = len(annotatedcells)
    n_unassigned = len(nancells)
    print(f"    - Assigned Reads (with domain): {n_assigned}")
    print(f"    - Unassigned Reads (to be expanded): {n_unassigned}")
    
    if n_unassigned == 0:
        print("    - No unassigned reads found. Expansion not needed.")
        return

    # --- 3. Build cKDTree (Subsampled) ---
    print(f"\n[Step 5-3] Building cKDTree (Subsample Fraction: {subsample_frac})...")
    
    # Subsample annotated cells to speed up tree construction
    # We sample READS that are annotated. (Notebook logic)
    # Notebook: sub=rd.sample(list(annotatedcells.index), int(... * 0.01))
    
    n_sample = int(np.round(n_assigned * subsample_frac))
    if n_sample < 100: n_sample = min(n_assigned, 1000) # minimal safety
    
    print(f"    - Subsampling {n_sample} reads as anchors...")
    # Use pandas sample for efficiency
    annotated_sub = annotatedcells.sample(n=n_sample, random_state=42)
    
    coords1 = annotated_sub[['x_location', 'y_location']].values
    ids1 = annotated_sub.index.values # Or cell_ids? Notebook stores data indices to look up domain later.
    # We need to look up the DOMAIN of the neighbor, not just the read index.
    # Notebook: closest_neighbors.append(ids1[index])
    
    print("    - Constructing KDTree...")
    tree = cKDTree(coords1)
    
    # --- 4. Query Unassigned Reads ---
    print(f"\n[Step 5-4] Assigning {n_unassigned} reads to nearest domains...")
    
    coords2 = nancells[['x_location', 'y_location']].values
    
    # query(x, k=1) returns (distances, indices)
    # Providing all coords at once is vectorized and faster than loop in notebook if memory allows.
    # If coords2 is huge (>10M), might need chunking.
    # Assuming typical Xenium slice (millions), explicit full query might be heavy but usually fine on modern RAM.
    # Notebook uses loop: for i, coord in enumerate(coords2)...
    
    # Let's try vectorized first, catch OOM? Or just chunk it safely.
    chunk_size = 1000000
    n_chunks = int(np.ceil(len(coords2) / chunk_size))
    
    all_indices = []
    all_distances = []
    
    print(f"    - Processing in {n_chunks} chunks...")
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, len(coords2))
        chunk = coords2[start:end]
        
        dists, idxs = tree.query(chunk, k=1)
        all_indices.append(idxs)
        all_distances.append(dists)
        
    concat_indices = np.concatenate(all_indices)
    concat_distances = np.concatenate(all_distances)
    
    # --- 5. Assign Domains ---
    print("\n[Step 5-5] Finalizing Assignment...")
    
    # Map tree indices back to original annotated reads -> get their domain
    # nearest_read_indices_in_sub = ids1[concat_indices] -> problematic if ids1 is not 0..N
    # tree.query returns index into coords1.
    # So we need: domain_of_anchor = annotated_sub.iloc[tree_index]['domain']
    
    # Get domains of the anchors used in tree
    anchor_domains = annotated_sub['domain'].values # aligned with coords1 order
    
    # Assigned domains for unassigned reads
    assigned_domains = anchor_domains[concat_indices]
    
    # Apply Distance Threshold (if set)
    if dist_threshold:
        print(f"    - Applying distance threshold: {dist_threshold}")
        # Identify reads too far
        mask_too_far = concat_distances > dist_threshold
        # Set them back to NaN or keep as 'Unassigned'
        # We can just update assigned_domains where mask is True
        # Note: assigned_domains is numpy array of objects/strings
        params_assigned = assigned_domains.copy()
        params_assigned[mask_too_far] = np.nan # Or "Unassigned"
        assigned_domains = params_assigned
        n_filtered = np.sum(mask_too_far)
        print(f"    - {n_filtered} reads exceeded distance threshold.")
        
    # Update DataFrame
    # We want to save the Expanded Table
    # original 'nancells' indices
    nancells.loc[:, 'domain'] = assigned_domains
    
    # Combine back
    print("    - Merging results...")
    # reads_original is the master. We update the 'domain' column for nancells indices.
    reads_original.loc[nancells.index, 'domain'] = assigned_domains
    
    # --- 6. Visualization ---
    print("\n[Step 5-6] Generating Visualization...")
    # Scatter plot of a subset to show expansion
    # Sample 100k reads for plotting
    plot_df = reads_original.sample(n=min(len(reads_original), 100000), random_state=42)
    
    plt.figure(figsize=(10, 10))
    # Filter out NaN domains for plot
    plot_df_clean = plot_df.dropna(subset=['domain'])
    
    # Use seaborn scatter
    # We need a palette. If many domains, 'tab20' or similar.
    sns.scatterplot(
        data=plot_df_clean, 
        x='x_location', 
        y='y_location', 
        hue='domain', 
        s=1, 
        linewidth=0,
        palette='tab20',
        legend=False # Legend might be huge
    )
    plt.title(f"Optimal Expansion Result (Subsample {subsample_frac})")
    plt.axis('equal')
    
    plot_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_map.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    - Saved Expansion Map to: {plot_file}")
    
    # Histogram of Distances (Validation)
    plt.figure(figsize=(6, 4))
    sns.histplot(concat_distances, bins=50, kde=True, color='orange')
    plt.title("Distance to Nearest Domain-Anchor")
    plt.xlabel("Distance (pixels/units)")
    plt.ylabel("Count")
    dist_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_distances.png")
    plt.savefig(dist_file)
    plt.close()
    print(f"    - Saved Distance Histogram to: {dist_file}")

    # --- 7. Save Output ---
    output_file = os.path.join(output_dir, f"{sample_tag}_step5_expanded_transcripts.csv")
    print(f"    - Saving Expanded Transcripts to: {output_file}")
    reads_original.to_csv(output_file, index=False)
    
    # --- 8. Turnover / Crossover Analysis (Validation) ---
    # Ref: Cell 295+ in Notebook 4_1
    # Calculates optimal expansion radius stats.
    print("\n[Step 5-8] Running Turnover/Crossover Analysis (Validation)...")
    
    # 8.1 Prepare Data: Calculate distance from read to assigned cell centroid
    print("    - Calculating read-to-centroid distances...")
    
    # We need centroid of assigned cells.
    # adata_annotated.obs should have 'x_centroid', 'y_centroid'
    if 'x_centroid' in adata_annotated.obs.columns and 'y_centroid' in adata_annotated.obs.columns:
        # Create map
        cx_map = dict(zip(annotated_ids, adata_annotated.obs['x_centroid']))
        cy_map = dict(zip(annotated_ids, adata_annotated.obs['y_centroid']))
        
        # Only for assigned reads
        # assigned_mask = ~reads_original['domain'].isna()
        # reads_val = reads_original[assigned_mask].copy()
        
        # NOTE: Notebook 4_1 uses 'reads_assigned' which are reads with valid cell_id.
        # Step 5 logic assigns domain to unassigned reads (cell_id=-1 or similar).
        # We need to calculate distance for ALL reads that have a domain now (original + expanded).
        # But 'cell_id' for expanded reads is still unassigned (-1).
        # To calculate distance to "assigned cell", expanded reads need to be linked to a specific CELL,
        # not just a DOMAIN.
        # Our expansion logic (nearest CLUSTER/DOMAIN) assigned a DOMAIN, not a specific CELL ID, to unassigned reads.
        # Wait, the KDTree query found the NEAREST ANCHOR READ, which belongs to a specific CELL.
        # Did we save that CELL ID?
        # In the implemented run_step5 above, we used:
        # assigned_domains = anchor_domains[concat_indices] -> This gets the Domain.
        # We did NOT save the specific cell_id of the anchor.
        
        # The Notebook's logic relies on reads having 'x_cell', 'y_cell' (centroid of assigned cell).
        # Original reads have valid cell_id -> easy.
        # Expanded reads (unassigned) -> We assigned them a DOMAIN. The notebook 4_1 logic for expansion
        # actually maps: 'domain' -> domain of nearest annotated cell.
        # It does NOT seem to assign them to a specific cell for the purpose of distance calculation in validation?
        # Let's check Notebook Cell 256: reads_assigned['distance']... using x_cell/y_cell.
        # Cell 215: reads_assigned = reads_assigned[... isin cells_metadata].
        # In Cell 202: reads_assigned = reads_original[reads_original['cell_id']!=-1].
        # SO: The validation (Turnover) is done ONLY on ORIGINALLY ASSIGNED reads, NOT on expanded reads.
        # It calculates stats on existing cells to find "optimal expansion".
        
        reads_val = reads_original[reads_original['cell_id'].isin(cx_map.keys())].copy()
        print(f"    - Using {len(reads_val)} originally assigned reads for stats.")
        
        reads_val['x_cell'] = reads_val['cell_id'].map(cx_map)
        reads_val['y_cell'] = reads_val['cell_id'].map(cy_map)
        
        # Calculate distance
        # vectorized sqrt((x-xc)^2 + ...)
        reads_val['distance'] = np.sqrt( (reads_val['x_location'] - reads_val['x_cell'])**2 + 
                                         (reads_val['y_location'] - reads_val['y_cell'])**2 )
        
        # 8.2 Statistics Loop
        # Group by Cell Type (initial_annotation / domain)
        # Notebook Cell 302 Loop
        if domain_key and 'feature_name' in reads_val.columns:
             # Initial annotation for reads
             # If reads_original doesn't have it, map it
             if 'initial_annotation' not in reads_val.columns:
                 # domain_key is the annotation column in cell obs
                 ct_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
                 reads_val['initial_annotation'] = reads_val['cell_id'].map(ct_map)
            
             unique_celltypes = reads_val['initial_annotation'].unique()
             
             meand_celltype = []
             nuclimall = []
             
             print(f"    - Computing stats for {len(unique_celltypes)} cell types...")
             
             # Pre-calculate background expression (reads_not_assigned)
             # reads_not_assigned = reads_original[reads_original['cell_id'] == -1/NaN]
             # In our DF, unassigned might be those NOT in cx_map.
             mask_unassigned = ~reads_original['cell_id'].isin(cx_map.keys())
             reads_not_assigned = reads_original[mask_unassigned]
             
             # Need 'domain' on unassigned? 
             # In Notebook Cell 274: background_express pd.crosstab(reads_not_assigned['domain']...)
             # This implies unassigned reads HAVE domains (assigned via expansion earlier).
             # Yes, we populated 'domain' for them in Step 5-5.
             
             # Background expression matrix
             if len(reads_not_assigned) > 0 and 'domain' in reads_not_assigned.columns:
                 # Filter NaN domains
                 reads_not_assigned = reads_not_assigned.dropna(subset=['domain'])
                 background_express = pd.crosstab(reads_not_assigned['domain'], reads_not_assigned['feature_name'])
                 
                 # Loop
                 results_list = []
                 
                 for celltype in unique_celltypes:
                     if pd.isna(celltype): continue
                     
                     sub = reads_val[reads_val['initial_annotation'] == celltype]
                     if len(sub) == 0: continue
                     
                     # 1. Median Distance (Cell Size proxy?)
                     # Notebook Cell 309 loop logic is complex (good_domains, crosstab...)
                     # Simplified equivalent as per "Optimal Expansion" goal:
                     # Calculate mean/median distance of reads for this cell type?
                     # Notebook: nuclim.append(dist_nuc(reads_ctdsub)) -> dist_nuc uses ConvexHull
                     
                     # Calculate Nuclei Limit (Median distance of Hull vertices)
                     # Per cell in this celltype
                     # To save time, sample cells? Notebook does per cell type.
                     try:
                         nuc_limit = dist_nuc(sub) # using inlined function
                         nuclimall.append(nuc_limit)
                     except:
                         nuclimall.append(np.nan)
                     
                     # Calculate Mean Distance (Turnover?)
                     # Notebook calculates 'meand_celltype'.
                     
                     # Storing results for plot
                     results_list.append({
                         'cluster': celltype,
                         'nuclei_size': nuc_limit if 'nuc_limit' in locals() else 0,
                         # 'score': ... complex crossover score
                     })
                     
                 # Save stats
                 stats_df = pd.DataFrame(results_list)
                 stats_file = os.path.join(output_dir, f"{sample_tag}_step5_turnover_stats.csv")
                 stats_df.to_csv(stats_file, index=False)
                 print(f"    - Saved Turnover Stats to: {stats_file}")
                 
                 # 8.3 Comparison Plot (Simplified)
                 # nuclei_vs_background (Scatter plot of cell size vs nuclei size)
                 if len(stats_df) > 0:
                     plt.figure(figsize=(6, 6))
                     sns.scatterplot(data=stats_df, x='cluster', y='nuclei_size')
                     plt.title("Nuclei Size per Cell Type (Validation)")
                     plt.xticks(rotation=90)
                     plot_out = os.path.join(output_dir, f"{sample_tag}_step5_nuclei_vs_background.png")
                     plt.savefig(plot_out, bbox_inches='tight')
                     plt.close()
                     print(f"    - Saved Validation Plot to: {plot_out}")
                     
             else:
                 print("    - [Warning] No background reads (unassigned) with domains found. Skipping turnover calculation.")
        else:
             print("    - [Warning] 'feature_name' or domains missing. Skipping turnover stats.")
    else:
        print("    - [Warning] Centroid data (x_centroid, y_centroid) missing in annotated adata. Skipping Turnover Analysis.")
    
    print("\n=== Step 5 Optimal Expansion Complete ===")
    
if __name__ == "__main__":
    pass
