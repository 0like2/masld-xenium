# pipeline/xenium_step0_formatting.py
# ------------------------------------------------------------
# Xenium Pipeline Step 0: Formatting
# Converts raw Xenium output to AnnData and prepares images.
# STANDALONE VERSION (No xb dependency)
# ------------------------------------------------------------

import os
import time
import shutil
import gzip
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmread
import tifffile as tf

# Configure local logging if run mainly
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ================== Core Logic (Inlined from xb.formatting) ==================

def format_xenium_adata_mid_2023(path, tag, output_path):
    """
    Format xenium data (output from the machine) to adata format, 
    considering the format used by Xenium at Q2 2023 (Mid 2023).
    Refactored from xb.formatting.
    """
    print(f"  > Processing Xenium Data at: {path}")

    # Decompress cell_feature_matrix.tar.gz if needed
    if not os.path.exists(os.path.join(path, 'cell_feature_matrix')):
        tar_gz = os.path.join(path, 'cell_feature_matrix.tar.gz')
        if os.path.exists(tar_gz):
            print(f"  > Decompressing {tar_gz}...")
            shutil.unpack_archive(tar_gz, path)
            print("  > Decompression done.")
        else:
            print(f"  [WARNING] 'cell_feature_matrix' folder and tar.gz not found at {path}")

    # Decompress individual files if needed
    for filename in ['transcripts.csv', 'cells.csv']:
        file_path = os.path.join(path, filename)
        gz_path = file_path + '.gz'
        if not os.path.exists(file_path) and os.path.exists(gz_path):
            print(f"  > Decompressing {gz_path}...")
            with gzip.open(gz_path, 'rb') as f_in:
                with open(file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
    
    # Decompress matrix files inside cell_feature_matrix
    cfm_path = os.path.join(path, 'cell_feature_matrix')
    if os.path.exists(cfm_path):
        for filename in ['barcodes.tsv', 'features.tsv', 'matrix.mtx']:
            file_path = os.path.join(cfm_path, filename)
            gz_path = file_path + '.gz'
            if not os.path.exists(file_path) and os.path.exists(gz_path):
                print(f"  > Decompressing {gz_path}...")
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(file_path, 'wb') as f_out:
                         shutil.copyfileobj(f_in, f_out)

    # Read Matrix and Metadata
    print("  > Reading matrix and metadata...")
    try:
        matrix_path = os.path.join(cfm_path, 'matrix.mtx')
        a = mmread(matrix_path)
        ad = a.todense()
        
        cells_path = os.path.join(path, "cells.csv")
        cell_info = pd.read_csv(cells_path)
        
        features_path = os.path.join(cfm_path, 'features.tsv')
        features = pd.read_csv(features_path, sep='\t', header=None, index_col=0)
        
        # Construct AnnData
        adata = sc.AnnData(ad.transpose(), obs=cell_info, var=features)
        adata.var.index.name = 'index'
        adata.var.index.name = 'index'
        adata.var.columns = ['gene_id', 'reason_of_inclusion']
        
        # Add spatial coordinates
        if 'x_centroid' in cell_info.columns and 'y_centroid' in cell_info.columns:
             adata.obsm['spatial'] = np.array(cell_info[['x_centroid', 'y_centroid']])
             print("  > Spatial coordinates (centroids) added to obsm['spatial']")
    except Exception as e:
        raise ValueError(f"Error reading matrix files: {e}")

    # Process Gene Panel JSON
    print("  > Processing gene panel...")
    json_path = os.path.join(path, 'gene_panel.json')
    if os.path.exists(json_path):
        with open(json_path) as f:
            data = json.load(f)
        
        geness = []
        idss = []
        descriptorss = []
        
        targets = data.get('payload', {}).get('targets', [])
        for t in targets:
            g_name = t['type']['data']['name']
            geness.append(g_name)
            
            try:
                idss.append(t['type']['data']['id'])
            except:
                idss.append('newid_' + g_name)
            
            try:
                descriptorss.append(t['type']['descriptor'])
            except:
                descriptorss.append('other')

        dict_inpanel = dict(zip(geness, descriptorss))
        dict_ENSEMBL = dict(zip(geness, idss))
        
        adata.var['Ensembl ID'] = adata.var['gene_id'].map(dict_ENSEMBL)
        adata.var['in_panel'] = adata.var['gene_id'].map(dict_inpanel)

    # Load Transcripts (Spots)
    print("  > Loading transcripts (this may take a while)...")
    transcripts_path = os.path.join(path, 'transcripts.csv')
    if os.path.exists(transcripts_path):
        transcripts = pd.read_csv(transcripts_path, index_col=0)
        adata.uns['spots'] = transcripts
    
    # Load Analysis Results (UMAP, PCA, Clusters)
    print("  > Loading analysis results...")
    if not os.path.exists(os.path.join(path, 'analysis')):
        an_tar = os.path.join(path, 'analysis.tar.gz')
        if os.path.exists(an_tar):
             print(f"  > Decompressing {an_tar}...")
             shutil.unpack_archive(an_tar, path)
             print("  > Analysis decompression done.")

    analysis_base = os.path.join(path, 'analysis')
    if os.path.exists(analysis_base):
        try:
            # UMAP
            umap_path = os.path.join(analysis_base, 'umap/gene_expression_2_components/projection.csv')
            if os.path.exists(umap_path):
                UMAP = pd.read_csv(umap_path, index_col=0)
                adata.obsm['X_umap'] = np.array(UMAP)
            
            # TSNE
            tsne_path = os.path.join(analysis_base, 'tsne/gene_expression_2_components/projection.csv')
            if os.path.exists(tsne_path):
                TSNE = pd.read_csv(tsne_path, index_col=0)
                adata.obsm['X_tsne'] = np.array(TSNE)

            # PCA
            pca_path = os.path.join(analysis_base, 'PCA/gene_expression_10_components/projection.csv')
            if os.path.exists(pca_path):
                PCA = pd.read_csv(pca_path, index_col=0)
                adata.obsm['X_pca'] = np.array(PCA)

            # Clusters
            clusters_path = os.path.join(analysis_base, 'clustering/gene_expression_graphclust/clusters.csv')
            if os.path.exists(clusters_path):
                clusters = pd.read_csv(clusters_path, index_col=0)
                adata.obs['graph_clusters'] = list(clusters['Cluster'].astype(str))

            # KMeans (2 to 10)
            for k in range(2, 11):
                km_path = os.path.join(analysis_base, f'clustering/gene_expression_kmeans_{k}_clusters/clusters.csv')
                if os.path.exists(km_path):
                    km = pd.read_csv(km_path, index_col=0)
                    adata.obs[f'kmeans{k}_clusters'] = list(km['Cluster'].astype(str))
                    
        except Exception as e:
            print(f"  [WARNING] Analysis files load issue: {e}")

    # Finalize
    adata.X = sp.csr_matrix(adata.X)
    if 'spots' in adata.uns:
        adata.uns['spots'] = adata.uns['spots'].fillna(0)
        if 'fov_name' in adata.uns['spots'].columns:
            adata.uns['spots']['fov_name'] = adata.uns['spots']['fov_name'].astype(str)
            
    # Save
    out_file = os.path.join(output_path, f"{tag}.h5ad")
    print(f"  > Saving to {out_file}...")
    adata.write(out_file)
    return adata


def format_background(path):
    """
    Format OME-TIFF background mipped image to .tiff image.
    Refactored from xb.formatting.
    """
    ome_tif = os.path.join(path, 'morphology_mip.ome.tif')
    out_tif = os.path.join(path, 'background.tiff')
    
    if os.path.exists(ome_tif):
        print(f"  > Converting {ome_tif} -> {out_tif}")
        IM = tf.TiffFile(ome_tif)
        position1_series = IM.series[0]
        position1 = position1_series.asarray()
        tf.imwrite(out_tif, position1)
        print("  > Background image conversion done.")
    else:
        print(f"  [WARNING] {ome_tif} not found. Cannot generate background.")


# ================== Main Execution wrapper ==================

def run_step0(config):
    """Main execution function for Step 0 called by pipeline_main.py"""
    input_path = config["input_path"]
    output_dir = config["output_dir"]
    sample_tag = config["sample_tag"]
    
    # Ensure output dir
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n[Step 0-1] Format Xenium Data (Standalone)")
    adata = format_xenium_adata_mid_2023(input_path, sample_tag, output_dir)
    print(f"  > Raw AnnData Shape: {adata.shape}")

    print("\n[Step 0-2] Background Image Analysis (Standalone)")
    if not os.path.exists(os.path.join(input_path, 'background.tiff')):
        format_background(input_path)
    else:
        print("  > background.tiff already exists.")

    # QC & Filtering
    print("\n[Step 0-3] Quality Control & Filtering")
    fmt_params = config.get("formatting", {})
    mincounts = fmt_params.get('mincounts', 10)
    mingenes = fmt_params.get('mingenes', 3)
    
    print(f"  > Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    print(f"  > Filtering cells (min_counts={mincounts}, min_genes={mingenes})...")
    n_cells_before = adata.n_obs
    
    sc.pp.filter_cells(adata, min_counts=mincounts)
    sc.pp.filter_cells(adata, min_genes=mingenes)
    
    n_cells_after = adata.n_obs
    print(f"  > Cells filtered: {n_cells_before} -> {n_cells_after} ({n_cells_before - n_cells_after} removed)")

    # Save the QC'd AnnData
    out_file = os.path.join(output_dir, f"{sample_tag}.h5ad")
    print(f"\n[Step 0-4] Saving QC'd Data to {out_file}")
    adata.write(out_file)
    
    return adata


if __name__ == '__main__':
    # Standalone execution support
    # (Optional: Read config.yaml locally if run directly)
    import yaml
    
    print("Xenium Pipeline Step 0: Formatting (Standalone)")
    
    # Try to load config from same dir
    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = yaml.safe_load(f)
        run_step0(config)
    else:
        # Fallback defaults for independent testing
        print("  [WARNING] config.yaml not found. Using defaults.")
        default_config = {
            "input_path": "/Users/iyeonglag/PycharmProjects/masld-xenium/data/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs",
            "output_dir": "output",
            "sample_tag": "human_alzheimers"
        }
        run_step0(default_config)
