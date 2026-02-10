# Step 0: Formatting (Ref: notebooks/0_formatting/)
# Converts raw Xenium machine output → AnnData (.h5ad)
#
# Flow:
#   0-1. Decompress raw Xenium files
#   0-2. Build AnnData from cell_feature_matrix
#   0-3. Flag & remove control probes
#   0-4. Parse gene_panel.json → Ensembl IDs
#   0-5. Load transcripts (parquet/csv) → sidecar
#   0-6. Import 10X analysis (UMAP, PCA, clusters)
#   0-7. Convert dtypes & save formatted .h5ad
#   0-8. Format background DAPI image
#   0-9. QC filtering (min_counts, min_genes)

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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def format_xenium_adata_mid_2023(path, tag, output_path, filter_nuclei_only=False):
    """
    Format xenium data (output from the machine) to adata format,
    considering the format used by Xenium at Q2 2023 (Mid 2023).
    """
    print(f"  > Processing Xenium Data at: {path}")

    # --- 0-1. Decompress cell_feature_matrix ---

    if not os.path.exists(os.path.join(path, 'cell_feature_matrix')):
        tar_gz = os.path.join(path, 'cell_feature_matrix.tar.gz')
        if os.path.exists(tar_gz):
            print(f"  > Decompressing {tar_gz}...")
            shutil.unpack_archive(tar_gz, path)
            print("  > Decompression done.")
        else:
            print(f"  [WARNING] 'cell_feature_matrix' folder and tar.gz not found at {path}")

    for filename in ['transcripts.csv', 'cells.csv']:
        file_path = os.path.join(path, filename)
        gz_path = file_path + '.gz'
        if not os.path.exists(file_path) and os.path.exists(gz_path):
            print(f"  > Decompressing {gz_path}...")
            with gzip.open(gz_path, 'rb') as f_in:
                with open(file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

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

    # --- 0-2. Build AnnData from cell_feature_matrix ---

    print("  > Reading matrix and metadata...")
    try:
        matrix_path = os.path.join(cfm_path, 'matrix.mtx')
        a = mmread(matrix_path)
        ad = a.todense()

        cells_path = os.path.join(path, "cells.csv")
        cell_info = pd.read_csv(cells_path)

        features_path = os.path.join(cfm_path, 'features.tsv')

        # Dynamic column handling: 3-col (id/name/reason) or 2-col (id/reason)
        try:
            features_df = pd.read_csv(features_path, header=None, sep='\t')

            if features_df.shape[1] == 3:
                features_df.columns = ['gene_id', 'gene_name', 'reason_of_inclusion']
                if features_df['gene_name'].is_unique:
                    features_df.index = features_df['gene_name']
                    features_df.index.name = 'index'
                else:
                     features_df.index = features_df['gene_id']
                     features_df.index.name = 'index'

            elif features_df.shape[1] == 2:
                features_df.columns = ['gene_id', 'reason_of_inclusion']
                features_df['gene_name'] = features_df['gene_id']
                features_df.index = features_df['gene_id']
                features_df.index.name = 'index'
            else:
                print(f"    [WARNING] Unexpected number of columns in features.tsv: {features_df.shape[1]}. Proceeding with default indexing.")
                features_df.index.name = 'index'

            features = features_df

        except Exception as e:
            print(f"    [ERROR] Failed to read/parse features.tsv: {e}")
            raise e

        adata = sc.AnnData(ad.transpose(), obs=cell_info, var=features)

        # --- 0-3. Flag & remove control probes ---

        print("  > Identifying and handling Negative Controls (QC)...")

        control_mask = (
            adata.var['gene_id'].str.contains('NegControlProbe_', case=False) |
            adata.var['gene_id'].str.contains('NegControlCodeword_', case=False) |
            adata.var['gene_id'].str.contains('antisense_', case=False) |
            adata.var['gene_id'].str.contains('BLANK', case=False)
        )

        adata.var['is_control'] = control_mask
        n_controls = control_mask.sum()
        print(f"  > Found {n_controls} negative control features.")

        # Calculate QC metrics before filtering to track noise levels per cell
        adata.obs['total_counts_raw'] = np.array(adata.X.sum(axis=1)).flatten()

        if n_controls > 0:
            control_genes = adata[:, control_mask]
            adata.obs['total_counts_control'] = np.array(control_genes.X.sum(axis=1)).flatten()

            with np.errstate(divide='ignore', invalid='ignore'):
                adata.obs['pct_counts_control'] = (adata.obs['total_counts_control'] / adata.obs['total_counts_raw']) * 100
            adata.obs['pct_counts_control'] = adata.obs['pct_counts_control'].fillna(0.0)

            total_control_reads = adata.obs['total_counts_control'].sum()
            total_raw_reads = adata.obs['total_counts_raw'].sum()
            if total_raw_reads > 0:
                global_control_pct = (total_control_reads / total_raw_reads) * 100
                print(f"  > Global % Control Reads: {global_control_pct:.2f}%")
        else:
            adata.obs['total_counts_control'] = 0.0
            adata.obs['pct_counts_control'] = 0.0

        # Filter out control probes
        print("  > Filtering out negative control probes from AnnData object...")
        original_shape = adata.shape
        adata = adata[:, ~adata.var['is_control']].copy()
        print(f"  > Filtered controls: {original_shape} -> {adata.shape}")

        if 'x_centroid' in cell_info.columns and 'y_centroid' in cell_info.columns:
             adata.obsm['spatial'] = np.array(cell_info[['x_centroid', 'y_centroid']])
             print("  > Spatial coordinates (centroids) added to obsm['spatial']")
    except Exception as e:
        raise ValueError(f"Error reading matrix files: {e}")

    # --- 0-4. Parse gene_panel.json → Ensembl IDs ---

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

    # --- 0-5. Load transcripts (parquet/csv) → sidecar ---

    print("  > Loading transcripts (this may take a while)...")
    transcripts_parquet = os.path.join(path, 'transcripts.parquet')
    transcripts_csv = os.path.join(path, 'transcripts.csv')
    transcripts_csv_gz = os.path.join(path, 'transcripts.csv.gz')

    transcripts = None
    if os.path.exists(transcripts_parquet):
        print(f"  > Reading transcripts from Parquet: {transcripts_parquet}")
        transcripts = pd.read_parquet(transcripts_parquet)
    elif os.path.exists(transcripts_csv):
        print(f"  > Reading transcripts from CSV: {transcripts_csv}")
        transcripts = pd.read_csv(transcripts_csv, low_memory=False)
    elif os.path.exists(transcripts_csv_gz):
        print(f"  > Reading transcripts from GZ CSV: {transcripts_csv_gz}")
        transcripts = pd.read_csv(transcripts_csv_gz, compression='gzip', low_memory=False)

    if transcripts is not None:
        print(f"  > Transcripts shape: {transcripts.shape}")

        if 'cell_id' not in transcripts.columns and transcripts.index.name == 'cell_id':
            transcripts = transcripts.reset_index()

        # Nuclei-only filtering (matches notebook 0_3 overlaps_nucleus=1)
        if filter_nuclei_only:
            if 'overlaps_nucleus' in transcripts.columns:
                n_before = len(transcripts)
                transcripts = transcripts[transcripts['overlaps_nucleus'] == 1]
                n_after = len(transcripts)
                print(f"  > Nuclei filtering: {n_before} -> {n_after} transcripts ({n_before - n_after} cytoplasmic removed)")
            else:
                print(f"  [WARNING] 'overlaps_nucleus' column not found. Skipping nuclei filter.")

        # Save as sidecar parquet (not adata.uns to avoid H5AD bloat)
        spots_path = os.path.join(output_path, f"{tag}_transcripts.parquet")
        print(f"  > Saving transcripts to independent file: {spots_path}")
        transcripts.to_parquet(spots_path)

        adata.uns['spots_path'] = spots_path
        print("  > Stored 'spots_path' in adata.uns")
    else:
        print("  [WARNING] No transcripts file found (checked parquet/csv/gz).")

    # --- 0-6. Import 10X analysis (UMAP, PCA, clusters) ---

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
            umap_path = os.path.join(analysis_base, 'umap/gene_expression_2_components/projection.csv')
            if os.path.exists(umap_path):
                UMAP = pd.read_csv(umap_path, index_col=0)
                adata.obsm['X_umap'] = np.array(UMAP)

            tsne_path = os.path.join(analysis_base, 'tsne/gene_expression_2_components/projection.csv')
            if os.path.exists(tsne_path):
                TSNE = pd.read_csv(tsne_path, index_col=0)
                adata.obsm['X_tsne'] = np.array(TSNE)

            pca_path = os.path.join(analysis_base, 'PCA/gene_expression_10_components/projection.csv')
            if os.path.exists(pca_path):
                PCA = pd.read_csv(pca_path, index_col=0)
                adata.obsm['X_pca'] = np.array(PCA)

            clusters_path = os.path.join(analysis_base, 'clustering/gene_expression_graphclust/clusters.csv')
            if os.path.exists(clusters_path):
                clusters = pd.read_csv(clusters_path, index_col=0)
                adata.obs['graph_clusters'] = list(clusters['Cluster'].astype(str))

            for k in range(2, 11):
                km_path = os.path.join(analysis_base, f'clustering/gene_expression_kmeans_{k}_clusters/clusters.csv')
                if os.path.exists(km_path):
                    km = pd.read_csv(km_path, index_col=0)
                    adata.obs[f'kmeans{k}_clusters'] = list(km['Cluster'].astype(str))

        except Exception as e:
            print(f"  [WARNING] Analysis files load issue: {e}")

    # --- 0-7. Convert dtypes & save formatted .h5ad ---

    adata.X = sp.csr_matrix(adata.X)

    if 'spots' in adata.uns:
        del adata.uns['spots']

    # Convert object-dtype columns to category/str for h5py compatibility
    for df in [adata.obs, adata.var]:
        for col in df.columns:
            if df[col].dtype == 'object':
                try:
                    df[col] = df[col].astype(str).astype('category')
                except:
                    print(f"  [WARNING] Could not convert column {col} to category. Converting to string.")
                    df[col] = df[col].astype(str)

    keys_to_remove = []
    for k, v in adata.uns.items():
        if isinstance(v, dict):
            continue
        if hasattr(v, 'dtype'):
             if v.dtype == 'object':
                 try:
                     adata.uns[k] = v.astype(str)
                 except:
                     pass

    out_file = os.path.join(output_path, f"{tag}.h5ad")
    print(f"  > Saving to {out_file}...")
    adata.write(out_file)
    return adata


def format_background(path):
    """Convert OME-TIFF background MIP image to plain .tiff."""
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


# --- Entry point ---

def run_step0(config):
    """Main execution function for Step 0 called by pipeline_main.py"""
    input_path = config["input_path"]
    output_dir = config["output_dir"]
    sample_tag = config["sample_tag"]

    os.makedirs(output_dir, exist_ok=True)

    fmt_params = config.get("formatting", {})
    filter_nuclei = fmt_params.get('filter_nuclei_only', False)

    print("\n[Step 0-1] Format Xenium Data")
    adata = format_xenium_adata_mid_2023(input_path, sample_tag, output_dir, filter_nuclei_only=filter_nuclei)
    print(f"  > Raw AnnData Shape: {adata.shape}")

    # --- 0-8. Format background DAPI image ---
    print("\n[Step 0-8] Background Image")
    if not os.path.exists(os.path.join(input_path, 'background.tiff')):
        format_background(input_path)
    else:
        print("  > background.tiff already exists.")

    # --- 0-9. QC filtering (min_counts, min_genes) ---

    print("\n[Step 0-9] Quality Control & Filtering")
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

    out_file = os.path.join(output_dir, f"{sample_tag}.h5ad")
    print(f"\n  > Saving QC'd data to {out_file}")
    adata.write(out_file)

    return adata


if __name__ == '__main__':
    import yaml

    print("Xenium Pipeline Step 0: Formatting (Standalone)")

    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = yaml.safe_load(f)
        run_step0(config)
    else:
        print("  [WARNING] config.yaml not found. Using defaults.")
        default_config = {
            "input_path": "/Users/iyeonglag/PycharmProjects/masld-xenium/data/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs",
            "output_dir": "output",
            "sample_tag": "human_alzheimers"
        }
        run_step0(default_config)
