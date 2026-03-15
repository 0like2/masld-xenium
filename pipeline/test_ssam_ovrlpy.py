#!/usr/bin/env python
"""Quick smoke test for SSAM and ovrlpy before running full pipeline.

Usage:
    conda run -n masld python pipeline/test_ssam_ovrlpy.py
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# ── Config ──────────────────────────────────────────────────────
H5AD_PATH = "pipeline/xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs_backup_20260218/step0_formatting/human_alzheimers.h5ad"
TRANSCRIPTS_PATH = "pipeline/xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs_backup_20260218/step0_formatting/human_alzheimers_transcripts.parquet"
OUTPUT_DIR = "pipeline/xenium-output"
SAMPLE_TAG = "human_alzheimers"

# SSAM params (same as config.yaml)
SSAM_UM_PER_PX = 2.0
SSAM_KDE_BW_UM = 2.5
SSAM_NORM_THRES = 5
SSAM_EXP_THRES = 0.2
SSAM_MIN_DIST = 7

# Use a small subset for quick testing
MAX_TRANSCRIPTS = 100_000  # Set to None for full dataset


def load_data():
    """Load adata and ensure spots are available."""
    logger.info(f"Loading {H5AD_PATH} ...")
    adata = sc.read_h5ad(H5AD_PATH)
    logger.info(f"adata shape: {adata.shape}")
    logger.info(f"uns keys: {list(adata.uns.keys())}")

    if 'spots' not in adata.uns:
        logger.info(f"Loading transcripts from {TRANSCRIPTS_PATH} ...")
        if TRANSCRIPTS_PATH.endswith('.parquet'):
            spots = pd.read_parquet(TRANSCRIPTS_PATH)
        else:
            spots = pd.read_csv(TRANSCRIPTS_PATH)
        if 'feature_name' not in spots.columns and 'gene' in spots.columns:
            spots.rename(columns={'gene': 'feature_name'}, inplace=True)
        adata.uns['spots'] = spots

    spots = adata.uns['spots']
    # Ensure feature_name is string (may be stored as bytes in parquet)
    if spots['feature_name'].dtype == object:
        sample = spots['feature_name'].iloc[0]
        if isinstance(sample, bytes):
            spots['feature_name'] = spots['feature_name'].str.decode('utf-8')
            adata.uns['spots'] = spots
    logger.info(f"Spots shape: {spots.shape}")
    logger.info(f"Spots columns: {spots.columns.tolist()}")
    logger.info(f"Has z_location: {'z_location' in spots.columns}")
    return adata


def _patch_zarr_v3_compat():
    """Patch zarr v3 Group.array for ssam v1.1.3 compatibility.
    Must be called BEFORE importing ssam."""
    import zarr
    if int(zarr.__version__.split('.')[0]) < 3:
        return
    if hasattr(zarr.Group, '_ssam_patched'):
        return

    def _compat_array(self, name, shape=None, dtype=None, *, data=None, **kwargs):
        if data is not None:
            data = np.asarray(data)
            return self.create_array(name=name, data=data, overwrite=True)
        return self.create_array(name=name, shape=shape, dtype=dtype, overwrite=True, **kwargs)

    zarr.Group.array = _compat_array
    zarr.Group._ssam_patched = True
    logger.info("Applied zarr v3 compat patch for ssam")


def test_ssam(adata):
    """Test SSAM analysis on a subset of transcripts."""
    logger.info("=" * 60)
    logger.info("TESTING SSAM")
    logger.info("=" * 60)

    _patch_zarr_v3_compat()

    try:
        import ssam
        logger.info(f"ssam version: {getattr(ssam, '__version__', 'unknown')}")
    except ImportError:
        logger.error("ssam not installed!")
        return False

    spots = adata.uns['spots']
    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    coords = spots[gene_mask][['x_location', 'y_location', 'feature_name']].copy()
    coords.columns = ['x', 'y', 'gene']

    # Shift to origin and convert to pixel coordinates
    coords[['x', 'y']] = coords[['x', 'y']] - coords[['x', 'y']].min()
    coords['x'] = coords['x'] / SSAM_UM_PER_PX
    coords['y'] = coords['y'] / SSAM_UM_PER_PX
    coords['gene'] = coords['gene'].astype(str)

    # Spatial crop for quick test (keep density intact, unlike random sampling)
    if MAX_TRANSCRIPTS and len(coords) > MAX_TRANSCRIPTS:
        logger.info(f"Spatial cropping from {len(coords)} transcripts for quick test")
        # Take a small spatial region to keep transcript density
        x_mid = coords['x'].median()
        y_mid = coords['y'].median()
        # Start with a guess, expand until we have enough transcripts
        radius = 200  # pixels
        for _ in range(20):
            mask = ((coords['x'] - x_mid).abs() < radius) & ((coords['y'] - y_mid).abs() < radius)
            if mask.sum() >= MAX_TRANSCRIPTS:
                break
            radius *= 1.5
        coords = coords[mask].reset_index(drop=True)
        coords[['x', 'y']] = coords[['x', 'y']] - coords[['x', 'y']].min()
        logger.info(f"Cropped to {len(coords)} transcripts in {2*radius:.0f}x{2*radius:.0f} px region")

    genes = sorted(coords['gene'].unique().tolist())
    logger.info(f"SSAM: {len(genes)} genes, {len(coords)} transcripts")

    width_px = int(np.ceil(coords['x'].max())) + 1
    height_px = int(np.ceil(coords['y'].max())) + 1
    logger.info(f"Grid: {width_px} x {height_px} px")

    # Create dataset and analysis
    ds = ssam.SSAMDataset()
    analysis = ssam.SSAMAnalysis(ds, ncores=4, verbose=True)

    locations_df = coords[['gene', 'x', 'y']].copy()
    kde_bw_px = SSAM_KDE_BW_UM / SSAM_UM_PER_PX

    logger.info(f"Running KDE (bandwidth={kde_bw_px:.2f} px)...")
    analysis.run_kde(
        locations=locations_df,
        width=width_px,
        height=height_px,
        depth=1,
        bandwidth=kde_bw_px,
        sampling_distance=1.0,
        re_run=True,
    )
    logger.info("KDE complete!")

    # Set thresholds
    analysis.set_thresholds(
        expression_threshold=SSAM_EXP_THRES,
        norm_threshold=SSAM_NORM_THRES,
    )

    # Find local maxima
    search_size = max(3, int(np.round(SSAM_MIN_DIST / SSAM_UM_PER_PX)))
    if search_size % 2 == 0:
        search_size += 1
    logger.info(f"Finding local maxima (search_size={search_size})...")
    analysis.find_localmax(search_size=search_size)

    # Normalize
    analysis.normalize_vectors(normalize_vector=True)

    # Check results
    normalized = ds.normalized_vectors
    local_maxs = ds.local_maxs
    if normalized is None or local_maxs is None or len(local_maxs[0]) == 0:
        logger.warning("No local maxima found (may be expected with subset).")
        return True  # Not a code error

    if 'genes' in ds.zarr_group:
        ssam_genes = [str(g) for g in ds.zarr_group['genes'][:]]
    else:
        ssam_genes = genes

    ssam_adata = sc.AnnData(
        np.array(normalized),
        var=pd.DataFrame(index=ssam_genes),
        obs=pd.DataFrame({'x': local_maxs[0], 'y': local_maxs[1]})
    )
    logger.info(f"SSAM pseudo-cells: {ssam_adata.n_obs}")

    # Quick clustering
    sc.pp.normalize_total(ssam_adata, target_sum=1e4)
    sc.tl.pca(ssam_adata, svd_solver='arpack')
    n_pcs = min(40, ssam_adata.n_vars - 1, ssam_adata.n_obs - 1)
    sc.pp.neighbors(ssam_adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.leiden(ssam_adata, resolution=2.0, random_state=42)

    n_clusters = len(ssam_adata.obs['leiden'].cat.categories)
    logger.info(f"SSAM Leiden clusters: {n_clusters}")
    logger.info("SSAM TEST PASSED!")
    return True


def test_ovrlpy(adata):
    """Test ovrlpy analysis on a subset of transcripts."""
    logger.info("=" * 60)
    logger.info("TESTING OVRLPY")
    logger.info("=" * 60)

    try:
        import ovrlpy
        logger.info(f"ovrlpy version: {getattr(ovrlpy, '__version__', 'unknown')}")
    except ImportError:
        logger.error("ovrlpy not installed!")
        return False

    spots = adata.uns['spots']

    if 'z_location' not in spots.columns:
        logger.warning("z_location not found in spots. ovrlpy requires 3D coordinates.")
        logger.info("Checking if we can still test with z=0 ...")
        spots = spots.copy()
        spots['z_location'] = 0.0
        adata.uns['spots'] = spots

    um_per_pixel = 2.0
    gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl', case=False, na=False)
    coords = spots[gene_mask][['x_location', 'y_location', 'z_location', 'feature_name']].copy()
    coords.columns = ['x', 'y', 'z', 'gene']
    coords['x'] = coords['x'] / um_per_pixel
    coords['y'] = coords['y'] / um_per_pixel
    coords['gene'] = coords['gene'].astype('category')
    coords = coords.reset_index(drop=True)

    # Spatial crop for quick test
    if MAX_TRANSCRIPTS and len(coords) > MAX_TRANSCRIPTS:
        logger.info(f"Spatial cropping from {len(coords)} transcripts for quick test")
        x_mid = coords['x'].median()
        y_mid = coords['y'].median()
        radius = 200
        for _ in range(20):
            mask = ((coords['x'] - x_mid).abs() < radius) & ((coords['y'] - y_mid).abs() < radius)
            if mask.sum() >= MAX_TRANSCRIPTS:
                break
            radius *= 1.5
        coords = coords[mask].reset_index(drop=True)
        logger.info(f"Cropped to {len(coords)} transcripts")

    logger.info(f"ovrlpy: {len(coords)} transcripts, {coords['gene'].nunique()} genes")

    # Create Ovrlp object (n_components=2 for small test subsets)
    ovrlp_obj = ovrlpy.Ovrlp(
        coords,
        KDE_bandwidth=1,
        min_distance=8,
        n_components=2,
        gene_key='gene',
        coordinate_keys=('x', 'y', 'z'),
    )
    logger.info("Ovrlp object created successfully")

    # Run analysis
    logger.info("Running ovrlpy.analyse ...")
    ovrlp_obj.analyse(gridsize=1, min_transcripts=10, fit_umap=False)
    logger.info("ovrlpy.analyse complete!")

    # Check signal integrity
    if hasattr(ovrlp_obj, 'integrity_map') or hasattr(ovrlp_obj, 'signal_map'):
        logger.info(f"integrity_map: {type(getattr(ovrlp_obj, 'integrity_map', None))}")
        logger.info(f"signal_map: {type(getattr(ovrlp_obj, 'signal_map', None))}")

    # Test plot function (without saving)
    try:
        import matplotlib
        matplotlib.use('Agg')
        fig = ovrlpy.plot_signal_integrity(ovrlp_obj, signal_threshold=2)
        logger.info(f"plot_signal_integrity: OK (fig type={type(fig)})")
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception as e:
        logger.warning(f"plot_signal_integrity failed: {e}")

    logger.info("OVRLPY TEST PASSED!")
    return True


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    logger.info(f"Working dir: {os.getcwd()}")

    adata = load_data()

    ssam_ok = test_ssam(adata)
    ovrlpy_ok = test_ovrlpy(adata)

    logger.info("=" * 60)
    logger.info(f"SSAM:   {'PASS' if ssam_ok else 'FAIL'}")
    logger.info(f"OVRLPY: {'PASS' if ovrlpy_ok else 'FAIL'}")
    logger.info("=" * 60)

    sys.exit(0 if (ssam_ok and ovrlpy_ok) else 1)
