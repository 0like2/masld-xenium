# Step 6: Segmentation Benchmark (Ref: notebooks/5_segmentation_benchmark/)
# Benchmarks multiple segmentation methods (nuclei, Cellpose, expansion, Baysor)
# by concatenating results, preprocessing, and computing comparison metrics.
#
# Flow:
#   6-1. Run Baysor (optional)
#   6-2. Load segmentation results (nuclei/cellpose/expansion/baysor)
#   6-3. Concatenate & preprocess (normalize/PCA/Leiden)
#   6-4. Annotation transfer (majority voting from reference)
#   6-5. Benchmark metrics (n_cells, median_reads, assigned_prop)
#   6-6. Visualizations
#     6-6a. UMAP per method
#     6-6b. Spatial scatter
#     6-6c. Cell type barplot
#     6-6d. Counts violin

import os
import json
import logging
import shutil
import multiprocessing
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from scipy.sparse import issparse
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from benchmark_utils import metrics

logger = logging.getLogger("Step6_Benchmark")
logging.basicConfig(level=logging.INFO)


# --- 6-1. Baysor data preparation (helper) ---

def prep_xenium_data_for_baysor(xenium_dir, out_dir, crop=False, coords=None):
    """Format Xenium transcripts into Baysor-compatible CSV (gene, x, y)."""
    out_dir = Path(out_dir)
    xenium_dir = Path(xenium_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Preparing Xenium data for Baysor in {out_dir}...")

    transcripts_path = xenium_dir / "transcripts.csv"
    if not transcripts_path.exists():
        transcripts_path = xenium_dir / "transcripts.parquet"
        if transcripts_path.exists():
            spots = pd.read_parquet(transcripts_path)
            spots = _decode_bytes_columns(spots)
        else:
            logger.error(f"Transcripts file not found in {xenium_dir}")
            return False
    else:
        spots = pd.read_csv(transcripts_path)

    logger.info(f"Loaded {len(spots)} transcripts.")

    col_map = {
        'feature_name': 'gene',
        'x_location': 'x',
        'y_location': 'y',
        'z_location': 'z'
    }

    spots_baysor = spots.rename(columns=col_map)

    req_cols = ['gene', 'x', 'y']
    if not all(col in spots_baysor.columns for col in req_cols):
        logger.error(f"Missing required columns for Baysor. Found: {spots_baysor.columns}")
        col_map_alt = {'feature_name': 'gene', 'x': 'x', 'y': 'y'}
        spots_baysor = spots.rename(columns=col_map_alt)
        if not all(col in spots_baysor.columns for col in req_cols):
            return False

    # Apply crop filter if coordinates are provided
    if crop and coords is not None:
        x_min, x_max, y_min, y_max = coords
        mask = spots_baysor['x'].between(x_min, x_max) & spots_baysor['y'].between(y_min, y_max)
        n_before = len(spots_baysor)
        spots_baysor = spots_baysor[mask]
        logger.info(f"Cropped transcripts: {n_before} -> {len(spots_baysor)} "
                     f"(region [{x_min:.0f}, {x_max:.0f}] x [{y_min:.0f}, {y_max:.0f}] µm)")

    spots_file = out_dir / "spots.csv"
    spots_baysor[req_cols].to_csv(spots_file, index=False)
    logger.info(f"Saved spots to {spots_file}")

    return True


# --- 6-1a. Crop / Tiling helpers ---

def find_densest_region(spots_df, size_um=2000, bin_um=50):
    """Find the densest square region using 2D histogram + prefix-sum sliding window.

    Parameters
    ----------
    spots_df : DataFrame with 'x' and 'y' columns (µm coordinates)
    size_um  : side length of the square crop region in µm
    bin_um   : bin size for the 2D histogram (smaller = more precise but slower)

    Returns
    -------
    (x_min, x_max, y_min, y_max) of the densest region
    """
    x = spots_df['x'].values
    y = spots_df['y'].values

    x_min_data, x_max_data = x.min(), x.max()
    y_min_data, y_max_data = y.min(), y.max()

    # Build 2D histogram
    x_bins = np.arange(x_min_data, x_max_data + bin_um, bin_um)
    y_bins = np.arange(y_min_data, y_max_data + bin_um, bin_um)
    hist, _, _ = np.histogram2d(x, y, bins=[x_bins, y_bins])

    # Prefix sum for O(1) window queries
    prefix = np.cumsum(np.cumsum(hist, axis=0), axis=1)

    window_bins = int(np.ceil(size_um / bin_um))
    nx, ny = hist.shape

    if window_bins >= nx or window_bins >= ny:
        logger.warning("Crop size >= tissue extent; using full tissue bounds.")
        return (x_min_data, x_max_data, y_min_data, y_max_data)

    best_count = -1
    best_i, best_j = 0, 0

    for i in range(nx - window_bins + 1):
        for j in range(ny - window_bins + 1):
            # Sum in rectangle [i, i+window_bins) x [j, j+window_bins)
            i2, j2 = i + window_bins - 1, j + window_bins - 1
            total = prefix[i2, j2]
            if i > 0:
                total -= prefix[i - 1, j2]
            if j > 0:
                total -= prefix[i2, j - 1]
            if i > 0 and j > 0:
                total += prefix[i - 1, j - 1]

            if total > best_count:
                best_count = total
                best_i, best_j = i, j

    crop_x_min = x_bins[best_i]
    crop_x_max = x_bins[min(best_i + window_bins, len(x_bins) - 1)]
    crop_y_min = y_bins[best_j]
    crop_y_max = y_bins[min(best_j + window_bins, len(y_bins) - 1)]

    logger.info(f"Densest {size_um}×{size_um} µm region: "
                f"x=[{crop_x_min:.0f}, {crop_x_max:.0f}], "
                f"y=[{crop_y_min:.0f}, {crop_y_max:.0f}] "
                f"({int(best_count)} molecules in {bin_um}µm bins)")

    return (crop_x_min, crop_x_max, crop_y_min, crop_y_max)


def _build_baysor_cmd(executable, params, spots_file, output_csv, prior_tif=None):
    """Build the Baysor CLI command list (shared by crop/tile/full modes)."""
    cmd = [executable, "run"]
    if "scale" in params:
        cmd.extend(["-s", str(params["scale"])])
    if "min_molecules_per_cell" in params:
        cmd.extend(["-m", str(params["min_molecules_per_cell"])])
    if "prior_segmentation_confidence" in params:
        cmd.extend(["--prior-segmentation-confidence",
                     str(params["prior_segmentation_confidence"])])
    cmd.extend(["-o", str(output_csv)])
    cmd.append(str(spots_file))
    if prior_tif and os.path.exists(prior_tif):
        cmd.append(str(prior_tif))
    return cmd


def _run_baysor_single(executable, params, spots_file, output_dir, tile_label="",
                       prior_tif=None, dry_run=False):
    """Run a single Baysor subprocess. Returns path to segmentation CSV or None."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    segmentation_csv = output_dir / "segmentation.csv"

    label = f"[{tile_label}] " if tile_label else ""

    # Cache check: skip re-run if output already exists
    if segmentation_csv.exists() and segmentation_csv.stat().st_size > 100:
        logger.info(f"{label}Baysor output already exists ({segmentation_csv}). Skipping re-run.")
        return str(segmentation_csv)

    cmd = _build_baysor_cmd(executable, params, spots_file, segmentation_csv,
                            prior_tif=prior_tif)

    logger.info(f"{label}Baysor command: {' '.join(cmd)}")

    if dry_run:
        logger.warning(f"{label}dry_run=True: writing minimal placeholder CSV.")
        if not segmentation_csv.exists():
            segmentation_csv.write_text("transcript_id,cell\n")
        return str(segmentation_csv)

    try:
        subprocess.run(cmd, check=True)
        logger.info(f"{label}Baysor completed successfully.")
        return str(segmentation_csv)
    except subprocess.CalledProcessError as e:
        logger.error(f"{label}Baysor failed (rc={e.returncode}): {e}")
        return None
    except FileNotFoundError:
        logger.error(f"{label}Baysor binary not found: {executable}")
        return None
    except Exception as e:
        logger.error(f"{label}Unexpected Baysor error: {e}")
        return None


def _decode_bytes_columns(df):
    """Decode any bytes-typed object columns to UTF-8 strings.

    Xenium parquet files may store string columns (feature_name, cell_id)
    as bytes, which causes downstream gene-name mismatches when compared
    with AnnData var_names loaded from h5ad or CSV sources.
    """
    for col in df.columns:
        if df[col].dtype == object and len(df) > 0 and isinstance(df[col].iloc[0], bytes):
            df[col] = df[col].str.decode('utf-8')
    return df


def _load_spots_df(xenium_dir):
    """Load transcripts as a DataFrame with gene/x/y columns (parquet preferred)."""
    xenium_dir = Path(xenium_dir)
    parquet_path = xenium_dir / "transcripts.parquet"
    csv_path = xenium_dir / "transcripts.csv"

    if parquet_path.exists():
        spots = pd.read_parquet(parquet_path)
    elif csv_path.exists():
        spots = pd.read_csv(csv_path)
    else:
        logger.error(f"No transcripts file found in {xenium_dir}")
        return None

    spots = _decode_bytes_columns(spots)

    col_map = {'feature_name': 'gene', 'x_location': 'x', 'y_location': 'y'}
    spots = spots.rename(columns=col_map)

    if not all(c in spots.columns for c in ['gene', 'x', 'y']):
        logger.error(f"Missing gene/x/y columns after rename. Found: {list(spots.columns)}")
        return None

    logger.info(f"Loaded {len(spots):,} transcripts from {xenium_dir}")
    return spots


def _run_baysor_crop(config, spots_df, baysor_out_dir, executable, params, dry_run):
    """Run Baysor on a cropped region (auto-detect or manual coords).

    Returns path to segmentation CSV or None.
    """
    crop_conf = config.get("baysor", {}).get("crop", {})
    size_um = crop_conf.get("size_um", 2000)
    manual_coords = crop_conf.get("coords")  # null or [x_min, x_max, y_min, y_max]

    # Determine crop coordinates
    if manual_coords is not None:
        coords = tuple(manual_coords)
        logger.info(f"Using manual crop coords: {coords}")
    else:
        coords = find_densest_region(spots_df, size_um=size_um)
        logger.info(f"Auto-detected densest region: {coords}")

    x_min, x_max, y_min, y_max = coords

    # Crop spots
    mask = spots_df['x'].between(x_min, x_max) & spots_df['y'].between(y_min, y_max)
    spots_crop = spots_df.loc[mask, ['gene', 'x', 'y']]
    logger.info(f"Cropped molecules: {len(spots_df):,} -> {len(spots_crop):,}")

    if len(spots_crop) == 0:
        logger.error("No molecules in crop region!")
        return None

    # Save cropped spots
    crop_dir = os.path.join(baysor_out_dir, "crop")
    os.makedirs(crop_dir, exist_ok=True)
    spots_file = os.path.join(crop_dir, "spots.csv")
    spots_crop.to_csv(spots_file, index=False)

    # Save crop metadata
    metadata = {
        "mode": "crop",
        "coords": list(coords),
        "size_um": size_um,
        "n_molecules_total": len(spots_df),
        "n_molecules_crop": len(spots_crop),
        "auto_detected": manual_coords is None,
    }
    meta_path = os.path.join(crop_dir, "crop_metadata.json")
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)
    logger.info(f"Saved crop metadata to {meta_path}")

    # Run Baysor (no prior TIF in crop mode — coordinate mismatch)
    return _run_baysor_single(
        executable, params, spots_file, crop_dir,
        tile_label="crop", prior_tif=None, dry_run=dry_run
    )


def compute_tile_coords(x_min, x_max, y_min, y_max, tile_size_um, overlap_um):
    """Compute tile coordinates with overlap for tiled Baysor execution.

    Returns list of dicts: {tile_id, bbox: (x0, x1, y0, y1), core: (x0, x1, y0, y1)}
    The 'core' region excludes overlap/2 from each side (for deduplication).
    """
    stride = tile_size_um - overlap_um
    half_ovlp = overlap_um / 2.0
    tiles = []
    tile_id = 0

    x = x_min
    while x < x_max:
        y = y_min
        x_end = min(x + tile_size_um, x_max)
        while y < y_max:
            y_end = min(y + tile_size_um, y_max)

            # Core region: shrink by half_ovlp except at tissue edges
            core_x0 = x if x == x_min else x + half_ovlp
            core_x1 = x_end if x_end >= x_max else x_end - half_ovlp
            core_y0 = y if y == y_min else y + half_ovlp
            core_y1 = y_end if y_end >= y_max else y_end - half_ovlp

            tiles.append({
                "tile_id": tile_id,
                "bbox": (x, x_end, y, y_end),
                "core": (core_x0, core_x1, core_y0, core_y1),
            })
            tile_id += 1
            y += stride
        x += stride

    return tiles


def _run_baysor_tile_worker(args):
    """Worker function for multiprocessing.Pool — runs Baysor on a single tile."""
    (executable, params, spots_file, tile_out_dir, tile_label, dry_run) = args
    return _run_baysor_single(
        executable, params, spots_file, tile_out_dir,
        tile_label=tile_label, prior_tif=None, dry_run=dry_run
    )


def _merge_tile_results(tile_results, baysor_out_dir):
    """Merge per-tile Baysor segmentation CSVs, deduplicating overlap regions.

    Parameters
    ----------
    tile_results : list of dicts with keys:
        csv_path, core (x0, x1, y0, y1), tile_id
    baysor_out_dir : output directory for merged CSV

    Returns
    -------
    Path to merged segmentation.csv
    """
    merged_parts = []
    cell_id_offset = 0

    for tr in tile_results:
        csv_path = tr["csv_path"]
        if csv_path is None or not os.path.exists(csv_path):
            logger.warning(f"Tile {tr['tile_id']}: missing segmentation CSV, skipping.")
            continue

        df = pd.read_csv(csv_path)
        if len(df) == 0:
            continue

        # Filter to core region (deduplicate overlaps)
        core_x0, core_x1, core_y0, core_y1 = tr["core"]
        if 'x' in df.columns and 'y' in df.columns:
            mask = (df['x'].between(core_x0, core_x1) &
                    df['y'].between(core_y0, core_y1))
            df = df[mask]

        # Offset cell IDs to ensure global uniqueness
        if 'cell' in df.columns:
            numeric_cells = pd.to_numeric(df['cell'], errors='coerce')
            valid = numeric_cells.notna() & (numeric_cells > 0)
            df.loc[valid, 'cell'] = (numeric_cells[valid] + cell_id_offset).astype(int)
            if valid.any():
                cell_id_offset = int(numeric_cells[valid].max()) + cell_id_offset

        merged_parts.append(df)
        logger.info(f"Tile {tr['tile_id']}: {len(df)} transcripts after core filter.")

    if not merged_parts:
        logger.error("No tile results to merge!")
        return None

    merged = pd.concat(merged_parts, ignore_index=True)
    out_path = os.path.join(baysor_out_dir, "segmentation.csv")
    merged.to_csv(out_path, index=False)
    logger.info(f"Merged {len(merged):,} transcripts from {len(merged_parts)} tiles -> {out_path}")
    return out_path


def _run_baysor_tiled(config, spots_df, baysor_out_dir, executable, params, dry_run):
    """Run Baysor in tiled parallel mode across the full tissue.

    Returns path to merged segmentation CSV or None.
    """
    tile_conf = config.get("baysor", {}).get("tiling", {})
    tile_size = tile_conf.get("tile_size_um", 2000)
    overlap = tile_conf.get("overlap_um", 200)
    n_workers = tile_conf.get("n_workers", 4)

    x_min, x_max = spots_df['x'].min(), spots_df['x'].max()
    y_min, y_max = spots_df['y'].min(), spots_df['y'].max()

    tiles = compute_tile_coords(x_min, x_max, y_min, y_max, tile_size, overlap)
    logger.info(f"Tiling: {len(tiles)} tiles ({tile_size}µm, {overlap}µm overlap, {n_workers} workers)")

    # Prepare per-tile spots files
    worker_args = []
    tile_meta = []
    for tile in tiles:
        tid = tile["tile_id"]
        bx0, bx1, by0, by1 = tile["bbox"]

        tile_spots = spots_df.loc[
            spots_df['x'].between(bx0, bx1) & spots_df['y'].between(by0, by1),
            ['gene', 'x', 'y']
        ]

        if len(tile_spots) == 0:
            logger.info(f"Tile {tid}: empty, skipping.")
            continue

        tile_dir = os.path.join(baysor_out_dir, f"tile_{tid:04d}")
        os.makedirs(tile_dir, exist_ok=True)
        spots_file = os.path.join(tile_dir, "spots.csv")
        tile_spots.to_csv(spots_file, index=False)

        worker_args.append((executable, params, spots_file, tile_dir,
                            f"tile_{tid:04d}", dry_run))
        tile_meta.append({
            "tile_id": tid,
            "core": tile["core"],
            "tile_dir": tile_dir,
        })

    if not worker_args:
        logger.error("No non-empty tiles found!")
        return None

    # Run tiles in parallel (or sequentially if n_workers=1)
    if n_workers <= 1 or dry_run:
        results = [_run_baysor_tile_worker(a) for a in worker_args]
    else:
        with multiprocessing.Pool(processes=n_workers) as pool:
            results = pool.map(_run_baysor_tile_worker, worker_args)

    # Assemble tile results for merge
    tile_results = []
    for meta, csv_path in zip(tile_meta, results):
        tile_results.append({
            "tile_id": meta["tile_id"],
            "core": meta["core"],
            "csv_path": csv_path,
        })

    return _merge_tile_results(tile_results, baysor_out_dir)


# --- 6-1. Baysor execution ---

def run_baysor(config, xenium_input_dir):
    """Execute Baysor segmentation with mode dispatch: crop > tiling > full.

    Mode priority:
      1. crop (default) — densest region only, fast (~5-10h)
      2. tiling — parallel tiles across full tissue (future)
      3. full — original whole-tissue single run (3d+)
    """
    logger.info("--- Starting Baysor Segmentation Phase ---")

    baysor_conf = config.get("baysor", {})
    if not baysor_conf.get("enabled", False):
        logger.info("Baysor disabled in config.")
        return None

    executable = baysor_conf.get("executable_path", "baysor")
    params = baysor_conf.get("params", {})
    dry_run = baysor_conf.get("dry_run", False)

    # Early check: skip data prep if Baysor binary is missing (unless dry_run)
    if not dry_run and not shutil.which(executable):
        logger.error(
            f"Baysor binary not found at '{executable}'. "
            "Install Baysor or set baysor.executable_path in config. "
            "To skip execution, set baysor.dry_run=True."
        )
        return None

    step6_out = config.get("output_dir", "xenium-output/step6_benchmark")
    baysor_out_dir = os.path.join(step6_out, "baysor_run")
    os.makedirs(baysor_out_dir, exist_ok=True)

    # Determine mode: crop > tiling > full
    crop_enabled = baysor_conf.get("crop", {}).get("enabled", False)
    tiling_enabled = baysor_conf.get("tiling", {}).get("enabled", False)

    if crop_enabled:
        mode = "crop"
    elif tiling_enabled:
        mode = "tiling"
    else:
        mode = "full"

    logger.info(f"Baysor mode: {mode}")

    # --- crop / tiling modes require transcript loading ---
    if mode in ("crop", "tiling"):
        spots_df = _load_spots_df(xenium_input_dir)
        if spots_df is None:
            logger.error("Failed to load transcripts for crop/tiling.")
            return None

        if mode == "crop":
            return _run_baysor_crop(config, spots_df, baysor_out_dir,
                                    executable, params, dry_run)
        else:  # tiling
            return _run_baysor_tiled(config, spots_df, baysor_out_dir,
                                     executable, params, dry_run)

    # --- full mode: original whole-tissue logic ---
    prep_dir = os.path.join(baysor_out_dir, "prep")
    if not prep_xenium_data_for_baysor(xenium_input_dir, prep_dir):
        logger.error("Failed to prepare data for Baysor.")
        return None

    # Prior segmentation TIF (from Step 3 or legacy input_advanced)
    prior_seg_path = baysor_conf.get("prior_segmentation_tif")
    if not prior_seg_path:
        val = config.get("benchmark", {}).get("input_advanced")
        if str(val).endswith('.tif'):
            prior_seg_path = val

    if prior_seg_path and str(prior_seg_path).endswith('.tif') and os.path.exists(prior_seg_path):
        logger.info(f"Using prior segmentation TIF: {prior_seg_path}")
    else:
        logger.warning(
            f"Prior segmentation TIF not valid or missing (Path: {prior_seg_path}). "
            "Running Baysor without prior."
        )
        prior_seg_path = None

    spots_file = os.path.join(prep_dir, "spots.csv")
    return _run_baysor_single(
        executable, params, spots_file, baysor_out_dir,
        tile_label="full", prior_tif=prior_seg_path, dry_run=dry_run
    )


# --- 6-2. Transcript aggregation (helper) ---

def load_transcripts_as_adata(csv_path, cell_col='cell', gene_col='gene',
                              current_method_name='unknown', ref_adata_path=None):
    """Load a transcripts CSV and aggregate into AnnData (cells x genes count matrix)."""
    logger.info(f"Loading transcripts from {csv_path} for method '{current_method_name}'...")
    try:
        df = pd.read_csv(csv_path)

        # Strip b'...' encoding artifacts from string columns (written by
        # parquet→CSV roundtrip where bytes were not decoded)
        for col in df.columns:
            if df[col].dtype == object and len(df) > 0:
                sample = df[col].iloc[0]
                if isinstance(sample, str) and sample.startswith("b'") and sample.endswith("'"):
                    df[col] = df[col].str[2:-1]

        # Resolve column names — prefer expansion assignment (closest_cell)
        # over nuclei (in_cell) over user-specified over defaults
        cols = df.columns
        actual_cell_col = None
        for candidate in ['closest_cell', 'in_cell', cell_col, 'cell_id', 'cell']:
            if candidate in cols:
                actual_cell_col = candidate
                break

        actual_gene_col = None
        if gene_col in cols:
            actual_gene_col = gene_col
        elif 'feature_name' in cols:
            actual_gene_col = 'feature_name'
        elif 'gene' in cols:
            actual_gene_col = 'gene'

        if not actual_cell_col or not actual_gene_col:
            logger.error(f"Missing required columns (cell/gene) in {csv_path}. Found: {list(cols)}")
            return None

        # Filter unassigned transcripts
        if pd.api.types.is_numeric_dtype(df[actual_cell_col]):
            df_assigned = df[df[actual_cell_col] > 0].copy()
        else:
            mask = (
                (df[actual_cell_col] != '0')
                & (df[actual_cell_col].astype(str).str.lower() != 'unassigned')
                & (df[actual_cell_col].astype(str).str.lower() != 'noise')
            )
            df_assigned = df[mask].copy()

        # Aggregate via crosstab (faster than loop + value_counts)
        counts = pd.crosstab(df_assigned[actual_cell_col], df_assigned[actual_gene_col])
        adata = ad.AnnData(counts)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        # Compute cell centroids from transcript coordinates for spatial plotting
        for xc, yc in [('x_location', 'y_location'), ('x_global_px', 'y_global_px'), ('x', 'y')]:
            if xc in df_assigned.columns and yc in df_assigned.columns:
                centroids = df_assigned.groupby(actual_cell_col)[[xc, yc]].mean()
                centroids.columns = ['x_centroid', 'y_centroid']
                adata.obs = adata.obs.join(centroids)
                logger.info(f"  Computed centroids from ({xc}, {yc}) for "
                            f"{adata.obs['x_centroid'].notna().sum()} cells.")
                break

        # Store raw spots for proportion_of_assigned_reads + Rand Index
        # Preserve coordinate columns if available for spatial matching
        keep_cols = [actual_cell_col, actual_gene_col]
        for coord_candidate in ['x_location', 'y_location', 'x_global_px', 'y_global_px', 'x', 'y']:
            if coord_candidate in cols:
                keep_cols.append(coord_candidate)
        adata.uns['spots'] = df[list(dict.fromkeys(keep_cols))].copy()

        return adata

    except Exception as e:
        logger.error(f"Error aggregating csv {csv_path}: {e}")
        return None


# --- 6-4. Annotation transfer ---

def annotate_by_majority_voting(adata_target, adata_reference,
                                cluster_key='leiden', ref_key='celltype',
                                n_neighbors=15):
    """Transfer cell type labels from reference to target via cluster-level majority voting.

    Matching notebook/paper methodology:
      1. For each cluster in adata_target, find cells that overlap with adata_reference
         (by obs_names / cell_id) and assign the most frequent reference label.
      2. If no direct cell overlap is found, fall back to kNN in PCA space.

    Adds 'celltype_majority' (per-cell) and 'celltype_cluster' (per-cluster consensus)."""

    if ref_key not in adata_reference.obs.columns:
        logger.error(f"Reference AnnData is missing the '{ref_key}' column.")
        return adata_target

    # --- Strategy 1: Cluster-level crosstab majority vote (notebook method) ---
    # Match cells between target and reference by obs_names or cell_id column
    ref_labels = adata_reference.obs[ref_key]

    # Try matching by obs_names first
    shared_cells = adata_target.obs_names.intersection(adata_reference.obs_names)

    # If obs_names don't match, try cell_id column
    if len(shared_cells) < 10 and 'cell_id' in adata_target.obs.columns and 'cell_id' in adata_reference.obs.columns:
        target_cids = adata_target.obs['cell_id'].astype(str)
        ref_cids = adata_reference.obs['cell_id'].astype(str)
        shared_cids = set(target_cids) & set(ref_cids)
        if len(shared_cids) > len(shared_cells):
            # Build mapping via cell_id
            ref_cid_to_label = dict(zip(ref_cids, adata_reference.obs[ref_key]))
            adata_target.obs['_ref_label'] = target_cids.map(ref_cid_to_label)
            shared_cells = adata_target.obs_names[adata_target.obs['_ref_label'].notna()]
            logger.info(f"Matched {len(shared_cells)} cells via cell_id column.")

    overlap_ratio = len(shared_cells) / max(adata_target.n_obs, 1)
    logger.info(f"Annotation transfer: {len(shared_cells)} overlapping cells "
                f"({overlap_ratio:.1%} of target).")

    if overlap_ratio > 0.01:
        # Cluster-level majority vote (matching notebook 5_1 / paper methodology)
        logger.info("Using cluster-level crosstab majority vote (notebook/paper method).")

        # Get reference labels for overlapping cells
        if '_ref_label' in adata_target.obs.columns:
            overlap_labels = adata_target.obs.loc[shared_cells, '_ref_label']
        else:
            overlap_labels = ref_labels.loc[shared_cells]

        # Build crosstab: cluster x cell_type counts
        overlap_clusters = adata_target.obs.loc[shared_cells, cluster_key]
        crosstab = pd.crosstab(overlap_clusters, overlap_labels)

        # For each cluster, assign the most frequent reference cell type
        cluster_map = {}
        for cluster_id in crosstab.index:
            winner = crosstab.loc[cluster_id].idxmax()
            cluster_map[cluster_id] = winner

        # Assign to ALL cells in each cluster (including non-overlapping ones like Baysor cells)
        adata_target.obs['celltype_cluster'] = (
            adata_target.obs[cluster_key].map(cluster_map).astype(str)
        )
        adata_target.obs['celltype_majority'] = adata_target.obs['celltype_cluster']

        # Confidence = proportion of overlapping cells in that cluster matching the winner
        def _cluster_confidence(row):
            cid = row[cluster_key]
            if cid in crosstab.index:
                total = crosstab.loc[cid].sum()
                return crosstab.loc[cid].max() / total if total > 0 else 0.0
            return 0.0

        adata_target.obs['celltype_confidence'] = adata_target.obs.apply(_cluster_confidence, axis=1)

        # Clean up temp column
        if '_ref_label' in adata_target.obs.columns:
            del adata_target.obs['_ref_label']

        logger.info(
            f"Annotation transfer complete (crosstab method). "
            f"{len(cluster_map)} clusters mapped to cell types."
        )
    else:
        # --- Strategy 2: kNN in gene-space-aligned PCA (when no cell overlap) ---
        logger.info("Low cell overlap — falling back to kNN in aligned PCA space.")
        from sklearn.neighbors import NearestNeighbors

        if '_ref_label' in adata_target.obs.columns:
            del adata_target.obs['_ref_label']

        shared_genes = sorted(adata_target.var_names.intersection(adata_reference.var_names))
        if len(shared_genes) < 10:
            logger.error(f"Too few shared genes ({len(shared_genes)}) for kNN transfer.")
            return adata_target

        logger.info(f"kNN transfer: {len(shared_genes)} shared genes, {n_neighbors} neighbours "
                     "(gene-space aligned PCA).")

        # Build reference PCA on shared gene space
        ref_sub = adata_reference[:, shared_genes].copy()

        # Filter out cells with missing cell type labels (NaN mixed with strings
        # causes TypeError in np.unique during kNN majority voting)
        valid_mask = ref_sub.obs[ref_key].notna()
        if hasattr(ref_sub.obs[ref_key], 'str'):
            valid_mask = valid_mask & (ref_sub.obs[ref_key].astype(str) != 'nan')
        n_dropped = (~valid_mask).sum()
        if n_dropped > 0:
            ref_sub = ref_sub[valid_mask].copy()
            logger.info(f"  Dropped {n_dropped} reference cells with missing '{ref_key}'.")

        if 'raw' in ref_sub.layers:
            ref_sub.X = ref_sub.layers['raw'].copy()
        if issparse(ref_sub.X):
            ref_sub.X = np.asarray(ref_sub.X.todense())
        sc.pp.normalize_total(ref_sub, target_sum=1e4)
        sc.pp.log1p(ref_sub)
        sc.pp.scale(ref_sub, max_value=10)
        sc.pp.pca(ref_sub)

        # Project target into reference PCA space using reference loadings
        tgt_sub = adata_target[:, shared_genes].copy()
        if 'raw' in tgt_sub.layers:
            tgt_sub.X = tgt_sub.layers['raw'].copy()
        if issparse(tgt_sub.X):
            tgt_sub.X = np.asarray(tgt_sub.X.todense())
        sc.pp.normalize_total(tgt_sub, target_sum=1e4)
        sc.pp.log1p(tgt_sub)

        # Scale target with reference statistics and project onto reference PCs
        tgt_X = np.array(tgt_sub.X)
        ref_mean = ref_sub.var['mean'].values
        ref_std = ref_sub.var['std'].values.copy()
        ref_std[ref_std == 0] = 1.0  # avoid division by zero
        tgt_X = np.clip((tgt_X - ref_mean) / ref_std, -10, 10)
        tgt_pca = tgt_X @ ref_sub.varm['PCs']

        # kNN in aligned PCA space
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
        nn.fit(ref_sub.obsm['X_pca'])
        distances, indices = nn.kneighbors(tgt_pca)

        ref_labels_arr = ref_sub.obs[ref_key].values
        cell_labels = []
        confidence_scores = []
        for idx_row in indices:
            neighbour_labels = ref_labels_arr[idx_row]
            values, counts = np.unique(neighbour_labels, return_counts=True)
            winner_idx = np.argmax(counts)
            cell_labels.append(values[winner_idx])
            confidence_scores.append(counts[winner_idx] / len(neighbour_labels))

        adata_target.obs['celltype_majority'] = cell_labels
        adata_target.obs['celltype_confidence'] = confidence_scores

        # Per-cluster consensus
        cluster_labels = []
        for cluster_id in adata_target.obs[cluster_key].unique():
            mask = adata_target.obs[cluster_key] == cluster_id
            cluster_celltypes = adata_target.obs.loc[mask, 'celltype_majority']
            values, counts = np.unique(cluster_celltypes, return_counts=True)
            winner = values[np.argmax(counts)]
            cluster_labels.append((cluster_id, winner))

        cluster_map = dict(cluster_labels)
        adata_target.obs['celltype_cluster'] = (
            adata_target.obs[cluster_key].map(cluster_map).astype(str)
        )

        logger.info("Annotation transfer complete (kNN in aligned PCA space).")

    return adata_target


# --- 6-3. Preprocessing ---

def preprocess_benchmark(adata, config):
    """Preprocess concatenated benchmark AnnData (normalize, PCA, Leiden).
    Parameters from config["benchmark"]["preprocessing"], defaults match xb/preprocessing.py."""
    pp = config.get("benchmark", {}).get("preprocessing", {})

    target_sum = pp.get("target_sum", 100)
    min_counts = pp.get("min_counts", 40)
    min_genes = pp.get("min_genes", 15)
    n_neighbors = pp.get("n_neighbors", 15)
    n_pcs = pp.get("n_pcs", 0)          # 0 = use ALL PCs
    umap_min_dist = pp.get("umap_min_dist", 0.1)
    resolution = pp.get("resolution", 1.0)
    scale = pp.get("scale", False)

    logger.info(
        f"Preprocessing params: target_sum={target_sum}, min_counts={min_counts}, "
        f"min_genes={min_genes}, n_neighbors={n_neighbors}, n_pcs={n_pcs}, "
        f"umap_min_dist={umap_min_dist}, resolution={resolution}, scale={scale}"
    )

    adata.layers['raw'] = adata.X.copy()

    # QC filtering
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata.raw = adata

    # Re-store raw layer after filtering (cell count may have changed)
    adata.layers['raw'] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # HVG selection (matching notebook's preprocess_adata defaults)
    hvg = pp.get("hvg", True)
    if hvg:
        sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
        logger.info(f"  HVG: {adata.var['highly_variable'].sum()} / {adata.shape[1]} genes selected")

    if scale:
        sc.pp.scale(adata)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=umap_min_dist)

    # Notebook labels this Louvain but actually uses Leiden
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

    # Auto-reduce resolution if too many clusters (>40 is usually unreadable)
    n_clusters = adata.obs['leiden'].nunique()
    if n_clusters > 40:
        reduced_res = resolution * 0.5
        logger.info(f"  Leiden produced {n_clusters} clusters (resolution={resolution}). "
                     f"Re-clustering with reduced resolution={reduced_res:.2f}")
        sc.tl.leiden(adata, resolution=reduced_res, key_added='leiden')
        n_clusters_new = adata.obs['leiden'].nunique()
        logger.info(f"  Re-clustered: {n_clusters} → {n_clusters_new} clusters")

    return adata


# --- 6-6. Visualizations ---

# --- 6-6-pre. QC plots after preprocessing ---
def _save_qc_histograms(adata, output_dir, min_counts=40, min_genes=15):
    """C3-1: Cell counts/genes histograms with filtering threshold lines."""
    try:
        raw_X = adata.layers['raw'] if 'raw' in adata.layers else adata.X
        if issparse(raw_X):
            total_counts = np.asarray(raw_X.sum(axis=1)).flatten()
            n_genes = np.asarray((raw_X > 0).sum(axis=1)).flatten()
        else:
            total_counts = np.sum(raw_X, axis=1)
            n_genes = np.sum(raw_X > 0, axis=1)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes[0].hist(total_counts, bins=50, edgecolor='black', alpha=0.7)
        axes[0].axvline(x=min_counts, color='red', linestyle='--', linewidth=1.5,
                        label=f'min_counts={min_counts}')
        axes[0].set_title('Total Counts per Cell')
        axes[0].set_xlabel('Counts')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()

        axes[1].hist(n_genes, bins=50, edgecolor='black', alpha=0.7, color='orange')
        axes[1].axvline(x=min_genes, color='red', linestyle='--', linewidth=1.5,
                        label=f'min_genes={min_genes}')
        axes[1].set_title('Genes Detected per Cell')
        axes[1].set_xlabel('Number of Genes')
        axes[1].set_ylabel('Frequency')
        axes[1].legend()

        fig.tight_layout()
        out_path = os.path.join(output_dir, "qc_counts_genes_histogram.png")
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved QC histograms to {out_path}")
    except Exception as e:
        logger.warning(f"QC histograms failed: {e}")
        plt.close('all')


def _save_hvg_plot(adata, output_dir):
    """C3-2: HVG selection plot."""
    try:
        if 'highly_variable' not in adata.var.columns:
            logger.info("No HVG column found; skipping HVG plot.")
            return
        sc.pl.highly_variable_genes(adata, show=False)
        out_path = os.path.join(output_dir, "hvg_selection_plot.png")
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved HVG plot to {out_path}")
    except Exception as e:
        logger.warning(f"HVG plot failed: {e}")
        plt.close('all')


def _save_pca_scree(adata, output_dir):
    """C3-3: PCA scree plot (variance explained)."""
    try:
        if 'pca' not in adata.uns:
            logger.info("No PCA results found; skipping scree plot.")
            return
        variance_ratio = adata.uns['pca']['variance_ratio']
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(range(1, len(variance_ratio) + 1), np.cumsum(variance_ratio), 'o-', markersize=3)
        ax.set_xlabel('PC')
        ax.set_ylabel('Cumulative Variance Explained')
        ax.set_title('PCA Scree Plot')
        ax.axhline(y=0.9, color='red', linestyle='--', alpha=0.5, label='90%')
        ax.legend()
        fig.tight_layout()
        out_path = os.path.join(output_dir, "pca_scree_plot.png")
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved PCA scree plot to {out_path}")
    except Exception as e:
        logger.warning(f"PCA scree plot failed: {e}")
        plt.close('all')


# --- 6-6a. UMAP per method ---
def _save_umap(adata, output_dir):
    """UMAP coloured by segmentation method, Leiden cluster, and cell type (if available)."""
    n_panels = 2
    ct_key = 'celltype_majority'
    has_ct = ct_key in adata.obs.columns
    if has_ct:
        n_panels = 3

    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 8))
    n_cells = adata.n_obs
    pt_size = max(3, min(20, 80000 / n_cells))  # auto-scale point size

    sc.pl.umap(adata, color='segmentation', ax=axes[0], show=False,
               title='Segmentation method', frameon=False, size=pt_size, alpha=0.7)
    # Suppress legend for Leiden if too many clusters (>20 makes legend unreadable)
    n_leiden = adata.obs['leiden'].nunique()
    sc.pl.umap(adata, color='leiden', ax=axes[1], show=False,
               title=f'Leiden clusters (n={n_leiden})', frameon=False,
               size=pt_size, alpha=0.7,
               legend_loc='on data' if n_leiden > 15 else 'right margin',
               legend_fontsize=5 if n_leiden > 15 else 8,
               legend_fontoutline=1)
    if has_ct:
        sc.pl.umap(adata, color=ct_key, ax=axes[2], show=False,
                    title='Cell type', frameon=False)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "umap_benchmark.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved UMAP plot to {out_path}")


# --- 6-6-dapi. DAPI / morphology background helpers ---

def _load_dapi_for_benchmark(config, max_size=2000):
    """Load and downsample DAPI/morphology image for benchmark visualizations.

    Returns (dapi_ds, step) or (None, 1) if unavailable.
    """
    import tifffile as tf

    dapi_path = config.get('dapi_image_path')
    if not dapi_path or not os.path.exists(dapi_path):
        return None, 1

    try:
        dapi = tf.imread(dapi_path)
        if dapi.ndim == 3:
            dapi = dapi[0]
        step = max(1, max(dapi.shape) // max_size)
        dapi_ds = dapi[::step, ::step]
        del dapi
        logger.info(f"Loaded DAPI for benchmark viz: {dapi_path} "
                    f"(downsampled {step}x -> {dapi_ds.shape})")
        return dapi_ds, step
    except Exception as e:
        logger.warning(f"Could not load DAPI image: {e}")
        return None, 1


def _load_masks_for_benchmark(config, dapi_step, crop_px=None, target_shape=None):
    """Load Step 3 mask TIFs and downsample to match DAPI.

    Parameters
    ----------
    crop_px : tuple or None
        (r0, r1, c0, c1) in full-res pixel coords.  When given the mask is
        first sliced to the crop region, *then* downsampled.
    target_shape : tuple or None
        (H, W) of the downsampled DAPI.  Masks are trimmed to match (avoids
        off-by-one from different crop+downsample paths).

    Returns (cellpose_masks_ds, nuclei_masks_ds) — either may be None.
    """
    import tifffile as tf

    def _trim(arr):
        if target_shape is not None:
            return arr[:target_shape[0], :target_shape[1]]
        return arr

    cellpose_masks_ds, nuclei_masks_ds = None, None

    cellpose_path = config.get('cellpose_mask_path')
    if cellpose_path and os.path.exists(cellpose_path):
        try:
            masks = tf.imread(cellpose_path)
            if crop_px is not None:
                r0, r1, c0, c1 = crop_px
                masks = masks[r0:r1, c0:c1]
            cellpose_masks_ds = _trim(masks[::dapi_step, ::dapi_step])
            del masks
            logger.info(f"Loaded cellpose masks: {cellpose_path} "
                        f"(downsampled {dapi_step}x -> {cellpose_masks_ds.shape})")
        except Exception as e:
            logger.warning(f"Could not load cellpose masks: {e}")

    nuclei_path = config.get('nuclei_mask_path')
    if nuclei_path and os.path.exists(nuclei_path):
        try:
            masks = tf.imread(nuclei_path)
            if crop_px is not None:
                r0, r1, c0, c1 = crop_px
                masks = masks[r0:r1, c0:c1]
            nuclei_masks_ds = _trim(masks[::dapi_step, ::dapi_step])
            del masks
            logger.info(f"Loaded nuclei masks: {nuclei_path} "
                        f"(downsampled {dapi_step}x -> {nuclei_masks_ds.shape})")
        except Exception as e:
            logger.warning(f"Could not load nuclei masks: {e}")

    return cellpose_masks_ds, nuclei_masks_ds


def _load_baysor_polygons(config):
    """Load Baysor polygon boundaries from GeoJSON.

    Returns list of (N,2) arrays (µm coords), or empty list.
    """
    step6_out = config.get("output_dir", "xenium-output/benchmark")
    poly_path = os.path.join(step6_out, "baysor_run", "crop",
                             "segmentation_polygons_2d.json")
    if not os.path.exists(poly_path):
        return []
    try:
        with open(poly_path) as f:
            geo = json.load(f)
        polys = []
        for feat in geo.get("features", []):
            coords = feat["geometry"]["coordinates"]
            # GeoJSON Polygon: coords[0] is the outer ring
            ring = np.array(coords[0])  # (N, 2) in µm
            polys.append(ring)
        logger.info(f"Loaded {len(polys)} Baysor polygons from {poly_path}")
        return polys
    except Exception as e:
        logger.warning(f"Could not load Baysor polygons: {e}")
        return []


def _crop_coords_um_to_px(crop_coords_um, um_per_pixel_inv):
    """Convert (x_min, x_max, y_min, y_max) µm → (r0, r1, c0, c1) full-res pixels.

    Note: Xenium coords are (x=col, y=row), so x→col, y→row.
    """
    x_min, x_max, y_min, y_max = crop_coords_um
    c0 = int(x_min * um_per_pixel_inv)
    c1 = int(x_max * um_per_pixel_inv)
    r0 = int(y_min * um_per_pixel_inv)
    r1 = int(y_max * um_per_pixel_inv)
    return (r0, r1, c0, c1)


def _draw_baysor_polygons(ax, polys_um, crop_um=None, px_per_um_ds=1.0):
    """Draw Baysor polygon outlines on a matplotlib axes.

    Parameters
    ----------
    polys_um : list of (N,2) arrays in µm (x, y)
    crop_um  : (x_min, x_max, y_min, y_max) or None
    px_per_um_ds : pixels-per-µm in the displayed (downsampled) image
    """
    from matplotlib.collections import LineCollection

    lines = []
    x_off = crop_um[0] if crop_um else 0
    y_off = crop_um[2] if crop_um else 0
    for ring in polys_um:
        xs = (ring[:, 0] - x_off) * px_per_um_ds
        ys = (ring[:, 1] - y_off) * px_per_um_ds
        pts = np.column_stack([xs, ys])
        lines.append(pts)
    lc = LineCollection(lines, colors='#FFD700', linewidths=0.4, alpha=0.8)
    ax.add_collection(lc)


def _save_dapi_mask_boundary_overlay(config, dapi_ds, dapi_step, output_dir,
                                     crop_coords_um=None):
    """DAPI background + mask boundary contours + Baysor polygons overlay.

    Panels: DAPI | cellpose expanded (cyan) | nuclei (orange) | Baysor (gold).
    When ``crop_coords_um`` is given, DAPI and masks are cropped to that region.
    """
    from skimage.segmentation import find_boundaries

    um_per_pixel_inv = config.get("resegmentation", {}).get("um_per_pixel_inv", 4.70588)

    # --- crop DAPI if needed ---
    crop_px = None
    if crop_coords_um is not None:
        crop_px = _crop_coords_um_to_px(crop_coords_um, um_per_pixel_inv)
        r0, r1, c0, c1 = crop_px
        dapi_ds = dapi_ds[r0 // dapi_step : r1 // dapi_step,
                          c0 // dapi_step : c1 // dapi_step]

    cellpose_masks_ds, nuclei_masks_ds = _load_masks_for_benchmark(
        config, dapi_step, crop_px=crop_px, target_shape=dapi_ds.shape[:2])

    baysor_polys = _load_baysor_polygons(config) if crop_coords_um else []

    if cellpose_masks_ds is None and nuclei_masks_ds is None and not baysor_polys:
        logger.info("No masks or Baysor polygons available; skipping mask boundary overlay.")
        return

    vmax = np.percentile(dapi_ds, 99.5)
    dapi_rgb = np.stack([dapi_ds / vmax] * 3, axis=-1).clip(0, 1)

    n_panels = (1 + (cellpose_masks_ds is not None)
                + (nuclei_masks_ds is not None)
                + (len(baysor_polys) > 0))
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 7))
    if n_panels == 1:
        axes = [axes]
    idx = 0

    # Panel 1: Raw DAPI
    axes[idx].imshow(dapi_ds, cmap='gray', vmax=vmax)
    axes[idx].set_title('DAPI')
    axes[idx].axis('off')
    idx += 1

    # Panel 2: DAPI + cellpose expanded boundaries (cyan)
    if cellpose_masks_ds is not None:
        boundaries = find_boundaries(cellpose_masks_ds, mode='thick')
        overlay = dapi_rgb.copy()
        overlay[boundaries, 0] = 0
        overlay[boundaries, 1] = 1
        overlay[boundaries, 2] = 1
        axes[idx].imshow(overlay)
        n_cells = len(np.unique(cellpose_masks_ds)) - (1 if 0 in cellpose_masks_ds else 0)
        axes[idx].set_title(f'DAPI + Cellpose ({n_cells:,} cells)')
        axes[idx].axis('off')
        idx += 1

    # Panel 3: DAPI + nuclei boundaries (orange)
    if nuclei_masks_ds is not None:
        nuc_boundaries = find_boundaries(nuclei_masks_ds, mode='thick')
        nuc_overlay = dapi_rgb.copy()
        nuc_overlay[nuc_boundaries, 0] = 1
        nuc_overlay[nuc_boundaries, 1] = 0.3
        nuc_overlay[nuc_boundaries, 2] = 0
        axes[idx].imshow(nuc_overlay)
        n_nuc = len(np.unique(nuclei_masks_ds)) - (1 if 0 in nuclei_masks_ds else 0)
        axes[idx].set_title(f'DAPI + Nuclei ({n_nuc:,})')
        axes[idx].axis('off')
        idx += 1

    # Panel 4: DAPI + Baysor polygon boundaries (gold)
    if baysor_polys:
        axes[idx].imshow(dapi_ds, cmap='gray', vmax=vmax)
        px_per_um_ds = um_per_pixel_inv / dapi_step
        _draw_baysor_polygons(axes[idx], baysor_polys,
                              crop_um=crop_coords_um, px_per_um_ds=px_per_um_ds)
        axes[idx].set_xlim(0, dapi_ds.shape[1])
        axes[idx].set_ylim(dapi_ds.shape[0], 0)
        axes[idx].set_title(f'DAPI + Baysor ({len(baysor_polys):,} cells)')
        axes[idx].axis('off')

    region_label = "crop region" if crop_coords_um else "full tissue"
    fig.suptitle(f'Step 6 — Segmentation mask overlay on DAPI ({region_label})', fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "dapi_mask_boundary_overlay.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DAPI mask boundary overlay: {out_path}")

    del cellpose_masks_ds, nuclei_masks_ds


def _save_dapi_mask_zoomed_roi(config, dapi_ds, dapi_step, output_dir,
                               crop_coords_um=None, roi_fraction=0.15):
    """Zoomed ROI: DAPI + mask boundaries + Baysor polygons + colored blend.

    4 panels: DAPI | dual mask boundaries (cyan+orange) | Baysor (gold) | blended.
    When ``crop_coords_um`` is given, operates on that region.
    """
    from skimage.segmentation import find_boundaries
    from skimage.color import label2rgb

    um_per_pixel_inv = config.get("resegmentation", {}).get("um_per_pixel_inv", 4.70588)

    # --- crop DAPI if needed ---
    crop_px = None
    if crop_coords_um is not None:
        crop_px = _crop_coords_um_to_px(crop_coords_um, um_per_pixel_inv)
        r0, r1, c0, c1 = crop_px
        dapi_ds = dapi_ds[r0 // dapi_step : r1 // dapi_step,
                          c0 // dapi_step : c1 // dapi_step]

    cellpose_masks_ds, nuclei_masks_ds = _load_masks_for_benchmark(
        config, dapi_step, crop_px=crop_px, target_shape=dapi_ds.shape[:2])

    baysor_polys = _load_baysor_polygons(config) if crop_coords_um else []

    if cellpose_masks_ds is None and not baysor_polys:
        logger.info("No cellpose mask or Baysor polygons; skipping mask zoomed ROI.")
        return

    h, w = dapi_ds.shape[:2]

    # Find center of dense region
    if cellpose_masks_ds is not None:
        ref_mask = nuclei_masks_ds if nuclei_masks_ds is not None else cellpose_masks_ds
        mask_present = np.argwhere(ref_mask > 0)
        if len(mask_present) == 0:
            logger.warning("No non-zero pixels in masks; skipping zoomed ROI.")
            return
        cy, cx = mask_present[len(mask_present) // 2]
    else:
        cy, cx = h // 2, w // 2

    # ROI crop
    rh, rw = int(h * roi_fraction), int(w * roi_fraction)
    r0, r1 = max(0, cy - rh // 2), min(h, cy + rh // 2)
    c0, c1 = max(0, cx - rw // 2), min(w, cx + rw // 2)

    dapi_crop = dapi_ds[r0:r1, c0:c1]

    vmax = np.percentile(dapi_ds, 99.5)
    dapi_rgb = np.stack([dapi_crop / vmax] * 3, axis=-1).clip(0, 1)

    has_masks = cellpose_masks_ds is not None
    has_baysor = len(baysor_polys) > 0
    n_panels = 1 + has_masks + has_baysor + has_masks  # dapi + boundaries + baysor + blend
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 7))
    if n_panels == 1:
        axes = [axes]
    idx = 0

    # Panel 1: DAPI zoomed
    axes[idx].imshow(dapi_crop, cmap='gray', vmax=vmax)
    axes[idx].set_title('DAPI (zoomed ROI)')
    axes[idx].axis('off')
    idx += 1

    # Panel 2: DAPI + dual mask boundaries
    if has_masks:
        mask_crop = cellpose_masks_ds[r0:r1, c0:c1]
        boundaries = find_boundaries(mask_crop, mode='thick')
        overlay = dapi_rgb.copy()
        overlay[boundaries, 0] = 0
        overlay[boundaries, 1] = 1
        overlay[boundaries, 2] = 1
        if nuclei_masks_ds is not None:
            nuc_crop = nuclei_masks_ds[r0:r1, c0:c1]
            nuc_bd = find_boundaries(nuc_crop, mode='thick')
            overlay[nuc_bd, 0] = 1
            overlay[nuc_bd, 1] = 0.5
            overlay[nuc_bd, 2] = 0
        axes[idx].imshow(overlay)
        axes[idx].set_title('Cellpose+Nuclei boundaries')
        axes[idx].axis('off')
        idx += 1

    # Panel 3: DAPI + Baysor polygons
    if has_baysor:
        axes[idx].imshow(dapi_crop, cmap='gray', vmax=vmax)
        px_per_um_ds = um_per_pixel_inv / dapi_step
        # Offset polygons to match the ROI crop within the (already-cropped) dapi_ds
        roi_x_off = c0 / px_per_um_ds  # µm offset of ROI within crop
        roi_y_off = r0 / px_per_um_ds
        crop_x_um = crop_coords_um[0] if crop_coords_um else 0
        crop_y_um = crop_coords_um[2] if crop_coords_um else 0
        total_x_off = crop_x_um + roi_x_off
        total_y_off = crop_y_um + roi_y_off
        shifted_crop_um = (total_x_off,
                           total_x_off + (c1 - c0) / px_per_um_ds,
                           total_y_off,
                           total_y_off + (r1 - r0) / px_per_um_ds)
        _draw_baysor_polygons(axes[idx], baysor_polys,
                              crop_um=shifted_crop_um, px_per_um_ds=px_per_um_ds)
        axes[idx].set_xlim(0, dapi_crop.shape[1])
        axes[idx].set_ylim(dapi_crop.shape[0], 0)
        axes[idx].set_title(f'Baysor ({len(baysor_polys):,} cells)')
        axes[idx].axis('off')
        idx += 1

    # Panel 4: DAPI + colored masks (blended)
    if has_masks:
        mask_crop = cellpose_masks_ds[r0:r1, c0:c1]
        mask_rgb = label2rgb(mask_crop, bg_label=0)
        blended = 0.5 * dapi_rgb + 0.5 * mask_rgb
        axes[idx].imshow(blended.clip(0, 1))
        axes[idx].set_title('Cellpose colored (blended)')
        axes[idx].axis('off')

    region_label = "crop region" if crop_coords_um else "full tissue"
    fig.suptitle(f'Step 6 — Zoomed ROI segmentation on DAPI ({region_label})', fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "dapi_mask_zoomed_roi.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DAPI mask zoomed ROI: {out_path}")

    del cellpose_masks_ds, nuclei_masks_ds


def _save_dapi_segmentation_comparison(adata, dapi_ds, output_dir):
    """DAPI background + per-method cell centroids overlay.

    Matching notebook pattern: plt.imshow(dapi) + plt.scatter(centroids).
    Creates one panel per segmentation method, all with shared DAPI background.
    """
    x_col, y_col = None, None
    for xc, yc in [('x_centroid', 'y_centroid'),
                    ('x_location', 'y_location'),
                    ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        logger.warning("No spatial coordinates for DAPI overlay.")
        return

    methods = adata.obs['segmentation'].unique()
    n_methods = len(methods)

    # Color cycle for methods
    method_colors = {'nuclei': '#00BFFF', 'cellpose': '#FF6347',
                     'expansion': '#32CD32', 'baysor': '#FFD700'}
    dapi_h, dapi_w = dapi_ds.shape
    vmax = np.percentile(dapi_ds, 99.5)

    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 7))
    if n_methods == 1:
        axes = [axes]

    for ax, method in zip(axes, methods):
        ax.imshow(dapi_ds, cmap='gray', vmax=vmax, extent=[0, dapi_w, dapi_h, 0])
        sub = adata[adata.obs['segmentation'] == method]
        x_vals = sub.obs[x_col].astype(float).values
        y_vals = sub.obs[y_col].astype(float).values

        # Scale coordinates to match downsampled DAPI
        x_range = max(x_vals.max() - x_vals.min(), 1)
        if x_range > 20000:  # pixel coords
            scale = dapi_w / (x_vals.max() - x_vals.min() + 1) if x_range > 0 else 1
            x_plot = (x_vals - x_vals.min()) * scale
            y_plot = (y_vals - y_vals.min()) * scale
        else:  # micron coords — use um_per_pixel_inv to get to pixels, then scale
            x_plot = x_vals * (dapi_w / (x_range + 1)) if x_range > 0 else x_vals
            y_plot = y_vals * (dapi_h / (y_vals.max() - y_vals.min() + 1)) if len(y_vals) > 0 else y_vals

        c = method_colors.get(method, '#FFFFFF')
        ax.scatter(x_plot, y_plot, s=0.2, alpha=0.4, c=c, rasterized=True)
        ax.set_title(f"{method} (n={len(sub):,})", fontsize=12, color=c, fontweight='bold')
        ax.axis('off')

    fig.suptitle("Segmentation methods on DAPI morphology", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "dapi_segmentation_comparison.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DAPI segmentation comparison: {out_path}")


def _save_dapi_zoomed_roi(adata, dapi_ds, output_dir, roi_fraction=0.2):
    """Zoomed ROI on DAPI comparing segmentation methods side by side.

    Picks a dense region and shows each method's cells overlaid on DAPI.
    """
    x_col, y_col = None, None
    for xc, yc in [('x_centroid', 'y_centroid'),
                    ('x_location', 'y_location'),
                    ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        return

    methods = adata.obs['segmentation'].unique()
    n_methods = len(methods)
    method_colors = {'nuclei': '#00BFFF', 'cellpose': '#FF6347',
                     'expansion': '#32CD32', 'baysor': '#FFD700'}

    dapi_h, dapi_w = dapi_ds.shape
    vmax = np.percentile(dapi_ds, 99.5)

    # Pick ROI center: find densest region from first method
    first_sub = adata[adata.obs['segmentation'] == methods[0]]
    x_vals = first_sub.obs[x_col].astype(float).values
    y_vals = first_sub.obs[y_col].astype(float).values

    # Scale to DAPI coords
    x_range = max(x_vals.max() - x_vals.min(), 1)
    if x_range > 20000:
        x_scale = dapi_w / (x_range + 1)
        y_scale = dapi_h / (max(y_vals.max() - y_vals.min(), 1) + 1)
        x_offset, y_offset = x_vals.min(), y_vals.min()
    else:
        x_scale = dapi_w / (x_range + 1)
        y_scale = dapi_h / (max(y_vals.max() - y_vals.min(), 1) + 1)
        x_offset, y_offset = x_vals.min(), y_vals.min()

    # Find median cell position as ROI center
    cx = int(np.median((x_vals - x_offset) * x_scale))
    cy = int(np.median((y_vals - y_offset) * y_scale))
    rw = int(dapi_w * roi_fraction)
    rh = int(dapi_h * roi_fraction)
    r0 = max(0, cy - rh // 2)
    r1 = min(dapi_h, cy + rh // 2)
    c0 = max(0, cx - rw // 2)
    c1 = min(dapi_w, cx + rw // 2)

    dapi_crop = dapi_ds[r0:r1, c0:c1]

    fig, axes = plt.subplots(1, n_methods + 1, figsize=(6 * (n_methods + 1), 6))

    # First panel: DAPI only
    axes[0].imshow(dapi_crop, cmap='gray', vmax=vmax)
    axes[0].set_title('DAPI (zoomed ROI)')
    axes[0].axis('off')

    # Per-method panels
    for i, method in enumerate(methods):
        ax = axes[i + 1]
        ax.imshow(dapi_crop, cmap='gray', vmax=vmax,
                  extent=[c0, c1, r1, r0])

        sub = adata[adata.obs['segmentation'] == method]
        xv = (sub.obs[x_col].astype(float).values - x_offset) * x_scale
        yv = (sub.obs[y_col].astype(float).values - y_offset) * y_scale

        # Filter to ROI
        in_roi = (xv >= c0) & (xv <= c1) & (yv >= r0) & (yv <= r1)
        c_color = method_colors.get(method, '#FFFFFF')
        ax.scatter(xv[in_roi], yv[in_roi], s=2, alpha=0.7, c=c_color, rasterized=True)
        ax.set_title(f"{method} ({in_roi.sum():,} cells)", color=c_color, fontweight='bold')
        ax.set_xlim(c0, c1)
        ax.set_ylim(r1, r0)
        ax.axis('off')

    fig.suptitle("Zoomed ROI — Segmentation on DAPI", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "dapi_zoomed_roi_comparison.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DAPI zoomed ROI comparison: {out_path}")


def _save_boundary_overlay(adata, config, dapi_ds, output_dir):
    """Overlay nucleus boundaries (from Xenium boundary files) on DAPI.

    Loads nucleus_boundaries.parquet and plots boundary polygons for a sample
    of cells, coloured by segmentation method assignment.
    """
    bd_path = config.get('nucleus_boundaries_path')
    if not bd_path or not os.path.exists(bd_path):
        logger.info("No nucleus boundary file available; skipping boundary overlay.")
        return

    try:
        if bd_path.endswith('.parquet'):
            bd_df = pd.read_parquet(bd_path)
        elif bd_path.endswith('.csv.gz'):
            bd_df = pd.read_csv(bd_path)
        else:
            bd_df = pd.read_csv(bd_path)
    except Exception as e:
        logger.warning(f"Failed to load boundaries: {e}")
        return

    # Expect columns: cell_id, vertex_x, vertex_y
    if 'vertex_x' not in bd_df.columns or 'vertex_y' not in bd_df.columns:
        logger.warning(f"Boundary file missing vertex_x/vertex_y columns: {list(bd_df.columns[:5])}")
        return

    dapi_h, dapi_w = dapi_ds.shape
    vmax = np.percentile(dapi_ds, 99.5)

    # Scale boundaries to downsampled DAPI
    vx = bd_df['vertex_x'].values
    vy = bd_df['vertex_y'].values
    vx_range = vx.max() - vx.min()
    vy_range = vy.max() - vy.min()

    # Boundaries are typically in pixel space
    if vx_range > 20000:
        bd_scale_x = dapi_w / (vx_range + 1)
        bd_scale_y = dapi_h / (vy_range + 1)
        vx_offset, vy_offset = vx.min(), vy.min()
    else:
        bd_scale_x = dapi_w / (vx_range + 1)
        bd_scale_y = dapi_h / (vy_range + 1)
        vx_offset, vy_offset = vx.min(), vy.min()

    # Sample cells to avoid overplotting (max 500 boundary polygons)
    unique_cells = bd_df['cell_id'].unique()
    if len(unique_cells) > 500:
        np.random.seed(42)
        sample_cells = np.random.choice(unique_cells, 500, replace=False)
        bd_sample = bd_df[bd_df['cell_id'].isin(sample_cells)]
    else:
        bd_sample = bd_df

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(dapi_ds, cmap='gray', vmax=vmax)

    for cell_id in bd_sample['cell_id'].unique():
        cell_bd = bd_sample[bd_sample['cell_id'] == cell_id]
        bx = (cell_bd['vertex_x'].values - vx_offset) * bd_scale_x
        by = (cell_bd['vertex_y'].values - vy_offset) * bd_scale_y
        ax.plot(bx, by, c='cyan', linewidth=0.3, alpha=0.6)

    ax.set_title(f'DAPI + nucleus boundaries ({len(bd_sample["cell_id"].unique())} cells)')
    ax.axis('off')
    fig.tight_layout()
    out_path = os.path.join(output_dir, "dapi_nucleus_boundaries.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved DAPI + nucleus boundaries: {out_path}")


# --- 6-6b. Spatial scatter ---
def _save_spatial_map(adata, output_dir):
    """Spatial scatter plot coloured by segmentation method.
    Requires x/y coordinates in adata.obs; skipped if not available."""
    x_col, y_col = None, None
    for xc, yc in [('x_centroid', 'y_centroid'),
                    ('x_location', 'y_location'),
                    ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        logger.warning(
            "Spatial coordinates not found in adata.obs. "
            "Skipping spatial segmentation map."
        )
        return

    methods = adata.obs['segmentation'].unique()
    n_methods = len(methods)
    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 6))
    if n_methods == 1:
        axes = [axes]

    # Determine common coordinate range for consistent axis limits
    all_x = adata.obs[x_col].dropna()
    all_y = adata.obs[y_col].dropna()
    x_lim = (all_x.quantile(0.001), all_x.quantile(0.999))
    y_lim = (all_y.quantile(0.001), all_y.quantile(0.999))

    for ax, method in zip(axes, methods):
        sub = adata[adata.obs['segmentation'] == method]
        ax.scatter(
            sub.obs[x_col].values, sub.obs[y_col].values,
            s=0.3, alpha=0.5, rasterized=True
        )
        ax.set_title(f"{method} (n={len(sub):,})")
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_aspect('equal')
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.invert_yaxis()

    fig.suptitle("Spatial map by segmentation method", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "spatial_segmentation_map.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved spatial segmentation map to {out_path}")


# --- 6-6c. Cell type barplot ---
def _save_celltype_barplot(adata, output_dir, ct_key='celltype_majority'):
    """Bar plot of cell type proportions across segmentation methods."""
    if ct_key not in adata.obs.columns:
        logger.warning(
            f"'{ct_key}' not in adata.obs -- skipping cell type frequency barplot."
        )
        return

    ct_counts = (
        adata.obs
        .groupby(['segmentation', ct_key])
        .size()
        .reset_index(name='count')
    )
    totals = ct_counts.groupby('segmentation')['count'].transform('sum')
    ct_counts['proportion'] = ct_counts['count'] / totals

    methods = ct_counts['segmentation'].unique()
    celltypes = ct_counts[ct_key].unique()
    n_ct = len(celltypes)

    x = np.arange(n_ct)
    width = 0.8 / len(methods)

    fig, ax = plt.subplots(figsize=(max(10, n_ct * 0.8), 6))
    for i, method in enumerate(methods):
        sub = ct_counts[ct_counts['segmentation'] == method].set_index(ct_key)
        heights = sub.reindex(celltypes, fill_value=0)['proportion']
        ax.bar(x + i * width, heights.values, width, label=method)

    ax.set_xticks(x + width * (len(methods) - 1) / 2)
    ax.set_xticklabels(celltypes, rotation=45, ha='right')
    ax.set_ylabel("Proportion")
    ax.set_title("Cell type frequency by segmentation method")
    ax.legend(title="Method")
    fig.tight_layout()
    out_path = os.path.join(output_dir, "celltype_frequency_barplot.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved cell type frequency barplot to {out_path}")


# --- 6-6d. Counts violin ---
def _save_counts_violin(adata, output_dir):
    """Violin plot of total counts per cell grouped by segmentation method."""
    if 'raw' not in adata.layers:
        logger.warning("adata.layers['raw'] not found -- skipping counts violin plot.")
        return

    raw_X = adata.layers['raw']
    if issparse(raw_X):
        total_counts = np.asarray(raw_X.sum(axis=1)).flatten()
    else:
        total_counts = np.sum(raw_X, axis=1)

    plot_df = pd.DataFrame({
        'total_counts': total_counts,
        'segmentation': adata.obs['segmentation'].values
    })

    methods = plot_df['segmentation'].unique()
    fig, ax = plt.subplots(figsize=(max(6, len(methods) * 2), 6))

    parts = ax.violinplot(
        [plot_df.loc[plot_df['segmentation'] == m, 'total_counts'].values
         for m in methods],
        positions=range(len(methods)),
        showmeans=True, showmedians=True
    )

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods)
    ax.set_ylabel("Total counts per cell")
    ax.set_title("Counts distribution by segmentation method")
    fig.tight_layout()
    out_path = os.path.join(output_dir, "counts_violin_by_method.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved counts violin plot to {out_path}")


# --- 6-6f. Reads Assigned vs NMP scatter (Extended Data Fig. 5e) ---
def _save_reads_vs_nmp_scatter(metrics_df, output_dir):
    """Scatter plot of proportion of assigned reads vs negative marker purity.

    Reproduces Extended Data Figure 5e from Salas et al. (Nature Methods, 2025):
    each segmentation method is a point, x = proportion of assigned reads,
    y = negative marker purity (NMP).  Both NMP-cells and NMP-reads variants
    are shown as separate subplots.
    """
    import seaborn as sns

    method_markers = {
        'nuclei': ('o', '#00BFFF', 'Xenium nuclei'),
        'cellpose': ('^', '#FF6347', 'Cellpose'),
        'expansion': ('s', '#32CD32', 'Expansion'),
        'baysor': ('D', '#FFD700', 'Baysor'),
    }

    # Parse benchmark_metrics into per-method rows
    if metrics_df is None or metrics_df.empty:
        logger.warning("No benchmark metrics available for reads-vs-NMP scatter.")
        return

    row = metrics_df.iloc[0]

    # Detect methods available in the CSV (columns: <method>_assigned_prop, <method>_nmp_cells, <method>_nmp_reads)
    methods_found = set()
    for col in metrics_df.columns:
        for suffix in ['_nmp_cells', '_nmp_reads', '_assigned_prop', '_median_reads']:
            if col.endswith(suffix):
                m = col[: -len(suffix)]
                methods_found.add(m)

    if not methods_found:
        logger.warning("No per-method NMP/assigned metrics found in benchmark_metrics.")
        return

    # Build per-method data
    records = []
    for m in sorted(methods_found):
        d = {'method': m}
        for k in ['nmp_cells', 'nmp_reads', 'assigned_prop', 'median_reads',
                   'n_cells', 'median_genes']:
            col = f'{m}_{k}'
            if col in row.index:
                val = row[col]
                try:
                    d[k] = float(val)
                except (ValueError, TypeError):
                    d[k] = np.nan
        records.append(d)

    df = pd.DataFrame(records)

    # Need at least assigned_prop OR median_reads + one NMP variant
    has_assigned_prop = 'assigned_prop' in df.columns and df['assigned_prop'].notna().any()
    has_nmp_cells = 'nmp_cells' in df.columns and df['nmp_cells'].notna().any()
    has_nmp_reads = 'nmp_reads' in df.columns and df['nmp_reads'].notna().any()

    if not (has_nmp_cells or has_nmp_reads):
        logger.info("No NMP data available; skipping reads-vs-NMP scatter.")
        return

    # Determine how many subplots we need
    panels = []
    if has_nmp_cells:
        panels.append(('nmp_cells', 'Negative Marker Purity (cells)'))
    if has_nmp_reads:
        panels.append(('nmp_reads', 'Negative Marker Purity (reads)'))

    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))
    if n_panels == 1:
        axes = [axes]

    for ax, (nmp_col, ylabel) in zip(axes, panels):
        # X-axis: prefer assigned_prop, fall back to median_reads
        if has_assigned_prop:
            x_col = 'assigned_prop'
            xlabel = 'Proportion of assigned reads'
        else:
            x_col = 'median_reads'
            xlabel = 'Median reads per cell'

        for _, mrow in df.iterrows():
            m = mrow['method']
            x_val = mrow.get(x_col, np.nan)
            y_val = mrow.get(nmp_col, np.nan)

            if pd.isna(x_val) or pd.isna(y_val):
                continue

            marker, color, label = method_markers.get(
                m, ('o', '#888888', m.capitalize()))

            ax.scatter(x_val, y_val, marker=marker, c=color, s=200,
                       edgecolors='black', linewidths=0.8, zorder=5,
                       label=label)

            # Annotate with method name + n_cells
            n_cells = mrow.get('n_cells', np.nan)
            annot = f"{label}"
            if not pd.isna(n_cells):
                annot += f"\n({int(n_cells):,} cells)"
            ax.annotate(annot, (x_val, y_val),
                        textcoords='offset points', xytext=(12, -5),
                        fontsize=8, fontstyle='italic',
                        arrowprops=dict(arrowstyle='->', color='gray',
                                        lw=0.5, connectionstyle='arc3,rad=0.2'))

        # Reference lines for paper expectations
        ax.axhline(y=0.8, color='red', linestyle='--', linewidth=1, alpha=0.5,
                    label='NMP = 0.8 (paper threshold)')

        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(f'Extended Data Fig. 5e: {ylabel}', fontsize=11, fontweight='bold')

        # Adjust y-axis to show relevant range
        y_vals = df[nmp_col].dropna()
        if len(y_vals) > 0:
            y_min = max(0, y_vals.min() - 0.05)
            y_max = min(1.0, y_vals.max() + 0.02)
            ax.set_ylim(y_min, y_max)

        ax.legend(loc='lower left', fontsize=8, framealpha=0.8)
        ax.grid(True, alpha=0.3, linestyle=':')

    fig.suptitle('Segmentation Benchmark: Reads Assigned vs. Negative Marker Purity\n'
                 '(Salas et al. Extended Data Fig. 5e style)',
                 fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "reads_vs_nmp_scatter.png")
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved reads-vs-NMP scatter (Ext. Data Fig. 5e): {out_path}")


# --- 6-6e. Spatial scatter coloured by cell type (notebook cell 27: map_of_clusters) ---
def _save_spatial_celltype_map(adata, output_dir, ct_key='celltype_majority'):
    """Spatial scatter per segmentation method, coloured by cell type (matching notebook)."""
    if ct_key not in adata.obs.columns:
        logger.info(f"'{ct_key}' not in adata.obs -- skipping spatial cell type map.")
        return

    x_col, y_col = None, None
    for xc, yc in [('x_centroid', 'y_centroid'),
                    ('x_location', 'y_location'),
                    ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        logger.warning("Spatial coordinates not found -- skipping spatial cell type map.")
        return

    methods = adata.obs['segmentation'].unique()
    celltypes = adata.obs[ct_key].unique()
    n_methods = len(methods)

    # Build a consistent color map across all panels
    cmap = plt.cm.get_cmap('tab20', len(celltypes))
    ct_colors = {ct: cmap(i) for i, ct in enumerate(sorted(celltypes))}

    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 6))
    if n_methods == 1:
        axes = [axes]

    for ax, method in zip(axes, methods):
        sub = adata[adata.obs['segmentation'] == method]
        for ct in sorted(celltypes):
            mask = sub.obs[ct_key] == ct
            if mask.sum() == 0:
                continue
            ax.scatter(
                sub.obs.loc[mask, x_col].values,
                sub.obs.loc[mask, y_col].values,
                s=0.3, alpha=0.5, color=ct_colors[ct],
                label=ct, rasterized=True
            )
        ax.set_title(f"{method} ({ct_key})")
        ax.set_aspect('equal')
        ax.invert_yaxis()

    # Single legend for all panels
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', fontsize=6, markerscale=5,
               bbox_to_anchor=(1.12, 0.5))
    fig.suptitle("Spatial map by cell type per segmentation", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "spatial_celltype_map.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved spatial cell type map to {out_path}")


# --- 6-5b. Rand Index between segmentation methods ---

def _compute_rand_index(adata_list, method_names, output_dir):
    """Compute pairwise Rand Index between segmentation methods using transcript-level assignments.

    Matches transcripts across methods by (x, y) coordinates, then calls
    metrics.rand_idx() on the aligned assignment matrix.
    """
    if len(adata_list) < 2:
        logger.info("Rand Index requires >=2 methods; skipping.")
        return

    # Build per-method transcript→cell DataFrames keyed by (x, y)
    dfs = {}
    skipped = []
    for ad, name in zip(adata_list, method_names):
        if 'spots' not in ad.uns:
            logger.warning(f"Rand Index: 'spots' missing for '{name}'; skipping this method. "
                           f"To include it, ensure adata.uns['spots'] contains transcript-level assignments.")
            skipped.append(name)
            continue
        spots = ad.uns['spots']
        # Find coordinate columns
        x_col = y_col = cell_col = None
        for xc, yc in [('x_location', 'y_location'), ('x_global_px', 'y_global_px'), ('x', 'y')]:
            if xc in spots.columns and yc in spots.columns:
                x_col, y_col = xc, yc
                break
        if x_col is None:
            logger.warning(f"Rand Index: no coordinate columns in spots for '{name}'; skipping this method.")
            continue
        cell_col = [c for c in spots.columns if c not in [x_col, y_col]][0]
        df = spots[[x_col, y_col, cell_col]].copy()
        df.columns = ['x', 'y', name]
        df['x'] = df['x'].round(3)
        df['y'] = df['y'].round(3)
        dfs[name] = df

    if len(dfs) < 2:
        logger.info(f"Rand Index requires >=2 methods with transcript-level data; skipping. "
                     f"Methods skipped (no spots): {skipped}")
        return

    # Merge all methods on coordinates
    names = list(dfs.keys())
    merged = dfs[names[0]]
    for n in names[1:]:
        merged = merged.merge(dfs[n], on=['x', 'y'], how='inner')

    if len(merged) < 100:
        logger.warning(f"Rand Index: only {len(merged)} matched transcripts; skipping.")
        return

    logger.info(f"Rand Index: {len(merged)} matched transcripts across {len(names)} methods.")
    assignments = merged[names].fillna('0').astype(str)  # unify types for ARI comparison

    try:
        rand_matrix = metrics.rand_idx(assignments)
        rand_matrix.index = names
        rand_matrix.columns = names

        # Save CSV
        csv_path = os.path.join(output_dir, "rand_index_matrix.csv")
        rand_matrix.to_csv(csv_path)
        logger.info(f"Saved Rand Index matrix to {csv_path}")

        # Heatmap
        fig, ax = plt.subplots(figsize=(max(6, len(names) * 1.5), max(5, len(names) * 1.2)))
        import seaborn as sns
        sns.heatmap(rand_matrix.astype(float), annot=True, fmt='.3f',
                    cmap='YlOrRd', vmin=0, vmax=1, square=True, ax=ax)
        ax.set_title("Adjusted Rand Index between Segmentation Methods")
        fig.tight_layout()
        fig_path = os.path.join(output_dir, "rand_index_heatmap.png")
        fig.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved Rand Index heatmap to {fig_path}")
    except Exception as e:
        logger.warning(f"Rand Index computation failed: {e}")
        plt.close('all')


# --- 6-1b. Spatial crop subsetting for fair Baysor comparison ---

def _load_baysor_crop_coords(config):
    """Load Baysor crop coordinates (in µm) from metadata JSON.

    Returns (x_min, x_max, y_min, y_max) in microns, or None if unavailable.
    """
    step6_out = config.get("output_dir", "xenium-output/benchmark")
    meta_path = os.path.join(step6_out, "baysor_run", "crop", "crop_metadata.json")
    if not os.path.exists(meta_path):
        return None
    try:
        with open(meta_path) as f:
            meta = json.load(f)
        coords = meta.get("coords")
        if coords and len(coords) == 4:
            logger.info(f"Loaded Baysor crop coords (µm): x=[{coords[0]:.0f}, {coords[1]:.0f}], "
                        f"y=[{coords[2]:.0f}, {coords[3]:.0f}]")
            return tuple(coords)
    except Exception as e:
        logger.warning(f"Could not load Baysor crop metadata: {e}")
    return None


def _subset_to_crop_region(adata, crop_coords_um, um_per_pixel_inv, method_name):
    """Subset an AnnData to cells within the Baysor crop region.

    Handles both micron- and pixel-space centroids by auto-detecting the scale.

    Parameters
    ----------
    adata : AnnData
    crop_coords_um : tuple (x_min, x_max, y_min, y_max) in microns
    um_per_pixel_inv : float, pixels per micron
    method_name : str, for logging

    Returns
    -------
    AnnData subset, or original if no spatial columns found.
    """
    x_col = y_col = None
    for xc, yc in [('x_centroid', 'y_centroid'), ('x_location', 'y_location'), ('x', 'y')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        logger.warning(f"[{method_name}] No spatial coordinates found; cannot subset to crop region.")
        return adata

    x_vals = adata.obs[x_col].astype(float)
    y_vals = adata.obs[y_col].astype(float)
    x_min_um, x_max_um, y_min_um, y_max_um = crop_coords_um

    # Auto-detect coordinate space: if max > 20000, coords are in pixels
    coord_range = max(x_vals.max() - x_vals.min(), y_vals.max() - y_vals.min())
    if coord_range > 20000:
        # Coordinates are in pixels — convert crop bounds to pixels
        x_min = x_min_um * um_per_pixel_inv
        x_max = x_max_um * um_per_pixel_inv
        y_min = y_min_um * um_per_pixel_inv
        y_max = y_max_um * um_per_pixel_inv
        coord_unit = "px"
    else:
        # Coordinates are in microns
        x_min, x_max = x_min_um, x_max_um
        y_min, y_max = y_min_um, y_max_um
        coord_unit = "µm"

    mask = (x_vals >= x_min) & (x_vals <= x_max) & (y_vals >= y_min) & (y_vals <= y_max)
    n_before = len(adata)
    n_after = mask.sum()

    adata_sub = adata[mask].copy()
    logger.info(f"[{method_name}] Subset to Baysor crop region ({coord_unit}): "
                f"{n_before:,} -> {n_after:,} cells "
                f"(x=[{x_min:.0f},{x_max:.0f}], y=[{y_min:.0f},{y_max:.0f}])")
    return adata_sub



# --- Benchmark Pipeline (reusable) ---

def _run_benchmark_pipeline(adata_list, config, output_dir, benchmark_label="main",
                            crop_coords_um=None):
    """Run the complete benchmark pipeline on a set of method AnnDatas.

    Flow: Rand Index -> coord scaling -> concat -> preprocess -> annotate -> metrics -> viz.
    All outputs are saved to ``output_dir``.

    Parameters
    ----------
    adata_list : list[AnnData]
        Each must have ``obs['segmentation']`` set to the method name.
    config : dict
        Full pipeline configuration.
    output_dir : str
        Directory for this benchmark run's outputs.
    benchmark_label : str
        Human-readable label for log messages (e.g. "full_tissue", "crop_region").
    crop_coords_um : tuple or None
        (x_min, x_max, y_min, y_max) in µm for crop region mask overlays.
    """
    os.makedirs(output_dir, exist_ok=True)

    method_names_in = [ad.obs['segmentation'].iloc[0] for ad in adata_list]
    logger.info(f"=== Benchmark [{benchmark_label}]: methods={method_names_in} -> {output_dir} ===")

    combined_output = os.path.join(output_dir, "benchmark_combined.h5ad")
    if os.path.exists(combined_output):
        logger.info(f"  [{benchmark_label}] [CACHE HIT] Loading: {combined_output}")
        adata = sc.read_h5ad(combined_output)
    else:
        # --- Rand Index (before concat, needs per-method spots) ---
        try:
            method_names = [ad.obs['segmentation'].iloc[0] for ad in adata_list]
            _compute_rand_index(adata_list, method_names, output_dir)
        except Exception as e:
            logger.warning(f"Rand Index step failed: {e}")

        # --- Coordinate scaling (microns -> pixels for spatial plots) ---
        um_per_pixel_inv = config.get("resegmentation", {}).get("um_per_pixel_inv", 4.70588)
        for ad_item in adata_list:
            for coord_col in ['x_centroid', 'y_centroid']:
                if coord_col in ad_item.obs.columns:
                    vals = ad_item.obs[coord_col].astype(float)
                    coord_range = vals.max() - vals.min()
                    # Heuristic: if range < 20000, coordinates are likely in um -> scale to pixels
                    if coord_range > 0 and coord_range < 20000:
                        ad_item.obs[coord_col] = vals * um_per_pixel_inv
                        logger.info(f"  Scaled {coord_col} um->px (x{um_per_pixel_inv}) "
                                    f"for method '{ad_item.obs['segmentation'].iloc[0]}'")

        # --- Concatenate & preprocess ---
        for ad_item in adata_list:
            method = ad_item.obs['segmentation'].iloc[0]
            logger.info(f"  [{method}] {ad_item.shape[0]} cells, {ad_item.shape[1]} genes "
                        f"(sample vars: {list(ad_item.var_names[:3])})")

        logger.info(f"[{benchmark_label}] Concatenating {len(adata_list)} datasets...")
        adata = sc.concat(adata_list)

        if adata.shape[0] == 0 or adata.shape[1] == 0:
            logger.error(
                f"[{benchmark_label}] Concatenation produced empty AnnData {adata.shape}. "
                "Likely gene-name mismatch between datasets (inner join). "
                "Check var_names across inputs."
            )
            return

        logger.info(f"  [{benchmark_label}] Concatenated shape: {adata.shape}")
        logger.info("Preprocessing and Clustering...")
        adata = preprocess_benchmark(adata, config)

        # --- QC visualizations (C3-1, C3-2, C3-3) ---
        pp = config.get("benchmark", {}).get("preprocessing", {})
        _save_qc_histograms(adata, output_dir,
                            min_counts=pp.get("min_counts", 40),
                            min_genes=pp.get("min_genes", 15))
        _save_hvg_plot(adata, output_dir)
        _save_pca_scree(adata, output_dir)

        # --- Annotation transfer ---
        ref_adata_path = config.get("benchmark", {}).get("reference_adata")
        if ref_adata_path and os.path.exists(ref_adata_path):
            logger.info(f"Loading reference AnnData for annotation transfer: {ref_adata_path}")
            adata_ref = sc.read(ref_adata_path)

            if 'X_pca' not in adata_ref.obsm:
                logger.info("Computing PCA on reference AnnData for annotation transfer...")
                sc.pp.pca(adata_ref)

            ref_key = config.get("benchmark", {}).get("ref_celltype_key", "celltype")
            # Fallback: if ref_key not in reference, try sc_reference.celltype_key
            if ref_key not in adata_ref.obs.columns:
                alt_key = config.get("sc_reference", {}).get("celltype_key", "subclass_label")
                if alt_key in adata_ref.obs.columns:
                    logger.info(f"ref_celltype_key '{ref_key}' not in reference; "
                                f"using '{alt_key}' instead.")
                    ref_key = alt_key
                else:
                    logger.warning(
                        f"Reference lacks both '{ref_key}' and '{alt_key}'. "
                        f"Available columns: {list(adata_ref.obs.columns[:10])}")
            pp_conf = config.get("benchmark", {}).get("preprocessing", {})
            n_neighbors = pp_conf.get("n_neighbors", 15)

            adata = annotate_by_majority_voting(
                adata, adata_ref,
                cluster_key='leiden', ref_key=ref_key,
                n_neighbors=n_neighbors
            )
        else:
            logger.info(
                "No reference AnnData provided (benchmark.reference_adata). "
                "Skipping annotation transfer."
            )

        adata.write_h5ad(combined_output)
        logger.info(f"Saved combined benchmark AnnData to {combined_output}")

    # --- Benchmark metrics ---
    # Convert layers['raw'] to dense ndarray -- sparse .sum() returns np.matrix
    # which causes shape alignment errors in metrics (matmul vs element-wise)
    if 'raw' in adata.layers and issparse(adata.layers['raw']):
        adata.layers['raw'] = np.asarray(adata.layers['raw'].todense())

    logger.info(f"[{benchmark_label}] Calculating Benchmark Metrics...")
    results = {}

    for seg_method in adata.obs['segmentation'].unique():
        subset = adata[adata.obs['segmentation'] == seg_method]
        n_cells = subset.shape[0]
        results[f'{seg_method}_n_cells'] = n_cells

        if 'raw' in subset.layers:
            results[f'{seg_method}_median_reads'] = metrics.median_reads_cells(subset)
            results[f'{seg_method}_median_genes'] = metrics.median_genes_cells(subset)
            results[f'{seg_method}_p5_reads'] = metrics.percentile_5th_reads_cells(subset)
            results[f'{seg_method}_p5_genes'] = metrics.percentile_5th_genes_cells(subset)
        else:
            logger.warning(
                f"layers['raw'] missing for '{seg_method}' -- "
                "skipping median_reads / median_genes metrics."
            )

        if 'spots' in subset.uns:
            results[f'{seg_method}_assigned_prop'] = metrics.proportion_of_assigned_reads(subset)

        # Clustering quality metrics (silhouette, Calinski-Harabasz, Davies-Bouldin)
        if 'X_pca' in subset.obsm and 'leiden' in subset.obs.columns:
            labels = subset.obs['leiden'].astype('category').cat.codes.values
            n_labels = len(np.unique(labels))
            if n_labels >= 2 and n_labels < subset.shape[0]:
                try:
                    X_pca = subset.obsm['X_pca']
                    results[f'{seg_method}_silhouette'] = silhouette_score(X_pca, labels)
                    results[f'{seg_method}_calinski_harabasz'] = calinski_harabasz_score(X_pca, labels)
                    results[f'{seg_method}_davies_bouldin'] = davies_bouldin_score(X_pca, labels)
                    logger.info(
                        f"  [{seg_method}] Silhouette={results[f'{seg_method}_silhouette']:.3f}, "
                        f"CH={results[f'{seg_method}_calinski_harabasz']:.1f}, "
                        f"DB={results[f'{seg_method}_davies_bouldin']:.3f}"
                    )
                except Exception as e:
                    logger.warning(f"  [{seg_method}] Clustering quality metrics failed: {e}")

    # --- NMP metrics (requires scRNA-seq reference with raw layer) ---
    sc_ref_path = config.get("comparison", {}).get("sc_reference_path")
    if not sc_ref_path:
        sc_ref_path = config.get("benchmark", {}).get("sc_reference_path")

    ct_key_for_nmp = 'celltype_majority' if 'celltype_majority' in adata.obs.columns else None
    if sc_ref_path and os.path.exists(sc_ref_path) and ct_key_for_nmp:
        logger.info(f"Computing NMP metrics using scRNA-seq reference: {sc_ref_path}")
        try:
            adata_sc = sc.read_h5ad(sc_ref_path)
            if 'raw' not in adata_sc.layers:
                adata_sc.layers['raw'] = adata_sc.X.copy()

            # Ensure reference has the cell type key expected by NMP metrics
            if ct_key_for_nmp not in adata_sc.obs.columns:
                for fallback_key in [
                    config.get("benchmark", {}).get("ref_celltype_key"),
                    config.get("sc_reference", {}).get("celltype_key"),
                    "subclass_label", "celltype", "cell_type",
                ]:
                    if fallback_key and fallback_key in adata_sc.obs.columns:
                        adata_sc.obs[ct_key_for_nmp] = adata_sc.obs[fallback_key]
                        logger.info(f"Mapped reference '{fallback_key}' -> "
                                    f"'{ct_key_for_nmp}' for NMP metrics.")
                        break
                else:
                    logger.warning(
                        f"Reference lacks cell type column for NMP "
                        f"(tried: celltype_majority, ref_celltype_key, "
                        f"sc_reference.celltype_key, subclass_label, celltype, cell_type). "
                        f"Available: {list(adata_sc.obs.columns[:10])}")
                    ct_key_for_nmp = None

            if not ct_key_for_nmp:
                logger.info("Skipping NMP: could not resolve cell type key in reference.")
            else:
                for seg_method in adata.obs['segmentation'].unique():
                    subset = adata[adata.obs['segmentation'] == seg_method].copy()
                    if 'raw' not in subset.layers:
                        continue
                    try:
                        nmp_cells = metrics.negative_marker_purity_cells(
                            subset, adata_sc, key=ct_key_for_nmp, pipeline_output=True)
                        results[f'{seg_method}_nmp_cells'] = nmp_cells
                        logger.info(f"  [{seg_method}] NMP (cells) = {nmp_cells}")
                    except Exception as e:
                        logger.warning(f"  [{seg_method}] NMP cells failed: {e}")

                    try:
                        nmp_reads = metrics.negative_marker_purity_reads(
                            subset, adata_sc, key=ct_key_for_nmp, pipeline_output=True)
                        results[f'{seg_method}_nmp_reads'] = nmp_reads
                        logger.info(f"  [{seg_method}] NMP (reads) = {nmp_reads}")
                    except Exception as e:
                        logger.warning(f"  [{seg_method}] NMP reads failed: {e}")
        except Exception as e:
            logger.warning(f"NMP computation failed: {e}")
    elif ct_key_for_nmp is None:
        logger.info("Skipping NMP: no cell type annotation available (run annotation transfer first).")
    else:
        logger.info("Skipping NMP: no scRNA-seq reference provided (comparison.sc_reference_path).")

    metrics_df = pd.DataFrame([results])
    metrics_file = os.path.join(output_dir, "benchmark_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    logger.info(f"Saved metrics to {metrics_file}")

    # --- Extended Data Fig. 5e: Reads Assigned vs NMP scatter ---
    try:
        _save_reads_vs_nmp_scatter(metrics_df, output_dir)
    except Exception as e:
        logger.warning(f"Reads-vs-NMP scatter (Ext. Data Fig. 5e) failed: {e}")
        plt.close('all')

    # --- Marker gene ranking per segmentation method ---
    try:
        sc.tl.rank_genes_groups(adata, groupby='segmentation', method='wilcoxon',
                                use_raw=True, key_added='rank_genes_segmentation')
        deg_path = os.path.join(output_dir, "marker_genes_by_segmentation.png")
        sc.pl.rank_genes_groups(adata, key='rank_genes_segmentation', show=False,
                                save=False, n_genes=10)
        plt.savefig(deg_path, dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved marker gene ranking plot: {deg_path}")
    except Exception as e:
        logger.warning(f"Marker gene ranking failed: {e}")
        plt.close('all')

    if 'celltype_majority' in adata.obs.columns:
        try:
            sc.tl.rank_genes_groups(adata, groupby='celltype_majority', method='wilcoxon',
                                    use_raw=True, key_added='rank_genes_celltype')
            deg_ct_path = os.path.join(output_dir, "marker_genes_by_celltype.png")
            n_ct = adata.obs['celltype_majority'].nunique()
            n_cols = min(4, n_ct)
            n_rows = (n_ct + n_cols - 1) // n_cols
            fig_w = max(12, n_cols * 4)
            fig_h = max(8, n_rows * 3.5)
            sc.settings.set_figure_params(figsize=(fig_w, fig_h))
            sc.pl.rank_genes_groups(adata, key='rank_genes_celltype', show=False,
                                    save=False, n_genes=5)
            plt.savefig(deg_ct_path, dpi=150, bbox_inches='tight')
            plt.close()
            sc.settings.set_figure_params(figsize=(4, 4))  # reset
            logger.info(f"Saved celltype marker gene plot: {deg_ct_path}")
        except Exception as e:
            logger.warning(f"Celltype marker gene ranking failed: {e}")
            plt.close('all')

    # C3-4: DEG dotplot per cluster
    try:
        if 'rank_genes_segmentation' in adata.uns:
            sc.pl.rank_genes_groups_dotplot(
                adata, key='rank_genes_segmentation', n_genes=5,
                show=False, save=False)
            dotplot_path = os.path.join(output_dir, "deg_dotplot_by_segmentation.png")
            plt.savefig(dotplot_path, dpi=150, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved DEG dotplot: {dotplot_path}")
    except Exception as e:
        logger.warning(f"DEG dotplot failed: {e}")
        plt.close('all')

    # --- Visualizations ---
    logger.info(f"[{benchmark_label}] Generating visualizations...")
    _save_umap(adata, output_dir)
    _save_spatial_map(adata, output_dir)
    _save_celltype_barplot(adata, output_dir, ct_key='celltype_majority')
    _save_counts_violin(adata, output_dir)
    _save_spatial_celltype_map(adata, output_dir)

    # --- DAPI / morphology overlays ---
    dapi_ds, _dapi_step = _load_dapi_for_benchmark(config)
    if dapi_ds is not None:
        logger.info(f"[{benchmark_label}] Generating DAPI overlay visualizations...")
        try:
            _save_dapi_segmentation_comparison(adata, dapi_ds, output_dir)
        except Exception as e:
            logger.warning(f"DAPI segmentation comparison failed: {e}")
            plt.close('all')
        try:
            _save_dapi_zoomed_roi(adata, dapi_ds, output_dir)
        except Exception as e:
            logger.warning(f"DAPI zoomed ROI failed: {e}")
            plt.close('all')
        try:
            _save_boundary_overlay(adata, config, dapi_ds, output_dir)
        except Exception as e:
            logger.warning(f"Boundary overlay failed: {e}")
            plt.close('all')
        # Mask boundary overlays (if Step 3 masks available)
        try:
            _save_dapi_mask_boundary_overlay(config, dapi_ds, _dapi_step, output_dir,
                                             crop_coords_um=crop_coords_um)
        except Exception as e:
            logger.warning(f"DAPI mask boundary overlay failed: {e}")
            plt.close('all')
        try:
            _save_dapi_mask_zoomed_roi(config, dapi_ds, _dapi_step, output_dir,
                                       crop_coords_um=crop_coords_um)
        except Exception as e:
            logger.warning(f"DAPI mask zoomed ROI failed: {e}")
            plt.close('all')
        del dapi_ds
    else:
        logger.info(f"[{benchmark_label}] No DAPI image available; skipping morphology overlays.")

    logger.info(f"Benchmark [{benchmark_label}] completed.")


# --- Entry Point ---

def run_step6(config):
    """Step 6 orchestrator: Baysor -> load methods -> benchmark(s).

    Supports two modes based on Baysor configuration:
      - **Dual mode** (Baysor crop + subset_other_methods): produces two benchmark sets:
        1. Full tissue benchmark (nuclei + cellpose + expansion, no Baysor)
        2. Crop comparison benchmark (all methods subset to Baysor crop region)
      - **Single mode** (Baysor disabled/full, or subset_other_methods=false):
        produces one benchmark with all available methods.
    """
    logger.info("Starting Step 6: Segmentation Benchmark")

    output_dir = config.get("output_dir", "xenium-output/benchmark")
    os.makedirs(output_dir, exist_ok=True)

    # --- 6-1. Run Baysor (optional) ---
    baysor_result_csv = None
    if config.get("baysor", {}).get("enabled", False):
        xenium_input_dir = config.get("input_dir")
        if xenium_input_dir:
            baysor_result_csv = run_baysor(config, xenium_input_dir)

    # --- 6-2. Load segmentation results ---
    logger.info("--- Loading Segmentation Results ---")

    input_nuclei = config.get("benchmark", {}).get("input_nuclei")
    input_advanced = config.get("benchmark", {}).get("input_advanced")
    input_expansion = config.get("benchmark", {}).get("input_expansion")

    adata_list = []

    if input_nuclei and os.path.exists(input_nuclei):
        logger.info(f"Loading Nuclei Data: {input_nuclei}")
        nuclei_adata = sc.read(input_nuclei)
        nuclei_adata.obs['segmentation'] = 'nuclei'
        adata_list.append(nuclei_adata)

    if input_advanced and os.path.exists(input_advanced) and str(input_advanced).endswith('.h5ad'):
        logger.info(f"Loading Cellpose Data: {input_advanced}")
        adv_adata = sc.read(input_advanced)
        adv_adata.obs['segmentation'] = 'cellpose'
        adata_list.append(adv_adata)

    if input_expansion and os.path.exists(input_expansion):
        logger.info(f"Loading Expansion Data (Step 5): {input_expansion}")
        exp_adata = load_transcripts_as_adata(
            input_expansion, current_method_name='expansion'
        )
        if exp_adata is not None:
            exp_adata.obs['segmentation'] = 'expansion'
            adata_list.append(exp_adata)

    if baysor_result_csv:
        baysor_adata = load_transcripts_as_adata(
            baysor_result_csv, cell_col='cell', gene_col='gene',
            current_method_name='baysor'
        )
        if baysor_adata is not None:
            baysor_adata.obs['segmentation'] = 'baysor'
            adata_list.append(baysor_adata)

    if not adata_list:
        logger.error("No valid input datasets found for benchmarking.")
        return

    # --- 6-3. Determine benchmark mode ---
    baysor_crop_enabled = config.get("baysor", {}).get("crop", {}).get("enabled", False)
    crop_coords_um = _load_baysor_crop_coords(config)
    has_baysor = any(ad.obs['segmentation'].iloc[0] == 'baysor' for ad in adata_list)
    subset_other = config.get("baysor", {}).get("crop", {}).get("subset_other_methods", True)

    if baysor_crop_enabled and has_baysor and crop_coords_um and subset_other:
        # ── DUAL BENCHMARK MODE ──
        # 1) Full tissue: nuclei + cellpose + expansion (no crop-mode Baysor)
        # 2) Crop comparison: all methods subset to Baysor crop region
        logger.info("=" * 60)
        logger.info("DUAL BENCHMARK MODE: Full tissue + Crop comparison")
        logger.info("=" * 60)

        # A) Full tissue benchmark (all methods except crop-mode Baysor)
        full_list = [ad.copy() for ad in adata_list
                     if ad.obs['segmentation'].iloc[0] != 'baysor']
        if full_list:
            logger.info("--- [1/2] Full tissue benchmark (without Baysor) ---")
            _run_benchmark_pipeline(full_list, config, output_dir, "full_tissue")
        else:
            logger.warning("No non-Baysor methods available for full tissue benchmark.")

        # B) Crop comparison benchmark (all methods in the Baysor crop region)
        crop_dir = os.path.join(output_dir, "crop_comparison")
        um_per_pixel_inv = config.get("resegmentation", {}).get("um_per_pixel_inv", 4.70588)
        x_min_um, x_max_um, y_min_um, y_max_um = crop_coords_um
        logger.info(f"--- [2/2] Crop comparison benchmark "
                    f"({x_max_um - x_min_um:.0f} x {y_max_um - y_min_um:.0f} um) ---")

        cropped_list = []
        for ad_item in adata_list:
            method = ad_item.obs['segmentation'].iloc[0]
            if method == 'baysor':
                cropped_list.append(ad_item)  # Baysor already cropped
            else:
                ad_cropped = _subset_to_crop_region(
                    ad_item, crop_coords_um, um_per_pixel_inv, method)
                if ad_cropped.shape[0] > 0:
                    cropped_list.append(ad_cropped)
                else:
                    logger.warning(f"[{method}] No cells in crop region — excluded.")

        if cropped_list:
            # Validate crop quality
            min_cells_for_benchmark = 100
            for ad_item in cropped_list:
                method = ad_item.obs['segmentation'].iloc[0]
                n_cells = ad_item.shape[0]
                if n_cells < min_cells_for_benchmark:
                    logger.warning(f"[{method}] Only {n_cells} cells in crop region — "
                                   f"results may be unreliable (min recommended: {min_cells_for_benchmark})")
                for ct_col in ['leiden', 'celltype_majority', 'Class', 'celltype']:
                    if ct_col in ad_item.obs.columns:
                        n_types = ad_item.obs[ct_col].nunique()
                        logger.info(f"[{method}] Crop region: {n_cells:,} cells, "
                                    f"{n_types} {ct_col} categories")
                        break

            _run_benchmark_pipeline(cropped_list, config, crop_dir, "crop_region",
                                    crop_coords_um=crop_coords_um)
        else:
            logger.warning("No methods have cells in the Baysor crop region.")

    else:
        # ── SINGLE BENCHMARK MODE ──
        # All available methods on full tissue (includes Baysor if full/tiled mode)
        logger.info("=" * 60)
        logger.info("SINGLE BENCHMARK MODE: all methods on full tissue")
        logger.info("=" * 60)
        _run_benchmark_pipeline(adata_list, config, output_dir, "main")

    logger.info("Step 6 Completed Successfully.")


if __name__ == "__main__":
    pass
