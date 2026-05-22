"""ROI-level data cropping helpers.

All cropping functions are cached via ``cache.cached`` so that repeated calls
(in the same or later pipeline runs) skip work as long as input file
fingerprints are unchanged.

Coordinate convention: ``bbox_um = (x_min, y_min, x_max, y_max)`` in µm.
The Xenium convention used by this project is ``um_per_pixel_inv = 4.70588``
(pixels per µm).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import cache
from .data_classes import ROIRecord, SampleData

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _bbox_to_pixel_slice(bbox_um: Tuple[float, float, float, float],
                          um_per_pixel_inv: float,
                          image_shape: Tuple[int, int]) -> Tuple[slice, slice]:
    """Return (row_slice, col_slice) clipped to image bounds. y → row, x → col."""
    x0, y0, x1, y1 = bbox_um
    H, W = image_shape
    row0 = max(0, int(np.floor(y0 * um_per_pixel_inv)))
    row1 = min(H, int(np.ceil(y1 * um_per_pixel_inv)))
    col0 = max(0, int(np.floor(x0 * um_per_pixel_inv)))
    col1 = min(W, int(np.ceil(x1 * um_per_pixel_inv)))
    return slice(row0, row1), slice(col0, col1)


def _decode_bytes_col(s: pd.Series) -> pd.Series:
    """Decode a pandas Series of bytes to str (Xenium parquet stores bytes)."""
    if len(s) == 0:
        return s
    sample = s.iloc[0]
    if isinstance(sample, (bytes, bytearray)):
        return s.str.decode("utf-8")
    return s.astype(str)


# ---------------------------------------------------------------------------
# DAPI / morphology image
# ---------------------------------------------------------------------------
def _crop_dapi_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    """Read a DAPI window from the OME-TIF for *roi*.

    Returns a dict with the cropped uint16 array, the µm bbox of the array,
    pixel coordinates, and a percentile-stretched float copy for display.
    """
    if sample.dapi_path is None or not sample.dapi_path.exists():
        return {"image": None, "bbox_um": tuple(roi.bbox_um), "pixel_bbox": None}

    import tifffile as tf
    with tf.TiffFile(sample.dapi_path) as t:
        page = t.pages[0]
        H, W = page.shape[-2:]
        rs, cs = _bbox_to_pixel_slice(roi.bbox_um, sample.um_per_pixel_inv, (H, W))
        img = page.asarray()
        if img.ndim == 3:
            img = img[0]
        crop = img[rs, cs]

    if crop.size == 0:
        logger.warning(f"crop_dapi: empty crop for {roi.roi_id} bbox={roi.bbox_um}")
        return {"image": None, "bbox_um": tuple(roi.bbox_um), "pixel_bbox": None}

    return {
        "image": np.asarray(crop, dtype=np.float32),
        "bbox_um": tuple(roi.bbox_um),
        "pixel_bbox": (rs.start, cs.start, rs.stop, cs.stop),
        "um_per_pixel_inv": sample.um_per_pixel_inv,
    }


crop_dapi = cache.cached(
    "raw",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_dapi_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.dapi_path or ""),
                                  tuple(roi.bbox_um),
                                  sample.um_per_pixel_inv],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_crop_dapi_impl)


# ---------------------------------------------------------------------------
# Transcripts (full parquet → ROI subset)
# ---------------------------------------------------------------------------
def _crop_transcripts_impl(sample: SampleData, roi: ROIRecord) -> pd.DataFrame:
    """Subset the full Xenium transcripts.parquet to the ROI bbox.

    Returns a DataFrame with at minimum: transcript_id, feature_name (str),
    cell_id (str, may be 'UNASSIGNED'), x_location, y_location, z_location,
    overlaps_nucleus, qv, nucleus_distance.
    """
    if sample.transcripts_path is None or not sample.transcripts_path.exists():
        return pd.DataFrame()

    cols = ["transcript_id", "cell_id", "overlaps_nucleus", "feature_name",
            "x_location", "y_location", "z_location", "qv", "nucleus_distance"]
    df = pd.read_parquet(sample.transcripts_path, columns=cols)

    x0, y0, x1, y1 = roi.bbox_um
    mask = ((df["x_location"] >= x0) & (df["x_location"] <= x1) &
            (df["y_location"] >= y0) & (df["y_location"] <= y1))
    df = df.loc[mask].copy()

    if "feature_name" in df.columns:
        df["feature_name"] = _decode_bytes_col(df["feature_name"])
    if "cell_id" in df.columns:
        df["cell_id"] = _decode_bytes_col(df["cell_id"])

    return df.reset_index(drop=True)


crop_transcripts = cache.cached(
    "raw",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_tx_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.transcripts_path or ""),
                                  tuple(roi.bbox_um)],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_crop_transcripts_impl)


# ---------------------------------------------------------------------------
# Label-mask cropping (returns label patch + pixel bbox)
# ---------------------------------------------------------------------------
def _crop_label_mask_impl(mask_path: str, roi_bbox: Tuple[float, float, float, float],
                           um_per_pixel_inv: float) -> Dict[str, Any]:
    if not mask_path or not Path(mask_path).exists():
        return {"mask": None, "pixel_bbox": None}
    import tifffile as tf

    # Try memmap first — for large untiled label TIFs this is ~100x faster
    # than asarray() because we only touch the bbox region's bytes.
    rs = cs = None
    crop = None
    try:
        m = tf.memmap(mask_path)
        if m.ndim == 3:
            m = m[0]
        H, W = m.shape
        rs, cs = _bbox_to_pixel_slice(roi_bbox, um_per_pixel_inv, (H, W))
        crop = np.asarray(m[rs, cs], dtype=np.int32)
        del m
    except Exception:
        # Fall back to full read.
        with tf.TiffFile(mask_path) as t:
            page = t.pages[0]
            H, W = page.shape[-2:]
            rs, cs = _bbox_to_pixel_slice(roi_bbox, um_per_pixel_inv, (H, W))
            full = page.asarray()
            if full.ndim == 3:
                full = full[0]
            crop = np.asarray(full[rs, cs], dtype=np.int32)

    return {
        "mask": crop,
        "pixel_bbox": (rs.start, cs.start, rs.stop, cs.stop),
        "um_per_pixel_inv": um_per_pixel_inv,
    }


crop_label_mask = cache.cached(
    "raw",
    name_fn=lambda mask_path, roi_bbox, um_per_pixel_inv:
        f"labelmask_{Path(mask_path).stem}_{int(roi_bbox[0])}_{int(roi_bbox[1])}",
    deps_fn=lambda mask_path, roi_bbox, um_per_pixel_inv:
        [mask_path, tuple(roi_bbox), um_per_pixel_inv],
)(_crop_label_mask_impl)


# ---------------------------------------------------------------------------
# P2R classification table
# ---------------------------------------------------------------------------
def _crop_p2r_impl(sample: SampleData, roi: ROIRecord) -> pd.DataFrame:
    """Per-spot P2R cluster assignment within the ROI.

    Preferred source: ``sample.p2r_h5ad_path`` (points2regions.h5ad), whose
    obs has spatial columns and per-spot ``p2r_name`` / ``p2r_compartment``.

    Fallback: ``sample.p2r_csv_path`` is per-cluster metadata only and has
    no spatial info — return empty DataFrame in that case (caller handles).
    """
    h5ad_path = sample.p2r_h5ad_path
    if h5ad_path is not None and Path(h5ad_path).exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(h5ad_path)
            obs = a.obs.copy()
            x_col = next((c for c in ("x_location", "x_um", "x", "x_centroid")
                          if c in obs.columns), None)
            y_col = next((c for c in ("y_location", "y_um", "y", "y_centroid")
                          if c in obs.columns), None)
            if x_col is None or y_col is None:
                logger.warning("crop_p2r: P2R h5ad has no spatial cols (%s)",
                               list(obs.columns)[:8])
                return pd.DataFrame()
            x = obs[x_col].astype(float)
            y = obs[y_col].astype(float)
            coord_range = max(x.max() - x.min(), y.max() - y.min())
            if coord_range > 20000:
                x = x / sample.um_per_pixel_inv
                y = y / sample.um_per_pixel_inv
            obs = obs.assign(_x_um=x, _y_um=y)
            x0, y0, x1, y1 = roi.bbox_um
            mask = ((obs["_x_um"] >= x0) & (obs["_x_um"] <= x1) &
                    (obs["_y_um"] >= y0) & (obs["_y_um"] <= y1))
            sub = obs.loc[mask].copy()
            keep = ["_x_um", "_y_um"] + [c for c in (
                "p2r_name", "p2r_compartment", "p2r_celltype",
                "compartment", "celltype", "leiden") if c in sub.columns]
            sub = sub[keep].rename(columns={"_x_um": "x_um", "_y_um": "y_um"})
            return sub.reset_index(drop=True)
        except Exception as e:
            logger.warning(f"crop_p2r: h5ad load failed ({e}); returning empty")
            return pd.DataFrame()

    # CSV-only path: no spatial info — return empty so renderer falls back.
    logger.warning("crop_p2r: only cluster-metadata CSV available; "
                   "P2R reference rendering requires p2r_h5ad_path.")
    return pd.DataFrame()


crop_p2r = cache.cached(
    "raw",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_p2r_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.p2r_csv_path or ""),
                                  str(sample.p2r_h5ad_path or ""),
                                  tuple(roi.bbox_um)],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_crop_p2r_impl)


# ---------------------------------------------------------------------------
# SSAM coordinate transform: obs[x,y] (px on 2 µm grid, origin-shifted) → µm
# ---------------------------------------------------------------------------
_SSAM_TRANSFORM_CACHE: Dict[str, Tuple[float, float, float]] = {}


def _ssam_obs_to_um(sample: SampleData, x_vals, y_vals):
    """Apply the (`um_per_px`, `x_min`, `y_min`) transform that step2 used
    when it built the SSAM AnnData. Result is in real µm.

    The offset is the per-sample minimum of the (panel-filtered) transcript
    coordinates — i.e. the same shift step2 applies before dividing by
    ``ssam_vf_um_per_px`` (default 2.0). We compute it once per sample and
    cache the tuple so subsequent ROIs are cheap.
    """
    key = sample.sample_tag
    cached = _SSAM_TRANSFORM_CACHE.get(key)
    if cached is None:
        um_per_px, x_min, y_min = _detect_ssam_transform(sample)
        _SSAM_TRANSFORM_CACHE[key] = (um_per_px, x_min, y_min)
    else:
        um_per_px, x_min, y_min = cached
    return x_vals * um_per_px + x_min, y_vals * um_per_px + y_min


def _detect_ssam_transform(sample: SampleData) -> Tuple[float, float, float]:
    """Determine ``(um_per_px, x_min, y_min)`` for the SSAM h5ad of *sample*.

    Strategy:
    1. Default ``um_per_px = 2.0`` (matches the step2 config and the SSAM
       library's default `ssam_vf_um_per_px`).
    2. ``x_min`` / ``y_min`` come from the transcript table — the
       ``transcripts.parquet`` columns ``x_location`` / ``y_location`` minimum.
    """
    um_per_px = 2.0
    x_min = 0.0
    y_min = 0.0
    try:
        if sample.transcripts_path and sample.transcripts_path.exists():
            tx = pd.read_parquet(sample.transcripts_path,
                                  columns=["x_location", "y_location"])
            x_min = float(tx["x_location"].min())
            y_min = float(tx["y_location"].min())
    except Exception as e:
        logger.debug(f"SSAM transform: transcripts read failed ({e}); "
                       f"using x_min=y_min=0")
    return um_per_px, x_min, y_min


# ---------------------------------------------------------------------------
# SSAM (h5ad → DataFrame of (x, y, leiden, celltype))
# ---------------------------------------------------------------------------
def _crop_ssam_impl(sample: SampleData, roi: ROIRecord) -> pd.DataFrame:
    if sample.ssam_h5ad_path is None or not sample.ssam_h5ad_path.exists():
        return pd.DataFrame()
    try:
        import anndata as ad
    except ImportError:
        logger.warning("anndata not available — SSAM crop skipped")
        return pd.DataFrame()
    a = ad.read_h5ad(sample.ssam_h5ad_path)
    obs = a.obs.copy()

    x_col = "x" if "x" in obs.columns else None
    y_col = "y" if "y" in obs.columns else None
    if x_col is None or y_col is None:
        return pd.DataFrame()
    x_vals = obs[x_col].astype(float)
    y_vals = obs[y_col].astype(float)

    # SSAM stores obs x/y in *vector-field pixels* on a downsampled grid,
    # shifted to (0,0). Real-µm coordinate is::
    #
    #     x_um = obs["x"] * um_per_px_grid + spots_x_min
    #     y_um = obs["y"] * um_per_px_grid + spots_y_min
    #
    # See ``_propagate_ssam_celltypes_to_spots`` in step2 (cf. lines 5677-5678).
    x_um, y_um = _ssam_obs_to_um(sample, x_vals, y_vals)
    obs = obs.assign(x_um=x_um, y_um=y_um)
    x0, y0, x1, y1 = roi.bbox_um
    mask = ((obs["x_um"] >= x0) & (obs["x_um"] <= x1) &
            (obs["y_um"] >= y0) & (obs["y_um"] <= y1))
    sub = obs.loc[mask].copy()
    cols = ["x_um", "y_um"] + [c for c in ("leiden", "leiden_assignment",
                                            "celltype", "ssam_celltype")
                                if c in sub.columns]
    return sub[cols].reset_index(drop=True)


crop_ssam = cache.cached(
    "raw",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_ssam_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.ssam_h5ad_path or ""),
                                  tuple(roi.bbox_um)],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_crop_ssam_impl)


# ---------------------------------------------------------------------------
# VSI / ovrlpy coherence map (npz cache)
# ---------------------------------------------------------------------------
def _crop_vsi_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    if sample.vsi_map_path is None or not sample.vsi_map_path.exists():
        return {"map": None, "extent_um": None}
    npz = np.load(sample.vsi_map_path, allow_pickle=True)
    keys = list(npz.files)
    arr = None
    for cand in ("vsi_map", "coherence_map", "integrity_map", "map"):
        if cand in keys:
            arr = npz[cand]
            break
    if arr is None:
        # fallback: first 2D float array
        for k in keys:
            v = npz[k]
            if hasattr(v, "ndim") and v.ndim == 2 and v.dtype.kind == "f":
                arr = v
                break
    if arr is None:
        logger.warning("crop_vsi: no recognizable map in %s (keys=%s)",
                       sample.vsi_map_path, keys)
        return {"map": None, "extent_um": None}

    # Determine grid pixel size: try common keys, else assume 2 µm/pixel.
    px_um = float(npz["um_per_pixel"]) if "um_per_pixel" in keys else 2.0
    H, W = arr.shape
    x_max_um = W * px_um
    y_max_um = H * px_um
    x0, y0, x1, y1 = roi.bbox_um
    rs = slice(max(0, int(y0 / px_um)), min(H, int(np.ceil(y1 / px_um))))
    cs = slice(max(0, int(x0 / px_um)), min(W, int(np.ceil(x1 / px_um))))
    crop = np.asarray(arr[rs, cs], dtype=np.float32)
    extent = (cs.start * px_um, cs.stop * px_um,
              rs.stop * px_um, rs.start * px_um)
    return {"map": crop, "extent_um": extent, "um_per_pixel_grid": px_um}


crop_vsi = cache.cached(
    "raw",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_vsi_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.vsi_map_path or ""),
                                  tuple(roi.bbox_um)],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_crop_vsi_impl)
