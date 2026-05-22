"""Per-method loaders.

Each ``load_method_*`` function returns a uniform dict:
    {
      "method":   str,
      "polygons": GeoDataFrame[cell_id, geometry] in µm  (may be empty),
      "raster":   ndarray int32 (H, W) on a fixed-resolution ROI grid (or None),
      "raster_extent_um": (x0, y0, x1, y1) of the raster,
      "transcripts": DataFrame[transcript_id, cell_id, assigned (bool),
                                x_location, y_location, feature_name],
      "centroids":   DataFrame[cell_id, x_centroid_um, y_centroid_um] (may be empty),
    }

Heavy operations (mask → polygon, expand_labels, geometry filtering) are
cached under stage='boundaries'.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import cache
from .cropping import (crop_dapi, crop_label_mask, crop_p2r, crop_ssam,
                        crop_transcripts, crop_vsi, _decode_bytes_col)
from .data_classes import ROIRecord, SampleData

logger = logging.getLogger(__name__)


# Resolution of the per-ROI raster grid used for difference masks (px / µm).
RASTER_PX_PER_UM = 2.0


def _empty_result(method: str) -> Dict[str, Any]:
    return {
        "method": method,
        "polygons": _empty_gdf(),
        "raster": None,
        "raster_extent_um": None,
        "transcripts": pd.DataFrame(columns=["transcript_id", "cell_id",
                                              "assigned", "x_location",
                                              "y_location", "feature_name"]),
        "centroids": pd.DataFrame(columns=["cell_id", "x_centroid_um",
                                            "y_centroid_um"]),
        "available": False,
    }


def _empty_gdf():
    try:
        import geopandas as gpd
        from shapely.geometry import Polygon
        return gpd.GeoDataFrame({"cell_id": pd.Series(dtype=str),
                                  "geometry": pd.Series(dtype=object)},
                                 geometry="geometry")
    except ImportError:
        return pd.DataFrame({"cell_id": [], "geometry": []})


def _make_raster_grid(roi: ROIRecord) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
    """Empty (H, W) int32 raster spanning the ROI at RASTER_PX_PER_UM."""
    x0, y0, x1, y1 = roi.bbox_um
    W = max(1, int(round((x1 - x0) * RASTER_PX_PER_UM)))
    H = max(1, int(round((y1 - y0) * RASTER_PX_PER_UM)))
    return np.zeros((H, W), dtype=np.int32), (x0, y0, x1, y1)


# ---------------------------------------------------------------------------
# 1. Xenium 10x default (nucleus boundaries parquet)
# ---------------------------------------------------------------------------
def _load_xenium_nucleus_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    method = "xenium_nucleus"
    nb_path = sample.xenium_nucleus_path
    if nb_path is None or not Path(nb_path).exists():
        logger.warning(f"[{method}] missing nucleus boundaries — skipping")
        return _empty_result(method)

    nb = pd.read_parquet(nb_path)
    if "cell_id" in nb.columns:
        nb["cell_id"] = _decode_bytes_col(nb["cell_id"])

    # Filter polygons that have *any* vertex inside the ROI bbox.
    x0, y0, x1, y1 = roi.bbox_um
    in_roi = ((nb["vertex_x"] >= x0) & (nb["vertex_x"] <= x1) &
              (nb["vertex_y"] >= y0) & (nb["vertex_y"] <= y1))
    cells_in_roi = nb.loc[in_roi, "cell_id"].unique()
    nb_sub = nb[nb["cell_id"].isin(cells_in_roi)]

    try:
        import geopandas as gpd
        from shapely.geometry import Polygon

        def _build_poly(g):
            xs = g["vertex_x"].values
            ys = g["vertex_y"].values
            if len(xs) < 3:
                return None
            return Polygon(zip(xs, ys))

        polys_obj = (nb_sub.groupby("cell_id")[["vertex_x", "vertex_y"]]
                            .apply(_build_poly))
        # Newer pandas can return a DataFrame from apply when the function
        # returns a non-scalar; coerce to a Series so reset_index has a
        # ``name=`` argument available.
        if isinstance(polys_obj, pd.DataFrame):
            if polys_obj.shape[1] >= 1:
                polys_obj = polys_obj.iloc[:, 0]
            else:
                polys_obj = pd.Series([], dtype=object, name="geometry")
        polys = polys_obj.rename("geometry").reset_index()
        polys = polys[polys["geometry"].notna()].reset_index(drop=True)
        gdf = gpd.GeoDataFrame(polys, geometry="geometry")
    except ImportError:
        gdf = pd.DataFrame({"cell_id": list(cells_in_roi), "geometry": [None] * len(cells_in_roi)})

    centroids = (nb_sub.groupby("cell_id")
                        .agg(x_centroid_um=("vertex_x", "mean"),
                             y_centroid_um=("vertex_y", "mean"))
                        .reset_index())

    raster, extent = _make_raster_grid(roi)
    if "geometry" in gdf.columns and len(gdf):
        labels = np.arange(1, len(gdf) + 1, dtype=np.int32)
        raster, n_drawn = _rasterize_polygons(
            gdf["geometry"].values, labels, raster.shape, extent, RASTER_PX_PER_UM)
        _sanity_check_raster(method, len(gdf), raster, n_drawn)
        gdf["raster_label"] = labels

    # Transcript assignment: use base transcripts.parquet 'cell_id' field.
    tx_full = crop_transcripts(sample, roi)
    if len(tx_full):
        tx = tx_full[["transcript_id", "cell_id", "feature_name",
                       "x_location", "y_location"]].copy()
        tx["assigned"] = tx["cell_id"].astype(str).str.upper() != "UNASSIGNED"
    else:
        tx = pd.DataFrame(columns=["transcript_id", "cell_id", "assigned",
                                    "x_location", "y_location", "feature_name"])

    return {
        "method": method,
        "polygons": gdf,
        "raster": raster,
        "raster_extent_um": extent,
        "transcripts": tx,
        "centroids": centroids,
        "available": True,
    }


load_xenium_nucleus = cache.cached(
    "boundaries",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_xenium_nucleus_{roi.roi_id}",
    deps_fn=lambda sample, roi: [str(sample.xenium_nucleus_path or ""),
                                  str(sample.transcripts_path or ""),
                                  tuple(roi.bbox_um)],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_load_xenium_nucleus_impl)


# ---------------------------------------------------------------------------
# 2. Cellpose nuclei (Step 3 nuclei_masks.tif)
# ---------------------------------------------------------------------------
def _label_mask_to_polygons_and_centroids(mask: np.ndarray,
                                            extent_um: Tuple[float, float, float, float],
                                            um_per_pixel_inv: float) -> Tuple[Any, pd.DataFrame, np.ndarray]:
    """Convert a label mask → per-cell polygons (µm) + centroids + ROI raster."""
    try:
        import geopandas as gpd
        from shapely.geometry import Polygon
        from skimage.measure import regionprops, find_contours
    except ImportError as e:
        logger.warning(f"polygon conversion deps missing: {e}")
        gdf_empty = pd.DataFrame({"cell_id": [], "geometry": []})
        return gdf_empty, pd.DataFrame(columns=["cell_id", "x_centroid_um", "y_centroid_um"]), \
            np.zeros((1, 1), dtype=np.int32)

    rows = []
    cents = []

    # extent_um is in µm coords (relative to whole tissue).
    x0_um, y0_um, x1_um, y1_um = extent_um

    # mask is at native resolution (1 px = 1/um_per_pixel_inv µm).
    px_to_um = 1.0 / um_per_pixel_inv

    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]
    if len(unique_labels) == 0:
        gdf = pd.DataFrame({"cell_id": [], "geometry": []})
        cents_df = pd.DataFrame(columns=["cell_id", "x_centroid_um", "y_centroid_um"])
        return gdf, cents_df, np.zeros((1, 1), dtype=np.int32)

    rp = {p.label: p for p in regionprops(mask)}
    for lbl in unique_labels:
        prop = rp.get(int(lbl))
        if prop is None:
            continue
        # Build polygon from contour of binary slice
        minr, minc, maxr, maxc = prop.bbox
        sub = (mask[minr:maxr, minc:maxc] == lbl).astype(np.uint8)
        if sub.size == 0:
            continue
        try:
            contours = find_contours(sub, level=0.5)
        except Exception:
            contours = []
        if not contours:
            continue
        # Largest contour
        contour = max(contours, key=len)
        # contour is (row, col) within sub; shift to global pixel, then µm
        global_rows = contour[:, 0] + minr
        global_cols = contour[:, 1] + minc
        xs_um = global_cols * px_to_um + x0_um  # only valid if mask was cropped at exact x0_um
        ys_um = global_rows * px_to_um + y0_um
        if len(xs_um) < 3:
            continue
        try:
            poly = Polygon(np.column_stack([xs_um, ys_um]))
            if not poly.is_valid:
                poly = poly.buffer(0)
        except Exception:
            continue
        cell_id = f"L{int(lbl)}"
        rows.append({"cell_id": cell_id, "geometry": poly,
                     "raster_label": int(lbl)})
        cy, cx = prop.centroid  # row, col
        cents.append({"cell_id": cell_id,
                       "x_centroid_um": cx * px_to_um + x0_um,
                       "y_centroid_um": cy * px_to_um + y0_um})

    try:
        import geopandas as gpd
        gdf = gpd.GeoDataFrame(rows, geometry="geometry") if rows else \
            gpd.GeoDataFrame({"cell_id": pd.Series(dtype=str),
                                "geometry": pd.Series(dtype=object)}, geometry="geometry")
    except ImportError:
        gdf = pd.DataFrame(rows) if rows else pd.DataFrame({"cell_id": [], "geometry": []})

    cents_df = pd.DataFrame(cents) if cents else pd.DataFrame(
        columns=["cell_id", "x_centroid_um", "y_centroid_um"])

    # ROI raster at RASTER_PX_PER_UM resolution (resample mask).
    return gdf, cents_df, mask  # caller resamples


def _resample_mask_to_roi_grid(mask: np.ndarray,
                                src_extent_um: Tuple[float, float, float, float],
                                dst_extent_um: Tuple[float, float, float, float],
                                src_px_per_um: float,
                                dst_px_per_um: float = RASTER_PX_PER_UM) -> np.ndarray:
    """Resample a label mask (cropped at src_extent_um) into the ROI raster grid."""
    x0d, y0d, x1d, y1d = dst_extent_um
    x0s, y0s, x1s, y1s = src_extent_um
    Wd = int(round((x1d - x0d) * dst_px_per_um))
    Hd = int(round((y1d - y0d) * dst_px_per_um))
    if Wd <= 0 or Hd <= 0:
        return np.zeros((1, 1), dtype=np.int32)
    out = np.zeros((Hd, Wd), dtype=np.int32)
    Hs, Ws = mask.shape
    # For each dst pixel, look up source pixel via nearest sampling.
    ys_dst = np.arange(Hd)
    xs_dst = np.arange(Wd)
    y_um = y0d + (ys_dst + 0.5) / dst_px_per_um
    x_um = x0d + (xs_dst + 0.5) / dst_px_per_um
    src_rows = ((y_um - y0s) * src_px_per_um).astype(np.int64)
    src_cols = ((x_um - x0s) * src_px_per_um).astype(np.int64)
    src_rows = np.clip(src_rows, 0, Hs - 1)
    src_cols = np.clip(src_cols, 0, Ws - 1)
    out[:, :] = mask[np.ix_(src_rows, src_cols)]
    return out


def _rasterize_polygons(geometries, labels, out_shape, extent_um, px_per_um):
    """Rasterize a list of shapely polygons into a (H, W) int32 label image.

    Single backend (skimage.draw.polygon) — no silent fallback so behavior is
    reproducible. The caller should sanity-check the returned raster against
    the polygon count via :func:`_sanity_check_raster`.
    """
    from skimage.draw import polygon as skpoly  # required dep

    H, W = out_shape
    raster = np.zeros((H, W), dtype=np.int32)
    x0, y0, _, _ = extent_um

    n_drawn = 0
    for geom, lbl in zip(geometries, labels):
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            rings = [list(geom.exterior.coords)]
        elif geom.geom_type == "MultiPolygon":
            rings = [list(p.exterior.coords) for p in geom.geoms]
        else:
            continue
        for ring in rings:
            xs, ys = zip(*ring)
            xs = np.asarray(xs, dtype=float)
            ys = np.asarray(ys, dtype=float)
            cols = (xs - x0) * px_per_um
            rows = (ys - y0) * px_per_um
            rr, cc = skpoly(rows, cols, shape=(H, W))
            if rr.size > 0:
                raster[rr, cc] = lbl
                n_drawn += 1
    return raster, n_drawn


def _sanity_check_raster(method: str, n_polygons: int, raster: np.ndarray,
                          n_drawn: int) -> None:
    """Loud warning if rasterization produced no pixels despite having input polygons.

    This catches silent failures (e.g., coordinate frame mismatch, empty
    contours, all polygons clipped out) instead of letting downstream
    difference-mask code silently report 100 % shared area.
    """
    n_filled = int((raster > 0).sum())
    if n_polygons > 0 and n_filled == 0:
        logger.error(f"[{method}] rasterize produced 0 filled pixels for "
                     f"{n_polygons} polygons (n_drawn={n_drawn}); "
                     "likely a coordinate frame mismatch — DOWNSTREAM DIFFERENCE "
                     "MASKS WILL BE WRONG")
    elif n_polygons > 0 and n_drawn == 0:
        logger.error(f"[{method}] no polygons were attempted in rasterizer "
                     "despite n_polygons>0")


def _build_method_result_from_label_mask(method: str,
                                          mask_path: Path,
                                          sample: SampleData,
                                          roi: ROIRecord,
                                          tx_assignment_col: Optional[str]) -> Dict[str, Any]:
    """Generic: load a Step 3 label TIF, build polygons / raster / transcripts."""
    if mask_path is None or not Path(mask_path).exists():
        logger.warning(f"[{method}] missing mask: {mask_path}")
        return _empty_result(method)

    crop = crop_label_mask(str(mask_path), tuple(roi.bbox_um), sample.um_per_pixel_inv)
    mask_crop = crop.get("mask")
    if mask_crop is None or mask_crop.size == 0:
        return _empty_result(method)

    px_to_um = 1.0 / sample.um_per_pixel_inv
    pix_bbox = crop.get("pixel_bbox")
    src_extent_um = (
        pix_bbox[1] * px_to_um, pix_bbox[0] * px_to_um,
        pix_bbox[3] * px_to_um, pix_bbox[2] * px_to_um,
    )

    gdf, cents_df, _ = _label_mask_to_polygons_and_centroids(
        mask_crop, src_extent_um, sample.um_per_pixel_inv)

    # Build ROI raster at RASTER_PX_PER_UM
    raster_extent = tuple(roi.bbox_um)
    raster = _resample_mask_to_roi_grid(
        mask_crop, src_extent_um, raster_extent,
        sample.um_per_pixel_inv, RASTER_PX_PER_UM)

    # Transcript assignment from Step 3 csv (if column provided)
    tx = pd.DataFrame()
    step3_tx_path = sample.step3_dir / f"{sample.sample_tag}_step3_transcripts_resegmented.csv" \
        if sample.step3_dir else None
    if tx_assignment_col and step3_tx_path and step3_tx_path.exists():
        cols = ["transcript_id", "feature_name", "x_location", "y_location",
                tx_assignment_col]
        try:
            df = pd.read_csv(step3_tx_path, usecols=cols, low_memory=False)
        except Exception:
            df = pd.read_csv(step3_tx_path, low_memory=False)
            df = df[[c for c in cols if c in df.columns]]
        x0, y0, x1, y1 = roi.bbox_um
        m = ((df["x_location"] >= x0) & (df["x_location"] <= x1) &
             (df["y_location"] >= y0) & (df["y_location"] <= y1))
        df = df.loc[m].copy()
        df["cell_id"] = df[tx_assignment_col].apply(
            lambda v: f"L{int(v)}" if pd.notna(v) and v not in (0, "0", "UNASSIGNED") else "UNASSIGNED")
        df["assigned"] = df["cell_id"] != "UNASSIGNED"
        tx = df[["transcript_id", "cell_id", "assigned",
                  "x_location", "y_location", "feature_name"]].reset_index(drop=True)
    else:
        # Fallback: derive assignment by sampling the raster at each tx position.
        tx_full = crop_transcripts(sample, roi)
        if len(tx_full):
            xs = tx_full["x_location"].values
            ys = tx_full["y_location"].values
            x0, y0, x1, y1 = roi.bbox_um
            cols = ((xs - x0) * RASTER_PX_PER_UM).astype(np.int64)
            rows = ((ys - y0) * RASTER_PX_PER_UM).astype(np.int64)
            cols = np.clip(cols, 0, raster.shape[1] - 1)
            rows = np.clip(rows, 0, raster.shape[0] - 1)
            labels = raster[rows, cols]
            tx = pd.DataFrame({
                "transcript_id": tx_full["transcript_id"].values,
                "cell_id": [f"L{int(l)}" if l > 0 else "UNASSIGNED" for l in labels],
                "feature_name": tx_full["feature_name"].values,
                "x_location": xs,
                "y_location": ys,
            })
            tx["assigned"] = tx["cell_id"] != "UNASSIGNED"

    return {
        "method": method,
        "polygons": gdf,
        "raster": raster,
        "raster_extent_um": raster_extent,
        "transcripts": tx,
        "centroids": cents_df,
        "available": True,
    }


def _load_cellpose_nuclei_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    mask_path = sample.step3_dir / f"{sample.sample_tag}_step3_nuclei_masks.tif" if sample.step3_dir else None
    return _build_method_result_from_label_mask(
        "cellpose_nuclei", mask_path, sample, roi, tx_assignment_col="in_cell")


load_cellpose_nuclei = cache.cached(
    "boundaries",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_cellpose_{roi.roi_id}",
    deps_fn=lambda sample, roi: [
        str((sample.step3_dir / f"{sample.sample_tag}_step3_nuclei_masks.tif")
            if sample.step3_dir else ""),
        str((sample.step3_dir / f"{sample.sample_tag}_step3_transcripts_resegmented.csv")
            if sample.step3_dir else ""),
        tuple(roi.bbox_um), sample.um_per_pixel_inv,
    ],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_load_cellpose_nuclei_impl)


# ---------------------------------------------------------------------------
# 3. Rigid expansion (Step 3 expand_labels d=400px)
# ---------------------------------------------------------------------------
def _load_rigid_expansion_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    mask_path = sample.step3_dir / f"{sample.sample_tag}_step3_resegmented_masks.tif" if sample.step3_dir else None
    return _build_method_result_from_label_mask(
        "rigid_expansion", mask_path, sample, roi, tx_assignment_col="closest_cell")


load_rigid_expansion = cache.cached(
    "boundaries",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_rigid_{roi.roi_id}",
    deps_fn=lambda sample, roi: [
        str((sample.step3_dir / f"{sample.sample_tag}_step3_resegmented_masks.tif")
            if sample.step3_dir else ""),
        str((sample.step3_dir / f"{sample.sample_tag}_step3_transcripts_resegmented.csv")
            if sample.step3_dir else ""),
        tuple(roi.bbox_um), sample.um_per_pixel_inv,
    ],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_load_rigid_expansion_impl)


# ---------------------------------------------------------------------------
# 4. Optimal expansion (apply expand_labels with optimal radius µm)
# ---------------------------------------------------------------------------
def _load_optimal_expansion_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    method = "optimal_expansion"
    if sample.step5_dir is None:
        return _empty_result(method)
    optimal_txt = sample.step5_dir / f"{sample.sample_tag}_step5_optimal_expansion.txt"
    if not optimal_txt.exists():
        logger.warning(f"[{method}] missing optimal_expansion.txt")
        return _empty_result(method)
    # Parse optimal radius (µm)
    radius_um = None
    try:
        with open(optimal_txt) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("optimal_expansion"):
                    parts = line.split()
                    radius_um = float(parts[-1])
                    break
    except Exception as e:
        logger.warning(f"[{method}] cannot parse {optimal_txt}: {e}")
    if radius_um is None or radius_um <= 0:
        return _empty_result(method)

    # Load nuclei mask crop, expand by radius_um pixels.
    nuclei_path = sample.step3_dir / f"{sample.sample_tag}_step3_nuclei_masks.tif" \
        if sample.step3_dir else None
    if nuclei_path is None or not nuclei_path.exists():
        return _empty_result(method)

    crop = crop_label_mask(str(nuclei_path), tuple(roi.bbox_um), sample.um_per_pixel_inv)
    mask_crop = crop.get("mask")
    if mask_crop is None or mask_crop.size == 0:
        return _empty_result(method)

    radius_px = int(round(radius_um * sample.um_per_pixel_inv))
    try:
        from skimage.segmentation import expand_labels
        expanded = expand_labels(mask_crop, distance=max(1, radius_px))
    except Exception as e:
        logger.warning(f"[{method}] expand_labels failed: {e}")
        expanded = mask_crop

    expanded = expanded.astype(np.int32)

    px_to_um = 1.0 / sample.um_per_pixel_inv
    pix_bbox = crop.get("pixel_bbox")
    src_extent_um = (
        pix_bbox[1] * px_to_um, pix_bbox[0] * px_to_um,
        pix_bbox[3] * px_to_um, pix_bbox[2] * px_to_um,
    )

    gdf, cents_df, _ = _label_mask_to_polygons_and_centroids(
        expanded, src_extent_um, sample.um_per_pixel_inv)

    raster = _resample_mask_to_roi_grid(
        expanded, src_extent_um, tuple(roi.bbox_um),
        sample.um_per_pixel_inv, RASTER_PX_PER_UM)

    # Transcript assignment by raster lookup.
    tx_full = crop_transcripts(sample, roi)
    if len(tx_full):
        xs = tx_full["x_location"].values
        ys = tx_full["y_location"].values
        x0, y0, x1, y1 = roi.bbox_um
        cols = np.clip(((xs - x0) * RASTER_PX_PER_UM).astype(np.int64), 0, raster.shape[1] - 1)
        rows = np.clip(((ys - y0) * RASTER_PX_PER_UM).astype(np.int64), 0, raster.shape[0] - 1)
        labels = raster[rows, cols]
        tx = pd.DataFrame({
            "transcript_id": tx_full["transcript_id"].values,
            "cell_id": [f"L{int(l)}" if l > 0 else "UNASSIGNED" for l in labels],
            "feature_name": tx_full["feature_name"].values,
            "x_location": xs,
            "y_location": ys,
        })
        tx["assigned"] = tx["cell_id"] != "UNASSIGNED"
    else:
        tx = pd.DataFrame(columns=["transcript_id", "cell_id", "assigned",
                                    "x_location", "y_location", "feature_name"])

    return {
        "method": method,
        "polygons": gdf,
        "raster": raster,
        "raster_extent_um": tuple(roi.bbox_um),
        "transcripts": tx,
        "centroids": cents_df,
        "extras": {"optimal_radius_um": radius_um},
        "available": True,
    }


load_optimal_expansion = cache.cached(
    "boundaries",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_optimal_{roi.roi_id}",
    deps_fn=lambda sample, roi: [
        str((sample.step3_dir / f"{sample.sample_tag}_step3_nuclei_masks.tif")
            if sample.step3_dir else ""),
        str((sample.step5_dir / f"{sample.sample_tag}_step5_optimal_expansion.txt")
            if sample.step5_dir else ""),
        str(sample.transcripts_path or ""),
        tuple(roi.bbox_um), sample.um_per_pixel_inv,
    ],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_load_optimal_expansion_impl)


# ---------------------------------------------------------------------------
# 5. Baysor (Step 6)
# ---------------------------------------------------------------------------
def _load_baysor_impl(sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    method = "baysor"
    if sample.step6_dir is None:
        return _empty_result(method)

    # Resolution priority (mirrors xenium_step6_segmentation_benchmark.py
    # _load_baysor_polygons): roi_targeted/ → crop/ → baysor_run/.
    roi_targeted_dir = sample.step6_dir / "baysor_run" / "roi_targeted"
    crop_dir = sample.step6_dir / "baysor_run" / "crop"
    bare_dir = sample.step6_dir / "baysor_run"

    baysor_dir = None
    mode = None
    for cand_dir, cand_mode in ((roi_targeted_dir, "roi_targeted"),
                                  (crop_dir, "crop"),
                                  (bare_dir, "full")):
        if (cand_dir / "segmentation_polygons_2d.json").exists() and \
           (cand_dir / "segmentation.csv").exists():
            baysor_dir = cand_dir
            mode = cand_mode
            break

    if baysor_dir is None:
        logger.warning(f"[{method}] missing Baysor outputs (no segmentation in "
                       f"roi_targeted/ or crop/ or baysor_run/)")
        return _empty_result(method)

    polys_path = baysor_dir / "segmentation_polygons_2d.json"
    seg_path = baysor_dir / "segmentation.csv"

    # roi_targeted mode: every manifest ROI is by construction inside its
    # padded region — but check the metadata to confirm THIS ROI was actually
    # processed (vs. pulled in from a different manifest). For crop mode we
    # still apply the bbox-inside-crop test.
    if mode == "roi_targeted":
        meta_path = baysor_dir / "roi_targeted_metadata.json"
        if meta_path.exists():
            try:
                with open(meta_path) as f:
                    meta = json.load(f)
                roi_ids_in_run = {rid
                                    for region in meta.get("regions", [])
                                    for rid in region.get("roi_ids", [])}
                if roi_ids_in_run and roi.roi_id not in roi_ids_in_run:
                    logger.info(f"[{method}] ROI {roi.roi_id} not in "
                                f"roi_targeted run (found {len(roi_ids_in_run)} "
                                "other ROIs) — skipping")
                    return _empty_result(method)
            except Exception as e:
                logger.warning(f"[{method}] roi_targeted_metadata read failed: {e}")
    elif mode == "crop":
        crop_meta_path = baysor_dir / "crop_metadata.json"
        if crop_meta_path.exists():
            try:
                with open(crop_meta_path) as f:
                    cm = json.load(f)
                cx0, cx1, cy0, cy1 = cm.get("coords", [None, None, None, None])
                if cx0 is not None:
                    rx0, ry0, rx1, ry1 = roi.bbox_um
                    if not (rx0 >= cx0 and rx1 <= cx1 and ry0 >= cy0 and ry1 <= cy1):
                        logger.warning(f"[{method}] ROI {roi.roi_id} outside Baysor crop "
                                       f"x=[{cx0:.0f},{cx1:.0f}] y=[{cy0:.0f},{cy1:.0f}] — skipping")
                        return _empty_result(method)
            except Exception as e:
                logger.warning(f"[{method}] failed to read crop_metadata: {e}")

    try:
        import geopandas as gpd
        from shapely.geometry import Polygon, shape
    except ImportError:
        return _empty_result(method)

    with open(polys_path) as f:
        gj = json.load(f)
    feats = gj.get("features", [])
    rows = []
    x0, y0, x1, y1 = roi.bbox_um
    for ft in feats:
        try:
            geom = shape(ft.get("geometry"))
        except Exception:
            continue
        if not geom.is_valid:
            geom = geom.buffer(0)
        bx0, by0, bx1, by1 = geom.bounds
        if bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1:
            continue
        rows.append({"cell_id": str(ft.get("id") or len(rows)),
                      "geometry": geom})
    gdf = gpd.GeoDataFrame(rows, geometry="geometry") if rows else \
        gpd.GeoDataFrame({"cell_id": pd.Series(dtype=str),
                            "geometry": pd.Series(dtype=object)},
                          geometry="geometry")
    if len(gdf):
        cents = gdf.copy()
        cents["x_centroid_um"] = cents.geometry.centroid.x
        cents["y_centroid_um"] = cents.geometry.centroid.y
        cents_df = cents[["cell_id", "x_centroid_um", "y_centroid_um"]]
    else:
        cents_df = pd.DataFrame(columns=["cell_id", "x_centroid_um", "y_centroid_um"])

    # Rasterize polygons into ROI grid for difference mask.
    raster, extent = _make_raster_grid(roi)
    if len(gdf):
        labels = np.arange(1, len(gdf) + 1, dtype=np.int32)
        raster, n_drawn = _rasterize_polygons(
            gdf["geometry"].values, labels, raster.shape, extent, RASTER_PX_PER_UM)
        _sanity_check_raster(method, len(gdf), raster, n_drawn)
        gdf = gdf.assign(raster_label=labels)
        cell_id_to_label = dict(zip(gdf["cell_id"].values, labels))
    else:
        cell_id_to_label = {}

    # Transcript assignment from Baysor segmentation.csv.
    seg = pd.read_csv(seg_path, low_memory=False)
    tx = pd.DataFrame()
    if {"x", "y", "molecule_id", "cell"}.issubset(seg.columns):
        m = ((seg["x"] >= x0) & (seg["x"] <= x1) &
             (seg["y"] >= y0) & (seg["y"] <= y1))
        sub = seg.loc[m].copy()
        # Match to xenium transcript_id by spatial nearest-neighbor (Baysor uses molecule_id).
        # For the difference framework we just use Baysor's molecule_id as transcript_id surrogate.
        tx = pd.DataFrame({
            "transcript_id": sub["molecule_id"].astype("int64"),
            "cell_id": sub["cell"].fillna("UNASSIGNED").astype(str),
            "feature_name": sub["gene"].astype(str) if "gene" in sub.columns else "",
            "x_location": sub["x"].astype(float),
            "y_location": sub["y"].astype(float),
        })
        tx["assigned"] = (tx["cell_id"] != "UNASSIGNED") & (tx["cell_id"].astype(str) != "")

    return {
        "method": method,
        "polygons": gdf,
        "raster": raster,
        "raster_extent_um": extent,
        "transcripts": tx,
        "centroids": cents_df,
        "available": True,
    }


load_baysor = cache.cached(
    "boundaries",
    name_fn=lambda sample, roi: f"{sample.sample_tag}_baysor_{roi.roi_id}",
    deps_fn=lambda sample, roi: [
        # roi_targeted artifacts (preferred when present)
        str((sample.step6_dir / "baysor_run" / "roi_targeted" / "segmentation_polygons_2d.json")
            if sample.step6_dir else ""),
        str((sample.step6_dir / "baysor_run" / "roi_targeted" / "segmentation.csv")
            if sample.step6_dir else ""),
        str((sample.step6_dir / "baysor_run" / "roi_targeted" / "roi_targeted_metadata.json")
            if sample.step6_dir else ""),
        # crop fallback
        str((sample.step6_dir / "baysor_run" / "crop" / "segmentation_polygons_2d.json")
            if sample.step6_dir else ""),
        str((sample.step6_dir / "baysor_run" / "crop" / "segmentation.csv")
            if sample.step6_dir else ""),
        tuple(roi.bbox_um),
    ],
    roi_id_fn=lambda sample, roi: roi.roi_id,
)(_load_baysor_impl)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------
METHOD_LOADERS = {
    "xenium_nucleus":     load_xenium_nucleus,
    "cellpose_nuclei":    load_cellpose_nuclei,
    "rigid_expansion":    load_rigid_expansion,
    "optimal_expansion":  load_optimal_expansion,
    "baysor":             load_baysor,
}


def load_method(method: str, sample: SampleData, roi: ROIRecord) -> Dict[str, Any]:
    fn = METHOD_LOADERS.get(method)
    if fn is None:
        raise ValueError(f"Unknown segmentation method: {method}")
    return fn(sample, roi)
