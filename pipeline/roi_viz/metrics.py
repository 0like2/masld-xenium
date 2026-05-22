"""Per-(ROI, method) metric computation and aggregation.

Each (ROI, method) tuple gets a single metric dict, persisted under
``_cache/metrics/{roi_id}_{method}_metrics.pkl``. ``aggregate_metrics`` walks
the cache and writes ``roi_metrics.csv`` atomically.

Metrics computed (all NaN if not derivable):

* ``n_cells``                  — number of polygons.
* ``n_transcripts``            — total tx in ROI.
* ``transcripts_per_cell``     — median.
* ``genes_per_cell``           — median (unique genes per cell).
* ``assigned_fraction``        — fraction of tx with non-UNASSIGNED cell_id.
* ``unassigned_fraction``      — 1 − assigned_fraction.
* ``mean_cell_area`` /
  ``median_cell_area``         — µm² from polygons.
* ``nuclei_per_cell_mean``     — for expansion methods, count of original
                                  Xenium nuclei polygons that fall inside
                                  each method polygon (mean).
* ``negative_marker_purity``   — Salas NMP: 1 − fraction of cells expressing
                                  cell-type-conflicting markers (placeholder
                                  using competing-marker exclusivity from
                                  ``markers.py``).
* ``mixed_marker_fraction``    — fraction of cells with mutually exclusive
                                  marker pairs both > threshold.
* ``ssam_agreement``           — fraction of method centroids whose closest
                                  SSAM celltype matches the polygon's
                                  cell-type label (if both available;
                                  NaN otherwise).
* ``mean_vsi`` / ``low_vsi_fraction`` — sampled from the ROI VSI map at
                                  cell-centroid pixels.
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import cache
from .cropping import crop_p2r, crop_ssam, crop_transcripts, crop_vsi
from .data_classes import ROIRecord, SampleData
from .loaders import load_method, RASTER_PX_PER_UM
from .markers import select_markers, DEFAULT_BRAIN_MARKERS
from .ssam_utils import (coarsen_celltype, method_celltype_per_polygon,
                            ssam_majority_in_polygon)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-method metrics
# ---------------------------------------------------------------------------
def _polygon_area_um2(g) -> float:
    if g is None or g.is_empty:
        return float("nan")
    try:
        return float(g.area)
    except Exception:
        return float("nan")


def _compute_metrics_impl(sample: SampleData, roi: ROIRecord, method: str) -> Dict[str, Any]:
    from .style import METHOD_GROUPS
    res = load_method(method, sample, roi)
    out: Dict[str, Any] = {
        "roi_id": roi.roi_id,
        "sample_id": sample.sample_id,
        "roi_class": roi.roi_class,
        "method": method,
        "method_group": METHOD_GROUPS.get(method, "unknown"),
        "available": bool(res["available"]),
    }
    if not res["available"]:
        return out

    polys = res["polygons"]
    tx = res["transcripts"]
    n_cells = len(polys)
    n_tx = int(len(tx))
    out["n_cells"] = n_cells
    out["n_transcripts"] = n_tx

    if "geometry" in polys.columns and n_cells:
        areas = polys["geometry"].apply(_polygon_area_um2).dropna().values
        out["mean_cell_area"] = float(np.nanmean(areas)) if areas.size else float("nan")
        out["median_cell_area"] = float(np.nanmedian(areas)) if areas.size else float("nan")
    else:
        out["mean_cell_area"] = float("nan")
        out["median_cell_area"] = float("nan")

    # Per-cell aggregates from transcripts
    if n_tx and "cell_id" in tx.columns:
        tx_assigned = tx[tx["assigned"]]
        out["assigned_fraction"] = float(tx["assigned"].mean()) if n_tx else float("nan")
        out["unassigned_fraction"] = 1.0 - out["assigned_fraction"]
        if len(tx_assigned):
            grp = tx_assigned.groupby("cell_id")
            out["transcripts_per_cell"] = float(grp.size().median())
            out["genes_per_cell"] = float(
                grp["feature_name"].nunique().median())
        else:
            out["transcripts_per_cell"] = 0.0
            out["genes_per_cell"] = 0.0
    else:
        out["assigned_fraction"] = float("nan")
        out["unassigned_fraction"] = float("nan")
        out["transcripts_per_cell"] = float("nan")
        out["genes_per_cell"] = float("nan")

    # Nuclei-per-cell (for expansion-style methods, intersect with xenium nuclei)
    out["nuclei_per_cell_mean"] = float("nan")
    if method in ("rigid_expansion", "optimal_expansion", "baysor"):
        try:
            xn = load_method("xenium_nucleus", sample, roi)
            if xn["available"] and "geometry" in xn["polygons"].columns and n_cells:
                import geopandas as gpd
                base_g = gpd.GeoDataFrame({"geometry": xn["polygons"]["geometry"]},
                                            geometry="geometry")
                method_g = gpd.GeoDataFrame({"geometry": polys["geometry"]},
                                              geometry="geometry")
                # spatial join (intersects)
                base_g["__base_idx"] = np.arange(len(base_g))
                method_g["__method_idx"] = np.arange(len(method_g))
                joined = gpd.sjoin(base_g, method_g,
                                    how="inner", predicate="intersects")
                if len(joined):
                    counts = joined.groupby("__method_idx").size().values
                    out["nuclei_per_cell_mean"] = float(np.nanmean(counts))
                else:
                    out["nuclei_per_cell_mean"] = 0.0
        except Exception as e:
            logger.warning(f"[{method}] nuclei_per_cell_mean failed: {e}")

    # Mixed-marker fraction (mutually exclusive marker pairs in same cell)
    out["mixed_marker_fraction"] = _mixed_marker_fraction(tx)
    out["negative_marker_purity"] = _negative_marker_purity(tx)

    # SSAM agreement / cell-type disagreement
    # `cell_type_disagreement_fraction` is computed first; ssam_agreement is
    # then defined as 1 - that value (when comparable). When SSAM data is
    # unavailable both are NaN.
    disagreement = _cell_type_disagreement_fraction(sample, roi, method, res)
    out["cell_type_disagreement_fraction"] = disagreement
    out["ssam_agreement"] = (1.0 - disagreement) if np.isfinite(disagreement) \
        else _ssam_agreement(sample, roi, res)

    # Boundary sharpness (architecture-relevant; NaN otherwise)
    out["boundary_sharpness"] = _boundary_sharpness(sample, roi)

    # VSI sampling at centroids
    out["mean_vsi"], out["low_vsi_fraction"] = _vsi_at_centroids(sample, roi, res)

    return out


compute_metrics = cache.cached(
    "metrics",
    name_fn=lambda sample, roi, method: f"{roi.roi_id}_{method}_metrics",
    deps_fn=lambda sample, roi, method: (
        [tuple(roi.bbox_um), method]
        + _method_source_paths_for_metrics(sample, method)
    ),
    roi_id_fn=lambda sample, roi, method: roi.roi_id,
)(_compute_metrics_impl)


def _method_source_paths_for_metrics(sample, method):
    """Mirror of diff_compute._method_source_paths so metric cache also
    invalidates when the underlying segmentation outputs change."""
    from .diff_compute import _method_source_paths
    return _method_source_paths(sample, method)


# ---------------------------------------------------------------------------
# Sub-metric helpers
# ---------------------------------------------------------------------------
_CONFLICT_PAIRS = [
    ("SLC17A7", "GAD1"),     # excitatory vs inhibitory
    ("MOBP", "GFAP"),         # oligo vs astro
    ("P2RY12", "GFAP"),       # microglia vs astro
    ("OLIG2", "RBFOX3"),      # oligo vs neuron
]


def _mixed_marker_fraction(tx: pd.DataFrame) -> float:
    if tx is None or len(tx) == 0 or "cell_id" not in tx.columns:
        return float("nan")
    sub = tx[tx["assigned"]]
    if len(sub) == 0:
        return float("nan")
    grp = sub.groupby("cell_id")["feature_name"].apply(set)
    n_cells = len(grp)
    if n_cells == 0:
        return float("nan")
    n_mixed = 0
    for genes in grp:
        for a, b in _CONFLICT_PAIRS:
            if a in genes and b in genes:
                n_mixed += 1
                break
    return n_mixed / n_cells


def _negative_marker_purity(tx: pd.DataFrame) -> float:
    """Approximate Salas NMP: fraction of cells whose dominant cell-type
    inferred from positive markers does NOT also express the competing
    marker. Returns NaN if no marker info usable.
    """
    if tx is None or len(tx) == 0:
        return float("nan")
    sub = tx[tx["assigned"]]
    if len(sub) == 0:
        return float("nan")
    grp = sub.groupby("cell_id")["feature_name"].apply(list)
    if len(grp) == 0:
        return float("nan")
    pure = 0
    counted = 0
    for genes in grp:
        gset = set(genes)
        # A cell's dominant marker assigns it to a class; competing classes
        # appearing in the same cell count as impure.
        cls_genes = []
        for cls, ms in DEFAULT_BRAIN_MARKERS.items():
            cls_genes.append((cls, [m for m in ms if m in gset]))
        cls_genes = [(c, gs) for c, gs in cls_genes if gs]
        if len(cls_genes) == 0:
            continue
        counted += 1
        # Pure if only one class's markers are present.
        if len(cls_genes) == 1:
            pure += 1
    return pure / counted if counted else float("nan")


def _ssam_agreement(sample: SampleData, roi: ROIRecord,
                     method_result: Dict[str, Any]) -> float:
    cents = method_result.get("centroids")
    if cents is None or len(cents) == 0:
        return float("nan")
    try:
        ssam_df = crop_ssam(sample, roi)
    except Exception:
        return float("nan")
    if len(ssam_df) == 0:
        return float("nan")
    # Use leiden cluster as proxy for celltype if no celltype col.
    label_col = next((c for c in ("celltype", "ssam_celltype",
                                     "leiden_assignment", "leiden")
                      if c in ssam_df.columns), None)
    if label_col is None:
        return float("nan")
    # For each centroid, find nearest SSAM point and check if a method-side
    # celltype label matches. We don't have method-side cell types yet so we
    # use a proxy: SSAM majority within a small radius vs centroid's cluster.
    # Without method-side cell types this is a one-sided diversity check.
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([ssam_df["x_um"].values,
                                      ssam_df["y_um"].values]))
    radius_um = 15.0
    counts = tree.query_ball_point(np.column_stack([cents["x_centroid_um"].values,
                                                       cents["y_centroid_um"].values]),
                                    r=radius_um)
    if len(counts) == 0:
        return float("nan")
    # Agreement proxy: fraction of cells whose neighborhood is dominated by
    # a single SSAM class (>= 60 %).
    agree = 0
    for idxs in counts:
        if len(idxs) == 0:
            continue
        labels = ssam_df.iloc[idxs][label_col].astype(str).values
        vc = pd.Series(labels).value_counts(normalize=True)
        if vc.iloc[0] >= 0.6:
            agree += 1
    return agree / len(counts)


def _cell_type_disagreement_fraction(sample: SampleData, roi: ROIRecord,
                                       method: str,
                                       res: Dict[str, Any]) -> float:
    """Fraction of method polygons whose dominant cell type disagrees with
    the SSAM majority cell type inside the same polygon. NaN if either side
    is undefined for this ROI.
    """
    polys = res.get("polygons")
    if polys is None or "geometry" not in polys.columns or len(polys) == 0:
        return float("nan")
    try:
        ssam_df = crop_ssam(sample, roi)
    except Exception:
        return float("nan")
    if len(ssam_df) == 0:
        return float("nan")
    label_col = next((c for c in ("celltype", "ssam_celltype",
                                     "leiden_assignment", "leiden")
                      if c in ssam_df.columns), None)
    if label_col is None:
        return float("nan")
    try:
        method_ct = method_celltype_per_polygon(
            sample, roi, method, polys, res.get("centroids", pd.DataFrame()),
            res.get("transcripts", pd.DataFrame()))
        ssam_ct = ssam_majority_in_polygon(ssam_df, polys, label_col)
    except Exception as e:
        logger.debug(f"disagreement compute failed for {method}: {e}")
        return float("nan")
    # Collapse both sides to the coarse class vocabulary (neuron / inhibitory /
    # astrocyte / oligodendrocyte / microglia / endothelial). SSAM and step6
    # use Allen-Brain subclass names (L6 IT / Sst / VLMC / …) while the
    # marker-fallback path can only emit coarse classes — the coarsening keeps
    # the comparison apples-to-apples.
    m = np.array([coarsen_celltype(v) for v in method_ct.astype(str).values])
    s = np.array([coarsen_celltype(v) for v in ssam_ct.astype(str).values])
    both = (m != "") & (s != "")
    if not both.any():
        return float("nan")
    return float((m[both] != s[both]).mean())


def _boundary_sharpness(sample: SampleData, roi: ROIRecord) -> float:
    """1 / (cross-over width in µm). 'Cross-over width' = total length along
    the ROI short axis where target / competing marker densities are within
    20 % of each other. Smaller width → sharper transition → higher score.
    NaN if no marker pair is available.
    """
    try:
        tx = crop_transcripts(sample, roi)
    except Exception:
        return float("nan")
    if len(tx) == 0:
        return float("nan")
    present = set(tx["feature_name"].astype(str).unique())
    pair = None
    for cls_a, cls_b in (("neuron", "astrocyte"),
                         ("neuron", "inhibitory"),
                         ("astrocyte", "oligodendrocyte")):
        ga = next((g for g in DEFAULT_BRAIN_MARKERS.get(cls_a, []) if g in present), None)
        gb = next((g for g in DEFAULT_BRAIN_MARKERS.get(cls_b, []) if g in present), None)
        if ga and gb:
            pair = (ga, gb); break
    if pair is None:
        return float("nan")
    g_a, g_b = pair
    x0, y0, x1, y1 = roi.bbox_um
    short = "y_location" if (y1 - y0) < (x1 - x0) else "x_location"
    short_lo = y0 if short == "y_location" else x0
    short_hi = y1 if short == "y_location" else x1
    bin_um = 5.0
    edges = np.arange(short_lo, short_hi + bin_um, bin_um)
    h_a, _ = np.histogram(tx[tx["feature_name"].astype(str) == g_a][short]
                           .astype(float).values, bins=edges)
    h_b, _ = np.histogram(tx[tx["feature_name"].astype(str) == g_b][short]
                           .astype(float).values, bins=edges)
    if h_a.sum() == 0 or h_b.sum() == 0:
        return float("nan")
    h_a = h_a / max(1, h_a.max())
    h_b = h_b / max(1, h_b.max())
    crossover = np.abs(h_a - h_b) < 0.20
    width_um = float(crossover.sum() * bin_um)
    return 1.0 / max(bin_um, width_um)


def _vsi_at_centroids(sample: SampleData, roi: ROIRecord,
                       method_result: Dict[str, Any]):
    cents = method_result.get("centroids")
    if cents is None or len(cents) == 0:
        return float("nan"), float("nan")
    try:
        vsi = crop_vsi(sample, roi)
    except Exception:
        return float("nan"), float("nan")
    if vsi.get("map") is None:
        return float("nan"), float("nan")
    arr = vsi["map"]
    extent = vsi["extent_um"]
    # extent: (x0, x1, y_max, y_min) per crop_vsi
    x0, x1, y1, y0 = extent  # caller convention; remap
    x_min = min(x0, x1)
    x_max = max(x0, x1)
    y_min = min(y0, y1)
    y_max = max(y0, y1)
    px_per_um_x = arr.shape[1] / max(1e-6, (x_max - x_min))
    px_per_um_y = arr.shape[0] / max(1e-6, (y_max - y_min))
    cols = ((cents["x_centroid_um"].values - x_min) * px_per_um_x).astype(int)
    rows = ((cents["y_centroid_um"].values - y_min) * px_per_um_y).astype(int)
    cols = np.clip(cols, 0, arr.shape[1] - 1)
    rows = np.clip(rows, 0, arr.shape[0] - 1)
    vals = arr[rows, cols]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(np.nanmean(vals)), float(np.mean(vals < 0.5))


# ---------------------------------------------------------------------------
# Aggregation → roi_metrics.csv
# ---------------------------------------------------------------------------
METRIC_COLUMNS = [
    "roi_id", "sample_id", "roi_class", "method", "method_group", "available",
    "n_cells", "n_transcripts", "transcripts_per_cell", "genes_per_cell",
    "assigned_fraction", "unassigned_fraction",
    "mean_cell_area", "median_cell_area", "nuclei_per_cell_mean",
    "negative_marker_purity", "mixed_marker_fraction",
    "ssam_agreement", "cell_type_disagreement_fraction", "boundary_sharpness",
    "mean_vsi", "low_vsi_fraction",
]


def _atomic_write_csv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.stem + ".", suffix=".csv", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as f:
            df.to_csv(f, index=False)
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def aggregate_metrics(sample: SampleData, rois: Iterable[ROIRecord],
                       methods: Sequence[str], outdir: Path) -> Path:
    """Compute (or load cached) metrics for every (ROI, method) and merge
    into ``outdir/roi_metrics.csv``. Existing rows for unchanged tuples are
    preserved.
    """
    rows: List[Dict[str, Any]] = []
    for roi in rois:
        for m in methods:
            rec = compute_metrics(sample, roi, m)
            row = {c: rec.get(c, np.nan) for c in METRIC_COLUMNS}
            rows.append(row)

    new_df = pd.DataFrame(rows, columns=METRIC_COLUMNS)
    out_csv = Path(outdir) / "roi_metrics.csv"

    if out_csv.exists():
        try:
            old_df = pd.read_csv(out_csv)
            new_keys = set(zip(new_df["roi_id"].astype(str),
                                 new_df["method"].astype(str)))
            old_keys = list(zip(old_df["roi_id"].astype(str),
                                  old_df["method"].astype(str)))
            keep_mask = [k not in new_keys for k in old_keys]
            keep = old_df.loc[keep_mask].copy()
            merged = pd.concat([keep, new_df], ignore_index=True)
        except Exception as e:
            logger.warning(f"failed to merge with existing roi_metrics.csv ({e}); overwriting")
            merged = new_df
    else:
        merged = new_df

    # Sort by roi_id, method for stable output.
    merged = merged.sort_values(["roi_id", "method"]).reset_index(drop=True)
    _atomic_write_csv(merged, out_csv)
    logger.info(f"roi_metrics.csv written: {len(merged)} rows → {out_csv}")
    return out_csv
