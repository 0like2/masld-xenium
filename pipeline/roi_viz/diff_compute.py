"""Heavy difference computations between baseline and method results.

These functions are cached at stage='diff' so that re-rendering the highlight
PNGs does not re-trigger the underlying numerical work.

Three primitives are produced per (baseline, method, ROI):

1. ``compute_difference_mask`` — categorical (H, W) int8 array:
       0 = unassigned in both (background)
       1 = baseline-only
       2 = method-only
       3 = shared
   This is the spatial set comparison used for the difference-mask render.

2. ``compute_reassignment`` — per-transcript category DataFrame with column
   ``category`` ∈ {unchanged, reassigned, newly, dropped}. This is the
   transcript-level comparison used for the reassignment scatter.

3. ``compute_boundary_overlap_meta`` — small dict of summary metrics
   (areas, IOU, transcript-category fractions) used by metric cards.

The colour composition of the difference mask into RGBA is in
``compose_difference_rgba`` — that is the *single imshow* trick used to
prevent matplotlib's alpha blending from making overlap regions look
black/grey (see plan §10.2).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors

from . import cache, style as st
from .data_classes import ROIRecord, SampleData
from .loaders import load_method, RASTER_PX_PER_UM


def _method_source_paths(sample: SampleData, method: str):
    """Return the on-disk source file paths whose mtime should drive cache
    invalidation for downstream artifacts (diff, metrics) that consume the
    output of ``load_method(method, ...)``.

    Without this, a re-run of (e.g.) Baysor in roi_targeted mode would NOT
    invalidate the diff/metrics caches because their fingerprints depend
    only on (bbox, method-name) and not on the method's source data.
    """
    paths = []
    if method == "xenium_nucleus":
        paths += [sample.xenium_nucleus_path, sample.transcripts_path]
    elif method in ("cellpose_nuclei", "rigid_expansion"):
        if sample.step3_dir:
            mname = ("nuclei_masks.tif" if method == "cellpose_nuclei"
                     else "resegmented_masks.tif")
            paths.append(sample.step3_dir / f"{sample.sample_tag}_step3_{mname}")
            paths.append(sample.step3_dir / f"{sample.sample_tag}_step3_transcripts_resegmented.csv")
    elif method == "optimal_expansion":
        if sample.step3_dir:
            paths.append(sample.step3_dir / f"{sample.sample_tag}_step3_nuclei_masks.tif")
        if sample.step5_dir:
            paths.append(sample.step5_dir / f"{sample.sample_tag}_step5_optimal_expansion.txt")
        paths.append(sample.transcripts_path)
    elif method == "baysor":
        if sample.step6_dir:
            base = sample.step6_dir / "baysor_run"
            paths += [
                base / "roi_targeted" / "segmentation.csv",
                base / "roi_targeted" / "segmentation_polygons_2d.json",
                base / "roi_targeted" / "roi_targeted_metadata.json",
                base / "crop" / "segmentation.csv",
                base / "crop" / "segmentation_polygons_2d.json",
                base / "segmentation.csv",
            ]
    return [str(p) for p in paths if p is not None]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cell-level set comparison
# ---------------------------------------------------------------------------
def _compute_difference_mask_impl(sample: SampleData, roi: ROIRecord,
                                    baseline: str, method: str) -> Dict[str, Any]:
    base = load_method(baseline, sample, roi)
    other = load_method(method, sample, roi)

    if not base["available"] or not other["available"]:
        return {"available": False, "reason": "method unavailable",
                 "category": None, "extent_um": None,
                 "area_summary": {}, "iou": None}

    A = base["raster"]
    B = other["raster"]
    if A is None or B is None:
        return {"available": False, "reason": "no raster",
                 "category": None, "extent_um": None,
                 "area_summary": {}, "iou": None}
    # Ensure same shape (raster_extent_um = roi.bbox_um for both via _make_raster_grid)
    if A.shape != B.shape:
        logger.error(f"raster shape mismatch {A.shape} vs {B.shape}")
        return {"available": False, "reason": "shape mismatch",
                 "category": None, "extent_um": None,
                 "area_summary": {}, "iou": None}

    a_filled = A > 0
    b_filled = B > 0
    cat = np.zeros(A.shape, dtype=np.int8)
    cat[a_filled & ~b_filled] = 1   # baseline-only
    cat[~a_filled & b_filled] = 2   # method-only
    cat[a_filled & b_filled] = 3    # shared

    total = cat.size
    n_a_only = int((cat == 1).sum())
    n_b_only = int((cat == 2).sum())
    n_shared = int((cat == 3).sum())
    n_bg = int((cat == 0).sum())
    union = n_a_only + n_b_only + n_shared
    iou = n_shared / union if union else 0.0

    return {
        "available": True,
        "category": cat,
        "extent_um": tuple(roi.bbox_um),
        "raster_px_per_um": RASTER_PX_PER_UM,
        "area_summary": {
            "baseline_only_frac": n_a_only / total,
            "method_only_frac":   n_b_only / total,
            "shared_frac":        n_shared / total,
            "background_frac":    n_bg / total,
            "iou":                iou,
        },
        "iou": iou,
        "baseline": baseline,
        "method": method,
    }


compute_difference_mask = cache.cached(
    "diff",
    name_fn=lambda sample, roi, baseline, method: f"{roi.roi_id}_diffmask_{baseline}__vs__{method}",
    deps_fn=lambda sample, roi, baseline, method: (
        [tuple(roi.bbox_um), baseline, method]
        + _method_source_paths(sample, baseline)
        + _method_source_paths(sample, method)
    ),
    roi_id_fn=lambda sample, roi, baseline, method: roi.roi_id,
)(_compute_difference_mask_impl)


# ---------------------------------------------------------------------------
# Transcript-level reassignment
# ---------------------------------------------------------------------------
REASSIGN_CATEGORIES = ("unchanged", "reassigned", "newly", "dropped")


def _compute_reassignment_impl(sample: SampleData, roi: ROIRecord,
                                baseline: str, method: str) -> Dict[str, Any]:
    base = load_method(baseline, sample, roi)
    other = load_method(method, sample, roi)

    if not base["available"] or not other["available"]:
        return {"available": False, "reason": "method unavailable",
                "frame": pd.DataFrame(), "fractions": {}}

    a = base["transcripts"]
    b = other["transcripts"]
    if a is None or b is None or len(a) == 0 or len(b) == 0:
        return {"available": False, "reason": "empty transcripts",
                "frame": pd.DataFrame(), "fractions": {}}

    # Inner-join on transcript_id (Baysor uses molecule_id which is in a
    # different namespace from the Xenium transcript_id, so for Baysor we
    # fall back to spatial nearest-neighbour join).
    a = a[["transcript_id", "cell_id", "assigned",
            "x_location", "y_location"]].copy()
    b = b[["transcript_id", "cell_id", "assigned",
            "x_location", "y_location"]].copy()

    if method == "baysor" or len(set(a["transcript_id"]) & set(b["transcript_id"])) == 0:
        # Spatial nearest-neighbour fallback
        from scipy.spatial import cKDTree
        tree = cKDTree(np.column_stack([b["x_location"].values,
                                          b["y_location"].values]))
        dists, idx = tree.query(np.column_stack([a["x_location"].values,
                                                   a["y_location"].values]),
                                 k=1)
        # Threshold: only accept matches within 1 µm
        valid = dists < 1.0
        joined = a.copy()
        joined["cell_id_b"] = "UNASSIGNED"
        joined.loc[valid, "cell_id_b"] = b["cell_id"].iloc[idx[valid]].values
        joined["assigned_b"] = joined["cell_id_b"] != "UNASSIGNED"
        joined = joined.rename(columns={"cell_id": "cell_id_a",
                                          "assigned": "assigned_a"})
    else:
        joined = a.merge(b, on="transcript_id", how="inner",
                         suffixes=("_a", "_b"))
        # Use baseline (a) coordinates for plotting
        joined["x_location"] = joined.get("x_location_a", joined.get("x_location"))
        joined["y_location"] = joined.get("y_location_a", joined.get("y_location"))

    cat = np.empty(len(joined), dtype=object)
    a_assigned = joined["assigned_a"].values
    b_assigned = joined["assigned_b"].values
    same_cell = (joined["cell_id_a"].astype(str).values ==
                  joined["cell_id_b"].astype(str).values)

    cat[a_assigned & b_assigned & same_cell] = "unchanged"
    cat[a_assigned & b_assigned & ~same_cell] = "reassigned"
    cat[~a_assigned & b_assigned] = "newly"
    cat[a_assigned & ~b_assigned] = "dropped"
    cat[~a_assigned & ~b_assigned] = "unchanged"
    joined["category"] = cat

    frac = {c: float((cat == c).sum() / len(cat)) for c in REASSIGN_CATEGORIES}
    return {
        "available": True,
        "frame": joined[["transcript_id", "x_location", "y_location",
                           "cell_id_a", "cell_id_b", "category"]].reset_index(drop=True),
        "fractions": frac,
        "n_total": len(joined),
        "baseline": baseline,
        "method": method,
    }


compute_reassignment = cache.cached(
    "diff",
    name_fn=lambda sample, roi, baseline, method: f"{roi.roi_id}_reassign_{baseline}__vs__{method}",
    deps_fn=lambda sample, roi, baseline, method: (
        [tuple(roi.bbox_um), baseline, method]
        + _method_source_paths(sample, baseline)
        + _method_source_paths(sample, method)
    ),
    roi_id_fn=lambda sample, roi, baseline, method: roi.roi_id,
)(_compute_reassignment_impl)


# ---------------------------------------------------------------------------
# RGBA compositing — key step that prevents alpha-stacking artifacts
# ---------------------------------------------------------------------------
def compose_difference_rgba(category: np.ndarray,
                              colors: Optional[Dict[int, Tuple[str, float]]] = None
                              ) -> np.ndarray:
    """Convert categorical (H, W) int8 into a single (H, W, 4) RGBA array.

    Default mapping:
        0 = transparent (background)
        1 = xenium_only red, alpha 0.65
        2 = method_only blue, alpha 0.65
        3 = shared gray, alpha 0.30
    """
    if colors is None:
        colors = {
            0: (st.HIGHLIGHT_COLORS["background"], 0.0),
            1: (st.HIGHLIGHT_COLORS["xenium_only"], 0.65),
            2: (st.HIGHLIGHT_COLORS["method_only"], 0.65),
            3: (st.HIGHLIGHT_COLORS["shared"], 0.30),
        }
    H, W = category.shape
    rgba = np.zeros((H, W, 4), dtype=np.float32)
    for k, (hex_color, alpha) in colors.items():
        mask = (category == k)
        if not mask.any():
            continue
        rgba[mask, :3] = mcolors.to_rgb(hex_color)
        rgba[mask, 3] = alpha
    return rgba
