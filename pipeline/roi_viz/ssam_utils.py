"""SSAM-based cell-type comparison helpers (architecture ROI add-ons).

Encapsulates the hybrid ``method-side cell type`` resolution used by
``render_architecture_addons``:

1. **Primary** — read ``step6_benchmark/benchmark_combined.h5ad`` (full tissue)
   and ``step6_benchmark/crop_comparison/benchmark_combined.h5ad`` (crop only)
   when present, look up cells inside the ROI bbox by ``segmentation`` value.
2. **Fallback** — derive a per-polygon majority cell type from transcript
   markers using ``markers.DEFAULT_BRAIN_MARKERS`` (gene → class reverse map).

Also provides:

- ``ssam_majority_in_polygon`` — for each cell polygon, the majority SSAM
  celltype label among SSAM points it contains.
- ``grid_majority`` — bin a ``(x, y, label)`` table into a uniform grid and
  return per-cell majority label as an int array (``-1`` = empty).
- ``stable_celltype_palette`` — deterministic celltype → hex color map.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .data_classes import ROIRecord, SampleData
from .markers import DEFAULT_BRAIN_MARKERS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Method ↔ step6 segmentation column mapping
# ---------------------------------------------------------------------------
METHOD_TO_STEP6_SEG = {
    "xenium_nucleus":    "nuclei",
    "cellpose_nuclei":   "cellpose",
    "rigid_expansion":   "expansion",
    "baysor":            "baysor",
    # `optimal_expansion` not in step6 — falls back to marker majority
}


def _gene_to_class_map() -> Dict[str, str]:
    """Reverse ``DEFAULT_BRAIN_MARKERS`` to gene → class. Excludes the
    ``nuclear`` / ``cytoplasmic`` compartment categories which are subcellular,
    not cell-type identities.
    """
    out: Dict[str, str] = {}
    for cls, genes in DEFAULT_BRAIN_MARKERS.items():
        if cls in ("nuclear", "cytoplasmic"):
            continue
        for g in genes:
            out.setdefault(g, cls)
    return out


_GENE_TO_CLASS = _gene_to_class_map()


# ---------------------------------------------------------------------------
# Granular SSAM/step6 celltype → coarse class collapse
# ---------------------------------------------------------------------------
# SSAM ``leiden_assignment`` and step6 ``celltype_majority`` use Allen-Brain
# subclass names (L6 IT / Sst / Pvalb / VLMC / OPC / Microglia-PVM / …) while
# the marker-majority fallback can only produce coarse classes (neuron /
# astrocyte / oligodendrocyte / microglia / endothelial). To keep
# disagreement comparable across all methods (including ``optimal_expansion``
# which has no step6 row) we collapse both sides to the coarse vocabulary.
_GRANULAR_TO_COARSE = {
    # Excitatory neuron subclasses
    "l2/3 it":   "neuron",
    "l4 it":     "neuron",
    "l5 it":     "neuron",
    "l5 et":     "neuron",
    "l5/6 np":   "neuron",
    "l6 it":     "neuron",
    "l6 it car3":"neuron",
    "l6 ct":     "neuron",
    "l6b":       "neuron",
    "car3":      "neuron",
    # Inhibitory subclasses
    "lamp5":     "inhibitory",
    "sncg":      "inhibitory",
    "vip":       "inhibitory",
    "sst":       "inhibitory",
    "sst chodl": "inhibitory",
    "pvalb":     "inhibitory",
    "pax6":      "inhibitory",
    "chandelier":"inhibitory",
    # Glia
    "astrocyte":      "astrocyte",
    "oligodendrocyte":"oligodendrocyte",
    "opc":            "oligodendrocyte",
    "microglia-pvm":  "microglia",
    "microglia":      "microglia",
    # Vascular / non-neuronal
    "endothelial":    "endothelial",
    "vlmc":           "endothelial",
    "pericyte":       "endothelial",
    "smc":            "endothelial",
    # Coarse marker classes already in canonical form
    "neuron":          "neuron",
    "inhibitory":      "inhibitory",
    "ad_pathology":    "ad_pathology",
}


def coarsen_celltype(label: str) -> str:
    """Collapse a granular celltype label to a coarse class. Empty strings
    pass through; unknown labels return the lower-cased original (so direct
    equality still works when both sides happen to use the same exact name).
    """
    if label is None:
        return ""
    s = str(label).strip().lower()
    if not s:
        return ""
    return _GRANULAR_TO_COARSE.get(s, s)


# ---------------------------------------------------------------------------
# Step6 cell-type loader
# ---------------------------------------------------------------------------
def load_step6_celltypes(sample: SampleData, roi: ROIRecord,
                          method: str) -> pd.DataFrame:
    """Load (x_centroid, y_centroid, celltype) for *method* inside ROI bbox.

    Tries the full-tissue h5ad first, then the crop-only h5ad (which is the
    only place ``baysor`` lives). Returns empty DataFrame if not available.
    """
    seg_value = METHOD_TO_STEP6_SEG.get(method)
    if seg_value is None:
        return pd.DataFrame(columns=["x", "y", "celltype"])

    paths: List[Any] = []
    p_full = getattr(sample, "step6_combined_h5ad", None)
    p_crop = getattr(sample, "step6_crop_combined_h5ad", None)
    if p_full is not None:
        paths.append(p_full)
    if p_crop is not None:
        paths.append(p_crop)

    if not paths:
        return pd.DataFrame(columns=["x", "y", "celltype"])

    try:
        import anndata as ad
    except ImportError:
        logger.warning("anndata unavailable — step6 celltype lookup skipped")
        return pd.DataFrame(columns=["x", "y", "celltype"])

    x0, y0, x1, y1 = roi.bbox_um
    for p in paths:
        try:
            if p is None or not (hasattr(p, "exists") and p.exists()):
                continue
            a = ad.read_h5ad(str(p), backed="r")
        except Exception as e:
            logger.debug(f"step6 h5ad read failed ({p}): {e}")
            continue
        obs = a.obs
        if "segmentation" not in obs.columns or "celltype_majority" not in obs.columns:
            continue
        if "x_centroid" not in obs.columns or "y_centroid" not in obs.columns:
            continue
        sub = obs[obs["segmentation"].astype(str) == seg_value]
        if len(sub) == 0:
            continue
        x = sub["x_centroid"].astype(float).values
        y = sub["y_centroid"].astype(float).values
        in_box = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        if not in_box.any():
            continue
        return pd.DataFrame({
            "x": x[in_box],
            "y": y[in_box],
            "celltype": sub["celltype_majority"].astype(str).values[in_box],
        }).reset_index(drop=True)

    return pd.DataFrame(columns=["x", "y", "celltype"])


# ---------------------------------------------------------------------------
# Marker-majority fallback
# ---------------------------------------------------------------------------
def marker_majority_celltype(tx: pd.DataFrame,
                              polygons: pd.DataFrame) -> pd.Series:
    """Return ``Series`` indexed by polygon ``cell_id`` (or polygon order if
    no cell_id) → marker-majority celltype string. Empty string when no
    marker transcript available.
    """
    if tx is None or len(tx) == 0 or "cell_id" not in tx.columns:
        # No transcript→cell mapping — return empties for every polygon.
        return _empty_polygon_series(polygons)

    sub = tx[tx["assigned"]] if "assigned" in tx.columns else tx
    if len(sub) == 0:
        return _empty_polygon_series(polygons)

    # Map gene → class
    cls = sub["feature_name"].astype(str).map(_GENE_TO_CLASS)
    sub2 = sub.assign(_class=cls).dropna(subset=["_class"])
    if len(sub2) == 0:
        return _empty_polygon_series(polygons)

    grp = sub2.groupby([sub2["cell_id"].astype(str), "_class"]).size()
    grp = grp.reset_index(name="n")
    # idxmax per cell
    out = (grp.sort_values(["cell_id", "n"], ascending=[True, False])
              .drop_duplicates("cell_id")
              .set_index("cell_id")["_class"])

    # Align to polygons
    if "cell_id" in polygons.columns:
        idx = polygons["cell_id"].astype(str)
    else:
        idx = pd.Series([str(i) for i in range(len(polygons))])
    return out.reindex(idx).fillna("").reset_index(drop=True)


def _empty_polygon_series(polygons: pd.DataFrame) -> pd.Series:
    n = len(polygons)
    return pd.Series([""] * n, name="celltype")


# ---------------------------------------------------------------------------
# Hybrid resolver — used by every architecture add-on
# ---------------------------------------------------------------------------
def method_celltype_per_polygon(sample: SampleData, roi: ROIRecord,
                                  method: str,
                                  polygons: pd.DataFrame,
                                  centroids: pd.DataFrame,
                                  tx: pd.DataFrame) -> pd.Series:
    """Return a Series aligned to *polygons* row order with each polygon's
    cell type. Uses step6 celltype if available, else marker-majority
    fallback. The Series uses string celltype values; empty string means
    unresolved.
    """
    n = len(polygons)
    if n == 0:
        return pd.Series([], dtype=str)

    s6 = load_step6_celltypes(sample, roi, method)
    used_step6 = False
    out = pd.Series([""] * n, name="celltype")

    if len(s6) and len(centroids):
        # Nearest-neighbour match per polygon centroid. Use a generous tolerance
        # (25 µm ≈ a typical cell radius) — step6 cell centroids come from a
        # different segmentation pass than this loader's polygons, so exact
        # equality is unrealistic.
        try:
            from scipy.spatial import cKDTree
            tree = cKDTree(np.column_stack([s6["x"].values, s6["y"].values]))
            cx = centroids["x_centroid_um"].astype(float).values
            cy = centroids["y_centroid_um"].astype(float).values
            d, idx = tree.query(np.column_stack([cx, cy]), distance_upper_bound=25.0)
            ok = np.isfinite(d) & (idx < len(s6))
            if ok.any():
                # Align centroid order → polygon order via cell_id (if both have it)
                ct = np.array([""] * len(centroids), dtype=object)
                ct[ok] = s6["celltype"].values[idx[ok]]
                # Map centroid cell_id → polygon order
                if "cell_id" in centroids.columns and "cell_id" in polygons.columns:
                    cid_to_ct = dict(zip(centroids["cell_id"].astype(str), ct))
                    out = pd.Series(
                        [cid_to_ct.get(c, "") for c in polygons["cell_id"].astype(str)],
                        name="celltype")
                else:
                    out = pd.Series(ct[:n] if len(ct) >= n else
                                    list(ct) + [""] * (n - len(ct)),
                                    name="celltype")
                used_step6 = bool(out.astype(str).str.len().sum() > 0)
        except Exception as e:
            logger.debug(f"step6 NN match failed for {method}: {e}")

    if not used_step6 or (out == "").mean() > 0.5:
        # Fallback (or supplement) with marker-majority for unfilled rows.
        marker_ct = marker_majority_celltype(tx, polygons)
        if not used_step6:
            logger.info(f"[ssam_utils] step6 celltype miss for {method} → marker fallback")
            out = marker_ct.astype(str)
        else:
            empty_mask = (out.astype(str) == "").values
            mc = marker_ct.astype(str).values
            out = out.copy()
            out.iloc[empty_mask] = mc[empty_mask]
    return out.astype(str)


# ---------------------------------------------------------------------------
# SSAM majority within each polygon
# ---------------------------------------------------------------------------
def ssam_majority_in_polygon(ssam_df: pd.DataFrame,
                              polygons: pd.DataFrame,
                              label_col: Optional[str] = None,
                              radius_um: float = 20.0) -> pd.Series:
    """For each polygon, return the majority SSAM label among nearby points.

    Uses a centroid + radius lookup (default 20 µm ≈ typical brain cell
    diameter) rather than a strict ``within`` join — nuclear-only polygons
    are tiny (~8 µm radius) and almost never *contain* a SSAM KDE sample,
    even though SSAM clearly assigns the area to a cell type.
    """
    n = len(polygons)
    if n == 0 or len(ssam_df) == 0:
        return pd.Series([""] * n, dtype=str)

    if label_col is None:
        for c in ("celltype", "ssam_celltype", "leiden_assignment", "leiden"):
            if c in ssam_df.columns:
                label_col = c
                break
    if label_col is None:
        return pd.Series([""] * n, dtype=str)

    if "geometry" not in polygons.columns:
        return pd.Series([""] * n, dtype=str)

    # Compute polygon centroids
    cx = np.zeros(n); cy = np.zeros(n); valid = np.zeros(n, dtype=bool)
    for i, geom in enumerate(polygons["geometry"].values):
        if geom is None or geom.is_empty:
            continue
        try:
            c = geom.centroid
            cx[i] = float(c.x); cy[i] = float(c.y); valid[i] = True
        except Exception:
            pass
    if not valid.any():
        return pd.Series([""] * n, dtype=str)

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return pd.Series([""] * n, dtype=str)

    tree = cKDTree(np.column_stack([ssam_df["x_um"].astype(float).values,
                                       ssam_df["y_um"].astype(float).values]))
    pts_xy = np.column_stack([cx, cy])
    neighbors = tree.query_ball_point(pts_xy, r=radius_um)
    labels = ssam_df[label_col].astype(str).values

    out = [""] * n
    for i in range(n):
        if not valid[i]:
            continue
        idxs = neighbors[i]
        if not idxs:
            continue
        ctr = Counter(labels[idxs])
        out[i] = ctr.most_common(1)[0][0]
    return pd.Series(out, name="ssam_celltype")


# ---------------------------------------------------------------------------
# Grid majority for heatmap
# ---------------------------------------------------------------------------
def grid_majority(df: pd.DataFrame, x_col: str, y_col: str, label_col: str,
                   bbox_um: Tuple[float, float, float, float],
                   grid_um: float = 30.0) -> Tuple[np.ndarray, List[str], Tuple[float, float, float, float]]:
    """Bin (x, y, label) into a grid_um × grid_um grid and return:

    * ``grid`` — int ndarray (H, W); ``-1`` means empty cell, otherwise index
      into the returned ``labels`` list of the majority label.
    * ``labels`` — list of unique label strings.
    * ``extent_um`` — (x0, x1, y1, y0) suitable for ``imshow(extent=...)`` with
      origin top-left (matching the rest of the package).
    """
    x0, y0, x1, y1 = bbox_um
    W = max(1, int(np.ceil((x1 - x0) / grid_um)))
    H = max(1, int(np.ceil((y1 - y0) / grid_um)))
    grid = np.full((H, W), -1, dtype=np.int32)
    labels: List[str] = []
    if df is None or len(df) == 0 or label_col not in df.columns:
        return grid, labels, (x0, x1, y1, y0)

    label_to_idx: Dict[str, int] = {}
    def _li(lbl: str) -> int:
        idx = label_to_idx.get(lbl)
        if idx is None:
            idx = len(labels)
            labels.append(lbl)
            label_to_idx[lbl] = idx
        return idx

    x = df[x_col].astype(float).values
    y = df[y_col].astype(float).values
    lab = df[label_col].astype(str).values

    cols = ((x - x0) / grid_um).astype(int)
    rows = ((y - y0) / grid_um).astype(int)
    in_box = (cols >= 0) & (cols < W) & (rows >= 0) & (rows < H) & (lab != "")
    cols = cols[in_box]; rows = rows[in_box]; lab = lab[in_box]
    if len(lab) == 0:
        return grid, labels, (x0, x1, y1, y0)

    # Per-cell Counter
    from collections import defaultdict
    bag: Dict[Tuple[int, int], Counter] = defaultdict(Counter)
    for r, c, l in zip(rows, cols, lab):
        bag[(int(r), int(c))][l] += 1
    for (r, c), ctr in bag.items():
        top_lbl, _ = ctr.most_common(1)[0]
        grid[r, c] = _li(top_lbl)

    return grid, labels, (x0, x1, y1, y0)


# ---------------------------------------------------------------------------
# Stable celltype color palette
# ---------------------------------------------------------------------------
_FALLBACK_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
]


def stable_celltype_palette(celltypes: Iterable[str]) -> Dict[str, str]:
    """Deterministic celltype string → hex color (sorted alphabetically)."""
    uniq = sorted({str(c) for c in celltypes if c is not None and str(c)})
    return {c: _FALLBACK_PALETTE[i % len(_FALLBACK_PALETTE)] for i, c in enumerate(uniq)}
