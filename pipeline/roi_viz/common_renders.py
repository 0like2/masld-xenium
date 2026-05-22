"""Common (non-class-specific) ROI rendering — C1–C4 + M1–M4.

All renderers:
- Build white-background figures (no transparent saves) per style §10.
- Use pre-cached ROI artifacts (DAPI / transcripts / per-method polygons).
- Skip work if the output PNG is newer than the cached artifacts (mtime cache).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon as MplPolygon, Rectangle

from . import cache
from . import style as st
from .cropping import crop_dapi, crop_transcripts, crop_p2r, crop_ssam, crop_vsi
from .data_classes import ROIRecord, SampleData
from .loaders import load_method
from .markers import marker_color_palette, select_markers

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dynamic scale bar (length = N × baseline cell diameter)
# ---------------------------------------------------------------------------
_BASELINE_AREA_CACHE: Dict[str, float] = {}


def _baseline_mean_cell_area(sample: SampleData, roi: ROIRecord,
                              baseline: str = "xenium_nucleus"):
    """Return mean polygon area (µm²) of *baseline* method on this ROI.

    Cached in-memory per (roi, baseline) so multiple renders for the same
    ROI don't repeat the geometry .area calculation."""
    key = f"{roi.roi_id}|{baseline}"
    if key in _BASELINE_AREA_CACHE:
        return _BASELINE_AREA_CACHE[key]
    try:
        res = load_method(baseline, sample, roi)
        if not res["available"] or not len(res["polygons"]):
            return None
        polys = res["polygons"]
        if "geometry" not in polys.columns:
            return None
        areas = [g.area for g in polys["geometry"]
                 if g is not None and not g.is_empty]
        if not areas:
            return None
        mean_a = float(np.nanmean(areas))
        _BASELINE_AREA_CACHE[key] = mean_a
        return mean_a
    except Exception:
        return None


def _add_dynamic_scale_bar(ax, roi: ROIRecord, sample: SampleData,
                            baseline: str = "xenium_nucleus",
                            n_cells: int = 5,
                            color: str = "white",
                            bbox_um: Optional[Any] = None) -> None:
    """Add a scale bar whose length equals *n_cells* × baseline cell diameter.

    Label format: ``~N cells (XXX µm)``. Falls back to fixed 50 µm when
    baseline cell area is unavailable.

    *bbox_um* defaults to ``roi.bbox_um``; pass an alternate bbox for the
    'wider context' renders that crop a larger area than the ROI itself.
    """
    if bbox_um is None:
        bbox_um = roi.bbox_um

    mean_area = _baseline_mean_cell_area(sample, roi, baseline)
    if mean_area is None or mean_area <= 0:
        st.add_scale_bar(ax, length_um=50, bbox_um=bbox_um, color=color)
        return

    diam = 2.0 * np.sqrt(mean_area / np.pi)
    raw_len = n_cells * diam
    # Round to a nice number
    if raw_len < 7:
        nice = 5
    elif raw_len < 15:
        nice = 10
    elif raw_len < 30:
        nice = 20
    elif raw_len < 75:
        nice = 50
    elif raw_len < 150:
        nice = 100
    else:
        nice = round(raw_len / 50) * 50

    actual_n = nice / diam
    _add_scale_bar_with_cell_label(ax, length_um=nice, bbox_um=bbox_um,
                                     n_cells=int(round(actual_n)), color=color)


def _add_scale_bar_with_cell_label(ax, length_um, bbox_um, n_cells, color):
    """style.add_scale_bar variant with 'N cells (X µm)' label."""
    from matplotlib import patheffects as mpe
    x0, y0, x1, y1 = bbox_um
    span_x = x1 - x0
    span_y = abs(y1 - y0)
    margin_x = span_x * 0.04
    margin_y = span_y * 0.04
    x_end = x1 - margin_x
    x_start = x_end - length_um
    y_pos = y1 - margin_y
    thickness = 2.0
    ax.plot([x_start, x_end], [y_pos, y_pos],
            color=color, linewidth=thickness, solid_capstyle="butt",
            path_effects=[mpe.withStroke(linewidth=thickness + 1.5,
                                            foreground="black")])
    label = f"~{int(n_cells)} cells ({int(length_um)} µm)"
    ax.text((x_start + x_end) / 2, y_pos - span_y * 0.02,
             label, ha="center", va="bottom", color=color, fontsize=7,
             path_effects=[mpe.withStroke(linewidth=2.0, foreground="black")])


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _draw_dapi(ax, dapi_dict: Dict[str, Any], roi: ROIRecord,
               alpha: float = 1.0, vmax_pct: float = 99.5) -> None:
    """Draw DAPI on *ax*. Always uses a black axes face so any alpha < 1
    blends DAPI signal toward black (preserving microscopy aesthetic) instead
    of toward the figure's white face (which would wash out transcript dots
    and boundary outlines on top)."""
    ax.set_facecolor("black")
    img = dapi_dict.get("image") if dapi_dict else None
    x0, y0, x1, y1 = roi.bbox_um
    if img is None or img.size == 0:
        ax.text(0.5, 0.5, "DAPI unavailable", transform=ax.transAxes,
                ha="center", va="center", color="#FFCC00")
        return
    vmin = float(np.percentile(img, 1))
    vmax = float(np.percentile(img, vmax_pct))
    if vmax <= vmin:
        vmax = vmin + 1.0
    ax.imshow(img, cmap=st.cmap_with_white_bad("gray"),
              vmin=vmin, vmax=vmax, alpha=alpha,
              extent=(x0, x1, y1, y0),
              interpolation="nearest")


def _draw_polygons_outline(ax, polygons,
                            color: str = "#00FFFF",
                            lw: float = 1.0,
                            alpha: float = 0.9,
                            halo: bool = True) -> None:
    """Draw polygon outlines with an optional white halo for visibility on
    dark backgrounds."""
    if polygons is None or len(polygons) == 0:
        return
    geoms = polygons["geometry"] if "geometry" in polygons.columns else None
    if geoms is None:
        return
    segments = []
    for g in geoms:
        if g is None or g.is_empty:
            continue
        if g.geom_type == "Polygon":
            xs, ys = zip(*list(g.exterior.coords))
            segments.append(np.column_stack([xs, ys]))
        elif g.geom_type == "MultiPolygon":
            for sub in g.geoms:
                xs, ys = zip(*list(sub.exterior.coords))
                segments.append(np.column_stack([xs, ys]))
    if not segments:
        return
    lc = LineCollection(segments, colors=color, linewidths=lw, alpha=alpha)
    if halo:
        lc.set_path_effects(st.boundary_path_effects(lw + 1.5))
    ax.add_collection(lc)


def _scatter_transcripts(ax, df: pd.DataFrame,
                          c: str = "#FFCC00", s: float = 0.4,
                          alpha: float = 0.6,
                          edgecolor: Optional[str] = None,
                          edge_lw: float = 0.0,
                          rasterized: bool = True) -> None:
    if df is None or len(df) == 0:
        return
    ax.scatter(df["x_location"], df["y_location"], s=s, c=c,
               alpha=alpha, edgecolors=edgecolor, linewidths=edge_lw,
               rasterized=rasterized)


# ---------------------------------------------------------------------------
# Output path helper
# ---------------------------------------------------------------------------
def _common_dir(outdir: Path, roi_id: str) -> Path:
    d = outdir / roi_id / "common"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _png_skip(out_path: Path, dep_paths: Sequence[Path]) -> bool:
    if cache.png_is_fresh(out_path, dep_paths):
        logger.debug(f"PNG fresh, skipping: {out_path.name}")
        return True
    return False


# ---------------------------------------------------------------------------
# C1: Context with ROI box
# ---------------------------------------------------------------------------
def render_context_with_box(sample: SampleData, roi: ROIRecord,
                              outdir: Path,
                              context_path: Optional[Path] = None) -> Path:
    """C1 — wide tissue overview (or step3 DAPI overview PNG) with ROI bbox."""
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_context.png"
    if context_path is None:
        # Reuse pre-rendered overview if available.
        candidates = [
            sample.step3_dir / f"{sample.sample_tag}_step3_dapi_overview.png" if sample.step3_dir else None,
            sample.output_dir / sample.input_path.name / "roi_test_output" / f"{sample.sample_tag}_step2_roi_overview_map.png",
        ]
        for c in candidates:
            if c and Path(c).exists():
                context_path = Path(c)
                break

    if _png_skip(out, [Path(context_path) if context_path else Path("")]):
        return out

    fig, ax = st.new_figure(figsize=(5, 5), dpi=180)
    if context_path and Path(context_path).exists():
        try:
            img = plt.imread(context_path)
            ax.imshow(img)
            # Note: we can't superimpose the ROI bbox without the image's
            # coordinate system. Annotate as text instead.
            ax.text(0.02, 0.98,
                    f"ROI: {roi.roi_id}\nclass: {roi.roi_class}\n"
                    f"bbox=({roi.bbox_um[0]:.0f},{roi.bbox_um[1]:.0f})-"
                    f"({roi.bbox_um[2]:.0f},{roi.bbox_um[3]:.0f}) µm",
                    transform=ax.transAxes, va="top", ha="left",
                    color="white", fontsize=7,
                    bbox=dict(facecolor="#1A1A1A", edgecolor="#444444",
                              alpha=0.85, pad=4))
        except Exception as e:
            logger.warning(f"context image read failed: {e}")
    else:
        ax.text(0.5, 0.5, "Context image\nunavailable",
                transform=ax.transAxes, ha="center", va="center")
    ax.axis("off")
    st.save_figure(fig, str(out), dpi=180)
    return out


# ---------------------------------------------------------------------------
# C2: DAPI crop
# ---------------------------------------------------------------------------
def render_dapi_crop(sample: SampleData, roi: ROIRecord, outdir: Path) -> Path:
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_dapi.png"
    if _png_skip(out, [sample.dapi_path] if sample.dapi_path else []):
        return out

    dapi = crop_dapi(sample, roi)
    fig, ax = st.new_figure(figsize=(5, 5))
    _draw_dapi(ax, dapi, roi)
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=f"{roi.roi_id} | DAPI")
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# C3: Reference crop (SSAM/P2R/VSI)
# ---------------------------------------------------------------------------
def render_reference_crop(sample: SampleData, roi: ROIRecord, outdir: Path) -> Optional[Path]:
    ref = (roi.primary_reference or "").lower()
    if not ref:
        # Pick by ROI class
        ref = {"architecture": "ssam", "compartment": "p2r",
               "low_vsi": "vsi", "fold_boundary": "vsi",
               "high_density": "p2r", "easy_control": "ssam"}.get(roi.roi_class, "ssam")

    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_reference_{ref}.png"
    if ref == "ssam":
        deps = [sample.ssam_h5ad_path] if sample.ssam_h5ad_path else []
    elif ref == "p2r":
        deps = [sample.p2r_csv_path] if sample.p2r_csv_path else []
    elif ref == "vsi":
        deps = [sample.vsi_map_path] if sample.vsi_map_path else []
    else:
        deps = []
    if _png_skip(out, deps):
        return out

    fig, ax = st.new_figure(figsize=(5, 5))
    title_suffix = ref.upper()
    drew = False

    if ref == "ssam":
        df = crop_ssam(sample, roi)
        if len(df):
            color_col = next((c for c in ("celltype", "ssam_celltype",
                                             "leiden_assignment", "leiden")
                              if c in df.columns), None)
            colors = df[color_col].astype(str) if color_col else "#888888"
            cats = pd.Categorical(colors).codes if color_col else None
            ax.scatter(df["x_um"], df["y_um"], c=cats if cats is not None else "#888",
                       cmap="tab20", s=4, alpha=0.85)
            drew = True

    elif ref == "p2r":
        df = crop_p2r(sample, roi)
        if len(df):
            cmp_col = next((c for c in ("p2r_compartment", "compartment",
                                          "p2r_celltype", "celltype",
                                          "p2r_name")
                            if c in df.columns), None)
            cats = pd.Categorical(df[cmp_col].astype(str)).codes if cmp_col else None
            x_col = next((c for c in ("x_um", "x_location", "x") if c in df.columns), None)
            y_col = next((c for c in ("y_um", "y_location", "y") if c in df.columns), None)
            if x_col and y_col:
                ax.scatter(df[x_col], df[y_col],
                           c=cats if cats is not None else "#3478F6",
                           cmap="tab10", s=2, alpha=0.6)
                drew = True

    elif ref == "vsi":
        d = crop_vsi(sample, roi)
        if d.get("map") is not None:
            ext = d["extent_um"]
            ax.imshow(d["map"], cmap="viridis", origin="upper",
                      extent=(ext[0], ext[1], ext[2], ext[3]),
                      vmin=0, vmax=1)
            drew = True

    if not drew:
        ax.text(0.5, 0.5, f"{title_suffix} reference\nunavailable",
                transform=ax.transAxes, ha="center", va="center",
                color="#888888")

    st.setup_axes(ax, bbox_um=roi.bbox_um, title=f"{roi.roi_id} | {title_suffix}")
    if drew:
        _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# C4: Metric card  (placeholder; populated by metrics.py later)
# ---------------------------------------------------------------------------
def render_metric_card(roi: ROIRecord, outdir: Path,
                         metrics_table: Optional[pd.DataFrame] = None) -> Path:
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_metric_card.png"
    fig, ax = st.new_figure(figsize=(6, 3.2))
    ax.axis("off")
    if metrics_table is None or metrics_table.empty:
        ax.text(0.5, 0.5,
                f"Metrics for {roi.roi_id}\n(populated after metrics stage)",
                ha="center", va="center", color="#666666")
    else:
        # Pivot a few key metrics × method into a small table.
        keep = ["method", "n_cells", "transcripts_per_cell", "genes_per_cell",
                "assigned_fraction", "mean_cell_area"]
        cols = [c for c in keep if c in metrics_table.columns]
        tbl = metrics_table[cols].copy()
        for c in tbl.columns:
            if tbl[c].dtype.kind == "f":
                tbl[c] = tbl[c].round(2)
        the_table = ax.table(cellText=tbl.values, colLabels=tbl.columns,
                              loc="center", cellLoc="center")
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1.0, 1.4)
        ax.set_title(f"{roi.roi_id} | metric card")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# M1: DAPI + boundary
# ---------------------------------------------------------------------------
def render_dapi_boundary(sample: SampleData, roi: ROIRecord, method: str,
                          outdir: Path) -> Path:
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_dapi_boundary_{method}.png"
    deps = []
    if sample.dapi_path: deps.append(sample.dapi_path)
    if _png_skip(out, deps):
        return out

    dapi = crop_dapi(sample, roi)
    res = load_method(method, sample, roi)

    fig, ax = st.new_figure(figsize=(5, 5))
    _draw_dapi(ax, dapi, roi, alpha=1.0)
    color = st.METHOD_COLORS.get(method, "#00FFFF")
    if res["available"]:
        _draw_polygons_outline(ax, res["polygons"], color=color, lw=0.9, alpha=0.95)
    else:
        ax.text(0.5, 0.05,
                f"{method} unavailable for this ROI",
                transform=ax.transAxes, ha="center", va="bottom",
                color="#FFCC00", fontsize=8,
                bbox=dict(facecolor="#1A1A1A", edgecolor="none", alpha=0.7, pad=3))
    title = f"{roi.roi_id} | DAPI + {st.METHOD_DISPLAY.get(method, method)}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# M2: All transcripts + boundary
# ---------------------------------------------------------------------------
def render_alltx_boundary(sample: SampleData, roi: ROIRecord, method: str,
                            outdir: Path) -> Path:
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_alltx_boundary_{method}.png"
    deps = []
    if sample.dapi_path: deps.append(sample.dapi_path)
    if sample.transcripts_path: deps.append(sample.transcripts_path)
    if _png_skip(out, deps):
        return out

    dapi = crop_dapi(sample, roi)
    tx = crop_transcripts(sample, roi)
    res = load_method(method, sample, roi)

    fig, ax = st.new_figure(figsize=(5, 5))
    _draw_dapi(ax, dapi, roi, alpha=0.6)
    if len(tx):
        _scatter_transcripts(ax, tx, c="#FFCC00", s=0.18, alpha=0.55,
                              edgecolor=None)
    color = st.METHOD_COLORS.get(method, "#00FFFF")
    if res["available"]:
        _draw_polygons_outline(ax, res["polygons"], color=color, lw=0.8, alpha=0.95)
    title = f"{roi.roi_id} | all tx + {st.METHOD_DISPLAY.get(method, method)} ({len(tx):,} tx)"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# M3: Marker transcripts + boundary
# ---------------------------------------------------------------------------
def render_markers_boundary(sample: SampleData, roi: ROIRecord, method: str,
                              outdir: Path,
                              markers: Optional[Sequence[str]] = None) -> Path:
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_markers_boundary_{method}.png"
    if _png_skip(out, [sample.dapi_path] if sample.dapi_path else []):
        return out

    dapi = crop_dapi(sample, roi)
    tx = crop_transcripts(sample, roi)
    res = load_method(method, sample, roi)

    if markers is None:
        mc = select_markers(roi, transcripts_in_roi=tx)
        markers = [g for g in (mc.target + mc.competing + mc.structural) if g]
    if not markers:
        markers = (tx["feature_name"].astype(str).value_counts()
                    .head(3).index.tolist()) if len(tx) else []

    palette = marker_color_palette(markers)

    fig, ax = st.new_figure(figsize=(5.4, 5))
    _draw_dapi(ax, dapi, roi, alpha=0.6)
    handles = []
    for gene in markers:
        sub = tx[tx["feature_name"].astype(str) == gene]
        if not len(sub):
            continue
        ax.scatter(sub["x_location"], sub["y_location"],
                   s=4.0, c=palette[gene], alpha=0.9,
                   edgecolors="white", linewidths=0.15,
                   rasterized=True, label=f"{gene} ({len(sub):,})")
        handles.append(gene)
    color = st.METHOD_COLORS.get(method, "#00FFFF")
    if res["available"]:
        _draw_polygons_outline(ax, res["polygons"], color=color, lw=0.7, alpha=0.85)
    if handles:
        ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
    title = f"{roi.roi_id} | markers ({', '.join(markers)}) + {st.METHOD_DISPLAY.get(method, method)}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# M4: Cell-type map (segmentation-based)
# ---------------------------------------------------------------------------
def render_celltype_map(sample: SampleData, roi: ROIRecord, method: str,
                          outdir: Path,
                          celltype_assignment: Optional[Dict[str, str]] = None) -> Path:
    """Polygon fill colored by cell type.

    *celltype_assignment* maps cell_id → celltype label. If None, falls back
    to coloring polygons by their raster_label (no celltype info available
    yet — populated when metrics stage produces cell-type tables).
    """
    out = _common_dir(outdir, roi.roi_id) / f"{roi.roi_id}_celltype_{method}.png"
    if _png_skip(out, [sample.dapi_path] if sample.dapi_path else []):
        return out

    dapi = crop_dapi(sample, roi)
    res = load_method(method, sample, roi)

    fig, ax = st.new_figure(figsize=(5.4, 5))
    _draw_dapi(ax, dapi, roi, alpha=0.4)

    if res["available"] and len(res["polygons"]):
        polys = res["polygons"]
        cell_ids = polys["cell_id"].astype(str).tolist()
        if celltype_assignment:
            ct = [celltype_assignment.get(cid, "unknown") for cid in cell_ids]
        else:
            ct = ["assigned"] * len(cell_ids)
        cats = pd.Categorical(ct)
        cmap = plt.get_cmap("tab20", max(2, len(cats.categories)))
        patches = []
        face_colors = []
        for poly, cat_code in zip(polys["geometry"].values, cats.codes):
            if poly is None or poly.is_empty:
                continue
            if poly.geom_type == "Polygon":
                xs, ys = zip(*list(poly.exterior.coords))
                patches.append(MplPolygon(np.column_stack([xs, ys]), closed=True))
                face_colors.append(cmap(cat_code))
            elif poly.geom_type == "MultiPolygon":
                for sub in poly.geoms:
                    xs, ys = zip(*list(sub.exterior.coords))
                    patches.append(MplPolygon(np.column_stack([xs, ys]), closed=True))
                    face_colors.append(cmap(cat_code))
        if patches:
            pc = PatchCollection(patches, facecolors=face_colors,
                                  edgecolors="black", linewidths=0.2,
                                  alpha=0.7)
            ax.add_collection(pc)
        # Legend (capped at 12 categories)
        cats_unique = list(cats.categories)
        if len(cats_unique) <= 12:
            for i, cat in enumerate(cats_unique):
                ax.scatter([], [], color=cmap(i), label=cat, s=20)
            ax.legend(loc="upper right", fontsize=7, framealpha=0.85)

    title = f"{roi.roi_id} | cell-type / {st.METHOD_DISPLAY.get(method, method)}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# Paper 4-bucket composite grid (rows = M1–M4, cols = methods)
# ---------------------------------------------------------------------------
def render_paper_4bucket_grid(sample: SampleData, roi: ROIRecord, outdir: Path,
                                methods: Sequence[str],
                                markers: Optional[Sequence[str]] = None) -> Path:
    """One composite figure with 4 rows (M1–M4) × N cols (methods).

    Each cell = small thumbnail of:
        row 0: DAPI + boundary (M1)
        row 1: all transcripts + boundary (M2)
        row 2: marker transcripts + boundary (M3)
        row 3: cell-type map (M4)

    Column header shows the paper bucket label `(N) ...` from
    ``style.METHOD_DISPLAY``.

    Useful for paper figure panels: lets a reviewer visually
    compare all 4–5 methods within a single ROI in one glance,
    instead of opening 20 separate PNGs.

    Reuses cached method results so repeated calls are near-zero cost.
    """
    out = outdir / roi.roi_id / f"{roi.roi_id}_paper_4bucket_grid.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    deps = []
    if sample.dapi_path: deps.append(sample.dapi_path)
    if sample.transcripts_path: deps.append(sample.transcripts_path)
    if _png_skip(out, deps):
        return out

    dapi = crop_dapi(sample, roi)
    tx = crop_transcripts(sample, roi)

    # Pick markers once (used for M3 row)
    if markers is None:
        from .markers import select_markers
        mc = select_markers(roi, transcripts_in_roi=tx)
        markers = [g for g in (mc.target + mc.competing + mc.structural) if g]
    if not markers and len(tx):
        markers = (tx["feature_name"].astype(str).value_counts()
                    .head(3).index.tolist())
    from .markers import marker_color_palette
    palette = marker_color_palette(markers)

    n_cols = len(methods)
    n_rows = 4
    # ~2 in per cell. Tight but readable.
    fig, axes = plt.subplots(n_rows, n_cols,
                                figsize=(2.0 * n_cols, 2.0 * n_rows + 0.5),
                                facecolor="white", dpi=160)
    if n_rows == 1: axes = axes[None, :]
    if n_cols == 1: axes = axes[:, None]

    row_labels = ["DAPI + boundary", "all-tx + boundary",
                    "markers + boundary", "celltype map"]

    for col_i, method in enumerate(methods):
        res = load_method(method, sample, roi)
        method_color = st.METHOD_COLORS.get(method, "#00FFFF")
        polys_avail = res["available"] and "geometry" in res["polygons"].columns \
            and len(res["polygons"]) > 0

        for row_i in range(n_rows):
            ax = axes[row_i, col_i]
            # Black axes facecolor → DAPI signal-less pixels stay black at any
            # alpha (instead of letting the figure's white face bleed through
            # and turning the panel grey). Required for transcript/marker dots
            # to keep contrast against a dark canvas across all rows.
            ax.set_facecolor("black")

            # Background: DAPI in rows 0-2, faint in row 3
            dapi_alpha = (1.0 if row_i == 0 else
                            0.6 if row_i in (1, 2) else 0.4)
            _draw_dapi(ax, dapi, roi, alpha=dapi_alpha)

            if row_i == 1 and len(tx):
                # all transcripts overlay
                _scatter_transcripts(ax, tx, c="#FFCC00", s=0.1,
                                       alpha=0.45, edgecolor=None)
            elif row_i == 2 and markers:
                # marker transcripts
                for gene in markers:
                    sub = tx[tx["feature_name"].astype(str) == gene]
                    if not len(sub):
                        continue
                    ax.scatter(sub["x_location"], sub["y_location"],
                                s=2.5, c=palette[gene], alpha=0.85,
                                edgecolors="white", linewidths=0.1,
                                rasterized=True)
            elif row_i == 3 and polys_avail:
                # filled cell-type map (single color = "assigned")
                from matplotlib.collections import PatchCollection
                from matplotlib.patches import Polygon as MplPolygon
                patches = []
                for g in res["polygons"]["geometry"]:
                    if g is None or g.is_empty:
                        continue
                    if g.geom_type == "Polygon":
                        xs, ys = zip(*list(g.exterior.coords))
                        patches.append(MplPolygon(np.column_stack([xs, ys]),
                                                     closed=True))
                if patches:
                    pc = PatchCollection(patches, facecolors=method_color,
                                            edgecolors="black", linewidths=0.15,
                                            alpha=0.55)
                    ax.add_collection(pc)

            # Boundary outline (rows 0-2)
            if row_i < 3 and polys_avail:
                _draw_polygons_outline(ax, res["polygons"],
                                         color=method_color, lw=0.5,
                                         alpha=0.95, halo=True)
            elif not polys_avail and row_i == 0:
                ax.text(0.5, 0.5, f"{method}\nunavailable",
                          transform=ax.transAxes,
                          ha="center", va="center", fontsize=7,
                          color="#FFCC00",
                          bbox=dict(facecolor="#1A1A1A", edgecolor="none",
                                      alpha=0.8, pad=2))

            st.setup_axes(ax, bbox_um=roi.bbox_um)

            # Column header (row 0 only)
            if row_i == 0:
                ax.set_title(st.METHOD_DISPLAY.get(method, method),
                                fontsize=8, pad=4)

    # Row labels via figure text (axes have axis off so set_ylabel hidden)
    for row_i, label in enumerate(row_labels):
        # vertical center of each row (approximate)
        y = 1.0 - (row_i + 0.5) / n_rows * 0.94 - 0.025
        fig.text(0.012, y, label, fontsize=9, fontweight="bold",
                  rotation=90, va="center", ha="center")

    # Single scale bar in bottom-right cell
    _add_dynamic_scale_bar(axes[-1, -1], roi, sample, color="white")

    fig.suptitle(
        f"{roi.roi_id} | Paper 4-bucket comparison "
        f"(class={roi.roi_class}, bbox=({roi.bbox_um[0]:.0f},{roi.bbox_um[1]:.0f})-"
        f"({roi.bbox_um[2]:.0f},{roi.bbox_um[3]:.0f}) µm)",
        fontsize=10, y=0.995)
    # Reserve left margin for row labels (~3% of width).
    fig.tight_layout(rect=(0.03, 0, 1, 0.97))
    st.save_figure(fig, str(out), dpi=160)
    return out


# ---------------------------------------------------------------------------
# Batch helper
# ---------------------------------------------------------------------------
def render_all_common(sample: SampleData, roi: ROIRecord, outdir: Path,
                        methods: Sequence[str],
                        celltype_assignment_by_method: Optional[Dict[str, Dict[str, str]]] = None,
                        ) -> Dict[str, Any]:
    """Render every C* and M* image for a single ROI; return path map."""
    paths: Dict[str, Any] = {}
    paths["context"] = render_context_with_box(sample, roi, outdir)
    paths["dapi"] = render_dapi_crop(sample, roi, outdir)
    paths["reference"] = render_reference_crop(sample, roi, outdir)
    paths["metric_card"] = render_metric_card(roi, outdir, metrics_table=None)

    paths["dapi_boundary"] = {}
    paths["alltx_boundary"] = {}
    paths["markers_boundary"] = {}
    paths["celltype"] = {}
    for m in methods:
        paths["dapi_boundary"][m] = render_dapi_boundary(sample, roi, m, outdir)
        paths["alltx_boundary"][m] = render_alltx_boundary(sample, roi, m, outdir)
        paths["markers_boundary"][m] = render_markers_boundary(sample, roi, m, outdir)
        ct_map = (celltype_assignment_by_method or {}).get(m)
        paths["celltype"][m] = render_celltype_map(sample, roi, m, outdir,
                                                     celltype_assignment=ct_map)

    # Paper 4-bucket composite — single figure summarizing all methods.
    paths["paper_grid"] = render_paper_4bucket_grid(sample, roi, outdir, methods)
    return paths
