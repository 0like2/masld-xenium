"""Render 3 PNGs for PPT slide 9 panel 'a. Low z-axis signal coherence'.

Produces zoomed-in figures (200 x 200 µm) inside ROI ovrlpy_low_coherence_1
for three methods: 10x default, Cellpose only, Cellpose + Baysor.

Each PNG (1600x1600 square):
- DAPI background (dim)
- Pink VSI < 0.5 mask (where low-VSI tissue is) — uses Otsu DAPI gate to
  exclude off-tissue pixels
- Transcripts colored by z (coolwarm)
- Method cell boundary contours
- Title / legend use PPT-consistent wording: "Low z-axis signal coherence"
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
os.chdir(Path(__file__).resolve().parent.parent)

from pipeline.roi_viz import cache
from pipeline.roi_viz import style as st
from pipeline.roi_viz.common_renders import _draw_dapi
from pipeline.roi_viz.cropping import crop_dapi, crop_transcripts, crop_vsi
from pipeline.roi_viz.data_classes import ROIRecord
from pipeline.roi_viz.loaders import load_method
from pipeline.roi_viz.pipeline import build_sample_data


METHODS_FOR_SLIDE = [
    ("xenium_nucleus",     "10x default"),
    ("optimal_expansion",  "Cellpose + expansion"),
    ("baysor",             "Cellpose + Baysor"),
]

ZOOM_OPTIONS = {
    "optA": (5174.0, 6164.0, 5254.0, 6244.0),  # 80x80 µm centered on (5214, 6204)
}
PARENT_BBOX  = (5114.0, 6104.0, 5614.0, 6604.0)
OUT_DIR = Path("xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/"
               "roi_test_output/roi_visualization/ovrlpy_low_coherence_1/low_vsi_slide9")


def _draw_polygons_outline(ax, polys, color, lw, alpha, bbox_um):
    """Draw polygon outlines clipped to bbox."""
    if polys is None or not len(polys):
        return
    from matplotlib.path import Path as MPath
    from matplotlib.collections import PathCollection
    x0, y0, x1, y1 = bbox_um
    paths = []
    for g in polys["geometry"]:
        if g is None or g.is_empty or g.geom_type != "Polygon":
            continue
        xs, ys = zip(*list(g.exterior.coords))
        # quick bbox cull
        if max(xs) < x0 or min(xs) > x1 or max(ys) < y0 or min(ys) > y1:
            continue
        verts = np.column_stack([xs, ys])
        codes = [MPath.MOVETO] + [MPath.LINETO] * (len(verts) - 1)
        paths.append(MPath(verts, codes))
    if not paths:
        return
    pc = PathCollection(paths, facecolors="none", edgecolors=color,
                        linewidths=lw, alpha=alpha, zorder=6)
    ax.add_collection(pc)


def _low_vsi_overlay(ax, vsi, dapi, bbox_um):
    """Pink VSI<0.5 mask, gated by Otsu DAPI to exclude background."""
    if vsi is None or vsi.get("map") is None:
        return False
    arr = vsi["map"]
    ext = vsi["extent_um"]  # (x0, x1, y0, y1)  origin upper → y0>y1
    low = (arr < 0.5) & np.isfinite(arr)
    dapi_img = dapi.get("image") if dapi else None
    if dapi_img is not None and dapi_img.size:
        from skimage.transform import resize as _resize
        dapi_rs = _resize(dapi_img, arr.shape, anti_aliasing=True, preserve_range=True)
        try:
            from skimage.filters import threshold_otsu
            from skimage.morphology import binary_closing, disk
            if float(dapi_rs.max()) > float(dapi_rs.min()):
                thr = threshold_otsu(dapi_rs)
            else:
                thr = float(dapi_rs.min())
            tissue = dapi_rs > thr
            tissue = binary_closing(tissue, disk(2))
        except Exception:
            tissue = dapi_rs > float(np.percentile(dapi_rs, 50))
        low = low & tissue
        try:
            from scipy.ndimage import binary_dilation
            low = binary_dilation(low, iterations=3) & tissue
        except Exception:
            pass
    if not low.any():
        return False
    rgba = np.zeros((*arr.shape, 4), dtype=np.float32)
    rgba[low, 0] = 1.0
    rgba[low, 1] = 0.05
    rgba[low, 2] = 0.65
    rgba[low, 3] = 0.55
    ax.imshow(rgba, extent=(ext[0], ext[1], ext[2], ext[3]),
              origin="upper", interpolation="nearest", zorder=3)
    return True


def render_one(sample, method_key, method_label, zoom_roi, parent_roi, out_path):
    res = load_method(method_key, sample, parent_roi)
    dapi = crop_dapi(sample, zoom_roi)
    tx = crop_transcripts(sample, zoom_roi)
    vsi = crop_vsi(sample, zoom_roi)

    fig, ax = plt.subplots(figsize=(8.0, 8.0), facecolor="black", dpi=200)
    ax.set_facecolor("black")

    _draw_dapi(ax, dapi, zoom_roi, alpha=0.55)

    low_drawn = _low_vsi_overlay(ax, vsi, dapi, zoom_roi.bbox_um)

    if tx is not None and len(tx):
        z_vals = tx["z_location"].values
        z_min = float(np.nanpercentile(z_vals, 1))
        z_max = float(np.nanpercentile(z_vals, 99))
        sc = ax.scatter(tx["x_location"], tx["y_location"],
                        c=z_vals, cmap="coolwarm",
                        vmin=z_min, vmax=z_max,
                        s=8.0, alpha=0.95, rasterized=True, edgecolors="none",
                        zorder=4)
        # compact inset colorbar (lower-right, big enough to read after shrink)
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(ax, width="4%", height="28%", loc="lower right",
                         bbox_to_anchor=(-0.02, 0.04, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(sc, cax=cax)
        cb.set_label("z (µm)", fontsize=15, color="white", fontweight="bold")
        cb.ax.tick_params(colors="white", labelsize=12)
        cb.outline.set_edgecolor("white")
        cb.outline.set_linewidth(1.4)

    polys = res.get("polygons") if res.get("available") else None
    method_color = st.METHOD_COLORS.get(method_key, "#FFD700")
    _draw_polygons_outline(ax, polys, color=method_color, lw=3.0, alpha=0.95,
                           bbox_um=zoom_roi.bbox_um)

    handles = []
    if low_drawn:
        from matplotlib.patches import Patch
        handles.append(Patch(facecolor=(1.0, 0.05, 0.65, 0.55),
                             edgecolor="none",
                             label="Low z-axis coherence"))
    handles.append(plt.Line2D([], [], color=method_color, lw=4.0,
                              label=method_label))
    leg = ax.legend(handles=handles, loc="upper right", fontsize=14,
                    framealpha=0.85, facecolor="black", edgecolor="white",
                    labelcolor="white", handlelength=1.4, handletextpad=0.6,
                    borderpad=0.4)
    leg.set_zorder(50)
    leg.get_frame().set_linewidth(1.4)
    for txt in leg.get_texts():
        txt.set_color("white")

    ax.set_title(f"{method_label}",
                 fontsize=24, color="white", pad=10, fontweight="bold")
    x0, y0, x1, y1 = zoom_roi.bbox_um
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("white")

    # scale bar (20 µm)
    bar_um = 20.0
    bx0 = x0 + (x1 - x0) * 0.05
    bx1 = bx0 + bar_um
    by  = y1 - (y1 - y0) * 0.05  # near bottom (remember y0>y1 inverted... actually y0<y1 in bbox)
    by  = y0 + (y1 - y0) * 0.92
    ax.plot([bx0, bx1], [by, by], color="white", lw=5.0, zorder=20)
    ax.text((bx0 + bx1) / 2, by - (y1 - y0) * 0.020, f"{int(bar_um)} µm",
            color="white", ha="center", va="bottom", fontsize=15,
            fontweight="bold", zorder=21)

    fig.tight_layout()
    fig.savefig(str(out_path), dpi=200, bbox_inches=None, pad_inches=0.05,
                facecolor="black", edgecolor="none")
    plt.close(fig)
    return out_path


def main():
    with open("pipeline/config.yaml") as f:
        cfg = yaml.safe_load(f)
    sample = build_sample_data(cfg)
    cache.configure(root=Path(OUT_DIR).parent.parent.parent / "roi_visualization",
                    force_stages=None, force_rois=None, version="0.2.0")

    parent_roi = ROIRecord(
        roi_id="ovrlpy_low_coherence_1", sample_id=sample.sample_id,
        roi_class="ovrlpy", bbox_um=PARENT_BBOX,
        source_map="ovrlpy", hypothesis="resegmentation_robustness",
        primary_reference=None,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for opt_tag, bbox in ZOOM_OPTIONS.items():
        zoom_roi = ROIRecord(
            roi_id=f"ovrlpy_low_coherence_1_zoom_{opt_tag}",
            sample_id=sample.sample_id, roi_class="ovrlpy", bbox_um=bbox,
            source_map="ovrlpy", hypothesis="resegmentation_robustness",
            primary_reference=None,
        )
        for key, label in METHODS_FOR_SLIDE:
            # legacy filename for option A (no suffix) so existing PPT refs stay valid
            if opt_tag == "optA":
                out = OUT_DIR / f"slide9_low_z_axis_coherence_{key}.png"
            else:
                out = OUT_DIR / f"slide9_low_z_axis_coherence_{key}_{opt_tag}.png"
            print(f"[{opt_tag}] Rendering {key} → {out}")
            render_one(sample, key, label, zoom_roi, parent_roi, out)

    print("Done. Files written under:", OUT_DIR)


if __name__ == "__main__":
    main()
