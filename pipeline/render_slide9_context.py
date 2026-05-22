"""Render a context overview that shows where the slide9 zoom region sits
inside the parent ROI ovrlpy_low_coherence_1.

Output: parent 500x500 µm DAPI + Baysor centroids + pink VSI<0.5 mask
+ yellow rectangle marking the slide9 zoom bbox.
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
from pipeline.roi_viz.cropping import crop_dapi, crop_vsi
from pipeline.roi_viz.data_classes import ROIRecord
from pipeline.roi_viz.loaders import load_method
from pipeline.roi_viz.pipeline import build_sample_data


PARENT_BBOX = (5114.0, 6104.0, 5614.0, 6604.0)
ZOOM_BOXES = {
    "A": ((5174.0, 6164.0, 5254.0, 6244.0), "#FFD700"),
}
OUT = Path("xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/"
           "roi_test_output/roi_visualization/ovrlpy_low_coherence_1/low_vsi_slide9/"
           "slide9_zoom_context.png")


def main():
    with open("pipeline/config.yaml") as f:
        cfg = yaml.safe_load(f)
    sample = build_sample_data(cfg)
    outdir = Path("xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/"
                  "roi_test_output/roi_visualization")
    cache.configure(root=outdir, force_stages=None, force_rois=None, version="0.2.0")

    parent = ROIRecord(
        roi_id="ovrlpy_low_coherence_1", sample_id=sample.sample_id,
        roi_class="ovrlpy", bbox_um=PARENT_BBOX,
        source_map="ovrlpy", hypothesis="resegmentation_robustness",
        primary_reference=None,
    )

    dapi = crop_dapi(sample, parent)
    vsi  = crop_vsi(sample, parent)
    res  = load_method("baysor", sample, parent)

    fig, ax = plt.subplots(figsize=(8.0, 8.0), facecolor="black", dpi=200)
    ax.set_facecolor("black")
    _draw_dapi(ax, dapi, parent, alpha=0.55)

    # pink VSI<0.5 overlay
    if vsi.get("map") is not None:
        arr = vsi["map"]; ext = vsi["extent_um"]
        low = (arr < 0.5) & np.isfinite(arr)
        dapi_img = dapi.get("image")
        if dapi_img is not None and dapi_img.size:
            from skimage.transform import resize as _resize
            from skimage.filters import threshold_otsu
            from skimage.morphology import binary_closing, disk
            from scipy.ndimage import binary_dilation
            dapi_rs = _resize(dapi_img, arr.shape, anti_aliasing=True, preserve_range=True)
            thr = threshold_otsu(dapi_rs) if dapi_rs.max() > dapi_rs.min() else dapi_rs.min()
            tissue = binary_closing(dapi_rs > thr, disk(2))
            low = binary_dilation(low & tissue, iterations=3) & tissue
        if low.any():
            rgba = np.zeros((*arr.shape, 4), dtype=np.float32)
            rgba[low, 0] = 1.0; rgba[low, 1] = 0.05; rgba[low, 2] = 0.65
            rgba[low, 3] = 0.55
            ax.imshow(rgba, extent=(ext[0], ext[1], ext[2], ext[3]),
                      origin="upper", interpolation="nearest", zorder=3)

    # baysor centroids
    c = res["centroids"]
    if c is not None and len(c):
        ax.scatter(c["x_centroid_um"], c["y_centroid_um"], s=4,
                   c=st.METHOD_COLORS["baysor"], alpha=0.85, zorder=5)

    # zoom rectangles
    from matplotlib.patches import Rectangle
    for tag, (bbox, col) in ZOOM_BOXES.items():
        zx0, zy0, zx1, zy1 = bbox
        ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0,
                                  fill=False, edgecolor=col,
                                  linewidth=3.0, zorder=20))
        ax.text(zx0 + 5, zy0 + 14, f"Zoom {tag}",
                color=col, fontsize=11, fontweight="bold", zorder=21)

    x0, y0, x1, y1 = PARENT_BBOX
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_color("white")
    ax.set_title("ovrlpy_low_coherence_1 (parent 500x500 µm)  -  Baysor centroids + low-VSI mask",
                 color="white", fontsize=11, pad=10)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT), dpi=200, bbox_inches=None, facecolor="black")
    plt.close(fig)
    print("written:", OUT)


if __name__ == "__main__":
    main()
