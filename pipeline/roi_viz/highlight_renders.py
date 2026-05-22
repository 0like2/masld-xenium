"""★ Difference-highlight renderers — user-priority module.

For each (baseline, method) pair within an ROI we generate three PNGs:

* ``boundary_overlap`` — DAPI + baseline boundary (yellow) + method boundary
  (cyan).
* ``difference_mask`` — DAPI + a *pre-composited* RGBA where each pixel
  belongs to one of {baseline-only, method-only, shared, background}. The
  RGBA is built once and ``imshow``-ed once, which avoids matplotlib's
  alpha-blending darkening overlapping regions (plan §10.2).
* ``reassigned_transcripts`` — DAPI + transcripts colored by category
  {unchanged, reassigned, newly, dropped}, plus a small inset bar of the
  category fractions.

Heavy computation lives in ``diff_compute``; this module only draws.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import cache, style as st
from .common_renders import _draw_dapi, _draw_polygons_outline, _png_skip
from .cropping import crop_dapi
from .data_classes import ROIRecord, SampleData
from .diff_compute import (_method_source_paths, compose_difference_rgba,
                            compute_difference_mask, compute_reassignment,
                            REASSIGN_CATEGORIES)
from .loaders import load_method


def _highlight_deps(sample: SampleData, baseline: str, method: str):
    """Dependency paths whose mtime drives PNG cache invalidation for any
    highlight render that compares (baseline, method). Includes DAPI plus
    every source file consumed by either method."""
    deps = []
    if sample.dapi_path:
        deps.append(sample.dapi_path)
    deps.extend(_method_source_paths(sample, baseline))
    deps.extend(_method_source_paths(sample, method))
    return deps


# Helpers moved to common_renders for shared use across all render modules.
from .common_renders import _add_dynamic_scale_bar  # noqa: E402

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------
def _highlights_dir(outdir: Path, roi_id: str) -> Path:
    d = outdir / roi_id / "highlights"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _unavailable_text(ax, msg: str) -> None:
    ax.text(0.5, 0.5, msg, transform=ax.transAxes,
            ha="center", va="center", color="#FFCC00", fontsize=9,
            bbox=dict(facecolor="#1A1A1A", edgecolor="#444444",
                      alpha=0.85, pad=4))


# ---------------------------------------------------------------------------
# Boundary overlap (yellow vs cyan)
# ---------------------------------------------------------------------------
def render_boundary_overlap(sample: SampleData, roi: ROIRecord,
                              baseline: str, method: str,
                              outdir: Path) -> Path:
    out = _highlights_dir(outdir, roi.roi_id) / \
        f"{roi.roi_id}_boundary_overlap_{baseline}__vs__{method}.png"
    if _png_skip(out, _highlight_deps(sample, baseline, method)):
        return out

    dapi = crop_dapi(sample, roi)
    base = load_method(baseline, sample, roi)
    other = load_method(method, sample, roi)

    fig, ax = st.new_figure(figsize=(5.4, 5.4))
    _draw_dapi(ax, dapi, roi, alpha=0.9)

    if base["available"] and other["available"]:
        _draw_polygons_outline(ax, base["polygons"],
                                 color=st.HIGHLIGHT_COLORS["baseline_boundary"],
                                 lw=1.0, alpha=0.95)
        _draw_polygons_outline(ax, other["polygons"],
                                 color=st.HIGHLIGHT_COLORS["method_boundary"],
                                 lw=1.0, alpha=0.95)
    else:
        _unavailable_text(ax,
                           f"Boundary overlap unavailable\n"
                           f"baseline={baseline} (avail={base['available']}) | "
                           f"method={method} (avail={other['available']})")

    # Legend swatches
    handles = [
        plt.Line2D([], [], color=st.HIGHLIGHT_COLORS["baseline_boundary"],
                    lw=2.0, label=f"baseline: {st.METHOD_DISPLAY.get(baseline, baseline)}"),
        plt.Line2D([], [], color=st.HIGHLIGHT_COLORS["method_boundary"],
                    lw=2.0, label=f"method: {st.METHOD_DISPLAY.get(method, method)}"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.85)

    title = f"{roi.roi_id} | boundary overlap : {baseline} vs {method}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample)
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# Difference mask (single-imshow RGBA)
# ---------------------------------------------------------------------------
def render_difference_mask(sample: SampleData, roi: ROIRecord,
                             baseline: str, method: str,
                             outdir: Path) -> Path:
    out = _highlights_dir(outdir, roi.roi_id) / \
        f"{roi.roi_id}_difference_mask_{baseline}__vs__{method}.png"
    if _png_skip(out, _highlight_deps(sample, baseline, method)):
        return out

    dapi = crop_dapi(sample, roi)
    diff = compute_difference_mask(sample, roi, baseline, method)

    fig, ax = st.new_figure(figsize=(5.6, 5.4))
    _draw_dapi(ax, dapi, roi, alpha=0.55)

    if diff["available"]:
        rgba = compose_difference_rgba(diff["category"])
        x0, y0, x1, y1 = diff["extent_um"]
        # ax.imshow with extent flips y already to match data coords.
        ax.imshow(rgba, extent=(x0, x1, y1, y0), interpolation="nearest")

        # Overlay baseline cell outlines so individual nuclei stay visible
        # inside large method-only blobs (e.g. rigid_expansion that fills
        # >60 % of the ROI). Without this, the diff mask collapses into a
        # single blue blob and per-cell detail is lost.
        try:
            base_res = load_method(baseline, sample, roi)
            if base_res["available"]:
                _draw_polygons_outline(
                    ax, base_res["polygons"],
                    color=st.HIGHLIGHT_COLORS["baseline_boundary"],
                    lw=0.5, alpha=0.85, halo=False)
        except Exception:
            pass
        # Area summary inset (text box in lower-left).
        a = diff["area_summary"]
        # Coverage of baseline by method (= shared / (baseline_only + shared))
        # — meaningful when one mask is much smaller than the other (e.g.
        # tiny nucleus vs huge expansion). Values close to 1.0 mean the
        # baseline is fully contained within the method's mask.
        # Numerical area summary (areas, IoU, coverage) lives in the
        # companion difference_mask_summary_{roi_id}.png — separate per-ROI
        # bar chart that compares all 4 methods. Spatial image stays clean.
    else:
        _unavailable_text(ax,
                           f"Difference mask unavailable\n"
                           f"reason: {diff.get('reason')}")

    # Legend
    handles = [
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color=st.HIGHLIGHT_COLORS["xenium_only"],
                    label=f"{st.METHOD_DISPLAY.get(baseline, baseline)}-only"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color=st.HIGHLIGHT_COLORS["method_only"],
                    label=f"{st.METHOD_DISPLAY.get(method, method)}-only"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color=st.HIGHLIGHT_COLORS["shared"], label="shared"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.85)

    title = f"{roi.roi_id} | difference mask : {baseline} vs {method}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample)
    st.save_figure(fig, str(out))
    return out


# ---------------------------------------------------------------------------
# Transcript reassignment (4 categories + mini bar)
# ---------------------------------------------------------------------------
def render_reassigned_transcripts(sample: SampleData, roi: ROIRecord,
                                     baseline: str, method: str,
                                     outdir: Path) -> Path:
    out = _highlights_dir(outdir, roi.roi_id) / \
        f"{roi.roi_id}_reassigned_transcripts_{baseline}__vs__{method}.png"
    if _png_skip(out, _highlight_deps(sample, baseline, method)):
        return out

    dapi = crop_dapi(sample, roi)
    re = compute_reassignment(sample, roi, baseline, method)
    base = load_method(baseline, sample, roi)
    other = load_method(method, sample, roi)

    cat_color = {
        "unchanged":  st.HIGHLIGHT_COLORS["tx_unchanged"],
        "reassigned": st.HIGHLIGHT_COLORS["tx_reassigned"],
        "newly":      st.HIGHLIGHT_COLORS["tx_newly"],
        "dropped":    st.HIGHLIGHT_COLORS["tx_dropped"],
    }

    # ---- 1) Spatial PNG (square, no inset) -------------------------------
    fig, ax = st.new_figure(figsize=(5.6, 5.4))
    _draw_dapi(ax, dapi, roi, alpha=0.55)

    if re["available"] and len(re["frame"]):
        df = re["frame"]
        for cat in ("unchanged", "dropped", "newly", "reassigned"):
            sub = df[df["category"] == cat]
            if not len(sub):
                continue
            s = 0.25 if cat == "unchanged" else 0.6
            alpha = 0.35 if cat == "unchanged" else 0.85
            ax.scatter(sub["x_location"], sub["y_location"],
                        c=cat_color[cat], s=s, alpha=alpha,
                        edgecolors="white" if cat != "unchanged" else None,
                        linewidths=0.05 if cat != "unchanged" else 0,
                        rasterized=True, label=f"{cat} ({len(sub):,})")
        if base["available"]:
            _draw_polygons_outline(ax, base["polygons"],
                                     color=st.HIGHLIGHT_COLORS["baseline_boundary"],
                                     lw=0.4, alpha=0.45, halo=False)
        if other["available"]:
            _draw_polygons_outline(ax, other["polygons"],
                                     color=st.HIGHLIGHT_COLORS["method_boundary"],
                                     lw=0.4, alpha=0.45, halo=False)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
    else:
        _unavailable_text(ax,
                           f"Transcript reassignment unavailable\n"
                           f"reason: {re.get('reason')}")

    title = f"{roi.roi_id} | tx reassignment : {baseline} vs {method}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample)
    st.save_figure(fig, str(out))

    # ---- 2) Bar-only PNG (companion file) -------------------------------
    bar_out = _highlights_dir(outdir, roi.roi_id) / \
        f"{roi.roi_id}_reassignment_bar_{baseline}__vs__{method}.png"
    if not _png_skip(bar_out, _highlight_deps(sample, baseline, method)):
        fig_b, bar_ax = plt.subplots(figsize=(4.0, 3.2), facecolor="white", dpi=200)
        bar_ax.set_facecolor("white")
        if re["available"] and len(re["frame"]):
            cats = list(REASSIGN_CATEGORIES)
            fracs = [re["fractions"].get(c, 0.0) for c in cats]
            colors = [cat_color[c] for c in cats]
            bar_ax.bar(cats, fracs, color=colors, edgecolor="black", linewidth=0.5)
            bar_ax.set_ylim(0, max(0.05, max(fracs) * 1.18))
            bar_ax.tick_params(axis="x", labelrotation=15, labelsize=10)
            bar_ax.tick_params(axis="y", labelsize=9)
            bar_ax.set_ylabel("fraction of transcripts", fontsize=10)
            bar_ax.set_title(
                f"{roi.roi_id} | reassignment fractions\n"
                f"baseline={baseline}  vs  method={method}\n"
                f"total tx joined = {re['n_total']:,}",
                fontsize=10)
            for i, (cat, f) in enumerate(zip(cats, fracs)):
                n = int(f * re["n_total"])
                bar_ax.text(i, f + 0.01, f"{n:,}\n({f*100:.1f}%)",
                              ha="center", va="bottom", fontsize=8)
            for spine in bar_ax.spines.values():
                spine.set_color("#888")
                spine.set_linewidth(0.5)
        else:
            bar_ax.text(0.5, 0.5,
                          f"unavailable\nreason: {re.get('reason')}",
                          transform=bar_ax.transAxes, ha="center", va="center")
            bar_ax.axis("off")
        fig_b.tight_layout()
        fig_b.savefig(str(bar_out), dpi=200, bbox_inches="tight",
                       pad_inches=0.05, facecolor="white", edgecolor="none")
        plt.close(fig_b)
    return out


# ---------------------------------------------------------------------------
# Batch helper
# ---------------------------------------------------------------------------
def render_all_highlights(sample: SampleData, roi: ROIRecord,
                            baseline: str, methods: Iterable[str],
                            outdir: Path) -> Dict[str, Dict[str, Path]]:
    paths: Dict[str, Dict[str, Path]] = {
        "boundary_overlap": {},
        "difference_mask":  {},
        "reassignment":     {},
    }
    method_list = [m for m in methods if m != baseline]
    for m in method_list:
        paths["boundary_overlap"][m] = render_boundary_overlap(
            sample, roi, baseline, m, outdir)
        paths["difference_mask"][m] = render_difference_mask(
            sample, roi, baseline, m, outdir)
        paths["reassignment"][m] = render_reassigned_transcripts(
            sample, roi, baseline, m, outdir)

    # Companion summary bar chart (separate from spatial PNGs).
    summary_path = render_difference_mask_summary(
        sample, roi, baseline, method_list, outdir)
    if summary_path is not None:
        paths["difference_mask_summary"] = summary_path
    return paths


def render_difference_mask_summary(sample: SampleData, roi: ROIRecord,
                                       baseline: str,
                                       methods: Iterable[str],
                                       outdir: Path):
    """Per-ROI bar chart comparing baseline vs each method on:
       baseline-only %, method-only %, shared %, IoU.
    Replaces the inline inset that used to live in difference_mask images."""
    out = _highlights_dir(outdir, roi.roi_id) / \
        f"{roi.roi_id}_difference_mask_summary.png"
    deps = [sample.dapi_path] if sample.dapi_path else []
    if _png_skip(out, deps):
        return out

    # Collect per-method numbers
    rows = []
    for m in methods:
        d = compute_difference_mask(sample, roi, baseline, m)
        if not d.get("available"):
            continue
        a = d["area_summary"]
        denom_a = a["baseline_only_frac"] + a["shared_frac"]
        denom_b = a["method_only_frac"] + a["shared_frac"]
        cov_a = (a["shared_frac"] / denom_a) if denom_a > 0 else float("nan")
        cov_b = (a["shared_frac"] / denom_b) if denom_b > 0 else float("nan")
        rows.append({
            "method": m,
            "baseline_only": a["baseline_only_frac"] * 100,
            "method_only":   a["method_only_frac"]   * 100,
            "shared":        a["shared_frac"]        * 100,
            "IoU":           a["iou"],
            "cov_a":         cov_a * 100 if cov_a == cov_a else 0,
            "cov_b":         cov_b * 100 if cov_b == cov_b else 0,
        })
    if not rows:
        return None

    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5),
                                facecolor="white", dpi=200,
                                gridspec_kw={"width_ratios": [3, 2]})

    # Left: grouped bar — baseline-only / method-only / shared %
    metrics = ["baseline_only", "method_only", "shared"]
    metric_colors = {
        "baseline_only": st.HIGHLIGHT_COLORS["xenium_only"],
        "method_only":   st.HIGHLIGHT_COLORS["method_only"],
        "shared":        st.HIGHLIGHT_COLORS["shared"],
    }
    n_methods = len(df)
    x = np.arange(n_methods)
    width = 0.27
    ax = axes[0]
    ax.set_facecolor("white")
    for i, met in enumerate(metrics):
        ax.bar(x + (i - 1) * width, df[met], width=width,
                color=metric_colors[met], edgecolor="black", linewidth=0.4,
                label=met.replace("_", "-"))
    ax.set_xticks(x)
    ax.set_xticklabels([st.METHOD_DISPLAY.get(m, m) for m in df["method"]],
                          rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("% of ROI area", fontsize=9)
    ax.set_title("Area composition (a-only / b-only / shared)", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    for s in ax.spines.values():
        s.set_color("#888")

    # Right: IoU + baseline⊂method coverage
    ax2 = axes[1]
    ax2.set_facecolor("white")
    width2 = 0.35
    ax2.bar(x - width2 / 2, df["IoU"] * 100, width=width2,
             color="#888", edgecolor="black", linewidth=0.4, label="IoU (%)")
    ax2.bar(x + width2 / 2, df["cov_a"], width=width2,
             color="#FFD700", edgecolor="black", linewidth=0.4,
             label="baseline ⊂ method (%)")
    ax2.set_xticks(x)
    ax2.set_xticklabels([st.METHOD_DISPLAY.get(m, m) for m in df["method"]],
                           rotation=15, ha="right", fontsize=8)
    ax2.set_ylabel("%", fontsize=9)
    ax2.set_title("Overlap quality (IoU, baseline coverage)", fontsize=10)
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
    for s in ax2.spines.values():
        s.set_color("#888")

    fig.suptitle(
        f"{roi.roi_id} | difference-mask summary  baseline = {baseline}",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(str(out), dpi=200, bbox_inches="tight", pad_inches=0.05,
                 facecolor="white", edgecolor="none")
    plt.close(fig)
    return out
