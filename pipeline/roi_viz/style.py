"""Single source of truth for colors, plot rcParams, scale bar, axes setup.

The aim is that any rendering function in this package uses these helpers and
constants instead of inline values, so that visual style is consistent across
all figures and easy to tweak in one place.
"""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib import patheffects as mpe
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------------------
# Color palette (single source — referenced everywhere)
# ---------------------------------------------------------------------------
HIGHLIGHT_COLORS = {
    "baseline_boundary":  "#FFD700",  # yellow (xenium baseline)
    "method_boundary":    "#00FFFF",  # cyan  (compared method)
    "xenium_only":        "#FF3B30",
    "method_only":        "#3478F6",
    "shared":             "#BFBFBF",
    "tx_unchanged":       "#BFBFBF",
    "tx_reassigned":      "#FF3B30",
    "tx_newly":           "#3478F6",
    "tx_dropped":         "#FF9500",
    "ssam_agree":         "#34C759",
    "ssam_disagree":      "#FF3B30",
    "background":         "#FFFFFF",
    "tissue_outline":     "#444444",
}

METHOD_COLORS = {
    "xenium_nucleus":     "#FFD700",   # gold
    "cellpose_nuclei":    "#FF6347",   # tomato
    "rigid_expansion":    "#32CD32",   # lime green
    "optimal_expansion":  "#A040FF",   # bright violet (moved off magenta so
                                         # the low-VSI pink mask never
                                         # collides with this method's fill)
    "baysor":             "#00BFFF",   # deep sky blue
}

METHOD_DISPLAY = {
    "xenium_nucleus":     "(1) 10x default (nucleus)",
    "cellpose_nuclei":    "(2) Cellpose only",
    "rigid_expansion":    "(3a) Cellpose + rigid expansion",
    "optimal_expansion":  "(3b) Cellpose + optimal expansion",
    "baysor":             "(4) Cellpose + Baysor",
}

# Paper Salas 2025 canonical 4-bucket framing.  (3) splits into rigid vs
# optimal sub-variants, both built on the same Cellpose nucleus prior.
METHOD_GROUPS = {
    "xenium_nucleus":     "(1) 10x default",
    "cellpose_nuclei":    "(2) Cellpose only",
    "rigid_expansion":    "(3) Cellpose + expansion",
    "optimal_expansion":  "(3) Cellpose + expansion",
    "baysor":             "(4) Cellpose + Baysor",
}


# ---------------------------------------------------------------------------
# Global rcParams
# ---------------------------------------------------------------------------
def apply_rc_params() -> None:
    """Apply consistent matplotlib rc params for the whole package."""
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#222222",
        "axes.linewidth": 0.6,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "image.interpolation": "nearest",
    })


# ---------------------------------------------------------------------------
# Colormaps with white "bad" pixel — prevents NaN→black background
# ---------------------------------------------------------------------------
def cmap_with_white_bad(name: str = "gray") -> mcolors.Colormap:
    """Return a copy of *name* colormap whose `bad` (NaN) color is white."""
    cmap = plt.get_cmap(name).copy()
    cmap.set_bad("white", alpha=1.0)
    return cmap


# ---------------------------------------------------------------------------
# Axes setup
# ---------------------------------------------------------------------------
def setup_axes(ax: Axes,
               bbox_um: Optional[Tuple[float, float, float, float]] = None,
               title: Optional[str] = None,
               keep_axis: bool = False,
               equal: bool = True,
               facecolor: Optional[str] = None) -> Axes:
    """Apply consistent axes config: equal aspect, optional bbox.

    *facecolor* is **only** applied when explicitly passed. By default this
    function preserves whatever face color the caller (or upstream helper
    such as ``_draw_dapi``) already set — important because DAPI imshow uses
    alpha < 1 in some rows and a black axes face is required to keep the
    composite from washing out to grey.
    """
    if facecolor is not None:
        ax.set_facecolor(facecolor)
    if bbox_um is not None:
        x0, y0, x1, y1 = bbox_um
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)  # invert y so image-style axes
    if equal:
        # adjustable="box" preserves the explicit limits set above instead of
        # overriding them to match the data range.
        ax.set_aspect("equal", adjustable="box")
    if not keep_axis:
        # IMPORTANT: ax.axis("off") would also hide the axes patch (face),
        # which makes alpha-blended imshow composite against the FIGURE
        # facecolor (typically white) instead of the AXES facecolor we set.
        # That turns "black" backgrounds gray. Instead hide ticks + spines
        # but keep the axes patch visible so the face color is honored.
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


# ---------------------------------------------------------------------------
# Scale bar
# ---------------------------------------------------------------------------
def add_scale_bar(ax: Axes,
                  length_um: float = 50.0,
                  bbox_um: Optional[Tuple[float, float, float, float]] = None,
                  loc: str = "lower right",
                  color: str = "white",
                  fontsize: int = 7,
                  thickness: float = 2.0) -> None:
    """Draw a fixed-length scale bar in µm."""
    if bbox_um is None:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if y0 > y1:
            y0, y1 = y1, y0
    else:
        x0, y0, x1, y1 = bbox_um
    span_x = x1 - x0
    span_y = abs(y1 - y0)
    margin_x = span_x * 0.04
    margin_y = span_y * 0.04

    if loc == "lower right":
        x_end = x1 - margin_x
        x_start = x_end - length_um
        y_pos = y1 - margin_y
        text_align = "right"
    elif loc == "lower left":
        x_start = x0 + margin_x
        x_end = x_start + length_um
        y_pos = y1 - margin_y
        text_align = "left"
    else:
        x_end = x1 - margin_x
        x_start = x_end - length_um
        y_pos = y1 - margin_y
        text_align = "right"

    ax.plot([x_start, x_end], [y_pos, y_pos],
            color=color, linewidth=thickness, solid_capstyle="butt",
            path_effects=[mpe.withStroke(linewidth=thickness + 1.5,
                                         foreground="black")])
    ax.text((x_start + x_end) / 2, y_pos - span_y * 0.02,
            f"{int(length_um)} µm",
            ha="center" if text_align == "right" else "left",
            va="bottom",
            color=color, fontsize=fontsize,
            path_effects=[mpe.withStroke(linewidth=2.0, foreground="black")])


# ---------------------------------------------------------------------------
# Boundary path effects (for high contrast on dark DAPI)
# ---------------------------------------------------------------------------
def boundary_path_effects(halo_lw: float = 2.0):
    """Return matplotlib path effects that add a soft white halo around lines."""
    return [mpe.Stroke(linewidth=halo_lw, foreground="white", alpha=0.65),
            mpe.Normal()]


# ---------------------------------------------------------------------------
# Categorical → RGBA composer (used by difference / reassignment renders)
# ---------------------------------------------------------------------------
def category_array_to_rgba(category: np.ndarray,
                            color_map: dict,
                            default: str = "background") -> np.ndarray:
    """Vectorize a 2D int category map into an (H, W, 4) RGBA array.

    *category*: 2D ndarray of small ints. Each integer key in *color_map* maps
    to a (color_hex, alpha) tuple. Unmapped values get the *default* entry.
    """
    h, w = category.shape
    rgba = np.zeros((h, w, 4), dtype=np.float32)
    base = color_map.get(default, ("#FFFFFF", 0.0))
    rgba[..., :3] = mcolors.to_rgb(base[0])
    rgba[..., 3] = base[1]
    for k, (hex_color, alpha) in color_map.items():
        if k == default:
            continue
        mask = (category == k)
        if not mask.any():
            continue
        rgba[mask, :3] = mcolors.to_rgb(hex_color)
        rgba[mask, 3] = alpha
    return rgba


# ---------------------------------------------------------------------------
# Convenience figure constructor
# ---------------------------------------------------------------------------
def new_figure(figsize: Tuple[float, float] = (6, 6),
               dpi: int = 200) -> Tuple[plt.Figure, Axes]:
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, facecolor="white")
    ax.set_facecolor("white")
    return fig, ax


def save_figure(fig: plt.Figure, path: str,
                dpi: int = 200,
                close: bool = True) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white", edgecolor="none")
    if close:
        plt.close(fig)
