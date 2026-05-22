"""ROI-class-specific add-on renderers.

Each ``render_*_addons`` function targets one ROI class and writes its PNGs
under ``{outdir}/{roi_id}/{class_name}/...``. When required reference data
is missing, the function logs a warning and skips just that image.

The add-on set per class is intentionally minimal — covering the most
diagnostic image from the plan. Less critical add-ons (3D inset, transect
plots) can be added incrementally without breaking callers.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import cache, style as st
from .common_renders import (_add_dynamic_scale_bar, _draw_dapi,
                                _draw_polygons_outline, _png_skip,
                                _scatter_transcripts)
from .cropping import (crop_dapi, crop_p2r, crop_ssam, crop_transcripts, crop_vsi)
from .data_classes import ROIRecord, SampleData
from .diff_compute import compute_difference_mask, compute_reassignment
from .loaders import load_method
from .markers import DEFAULT_BRAIN_MARKERS, marker_color_palette, select_markers
from .ssam_utils import (coarsen_celltype, grid_majority,
                          method_celltype_per_polygon,
                          ssam_majority_in_polygon, stable_celltype_palette)

logger = logging.getLogger(__name__)


def _class_dir(outdir: Path, roi_id: str, roi_class: str) -> Path:
    d = outdir / roi_id / roi_class
    d.mkdir(parents=True, exist_ok=True)
    return d


def _save(fig, path: Path) -> Path:
    st.save_figure(fig, str(path))
    return path


# ---------------------------------------------------------------------------
# 1. easy_control: clean boundary-only + cell area map
# ---------------------------------------------------------------------------
def render_easy_control_addons(sample: SampleData, roi: ROIRecord,
                                  baseline: str, methods: Sequence[str],
                                  outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "easy_control")

    for m in methods:
        out = out_dir / f"{roi.roi_id}_boundary_only_{m}.png"
        if not _png_skip(out, []):
            res = load_method(m, sample, roi)
            fig, ax = st.new_figure(figsize=(4.5, 4.5))
            color = st.METHOD_COLORS.get(m, "#00FFFF")
            if res["available"]:
                _draw_polygons_outline(ax, res["polygons"], color=color,
                                         lw=0.8, alpha=0.95, halo=False)
            st.setup_axes(ax, bbox_um=roi.bbox_um,
                            title=f"{roi.roi_id} | clean boundary | {m}")
            _add_dynamic_scale_bar(ax, roi, sample, color="black")
            out_paths[f"boundary_only/{m}"] = _save(fig, out)

    # Cell-area map (per method)
    for m in methods:
        out = out_dir / f"{roi.roi_id}_cell_area_{m}.png"
        if _png_skip(out, []):
            continue
        res = load_method(m, sample, roi)
        fig, ax = st.new_figure(figsize=(5, 5))
        dapi = crop_dapi(sample, roi)
        _draw_dapi(ax, dapi, roi, alpha=0.4)
        if res["available"] and "geometry" in res["polygons"].columns and len(res["polygons"]):
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Polygon as MplPolygon
            areas = []
            patches = []
            for g in res["polygons"]["geometry"]:
                if g is None or g.is_empty:
                    continue
                if g.geom_type == "Polygon":
                    xs, ys = zip(*list(g.exterior.coords))
                    patches.append(MplPolygon(np.column_stack([xs, ys]), closed=True))
                    areas.append(g.area)
            if patches:
                pc = PatchCollection(patches, cmap="viridis", alpha=0.7,
                                       edgecolors="black", linewidths=0.2)
                pc.set_array(np.array(areas))
                ax.add_collection(pc)
                cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
                cb.set_label("cell area (µm²)", fontsize=8)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | cell area | {m}")
        _add_dynamic_scale_bar(ax, roi, sample, color="white")
        out_paths[f"cell_area/{m}"] = _save(fig, out)
    return out_paths


# ---------------------------------------------------------------------------
# 2. architecture: SSAM side-by-side + agreement
# ---------------------------------------------------------------------------
def render_architecture_addons(sample: SampleData, roi: ROIRecord,
                                  baseline: str, methods: Sequence[str],
                                  outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "architecture")
    ssam = crop_ssam(sample, roi)
    if len(ssam) == 0:
        logger.warning(f"[{roi.roi_id}] SSAM unavailable for architecture add-ons")
        return out_paths

    label_col = next((c for c in ("celltype", "ssam_celltype",
                                       "leiden_assignment", "leiden")
                       if c in ssam.columns), None)
    if label_col is None:
        return out_paths
    cats = pd.Categorical(ssam[label_col].astype(str))
    cmap = plt.get_cmap("tab20", max(2, len(cats.categories)))

    # SSAM crop image (always)
    out = out_dir / f"{roi.roi_id}_ssam_crop.png"
    if not _png_skip(out, [sample.ssam_h5ad_path] if sample.ssam_h5ad_path else []):
        fig, ax = st.new_figure(figsize=(5, 5))
        ax.scatter(ssam["x_um"], ssam["y_um"], c=cats.codes, cmap=cmap, s=4, alpha=0.8)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | SSAM celltypes ({label_col})")
        _add_dynamic_scale_bar(ax, roi, sample, color="black")
        out_paths["ssam_crop"] = _save(fig, out)

    # SSAM vs each method centroid alignment
    for m in methods:
        out = out_dir / f"{roi.roi_id}_ssam_alignment_{m}.png"
        if _png_skip(out, []):
            continue
        res = load_method(m, sample, roi)
        fig, ax = st.new_figure(figsize=(5.4, 5))
        ax.scatter(ssam["x_um"], ssam["y_um"], c=cats.codes, cmap=cmap,
                    s=2, alpha=0.5)
        if res["available"] and len(res["centroids"]):
            ax.scatter(res["centroids"]["x_centroid_um"],
                        res["centroids"]["y_centroid_um"],
                        s=12, facecolors="none",
                        edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                        linewidths=0.6)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | SSAM + {m} centroids")
        _add_dynamic_scale_bar(ax, roi, sample, color="black")
        out_paths[f"ssam_alignment/{m}"] = _save(fig, out)

    # ---------------------------------------------------------------
    # A1–A5: SSAM cell-type transition diagnostics (per method)
    # ---------------------------------------------------------------
    # Stable celltype color palette: union of SSAM celltypes and any step6
    # celltypes resolved across methods, sorted alphabetically. Built once.
    palette_seed = list(ssam[label_col].astype(str).unique())
    for m in methods:
        try:
            res_m = load_method(m, sample, roi)
            ct = method_celltype_per_polygon(
                sample, roi, m, res_m["polygons"],
                res_m["centroids"], res_m["transcripts"])
            palette_seed.extend([c for c in ct.unique() if c])
        except Exception:
            pass
    celltype_palette = stable_celltype_palette(palette_seed)

    # Marker pair for marker_exclusivity / boundary_transect (architecture
    # panel = neuron / inhibitory / astrocyte). Pick first two present.
    arch_pair = _arch_marker_pair(sample, roi)

    for m in methods:
        try:
            res_m = load_method(m, sample, roi)
        except Exception as e:
            logger.warning(f"[{roi.roi_id}] load_method({m}) failed: {e}")
            continue
        if not res_m.get("available") or "geometry" not in res_m["polygons"].columns \
           or len(res_m["polygons"]) == 0:
            continue

        # Compute once per method (used by A1, A3, A4)
        method_ct = method_celltype_per_polygon(
            sample, roi, m, res_m["polygons"], res_m["centroids"],
            res_m["transcripts"])
        ssam_ct = ssam_majority_in_polygon(ssam, res_m["polygons"], label_col)

        # Detect classifier-bias once per (ROI, method).  All three classifier-
        # dependent add-ons surface this in their title so the reader doesn't
        # over-interpret a misleading agree/disagree pattern.
        caveat = _classifier_bias_note(method_ct, ssam, label_col,
                                            crop_transcripts(sample, roi))

        # ---- A1 ssam_disagreement -----------------------------------
        out = out_dir / f"{roi.roi_id}_ssam_disagreement_{m}.png"
        if not _png_skip(out, []):
            out_paths[f"ssam_disagreement/{m}"] = _save_ssam_disagreement(
                sample, roi, m, res_m, method_ct, ssam_ct, out, caveat)
        # ---- A2 marker_exclusivity ----------------------------------
        out = out_dir / f"{roi.roi_id}_marker_exclusivity_{m}.png"
        if not _png_skip(out, []) and arch_pair:
            out_paths[f"marker_exclusivity/{m}"] = _save_marker_exclusivity(
                sample, roi, m, res_m, arch_pair, out)
        # ---- A3 ssam_agreement_heatmap -------------------------------
        out = out_dir / f"{roi.roi_id}_ssam_agreement_heatmap_{m}.png"
        if not _png_skip(out, []):
            out_paths[f"ssam_agreement_heatmap/{m}"] = _save_ssam_agreement_heatmap(
                sample, roi, m, ssam, label_col, res_m, method_ct, out, caveat)
        # ---- A4 ssam_sidebyside --------------------------------------
        out = out_dir / f"{roi.roi_id}_ssam_sidebyside_{m}.png"
        if not _png_skip(out, []):
            out_paths[f"ssam_sidebyside/{m}"] = _save_ssam_sidebyside(
                sample, roi, m, ssam, label_col, res_m, method_ct,
                celltype_palette, out, caveat)
        # ---- A5 boundary_transect ------------------------------------
        if arch_pair:
            out = out_dir / f"{roi.roi_id}_boundary_transect_{m}.png"
            if not _png_skip(out, []):
                out_paths[f"boundary_transect/{m}"] = _save_boundary_transect(
                    sample, roi, m, arch_pair, out)

    # ------------------------------------------------------------------
    # P1 chimera-cell highlight + per-method panel_compare
    # ------------------------------------------------------------------
    # For each cell: chimera_score = 1 - (max_celltype_count / total_count)
    # of SSAM dots inside the cell.  Pool across methods, use shared q75/q50
    # threshold so the red intensity is directly comparable across methods.
    try:
        import geopandas as gpd
    except Exception:
        gpd = None
    if gpd is not None and len(ssam):
        ssam_gdf = gpd.GeoDataFrame(
            {"label": ssam[label_col].astype(str).values},
            geometry=gpd.points_from_xy(ssam["x_um"], ssam["y_um"]),
            crs=None)

        method_chim: Dict[str, "pd.DataFrame"] = {}
        for m in methods:
            try:
                res_m = load_method(m, sample, roi)
            except Exception:
                continue
            if not res_m.get("available") or "geometry" not in res_m["polygons"].columns \
               or len(res_m["polygons"]) == 0:
                continue
            cells_gdf = res_m["polygons"][res_m["polygons"].geometry.notna()].copy()
            cells_gdf = cells_gdf.reset_index(drop=True)
            cells_gdf["cell_idx"] = cells_gdf.index

            chimera = pd.Series(np.nan, index=cells_gdf.index)
            j = gpd.sjoin(ssam_gdf, cells_gdf[["cell_idx", "geometry"]],
                            predicate="within", how="inner")
            if len(j):
                grp = j.groupby("cell_idx")["label"]
                # max single-celltype count per cell vs total dot count
                def _mix(x):
                    vc = x.value_counts()
                    if len(vc) == 0 or vc.sum() == 0:
                        return np.nan
                    return 1.0 - vc.iloc[0] / vc.sum()
                chimera_vals = grp.apply(_mix)
                chimera.loc[chimera_vals.index] = chimera_vals.values
            cells_gdf["chimera"] = chimera
            method_chim[m] = cells_gdf

        # Pooled thresholds
        if method_chim:
            pooled = pd.concat([cdf["chimera"] for cdf in method_chim.values()],
                                  ignore_index=True).dropna()
            if len(pooled):
                gq75 = float(pooled.quantile(0.75))
                gq50 = float(pooled.quantile(0.50))
            else:
                gq75 = gq50 = float("nan")

            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Polygon as MplPolygon, Patch

            # ---- chimera_cell_highlight_{method}.png × 5
            for m in methods:
                if m not in method_chim:
                    continue
                out = out_dir / f"{roi.roi_id}_chimera_cell_highlight_{m}.png"
                if _png_skip(out, []):
                    continue
                cells_gdf = method_chim[m]
                bad_p, mid_p, ok_p = [], [], []
                for _, row in cells_gdf.iterrows():
                    g = row["geometry"]
                    if g is None or g.is_empty or g.geom_type != "Polygon":
                        continue
                    xs, ys = zip(*list(g.exterior.coords))
                    poly = MplPolygon(np.column_stack([xs, ys]), closed=True)
                    cm = row["chimera"]
                    if np.isnan(cm):
                        ok_p.append(poly)
                    elif not np.isnan(gq75) and cm >= gq75:
                        bad_p.append(poly)
                    elif not np.isnan(gq50) and cm >= gq50:
                        mid_p.append(poly)
                    else:
                        ok_p.append(poly)
                fig, ax = plt.subplots(figsize=(5.6, 5.4),
                                            facecolor="black", dpi=180)
                # SSAM dots first (faint reference)
                ax.scatter(ssam["x_um"], ssam["y_um"], c=cats.codes,
                            cmap=cmap, s=2, alpha=0.45)
                if ok_p:
                    ax.add_collection(PatchCollection(
                        ok_p, facecolors="none",
                        edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                        linewidths=0.5, alpha=0.85))
                if mid_p:
                    ax.add_collection(PatchCollection(
                        mid_p, facecolors="#FFA50066",
                        edgecolors="#FFA500", linewidths=1.0, alpha=0.85))
                if bad_p:
                    ax.add_collection(PatchCollection(
                        bad_p, facecolors="#FF000099",
                        edgecolors="#FF0000", linewidths=1.4, alpha=0.95))

                valid_m = cells_gdf["chimera"].dropna()
                pct_bad = (100.0 * (valid_m >= gq75).sum() / len(valid_m)
                            if len(valid_m) and not np.isnan(gq75) else float("nan"))
                handles = [
                    Patch(facecolor="#FF000099", edgecolor="#FF0000",
                          label=(f"top-25% pooled (≥{gq75:.2f}) this method: {pct_bad:.0f}%"
                                  if not np.isnan(gq75) else "top-25% pooled")),
                    Patch(facecolor="#FFA50066", edgecolor="#FFA500",
                          label=f"50-75% pooled (≥{gq50:.2f})" if not np.isnan(gq50) else "mid"),
                    plt.Line2D([], [], color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                lw=1.5, label="bottom 50% (pure)"),
                ]
                leg = ax.legend(handles=handles, loc="upper right", fontsize=7,
                                  framealpha=0.85, facecolor="black",
                                  edgecolor="white")
                for t in leg.get_texts():
                    t.set_color("white")
                ax.set_title(
                    f"{roi.roi_id} | chimera cells | "
                    f"{st.METHOD_DISPLAY.get(m, m)}\n"
                    f"score = 1 − (majority SSAM celltype fraction in cell), "
                    f"red = pooled top-25%",
                    fontsize=8, color="white")
                st.setup_axes(ax, bbox_um=roi.bbox_um)
                _add_dynamic_scale_bar(ax, roi, sample, color="white")
                fig.tight_layout()
                fig.savefig(str(out), dpi=180, bbox_inches="tight",
                              pad_inches=0.05, facecolor="black")
                plt.close(fig)
                out_paths[f"chimera_cell_highlight/{m}"] = out

            # ---- architecture_panel_compare.png (2×3)
            out_panel = out_dir / f"{roi.roi_id}_architecture_panel_compare.png"
            if not _png_skip(out_panel, []):
                n = len(methods)
                n_cols = 3 if n > 3 else n
                n_rows = (n + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols,
                                              figsize=(5.0 * n_cols, 5.0 * n_rows + 0.6),
                                              facecolor="black", dpi=160)
                if n_rows == 1 and n_cols == 1:
                    axes = np.array([[axes]])
                elif n_rows == 1:
                    axes = axes[None, :]
                elif n_cols == 1:
                    axes = axes[:, None]
                for i, m in enumerate(methods):
                    r, c = i // n_cols, i % n_cols
                    ax = axes[r, c]
                    ax.set_facecolor("black")
                    ax.scatter(ssam["x_um"], ssam["y_um"], c=cats.codes,
                                cmap=cmap, s=3, alpha=0.6)
                    if m in method_chim:
                        cells_gdf = method_chim[m]
                        bad_p, ok_p = [], []
                        for _, row in cells_gdf.iterrows():
                            g = row["geometry"]
                            if g is None or g.is_empty or g.geom_type != "Polygon":
                                continue
                            xs, ys = zip(*list(g.exterior.coords))
                            poly = MplPolygon(np.column_stack([xs, ys]), closed=True)
                            if not np.isnan(row["chimera"]) and \
                               not np.isnan(gq75) and row["chimera"] >= gq75:
                                bad_p.append(poly)
                            else:
                                ok_p.append(poly)
                        if ok_p:
                            ax.add_collection(PatchCollection(
                                ok_p, facecolors="none",
                                edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                                linewidths=0.5, alpha=0.90))
                        if bad_p:
                            ax.add_collection(PatchCollection(
                                bad_p, facecolors="#FF000099",
                                edgecolors="#FF0000", linewidths=1.0, alpha=0.95))
                        valid_m = cells_gdf["chimera"].dropna()
                        pure_pct = (100.0 * (valid_m < (gq50 if not np.isnan(gq50) else 0.2)).sum()
                                       / max(1, len(valid_m)))
                        ax.set_title(
                            f"{st.METHOD_DISPLAY.get(m, m)}\n"
                            f"pure cells: {pure_pct:.0f}% | "
                            f"chimera (top-25%): "
                            f"{100*(valid_m>=gq75).sum()/max(1,len(valid_m)):.0f}%",
                            fontsize=10, color="white")
                    else:
                        ax.set_title(f"{st.METHOD_DISPLAY.get(m, m)} (unavailable)",
                                        fontsize=10, color="white")
                    st.setup_axes(ax, bbox_um=roi.bbox_um)
                for j in range(n, n_rows * n_cols):
                    r, c = j // n_cols, j % n_cols
                    axes[r, c].set_visible(False)
                fig.suptitle(
                    f"{roi.roi_id} | architecture per-method comparison "
                    f"(SSAM celltype dots = colored by celltype, "
                    f"red fill = chimera cell)",
                    fontsize=11, color="white")
                fig.tight_layout(rect=(0, 0, 1, 0.95))
                fig.savefig(str(out_panel), dpi=160, bbox_inches="tight",
                              pad_inches=0.05, facecolor="black")
                plt.close(fig)
                out_paths["architecture_panel_compare"] = out_panel

    return out_paths


def _classifier_bias_note(method_ct: pd.Series,
                            ssam: pd.DataFrame, ssam_label_col: Optional[str],
                            tx: Optional[pd.DataFrame]) -> Optional[str]:
    """Return a short warning string when the segmentation-side classifier
    looks marker-majority biased (= falls back to coarse marker counts and
    one gene's transcript count dominates the cell-type call).  Returns
    ``None`` when the comparison appears reliable.

    Trigger: method-side astrocyte fraction > 50% AND
    (SSAM-side astrocyte fraction < 20% OR CLU > 40% of marker tx).
    """
    if method_ct is None or len(method_ct) == 0:
        return None
    nonempty = method_ct[method_ct.astype(str).str.strip() != ""]
    if len(nonempty) == 0:
        return None
    method_astro = float((nonempty.astype(str).str.lower() == "astrocyte").mean())

    if ssam is not None and len(ssam) and ssam_label_col in ssam.columns:
        s_lbl = ssam[ssam_label_col].astype(str).str.lower()
        ssam_astro = float((s_lbl == "astrocyte").mean())
    else:
        ssam_astro = float("nan")

    clu_frac = float("nan")
    if tx is not None and len(tx) and "feature_name" in tx.columns:
        gname = tx["feature_name"].astype(str)
        marker_genes: set = set()
        for cls, gs in DEFAULT_BRAIN_MARKERS.items():
            if cls in ("nuclear", "cytoplasmic"):
                continue
            marker_genes.update(gs)
        n_marker = int(gname.isin(marker_genes).sum())
        n_clu = int((gname == "CLU").sum())
        if n_marker > 0:
            clu_frac = n_clu / n_marker

    suspicious = method_astro > 0.50 and (
        (not np.isnan(ssam_astro) and ssam_astro < 0.20) or
        (not np.isnan(clu_frac) and clu_frac > 0.40)
    )
    if not suspicious:
        return None
    parts = [f"⚠ classifier bias: {100*method_astro:.0f}% labeled astrocyte"]
    if not np.isnan(clu_frac):
        parts.append(f"CLU = {100*clu_frac:.0f}% of marker tx")
    if not np.isnan(ssam_astro):
        parts.append(f"SSAM-side astrocyte = {100*ssam_astro:.0f}%")
    return "  |  ".join(parts)


def _arch_marker_pair(sample: SampleData, roi: ROIRecord
                       ) -> Optional[Tuple[str, str]]:
    """Pick a (target, competing) gene pair for architecture marker views.

    Strategy: take the architecture class panel (neuron / astrocyte /
    inhibitory), keep only genes present in the ROI's transcript table, and
    return the first two distinct classes.
    """
    try:
        tx = crop_transcripts(sample, roi)
    except Exception:
        return None
    if len(tx) == 0:
        return None
    present = set(tx["feature_name"].astype(str).unique())
    classes = ["neuron", "astrocyte", "inhibitory", "oligodendrocyte"]
    picks: list = []
    for cls in classes:
        for g in DEFAULT_BRAIN_MARKERS.get(cls, []):
            if g in present:
                picks.append(g)
                break
        if len(picks) >= 2:
            break
    if len(picks) < 2:
        return None
    return picks[0], picks[1]


def _polygon_patches(polygons: pd.DataFrame):
    """Yield (idx, MplPolygon) for valid polygon rows."""
    from matplotlib.patches import Polygon as MplPolygon
    for idx, geom in enumerate(polygons["geometry"]):
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type != "Polygon":
            continue
        xs, ys = zip(*list(geom.exterior.coords))
        yield idx, MplPolygon(np.column_stack([xs, ys]), closed=True)


# ---------------------------------------------------------------------------
# A1 — SSAM disagreement
# ---------------------------------------------------------------------------
def _save_ssam_disagreement(sample: SampleData, roi: ROIRecord, method: str,
                              res_m: Dict[str, Any],
                              method_ct: pd.Series, ssam_ct: pd.Series,
                              out: Path,
                              caveat: Optional[str] = None) -> Path:
    from matplotlib.collections import PatchCollection
    polys = res_m["polygons"]
    fig, ax = st.new_figure(figsize=(5.4, 5))
    dapi = crop_dapi(sample, roi)
    _draw_dapi(ax, dapi, roi, alpha=0.45)

    match_patches, mismatch_patches, gray_patches = [], [], []
    n_match = n_mismatch = n_gray = 0
    for idx, patch in _polygon_patches(polys):
        m_lbl = coarsen_celltype(method_ct.iloc[idx]) if idx < len(method_ct) else ""
        s_lbl = coarsen_celltype(ssam_ct.iloc[idx]) if idx < len(ssam_ct) else ""
        if not m_lbl or not s_lbl:
            gray_patches.append(patch); n_gray += 1
        elif m_lbl == s_lbl:
            match_patches.append(patch); n_match += 1
        else:
            mismatch_patches.append(patch); n_mismatch += 1

    if match_patches:
        ax.add_collection(PatchCollection(
            match_patches, facecolor=st.HIGHLIGHT_COLORS["ssam_agree"],
            edgecolor="#1F8E3C", linewidths=0.4, alpha=0.55))
    if mismatch_patches:
        ax.add_collection(PatchCollection(
            mismatch_patches, facecolor=st.HIGHLIGHT_COLORS["ssam_disagree"],
            edgecolor="#A91500", linewidths=0.4, alpha=0.65))
    if gray_patches:
        ax.add_collection(PatchCollection(
            gray_patches, facecolor="none",
            edgecolor="#888888", linewidths=0.4, alpha=0.7))

    handles = [
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color=st.HIGHLIGHT_COLORS["ssam_agree"],
                    label=f"agree ({n_match})"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color=st.HIGHLIGHT_COLORS["ssam_disagree"],
                    label=f"disagree ({n_mismatch})"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color="#888888", label=f"undefined ({n_gray})"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.85)
    title = (f"{roi.roi_id} | SSAM ↔ {st.METHOD_DISPLAY.get(method, method)} | "
             f"disagreement")
    if caveat:
        title += f"\n{caveat}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# A2 — marker exclusivity
# ---------------------------------------------------------------------------
def _save_marker_exclusivity(sample: SampleData, roi: ROIRecord, method: str,
                               res_m: Dict[str, Any],
                               arch_pair: Tuple[str, str], out: Path) -> Path:
    from matplotlib.collections import PatchCollection
    target_g, competing_g = arch_pair
    polys = res_m["polygons"]
    tx = res_m.get("transcripts", pd.DataFrame())
    fig, ax = st.new_figure(figsize=(5.4, 5))
    dapi = crop_dapi(sample, roi)
    _draw_dapi(ax, dapi, roi, alpha=0.4)

    # Per-cell counts of target / competing
    if "cell_id" in tx.columns and "feature_name" in tx.columns:
        sub = tx[tx["feature_name"].isin([target_g, competing_g])]
        sub = sub[sub["assigned"]] if "assigned" in sub.columns else sub
        if len(sub):
            counts = (sub.groupby([sub["cell_id"].astype(str),
                                      "feature_name"]).size().unstack(fill_value=0))
        else:
            counts = pd.DataFrame()
    else:
        counts = pd.DataFrame()

    mixed = set()
    target_only = set()
    competing_only = set()
    if len(counts):
        for cid, row in counts.iterrows():
            t = int(row.get(target_g, 0))
            c = int(row.get(competing_g, 0))
            if t >= 2 and c >= 2:
                mixed.add(str(cid))
            elif t >= 2:
                target_only.add(str(cid))
            elif c >= 2:
                competing_only.add(str(cid))

    bg, mix_p, t_p, c_p = [], [], [], []
    for idx, patch in _polygon_patches(polys):
        cid = str(polys["cell_id"].iloc[idx]) if "cell_id" in polys.columns else str(idx)
        if cid in mixed:
            mix_p.append(patch)
        elif cid in target_only:
            t_p.append(patch)
        elif cid in competing_only:
            c_p.append(patch)
        else:
            bg.append(patch)
    if bg:
        ax.add_collection(PatchCollection(
            bg, facecolor="none", edgecolor="#888888", linewidths=0.3, alpha=0.5))
    if t_p:
        ax.add_collection(PatchCollection(
            t_p, facecolor="#FF3B30", edgecolor="#A91500",
            linewidths=0.3, alpha=0.45))
    if c_p:
        ax.add_collection(PatchCollection(
            c_p, facecolor="#34C759", edgecolor="#1F8E3C",
            linewidths=0.3, alpha=0.45))
    if mix_p:
        ax.add_collection(PatchCollection(
            mix_p, facecolor="#FFD60A", edgecolor="#7A5C00",
            linewidths=0.6, alpha=0.85))

    handles = [
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color="#FF3B30", label=f"{target_g} only ({len(target_only)})"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color="#34C759", label=f"{competing_g} only ({len(competing_only)})"),
        plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                    color="#FFD60A", label=f"mixed ({len(mixed)})"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.85)
    st.setup_axes(ax, bbox_um=roi.bbox_um,
                    title=f"{roi.roi_id} | marker exclusivity | "
                          f"{target_g} vs {competing_g} | "
                          f"{st.METHOD_DISPLAY.get(method, method)}")
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# A3 — agreement heatmap (grid majority)
# ---------------------------------------------------------------------------
def _save_ssam_agreement_heatmap(sample: SampleData, roi: ROIRecord,
                                    method: str, ssam: pd.DataFrame,
                                    label_col: str, res_m: Dict[str, Any],
                                    method_ct: pd.Series, out: Path,
                                    caveat: Optional[str] = None) -> Path:
    # SSAM grid
    ssam_grid, ssam_labels, ext = grid_majority(
        ssam, "x_um", "y_um", label_col, roi.bbox_um, grid_um=30.0)
    # Method grid: build (x, y, celltype) from centroids + method_ct
    cents = res_m["centroids"]
    if len(cents) and len(method_ct) == len(res_m["polygons"]):
        # Need centroid → polygon idx alignment. Both come from same loader so
        # we assume row order matches when cell_id matches. Build a per-centroid
        # celltype via cell_id mapping if available.
        if "cell_id" in cents.columns and "cell_id" in res_m["polygons"].columns:
            ct_by_cid = dict(zip(
                res_m["polygons"]["cell_id"].astype(str),
                method_ct.astype(str).values))
            ct_per_centroid = [ct_by_cid.get(str(c), "")
                                for c in cents["cell_id"].astype(str)]
        else:
            ct_per_centroid = list(method_ct.astype(str).values[:len(cents)])
        method_df = pd.DataFrame({
            "x": cents["x_centroid_um"].astype(float).values,
            "y": cents["y_centroid_um"].astype(float).values,
            "celltype": ct_per_centroid,
        })
    else:
        method_df = pd.DataFrame(columns=["x", "y", "celltype"])
    method_grid, method_labels, _ = grid_majority(
        method_df, "x", "y", "celltype", roi.bbox_um, grid_um=30.0)

    # Translate label indices back to strings for cross-grid comparison.
    def _to_str(grid: np.ndarray, labels: list) -> np.ndarray:
        out = np.full(grid.shape, "", dtype=object)
        for i, lbl in enumerate(labels):
            out[grid == i] = lbl
        return out
    s_str = _to_str(ssam_grid, ssam_labels)
    m_str = _to_str(method_grid, method_labels)
    # Collapse SSAM ("L6 IT") and method ("oligodendrocyte" / "Oligodendrocyte")
    # labels to the same coarse vocabulary so the heatmap is meaningful.
    s_str = np.vectorize(coarsen_celltype, otypes=[object])(s_str)
    m_str = np.vectorize(coarsen_celltype, otypes=[object])(m_str)
    s_str = s_str.astype("U64"); m_str = m_str.astype("U64")
    both = (s_str != "") & (m_str != "")
    agree = (s_str == m_str) & both
    score = np.full(s_str.shape, np.nan, dtype=float)
    score[both] = agree[both].astype(float)

    fig, ax = st.new_figure(figsize=(5.6, 5))
    dapi = crop_dapi(sample, roi)
    _draw_dapi(ax, dapi, roi, alpha=0.55)
    cmap = plt.get_cmap("RdYlGn")
    cmap = cmap.copy(); cmap.set_bad(color=(0, 0, 0, 0))
    masked = np.ma.masked_invalid(score)
    ax.imshow(masked, cmap=cmap, vmin=0, vmax=1,
                extent=ext, alpha=0.65, interpolation="nearest")
    pct = 100.0 * np.nansum(score) / max(1, np.sum(np.isfinite(score)))
    title = (f"{roi.roi_id} | SSAM↔{st.METHOD_DISPLAY.get(method, method)} | "
             f"agreement heatmap (30 µm bins, mean = {pct:.0f}%)")
    if caveat:
        title += f"\n{caveat}"
    st.setup_axes(ax, bbox_um=roi.bbox_um, title=title)
    _add_dynamic_scale_bar(ax, roi, sample, color="white")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("agreement", fontsize=8)
    return _save(fig, out)


# ---------------------------------------------------------------------------
# A4 — side-by-side SSAM scatter vs method polygon fill
# ---------------------------------------------------------------------------
def _save_ssam_sidebyside(sample: SampleData, roi: ROIRecord, method: str,
                            ssam: pd.DataFrame, label_col: str,
                            res_m: Dict[str, Any], method_ct: pd.Series,
                            palette: Dict[str, str], out: Path,
                            caveat: Optional[str] = None) -> Path:
    from matplotlib.collections import PatchCollection
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.4), dpi=200,
                                facecolor="white")
    ax_l, ax_r = axes

    # Left: SSAM scatter (stable celltype palette)
    s_lbl = ssam[label_col].astype(str).values
    s_colors = [palette.get(l, "#CCCCCC") for l in s_lbl]
    ax_l.scatter(ssam["x_um"], ssam["y_um"], c=s_colors, s=4, alpha=0.85,
                  edgecolors="none")
    st.setup_axes(ax_l, bbox_um=roi.bbox_um,
                    title=f"{roi.roi_id} | SSAM celltype")
    _add_dynamic_scale_bar(ax_l, roi, sample, color="black")

    # Right: method polygon fill
    polys = res_m["polygons"]
    grouped: Dict[str, list] = {}
    for idx, patch in _polygon_patches(polys):
        lbl = str(method_ct.iloc[idx]) if idx < len(method_ct) else ""
        grouped.setdefault(lbl, []).append(patch)
    dapi = crop_dapi(sample, roi)
    _draw_dapi(ax_r, dapi, roi, alpha=0.35)
    for lbl, patches in grouped.items():
        color = palette.get(lbl, "#CCCCCC") if lbl else "#888888"
        ax_r.add_collection(PatchCollection(
            patches, facecolor=color, edgecolor="black",
            linewidths=0.3, alpha=0.7 if lbl else 0.3))
    r_title = f"{roi.roi_id} | {st.METHOD_DISPLAY.get(method, method)} | celltype"
    if caveat:
        r_title += f"\n{caveat}"
    st.setup_axes(ax_r, bbox_um=roi.bbox_um, title=r_title)
    _add_dynamic_scale_bar(ax_r, roi, sample, color="white")

    # Shared compact legend (top of figure)
    handles = [plt.Line2D([], [], marker="s", linestyle="", markersize=10,
                            color=palette[k], label=k)
                for k in sorted(palette.keys())[:14]]
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=min(7, len(handles)),
                    fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(str(out), dpi=200, bbox_inches="tight", pad_inches=0.05,
                  facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# A5 — boundary transect (line profile)
# ---------------------------------------------------------------------------
def _save_boundary_transect(sample: SampleData, roi: ROIRecord, method: str,
                              arch_pair: Tuple[str, str], out: Path) -> Path:
    target_g, competing_g = arch_pair
    tx = crop_transcripts(sample, roi)
    if len(tx) == 0:
        return out
    x0, y0, x1, y1 = roi.bbox_um
    width = x1 - x0
    height = y1 - y0
    short = "y_location" if height < width else "x_location"
    short_lo = y0 if short == "y_location" else x0
    short_hi = y1 if short == "y_location" else x1
    bin_um = 5.0
    edges = np.arange(short_lo, short_hi + bin_um, bin_um)

    sub_t = tx[tx["feature_name"].astype(str) == target_g]
    sub_c = tx[tx["feature_name"].astype(str) == competing_g]
    h_t, _ = np.histogram(sub_t[short].astype(float).values, bins=edges)
    h_c, _ = np.histogram(sub_c[short].astype(float).values, bins=edges)

    # Smooth a little (3-bin moving average) for readability
    if len(h_t) >= 3:
        kernel = np.ones(3) / 3.0
        h_t_s = np.convolve(h_t, kernel, mode="same")
        h_c_s = np.convolve(h_c, kernel, mode="same")
    else:
        h_t_s, h_c_s = h_t, h_c

    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = st.new_figure(figsize=(6.0, 3.6))
    ax.plot(centers, h_t_s, color="#FF3B30", lw=1.4, label=f"{target_g}")
    ax.plot(centers, h_c_s, color="#34C759", lw=1.4, label=f"{competing_g}")
    ax.fill_between(centers, h_t_s, h_c_s,
                     where=h_t_s > h_c_s, alpha=0.15, color="#FF3B30")
    ax.fill_between(centers, h_t_s, h_c_s,
                     where=h_c_s > h_t_s, alpha=0.15, color="#34C759")
    ax.set_xlabel(f"{short.replace('_location','')} (µm)")
    ax.set_ylabel("transcripts / 5 µm")
    ax.set_title(f"{roi.roi_id} | boundary transect | "
                  f"{st.METHOD_DISPLAY.get(method, method)}")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    st.setup_axes(ax, keep_axis=True, equal=False)
    return _save(fig, out)


# ---------------------------------------------------------------------------
# 3. compartment: nuclear vs cytoplasmic markers + nucleus-distance plot
# ---------------------------------------------------------------------------
def render_compartment_addons(sample: SampleData, roi: ROIRecord,
                                  baseline: str, methods: Sequence[str],
                                  outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "compartment")
    tx = crop_transcripts(sample, roi)
    if len(tx) == 0:
        return out_paths
    dapi = crop_dapi(sample, roi)

    nuclear = [g for g in DEFAULT_BRAIN_MARKERS["nuclear"]
               if g in tx["feature_name"].astype(str).unique()]
    cyto = [g for g in DEFAULT_BRAIN_MARKERS["cytoplasmic"]
            if g in tx["feature_name"].astype(str).unique()]

    for m in methods:
        # Nuclear marker map
        if nuclear:
            out = out_dir / f"{roi.roi_id}_nuclear_marker_{m}.png"
            if not _png_skip(out, []):
                fig, ax = st.new_figure(figsize=(5.2, 5))
                _draw_dapi(ax, dapi, roi, alpha=0.55)
                pal = marker_color_palette(nuclear)
                for g in nuclear:
                    sub = tx[tx["feature_name"].astype(str) == g]
                    ax.scatter(sub["x_location"], sub["y_location"],
                                s=4, c=pal[g], alpha=0.85,
                                edgecolors="white", linewidths=0.15,
                                label=f"{g} ({len(sub):,})")
                res = load_method(m, sample, roi)
                if res["available"]:
                    _draw_polygons_outline(ax, res["polygons"],
                                             color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                             lw=0.7, alpha=0.85)
                ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
                st.setup_axes(ax, bbox_um=roi.bbox_um,
                                title=f"{roi.roi_id} | nuclear markers | {m}")
                _add_dynamic_scale_bar(ax, roi, sample, color="white")
                out_paths[f"nuclear_marker/{m}"] = _save(fig, out)
        # Cytoplasmic marker map
        if cyto:
            out = out_dir / f"{roi.roi_id}_cytoplasmic_marker_{m}.png"
            if not _png_skip(out, []):
                fig, ax = st.new_figure(figsize=(5.2, 5))
                _draw_dapi(ax, dapi, roi, alpha=0.55)
                pal = marker_color_palette(cyto)
                for g in cyto:
                    sub = tx[tx["feature_name"].astype(str) == g]
                    ax.scatter(sub["x_location"], sub["y_location"],
                                s=4, c=pal[g], alpha=0.85,
                                edgecolors="white", linewidths=0.15,
                                label=f"{g} ({len(sub):,})")
                res = load_method(m, sample, roi)
                if res["available"]:
                    _draw_polygons_outline(ax, res["polygons"],
                                             color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                             lw=0.7, alpha=0.85)
                ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
                st.setup_axes(ax, bbox_um=roi.bbox_um,
                                title=f"{roi.roi_id} | cytoplasmic markers | {m}")
                _add_dynamic_scale_bar(ax, roi, sample, color="white")
                out_paths[f"cytoplasmic_marker/{m}"] = _save(fig, out)

    # ------------------------------------------------------------------
    # P1/P2 compartment add-ons (plan/roi_visualization §5.3)
    # ------------------------------------------------------------------
    NUC_COLOR = "#1E90FF"   # dodger blue — nuclear marker dot
    CYTO_COLOR = "#FF8C00"  # dark orange — cytoplasmic marker dot

    nuclear_tx = tx[tx["feature_name"].astype(str).isin(nuclear)] if nuclear else tx.iloc[0:0]
    cyto_tx = tx[tx["feature_name"].astype(str).isin(cyto)] if cyto else tx.iloc[0:0]

    nucleus_res = load_method("xenium_nucleus", sample, roi)
    nuc_polys = nucleus_res["polygons"] if nucleus_res.get("available") else None

    def _per_method_compartment_stats(cell_polys: pd.DataFrame) -> Tuple[float, float]:
        """Return (nuclear-in-nucleus%, cyto-extranuclear%) for cells in this
        method, using per-tx ``overlaps_nucleus`` flag.  Higher = better
        compartment fidelity."""
        if "overlaps_nucleus" not in tx.columns or cell_polys is None or len(cell_polys) == 0:
            return float("nan"), float("nan")
        try:
            import geopandas as gpd
            from shapely.geometry import Point
        except Exception:
            return float("nan"), float("nan")
        nuc_in = float("nan")
        cyto_out = float("nan")
        cells_gdf = cell_polys[cell_polys.geometry.notna()].copy()
        if len(cells_gdf) == 0:
            return nuc_in, cyto_out
        if len(nuclear_tx):
            n_gdf = gpd.GeoDataFrame(
                {"on": nuclear_tx["overlaps_nucleus"].fillna(0).astype(int).values},
                geometry=gpd.points_from_xy(nuclear_tx["x_location"], nuclear_tx["y_location"]),
                crs=None)
            j = gpd.sjoin(n_gdf, cells_gdf[["geometry"]], predicate="within", how="inner")
            if len(j):
                nuc_in = 100.0 * (j["on"] == 1).sum() / len(j)
        if len(cyto_tx):
            c_gdf = gpd.GeoDataFrame(
                {"on": cyto_tx["overlaps_nucleus"].fillna(0).astype(int).values},
                geometry=gpd.points_from_xy(cyto_tx["x_location"], cyto_tx["y_location"]),
                crs=None)
            j = gpd.sjoin(c_gdf, cells_gdf[["geometry"]], predicate="within", how="inner")
            if len(j):
                cyto_out = 100.0 * (j["on"] == 0).sum() / len(j)
        return nuc_in, cyto_out

    # Pre-compute per-method stats once (re-used by overlay + panel_compare)
    method_stats: Dict[str, Tuple[float, float]] = {}
    for m in methods:
        try:
            res_m = load_method(m, sample, roi)
            method_stats[m] = _per_method_compartment_stats(
                res_m["polygons"] if res_m.get("available") else None)
        except Exception as e:
            logger.debug(f"compartment stats {m}: {e}")
            method_stats[m] = (float("nan"), float("nan"))

    # P1.a: compartment_overlay_{method}.png × 5
    if (len(nuclear_tx) + len(cyto_tx)) > 0:
        for m in methods:
            out = out_dir / f"{roi.roi_id}_compartment_overlay_{m}.png"
            if _png_skip(out, []):
                continue
            fig, ax = plt.subplots(figsize=(5.6, 5.4), facecolor="black", dpi=180)
            _draw_dapi(ax, dapi, roi, alpha=0.55)
            if len(nuclear_tx):
                ax.scatter(nuclear_tx["x_location"], nuclear_tx["y_location"],
                            s=5, c=NUC_COLOR, alpha=0.90,
                            edgecolors="white", linewidths=0.20,
                            label=f"nuclear markers ({len(nuclear_tx):,})")
            if len(cyto_tx):
                ax.scatter(cyto_tx["x_location"], cyto_tx["y_location"],
                            s=4, c=CYTO_COLOR, alpha=0.80,
                            edgecolors="white", linewidths=0.15,
                            label=f"cyto markers ({len(cyto_tx):,})")
            if nuc_polys is not None:
                _draw_polygons_outline(ax, nuc_polys, color="#FFFFFF",
                                         lw=0.4, alpha=0.55, halo=False)
            res_m = load_method(m, sample, roi)
            if res_m["available"]:
                _draw_polygons_outline(ax, res_m["polygons"],
                                         color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                         lw=1.1, alpha=0.95, halo=False)
            nuc_in, cyto_out = method_stats.get(m, (float("nan"), float("nan")))
            nuc_str = (f"{nuc_in:.0f}%" if not np.isnan(nuc_in)
                          else "n/a (NEAT1/MALAT1 absent)")
            cyto_str = (f"{cyto_out:.0f}%" if not np.isnan(cyto_out)
                           else "n/a")
            ax.set_title(
                f"{roi.roi_id} | compartment overlay | {st.METHOD_DISPLAY.get(m, m)}\n"
                f"nuclear-in-nucleus: {nuc_str} | cyto-extranuclear: {cyto_str}",
                fontsize=9, color="white")
            leg = ax.legend(loc="upper right", fontsize=7, framealpha=0.85,
                              facecolor="black", edgecolor="white")
            for t in leg.get_texts():
                t.set_color("white")
            st.setup_axes(ax, bbox_um=roi.bbox_um)
            _add_dynamic_scale_bar(ax, roi, sample, color="white")
            fig.tight_layout()
            fig.savefig(str(out), dpi=180, bbox_inches="tight",
                          pad_inches=0.05, facecolor="black")
            plt.close(fig)
            out_paths[f"compartment_overlay/{m}"] = out

    # P1.b: compartment_panel_compare.png — multi-panel grid for direct
    # method comparison.  2×3 layout (5 panels + 1 empty) gives larger
    # individual panels than the original 1×5.
    if (len(nuclear_tx) + len(cyto_tx)) > 0:
        out_panel = out_dir / f"{roi.roi_id}_compartment_panel_compare.png"
        if not _png_skip(out_panel, []):
            n = len(methods)
            n_cols = 3 if n > 3 else n
            n_rows = (n + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols,
                                          figsize=(5.2 * n_cols, 5.2 * n_rows + 0.6),
                                          facecolor="black", dpi=160)
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes[None, :]
            elif n_cols == 1:
                axes = axes[:, None]
            for i, m in enumerate(methods):
                r, c = i // n_cols, i % n_cols
                ax = axes[r, c]
                _draw_dapi(ax, dapi, roi, alpha=0.55)
                if len(nuclear_tx):
                    ax.scatter(nuclear_tx["x_location"], nuclear_tx["y_location"],
                                s=4, c=NUC_COLOR, alpha=0.90,
                                edgecolors="white", linewidths=0.15)
                if len(cyto_tx):
                    ax.scatter(cyto_tx["x_location"], cyto_tx["y_location"],
                                s=3, c=CYTO_COLOR, alpha=0.80,
                                edgecolors="white", linewidths=0.10)
                if nuc_polys is not None:
                    _draw_polygons_outline(ax, nuc_polys, color="#FFFFFF",
                                             lw=0.3, alpha=0.50, halo=False)
                res_m = load_method(m, sample, roi)
                if res_m["available"]:
                    _draw_polygons_outline(ax, res_m["polygons"],
                                             color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                             lw=0.9, alpha=0.95, halo=False)
                nuc_in, cyto_out = method_stats.get(m, (float("nan"), float("nan")))
                ax.set_title(
                    f"{st.METHOD_DISPLAY.get(m, m)}\n"
                    f"nuc-in-nuc: {nuc_in:.0f}% | cyto-extra: {cyto_out:.0f}%",
                    fontsize=10, color="white")
                st.setup_axes(ax, bbox_um=roi.bbox_um)
            for j in range(n, n_rows * n_cols):
                r, c = j // n_cols, j % n_cols
                axes[r, c].set_visible(False)
            warn = ("" if (len(nuclear_tx) and len(cyto_tx))
                     else "  ⚠ only one marker class present in panel — comparison is one-sided")
            fig.suptitle(
                f"{roi.roi_id} | compartment per-method comparison "
                f"(blue=nuclear marker, orange=cyto marker, white=nucleus boundary)"
                f"{warn}",
                fontsize=11, color="white")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(str(out_panel), dpi=180, bbox_inches="tight",
                          pad_inches=0.05, facecolor="black")
            plt.close(fig)
            out_paths["compartment_panel_compare"] = out_panel

    # P2.a: nucleus_distance_split — split histogram by marker class.
    if "nucleus_distance" in tx.columns and (len(nuclear_tx) + len(cyto_tx)) > 0:
        out = out_dir / f"{roi.roi_id}_nucleus_distance_split.png"
        if not _png_skip(out, []):
            fig, ax = st.new_figure(figsize=(6.0, 3.6))
            n_d = nuclear_tx["nucleus_distance"].dropna() if len(nuclear_tx) else pd.Series([], dtype=float)
            c_d = cyto_tx["nucleus_distance"].dropna() if len(cyto_tx) else pd.Series([], dtype=float)
            bins = np.linspace(0, max(40, float(tx["nucleus_distance"].quantile(0.99))), 40)
            if len(n_d):
                ax.hist(n_d, bins=bins, color=NUC_COLOR, alpha=0.65,
                          density=True, label=f"nuclear (n={len(n_d):,}, median={n_d.median():.1f}µm)")
            if len(c_d):
                ax.hist(c_d, bins=bins, color=CYTO_COLOR, alpha=0.55,
                          density=True, label=f"cyto (n={len(c_d):,}, median={c_d.median():.1f}µm)")
            ax.axvline(5, color="red", linestyle="--", lw=0.8)
            ax.axvline(10, color="orange", linestyle="--", lw=0.8)
            ax.set_xlabel("nucleus distance (µm)")
            ax.set_ylabel("density")
            sep = ""
            if len(n_d) and len(c_d):
                shift = c_d.median() - n_d.median()
                sep = f"  Δmedian = {shift:+.1f}µm (cyto − nuclear)"
            elif len(c_d):
                sep = "  ⚠ nuclear class absent (NEAT1/MALAT1 not in panel)"
            elif len(n_d):
                sep = "  ⚠ cyto class absent"
            ax.set_title(f"{roi.roi_id} | tx → nucleus distance, by marker class{sep}",
                            fontsize=9)
            ax.legend(fontsize=7)
            st.setup_axes(ax, keep_axis=True, equal=False)
            out_paths["nucleus_distance_split"] = _save(fig, out)

    # P2.b: compartment_mismatch_{method}.png — highlight cells with poor
    # compartment fidelity (high nuc-marker outside nucleus / high
    # cyto-marker stuck inside nucleus).
    # Threshold is computed from the POOLED distribution of all 5 methods'
    # per-cell mismatch scores so that "red" means the same thing in every
    # PNG → cross-method comparison is direct.
    if (len(nuclear_tx) + len(cyto_tx)) > 0 and "overlaps_nucleus" in tx.columns:
        try:
            import geopandas as gpd
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Polygon as MplPolygon
        except Exception:
            gpd = None
        if gpd is not None:
            n_gdf = gpd.GeoDataFrame(
                {"on": nuclear_tx["overlaps_nucleus"].fillna(0).astype(int).values},
                geometry=gpd.points_from_xy(nuclear_tx["x_location"], nuclear_tx["y_location"]),
                crs=None) if len(nuclear_tx) else None
            c_gdf = gpd.GeoDataFrame(
                {"on": cyto_tx["overlaps_nucleus"].fillna(0).astype(int).values},
                geometry=gpd.points_from_xy(cyto_tx["x_location"], cyto_tx["y_location"]),
                crs=None) if len(cyto_tx) else None

            # Pre-compute per-cell mismatch for every method.
            method_cells: Dict[str, "gpd.GeoDataFrame"] = {}
            for m in methods:
                res_m = load_method(m, sample, roi)
                if not res_m.get("available") or len(res_m["polygons"]) == 0:
                    continue
                cells_gdf = res_m["polygons"][res_m["polygons"].geometry.notna()].copy()
                cells_gdf = cells_gdf.reset_index(drop=True)
                cells_gdf["cell_idx"] = cells_gdf.index
                nuc_extranuc = pd.Series(np.nan, index=cells_gdf.index)
                if n_gdf is not None and len(n_gdf):
                    j = gpd.sjoin(n_gdf, cells_gdf[["cell_idx", "geometry"]],
                                    predicate="within", how="inner")
                    if len(j):
                        grp = j.groupby("cell_idx")["on"].agg(
                            ['count', lambda x: (x == 0).sum()])
                        nuc_extranuc.loc[grp.index] = grp.iloc[:, 1] / grp["count"]
                cyto_intranuc = pd.Series(np.nan, index=cells_gdf.index)
                if c_gdf is not None and len(c_gdf):
                    j = gpd.sjoin(c_gdf, cells_gdf[["cell_idx", "geometry"]],
                                    predicate="within", how="inner")
                    if len(j):
                        grp = j.groupby("cell_idx")["on"].agg(
                            ['count', lambda x: (x == 1).sum()])
                        cyto_intranuc.loc[grp.index] = grp.iloc[:, 1] / grp["count"]
                stack = np.vstack([nuc_extranuc.values, cyto_intranuc.values])
                with np.errstate(invalid="ignore"):
                    cells_gdf["mismatch"] = np.nanmean(stack, axis=0)
                method_cells[m] = cells_gdf

            # Global thresholds from the pooled per-cell mismatch values.
            pooled = pd.concat([cdf["mismatch"] for cdf in method_cells.values()],
                                  ignore_index=True).dropna()
            if len(pooled):
                gq75 = float(pooled.quantile(0.75))
                gq50 = float(pooled.quantile(0.50))
            else:
                gq75 = gq50 = float("nan")

            for m in methods:
                if m not in method_cells:
                    continue
                out = out_dir / f"{roi.roi_id}_compartment_mismatch_{m}.png"
                if _png_skip(out, []):
                    continue
                cells_gdf = method_cells[m]

                # Render: red fill = pooled-top-quartile, orange = pooled-top-half,
                # method-color outline = bottom half.
                fig, ax = plt.subplots(figsize=(5.6, 5.4), facecolor="black", dpi=180)
                _draw_dapi(ax, dapi, roi, alpha=0.55)

                q75 = gq75
                q50 = gq50

                bad_p, mid_p, ok_p = [], [], []
                for _, row in cells_gdf.iterrows():
                    g = row["geometry"]
                    if g is None or g.is_empty or g.geom_type != "Polygon":
                        continue
                    xs, ys = zip(*list(g.exterior.coords))
                    poly = MplPolygon(np.column_stack([xs, ys]), closed=True)
                    mm = row["mismatch"]
                    if np.isnan(mm):
                        ok_p.append(poly)
                    elif not np.isnan(q75) and mm >= q75:
                        bad_p.append(poly)
                    elif not np.isnan(q50) and mm >= q50:
                        mid_p.append(poly)
                    else:
                        ok_p.append(poly)
                if ok_p:
                    ax.add_collection(PatchCollection(
                        ok_p, facecolors="none",
                        edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                        linewidths=0.5, alpha=0.85))
                if mid_p:
                    ax.add_collection(PatchCollection(
                        mid_p, facecolors="#FFA50066",
                        edgecolors="#FFA500", linewidths=1.0, alpha=0.85))
                if bad_p:
                    ax.add_collection(PatchCollection(
                        bad_p, facecolors="#FF000099",
                        edgecolors="#FF0000", linewidths=1.4, alpha=0.95))

                # Reference markers (small) so user sees what drove the score
                if len(nuclear_tx):
                    ax.scatter(nuclear_tx["x_location"], nuclear_tx["y_location"],
                                s=2, c=NUC_COLOR, alpha=0.55,
                                edgecolors="none")
                if len(cyto_tx):
                    ax.scatter(cyto_tx["x_location"], cyto_tx["y_location"],
                                s=1.5, c=CYTO_COLOR, alpha=0.50,
                                edgecolors="none")

                # Per-method % of cells in each bucket → useful in title.
                valid_m = cells_gdf["mismatch"].dropna()
                pct_bad = (100.0 * (valid_m >= q75).sum() / len(valid_m)
                           if len(valid_m) and not np.isnan(q75) else float("nan"))
                from matplotlib.patches import Patch
                handles = [
                    Patch(facecolor="#FF000099", edgecolor="#FF0000",
                          label=(f"top-25% pooled (≥{q75:.2f})  this method: {pct_bad:.0f}%"
                                  if not np.isnan(q75) else "top-25% pooled")),
                    Patch(facecolor="#FFA50066", edgecolor="#FFA500",
                          label=f"50-75% pooled (≥{q50:.2f})" if not np.isnan(q50) else "mid"),
                    plt.Line2D([], [], color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                lw=1.5, label="bottom 50% (good)"),
                ]
                leg = ax.legend(handles=handles, loc="upper right", fontsize=7,
                                framealpha=0.85, facecolor="black",
                                edgecolor="white")
                for t in leg.get_texts():
                    t.set_color("white")
                if len(nuclear_tx) and len(cyto_tx):
                    basis = "½(nuclear-extranuclear + cyto-intranuclear)"
                elif len(cyto_tx):
                    basis = "cyto-marker fraction stuck inside nucleus only ⚠ nuclear class absent"
                else:
                    basis = "nuclear-marker fraction outside nucleus ⚠ cyto class absent"
                ax.set_title(
                    f"{roi.roi_id} | compartment mismatch | "
                    f"{st.METHOD_DISPLAY.get(m, m)}\n"
                    f"score = {basis}  |  red = pooled top-25% across all 5 methods",
                    fontsize=8, color="white")
                st.setup_axes(ax, bbox_um=roi.bbox_um)
                _add_dynamic_scale_bar(ax, roi, sample, color="white")
                fig.tight_layout()
                fig.savefig(str(out), dpi=180, bbox_inches="tight",
                              pad_inches=0.05, facecolor="black")
                plt.close(fig)
                out_paths[f"compartment_mismatch/{m}"] = out

    return out_paths


# ---------------------------------------------------------------------------
# 4. low_vsi: whole-tissue VSI map with ROI box + z-colored scatter
# ---------------------------------------------------------------------------
def render_low_vsi_addons(sample: SampleData, roi: ROIRecord,
                              baseline: str, methods: Sequence[str],
                              outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "low_vsi")

    # ROI VSI crop
    vsi = crop_vsi(sample, roi)
    out = out_dir / f"{roi.roi_id}_vsi_crop.png"
    if vsi.get("map") is not None and not _png_skip(out, []):
        fig, ax = st.new_figure(figsize=(5, 5))
        ext = vsi["extent_um"]
        im = ax.imshow(vsi["map"], cmap="viridis",
                         extent=(ext[0], ext[1], ext[2], ext[3]),
                         vmin=0, vmax=1, interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label("VSI / coherence", fontsize=8)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | VSI crop")
        _add_dynamic_scale_bar(ax, roi, sample, color="white")
        out_paths["vsi_crop"] = _save(fig, out)

    # z-colored transcript scatter
    tx = crop_transcripts(sample, roi)
    if len(tx) and "z_location" in tx.columns:
        out = out_dir / f"{roi.roi_id}_zscatter.png"
        if not _png_skip(out, []):
            fig, ax = st.new_figure(figsize=(5.4, 5))
            sc = ax.scatter(tx["x_location"], tx["y_location"],
                              c=tx["z_location"], cmap="coolwarm",
                              s=0.5, alpha=0.7, rasterized=True)
            cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label("z (µm)", fontsize=8)
            st.setup_axes(ax, bbox_um=roi.bbox_um,
                            title=f"{roi.roi_id} | z-colored transcripts")
            _add_dynamic_scale_bar(ax, roi, sample, color="black")
            out_paths["zscatter"] = _save(fig, out)
        # x-z side view
        out_xz = out_dir / f"{roi.roi_id}_sideview_xz.png"
        if not _png_skip(out_xz, []):
            fig, ax = st.new_figure(figsize=(5.4, 3))
            ax.scatter(tx["x_location"], tx["z_location"],
                        s=0.3, c="#3478F6", alpha=0.5, rasterized=True)
            ax.set_xlabel("x (µm)"); ax.set_ylabel("z (µm)")
            ax.set_title(f"{roi.roi_id} | x-z side view")
            st.setup_axes(ax, keep_axis=True, equal=False)
            out_paths["sideview_xz"] = _save(fig, out_xz)

    # ★ low-VSI overlay × per method — directly answers "어디가 low-VSI 이고
    # 그 영역에 각 method 의 cell 이 어떻게 대응했나" (paper §5.4 진단).
    vsi_overlay_stats = {}  # method -> (n_in_low, n_total)
    if vsi.get("map") is not None:
        ext = vsi["extent_um"]
        vsi_arr = vsi["map"]
        for m in methods:
            out_m = out_dir / f"{roi.roi_id}_vsi_overlay_{m}.png"
            res = load_method(m, sample, roi)

            # Always compute stats (used by summary chart even if PNG skipped)
            low = (vsi_arr < 0.5) & np.isfinite(vsi_arr)
            if res["available"] and len(res["centroids"]):
                cents = res["centroids"]
                px_um = (ext[1] - ext[0]) / vsi_arr.shape[1]
                py_um = abs(ext[2] - ext[3]) / vsi_arr.shape[0]
                cols = ((cents["x_centroid_um"] - ext[0]) / px_um
                        ).astype(int).clip(0, vsi_arr.shape[1] - 1)
                rows = ((cents["y_centroid_um"] - min(ext[2], ext[3]))
                        / py_um).astype(int).clip(0, vsi_arr.shape[0] - 1)
                in_low = low[rows, cols]
                vsi_overlay_stats[m] = (int(in_low.sum()), int(len(cents)))

            if _png_skip(out_m, []):
                continue
            fig, ax = st.new_figure(figsize=(5.4, 5))
            dapi = crop_dapi(sample, roi)
            _draw_dapi(ax, dapi, roi, alpha=0.7)

            # Magenta semi-transparent overlay where VSI < 0.5 (= problematic).
            rgba = np.zeros((*vsi_arr.shape, 4), dtype=np.float32)
            rgba[low, 0] = 1.0
            rgba[low, 1] = 0.2
            rgba[low, 2] = 0.85
            rgba[low, 3] = 0.55
            ax.imshow(rgba, extent=(ext[0], ext[1], ext[2], ext[3]),
                       interpolation="nearest")

            if res["available"]:
                _draw_polygons_outline(
                    ax, res["polygons"],
                    color=st.METHOD_COLORS.get(m, "#00FFFF"),
                    lw=0.7, alpha=0.95, halo=True)
                # Red X markers for centroids in low-VSI (visual highlight).
                if m in vsi_overlay_stats and vsi_overlay_stats[m][0]:
                    cents = res["centroids"]
                    px_um = (ext[1] - ext[0]) / vsi_arr.shape[1]
                    py_um = abs(ext[2] - ext[3]) / vsi_arr.shape[0]
                    cols = ((cents["x_centroid_um"] - ext[0]) / px_um
                            ).astype(int).clip(0, vsi_arr.shape[1] - 1)
                    rows = ((cents["y_centroid_um"] - min(ext[2], ext[3]))
                            / py_um).astype(int).clip(0, vsi_arr.shape[0] - 1)
                    in_low_arr = low[rows, cols]
                    bad = cents.iloc[in_low_arr.values] if hasattr(in_low_arr, 'values') else cents[in_low_arr]
                    ax.scatter(bad["x_centroid_um"], bad["y_centroid_um"],
                                s=20, marker="x", c="#FF3B30",
                                linewidths=1.4, alpha=0.95)

            # Spatial image legend (no numerical inset).
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=(1.0, 0.2, 0.85, 0.55),
                              edgecolor="none", label="VSI < 0.5 (low-VSI)"),
                        plt.Line2D([], [], color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                    lw=2.0, label=f"{st.METHOD_DISPLAY.get(m, m)} boundary")]
            ax.legend(handles=handles, loc="upper right", fontsize=7,
                       framealpha=0.85)

            st.setup_axes(ax, bbox_um=roi.bbox_um,
                           title=f"{roi.roi_id} | low-VSI overlay | {m}")
            _add_dynamic_scale_bar(ax, roi, sample, color="white")
            out_paths[f"vsi_overlay/{m}"] = _save(fig, out_m)

        # ----- Companion summary bar chart (separate from spatial PNGs) -----
        if vsi_overlay_stats:
            out_sum = out_dir / f"{roi.roi_id}_vsi_overlay_summary.png"
            if not _png_skip(out_sum, []):
                fig, ax = st.new_figure(figsize=(6.5, 4.0))
                ax.set_facecolor("white")
                method_keys = [m for m in methods if m in vsi_overlay_stats]
                fracs = [100 * vsi_overlay_stats[m][0] / max(1, vsi_overlay_stats[m][1])
                         for m in method_keys]
                colors = [st.METHOD_COLORS.get(m, "#888") for m in method_keys]
                labels = [st.METHOD_DISPLAY.get(m, m) for m in method_keys]
                ax.bar(range(len(method_keys)), fracs, color=colors,
                        edgecolor="black", linewidth=0.5)
                ax.set_xticks(range(len(method_keys)))
                ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
                ax.set_ylabel(
                    "% cells with centroid in low-VSI region\n(VSI < 0.5)",
                    fontsize=9)
                ax.set_title(
                    f"{roi.roi_id} | per-method low-VSI invasion\n"
                    "(higher = more cells placed in problematic vertical-overlap area)",
                    fontsize=10)
                for i, (m, f) in enumerate(zip(method_keys, fracs)):
                    n_in, n_t = vsi_overlay_stats[m]
                    ax.text(i, f + 0.5, f"{n_in}/{n_t}\n({f:.1f}%)",
                              ha="center", va="bottom", fontsize=7)
                ax.set_ylim(0, max(20, max(fracs) * 1.25) if fracs else 20)
                for s in ax.spines.values():
                    s.set_color("#888")
                fig.tight_layout()
                fig.savefig(str(out_sum), dpi=200, bbox_inches="tight",
                              pad_inches=0.05, facecolor="white")
                plt.close(fig)
                out_paths["vsi_overlay_summary"] = out_sum

    # ★ z-color + boundary diagnostic — paper §5.4 의 진짜 검증.
    # 같은 cell 안 transcript 의 z 값이 단일 층(정상)인지, 위·아래 섞임
    # (vertical overlap 으로 인한 false-merge)인지 시각적으로 보여줌.
    z_mixing_per_method = {}  # method -> (n_total_cells, n_mixing_cells)
    if len(tx) and "z_location" in tx.columns:
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon as MplPolygon

        z_min = float(tx["z_location"].quantile(0.02))
        z_max = float(tx["z_location"].quantile(0.98))

        # ----- Determine z-range threshold from BASELINE (xenium_nucleus) -----
        # Strategy: compute per-cell z_range distribution of the baseline
        # method, take its 90th percentile as a "single-layer max z-spread"
        # bound, then flag cells in OTHER methods that exceed this bound.
        # This makes the threshold:
        #   • data-aware (fits actual tissue thickness),
        #   • cross-method comparable (same threshold for all),
        #   • baseline-relative (over-baseline cells = potential false-merge).
        baseline_method = "xenium_nucleus"
        Z_RANGE_THRESHOLD = None
        try:
            base_res = load_method(baseline_method, sample, roi)
            if base_res["available"] and base_res["raster"] is not None:
                br = base_res["raster"]
                bx0, by0, bx1, by1 = base_res["raster_extent_um"]
                Hb, Wb = br.shape
                bppx = Wb / max(1e-6, bx1 - bx0)
                bppy = Hb / max(1e-6, by1 - by0)
                bcols = ((tx["x_location"].values - bx0) * bppx).astype(int)
                brows = ((tx["y_location"].values - by0) * bppy).astype(int)
                bin = ((bcols >= 0) & (bcols < Wb) &
                       (brows >= 0) & (brows < Hb))
                blab = np.full(len(tx), -1, dtype=np.int64)
                blab[bin] = br[brows[bin], bcols[bin]]
                bdf = pd.DataFrame({"cell": blab, "z": tx["z_location"].values})
                bdf = bdf[bdf["cell"] > 0]
                if len(bdf):
                    base_z_ranges = (bdf.groupby("cell")["z"]
                                       .agg(lambda v: float(v.max() - v.min())))
                    if len(base_z_ranges):
                        Z_RANGE_THRESHOLD = float(np.percentile(base_z_ranges, 90))
        except Exception:
            pass
        if Z_RANGE_THRESHOLD is None or Z_RANGE_THRESHOLD < 1.0:
            Z_RANGE_THRESHOLD = 5.0  # fallback

        for m in methods:
            out_m = out_dir / f"{roi.roi_id}_z_color_boundary_{m}.png"
            if _png_skip(out_m, []):
                continue
            res = load_method(m, sample, roi)
            if not res["available"] or res["raster"] is None:
                continue

            # Compute per-cell z-range using the per-method raster
            raster = res["raster"]
            r_x0, r_y0, r_x1, r_y1 = res["raster_extent_um"]
            H, W = raster.shape
            ppx = W / max(1e-6, (r_x1 - r_x0))
            ppy = H / max(1e-6, (r_y1 - r_y0))
            cols = ((tx["x_location"].values - r_x0) * ppx).astype(int)
            rows = ((tx["y_location"].values - r_y0) * ppy).astype(int)
            inside = (cols >= 0) & (cols < W) & (rows >= 0) & (rows < H)
            cell_label = np.full(len(tx), -1, dtype=np.int64)
            cell_label[inside] = raster[rows[inside], cols[inside]]
            df = pd.DataFrame({
                "cell": cell_label,
                "z": tx["z_location"].values,
            })
            df = df[df["cell"] > 0]
            if len(df) == 0:
                continue
            z_stats = df.groupby("cell")["z"].agg([
                ("z_range", lambda v: float(v.max() - v.min())),
                ("z_std",   "std"),
                ("n",       "count"),
            ]).reset_index()
            n_total = int(len(z_stats))
            # Use baseline-derived threshold (top 10 % of 10x cell z_ranges).
            # A method's cell flagged as z-mixing iff its intra-cell z_range
            # exceeds what 90 % of baseline single-layer cells exhibit.
            mixing_mask = z_stats["z_range"] > Z_RANGE_THRESHOLD
            n_mixing = int(mixing_mask.sum())
            mixing_labels = set(z_stats.loc[mixing_mask, "cell"].astype(int))
            # Save full state for the multi-panel comparison plot below.
            z_mixing_per_method[m] = {
                "n_total": n_total,
                "n_mixing": n_mixing,
                "mixing_labels": mixing_labels,
                "polygons": res["polygons"],
            }

            # ----- draw -----
            fig, ax = st.new_figure(figsize=(5.6, 5.2))
            dapi = crop_dapi(sample, roi)
            _draw_dapi(ax, dapi, roi, alpha=0.30)
            sc = ax.scatter(tx["x_location"], tx["y_location"],
                              c=tx["z_location"], cmap="coolwarm",
                              vmin=z_min, vmax=z_max,
                              s=0.6, alpha=0.85, rasterized=True,
                              edgecolors="none")
            cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
            cb.set_label("z (µm)", fontsize=8)

            # Boundary in two passes: clean cells in method color (thin),
            # z-mixing cells in red (thick) — visually highlights suspicious.
            polys = res["polygons"]
            if "raster_label" in polys.columns:
                clean_polys = polys[~polys["raster_label"].isin(mixing_labels)]
                mix_polys = polys[polys["raster_label"].isin(mixing_labels)]
                _draw_polygons_outline(ax, clean_polys,
                                         color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                         lw=0.5, alpha=0.85, halo=False)
                _draw_polygons_outline(ax, mix_polys,
                                         color="#FF1493",  # deep pink for mixing
                                         lw=1.4, alpha=0.95, halo=True)
            else:
                _draw_polygons_outline(ax, polys,
                                         color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                         lw=0.6, alpha=0.85, halo=False)

            # No text inset on spatial image — numerical summary moved to
            # separate z_mixing_comparison.png (per-ROI bar chart).
            frac = 100 * n_mixing / max(1, n_total)

            # Legend
            handles = [
                plt.Line2D([], [], color=st.METHOD_COLORS.get(m, "#00FFFF"),
                            lw=2.0,
                            label=f"{st.METHOD_DISPLAY.get(m, m)} (clean)"),
                plt.Line2D([], [], color="#FF1493",
                            lw=2.0,
                            label=f"z-mixing (range > {Z_RANGE_THRESHOLD:.0f} µm)"),
            ]
            ax.legend(handles=handles, loc="upper right", fontsize=7,
                       framealpha=0.85)
            st.setup_axes(ax, bbox_um=roi.bbox_um,
                           title=f"{roi.roi_id} | z-color + boundary | "
                                 f"{st.METHOD_DISPLAY.get(m, m)}")
            _add_dynamic_scale_bar(ax, roi, sample, color="white")
            out_paths[f"z_color_boundary/{m}"] = _save(fig, out_m)

        # ----- Cross-method comparison bar chart -----
        if z_mixing_per_method:
            out_cmp = out_dir / f"{roi.roi_id}_z_mixing_comparison.png"
            if not _png_skip(out_cmp, []):
                fig, ax = st.new_figure(figsize=(6.5, 4.0))
                ax.set_facecolor("white")
                method_keys = [m for m in methods if m in z_mixing_per_method]
                fracs = [
                    100 * z_mixing_per_method[m]["n_mixing"]
                    / max(1, z_mixing_per_method[m]["n_total"])
                    for m in method_keys
                ]
                colors = [st.METHOD_COLORS.get(m, "#888") for m in method_keys]
                labels = [st.METHOD_DISPLAY.get(m, m) for m in method_keys]
                ax.bar(range(len(method_keys)), fracs, color=colors,
                        edgecolor="black", linewidth=0.5)
                ax.set_xticks(range(len(method_keys)))
                ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
                ax.set_ylabel(
                    f"% cells with z-mixing\n"
                    f"(intra-cell z_range > {Z_RANGE_THRESHOLD:.1f} µm,\n"
                    f"baseline 10x p90)",
                    fontsize=9)
                ax.set_title(
                    f"{roi.roi_id} | per-method z-mixing comparison\n"
                    f"(higher = more vertical-overlap false-merges)",
                    fontsize=10)
                for i, m in enumerate(method_keys):
                    n_t = z_mixing_per_method[m]["n_total"]
                    n_m = z_mixing_per_method[m]["n_mixing"]
                    f = fracs[i]
                    ax.text(i, f + 0.5, f"{n_m}/{n_t}\n({f:.1f}%)",
                              ha="center", va="bottom", fontsize=7)
                ax.set_ylim(0, max(20, max(fracs) * 1.25))
                for s in ax.spines.values():
                    s.set_color("#888")
                fig.tight_layout()
                fig.savefig(str(out_cmp), dpi=200, bbox_inches="tight",
                              pad_inches=0.05, facecolor="white")
                plt.close(fig)
                out_paths["z_mixing_comparison"] = out_cmp

        # ----- ★ Single overlay: z-mixing cells of every method, color-coded ----
        # 사용자 요청: 한 그림 / 잡음 제거 / low-VSI 영역(분홍) 같이 보이되
        # 배경은 제외 / (1)(2) 가 (3)(4) 밑에 묻히지 않게.
        if z_mixing_per_method:
            method_keys = [m for m in methods if m in z_mixing_per_method]
            if method_keys:
                out_panel = out_dir / f"{roi.roi_id}_z_mixing_panel_compare.png"
                if not _png_skip(out_panel, []):
                    fig, ax = plt.subplots(figsize=(8.0, 8.0),
                                              facecolor="black", dpi=200)
                    dapi = crop_dapi(sample, roi)
                    _draw_dapi(ax, dapi, roi, alpha=0.60)

                    # ---- low-VSI pink overlay (VSI<0.5).  NaN pixels = no
                    # ovrlpy support there → naturally excludes off-tissue
                    # background.  Optional very-loose DAPI gate to drop pixels
                    # that are also below noise floor.
                    has_low_vsi = False
                    if vsi.get("map") is not None:
                        try:
                            vsi_arr = vsi["map"]
                            ext = vsi["extent_um"]
                            low_mask = (vsi_arr < 0.5) & np.isfinite(vsi_arr)
                            dapi_img = dapi.get("image") if dapi else None
                            if dapi_img is not None and dapi_img.size:
                                from skimage.transform import resize as _resize
                                dapi_rs = _resize(dapi_img, vsi_arr.shape,
                                                     anti_aliasing=True,
                                                     preserve_range=True)
                                # Otsu threshold splits bimodal tissue vs.
                                # off-tissue background reliably (works on
                                # both dense ROI1 and sparse-tissue ROI2).
                                # Light morphological closing fills small
                                # holes inside tissue.
                                try:
                                    from skimage.filters import threshold_otsu
                                    from skimage.morphology import (
                                        binary_closing, disk)
                                    if float(dapi_rs.max()) > float(dapi_rs.min()):
                                        thr = threshold_otsu(dapi_rs)
                                    else:
                                        thr = float(dapi_rs.min())
                                    tissue = dapi_rs > thr
                                    # Bridge inter-cell gaps so pixel-scale
                                    # low-VSI INSIDE tissue isn't lost, but
                                    # don't remove_small_objects (would also
                                    # wipe pixel-scale tissue islands).
                                    tissue = binary_closing(tissue, disk(2))
                                except Exception:
                                    p50 = float(np.percentile(dapi_rs, 50))
                                    tissue = dapi_rs > p50
                                low_mask = low_mask & tissue
                                # Dilate inside tissue so single-pixel low-VSI
                                # hotspots become visually prominent blobs
                                # (~6µm radius) without bleeding into the
                                # excluded background.
                                try:
                                    from scipy.ndimage import binary_dilation
                                    low_mask = (binary_dilation(
                                        low_mask, iterations=3) & tissue)
                                except Exception:
                                    pass
                            if low_mask.any():
                                # ★ Vivid hot pink at higher alpha so it
                                # visually dominates the cell fills (which
                                # use alpha=0.32).  Now safe to use pink —
                                # (3b) optimal_expansion was moved to violet
                                # in style.py so colors don't collide.
                                rgba = np.zeros((*vsi_arr.shape, 4),
                                                  dtype=np.float32)
                                rgba[low_mask, 0] = 1.0
                                rgba[low_mask, 1] = 0.05
                                rgba[low_mask, 2] = 0.65
                                rgba[low_mask, 3] = 0.85
                                ax.imshow(rgba,
                                            extent=(ext[0], ext[1],
                                                     ext[2], ext[3]),
                                            origin="upper",
                                            interpolation="nearest",
                                            zorder=2)
                                has_low_vsi = True
                        except Exception as e:
                            logger.debug(f"low-VSI overlay skipped: {e}")

                    from matplotlib.collections import PatchCollection
                    from matplotlib.patches import Polygon as MplPolygon

                    # Draw order: largest fill first → smallest outline last so
                    # nucleus-only methods (1)/(2) remain visible on top.
                    # (2) Cellpose drawn dashed so it differentiates from (1)
                    # 10x when both target the same nucleus.
                    fill_methods = {"baysor", "rigid_expansion",
                                       "optimal_expansion"}
                    line_styles = {"cellpose_nuclei": "--",
                                       "xenium_nucleus": "-"}
                    draw_order = ["baysor", "rigid_expansion",
                                    "optimal_expansion", "xenium_nucleus",
                                    "cellpose_nuclei"]
                    draw_order = [m for m in draw_order if m in method_keys]

                    legend_entries = []  # (artist_kind, color, label, style)
                    for zo, m in enumerate(draw_order):
                        polys = z_mixing_per_method[m]["polygons"]
                        mixing_labels = z_mixing_per_method[m]["mixing_labels"]
                        if "raster_label" not in polys.columns or not len(polys):
                            continue
                        mix_patches = []
                        for _, row in polys.iterrows():
                            g = row["geometry"]
                            if g is None or g.is_empty:
                                continue
                            if g.geom_type != "Polygon":
                                continue
                            if int(row["raster_label"]) not in mixing_labels:
                                continue
                            xs, ys = zip(*list(g.exterior.coords))
                            mix_patches.append(MplPolygon(
                                np.column_stack([xs, ys]), closed=True))
                        color = st.METHOD_COLORS.get(m, "#FFFFFF")
                        ls = line_styles.get(m, "-")
                        if mix_patches:
                            if m in fill_methods:
                                pc = PatchCollection(mix_patches,
                                                          facecolors=color,
                                                          edgecolors=color,
                                                          linewidths=0.5,
                                                          alpha=0.32,
                                                          zorder=3 + zo)
                                kind = "fill"
                            else:
                                # Nucleus-only — outline only, thick & opaque so
                                # it punches through the larger fills above.
                                pc = PatchCollection(mix_patches,
                                                          facecolors="none",
                                                          edgecolors=color,
                                                          linewidths=2.4,
                                                          linestyles=ls,
                                                          alpha=1.0,
                                                          zorder=10 + zo)
                                kind = "line"
                            ax.add_collection(pc)
                        else:
                            kind = "fill"
                        n_t = z_mixing_per_method[m]["n_total"]
                        n_m = z_mixing_per_method[m]["n_mixing"]
                        frac = 100 * n_m / max(1, n_t)
                        legend_entries.append(
                            (kind, color, ls,
                             f"{st.METHOD_DISPLAY.get(m, m)}  "
                             f"{n_m}/{n_t} ({frac:.1f}%)"))

                    from matplotlib.patches import Patch
                    handles = []
                    if has_low_vsi:
                        handles.append(Patch(facecolor=(1.0, 0.05, 0.65, 0.85),
                                                edgecolor="none",
                                                label="VSI < 0.5 (low-VSI)"))
                    for kind, c, ls_, lbl in legend_entries:
                        if kind == "line":
                            handles.append(plt.Line2D([], [], color=c, lw=2.4,
                                                          linestyle=ls_,
                                                          label=lbl))
                        else:
                            handles.append(Patch(facecolor=c, edgecolor=c,
                                                    alpha=0.55, label=lbl))
                    if handles:
                        leg = ax.legend(handles=handles, loc="upper right",
                                            fontsize=8, framealpha=0.96,
                                            facecolor="black",
                                            edgecolor="white",
                                            labelcolor="white")
                        leg.set_zorder(50)
                        for txt in leg.get_texts():
                            txt.set_color("white")

                    ax.set_title(
                        f"{roi.roi_id} | z-mixing cells overlay "
                        f"(z_range > {Z_RANGE_THRESHOLD:.1f} µm) "
                        f"+ low-VSI tissue mask",
                        fontsize=11, color="white")
                    st.setup_axes(ax, bbox_um=roi.bbox_um)
                    fig.tight_layout()
                    fig.savefig(str(out_panel), dpi=200, bbox_inches="tight",
                                  pad_inches=0.05, facecolor="black")
                    plt.close(fig)
                    out_paths["z_mixing_panel_compare"] = out_panel

    return out_paths


# ---------------------------------------------------------------------------
# 5. high_density: nuclei centroids + cell-area map + merge/split highlight
# ---------------------------------------------------------------------------
def render_high_density_addons(sample: SampleData, roi: ROIRecord,
                                  baseline: str, methods: Sequence[str],
                                  outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "high_density")

    # Nucleus centroid map (xenium baseline)
    out = out_dir / f"{roi.roi_id}_nucleus_centroids.png"
    if not _png_skip(out, []):
        dapi = crop_dapi(sample, roi)
        base = load_method(baseline, sample, roi)
        fig, ax = st.new_figure(figsize=(8.0, 8.0))
        _draw_dapi(ax, dapi, roi, alpha=0.7)
        if base["available"] and len(base["centroids"]):
            ax.scatter(base["centroids"]["x_centroid_um"],
                        base["centroids"]["y_centroid_um"],
                        s=8, c="#FFD700", edgecolors="black",
                        linewidths=0.3, alpha=0.95)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | nucleus centroids ({baseline})")
        _add_dynamic_scale_bar(ax, roi, sample, color="white")
        fig.savefig(str(out), dpi=200, bbox_inches=None, pad_inches=0.05,
                       facecolor="white", edgecolor="none")
        plt.close(fig)
        out_paths["nucleus_centroids"] = out

    # Merge / split highlight: per non-baseline method
    base = load_method(baseline, sample, roi)
    for m in methods:
        if m == baseline:
            continue
        out = out_dir / f"{roi.roi_id}_merge_split_{m}.png"
        if _png_skip(out, []):
            continue
        re = compute_reassignment(sample, roi, baseline, m)
        if not re["available"] or not len(re["frame"]):
            continue
        df = re["frame"]
        # Identify merges (multiple baseline cells → one method cell) and
        # splits (one baseline cell → multiple method cells).
        merge_groups = (df.groupby("cell_id_b")["cell_id_a"]
                          .nunique().reset_index(name="n_a"))
        merge_cells = set(merge_groups[merge_groups["n_a"] > 1]["cell_id_b"])
        split_groups = (df.groupby("cell_id_a")["cell_id_b"]
                          .nunique().reset_index(name="n_b"))
        split_cells = set(split_groups[split_groups["n_b"] > 1]["cell_id_a"])

        df["highlight"] = "neutral"
        df.loc[df["cell_id_b"].isin(merge_cells), "highlight"] = "merged"
        df.loc[df["cell_id_a"].isin(split_cells), "highlight"] = "split"

        fig, ax = st.new_figure(figsize=(8.0, 8.0))
        dapi = crop_dapi(sample, roi)
        _draw_dapi(ax, dapi, roi, alpha=0.45)
        for cat, color in [("neutral", "#BFBFBF"),
                            ("merged", "#FF3B30"),
                            ("split", "#34C759")]:
            sub = df[df["highlight"] == cat]
            if not len(sub):
                continue
            ax.scatter(sub["x_location"], sub["y_location"],
                        s=0.4 if cat == "neutral" else 0.7,
                        alpha=0.4 if cat == "neutral" else 0.85,
                        c=color, edgecolors=None,
                        rasterized=True, label=f"{cat} ({len(sub):,})")
        ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
        st.setup_axes(ax, bbox_um=roi.bbox_um,
                        title=f"{roi.roi_id} | merge/split | {baseline} vs {m}")
        _add_dynamic_scale_bar(ax, roi, sample, color="white")
        fig.savefig(str(out), dpi=200, bbox_inches=None, pad_inches=0.05,
                       facecolor="white", edgecolor="none")
        plt.close(fig)
        out_paths[f"merge_split/{m}"] = out
    return out_paths


# ---------------------------------------------------------------------------
# 6. fold_boundary: wider context + structure marker map
# ---------------------------------------------------------------------------
def render_fold_boundary_addons(sample: SampleData, roi: ROIRecord,
                                    baseline: str, methods: Sequence[str],
                                    outdir: Path) -> Dict[str, Path]:
    out_paths: Dict[str, Path] = {}
    out_dir = _class_dir(outdir, roi.roi_id, "fold_boundary")

    # Wider context: bbox grown 2x
    x0, y0, x1, y1 = roi.bbox_um
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    half = max(x1 - x0, y1 - y0)
    wide = ROIRecord(
        roi_id=roi.roi_id + "__wide", sample_id=roi.sample_id,
        roi_class=roi.roi_class,
        bbox_um=(cx - half, cy - half, cx + half, cy + half),
        source_map=roi.source_map, hypothesis=roi.hypothesis,
        primary_reference=roi.primary_reference,
    )
    out = out_dir / f"{roi.roi_id}_wider_context.png"
    if not _png_skip(out, [sample.dapi_path] if sample.dapi_path else []):
        dapi = crop_dapi(sample, wide)
        fig, ax = st.new_figure(figsize=(5, 5))
        _draw_dapi(ax, dapi, wide)
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                  fill=False, edgecolor="#FFD700", linewidth=1.4))
        st.setup_axes(ax, bbox_um=wide.bbox_um,
                        title=f"{roi.roi_id} | wider context (yellow = ROI)")
        _add_dynamic_scale_bar(ax, roi, sample, color="white",
                                  bbox_um=wide.bbox_um)
        out_paths["wider_context"] = _save(fig, out)

    # Structure marker (ad_pathology / endothelial)
    tx = crop_transcripts(sample, roi)
    if len(tx):
        struct_genes = (DEFAULT_BRAIN_MARKERS["ad_pathology"] +
                         DEFAULT_BRAIN_MARKERS["endothelial"])
        present = [g for g in struct_genes
                    if g in tx["feature_name"].astype(str).unique()]
        if present:
            out = out_dir / f"{roi.roi_id}_structure_marker.png"
            if not _png_skip(out, []):
                pal = marker_color_palette(present)
                fig, ax = st.new_figure(figsize=(5.2, 5))
                dapi = crop_dapi(sample, roi)
                _draw_dapi(ax, dapi, roi, alpha=0.55)
                for g in present:
                    sub = tx[tx["feature_name"].astype(str) == g]
                    ax.scatter(sub["x_location"], sub["y_location"],
                                s=4, c=pal[g], alpha=0.85,
                                edgecolors="white", linewidths=0.15,
                                label=f"{g} ({len(sub):,})")
                ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
                st.setup_axes(ax, bbox_um=roi.bbox_um,
                                title=f"{roi.roi_id} | structure markers")
                _add_dynamic_scale_bar(ax, roi, sample, color="white")
                out_paths["structure_marker"] = _save(fig, out)

    # ------------------------------------------------------------------
    # P1/P2 boundary overflow add-ons (plan/roi_visualization §5.6)
    # ------------------------------------------------------------------
    # Build a tissue mask (DAPI Otsu) in ROI coordinates so we can:
    #  - shade the *off-tissue* area in magenta (the "problem region",
    #    matching the low-VSI pink convention),
    #  - flag cells whose centroid lies off-tissue as overflow (red fill),
    #  - quantify per-method overflow fraction.
    dapi_full = crop_dapi(sample, roi)
    dapi_img = dapi_full.get("image") if dapi_full else None
    tissue_mask = None
    tx0, ty0, tx1, ty1 = roi.bbox_um  # extent of tissue mask
    if dapi_img is not None and dapi_img.size:
        try:
            from skimage.filters import threshold_otsu
            from skimage.morphology import binary_closing, binary_opening, disk
            thr = threshold_otsu(dapi_img) if (dapi_img.max() > dapi_img.min()) else dapi_img.min()
            tm = dapi_img > thr
            tm = binary_closing(tm, disk(4))
            tm = binary_opening(tm, disk(2))
            tissue_mask = tm  # (H, W) boolean
        except Exception as e:
            logger.debug(f"tissue mask failed: {e}")

    def _is_off_tissue(cx_um: float, cy_um: float) -> Optional[bool]:
        if tissue_mask is None:
            return None
        H, W = tissue_mask.shape
        # ROI bbox spans (tx0,ty0)-(tx1,ty1).  DAPI image's pixel (0,0) is at
        # (tx0, ty0); (W-1, H-1) at (tx1, ty1). Standard image origin "upper".
        col = int((cx_um - tx0) / max(tx1 - tx0, 1e-6) * W)
        row = int((cy_um - ty0) / max(ty1 - ty0, 1e-6) * H)
        col = min(max(col, 0), W - 1)
        row = min(max(row, 0), H - 1)
        return not bool(tissue_mask[row, col])

    # P1: per-method overflow stats + render
    method_overflow: Dict[str, Tuple[int, int]] = {}  # m -> (n_overflow, n_total)
    method_cells_with_off: Dict[str, pd.DataFrame] = {}
    for m in methods:
        try:
            res_m = load_method(m, sample, roi)
        except Exception:
            continue
        if not res_m.get("available") or len(res_m["centroids"]) == 0:
            continue
        cents = res_m["centroids"].copy()
        cents["off_tissue"] = [
            _is_off_tissue(float(x), float(y))
            for x, y in zip(cents["x_centroid_um"], cents["y_centroid_um"])
        ]
        n_off = int((cents["off_tissue"] == True).sum())
        n_tot = int(cents["off_tissue"].notna().sum())
        method_overflow[m] = (n_off, n_tot)
        method_cells_with_off[m] = cents

    OFF_COLOR_RGBA = (1.0, 0.08, 0.6, 0.32)  # magenta (problem-region convention)

    def _shade_off_tissue(ax) -> None:
        if tissue_mask is None:
            return
        H, W = tissue_mask.shape
        rgba = np.zeros((H, W, 4), dtype=np.float32)
        off = ~tissue_mask
        rgba[off, 0] = OFF_COLOR_RGBA[0]
        rgba[off, 1] = OFF_COLOR_RGBA[1]
        rgba[off, 2] = OFF_COLOR_RGBA[2]
        rgba[off, 3] = OFF_COLOR_RGBA[3]
        ax.imshow(rgba, extent=(tx0, tx1, ty1, ty0),
                    origin="upper", interpolation="nearest", zorder=2)

    # ----- boundary_overflow_{method}.png × 5
    if tissue_mask is not None and method_overflow:
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon as MplPolygon, Patch
        struct_genes = (DEFAULT_BRAIN_MARKERS["ad_pathology"] +
                          DEFAULT_BRAIN_MARKERS["endothelial"])
        struct_present = [g for g in struct_genes
                            if g in tx["feature_name"].astype(str).unique()]
        for m in methods:
            if m not in method_overflow:
                continue
            out = out_dir / f"{roi.roi_id}_boundary_overflow_{m}.png"
            if _png_skip(out, []):
                continue
            res_m = load_method(m, sample, roi)
            cents = method_cells_with_off[m]

            fig, ax = plt.subplots(figsize=(5.6, 5.4),
                                        facecolor="black", dpi=180)
            _draw_dapi(ax, dapi_full, roi, alpha=0.55)
            _shade_off_tissue(ax)

            # Structure markers (small dots so they don't overpower)
            for g in struct_present[:3]:
                sub = tx[tx["feature_name"].astype(str) == g]
                ax.scatter(sub["x_location"], sub["y_location"],
                            s=2.5, c="#FFD700", alpha=0.65,
                            edgecolors="none")

            # Cell polygons: red fill if any vertex outside tissue OR
            # centroid off-tissue.  Otherwise method-color outline.
            polys = res_m.get("polygons")
            if polys is not None and len(polys):
                bad_p, ok_p = [], []
                off_set = set(cents.index[cents["off_tissue"] == True])
                for idx in range(len(polys)):
                    g = polys.iloc[idx]["geometry"]
                    if g is None or g.is_empty or g.geom_type != "Polygon":
                        continue
                    xs, ys = zip(*list(g.exterior.coords))
                    poly = MplPolygon(np.column_stack([xs, ys]), closed=True)
                    if idx in off_set:
                        bad_p.append(poly)
                    else:
                        ok_p.append(poly)
                if ok_p:
                    ax.add_collection(PatchCollection(
                        ok_p, facecolors="none",
                        edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                        linewidths=0.6, alpha=0.85, zorder=4))
                if bad_p:
                    ax.add_collection(PatchCollection(
                        bad_p, facecolors="#FF000099",
                        edgecolors="#FF0000",
                        linewidths=1.4, alpha=0.95, zorder=5))

            n_off, n_tot = method_overflow[m]
            ovf = 100.0 * n_off / max(1, n_tot)
            handles = [
                Patch(facecolor=OFF_COLOR_RGBA, edgecolor="none",
                      label="off-tissue (DAPI Otsu)"),
                Patch(facecolor="#FF000099", edgecolor="#FF0000",
                      label=f"overflow cell  {n_off}/{n_tot} ({ovf:.1f}%)"),
                plt.Line2D([], [], color=st.METHOD_COLORS.get(m, "#00FFFF"),
                            lw=2.0, label=f"{st.METHOD_DISPLAY.get(m, m)} boundary"),
                plt.Line2D([], [], marker="o", linestyle="", markersize=4,
                            color="#FFD700", label="structure marker"),
            ]
            leg = ax.legend(handles=handles, loc="upper right", fontsize=7,
                              framealpha=0.85, facecolor="black",
                              edgecolor="white")
            for t in leg.get_texts():
                t.set_color("white")
            ax.set_title(
                f"{roi.roi_id} | tissue overflow | "
                f"{st.METHOD_DISPLAY.get(m, m)}\n"
                f"red = cell centroid off-tissue ({ovf:.1f}%)",
                fontsize=9, color="white")
            st.setup_axes(ax, bbox_um=roi.bbox_um)
            _add_dynamic_scale_bar(ax, roi, sample, color="white")
            fig.tight_layout()
            fig.savefig(str(out), dpi=180, bbox_inches="tight",
                          pad_inches=0.05, facecolor="black")
            plt.close(fig)
            out_paths[f"boundary_overflow/{m}"] = out

    # ----- boundary_panel_compare.png (2×3)
    if tissue_mask is not None and method_overflow:
        out_panel = out_dir / f"{roi.roi_id}_boundary_panel_compare.png"
        if not _png_skip(out_panel, []):
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Polygon as MplPolygon
            n = len(methods)
            n_cols = 3 if n > 3 else n
            n_rows = (n + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols,
                                          figsize=(5.0 * n_cols, 5.0 * n_rows + 0.6),
                                          facecolor="black", dpi=160)
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes[None, :]
            elif n_cols == 1:
                axes = axes[:, None]
            for i, m in enumerate(methods):
                r, c = i // n_cols, i % n_cols
                ax = axes[r, c]
                _draw_dapi(ax, dapi_full, roi, alpha=0.55)
                _shade_off_tissue(ax)
                if m not in method_overflow:
                    ax.set_title(f"{st.METHOD_DISPLAY.get(m, m)}\n(unavailable)",
                                    fontsize=10, color="white")
                    st.setup_axes(ax, bbox_um=roi.bbox_um)
                    continue
                res_m = load_method(m, sample, roi)
                cents = method_cells_with_off[m]
                off_set = set(cents.index[cents["off_tissue"] == True])
                bad_p, ok_p = [], []
                polys = res_m.get("polygons")
                if polys is not None and len(polys):
                    for idx in range(len(polys)):
                        g = polys.iloc[idx]["geometry"]
                        if g is None or g.is_empty or g.geom_type != "Polygon":
                            continue
                        xs, ys = zip(*list(g.exterior.coords))
                        poly = MplPolygon(np.column_stack([xs, ys]), closed=True)
                        (bad_p if idx in off_set else ok_p).append(poly)
                if ok_p:
                    ax.add_collection(PatchCollection(
                        ok_p, facecolors="none",
                        edgecolors=st.METHOD_COLORS.get(m, "#00FFFF"),
                        linewidths=0.5, alpha=0.85, zorder=4))
                if bad_p:
                    ax.add_collection(PatchCollection(
                        bad_p, facecolors="#FF000099",
                        edgecolors="#FF0000", linewidths=1.0,
                        alpha=0.95, zorder=5))
                n_off, n_tot = method_overflow[m]
                ovf = 100.0 * n_off / max(1, n_tot)
                ax.set_title(
                    f"{st.METHOD_DISPLAY.get(m, m)}\n"
                    f"overflow: {n_off}/{n_tot} ({ovf:.1f}%)",
                    fontsize=10, color="white")
                st.setup_axes(ax, bbox_um=roi.bbox_um)
            for j in range(n, n_rows * n_cols):
                r, c = j // n_cols, j % n_cols
                axes[r, c].set_visible(False)
            fig.suptitle(
                f"{roi.roi_id} | tissue overflow per-method "
                f"(magenta=off-tissue, red=overflow cell)",
                fontsize=11, color="white")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(str(out_panel), dpi=160, bbox_inches="tight",
                          pad_inches=0.05, facecolor="black")
            plt.close(fig)
            out_paths["boundary_panel_compare"] = out_panel

    # ----- boundary_overflow_fraction_chart.png (P2)
    if method_overflow:
        out = out_dir / f"{roi.roi_id}_boundary_overflow_fraction_chart.png"
        if not _png_skip(out, []):
            keys = [m for m in methods if m in method_overflow]
            ovfs = []
            for m in keys:
                n_off, n_tot = method_overflow[m]
                ovfs.append(100.0 * n_off / max(1, n_tot))
            fig, ax = st.new_figure(figsize=(6.5, 3.5))
            colors = [st.METHOD_COLORS.get(m, "#00FFFF") for m in keys]
            bars = ax.bar(range(len(keys)), ovfs, color=colors,
                            edgecolor="black", linewidth=0.4)
            for b, v in zip(bars, ovfs):
                ax.text(b.get_x() + b.get_width() / 2,
                          v + max(ovfs) * 0.02, f"{v:.1f}%",
                          ha="center", va="bottom", fontsize=8)
            ax.set_xticks(range(len(keys)))
            ax.set_xticklabels([st.METHOD_DISPLAY.get(k, k) for k in keys],
                                  rotation=20, ha="right", fontsize=7)
            ax.set_ylabel("cells with centroid off-tissue (%)")
            ax.set_title(f"{roi.roi_id} | tissue overflow fraction per method")
            st.setup_axes(ax, keep_axis=True, equal=False)
            out_paths["boundary_overflow_fraction_chart"] = _save(fig, out)

    # ----- structure_marker_panel_compare.png (P2 optional)
    if (struct_genes := (DEFAULT_BRAIN_MARKERS["ad_pathology"] +
                            DEFAULT_BRAIN_MARKERS["endothelial"])) and len(tx):
        struct_present = [g for g in struct_genes
                            if g in tx["feature_name"].astype(str).unique()]
        if struct_present:
            out_panel = out_dir / f"{roi.roi_id}_structure_marker_panel_compare.png"
            if not _png_skip(out_panel, []):
                pal = marker_color_palette(struct_present)
                n = len(methods)
                n_cols = 3 if n > 3 else n
                n_rows = (n + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols,
                                              figsize=(5.0 * n_cols, 5.0 * n_rows + 0.6),
                                              facecolor="black", dpi=160)
                if n_rows == 1 and n_cols == 1:
                    axes = np.array([[axes]])
                elif n_rows == 1:
                    axes = axes[None, :]
                elif n_cols == 1:
                    axes = axes[:, None]
                for i, m in enumerate(methods):
                    r, c = i // n_cols, i % n_cols
                    ax = axes[r, c]
                    _draw_dapi(ax, dapi_full, roi, alpha=0.55)
                    if tissue_mask is not None:
                        _shade_off_tissue(ax)
                    for g in struct_present[:4]:
                        sub = tx[tx["feature_name"].astype(str) == g]
                        ax.scatter(sub["x_location"], sub["y_location"],
                                    s=4, c=pal[g], alpha=0.85,
                                    edgecolors="white", linewidths=0.15,
                                    label=f"{g} ({len(sub):,})" if i == 0 else None)
                    try:
                        res_m = load_method(m, sample, roi)
                        if res_m.get("available"):
                            _draw_polygons_outline(
                                ax, res_m["polygons"],
                                color=st.METHOD_COLORS.get(m, "#00FFFF"),
                                lw=0.8, alpha=0.95, halo=False)
                    except Exception:
                        pass
                    ax.set_title(f"{st.METHOD_DISPLAY.get(m, m)}",
                                    fontsize=10, color="white")
                    st.setup_axes(ax, bbox_um=roi.bbox_um)
                for j in range(n, n_rows * n_cols):
                    r, c = j // n_cols, j % n_cols
                    axes[r, c].set_visible(False)
                if axes[0, 0].get_visible():
                    leg = axes[0, 0].legend(loc="upper right", fontsize=7,
                                                  framealpha=0.85,
                                                  facecolor="black",
                                                  edgecolor="white")
                    for t in leg.get_texts():
                        t.set_color("white")
                fig.suptitle(
                    f"{roi.roi_id} | structure marker per-method "
                    f"(magenta=off-tissue)",
                    fontsize=11, color="white")
                fig.tight_layout(rect=(0, 0, 1, 0.95))
                fig.savefig(str(out_panel), dpi=160, bbox_inches="tight",
                              pad_inches=0.05, facecolor="black")
                plt.close(fig)
                out_paths["structure_marker_panel_compare"] = out_panel

    return out_paths


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
CLASS_DISPATCH = {
    "easy_control":  render_easy_control_addons,
    "architecture":  render_architecture_addons,
    "compartment":   render_compartment_addons,
    "low_vsi":       render_low_vsi_addons,
    "high_density":  render_high_density_addons,
    "fold_boundary": render_fold_boundary_addons,
}


def render_class_addons(sample: SampleData, roi: ROIRecord,
                          baseline: str, methods: Sequence[str],
                          outdir: Path) -> Dict[str, Path]:
    fn = CLASS_DISPATCH.get(roi.roi_class)
    if fn is None:
        logger.info(f"[{roi.roi_id}] no class-specific add-ons for class={roi.roi_class}")
        return {}
    try:
        return fn(sample, roi, baseline, methods, outdir)
    except Exception as e:
        logger.exception(f"[{roi.roi_id}] class add-ons failed: {e}")
        return {}
