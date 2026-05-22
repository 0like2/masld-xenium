"""Marker gene selection per ROI class.

The selection is intentionally simple: pick a small set from a curated brain
panel based on ROI class. The user can override per-ROI markers via
``roi.notes`` JSON (`{"markers": ["GFAP","SLC17A7"]}`) or via a config block.
"""

from __future__ import annotations

import json
import logging
from typing import Dict, List, Optional, Sequence

import pandas as pd

from .data_classes import ROIRecord, MarkerConfig

logger = logging.getLogger(__name__)


# Curated brain markers (Alzheimer's panel). Cross-referenced with the gene
# panel in this dataset and Step 4 marker analyses.
DEFAULT_BRAIN_MARKERS = {
    "neuron":          ["SLC17A7", "RBFOX3", "NRGN"],
    "astrocyte":       ["GFAP", "AQP4", "CLU"],
    "oligodendrocyte": ["MOBP", "OLIG2", "PLP1"],
    "microglia":       ["P2RY12", "C3", "PTPRC"],
    "endothelial":     ["FLT1", "CLDN5", "PECAM1"],
    "inhibitory":      ["GAD1", "GAD2"],
    "ad_pathology":    ["APP", "APOE"],
    "nuclear":         ["NEAT1", "MALAT1"],
    "cytoplasmic":     ["MTRNR2L12", "EEF1A1", "FTH1"],
}


CLASS_TO_PANEL = {
    "high_density":  ["neuron", "astrocyte", "oligodendrocyte"],
    "compartment":   ["nuclear", "cytoplasmic", "neuron"],
    "low_vsi":       ["neuron", "oligodendrocyte", "astrocyte"],
    "architecture":  ["neuron", "astrocyte", "inhibitory"],
    "fold_boundary": ["endothelial", "astrocyte", "ad_pathology"],
    "easy_control":  ["neuron", "astrocyte", "oligodendrocyte"],
    "unknown":       ["neuron", "astrocyte", "oligodendrocyte"],
}


def select_markers(roi: ROIRecord,
                    transcripts_in_roi: Optional[pd.DataFrame] = None,
                    n_per_group: int = 1,
                    user_override: Optional[Sequence[str]] = None) -> MarkerConfig:
    """Choose marker genes for a ROI.

    Strategy:
    1. If *user_override* is provided, use it (target=first, competing=second,
       structural=third).
    2. Else if ``roi.notes`` contains JSON {"markers": [...]}, use those.
    3. Else use ``CLASS_TO_PANEL[roi.roi_class]`` to pick markers from
       ``DEFAULT_BRAIN_MARKERS``.
    4. If *transcripts_in_roi* is provided, filter to genes that actually
       appear in the ROI (sorted by abundance).
    """
    chosen = list(user_override or [])

    if not chosen and roi.notes:
        try:
            obj = json.loads(roi.notes)
            if isinstance(obj, dict) and "markers" in obj:
                chosen = list(obj["markers"])
        except Exception:
            pass

    if not chosen:
        groups = CLASS_TO_PANEL.get(roi.roi_class, CLASS_TO_PANEL["unknown"])
        for g in groups:
            chosen.extend(DEFAULT_BRAIN_MARKERS.get(g, [])[:n_per_group])

    if transcripts_in_roi is not None and len(transcripts_in_roi):
        present = set(transcripts_in_roi["feature_name"].astype(str).unique())
        kept = [g for g in chosen if g in present]
        if not kept:
            # Fallback: top 3 most abundant in ROI.
            top = (transcripts_in_roi["feature_name"].astype(str)
                   .value_counts().head(3).index.tolist())
            kept = top
        chosen = kept

    chosen = list(dict.fromkeys(chosen))[:6]   # dedupe + cap

    nuclear = [g for g in DEFAULT_BRAIN_MARKERS["nuclear"] if g in chosen]
    cyto = [g for g in DEFAULT_BRAIN_MARKERS["cytoplasmic"] if g in chosen]

    return MarkerConfig(
        target=chosen[:1],
        competing=chosen[1:2],
        structural=chosen[2:3],
        nuclear=nuclear,
        cytoplasmic=cyto,
    )


def marker_color_palette(markers: Sequence[str]) -> Dict[str, str]:
    """Return a stable color dict for a small set of markers."""
    base = ["#FF3B30", "#34C759", "#3478F6", "#FF9500",
            "#AF52DE", "#FF2D92", "#5AC8FA", "#FFCC00"]
    return {g: base[i % len(base)] for i, g in enumerate(markers)}
