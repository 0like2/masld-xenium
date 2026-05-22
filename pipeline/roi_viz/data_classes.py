"""Dataclasses shared across the roi_viz pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


SOURCE_TO_CLASS = {
    "density":  "high_density",
    "boundary": "fold_boundary",
    "ovrlpy":   "low_vsi",
    "p2r":      "compartment",
    "ssam":     "architecture",
    "easy":     "easy_control",
    "easy_control": "easy_control",
}


@dataclass
class ROIRecord:
    """Canonical ROI record consumed by the pipeline."""
    roi_id: str
    sample_id: str
    roi_class: str
    bbox_um: Tuple[float, float, float, float]   # (x_min, y_min, x_max, y_max)
    source_map: str = "unknown"
    hypothesis: str = "unknown"
    selection_reason: str = ""
    primary_reference: Optional[str] = None
    secondary_reference: Optional[str] = None
    priority: str = "medium"
    notes: str = ""

    @classmethod
    def from_manifest_entry(cls, raw: Dict[str, Any], sample_id: str) -> "ROIRecord":
        bbox = raw.get("bbox_um") or raw.get("bbox")
        if bbox is None or len(bbox) != 4:
            raise ValueError(f"Invalid bbox in ROI manifest entry: {raw}")
        bbox_t = tuple(float(v) for v in bbox)

        roi_id = str(raw.get("roi_id") or raw.get("id") or "roi_unknown")
        source = str(raw.get("source_map") or "unknown")
        crit = str(raw.get("criterion") or "")
        hyp = str(raw.get("hypothesis") or "unknown")

        # Map source_map → roi_class. Fall back to criterion or string match.
        roi_class = SOURCE_TO_CLASS.get(source.lower())
        if roi_class is None:
            for key, cls_name in SOURCE_TO_CLASS.items():
                if key in roi_id.lower() or key in crit.lower() or key in hyp.lower():
                    roi_class = cls_name
                    break
        if roi_class is None:
            roi_class = "unknown"

        primary_ref = {
            "ssam": "ssam", "p2r": "p2r", "ovrlpy": "vsi",
            "boundary": "vsi", "density": "p2r", "easy": "ssam",
        }.get(source.lower())

        return cls(
            roi_id=roi_id,
            sample_id=sample_id,
            roi_class=roi_class,
            bbox_um=bbox_t,
            source_map=source,
            hypothesis=hyp,
            selection_reason=crit or hyp,
            primary_reference=primary_ref,
            priority=str(raw.get("priority") or "medium"),
            notes=str(raw.get("notes") or ""),
        )


@dataclass
class SegmentationMethod:
    """Per-method handle resolved at config time."""
    key: str                    # 'xenium_nucleus', 'cellpose_nuclei', ...
    display_name: str
    boundary_kind: str          # 'polygon' | 'label_map' | 'centroid_radius'
    paths: Dict[str, str] = field(default_factory=dict)
    params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class MarkerConfig:
    """Marker selection per ROI class."""
    target: List[str] = field(default_factory=list)
    competing: List[str] = field(default_factory=list)
    structural: List[str] = field(default_factory=list)
    nuclear: List[str] = field(default_factory=list)
    cytoplasmic: List[str] = field(default_factory=list)


@dataclass
class PlotConfig:
    """Style knobs (kept in style.py for default values)."""
    figsize: Tuple[float, float] = (6.0, 6.0)
    dpi: int = 200
    boundary_lw: float = 1.0
    boundary_alpha: float = 0.9
    transcript_size: float = 0.3
    transcript_alpha: float = 0.6
    marker_size: float = 4.0
    scale_bar_um: float = 50.0
    facecolor: str = "white"
    save_kwargs: Dict[str, Any] = field(default_factory=lambda: {
        "bbox_inches": "tight",
        "pad_inches": 0.05,
        "facecolor": "white",
        "edgecolor": "none",
    })


@dataclass
class SampleData:
    """Lazy bundle of sample-level data sources.

    Concrete loaders only resolve a path each; downstream modules call
    cropping/loaders to fetch the actual ROI-cropped arrays via cache.
    """
    sample_id: str
    sample_tag: str
    input_path: Path
    output_dir: Path
    dapi_path: Optional[Path] = None
    transcripts_path: Optional[Path] = None
    xenium_nucleus_path: Optional[Path] = None
    step2_dir: Optional[Path] = None
    step3_dir: Optional[Path] = None
    step5_dir: Optional[Path] = None
    step6_dir: Optional[Path] = None
    ssam_h5ad_path: Optional[Path] = None
    p2r_csv_path: Optional[Path] = None
    p2r_h5ad_path: Optional[Path] = None
    vsi_map_path: Optional[Path] = None
    step6_combined_h5ad: Optional[Path] = None
    step6_crop_combined_h5ad: Optional[Path] = None
    um_per_pixel_inv: float = 4.70588
    extras: Dict[str, Any] = field(default_factory=dict)
