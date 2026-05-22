"""End-to-end ROI visualization orchestrator.

Loads the Step-2 ROI manifest, builds a :class:`SampleData` from project
config, then for each ROI:

1. Renders C1–C4 + M1–M4 common images for every requested method.
2. Renders the highlight set (boundary overlap, difference mask,
   transcript reassignment) for *each* method against the chosen baseline.
3. Renders ROI-class-specific add-ons (where data is available).
4. Computes per-(ROI, method) metrics and aggregates them into a single
   ``roi_metrics.csv``.

Cache invalidation is plumbed through ``cache.configure``: pass
``--force``, ``--force-cache=stage``, ``--force-roi=ID``.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from . import __version__, cache
from . import style as st
from .class_renders import render_class_addons
from .common_renders import render_all_common
from .data_classes import ROIRecord, SampleData
from .highlight_renders import render_all_highlights
from .metrics import aggregate_metrics

logger = logging.getLogger(__name__)


DEFAULT_METHODS = (
    "xenium_nucleus",
    "cellpose_nuclei",
    "rigid_expansion",
    "optimal_expansion",
    "baysor",
)


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------
def load_manifest(path: Path) -> List[Dict[str, Any]]:
    with open(path) as f:
        d = json.load(f)
    rois = d.get("rois") or d.get("roi_list") or []
    return list(rois)


def build_roi_records(manifest_rois: Iterable[Dict[str, Any]],
                        sample_id: str) -> List[ROIRecord]:
    out = []
    for raw in manifest_rois:
        try:
            out.append(ROIRecord.from_manifest_entry(raw, sample_id))
        except Exception as e:
            logger.warning(f"skipping ROI entry {raw.get('roi_id')}: {e}")
    return out


# ---------------------------------------------------------------------------
# SampleData construction from a project config
# ---------------------------------------------------------------------------
def build_sample_data(config: Dict[str, Any]) -> SampleData:
    """Build a :class:`SampleData` from a parsed config dict.

    Expected keys:
    - ``input_path``      Xenium output bundle
    - ``output_dir``      project output root
    - ``sample_tag``      naming prefix (e.g. 'human_alzheimers')
    """
    in_path = Path(config["input_path"])
    out_dir_root = Path(config["output_dir"])
    sample_tag = config["sample_tag"]
    sample_id = config.get("sample_id", sample_tag)

    # The pipeline writes outputs under
    #   {output_dir}/{input_path.name}/...
    out_per_sample = out_dir_root / in_path.name
    if not out_per_sample.exists():
        # Fallback to first existing per-sample dir.
        candidates = [d for d in out_dir_root.iterdir() if d.is_dir()] \
            if out_dir_root.exists() else []
        if candidates:
            out_per_sample = candidates[0]

    return SampleData(
        sample_id=sample_id,
        sample_tag=sample_tag,
        input_path=in_path,
        output_dir=out_per_sample,
        dapi_path=in_path / "morphology_focus.ome.tif",
        transcripts_path=in_path / "transcripts.parquet",
        xenium_nucleus_path=in_path / "nucleus_boundaries.parquet",
        step2_dir=out_per_sample / "step2_segmentation_free",
        step3_dir=out_per_sample / "step3_resegmentation",
        step5_dir=out_per_sample / "step5_optimal_expansion",
        step6_dir=out_per_sample / "step6_benchmark",
        ssam_h5ad_path=out_per_sample / "step2_segmentation_free"
                        / f"{sample_tag}_step2_ssam.h5ad",
        p2r_csv_path=out_per_sample / "step2_segmentation_free"
                      / f"{sample_tag}_step2_p2r_classification.csv",
        p2r_h5ad_path=out_per_sample / "step2_segmentation_free"
                       / f"{sample_tag}_step2_points2regions.h5ad",
        vsi_map_path=out_per_sample / "step2_segmentation_free"
                      / f"{sample_tag}_step2_ovrlpy_cache.npz",
        step6_combined_h5ad=out_per_sample / "step6_benchmark"
                              / "benchmark_combined.h5ad",
        step6_crop_combined_h5ad=out_per_sample / "step6_benchmark"
                                   / "crop_comparison" / "benchmark_combined.h5ad",
    )


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------
def run_roi_pipeline(config: Dict[str, Any],
                       manifest_path: Path,
                       methods: Sequence[str] = DEFAULT_METHODS,
                       baseline: str = "xenium_nucleus",
                       outdir: Optional[Path] = None,
                       roi_ids: Optional[Sequence[str]] = None,
                       force_stages: Optional[Sequence[str]] = None,
                       force_rois: Optional[Sequence[str]] = None,
                       skip_metrics: bool = False) -> Dict[str, Any]:
    sample = build_sample_data(config)

    if outdir is None:
        outdir = sample.output_dir / "roi_test_output" / "roi_visualization"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cache.configure(root=outdir,
                     force_stages=force_stages,
                     force_rois=force_rois,
                     version=__version__)

    st.apply_rc_params()

    # Logging to file
    log_path = outdir / "roi_visualization.log"
    fh = logging.FileHandler(log_path, mode="a")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
    root_logger = logging.getLogger()
    root_logger.addHandler(fh)

    logger.info(f"roi_viz run start | version={__version__} | sample={sample.sample_tag}")
    logger.info(f"manifest={manifest_path} | baseline={baseline} | methods={list(methods)}")
    logger.info(f"outdir={outdir} | force_stages={force_stages} | force_rois={force_rois}")

    rois_raw = load_manifest(manifest_path)
    rois = build_roi_records(rois_raw, sample.sample_id)
    if roi_ids:
        rois = [r for r in rois if r.roi_id in set(roi_ids)]
    logger.info(f"running {len(rois)} ROI(s)")

    summary: Dict[str, Any] = {"rois": {}, "outdir": str(outdir),
                                 "log": str(log_path)}

    for roi in rois:
        roi_summary: Dict[str, Any] = {"roi_class": roi.roi_class,
                                         "bbox_um": list(roi.bbox_um)}
        try:
            common_paths = render_all_common(sample, roi, outdir, methods)
            highlight_paths = render_all_highlights(
                sample, roi, baseline=baseline, methods=methods, outdir=outdir)
            class_paths = render_class_addons(
                sample, roi, baseline=baseline, methods=methods, outdir=outdir)
            def _flatten(d):
                """Render a {key: Path | dict[k,Path]} mapping into all-strings."""
                out = {}
                for k, v in d.items():
                    if isinstance(v, dict):
                        out[k] = {kk: str(vv) for kk, vv in v.items()}
                    else:
                        out[k] = str(v)
                return out
            roi_summary["common"] = _flatten(common_paths)
            roi_summary["highlights"] = _flatten(highlight_paths)
            roi_summary["class_addons"] = _flatten(class_paths)
        except Exception as e:
            logger.exception(f"render failed for {roi.roi_id}: {e}")
            roi_summary["error"] = str(e)
        summary["rois"][roi.roi_id] = roi_summary

    # Aggregate metrics
    if not skip_metrics and rois:
        try:
            csv_path = aggregate_metrics(sample, rois, methods, outdir)
            summary["roi_metrics_csv"] = str(csv_path)
        except Exception as e:
            logger.exception(f"metric aggregation failed: {e}")
            summary["metric_error"] = str(e)

    logger.info("roi_viz run complete")
    root_logger.removeHandler(fh)
    fh.close()
    return summary
