#!/usr/bin/env python
"""CLI entry point for the ROI visualization pipeline.

Example:
    python pipeline/run_roi_visualization.py \
        --config pipeline/config.yaml \
        --roi-manifest xenium-output/.../roi_test_output/human_alzheimers_step2_rois.json \
        --methods xenium_nucleus,cellpose_nuclei,rigid_expansion,optimal_expansion,baysor \
        --baseline xenium_nucleus

Cache controls:
    --force                   ignore all caches (full rebuild)
    --force-cache=raw         invalidate a single stage
                              (raw|boundaries|diff|metrics|render)
    --force-roi=ID[,ID...]    only invalidate caches for given ROI(s)
    --force-render            re-render PNGs from cached artifacts (style tweaks)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

# Allow `python pipeline/run_roi_visualization.py` as well as
# `python -m pipeline.run_roi_visualization`.
_HERE = Path(__file__).resolve().parent
if str(_HERE.parent) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))

from pipeline.roi_viz.pipeline import (DEFAULT_METHODS, run_roi_pipeline)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, type=Path,
                    help="Path to pipeline config.yaml")
    p.add_argument("--roi-manifest", required=True, type=Path,
                    help="Path to step2_rois.json")
    p.add_argument("--methods", default=",".join(DEFAULT_METHODS),
                    help="Comma-separated method keys (default: all 5)")
    p.add_argument("--baseline", default="xenium_nucleus",
                    help="Baseline method for highlight comparisons")
    p.add_argument("--outdir", type=Path, default=None,
                    help="Output directory (default: <sample>/roi_test_output/roi_visualization)")
    p.add_argument("--roi-ids", default=None,
                    help="Optional comma-separated subset of roi_id values")

    p.add_argument("--force", action="store_true",
                    help="Invalidate all caches (full rebuild)")
    p.add_argument("--force-cache", default="",
                    help="Comma-separated stages to invalidate "
                         "(raw|boundaries|diff|metrics|render)")
    p.add_argument("--force-roi", default="",
                    help="Comma-separated ROI ids to invalidate")
    p.add_argument("--force-render", action="store_true",
                    help="Shortcut for --force-cache=render")

    p.add_argument("--skip-metrics", action="store_true",
                    help="Skip metric aggregation (faster smoke runs)")

    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

    with open(args.config) as f:
        config = yaml.safe_load(f)

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    roi_ids = [r.strip() for r in args.roi_ids.split(",")] if args.roi_ids else None
    force_stages = []
    if args.force:
        force_stages.append("all")
    if args.force_render:
        force_stages.append("render")
    if args.force_cache:
        force_stages.extend([s.strip() for s in args.force_cache.split(",") if s.strip()])
    force_rois = [r.strip() for r in args.force_roi.split(",")] if args.force_roi else None

    summary = run_roi_pipeline(
        config=config,
        manifest_path=args.roi_manifest,
        methods=methods,
        baseline=args.baseline,
        outdir=args.outdir,
        roi_ids=roi_ids,
        force_stages=force_stages or None,
        force_rois=force_rois,
        skip_metrics=args.skip_metrics,
    )

    print()
    print(f"ROIs processed: {len(summary['rois'])}")
    print(f"Output: {summary['outdir']}")
    print(f"Log:    {summary['log']}")
    if "roi_metrics_csv" in summary:
        print(f"Metrics: {summary['roi_metrics_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
