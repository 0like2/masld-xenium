"""Incrementally run Baysor on extra ROIs (e.g. ssam_transition_*).

Existing ``step6_benchmark/baysor_run/roi_targeted/`` is preserved. The
script:

1. Reads the augmented manifest (with new SSAM ROIs).
2. Looks at existing ``roi_targeted_metadata.json`` to figure out which
   ROI ids have *already* been processed.
3. For every ROI in the manifest that is **not** yet processed, runs
   Baysor in its own ``region_<roi_id>/`` directory using the same
   helpers as ``xenium_step6_segmentation_benchmark._run_baysor_roi_targeted``.
4. After all new regions finish, re-merges ``segmentation.csv`` and
   ``segmentation_polygons_2d.json`` from ALL regions (existing + new)
   so the consolidated outputs include every ROI.
5. Writes the updated metadata.

Usage::

    cd pipeline
    python -m roi_viz.baysor_incremental \\
        --config config.yaml \\
        --manifest ../xenium-output/.../roi_test_output/human_alzheimers_step2_rois_with_ssam.json \\
        --only-roi-ids ssam_transition_1,ssam_transition_2,ssam_transition_3,ssam_transition_4
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import yaml

# Ensure pipeline/ is importable so the step6 module can resolve `utils.*`.
_PIPELINE_DIR = Path(__file__).resolve().parent.parent
if str(_PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(_PIPELINE_DIR))

# Heavy imports gated behind sys.path tweak above.
from xenium_step6_segmentation_benchmark import (    # type: ignore
    _compute_roi_baysor_regions,
    _crop_nuclei_mask_for_region,
    _load_spots_df,
    _merge_roi_polygons,
    _merge_roi_results,
    _run_baysor_tile_worker,
)

logger = logging.getLogger(__name__)


def _load_yaml(path: str) -> Dict[str, Any]:
    with open(path) as f:
        return yaml.safe_load(f)


def _resolve_prior_tif(config: Dict[str, Any], step6_out: Path,
                         per_sample_root: Path) -> Optional[Path]:
    sample_tag = config.get("sample_tag", "human_alzheimers")
    bcfg = config.get("baysor", {}) or {}
    prior_full_tif = bcfg.get("prior_segmentation_tif")
    if prior_full_tif and Path(prior_full_tif).exists():
        return Path(prior_full_tif)
    cand = per_sample_root / "step3_resegmentation" / f"{sample_tag}_step3_nuclei_masks.tif"
    if cand.exists():
        return cand
    return None


def _existing_region_meta(metadata_path: Path) -> List[Dict[str, Any]]:
    if not metadata_path.exists():
        return []
    with open(metadata_path) as f:
        meta = json.load(f)
    return list(meta.get("regions", []))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--manifest", required=True,
                     help="ROI manifest (the *augmented* one with new ROIs)")
    ap.add_argument("--only-roi-ids", default="",
                     help="Comma-separated subset; leave empty to process every "
                          "ROI in the manifest that's not already in metadata")
    ap.add_argument("--padding-um", type=float, default=None,
                     help="Override baysor.roi_targeted.padding_um (default 250)")
    ap.add_argument("--n-workers", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                          format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")

    cfg = _load_yaml(args.config)
    sample_tag = cfg.get("sample_tag", "human_alzheimers")
    output_dir = cfg.get("output_dir", "xenium-output")
    in_path = Path(cfg.get("input_path", ""))

    per_sample_root = Path(output_dir) / in_path.name
    step6_out = per_sample_root / "step6_benchmark"
    baysor_out_dir = step6_out / "baysor_run"
    roi_out_dir = baysor_out_dir / "roi_targeted"
    roi_out_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = roi_out_dir / "roi_targeted_metadata.json"

    # 1. Load augmented manifest
    with open(args.manifest) as f:
        manifest = json.load(f)
    rois = manifest.get("rois", [])

    # 2. Decide which ROI ids to process
    existing_regions = _existing_region_meta(metadata_path)
    existing_roi_ids = {rid
                         for r in existing_regions
                         for rid in r.get("roi_ids", [])}
    if args.only_roi_ids.strip():
        wanted = {x.strip() for x in args.only_roi_ids.split(",") if x.strip()}
        rois_to_run = [r for r in rois if r["roi_id"] in wanted
                        and r["roi_id"] not in existing_roi_ids]
        if not rois_to_run:
            logger.warning(f"All requested ROI ids already processed or not in "
                           f"manifest. Wanted={sorted(wanted)} "
                           f"existing={sorted(existing_roi_ids)}")
    else:
        rois_to_run = [r for r in rois if r["roi_id"] not in existing_roi_ids]

    logger.info(f"Existing ROI ids in metadata: {sorted(existing_roi_ids)}")
    logger.info(f"New ROIs to run baysor on:   {[r['roi_id'] for r in rois_to_run]}")

    if not rois_to_run:
        logger.info("Nothing to do.")
        return

    # 3. Baysor parameters from config
    bcfg = cfg.get("baysor", {}) or {}
    rt_conf = bcfg.get("roi_targeted", {}) or {}
    padding_um = args.padding_um if args.padding_um is not None \
        else rt_conf.get("padding_um", 250)
    params = bcfg.get("params", {}) or {}
    executable = bcfg.get("executable_path", "baysor")
    dry_run = args.dry_run or bcfg.get("dry_run", False)

    # 4. Spots & tissue bounds
    spots_df = _load_spots_df(in_path)
    if spots_df is None:
        logger.error("Failed to load transcripts — cannot run baysor.")
        sys.exit(2)
    tissue_bounds = (spots_df["x"].min(), spots_df["x"].max(),
                     spots_df["y"].min(), spots_df["y"].max())

    # 5. Compute regions for the *new* ROIs only (no merge with existing —
    # existing regions are already on disk; they get appended at the merge
    # step, not re-processed).
    new_regions = _compute_roi_baysor_regions(
        rois_to_run, padding_um, tissue_bounds,
        merge_overlapping=rt_conf.get("merge_overlapping", True))
    if not new_regions:
        logger.error("No valid baysor regions for the new ROIs.")
        sys.exit(2)

    # 6. Resolve Cellpose prior TIF
    prior_full_tif = _resolve_prior_tif(cfg, step6_out, per_sample_root)
    um_per_pixel_inv = (cfg.get("resegmentation", {}) or {}).get(
        "um_per_pixel_inv", 4.70588)
    if prior_full_tif and prior_full_tif.exists():
        logger.info(f"Cellpose prior TIF: {prior_full_tif}")
    else:
        logger.warning("Cellpose prior TIF not found — Baysor runs WITHOUT prior.")
        prior_full_tif = None

    # 7. Prepare per-region input + run
    worker_args = []
    new_region_meta: List[Dict[str, Any]] = []

    for region in new_regions:
        rid = region["region_id"]
        px0, px1, py0, py1 = region["padded_bbox"]

        region_spots = spots_df.loc[
            spots_df["x"].between(px0, px1) & spots_df["y"].between(py0, py1),
            ["gene", "x", "y"],
        ]
        if len(region_spots) == 0:
            logger.warning(f"Region {rid}: 0 molecules — skipping.")
            continue

        region_dir = roi_out_dir / f"region_{rid}"
        region_dir.mkdir(parents=True, exist_ok=True)
        spots_file = region_dir / "spots.csv"
        region_spots.to_csv(spots_file, index=False)

        region_prior_tif = None
        if prior_full_tif is not None:
            region_prior_tif = region_dir / "prior_nuclei_mask.tif"
            try:
                _crop_nuclei_mask_for_region(
                    str(prior_full_tif), region["padded_bbox"],
                    um_per_pixel_inv, str(region_prior_tif))
                logger.info(f"Region {rid}: prior cropped "
                              f"({region_prior_tif.stat().st_size/1024:.0f} KB)")
            except Exception as e:
                logger.warning(f"Region {rid}: prior crop failed ({e}); "
                                  "running without prior.")
                region_prior_tif = None

        logger.info(f"Region {rid}: {len(region_spots):,} molecules in "
                       f"[{px0:.0f},{px1:.0f}]x[{py0:.0f},{py1:.0f}] µm")

        worker_args.append((executable, params, str(spots_file), str(region_dir),
                              f"roi_{rid}", dry_run,
                              str(region_prior_tif) if region_prior_tif else None))
        new_region_meta.append({
            "region_id": rid,
            "core_regions": region["core_regions"],
            "roi_ids": region["roi_ids"],
            "region_dir": str(region_dir),
            "n_molecules": int(len(region_spots)),
            "padded_bbox": list(region["padded_bbox"]),
            "prior_tif": str(region_prior_tif) if region_prior_tif else None,
        })

    if not worker_args:
        logger.error("No non-empty new regions.")
        sys.exit(2)

    # 8. Run baysor (sequential by default — workers can be spawned via -n)
    if args.n_workers <= 1 or dry_run:
        results = [_run_baysor_tile_worker(a) for a in worker_args]
    else:
        import multiprocessing
        with multiprocessing.Pool(processes=args.n_workers) as pool:
            results = pool.map(_run_baysor_tile_worker, worker_args)

    new_region_results = [{
        "region_id": meta["region_id"],
        "core_regions": meta["core_regions"],
        "csv_path": csv_path,
    } for meta, csv_path in zip(new_region_meta, results)]

    # 9. Build merge inputs from EXISTING + NEW regions, then re-merge.
    existing_region_results = []
    for r in existing_regions:
        rid = r["region_id"]
        region_csv = Path(r.get("region_dir", str(roi_out_dir / f"region_{rid}"))) \
            / "segmentation.csv"
        if not region_csv.exists():
            # legacy region dirs put segmentation_<label>.csv files; try fallback
            legacy = list(Path(r.get("region_dir",
                                       str(roi_out_dir / f"region_{rid}")))
                              .glob("segmentation*.csv"))
            if legacy:
                region_csv = legacy[0]
        existing_region_results.append({
            "region_id": rid,
            "core_regions": [tuple(c) for c in r.get("core_regions", [])],
            "csv_path": str(region_csv) if region_csv.exists() else None,
        })

    all_region_results = existing_region_results + new_region_results

    merged_csv = _merge_roi_results(all_region_results, str(roi_out_dir))
    _merge_roi_polygons(all_region_results, str(roi_out_dir))

    # 10. Update metadata
    all_region_meta = list(existing_regions) + new_region_meta
    metadata = {
        "mode": "roi_targeted",
        "padding_um": padding_um,
        "n_workers": args.n_workers,
        "merge_overlapping": rt_conf.get("merge_overlapping", True),
        "n_rois": len(rois),
        "n_regions": len(all_region_meta),
        "regions": all_region_meta,
        "n_molecules_total": int(len(spots_df)),
    }
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    logger.info(f"Updated metadata: {metadata_path} ({len(all_region_meta)} regions)")
    logger.info(f"Merged segmentation.csv: {merged_csv}")


if __name__ == "__main__":
    main()
