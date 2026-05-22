"""Augment an existing Step-2 ROI manifest with SSAM transition-zone ROIs.

This calls ``_select_ssam_rois`` (which now has a per-spot Shannon-entropy
fallback) to derive architecture ROIs from the SSAM h5ad without re-running
Step 2. It writes a new manifest JSON with the additional ROIs appended.

Usage::

    python -m pipeline.roi_viz.augment_manifest_with_ssam \
        --config pipeline/config.yaml \
        --manifest xenium-output/.../roi_test_output/human_alzheimers_step2_rois.json \
        --out      xenium-output/.../roi_test_output/human_alzheimers_step2_rois_with_ssam.json \
        --n-rois 4
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(path: str) -> Dict[str, Any]:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    # Resolve relative paths against the repo root (cwd) to match the rest of
    # the pipeline.
    cfg.setdefault("output_dir", "xenium-output")
    return cfg


def _pick_separated_points(candidates: np.ndarray, n: int,
                              min_sep: float) -> List:
    """Greedy: pick top ``n`` points keeping pair-wise separation ≥ ``min_sep``."""
    selected: List = []
    for cand in candidates:
        ok = True
        for s in selected:
            if (cand[0] - s[0]) ** 2 + (cand[1] - s[1]) ** 2 < min_sep ** 2:
                ok = False; break
        if ok:
            selected.append(tuple(map(float, cand)))
        if len(selected) >= n:
            break
    return selected


def _select_ssam_transition_rois(ssam_h5ad_path: Path,
                                    transcripts_path: Path,
                                    n_rois: int, window_um: float,
                                    min_sep_um: float,
                                    ssam_um_per_px: float = 2.0) -> List[Dict[str, Any]]:
    """Per-spot Shannon-entropy based SSAM transition-zone selection.

    Self-contained reimplementation that mirrors the fallback path inside
    ``_select_ssam_rois`` — does not require importing the step2 module.
    """
    import anndata as ad
    from collections import Counter, defaultdict
    from scipy.ndimage import uniform_filter
    from scipy.stats import entropy as _entropy

    a = ad.read_h5ad(str(ssam_h5ad_path))
    obs = a.obs
    label_col = next((c for c in ("celltype", "ssam_celltype",
                                       "leiden_assignment", "leiden")
                       if c in obs.columns), None)
    if label_col is None:
        logger.warning("SSAM h5ad missing celltype/leiden column.")
        return []
    if "x" not in obs.columns or "y" not in obs.columns:
        logger.warning("SSAM h5ad missing x/y columns.")
        return []
    sx = obs["x"].astype(float).values
    sy = obs["y"].astype(float).values
    # SSAM stores obs x/y as VF-grid pixels (1 px = 2 µm) shifted to (0,0).
    # Real µm = obs * um_per_px + spots_x_min. Pull spots min from
    # transcripts.parquet to recover the offset that step2 applied.
    x_off = y_off = 0.0
    try:
        tx = pd.read_parquet(transcripts_path,
                              columns=["x_location", "y_location"])
        x_off = float(tx["x_location"].min())
        y_off = float(tx["y_location"].min())
    except Exception as e:
        logger.warning(f"transcripts read failed ({e}) — assuming offset 0")
    sx = sx * ssam_um_per_px + x_off
    sy = sy * ssam_um_per_px + y_off
    slabel = obs[label_col].astype(str).values

    bin_um = 50.0
    x_bins = np.arange(sx.min(), sx.max() + bin_um, bin_um)
    y_bins = np.arange(sy.min(), sy.max() + bin_um, bin_um)
    xi = np.clip(np.digitize(sx, x_bins) - 1, 0, len(x_bins) - 2)
    yi = np.clip(np.digitize(sy, y_bins) - 1, 0, len(y_bins) - 2)
    ent_grid = np.zeros((len(x_bins) - 1, len(y_bins) - 1))
    cnt_grid = np.zeros_like(ent_grid)
    bag = defaultdict(Counter)
    for ix, iy, lb in zip(xi, yi, slabel):
        bag[(int(ix), int(iy))][lb] += 1
    for (ix, iy), ctr in bag.items():
        vals = np.array(list(ctr.values()), dtype=float)
        if vals.sum() < 5:
            continue
        p = vals / vals.sum()
        ent_grid[ix, iy] = float(_entropy(p, base=2))
        cnt_grid[ix, iy] = float(vals.sum())
    ent_smooth = uniform_filter(ent_grid, size=3)
    ent_smooth[cnt_grid < 5] = 0.0

    n_candidates = n_rois * 20
    flat = ent_smooth.ravel()
    top_idxs = np.argsort(flat)[-n_candidates:][::-1]
    ix_top, iy_top = np.unravel_index(top_idxs, ent_smooth.shape)
    cx_um = x_bins[ix_top] + bin_um / 2
    cy_um = y_bins[iy_top] + bin_um / 2
    cands = np.column_stack([cx_um, cy_um])
    selected = _pick_separated_points(cands, n=n_rois, min_sep=min_sep_um)
    rois: List[Dict[str, Any]] = []
    for cx, cy in selected:
        bi = min(int((cx - x_bins[0]) / bin_um), ent_smooth.shape[0] - 1)
        bj = min(int((cy - y_bins[0]) / bin_um), ent_smooth.shape[1] - 1)
        rois.append({
            "center_um": [float(cx), float(cy)],
            "score": float(ent_smooth[bi, bj]),
            "metadata": {"shannon_entropy": float(ent_smooth[bi, bj]),
                          "n_spots": int(cnt_grid[bi, bj])},
        })
    return rois


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to pipeline/config.yaml")
    ap.add_argument("--manifest", required=True, help="Existing manifest JSON")
    ap.add_argument("--out", required=True, help="Output augmented manifest JSON")
    ap.add_argument("--n-rois", type=int, default=4,
                     help="How many architecture ROIs to add")
    ap.add_argument("--window-um", type=float, default=250.0)
    ap.add_argument("--min-sep-um", type=float, default=200.0)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                          format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")

    cfg = _load_config(args.config)
    sample_tag = cfg.get("sample_tag", "sample")
    output_dir = cfg.get("output_dir", "xenium-output")

    in_path = Path(cfg.get("input_path", ""))
    per_sample = Path(output_dir) / in_path.name
    ssam_h5ad_path = (per_sample / "step2_segmentation_free"
                       / f"{sample_tag}_step2_ssam.h5ad")
    transcripts_path = in_path / "transcripts.parquet"
    if not ssam_h5ad_path.exists():
        logger.error(f"SSAM h5ad not found: {ssam_h5ad_path}")
        raise SystemExit(2)

    ssam_rois = _select_ssam_transition_rois(
        ssam_h5ad_path,
        transcripts_path=transcripts_path,
        n_rois=args.n_rois,
        window_um=args.window_um,
        min_sep_um=args.min_sep_um)
    if not ssam_rois:
        logger.warning("No SSAM ROIs derived — output will equal input manifest.")

    # Load existing manifest
    with open(args.manifest) as f:
        manifest = json.load(f)
    existing = list(manifest.get("rois", []))

    # Convert ROI dicts (center_um → bbox_um) and assign IDs.
    half = args.window_um
    next_id = 1
    for r in ssam_rois:
        cx, cy = r["center_um"]
        bbox = [cx - half, cy - half, cx + half, cy + half]
        existing.append({
            "roi_id": f"ssam_transition_{next_id}",
            "source_map": "ssam",
            "criterion": "transition_zone",
            "hypothesis": "architecture_preservation",
            "bbox_um": bbox,
            "center_um": [float(cx), float(cy)],
            "score": float(r.get("score", 0.0)),
            "priority": "high",
            "metadata": r.get("metadata", {}),
            "notes": "auto-derived (Shannon entropy of SSAM clusters)",
        })
        next_id += 1

    manifest["rois"] = existing
    manifest.setdefault("provenance", {})
    manifest["provenance"]["augmented_with_ssam"] = {
        "n_added": len(ssam_rois),
        "window_um": args.window_um,
        "method": "per-spot Shannon entropy fallback",
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info(f"Augmented manifest written: {out_path}  "
                f"({len(existing)} total ROIs, +{len(ssam_rois)} SSAM)")


if __name__ == "__main__":
    main()
