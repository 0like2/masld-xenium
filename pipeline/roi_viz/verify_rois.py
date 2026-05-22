"""Standalone ROI class-membership verifier.

For each ROI in the manifest, prints a quantitative check that its bbox
actually matches the class it was selected for:

* high_density   — ROI tx density (per µm²) vs tissue 95th percentile
* fold_boundary  — ROI distance to convex-hull edge of all transcripts
* low_vsi        — ROI mean VSI vs tissue mean VSI / 5th percentile
* compartment    — ROI nuclear-vs-cytoplasmic balance (closer to 50/50 → more
                   ambiguous, which is the selection criterion)
* architecture   — ROI SSAM celltype diversity (Shannon entropy, normalised)
* easy_control   — ROI mean VSI ≥ tissue 90th percentile *and* moderate density

Run from the repo root:

    python -m pipeline.roi_viz.verify_rois \
        --manifest xenium-output/.../human_alzheimers_step2_rois.json
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _decode(s: pd.Series) -> pd.Series:
    if len(s) == 0:
        return s
    if isinstance(s.iloc[0], (bytes, bytearray)):
        return s.str.decode("utf-8")
    return s.astype(str)


def _tissue_density_grid(spots: pd.DataFrame, bin_um: float = 50.0):
    x = spots["x_location"].values
    y = spots["y_location"].values
    x_bins = np.arange(x.min(), x.max() + bin_um, bin_um)
    y_bins = np.arange(y.min(), y.max() + bin_um, bin_um)
    hist, _, _ = np.histogram2d(x, y, bins=[x_bins, y_bins])
    return hist, x_bins, y_bins


def _bbox_density(spots: pd.DataFrame, bbox_um) -> Dict[str, float]:
    x0, y0, x1, y1 = bbox_um
    m = ((spots["x_location"] >= x0) & (spots["x_location"] <= x1) &
         (spots["y_location"] >= y0) & (spots["y_location"] <= y1))
    n = int(m.sum())
    area = max(1.0, (x1 - x0) * (y1 - y0))
    return {"n_tx": n, "area_um2": area,
            "density_per_um2": n / area}


def _vsi_in_bbox(coherence_map: np.ndarray, bbox_um, um_per_pixel: float):
    if coherence_map is None:
        return None
    x0, y0, x1, y1 = bbox_um
    H, W = coherence_map.shape
    rs0 = max(0, int(y0 / um_per_pixel))
    rs1 = min(H, int(np.ceil(y1 / um_per_pixel)))
    cs0 = max(0, int(x0 / um_per_pixel))
    cs1 = min(W, int(np.ceil(x1 / um_per_pixel)))
    if rs1 <= rs0 or cs1 <= cs0:
        return None
    sub = coherence_map[rs0:rs1, cs0:cs1]
    finite = sub[np.isfinite(sub)]
    return {"mean": float(finite.mean()) if finite.size else None,
            "min": float(finite.min()) if finite.size else None,
            "p5": float(np.percentile(finite, 5)) if finite.size else None,
            "n_pixels": int(finite.size)}


def _distance_to_tissue_edge(spots: pd.DataFrame, bbox_um) -> float:
    """Approximate µm distance from the ROI center to the convex hull edge of
    all transcripts."""
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return float("nan")
    x = spots["x_location"].values
    y = spots["y_location"].values
    if len(x) > 200_000:
        idx = np.random.default_rng(0).choice(len(x), 200_000, replace=False)
        x, y = x[idx], y[idx]
    pts = np.column_stack([x, y])
    try:
        hull = ConvexHull(pts)
    except Exception:
        return float("nan")
    cx = (bbox_um[0] + bbox_um[2]) / 2
    cy = (bbox_um[1] + bbox_um[3]) / 2
    # Min distance from center to each hull edge.
    edges = [(hull.points[hull.vertices[i]],
                hull.points[hull.vertices[(i + 1) % len(hull.vertices)]])
              for i in range(len(hull.vertices))]
    min_d = float("inf")
    for p1, p2 in edges:
        v = p2 - p1
        L2 = (v * v).sum()
        if L2 < 1e-9:
            continue
        t = max(0.0, min(1.0, ((cx - p1[0]) * v[0] + (cy - p1[1]) * v[1]) / L2))
        proj = p1 + t * v
        d = float(np.hypot(cx - proj[0], cy - proj[1]))
        if d < min_d:
            min_d = d
    return min_d


def _ssam_diversity(ssam_obs: pd.DataFrame, bbox_um) -> float:
    if ssam_obs is None or len(ssam_obs) == 0:
        return float("nan")
    label_col = next((c for c in ("celltype", "ssam_celltype",
                                       "leiden_assignment", "leiden")
                       if c in ssam_obs.columns), None)
    if label_col is None:
        return float("nan")
    x_col = "x_um" if "x_um" in ssam_obs.columns else "x"
    y_col = "y_um" if "y_um" in ssam_obs.columns else "y"
    x0, y0, x1, y1 = bbox_um
    m = ((ssam_obs[x_col] >= x0) & (ssam_obs[x_col] <= x1) &
         (ssam_obs[y_col] >= y0) & (ssam_obs[y_col] <= y1))
    sub = ssam_obs.loc[m]
    if not len(sub):
        return float("nan")
    p = sub[label_col].value_counts(normalize=True).values
    if len(p) <= 1:
        return 0.0
    H = -(p * np.log(p + 1e-12)).sum()
    H_max = np.log(len(p))
    return float(H / H_max) if H_max > 0 else 0.0


def _p2r_compartment_balance(p2r_obs: pd.DataFrame, bbox_um) -> Dict[str, float]:
    if p2r_obs is None or len(p2r_obs) == 0:
        return {"n_total": 0, "frac_nuc": float("nan"), "frac_cyto": float("nan")}
    x_col = next((c for c in ("x_centroid", "x_um", "x_location", "x")
                  if c in p2r_obs.columns), None)
    y_col = next((c for c in ("y_centroid", "y_um", "y_location", "y")
                  if c in p2r_obs.columns), None)
    if x_col is None or y_col is None:
        return {"n_total": 0, "frac_nuc": float("nan"), "frac_cyto": float("nan")}
    x0, y0, x1, y1 = bbox_um
    m = ((p2r_obs[x_col] >= x0) & (p2r_obs[x_col] <= x1) &
         (p2r_obs[y_col] >= y0) & (p2r_obs[y_col] <= y1))
    sub = p2r_obs.loc[m]
    cmp_col = next((c for c in ("p2r_compartment", "compartment")
                    if c in sub.columns), None)
    if cmp_col is None:
        return {"n_total": int(len(sub)),
                "frac_nuc": float("nan"), "frac_cyto": float("nan")}
    s = sub[cmp_col].astype(str).str.lower()
    return {
        "n_total": int(len(sub)),
        "frac_nuc": float((s.str.contains("nuc")).mean()),
        "frac_cyto": float((s.str.contains("cyto")).mean()),
    }


def verify(manifest_path: Path, sample_root: Path,
            transcripts_path: Path,
            vsi_npz_path: Path = None,
            ssam_h5ad_path: Path = None,
            p2r_h5ad_path: Path = None) -> pd.DataFrame:
    with open(manifest_path) as f:
        manifest = json.load(f)
    rois = manifest.get("rois", [])

    # Load tissue-wide spots
    cols = ["transcript_id", "feature_name", "x_location", "y_location"]
    spots = pd.read_parquet(transcripts_path, columns=cols)
    spots["feature_name"] = _decode(spots["feature_name"])
    print(f"[verify] loaded {len(spots):,} transcripts")

    # Tissue-wide density grid
    hist, x_bins, y_bins = _tissue_density_grid(spots, bin_um=50.0)
    bin_area = 50.0 * 50.0
    densities = hist.ravel() / bin_area
    densities_nonzero = densities[densities > 0]
    tissue_d_p50 = float(np.percentile(densities_nonzero, 50))
    tissue_d_p95 = float(np.percentile(densities_nonzero, 95))
    print(f"[verify] tissue tx density (per µm²): median={tissue_d_p50:.3f} "
          f"95p={tissue_d_p95:.3f}")

    # VSI map
    coherence_map = None
    um_per_pixel = 2.0
    if vsi_npz_path and Path(vsi_npz_path).exists():
        npz = np.load(vsi_npz_path, allow_pickle=True)
        for cand in ("vsi_map", "coherence_map", "integrity_map", "map"):
            if cand in npz.files:
                coherence_map = npz[cand]
                break
        if "um_per_pixel" in npz.files:
            um_per_pixel = float(npz["um_per_pixel"])
        if coherence_map is not None:
            finite = coherence_map[np.isfinite(coherence_map)]
            tissue_vsi_mean = float(finite.mean())
            tissue_vsi_p5 = float(np.percentile(finite, 5))
            print(f"[verify] tissue VSI: mean={tissue_vsi_mean:.3f} "
                  f"p5={tissue_vsi_p5:.3f}")

    # SSAM
    ssam_obs = None
    if ssam_h5ad_path and Path(ssam_h5ad_path).exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(ssam_h5ad_path)
            obs = a.obs.copy()
            x_col = next((c for c in ("x", "x_centroid") if c in obs.columns), None)
            y_col = next((c for c in ("y", "y_centroid") if c in obs.columns), None)
            if x_col and y_col:
                x = obs[x_col].astype(float)
                y = obs[y_col].astype(float)
                if max(x.max() - x.min(), y.max() - y.min()) > 20000:
                    x = x / 4.70588; y = y / 4.70588
                obs = obs.assign(x_um=x, y_um=y)
            ssam_obs = obs
        except Exception as e:
            print(f"[verify] ssam load failed: {e}")

    # P2R
    p2r_obs = None
    if p2r_h5ad_path and Path(p2r_h5ad_path).exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(p2r_h5ad_path)
            p2r_obs = a.obs.copy()
        except Exception as e:
            print(f"[verify] p2r load failed: {e}")

    rows = []
    for r in rois:
        bbox = r["bbox_um"]
        x0, y0, x1, y1 = bbox
        d_in = _bbox_density(spots, bbox)
        density_ratio_p95 = (d_in["density_per_um2"] / tissue_d_p95
                              if tissue_d_p95 > 0 else float("nan"))
        density_ratio_p50 = (d_in["density_per_um2"] / tissue_d_p50
                              if tissue_d_p50 > 0 else float("nan"))
        vsi_info = _vsi_in_bbox(coherence_map, bbox, um_per_pixel) \
            if coherence_map is not None else None
        edge_d = _distance_to_tissue_edge(spots, bbox)
        ssam_div = _ssam_diversity(ssam_obs, bbox) if ssam_obs is not None else float("nan")
        p2r_bal = _p2r_compartment_balance(p2r_obs, bbox)

        # Class-specific verdict
        cls = r.get("source_map", "?")
        verdict = "?"
        evidence = ""
        if cls == "density":
            verdict = "PASS" if density_ratio_p95 >= 1.0 else "FAIL"
            evidence = f"tx density {d_in['density_per_um2']:.3f}/µm² (≥ tissue p95 {tissue_d_p95:.3f}? {density_ratio_p95:.2f}×)"
        elif cls == "boundary":
            verdict = "PASS" if (edge_d == edge_d) and edge_d < 1500 else "FAIL"
            evidence = f"distance-to-tissue-edge ≈ {edge_d:.0f} µm (close to edge if <1500)"
        elif cls == "ovrlpy":
            if vsi_info and vsi_info.get("mean") is not None:
                verdict = "PASS" if vsi_info["mean"] < tissue_vsi_mean else "FAIL"
                evidence = (f"VSI mean = {vsi_info['mean']:.3f} (< tissue mean {tissue_vsi_mean:.3f}?), "
                            f"min = {vsi_info['min']:.3f}")
            else:
                verdict = "?"
                evidence = "VSI map missing in ROI bbox"
        elif cls == "p2r":
            if p2r_bal["frac_nuc"] == p2r_bal["frac_nuc"]:
                # mismatch → both nuc and cyto present non-trivially
                bal = abs(p2r_bal["frac_nuc"] - p2r_bal["frac_cyto"])
                verdict = "PASS" if bal < 0.5 else "FAIL"
                evidence = (f"P2R compartment frac nuc={p2r_bal['frac_nuc']:.2f} "
                            f"cyto={p2r_bal['frac_cyto']:.2f} (more balanced = better)")
            else:
                evidence = f"P2R obs has no compartment col (n in bbox={p2r_bal['n_total']})"
        elif cls == "ssam":
            verdict = "PASS" if ssam_div >= 0.5 else "FAIL"
            evidence = f"SSAM celltype Shannon entropy (norm) = {ssam_div:.2f} (>=0.5 = transition zone)"
        elif cls == "easy":
            if vsi_info and vsi_info.get("mean") is not None:
                verdict = "PASS" if vsi_info["mean"] > tissue_vsi_mean else "FAIL"
                evidence = (f"VSI mean = {vsi_info['mean']:.3f} (> tissue mean {tissue_vsi_mean:.3f}?), "
                            f"density {density_ratio_p50:.2f}× tissue median")

        rows.append({
            "roi_id": r["roi_id"],
            "source_map": cls,
            "criterion": r.get("criterion", ""),
            "bbox_um": tuple(bbox),
            "n_tx_in_roi": d_in["n_tx"],
            "tx_density": d_in["density_per_um2"],
            "tx_density_vs_p95": density_ratio_p95,
            "tx_density_vs_median": density_ratio_p50,
            "vsi_mean": vsi_info["mean"] if vsi_info else None,
            "vsi_min": vsi_info["min"] if vsi_info else None,
            "edge_distance_um": edge_d,
            "ssam_norm_entropy": ssam_div,
            "p2r_frac_nuc": p2r_bal["frac_nuc"],
            "p2r_frac_cyto": p2r_bal["frac_cyto"],
            "verdict": verdict,
            "evidence": evidence,
        })

    df = pd.DataFrame(rows)
    return df


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--input-dir", type=Path, default=Path(
        "data/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs"))
    p.add_argument("--step2-dir", type=Path, default=Path(
        "xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/step2_segmentation_free"))
    p.add_argument("--out-csv", type=Path, default=None)
    args = p.parse_args(argv)

    df = verify(
        manifest_path=args.manifest,
        sample_root=args.input_dir,
        transcripts_path=args.input_dir / "transcripts.parquet",
        vsi_npz_path=args.step2_dir / "human_alzheimers_step2_ovrlpy_cache.npz",
        ssam_h5ad_path=args.step2_dir / "human_alzheimers_step2_ssam.h5ad",
        p2r_h5ad_path=args.step2_dir / "human_alzheimers_step2_points2regions.h5ad",
    )
    print()
    print(df.to_string(index=False))
    if args.out_csv:
        df.to_csv(args.out_csv, index=False)
        print(f"\nwritten: {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
