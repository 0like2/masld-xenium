#!/usr/bin/env python
"""Re-run SSAM analysis only (skip P2R/ovrlpy/distance which are unchanged).

Usage:
    conda activate masld
    python pipeline/rerun_ssam_only.py
"""

import os
import sys
import yaml
import logging
import time
import glob

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

import scanpy as sc
import xenium_step2_segmentation_free_analysis as step2

LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


def main():
    config_path = os.path.join(current_dir, "config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    sample_tag = config["sample_tag"]
    input_path = config["input_path"]
    sample_name = os.path.basename(input_path.rstrip(os.sep))

    repo_root = os.path.dirname(current_dir)
    base_output = os.path.join(repo_root, config.get("output_dir", "xenium-output"))
    sample_output_dir = os.path.join(base_output, sample_name)
    step2_dir = os.path.join(sample_output_dir, "step2_segmentation_free")

    # File logging
    log_dir = os.path.join(step2_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    log_path = os.path.join(log_dir, f"ssam_rerun_{timestamp}.log")

    fh = logging.FileHandler(log_path, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(LOG_FORMAT))
    logging.getLogger().addHandler(fh)

    # Done marker
    done_marker = os.path.join(step2_dir, "SSAM_RERUN_DONE.marker")
    if os.path.exists(done_marker):
        os.remove(done_marker)

    logger.info(f"=== SSAM-only re-run | sample={sample_tag} ===")
    logger.info(f"Log: {log_path}")

    # Load existing step2 checkpoint (has P2R/distance results already)
    adata_path = os.path.join(step2_dir, f"{sample_tag}_step2_points2regions.h5ad")
    if not os.path.exists(adata_path):
        logger.error(f"Step 2 checkpoint not found: {adata_path}")
        sys.exit(1)

    logger.info(f"Loading existing Step 2 adata: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    logger.info(f"  adata: {adata.n_obs} cells, {adata.n_vars} genes")

    # Remove old SSAM output files (preserve ssam_cache dir + h5ad cache)
    old_ssam = glob.glob(os.path.join(step2_dir, f"*ssam*"))
    if old_ssam:
        logger.info(f"Removing old SSAM plot files (preserving cache)...")
        for f_path in old_ssam:
            if os.path.isdir(f_path):
                continue  # Skip ssam_cache directory
            if f_path.endswith('_ssam.h5ad'):
                continue  # Preserve h5ad cache for faster reruns
            os.remove(f_path)
            logger.info(f"  Removed: {os.path.basename(f_path)}")

    # Build config for step2
    step2_config = config.copy()
    step2_config["output_dir"] = step2_dir

    # Fix sc_reference path: config uses relative "data/scRNAseq", resolve to absolute
    sc_ref = step2_config.get('sc_reference', {})
    ref_dest = sc_ref.get('dest_dir', '')
    if ref_dest and not os.path.isabs(ref_dest):
        # Try both repo_root-relative and output-relative paths
        for candidate in [
            os.path.join(base_output, ref_dest),  # xenium-output/data/scRNAseq
            os.path.join(repo_root, ref_dest),     # repo/data/scRNAseq
        ]:
            if os.path.isdir(candidate):
                step2_config['sc_reference']['dest_dir'] = candidate
                logger.info(f"Resolved sc_reference path: {candidate}")
                break

    # Run SSAM only
    t0 = time.time()
    try:
        adata = step2.run_ssam_analysis(adata, step2_config, step2_dir, sample_tag)
        elapsed = time.time() - t0
        logger.info(f"SSAM analysis completed in {elapsed/60:.1f} min")

        # Save updated adata
        adata.write_h5ad(adata_path)
        logger.info(f"Updated adata saved to {adata_path}")

        # Done marker
        with open(done_marker, 'w') as f:
            f.write(f"SSAM re-run completed at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Elapsed: {elapsed/60:.1f} min\n")
            f.write(f"Log: {log_path}\n")

        logger.info(f"DONE marker: {done_marker}")
        print(f"\n✓ SSAM re-run completed ({elapsed/60:.1f} min)")
        print(f"  Output: {step2_dir}")
        print(f"  Log: {log_path}")
        print(f"  Marker: {done_marker}")

    except Exception as e:
        elapsed = time.time() - t0
        logger.error(f"SSAM re-run FAILED after {elapsed/60:.1f} min: {e}", exc_info=True)
        # Write failure marker
        with open(done_marker, 'w') as f:
            f.write(f"SSAM re-run FAILED at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Error: {e}\n")
            f.write(f"Log: {log_path}\n")
        print(f"\n✗ SSAM re-run FAILED ({elapsed/60:.1f} min). Check log: {log_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
