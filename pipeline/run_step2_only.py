#!/usr/bin/env python
"""Run Step 2 (Segmentation-Free Analysis) only, with file logging."""

import os
import sys
import yaml
import logging
import time

# Ensure pipeline modules are importable
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

import xenium_step2_segmentation_free_analysis as step2

LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


def main():
    # Load config
    config_path = os.path.join(current_dir, "config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    sample_tag = config["sample_tag"]
    input_path = config["input_path"]
    sample_name = os.path.basename(input_path.rstrip(os.sep))

    # Use absolute paths from repo root
    repo_root = os.path.dirname(current_dir)
    base_output = os.path.join(repo_root, config.get("output_dir", "xenium-output"))
    sample_output_dir = os.path.join(base_output, sample_name)

    step2_dir = os.path.join(sample_output_dir, "step2_segmentation_free")
    os.makedirs(step2_dir, exist_ok=True)

    # File logging
    log_dir = os.path.join(step2_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    log_path = os.path.join(log_dir, f"step2_run_{timestamp}.log")

    fh = logging.FileHandler(log_path, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(LOG_FORMAT))
    logging.getLogger().addHandler(fh)

    logger.info(f"Step 2 standalone run | sample={sample_tag}")
    logger.info(f"Log file: {log_path}")

    # Build step2 config
    step1_adata = os.path.join(
        sample_output_dir, "step1_exploration",
        f"{sample_tag}_step1_exploration.h5ad"
    )
    if not os.path.exists(step1_adata):
        logger.error(f"Step 1 output not found: {step1_adata}")
        sys.exit(1)

    step2_config = config.copy()
    step2_config["output_dir"] = step2_dir
    step2_config["previous_step_adata_path"] = step1_adata

    logger.info(f"Input adata: {step1_adata}")
    logger.info(f"Output dir:  {step2_dir}")

    # Run
    try:
        step2.run_step2(step2_config)
        logger.info("Step 2 completed successfully.")
        print(f"\nStep 2 completed. Output: {step2_dir}")
        print(f"Log: {log_path}")
    except Exception as e:
        logger.error(f"Step 2 failed: {e}", exc_info=True)
        print(f"\nStep 2 FAILED. Check log: {log_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
