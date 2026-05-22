"""Automatic scRNA-seq reference downloader for the Xenium pipeline.

Downloads a matching scRNA-seq reference .h5ad when the user has not
provided one manually.  The download is skipped when:
  - ``sc_reference.auto_download`` is ``false`` in the config, or
  - the destination file already exists on disk, or
  - the user has already set a manual path in ``comparison.sc_reference_path``.
"""

import logging
import os
import urllib.request

logger = logging.getLogger(__name__)

# Default URL for the SEA-AD Allen Institute snRNA-seq reference
_DEFAULT_URL = (
    "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/"
    "MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad"
)


def ensure_reference(config: dict, base_output_dir: str) -> str | None:
    """Download scRNA-seq reference if needed and return its local path.

    Parameters
    ----------
    config : dict
        Full pipeline config (as loaded from ``config.yaml``).
    base_output_dir : str
        Root output directory (e.g. ``xenium-output``).  The reference file
        is stored under ``<base_output_dir>/<dest_dir>/<filename>``.

    Returns
    -------
    str or None
        Absolute path to the reference ``.h5ad`` on success, ``None`` when
        the download is disabled or fails.
    """
    ref_cfg = config.get("sc_reference", {})
    if not ref_cfg:
        return None

    if not ref_cfg.get("auto_download", False):
        logger.info("sc_reference.auto_download is disabled – skipping.")
        return None

    # If the user already provided a manual path that exists, respect it
    manual = config.get("comparison", {}).get("sc_reference_path")
    if manual and os.path.isfile(manual):
        logger.info("Using user-provided sc_reference_path: %s", manual)
        return manual

    url = ref_cfg.get("url", _DEFAULT_URL)
    dest_dir = ref_cfg.get("dest_dir", "data/scRNAseq")
    filename = os.path.basename(url)

    dest_folder = os.path.join(base_output_dir, dest_dir)
    os.makedirs(dest_folder, exist_ok=True)
    dest_path = os.path.join(dest_folder, filename)

    if os.path.isfile(dest_path):
        logger.info("scRNA-seq reference already exists: %s", dest_path)
        return dest_path

    logger.info("Downloading scRNA-seq reference from %s ...", url)
    try:
        urllib.request.urlretrieve(url, dest_path)
        logger.info("Download complete: %s", dest_path)
        return dest_path
    except Exception as exc:
        logger.warning(
            "Failed to download scRNA-seq reference (%s). "
            "Pipeline will continue without it. Error: %s",
            url,
            exc,
        )
        # Clean up partial file
        if os.path.exists(dest_path):
            os.remove(dest_path)
        return None
