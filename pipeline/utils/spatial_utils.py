# Shared spatial utility functions for the Xenium pipeline.

import numpy as np

# Pixel-to-micrometer conversion factors by technology platform.
# Factor = pixels per micrometer; divide pixel distances by this to get µm.
PIXEL_TO_UM_FACTORS = {
    'xenium': 4.70588,
    'cosmx': 8.3333,
    'vizgen': 9.20586,
    'merfish': 9.28,
    'hybriss': 3.11,
    'resolvedbio': 7.24,
}


def pixels_to_um(pixels, technology='xenium'):
    """Convert pixel distances to micrometers for a given technology.

    Parameters
    ----------
    pixels : float or array-like
        Distance(s) in pixels.
    technology : str
        Technology name (case-insensitive). One of: xenium, cosmx, vizgen,
        merfish, hybriss, resolvedbio.

    Returns
    -------
    float or ndarray
        Distance(s) in micrometers.
    """
    factor = PIXEL_TO_UM_FACTORS.get(technology.lower(), 1.0)
    return np.asarray(pixels) / factor


def um_to_pixels(um, technology='xenium'):
    """Convert micrometer distances to pixels for a given technology.

    Parameters
    ----------
    um : float or array-like
        Distance(s) in micrometers.
    technology : str
        Technology name (case-insensitive).

    Returns
    -------
    float or ndarray
        Distance(s) in pixels.
    """
    factor = PIXEL_TO_UM_FACTORS.get(technology.lower(), 1.0)
    return np.asarray(um) * factor
