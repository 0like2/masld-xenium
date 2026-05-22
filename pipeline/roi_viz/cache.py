"""Disk cache layer for the roi_viz pipeline.

Goals:
- Heavy computations (ROI crop, label-map → polygon, difference mask, metric)
  are persisted to disk so partial edits do not re-trigger the whole pipeline.
- Cache fingerprints depend on input file mtimes/sizes/hashes, so external
  changes invalidate downstream caches automatically.
- Stage names ('raw', 'boundaries', 'diff', 'metrics', 'render') let the CLI
  selectively force a single stage.

Storage:
- pickled .pkl for arbitrary Python objects (DataFrame, dict, ndarray bundles)
- the caller is encouraged to pre-cast to a pickle-friendly form if speed
  matters; for very large arrays use np.savez_compressed elsewhere.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Iterable, List, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Globals (set up via configure())
# ---------------------------------------------------------------------------
_ROOT: Optional[Path] = None
_FORCE_STAGES: Set[str] = set()           # stages whose cache is bypassed
_FORCE_ROIS: Set[str] = set()             # if non-empty, only invalidate these ROIs
_VERSION: str = "0.1.0"
_DISABLE: bool = False


def configure(root: Path,
              force_stages: Optional[Iterable[str]] = None,
              force_rois: Optional[Iterable[str]] = None,
              version: Optional[str] = None,
              disable: bool = False) -> None:
    """Set global cache configuration. Call once at pipeline entry."""
    global _ROOT, _FORCE_STAGES, _FORCE_ROIS, _VERSION, _DISABLE
    _ROOT = Path(root)
    _ROOT.mkdir(parents=True, exist_ok=True)
    _FORCE_STAGES = set(force_stages or [])
    _FORCE_ROIS = set(force_rois or [])
    if version:
        _VERSION = version
    _DISABLE = bool(disable)
    logger.info(f"Cache configured: root={_ROOT}, force_stages={_FORCE_STAGES}, "
                f"force_rois={_FORCE_ROIS}, disable={_DISABLE}")


def get_root() -> Path:
    if _ROOT is None:
        raise RuntimeError("cache.configure() must be called first")
    return _ROOT


# ---------------------------------------------------------------------------
# Fingerprint computation
# ---------------------------------------------------------------------------
def _file_fingerprint(path: Path) -> str:
    p = Path(path)
    if not p.exists():
        return f"missing:{p}"
    st = p.stat()
    h = hashlib.blake2b(digest_size=12)
    # mtime + size are usually enough; add a small content hash for safety.
    h.update(str(st.st_size).encode())
    h.update(str(int(st.st_mtime_ns // 1_000_000)).encode())
    try:
        with open(p, "rb") as f:
            h.update(f.read(1 << 20))  # first 1 MB
    except OSError:
        pass
    return h.hexdigest()


def compute_fingerprint(deps: Sequence[Any]) -> str:
    """Hash a list of dependencies (file paths, scalars, dicts) into a short key."""
    h = hashlib.blake2b(digest_size=16)
    h.update(_VERSION.encode())
    for d in deps:
        if isinstance(d, (str, Path)) and Path(d).exists():
            h.update(_file_fingerprint(Path(d)).encode())
        elif isinstance(d, (str, Path)):
            h.update(f"path:{d}".encode())
        else:
            try:
                h.update(json.dumps(d, sort_keys=True, default=str).encode())
            except Exception:
                h.update(repr(d).encode())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Cache directory layout
# ---------------------------------------------------------------------------
def stage_dir(stage: str) -> Path:
    d = get_root() / "_cache" / stage
    d.mkdir(parents=True, exist_ok=True)
    return d


def cache_path(stage: str, name: str, fingerprint: str,
               ext: str = "pkl") -> Path:
    return stage_dir(stage) / f"{name}__{fingerprint[:16]}.{ext}"


def _atomic_write_pickle(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.stem + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def _force_invalidate(stage: str, roi_id: Optional[str]) -> bool:
    if _DISABLE:
        return True
    if "all" in _FORCE_STAGES or stage in _FORCE_STAGES:
        if not _FORCE_ROIS:
            return True
        if roi_id and roi_id in _FORCE_ROIS:
            return True
    return False


# ---------------------------------------------------------------------------
# Decorator
# ---------------------------------------------------------------------------
def cached(stage: str, name_fn: Callable[..., str],
           deps_fn: Callable[..., Sequence[Any]],
           roi_id_fn: Optional[Callable[..., Optional[str]]] = None,
           ext: str = "pkl") -> Callable:
    """Disk-cache decorator.

    Parameters
    ----------
    stage : str
        Logical stage name ('raw', 'boundaries', 'diff', 'metrics', 'render').
    name_fn : callable
        Returns the human-readable cache base name from the wrapped function's
        args (used in the on-disk filename, prefixing the fingerprint).
    deps_fn : callable
        Returns the list of dependencies (file paths, scalars, dicts) that
        contribute to the fingerprint.
    roi_id_fn : callable, optional
        Returns the ROI id from args, used for ``--force-roi`` filtering.
    ext : str
        File extension (default 'pkl').
    """

    def decorator(fn: Callable) -> Callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            try:
                name = name_fn(*args, **kwargs)
                deps = list(deps_fn(*args, **kwargs))
                roi_id = roi_id_fn(*args, **kwargs) if roi_id_fn else None
            except Exception as e:
                logger.warning(f"cache key computation failed for {fn.__name__}: {e}; "
                               "executing without cache")
                return fn(*args, **kwargs)

            fp = compute_fingerprint(deps + [fn.__name__])
            path = cache_path(stage, name, fp, ext=ext)

            if not _force_invalidate(stage, roi_id) and path.exists():
                try:
                    with open(path, "rb") as f:
                        obj = pickle.load(f)
                    logger.debug(f"cache hit ({stage}): {path.name}")
                    return obj
                except Exception as e:
                    logger.warning(f"cache load failed for {path}: {e}; recomputing")

            obj = fn(*args, **kwargs)
            try:
                _atomic_write_pickle(path, obj)
                logger.debug(f"cache write ({stage}): {path.name}")
            except Exception as e:
                logger.warning(f"cache write failed for {path}: {e}")
            return obj

        return wrapper

    return decorator


# ---------------------------------------------------------------------------
# PNG output cache (mtime-based, no fingerprint file)
# ---------------------------------------------------------------------------
def png_is_fresh(out_path: Path, dep_paths: Sequence[Path]) -> bool:
    """Return True iff *out_path* exists and is newer than every dep.

    When any non-render stage is in ``_FORCE_STAGES`` (e.g.
    ``--force-cache=diff``), we *also* require the PNG to be newer than
    every cached artifact in that stage's directory — so invalidating an
    upstream stage cascades into PNG regeneration.
    """
    out = Path(out_path)
    if not out.exists():
        return False
    if "render" in _FORCE_STAGES or "all" in _FORCE_STAGES or _DISABLE:
        return False
    out_m = out.stat().st_mtime
    for d in dep_paths:
        d = Path(d)
        if d.exists() and d.stat().st_mtime > out_m:
            return False
    # Cascade only when an upstream stage was explicitly forced.
    upstream_forced = _FORCE_STAGES - {"render"}
    if upstream_forced and _ROOT is not None:
        cache_root = _ROOT / "_cache"
        for stage in upstream_forced:
            stage_path = cache_root / stage
            if not stage_path.exists():
                continue
            try:
                newest = max((p.stat().st_mtime for p in stage_path.iterdir()
                                if p.is_file()), default=0.0)
            except OSError:
                newest = 0.0
            if newest > out_m:
                return False
    return True
