"""Download and local caching of NASA's Besselian elements dataset."""

from __future__ import annotations

import shutil
from pathlib import Path

import requests
from platformdirs import user_cache_dir

BESSELIAN_ELEMENTS_URL = (
    "https://eclipse.gsfc.nasa.gov/eclipse_besselian_from_mysqldump2.csv"
)

APP_NAME = "besselian2shape"


def get_cache_dir() -> Path:
    """Return the local cache directory, creating it if necessary."""
    cache_dir = Path(user_cache_dir(APP_NAME))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def get_cached_csv_path() -> Path:
    """Return the path where the Besselian elements CSV is (or will be) cached."""
    return get_cache_dir() / "eclipse_besselian_from_mysqldump2.csv"


def download_besselian_csv(force: bool = False) -> Path:
    """Ensure the NASA Besselian elements CSV is cached locally, downloading it
    if it is not already present (or if `force` is True).

    Returns the path to the cached file.
    """
    dest = get_cached_csv_path()
    if dest.exists() and not force:
        return dest

    response = requests.get(BESSELIAN_ELEMENTS_URL, stream=True, timeout=60)
    response.raise_for_status()

    tmp_path = dest.with_suffix(dest.suffix + ".part")
    try:
        with open(tmp_path, "wb") as f:
            shutil.copyfileobj(response.raw, f)
        tmp_path.replace(dest)
    finally:
        tmp_path.unlink(missing_ok=True)

    return dest
