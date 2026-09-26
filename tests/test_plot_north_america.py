"""Generates the 2024-04-08 total eclipse shapefiles and plots them over a
North America basemap, reading geometry back from the .shp files (not from
in-memory objects), as a smoke test that the shapefiles are valid and
geographically sane.
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shapefile
from pyproj import CRS

from besselian2shape.export import generate_eclipse_shapefiles

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"

# Covers the contiguous US plus enough margin north and south to show the
# eclipse paths' full extent, in (lon, lat) degrees.
_MAP_EXTENT = (-125, -65, 0, 75)

# USA Contiguous Albers Equal Area Conic (USGS version), the map's CRS.
_MAP_CRS = ccrs.Projection(CRS.from_user_input("ESRI:102039"))

# ESRI:102039 is officially registered for CONUS only (area_of_use roughly
# 24-49N), and cartopy's generic Projection wrapper uses that area_of_use as
# the projection's clipping domain -- so without this override, ax.set_extent
# and every drawn feature get silently clipped to that box regardless of
# _MAP_EXTENT, even though the underlying Albers projection math is well
# behaved well outside it. Recompute the domain from _MAP_EXTENT instead.
_lon_min, _lon_max, _lat_min, _lat_max = _MAP_EXTENT
_extent_corners = _MAP_CRS.transform_points(
    ccrs.PlateCarree(),
    np.array([_lon_min, _lon_min, _lon_max, _lon_max]),
    np.array([_lat_min, _lat_max, _lat_max, _lat_min]),
)
_MAP_CRS.bounds = (
    _extent_corners[:, 0].min(), _extent_corners[:, 0].max(),
    _extent_corners[:, 1].min(), _extent_corners[:, 1].max(),
)

# Shapefile coordinates are geographic (WGS84 lon/lat); this tells cartopy
# how to reproject them onto the map's CRS.
_DATA_CRS = ccrs.PlateCarree()


def _plot_shape(ax, reader: shapefile.Reader, **kwargs):
    for shape in reader.shapes():
        parts = list(shape.parts) + [len(shape.points)]
        for start, end in zip(parts, parts[1:]):
            xs, ys = zip(*shape.points[start:end])
            ax.plot(xs, ys, transform=_DATA_CRS, **kwargs)


def test_plot_2024_total_eclipse_shapefiles_on_north_america_map():
    # Written alongside the other eclipse_YYYY-MM-DD/ output directories this
    # project's CLI produces (see .gitignore), so the plot persists for
    # inspection after the test run instead of being deleted with tmp_path.
    OUTPUT_DIR = Path(__file__).parent.parent / "eclipse_2024-04-08"
    written = generate_eclipse_shapefiles(
        2024, 4, 8, OUTPUT_DIR, csv_path=FIXTURE, step_minutes=1.0,
        penumbral_resolution_deg=0.3, penumbral_step_minutes=2.0,
    )
    assert set(written) == {"penumbral_path", "central_line", "umbral_path"}

    fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={"projection": _MAP_CRS})
    plt.tight_layout()
    print(_MAP_EXTENT)
    ax.set_extent(_MAP_EXTENT, crs=_DATA_CRS)
    ax.coastlines(resolution="110m")
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES.with_scale("110m"), linewidth=0.3)

    with shapefile.Reader(str(written["penumbral_path"])) as r:
        _plot_shape(ax, r, color="goldenrod", linewidth=1, label="Penumbral path")
    with shapefile.Reader(str(written["umbral_path"])) as r:
        _plot_shape(ax, r, color="black", linewidth=1.5, label="Umbral path")
    with shapefile.Reader(str(written["central_line"])) as r:
        _plot_shape(ax, r, color="red", linewidth=1, label="Central line")

    ax.set_title("2024-04-08 Total Solar Eclipse")
    ax.legend(loc="lower left")

    out_path = OUTPUT_DIR / "eclipse_2024-04-08_north_america.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    assert out_path.exists()
    assert out_path.stat().st_size > 0

    # The umbral path (path of totality) crosses the US mainland before
    # continuing on into the North Atlantic, so only part of it -- not
    # necessarily all of it -- need fall within the North America extent.
    lon_min, lon_max, lat_min, lat_max = _MAP_EXTENT
    with shapefile.Reader(str(written["umbral_path"])) as r:
        umbral_points = r.shape(0).points
    assert any(
        lon_min < lon < lon_max and lat_min < lat < lat_max
        for lon, lat in umbral_points
    )


def test_plot_2017_total_eclipse_shapefiles_on_north_america_map():
    # Written alongside the other eclipse_YYYY-MM-DD/ output directories this
    # project's CLI produces (see .gitignore), so the plot persists for
    # inspection after the test run instead of being deleted with tmp_path.
    OUTPUT_DIR = Path(__file__).parent.parent / "eclipse_2017-08-21"
    written = generate_eclipse_shapefiles(
        2017, 8, 21, OUTPUT_DIR, csv_path=FIXTURE, step_minutes=1.0,
        penumbral_resolution_deg=0.3, penumbral_step_minutes=2.0,
    )
    assert set(written) == {"penumbral_path", "central_line", "umbral_path"}

    fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={"projection": _MAP_CRS})
    plt.tight_layout()

    ax.set_extent(_MAP_EXTENT, crs=_DATA_CRS)
    ax.coastlines(resolution="110m")
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES.with_scale("110m"), linewidth=0.3)

    with shapefile.Reader(str(written["penumbral_path"])) as r:
        _plot_shape(ax, r, color="goldenrod", linewidth=1, label="Penumbral path")
    with shapefile.Reader(str(written["umbral_path"])) as r:
        _plot_shape(ax, r, color="black", linewidth=1.5, label="Umbral path")
    with shapefile.Reader(str(written["central_line"])) as r:
        _plot_shape(ax, r, color="red", linewidth=1, label="Central line")

    ax.set_title("2017-08-21 Total Solar Eclipse")
    ax.legend(loc="lower left")

    out_path = OUTPUT_DIR / "eclipse_2017-08-17_north_america.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    assert out_path.exists()
    assert out_path.stat().st_size > 0

    # The umbral path (path of totality) crosses the US mainland before
    # continuing on into the North Atlantic, so only part of it -- not
    # necessarily all of it -- need fall within the North America extent.
    lon_min, lon_max, lat_min, lat_max = _MAP_EXTENT
    with shapefile.Reader(str(written["umbral_path"])) as r:
        umbral_points = r.shape(0).points
    assert any(
        lon_min < lon < lon_max and lat_min < lat < lat_max
        for lon, lat in umbral_points
    )
