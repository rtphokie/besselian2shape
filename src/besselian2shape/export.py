"""Top-level API: generate ESRI shapefiles or KML/KMZ for a solar eclipse
by date.
"""

from __future__ import annotations

from pathlib import Path

from . import kml as kml_writer
from .elements import BesselianElements, find_by_date
from .geometry import central_line, is_central, limits_to_ring, umbral_limits
from .raster import penumbral_coverage_rings
from .shapefiles import (
    write_polygon_shapefile,
    write_polygon_shapefile_with_holes,
    write_polyline_shapefile,
)

# Grid resolution and time step for the penumbral visibility raster scan
# (see raster.py). The path/umbral computation's own `step_minutes` isn't
# reused here -- it controls sampling density along a 1-D analytic curve,
# a different kind of parameter from a 2-D grid's resolution.
_PENUMBRAL_RESOLUTION_DEG = 0.25
_PENUMBRAL_STEP_MINUTES = 1.5


class _Layers:
    def __init__(
        self,
        e: BesselianElements,
        step_minutes: float,
        penumbral_resolution_deg: float,
        penumbral_step_minutes: float,
    ):
        self.penumbral_polygons = penumbral_coverage_rings(
            e,
            lambda s: (s.l1, s.l1p),
            resolution_deg=penumbral_resolution_deg,
            step_minutes=penumbral_step_minutes,
        )
        self.central_line: list[tuple[float, float]] = []
        self.umbral_ring: list[tuple[float, float]] = []
        if is_central(e):
            self.central_line = central_line(e, step_minutes)
            self.umbral_ring = limits_to_ring(umbral_limits(e, step_minutes))


def _load_layers(
    year: int,
    month: int,
    day: int,
    csv_path: Path | None,
    step_minutes: float,
    elements: BesselianElements | None,
    penumbral_resolution_deg: float,
    penumbral_step_minutes: float,
) -> tuple[BesselianElements, _Layers]:
    e = elements if elements is not None else find_by_date(year, month, day, csv_path=csv_path)
    layers = _Layers(e, step_minutes, penumbral_resolution_deg, penumbral_step_minutes)
    if not layers.penumbral_polygons:
        raise RuntimeError(
            f"Could not compute a penumbral path for {year:04d}-{month:02d}-{day:02d}"
        )
    return e, layers


def generate_eclipse_shapefiles(
    year: int,
    month: int,
    day: int,
    output_dir: Path | str,
    csv_path: Path | None = None,
    step_minutes: float = 0.5,
    elements: BesselianElements | None = None,
    penumbral_resolution_deg: float = _PENUMBRAL_RESOLUTION_DEG,
    penumbral_step_minutes: float = _PENUMBRAL_STEP_MINUTES,
) -> dict[str, Path]:
    """Generate ESRI shapefiles for the solar eclipse on the given
    (proleptic Gregorian, astronomical-numbered) calendar date.

    Writes into `output_dir`:
      - "penumbral_path.shp": polygon (possibly multi-part, e.g. a main
        region plus separate disjoint loops for a shadow path that passes
        close to a pole) of the region on Earth from which at least a
        partial eclipse is visible (all eclipse types).
      - "central_line.shp": polyline of the path of totality/annularity's
        central line (total/annular/hybrid eclipses only).
      - "umbral_path.shp": polygon of the path of totality/annularity
        (total/annular/hybrid eclipses only).

    Each shapefile is written as a single-feature .shp/.shx/.dbf/.prj set
    in WGS84 geographic coordinates. `step_minutes` controls how finely
    the umbral/central-line path is sampled along its length.
    `penumbral_resolution_deg`/`penumbral_step_minutes` control the
    coverage-grid resolution and time step used to compute the penumbral
    boundary (see `raster.penumbral_coverage_rings`); the defaults are a
    reasonable accuracy/runtime balance, coarsened here mainly for fast
    tests.

    If `elements` is given, it is used directly instead of looking `year`,
    `month`, `day` up in the (cached) CSV -- useful when generating more
    than one output format for the same eclipse, to avoid re-parsing the
    dataset each time.

    Raises LookupError if no eclipse occurred on that date. Returns a dict
    mapping each layer name written to its .shp path.
    """
    e, layers = _load_layers(
        year, month, day, csv_path, step_minutes, elements,
        penumbral_resolution_deg, penumbral_step_minutes,
    )
    output_dir = Path(output_dir)

    written: dict[str, Path] = {
        "penumbral_path": write_polygon_shapefile_with_holes(
            output_dir / "penumbral_path.shp", layers.penumbral_polygons, e
        )
    }

    if layers.central_line or layers.umbral_ring:
        path_fields = [("duration", "C", 8, 0), ("path_km", "F", 10, 2)]
        path_record = {"duration": e.central_duration, "path_km": e.path_width}

        if layers.central_line:
            written["central_line"] = write_polyline_shapefile(
                output_dir / "central_line.shp",
                layers.central_line,
                e,
                extra_fields=path_fields,
                extra_record=path_record,
            )
        if layers.umbral_ring:
            written["umbral_path"] = write_polygon_shapefile(
                output_dir / "umbral_path.shp",
                layers.umbral_ring,
                e,
                extra_fields=path_fields,
                extra_record=path_record,
            )

    return written


def generate_eclipse_kml(
    year: int,
    month: int,
    day: int,
    output_path: Path | str,
    csv_path: Path | None = None,
    step_minutes: float = 0.5,
    elements: BesselianElements | None = None,
    penumbral_resolution_deg: float = _PENUMBRAL_RESOLUTION_DEG,
    penumbral_step_minutes: float = _PENUMBRAL_STEP_MINUTES,
) -> Path:
    """Generate a single KML (or, if `output_path` ends in .kmz, KMZ) file
    for the solar eclipse on the given date, containing the same layers as
    `generate_eclipse_shapefiles` as separate Placemarks in one Document:
    the penumbral path (all eclipse types), and, for total/annular/hybrid
    eclipses, the central line and umbral path.

    See `generate_eclipse_shapefiles` for the other parameters.

    Raises LookupError if no eclipse occurred on that date.
    """
    e, layers = _load_layers(
        year, month, day, csv_path, step_minutes, elements,
        penumbral_resolution_deg, penumbral_step_minutes,
    )
    output_path = Path(output_path)

    doc = kml_writer.build_document(e)
    kml_writer.add_penumbral_path(doc, e, layers.penumbral_polygons)
    if layers.central_line:
        kml_writer.add_central_line(doc, e, layers.central_line)
    if layers.umbral_ring:
        kml_writer.add_umbral_path(doc, e, layers.umbral_ring)

    return kml_writer.write_kml(doc, output_path)
