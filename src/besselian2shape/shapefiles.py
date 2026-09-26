"""Low-level ESRI shapefile writing for eclipse geometry."""

from __future__ import annotations

from pathlib import Path

import shapefile as pyshp

from .elements import BesselianElements

_PRJ_WGS84 = (
    'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
    'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
    'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'
)


def _write_prj(path: Path) -> None:
    path.write_text(_PRJ_WGS84)


def _signed_area(ring_lonlat: list[list[float]]) -> float:
    total = 0.0
    for (x1, y1), (x2, y2) in zip(ring_lonlat, ring_lonlat[1:]):
        total += x1 * y2 - x2 * y1
    return total / 2.0


def _to_wound_lonlat(
    ring_latlon: list[tuple[float, float]], clockwise: bool = True
) -> list[list[float]]:
    """Convert a (lat, lon) ring to [lon, lat] pairs, reordered clockwise or
    counter-clockwise as requested. The shapefile spec requires exterior
    (outer) polygon rings to be clockwise and interior (hole) rings
    counter-clockwise -- readers use that winding, not anything else about
    the ring, to tell which is which.
    """
    lonlat = [[lon, lat] for lat, lon in ring_latlon]
    is_clockwise = _signed_area(lonlat) < 0
    if is_clockwise != clockwise:
        lonlat.reverse()
    return lonlat


def _common_fields() -> list[tuple]:
    return [
        ("date", "C", 10, 0),
        ("ecl_type", "C", 4, 0),
        ("saros", "N", 6, 0),
        ("gamma", "F", 10, 6),
        ("magnitude", "F", 10, 6),
    ]


def _common_record(e: BesselianElements) -> dict:
    return {
        "date": f"{e.year:04d}-{e.month:02d}-{e.day:02d}",
        "ecl_type": e.eclipse_type,
        "saros": e.saros,
        "gamma": e.gamma,
        "magnitude": e.magnitude,
    }


def write_polygon_shapefile(
    path: Path,
    rings_latlon: list[tuple[float, float]] | list[list[tuple[float, float]]],
    e: BesselianElements,
    extra_fields: list[tuple] | None = None,
    extra_record: dict | None = None,
) -> Path:
    """Write a single-feature (possibly multi-part) polygon shapefile from
    one closed (lat, lon) ring, or a list of them (e.g. a main region plus
    separate disjoint loops, as can happen near a pole).
    """
    if rings_latlon and isinstance(rings_latlon[0], tuple):
        rings_latlon = [rings_latlon]  # a single ring, not a list of rings

    path.parent.mkdir(parents=True, exist_ok=True)
    with pyshp.Writer(str(path), shapeType=pyshp.POLYGON) as w:
        for spec in _common_fields() + (extra_fields or []):
            w.field(*spec)
        w.poly([_to_wound_lonlat(ring) for ring in rings_latlon])
        w.record(**{**_common_record(e), **(extra_record or {})})
    _write_prj(path.with_suffix(".prj"))
    return path


def write_polygon_shapefile_with_holes(
    path: Path,
    polygons: list[
        tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]
    ],
    e: BesselianElements,
    extra_fields: list[tuple] | None = None,
    extra_record: dict | None = None,
) -> Path:
    """Write a single-feature (possibly multi-part) polygon shapefile from
    a list of (outer_ring, [hole_ring, ...]) polygons -- e.g. a main
    region plus separate disjoint loops and/or gaps, as can happen near a
    pole. Outer rings are wound clockwise and holes counter-clockwise, per
    the shapefile spec, so readers correctly render holes as gaps rather
    than extra filled regions.
    """
    parts: list[list[list[float]]] = []
    for outer, holes in polygons:
        parts.append(_to_wound_lonlat(outer, clockwise=True))
        for hole in holes:
            parts.append(_to_wound_lonlat(hole, clockwise=False))

    path.parent.mkdir(parents=True, exist_ok=True)
    with pyshp.Writer(str(path), shapeType=pyshp.POLYGON) as w:
        for spec in _common_fields() + (extra_fields or []):
            w.field(*spec)
        w.poly(parts)
        w.record(**{**_common_record(e), **(extra_record or {})})
    _write_prj(path.with_suffix(".prj"))
    return path


def write_polyline_shapefile(
    path: Path,
    line_latlon: list[tuple[float, float]],
    e: BesselianElements,
    extra_fields: list[tuple] | None = None,
    extra_record: dict | None = None,
) -> Path:
    """Write a single-feature polyline shapefile from an ordered (lat, lon) path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with pyshp.Writer(str(path), shapeType=pyshp.POLYLINE) as w:
        for spec in _common_fields() + (extra_fields or []):
            w.field(*spec)
        w.line([[[lon, lat] for lat, lon in line_latlon]])
        w.record(**{**_common_record(e), **(extra_record or {})})
    _write_prj(path.with_suffix(".prj"))
    return path
