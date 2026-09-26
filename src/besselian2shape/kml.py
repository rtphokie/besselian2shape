"""KML/KMZ writing for eclipse geometry.

Builds directly on the same (lat, lon) geometry produced for shapefiles
(see `geometry.py`) rather than converting shapefiles after the fact --
there's no intermediate format to round-trip through, and no dependency
on GDAL/ogr2ogr. KMZ is just the KML document zip-compressed, per the
KML spec.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from .elements import BesselianElements

_KML_NS = "http://www.opengis.net/kml/2.2"

# KML color format is aabbggrr (alpha, blue, green, red), not rrggbb.
_CENTRAL_LINE_COLOR = "ff0000ff"  # opaque red
_UMBRAL_FILL_COLOR = "7f000000"  # ~50% black
_PENUMBRAL_FILL_COLOR = "4c00ffff"  # ~30% yellow

_CENTRAL_STYLE_ID = "centralLineStyle"
_UMBRAL_STYLE_ID = "umbralPathStyle"
_PENUMBRAL_STYLE_ID = "penumbralPathStyle"


def _coordinates_text(points_latlon: list[tuple[float, float]]) -> str:
    return " ".join(f"{lon:.6f},{lat:.6f},0" for lat, lon in points_latlon)


def _extended_data(fields: dict[str, str]) -> ET.Element:
    ext = ET.Element("ExtendedData")
    for name, value in fields.items():
        data = ET.SubElement(ext, "Data", {"name": name})
        ET.SubElement(data, "value").text = str(value)
    return ext


def _line_style(style_id: str, color: str, width: float) -> ET.Element:
    style = ET.Element("Style", {"id": style_id})
    line = ET.SubElement(style, "LineStyle")
    ET.SubElement(line, "color").text = color
    ET.SubElement(line, "width").text = str(width)
    return style


def _polygon_style(style_id: str, fill_color: str) -> ET.Element:
    style = ET.Element("Style", {"id": style_id})
    line = ET.SubElement(style, "LineStyle")
    ET.SubElement(line, "color").text = "ff" + fill_color[2:]  # opaque outline, same hue
    poly = ET.SubElement(style, "PolyStyle")
    ET.SubElement(poly, "color").text = fill_color
    ET.SubElement(poly, "fill").text = "1"
    ET.SubElement(poly, "outline").text = "1"
    return style


def _line_placemark(
    name: str, points_latlon: list[tuple[float, float]], style_id: str, fields: dict[str, str]
) -> ET.Element:
    placemark = ET.Element("Placemark")
    ET.SubElement(placemark, "name").text = name
    ET.SubElement(placemark, "styleUrl").text = f"#{style_id}"
    placemark.append(_extended_data(fields))
    line_string = ET.SubElement(placemark, "LineString")
    ET.SubElement(line_string, "tessellate").text = "1"
    ET.SubElement(line_string, "coordinates").text = _coordinates_text(points_latlon)
    return placemark


def _polygon_geometry(
    outer_latlon: list[tuple[float, float]],
    holes_latlon: list[list[tuple[float, float]]] | None = None,
) -> ET.Element:
    polygon = ET.Element("Polygon")
    ET.SubElement(polygon, "tessellate").text = "1"
    outer = ET.SubElement(polygon, "outerBoundaryIs")
    ring = ET.SubElement(outer, "LinearRing")
    ET.SubElement(ring, "coordinates").text = _coordinates_text(outer_latlon)
    for hole in holes_latlon or []:
        inner = ET.SubElement(polygon, "innerBoundaryIs")
        hole_ring = ET.SubElement(inner, "LinearRing")
        ET.SubElement(hole_ring, "coordinates").text = _coordinates_text(hole)
    return polygon


def _polygon_placemark(
    name: str,
    rings_latlon: list[tuple[float, float]] | list[list[tuple[float, float]]],
    style_id: str,
    fields: dict[str, str],
) -> ET.Element:
    """Placemark for one polygon, or (as a <MultiGeometry> of several
    <Polygon>s) a list of separate closed rings -- e.g. a main region plus
    separate disjoint loops, as can happen near a pole.
    """
    if rings_latlon and isinstance(rings_latlon[0], tuple):
        rings_latlon = [rings_latlon]  # a single ring, not a list of rings

    placemark = ET.Element("Placemark")
    ET.SubElement(placemark, "name").text = name
    ET.SubElement(placemark, "styleUrl").text = f"#{style_id}"
    placemark.append(_extended_data(fields))
    if len(rings_latlon) == 1:
        placemark.append(_polygon_geometry(rings_latlon[0]))
    else:
        multi = ET.SubElement(placemark, "MultiGeometry")
        for ring in rings_latlon:
            multi.append(_polygon_geometry(ring))
    return placemark


def _polygon_placemark_with_holes(
    name: str,
    polygons: list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]],
    style_id: str,
    fields: dict[str, str],
) -> ET.Element:
    """Placemark for one polygon, or (as a <MultiGeometry> of several
    <Polygon>s) a list of (outer_ring, [hole_ring, ...]) polygons -- each
    hole becomes an <innerBoundaryIs> nested within its own outer ring's
    <Polygon>, so viewers render it as a gap rather than an extra filled
    region.
    """
    placemark = ET.Element("Placemark")
    ET.SubElement(placemark, "name").text = name
    ET.SubElement(placemark, "styleUrl").text = f"#{style_id}"
    placemark.append(_extended_data(fields))
    if len(polygons) == 1:
        outer, holes = polygons[0]
        placemark.append(_polygon_geometry(outer, holes))
    else:
        multi = ET.SubElement(placemark, "MultiGeometry")
        for outer, holes in polygons:
            multi.append(_polygon_geometry(outer, holes))
    return placemark


def _common_fields(e: BesselianElements) -> dict[str, str]:
    return {
        "date": f"{e.year:04d}-{e.month:02d}-{e.day:02d}",
        "eclipse_type": e.eclipse_type,
        "saros": str(e.saros),
        "gamma": f"{e.gamma:.6f}",
        "magnitude": f"{e.magnitude:.6f}",
    }


def build_document(e: BesselianElements) -> ET.Element:
    """Build the root <kml> element with styles and an (initially empty)
    <Document>, ready for placemarks to be appended.
    """
    kml = ET.Element("kml", {"xmlns": _KML_NS})
    document = ET.SubElement(kml, "Document")
    ET.SubElement(document, "name").text = (
        f"Solar Eclipse {e.year:04d}-{e.month:02d}-{e.day:02d}"
    )
    document.append(_line_style(_CENTRAL_STYLE_ID, _CENTRAL_LINE_COLOR, 3))
    document.append(_polygon_style(_UMBRAL_STYLE_ID, _UMBRAL_FILL_COLOR))
    document.append(_polygon_style(_PENUMBRAL_STYLE_ID, _PENUMBRAL_FILL_COLOR))
    return kml


def document_element(kml: ET.Element) -> ET.Element:
    return kml.find("Document")


def add_central_line(kml: ET.Element, e: BesselianElements, points_latlon: list[tuple[float, float]]) -> None:
    fields = {**_common_fields(e), "duration": e.central_duration, "path_km": f"{e.path_width:.2f}"}
    document_element(kml).append(
        _line_placemark("Central Line", points_latlon, _CENTRAL_STYLE_ID, fields)
    )


def add_umbral_path(kml: ET.Element, e: BesselianElements, ring_latlon: list[tuple[float, float]]) -> None:
    fields = {**_common_fields(e), "duration": e.central_duration, "path_km": f"{e.path_width:.2f}"}
    document_element(kml).append(
        _polygon_placemark("Umbral Path", ring_latlon, _UMBRAL_STYLE_ID, fields)
    )


def add_penumbral_path(
    kml: ET.Element,
    e: BesselianElements,
    polygons: list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]],
) -> None:
    document_element(kml).append(
        _polygon_placemark_with_holes(
            "Penumbral Path", polygons, _PENUMBRAL_STYLE_ID, _common_fields(e)
        )
    )


def write_kml(kml: ET.Element, path: Path) -> Path:
    """Write the KML tree to `path`. If `path` ends in .kmz, it is written
    as a zip-compressed KMZ archive (containing "doc.kml") instead.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(kml, space="  ")
    xml_bytes = b'<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(kml, encoding="UTF-8")

    if path.suffix.lower() == ".kmz":
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("doc.kml", xml_bytes)
    else:
        path.write_bytes(xml_bytes)
    return path
