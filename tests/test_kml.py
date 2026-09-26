import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pytest

from besselian2shape.export import generate_eclipse_kml

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"
NS = {"k": "http://www.opengis.net/kml/2.2"}

# Coarse/fast settings for the penumbral raster scan in tests.
FAST_RASTER = {"penumbral_resolution_deg": 3.0, "penumbral_step_minutes": 20.0}


def _placemark_names(root):
    return [p.find("k:name", NS).text for p in root.findall(".//k:Placemark", NS)]


def test_generate_kml_for_total_eclipse(tmp_path):
    path = generate_eclipse_kml(
        2024, 4, 8, tmp_path / "eclipse.kml", csv_path=FIXTURE, step_minutes=2, **FAST_RASTER
    )

    assert path.exists()
    root = ET.parse(path).getroot()
    assert _placemark_names(root) == ["Penumbral Path", "Central Line", "Umbral Path"]

    central = root.find(".//k:Placemark[k:name='Central Line']", NS)
    coords = central.find(".//k:coordinates", NS).text.split()
    assert len(coords) > 50
    lon, lat, alt = coords[0].split(",")
    float(lon), float(lat), float(alt)  # parses as numbers

    fields = {
        d.get("name"): d.find("k:value", NS).text
        for d in central.find("k:ExtendedData", NS).findall("k:Data", NS)
    }
    assert fields["date"] == "2024-04-08"
    assert fields["eclipse_type"] == "T"


def test_generate_kml_for_partial_eclipse_has_only_penumbral_layer(tmp_path):
    path = generate_eclipse_kml(
        2025, 3, 29, tmp_path / "eclipse.kml", csv_path=FIXTURE, step_minutes=2, **FAST_RASTER
    )

    root = ET.parse(path).getroot()
    assert _placemark_names(root) == ["Penumbral Path"]


def test_generate_kmz_is_a_valid_zip_containing_doc_kml(tmp_path):
    path = generate_eclipse_kml(
        2024, 4, 8, tmp_path / "eclipse.kmz", csv_path=FIXTURE, step_minutes=2, **FAST_RASTER
    )

    assert path.suffix == ".kmz"
    with zipfile.ZipFile(path) as zf:
        assert zf.namelist() == ["doc.kml"]
        inner_root = ET.fromstring(zf.read("doc.kml"))
    assert _placemark_names(inner_root) == ["Penumbral Path", "Central Line", "Umbral Path"]


def test_generate_kml_matches_shapefile_geometry_point_counts(tmp_path):
    from besselian2shape.export import generate_eclipse_shapefiles
    import shapefile

    shp_result = generate_eclipse_shapefiles(
        2024, 4, 8, tmp_path / "shp", csv_path=FIXTURE, step_minutes=2, **FAST_RASTER
    )
    kml_path = generate_eclipse_kml(
        2024, 4, 8, tmp_path / "eclipse.kml", csv_path=FIXTURE, step_minutes=2, **FAST_RASTER
    )
    root = ET.parse(kml_path).getroot()

    for layer, placemark_name in [
        ("central_line", "Central Line"),
        ("umbral_path", "Umbral Path"),
        ("penumbral_path", "Penumbral Path"),
    ]:
        shp_points = shapefile.Reader(str(shp_result[layer])).shape(0).points
        placemark = root.find(f".//k:Placemark[k:name='{placemark_name}']", NS)
        # Sums across every <coordinates> element (a polygon may have
        # multiple parts, as a MultiGeometry of several <Polygon>s in KML).
        kml_points = sum(
            len(el.text.split()) for el in placemark.findall(".//k:coordinates", NS)
        )
        assert len(shp_points) == kml_points


def test_generate_kml_raises_for_missing_date(tmp_path):
    with pytest.raises(LookupError):
        generate_eclipse_kml(2024, 4, 9, tmp_path / "eclipse.kml", csv_path=FIXTURE, **FAST_RASTER)
