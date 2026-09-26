from pathlib import Path

import shapefile
import pytest

from besselian2shape.export import generate_eclipse_shapefiles

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"

# Coarse/fast settings for the penumbral raster scan in tests.
FAST_RASTER = {"penumbral_resolution_deg": 3.0, "penumbral_step_minutes": 20.0}


def test_generate_shapefiles_for_total_eclipse(tmp_path):
    result = generate_eclipse_shapefiles(2024, 4, 8, tmp_path, csv_path=FIXTURE, step_minutes=2, **FAST_RASTER)

    assert set(result) == {"penumbral_path", "central_line", "umbral_path"}
    for name, shp_path in result.items():
        assert shp_path == tmp_path / f"{name}.shp"
        assert shp_path.exists()
        assert shp_path.with_suffix(".shx").exists()
        assert shp_path.with_suffix(".dbf").exists()
        assert shp_path.with_suffix(".prj").exists()

    central = shapefile.Reader(str(result["central_line"]))
    assert central.shapeType == shapefile.POLYLINE
    record = central.record(0)
    assert record["date"] == "2024-04-08"
    assert record["ecl_type"] == "T"

    umbral = shapefile.Reader(str(result["umbral_path"]))
    assert umbral.shapeType == shapefile.POLYGON
    shape = umbral.shape(0)
    assert shape.points[0] == pytest.approx(shape.points[-1])


def test_generate_shapefiles_for_partial_eclipse(tmp_path):
    result = generate_eclipse_shapefiles(2025, 3, 29, tmp_path, csv_path=FIXTURE, step_minutes=2, **FAST_RASTER)

    assert set(result) == {"penumbral_path"}
    penumbral = shapefile.Reader(str(result["penumbral_path"]))
    assert penumbral.shapeType == shapefile.POLYGON
    record = penumbral.record(0)
    assert record["ecl_type"] == "P"


def test_generate_shapefiles_raises_for_missing_date(tmp_path):
    with pytest.raises(LookupError):
        generate_eclipse_shapefiles(2024, 4, 9, tmp_path, csv_path=FIXTURE, **FAST_RASTER)
