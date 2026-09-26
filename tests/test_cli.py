from pathlib import Path

import pytest

from besselian2shape import cli

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"

# Coarse/fast settings for the penumbral raster scan, so tests don't pay
# for production-quality resolution.
FAST_RASTER = ["--penumbral-resolution", "3", "--penumbral-step-minutes", "20"]


@pytest.fixture(autouse=True)
def _use_fixture_csv(monkeypatch):
    # Point find_by_date (used inside cli.main) at the small test fixture
    # instead of the full cached dataset, and skip the network entirely.
    monkeypatch.setattr(cli, "find_by_date", lambda y, m, d: _find(y, m, d))
    yield


def _find(year, month, day):
    from besselian2shape.elements import find_by_date

    return find_by_date(year, month, day, csv_path=FIXTURE)


def test_cli_generates_shapefiles_by_default(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rc = cli.main(["2024", "4", "8", *FAST_RASTER])

    assert rc == 0
    out_dir = tmp_path / "eclipse_2024-04-08"
    assert (out_dir / "penumbral_path.shp").exists()
    assert (out_dir / "central_line.shp").exists()
    assert (out_dir / "umbral_path.shp").exists()
    assert not (out_dir / "eclipse.kml").exists()


def test_cli_custom_output_dir_and_kml_format(tmp_path):
    out_dir = tmp_path / "custom"
    rc = cli.main(["2024", "4", "8", "-o", str(out_dir), "-f", "kml", *FAST_RASTER])

    assert rc == 0
    assert (out_dir / "eclipse.kml").exists()
    assert not (out_dir / "penumbral_path.shp").exists()


def test_cli_all_formats(tmp_path):
    out_dir = tmp_path / "all"
    rc = cli.main(["2024", "4", "8", "-o", str(out_dir), "-f", "all", *FAST_RASTER])

    assert rc == 0
    assert (out_dir / "penumbral_path.shp").exists()
    assert (out_dir / "eclipse.kml").exists()
    assert (out_dir / "eclipse.kmz").exists()


def test_cli_repeated_format_flags(tmp_path):
    out_dir = tmp_path / "repeated"
    rc = cli.main(["2024", "4", "8", "-o", str(out_dir), "-f", "kml", "-f", "kmz", *FAST_RASTER])

    assert rc == 0
    assert (out_dir / "eclipse.kml").exists()
    assert (out_dir / "eclipse.kmz").exists()
    assert not (out_dir / "penumbral_path.shp").exists()


def test_cli_partial_eclipse_only_writes_penumbral_layer(tmp_path):
    out_dir = tmp_path / "partial"
    rc = cli.main(["2025", "3", "29", "-o", str(out_dir), "-f", "shp", *FAST_RASTER])

    assert rc == 0
    assert (out_dir / "penumbral_path.shp").exists()
    assert not (out_dir / "central_line.shp").exists()
    assert not (out_dir / "umbral_path.shp").exists()


def test_cli_negative_astronomical_year_is_parsed(tmp_path, capsys):
    out_dir = tmp_path / "ancient"
    rc = cli.main(["-1999", "6", "12", "-o", str(out_dir), "-f", "kml", *FAST_RASTER])

    assert rc == 0
    assert (out_dir / "eclipse.kml").exists()


def test_cli_missing_eclipse_date_returns_error(tmp_path, capsys):
    rc = cli.main(["2024", "4", "9", "-o", str(tmp_path / "x")])

    assert rc == 1
    captured = capsys.readouterr()
    assert "No solar eclipse found" in captured.err
