from pathlib import Path

import pytest

from besselian2shape import BesselianElements, find_by_date, load_all

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"


def test_load_all_parses_every_row():
    elements = load_all(FIXTURE)

    assert len(elements) == 5
    assert all(isinstance(e, BesselianElements) for e in elements)


def test_load_all_casts_field_types():
    elements = load_all(FIXTURE)
    great_american_eclipse = next(e for e in elements if e.year == 2017)

    assert great_american_eclipse.month == 8
    assert great_american_eclipse.day == 21
    assert isinstance(great_american_eclipse.year, int)
    assert isinstance(great_american_eclipse.gamma, float)
    assert isinstance(great_american_eclipse.eclipse_type, str)
    assert isinstance(great_american_eclipse.cat_no, int)


def test_find_by_date_returns_matching_eclipse():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)

    assert e.eclipse_type == "T"
    assert e.gamma == pytest.approx(0.34314)
    assert e.path_width == pytest.approx(197.5)
    assert e.central_duration == "04m28s"
    assert e.x0 == pytest.approx(-0.318244)


def test_find_by_date_raises_when_no_eclipse():
    with pytest.raises(LookupError):
        find_by_date(2024, 4, 9, csv_path=FIXTURE)


def test_find_by_date_2017_eclipse():
    e = find_by_date(2017, 8, 21, csv_path=FIXTURE)

    assert e.lat_ge == "37.0N"
    assert e.lng_ge == "87.7W"
    assert e.saros == 145
