"""Parsing of NASA's Besselian elements CSV and lookup by date."""

from __future__ import annotations

import csv
from dataclasses import dataclass, fields
from pathlib import Path

from .cache import download_besselian_csv

_INT_FIELDS = {
    "year", "month", "day", "luna_num", "saros", "cat_no", "canon_plate",
    "etype", "PNS", "UNS", "NCN", "nSer", "nSeq", "nJLE",
}

_FLOAT_FIELDS = {
    "dt", "gamma", "magnitude", "lat_dd_ge", "lng_dd_ge", "sun_alt", "sun_azm",
    "path_width", "duration_secs", "julian_date", "t0",
    "x0", "x1", "x2", "x3", "y0", "y1", "y2", "y3", "d0", "d1", "d2",
    "mu0", "mu1", "mu2", "l10", "l11", "l12", "l20", "l21", "l22",
    "tan_f1", "tan_f2", "tmin", "tmax",
}


@dataclass(frozen=True, slots=True)
class BesselianElements:
    """Besselian elements and metadata for a single solar eclipse.

    Field names and units follow NASA's Five Millennium Canon of Solar
    Eclipses mysqldump export: x/y are the coordinates of the shadow axis
    on the fundamental plane, d/mu describe the axis's declination and
    hour angle, l1/l2 are the penumbral/umbral radii, and tan_f1/tan_f2
    are the penumbral/umbral cone angles. All are polynomials in time (in
    hours) measured from t0 (TDT), valid over [tmin, tmax].
    """

    year: int
    month: int
    day: int
    td_ge: str
    dt: float
    luna_num: int
    saros: int
    eclipse_type: str
    gamma: float
    magnitude: float
    lat_ge: str
    lng_ge: str
    lat_dd_ge: float
    lng_dd_ge: float
    sun_alt: float
    sun_azm: float
    path_width: float
    central_duration: str
    duration_secs: float
    cat_no: int
    canon_plate: int
    julian_date: float
    t0: float
    x0: float
    x1: float
    x2: float
    x3: float
    y0: float
    y1: float
    y2: float
    y3: float
    d0: float
    d1: float
    d2: float
    mu0: float
    mu1: float
    mu2: float
    l10: float
    l11: float
    l12: float
    l20: float
    l21: float
    l22: float
    tan_f1: float
    tan_f2: float
    tmin: float
    tmax: float
    etype: int
    PNS: int
    UNS: int
    NCN: int
    nSer: int
    nSeq: int
    nJLE: int

    @classmethod
    def from_row(cls, row: dict[str, str]) -> "BesselianElements":
        kwargs = {}
        for f in fields(cls):
            raw = row[f.name]
            if f.name in _INT_FIELDS:
                kwargs[f.name] = int(float(raw))
            elif f.name in _FLOAT_FIELDS:
                kwargs[f.name] = float(raw)
            else:
                kwargs[f.name] = raw
        return cls(**kwargs)


def load_all(csv_path: Path | None = None) -> list[BesselianElements]:
    """Load every eclipse's Besselian elements from the (cached) CSV.

    Downloads and caches the CSV first if `csv_path` is not given and no
    cached copy exists yet.
    """
    path = csv_path if csv_path is not None else download_besselian_csv()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [BesselianElements.from_row(row) for row in reader]


def find_by_date(
    year: int, month: int, day: int, csv_path: Path | None = None
) -> BesselianElements:
    """Return the Besselian elements for the solar eclipse on the given
    (proleptic Gregorian, astronomical-numbered) calendar date.

    Use astronomical year numbering for BCE dates (e.g. 1 BCE is year 0,
    2 BCE is year -1). Raises LookupError if no eclipse occurred on that
    date in the dataset.
    """
    matches = [
        e for e in load_all(csv_path)
        if e.year == year and e.month == month and e.day == day
    ]
    if not matches:
        raise LookupError(f"No solar eclipse found for {year:04d}-{month:02d}-{day:02d}")
    if len(matches) > 1:
        raise LookupError(
            f"Multiple eclipses found for {year:04d}-{month:02d}-{day:02d}: {matches}"
        )
    return matches[0]
