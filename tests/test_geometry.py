import math
import random
from pathlib import Path

import pytest

from besselian2shape.elements import find_by_date
from besselian2shape.geometry import (
    _B,
    central_line,
    central_line_point,
    invert,
    is_central,
    limits_to_ring,
    umbral_limits,
)

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"


def _forward(lat, lon, d, mu):
    """Reference forward transform (geographic -> fundamental plane),
    independent of the code under test, used to validate invert()."""
    u = math.atan(_B * math.tan(math.radians(lat)))
    rho_sin = _B * math.sin(u)
    rho_cos = math.cos(u)
    h_angle = math.radians(math.degrees(mu) + lon)
    xi = rho_cos * math.sin(h_angle)
    eta = rho_sin * math.cos(d) - rho_cos * math.cos(h_angle) * math.sin(d)
    zeta = rho_sin * math.sin(d) + rho_cos * math.cos(h_angle) * math.cos(d)
    return xi, eta, zeta


def test_invert_round_trips_forward_transform():
    rng = random.Random(42)
    tested = 0
    for _ in range(500):
        lat = rng.uniform(-80, 80)
        lon = rng.uniform(-179, 179)
        d = math.radians(rng.uniform(-23, 23))
        mu = math.radians(rng.uniform(0, 360))
        xi, eta, zeta = _forward(lat, lon, d, mu)
        if zeta <= 0:
            continue
        tested += 1
        result = invert(xi, eta, d, mu)
        assert result is not None
        lat2, lon2 = result
        assert lat2 == pytest.approx(lat, abs=1e-3)
        dlon = ((lon2 - lon) + 180) % 360 - 180
        assert dlon == pytest.approx(0.0, abs=1e-3)
    assert tested > 100


def test_invert_returns_none_for_far_side():
    # (x, y) with x^2 + y^2 > 1 cannot be on the Earth's disk at all.
    assert invert(2.0, 2.0, 0.0, 0.0) is None


def test_total_eclipse_2024_is_central_with_all_layers():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    assert is_central(e)

    line = central_line(e, step_minutes=2)
    assert len(line) > 50
    # Path crosses the Pacific to the North Atlantic, roughly equator to ~49N.
    lats = [p[0] for p in line]
    lons = [p[1] for p in line]
    assert min(lats) < 0 < max(lats) < 55
    assert min(lons) < -150
    assert max(lons) > -40

    u_ring = limits_to_ring(umbral_limits(e, step_minutes=2))
    assert len(u_ring) > 50
    assert u_ring[0] == u_ring[-1]


def test_total_eclipse_central_line_near_known_greatest_eclipse():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    h, m, s = (int(x) for x in e.td_ge.split(":"))
    t_ge = (h + m / 60 + s / 3600) - e.t0

    point = central_line_point(e, t_ge)
    assert point is not None
    lat, lon = point
    # Greatest eclipse and the central line at that same clock time are close
    # but not identical points; a loose tolerance just guards against gross
    # sign/formula errors.
    assert lat == pytest.approx(e.lat_dd_ge, abs=1.0)
    assert lon == pytest.approx(e.lng_dd_ge, abs=1.0)


def test_partial_eclipse_2025_is_not_central():
    e = find_by_date(2025, 3, 29, csv_path=FIXTURE)
    assert not is_central(e)
    assert central_line(e) == []
    assert limits_to_ring(umbral_limits(e)) == []
