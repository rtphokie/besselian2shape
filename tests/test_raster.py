import math
from pathlib import Path

import numpy as np
import pytest

from besselian2shape.elements import find_by_date
from besselian2shape.raster import _split_at_antimeridian, _trace_contours, penumbral_coverage_rings

FIXTURE = Path(__file__).parent / "fixtures" / "sample_besselian.csv"


def _dist_km(p1, p2):
    lat1, lon1 = p1
    lat2, lon2 = p2
    dlon = (lon2 - lon1 + 180) % 360 - 180
    return math.hypot((lat2 - lat1) * 111, dlon * 111 * math.cos(math.radians((lat1 + lat2) / 2)))


def test_trace_contours_simple_square():
    # A 4x4 grid with a 2x2 True block in the middle should trace one
    # closed ring around it.
    lats = np.array([0.0, 1.0, 2.0, 3.0])
    lons = np.array([0.0, 1.0, 2.0, 3.0])
    covered = np.array(
        [
            [False, False, False, False],
            [False, True, True, False],
            [False, True, True, False],
            [False, False, False, False],
        ]
    )
    rings = _trace_contours(lats, lons, covered)
    assert len(rings) == 1
    ring = rings[0]
    assert ring[0] == ring[-1]
    ring_lats = [p[0] for p in ring]
    ring_lons = [p[1] for p in ring]
    assert min(ring_lats) == pytest.approx(0.5)
    assert max(ring_lats) == pytest.approx(2.5)
    assert min(ring_lons) == pytest.approx(0.5)
    assert max(ring_lons) == pytest.approx(2.5)


def test_trace_contours_periodic_connects_across_seam():
    # A True region spanning the last and first longitude columns (270 and
    # 0, adjacent under wraparound) should trace as one connected ring,
    # not be cut into pieces at that seam.
    lats = np.array([0.0, 1.0, 2.0])
    lons = np.array([0.0, 90.0, 180.0, 270.0])  # periodic: 270 adjacent to 0
    covered = np.array(
        [
            [False, False, False, False],
            [True, False, False, True],
            [False, False, False, False],
        ]
    )
    rings = _trace_contours(lats, lons, covered, periodic_lon=True)
    assert len(rings) == 1
    ring = rings[0]
    assert ring[0] == ring[-1]
    assert len(ring) > 4  # more than a single isolated cell's own small ring

    # Without wraparound, the same grid must split into two separate rings
    # (one per isolated True cell), confirming periodic mode is doing
    # something different -- not that it happens to produce one ring by
    # coincidence of this particular grid.
    rings_flat = _trace_contours(lats, lons, covered, periodic_lon=False)
    assert len(rings_flat) == 2


def test_penumbral_coverage_rings_total_eclipse_matches_reference_bounds():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    polygons = penumbral_coverage_rings(
        e, lambda s: (s.l1, s.l1p), resolution_deg=0.3, step_minutes=2.0
    )
    assert len(polygons) == 1
    outer, holes = polygons[0]
    assert outer[0] == outer[-1]
    assert holes == []

    # Reference (Xavier Jubier KMZ) Penumbra Northern/Southern Limit
    # endpoints for this eclipse; our raster boundary should pass close by.
    reference_points = [
        (33.53159, -177.05044),
        (89.36876, 72.32444),
        (-38.74931, -151.56625),
        (16.78002, -27.89899),
    ]
    for ref in reference_points:
        nearest = min(outer, key=lambda p: _dist_km(p, ref))
        assert _dist_km(nearest, ref) < 100


def test_penumbral_coverage_rings_partial_eclipse_no_freeze():
    e = find_by_date(2025, 3, 29, csv_path=FIXTURE)
    polygons = penumbral_coverage_rings(
        e, lambda s: (s.l1, s.l1p), resolution_deg=0.5, step_minutes=3.0
    )
    assert len(polygons) >= 1
    for outer, holes in polygons:
        assert outer[0] == outer[-1]
        assert len(outer) > 3
        for hole in holes:
            assert hole[0] == hole[-1]


def test_penumbral_coverage_rings_empty_for_out_of_range_radius():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    polygons = penumbral_coverage_rings(
        e, lambda s: (0.0, 0.0), resolution_deg=1.0, step_minutes=5.0
    )
    assert polygons == []


def _max_lon_jump(ring):
    n = len(ring) - 1
    return max(abs(ring[j + 1][1] - ring[j][1]) for j in range(n))


def test_split_at_antimeridian_no_crossing_is_unchanged():
    ring = [(0.0, 10.0), (1.0, 20.0), (0.0, 30.0), (0.0, 10.0)]
    assert _split_at_antimeridian(ring) == [ring]


def test_split_at_antimeridian_pole_enclosing_ring_has_no_real_jump():
    # A ring circling the north pole at ~80N, crossing the seam once.
    ring = [
        (80.0, -179.0),
        (80.0, -90.0),
        (80.0, 0.0),
        (80.0, 90.0),
        (80.0, 179.0),
        (80.0, -179.0),
    ]
    pieces = _split_at_antimeridian(ring)
    assert len(pieces) == 1
    piece = pieces[0]
    assert piece[0] == piece[-1]
    # Every remaining longitude jump must be confined to the pole itself
    # (lat == 90), which is geographically a single point regardless of
    # its stated longitude -- not a real jump across the map.
    n = len(piece) - 1
    for j in range(n):
        jump = abs(piece[j + 1][1] - piece[j][1])
        if jump > 170:
            assert piece[j][0] == 90.0 and piece[j + 1][0] == 90.0


def test_split_at_antimeridian_simple_excursion_splits_into_two_pieces():
    # A ring that dips across the antimeridian once and comes back,
    # crossing it exactly twice.
    ring = [
        (0.0, 170.0),
        (5.0, 179.0),
        (5.0, -179.0),  # crossing 1 (into the "east of seam" arc)
        (0.0, -170.0),
        (-5.0, -179.0),
        (-5.0, 179.0),  # crossing 2 (back to the "west of seam" arc)
        (0.0, 170.0),
    ]
    pieces = _split_at_antimeridian(ring)
    assert len(pieces) == 2
    for piece in pieces:
        assert piece[0] == piece[-1]
        assert _max_lon_jump(piece) < 20  # no spurious near-360 jump


def test_penumbral_coverage_rings_no_uncut_antimeridian_jump():
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    polygons = penumbral_coverage_rings(
        e, lambda s: (s.l1, s.l1p), resolution_deg=0.3, step_minutes=2.0
    )
    for outer, holes in polygons:
        for ring in [outer, *holes]:
            n = len(ring) - 1
            for j in range(n):
                jump = abs(ring[j + 1][1] - ring[j][1])
                if jump > 170:
                    assert ring[j][0] == 90.0 and ring[j + 1][0] == 90.0


def test_penumbral_coverage_rings_annular_eclipse_has_a_hole():
    # This eclipse's penumbral boundary has a genuine uncovered gap near
    # the antimeridian, plus a south-polar-cap hole from the
    # envelope-ceiling clip (see `_envelope_lat_bounds`). At overly coarse
    # resolution these two holes can fragment into extra disconnected
    # pieces, so this test uses a resolution known to be stable.
    e = find_by_date(2024, 10, 2, csv_path=FIXTURE)
    polygons = penumbral_coverage_rings(
        e, lambda s: (s.l1, s.l1p), resolution_deg=0.35, step_minutes=2.0
    )
    assert len(polygons) == 1
    outer, holes = polygons[0]
    assert outer[0] == outer[-1]
    assert len(holes) == 2
    for hole in holes:
        assert hole[0] == hole[-1]
        assert len(hole) > 3


def test_group_rings_with_holes_nests_inner_ring():
    from besselian2shape.raster import _group_rings_with_holes

    outer = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)]
    hole = [(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0), (4.0, 4.0)]
    island = [(20.0, 20.0), (20.0, 22.0), (22.0, 22.0), (22.0, 20.0), (20.0, 20.0)]

    groups = _group_rings_with_holes([outer, hole, island])
    assert len(groups) == 2
    by_first_point = {g[0][0]: g for g in groups}
    assert by_first_point[outer[0]][1] == [hole]
    assert by_first_point[island[0]][1] == []


def test_penumbral_coverage_capped_at_envelope_ceiling():
    # Coverage is capped at the envelope formula's own max latitude reach
    # (a standard, validated definition) as a sanity ceiling, so a narrow
    # high-latitude band cannot spuriously extend the raster boundary all
    # the way to the pole.
    from besselian2shape.raster import _envelope_lat_bounds, _ENVELOPE_CEILING_MARGIN_DEG

    e = find_by_date(2017, 8, 21, csv_path=FIXTURE)
    radius_fn = lambda s: (s.l1, s.l1p)
    _, envelope_max_lat = _envelope_lat_bounds(e, radius_fn)

    polygons = penumbral_coverage_rings(e, radius_fn, resolution_deg=0.3, step_minutes=2.0)
    assert len(polygons) == 1
    outer, holes = polygons[0]
    # The outer ring's max is always ~90 (the synthetic pole-closing
    # point), but any polar-cap hole excising over-extension must start
    # at or below the envelope ceiling plus its margin.
    for hole in holes:
        hole_min_lat = min(p[0] for p in hole)
        assert hole_min_lat <= envelope_max_lat + _ENVELOPE_CEILING_MARGIN_DEG + 1.0


def test_penumbral_coverage_2024_still_reaches_near_pole():
    # A genuine near-pole reach must not be clipped away by the
    # envelope-ceiling guard.
    e = find_by_date(2024, 4, 8, csv_path=FIXTURE)
    radius_fn = lambda s: (s.l1, s.l1p)
    polygons = penumbral_coverage_rings(e, radius_fn, resolution_deg=0.3, step_minutes=2.0)
    assert len(polygons) == 1
    outer, holes = polygons[0]
    assert holes == []
    assert max(p[0] for p in outer) >= 85.0
