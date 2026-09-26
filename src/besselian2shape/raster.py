"""Raster (coverage-grid) computation of the penumbral visibility boundary.

The region on Earth from which at least a partial eclipse is visible does
not, in general, reduce to a single closed curve split into "northern" and
"southern" limits -- for a shadow whose path passes close to a pole, the
true boundary (as professional eclipse cartography tools compute it) is a
main lens-shaped region plus separate closed loops near the contact
points, stitched from several different curve types (limits, "sun
rise/set" curves, "maximum on horizon" curves).

Rather than classify and stitch those curve types analytically, this
module answers a simpler, robust question directly for a grid of
candidate points: is this point ever, for some instant in the eclipse,
both within the penumbral shadow circle and on Earth's sunlit hemisphere?
Marching squares then extracts whatever boundary curve(s) -- one loop, or
several -- actually result, with no special-casing for topology.
"""

from __future__ import annotations

import math

import numpy as np

from .elements import BesselianElements
from .geometry import ShadowState, envelope_points, sample_times, shadow_state

_B = 0.99664719


def _coverage_mask(
    e: BesselianElements,
    radius_fn,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
    resolution_deg: float,
    step_minutes: float,
    periodic_lon: bool = False,
    seam_lon: float = -180.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lats, lons, covered): a regular grid over the given bounds, and a
    boolean array of shape (len(lats), len(lons)) marking every point ever
    within the shadow circle on Earth's sunlit hemisphere, at any sampled
    instant across the eclipse's tmin..tmax.

    If `periodic_lon`, the grid spans the full globe as a periodic ring in
    longitude (excluding the duplicate column one full turn past `seam_lon`)
    rather than a flat sheet -- needed whenever the covered region wraps
    around the antimeridian, which a flat sheet would otherwise cut
    through, corrupting the traced contour there. `seam_lon` places that
    wraparound seam wherever the caller chooses (default -180, the
    antimeridian itself) -- see `_choose_seam_lon`.
    """
    lats = np.arange(max(lat_min, -90.0), min(lat_max, 90.0) + resolution_deg / 2, resolution_deg)
    lats = lats[lats <= 90.0]
    if periodic_lon:
        lons = np.arange(seam_lon, seam_lon + 360.0, resolution_deg)
        lons = (lons + 180.0) % 360.0 - 180.0
    else:
        lons = np.arange(lon_min, lon_max + resolution_deg / 2, resolution_deg)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    lat_rad = np.radians(lat_grid)
    u = np.arctan(_B * np.tan(lat_rad))
    rho_sin = _B * np.sin(u)
    rho_cos = np.cos(u)
    lon_rad = np.radians(lon_grid)

    covered = np.zeros(lat_grid.shape, dtype=bool)
    for t in sample_times(e.tmin, e.tmax, step_minutes):
        s = shadow_state(e, t)
        radius, _ = radius_fn(s)
        if radius <= 0:
            continue
        h_angle = s.mu + lon_rad
        cos_h = np.cos(h_angle)
        xi = rho_cos * np.sin(h_angle)
        eta = rho_sin * math.cos(s.d) - rho_cos * cos_h * math.sin(s.d)
        zeta = rho_sin * math.sin(s.d) + rho_cos * cos_h * math.cos(s.d)
        dist_sq = (xi - s.x) ** 2 + (eta - s.y) ** 2
        covered |= (zeta > 0) & (dist_sq <= radius * radius)

    return lats, lons, covered


# Margin added to the envelope-formula latitude ceiling (see
# `_envelope_lat_bounds`) before clipping raster coverage to it: the
# envelope formula itself isn't valid right at a start/end tangent point
# or other cusp (see geometry.py), where the true boundary can reach
# slightly further poleward than the smooth envelope alone suggests.
_ENVELOPE_CEILING_MARGIN_DEG = 3.0


def _envelope_lat_bounds(
    e: BesselianElements, radius_fn, step_minutes: float = 1.0
) -> tuple[float, float]:
    """(min_lat, max_lat): the latitude range the standard analytic
    envelope formula (geometry.envelope_points -- the validated, standard
    definition of a shadow circle's northern/southern limit) ever reaches
    for this eclipse's shadow circle. Used as a sanity ceiling on the
    raster scan: coverage claimed well beyond what the envelope itself
    ever supports is treated as suspect rather than trusted outright (see
    `penumbral_coverage_rings`).
    """
    min_lat, max_lat = 90.0, -90.0
    for t in sample_times(e.tmin, e.tmax, step_minutes):
        s = shadow_state(e, t)
        radius, radius_deriv = radius_fn(s)
        if radius <= 0:
            continue
        for p in envelope_points(s, radius, radius_deriv):
            if p is not None:
                min_lat = min(min_lat, p[0])
                max_lat = max(max_lat, p[0])
    return min_lat, max_lat


def _find_bounds(
    e: BesselianElements, radius_fn, coarse_resolution_deg: float, coarse_step_minutes: float
) -> tuple[float, float, float, float, bool] | None:
    """A padded (lat_min, lat_max, lon_min, lon_max, periodic_lon) bounding
    box containing all coverage, found with a fast, coarse global scan.
    None if the shadow never touches Earth's visible surface at all.

    `periodic_lon` is True whenever coverage reaches either edge of the
    -180..180 range, since that means it may wrap around the antimeridian
    -- ambiguous from a coarse scan alone (it can't distinguish "wraps
    around" from "just happens to reach the edge"), so the full-globe
    periodic grid is used to be safe, at the cost of a larger fine pass.
    """
    lats, lons, covered = _coverage_mask(
        e, radius_fn, -90.0, 90.0, -180.0, 180.0, coarse_resolution_deg, coarse_step_minutes
    )
    if not covered.any():
        return None
    rows = np.any(covered, axis=1)
    cols = np.any(covered, axis=0)
    lat_min = lats[rows].min() - 2 * coarse_resolution_deg
    lat_max = lats[rows].max() + 2 * coarse_resolution_deg

    if bool(cols[0]) or bool(cols[-1]):
        return max(lat_min, -90.0), min(lat_max, 90.0), -180.0, 180.0, True

    lon_min = lons[cols].min() - 2 * coarse_resolution_deg
    lon_max = lons[cols].max() + 2 * coarse_resolution_deg
    return max(lat_min, -90.0), min(lat_max, 90.0), lon_min, lon_max, False


# Marching-squares case table: for a cell with corners TL/TR/BR/BL (each
# True if "covered"), which pairs of edges the boundary crosses. Edges are
# named by compass position within the cell (top/right/bottom/left); with
# boolean corner data, a crossed edge's contour point is always its
# midpoint (linear interpolation between a 0 and a 1, at threshold 0.5).
# Cases 5 and 10 are the ambiguous "saddle" configurations; one of the two
# valid diagonal resolutions is picked arbitrarily but consistently.
_CASE_EDGES: dict[int, list[tuple[str, str]]] = {
    0: [],
    1: [("left", "top")],
    2: [("top", "right")],
    3: [("left", "right")],
    4: [("right", "bottom")],
    5: [("left", "top"), ("right", "bottom")],
    6: [("top", "bottom")],
    7: [("left", "bottom")],
    8: [("bottom", "left")],
    9: [("top", "bottom")],
    10: [("top", "right"), ("bottom", "left")],
    11: [("right", "bottom")],
    12: [("left", "right")],
    13: [("top", "right")],
    14: [("left", "top")],
    15: [],
}


def _trace_contours(
    lats: np.ndarray, lons: np.ndarray, covered: np.ndarray, periodic_lon: bool = False
) -> list[list[tuple[float, float]]]:
    """Closed polygon ring(s) bounding the True region of `covered`, via
    marching squares. Each ring is a list of (lat, lon) points, closed
    (first point repeated at the end).

    If `periodic_lon`, `lons` is treated as wrapping around (column n_lon-1
    adjacent to column 0) rather than as a flat sheet with a hard edge at
    its last column -- needed for a full-globe grid, so a region that
    happens to cross the antimeridian isn't cut in two there.
    """
    n_lat, n_lon = covered.shape

    def edge_point(edge: tuple) -> tuple[float, float]:
        kind, a, b = edge
        if kind == "H":
            # Average the two columns' longitudes, unwrapped relative to
            # each other first: every longitude here is stored wrapped
            # into [-180, 180], and the one pair of geographically
            # adjacent columns whose true (unwrapped) longitudes straddle
            # +/-180 would otherwise average to the antipodal point.
            lon_a = float(lons[b])
            lon_b = float(lons[(b + 1) % n_lon])
            if lon_b - lon_a > 180:
                lon_b -= 360
            elif lon_a - lon_b > 180:
                lon_b += 360
            lon_mid = (lon_a + lon_b) / 2
            return float(lats[a]), (lon_mid + 180.0) % 360.0 - 180.0
        return float((lats[a] + lats[a + 1]) / 2), float(lons[b])

    adjacency: dict[tuple, list[tuple]] = {}

    def connect(e1: tuple, e2: tuple) -> None:
        adjacency.setdefault(e1, []).append(e2)
        adjacency.setdefault(e2, []).append(e1)

    j_span = range(n_lon) if periodic_lon else range(n_lon - 1)
    for i in range(n_lat - 1):
        for j in j_span:
            j_next = (j + 1) % n_lon if periodic_lon else j + 1
            bl, br, tl, tr = (
                covered[i, j],
                covered[i, j_next],
                covered[i + 1, j],
                covered[i + 1, j_next],
            )
            case = tl * 1 + tr * 2 + br * 4 + bl * 8
            pairs = _CASE_EDGES[case]
            if not pairs:
                continue
            edges = {
                "top": ("H", i + 1, j),
                "bottom": ("H", i, j),
                "left": ("V", i, j),
                "right": ("V", i, j_next),
            }
            for a, b in pairs:
                connect(edges[a], edges[b])

    rings: list[list[tuple[float, float]]] = []
    visited: set[tuple] = set()
    for start in list(adjacency):
        if start in visited:
            continue
        ring_edges = [start]
        visited.add(start)
        prev, current = None, start
        while True:
            neighbors = [n for n in adjacency[current] if n != prev]
            nxt = neighbors[0] if neighbors else None
            if nxt is None or nxt == start:
                break
            ring_edges.append(nxt)
            visited.add(nxt)
            prev, current = current, nxt

        if len(ring_edges) < 3:
            continue
        ring = [edge_point(edge) for edge in ring_edges]
        ring.append(ring[0])
        rings.append(ring)

    return rings


def _crossings(pts: list[tuple[float, float]]) -> list[tuple[int, int]]:
    """(index, direction) for every antimeridian crossing in a point
    sequence treated as a closed loop -- direction +1 for a "positive to
    negative" jump (continuing to increasing longitude, wrapping past
    +180 to -180), -1 for the reverse.
    """
    n = len(pts)
    result = []
    for i in range(n):
        lon1, lon2 = pts[i][1], pts[(i + 1) % n][1]
        if lon2 - lon1 > 180:
            result.append((i, -1))
        elif lon1 - lon2 > 180:
            result.append((i, 1))
    return result


def _close_directly(arc: list[tuple[float, float]]) -> list[tuple[float, float]]:
    closed = list(arc)
    closed.append(closed[0])
    return closed


def _close_through_pole(chain: list[tuple[float, float]]) -> list[tuple[float, float]]:
    pole_lat = 90.0 if sum(p[0] for p in chain) >= 0 else -90.0
    closed = chain + [(pole_lat, chain[-1][1]), (pole_lat, chain[0][1])]
    closed.append(closed[0])
    return closed


def _split_at_antimeridian(
    ring: list[tuple[float, float]],
) -> list[list[tuple[float, float]]]:
    """Split a closed (lat, lon) ring into one or more pieces that never
    cross longitude +/-180 without being cut there -- KML and shapefile
    polygon rings are just flat sequences of vertices with no built-in
    understanding that longitude wraps around, so a ring that crosses the
    antimeridian uncut is commonly mis-rendered (a wraparound band across
    the whole map) by viewers that just connect consecutive vertices.

    A ring's crossings can be more than the minimum needed for its net
    effect: e.g. a boundary that encircles a pole (an odd net number of
    crossings of any fixed meridian) can *also* have a small excursion
    that dips across the antimeridian and back nearby, adding one
    "cancelling" pair of crossings (opposite directions, adjacent in
    crossing order) on top of the net crossing that really needs pole
    routing. Such adjacent opposite-direction pairs are repeatedly
    extracted as their own standalone, directly-closed pieces (the same
    as a simple 2-crossing dip-and-return) until at most one crossing is
    left; if one remains, that final loop is closed by routing through
    whichever pole it encircles. A pattern that doesn't resolve this way
    (e.g. adjacent crossings that never alternate direction) is left
    unsplit as a safe fallback.
    """
    pts = ring[:-1]
    if not _crossings(pts):
        return [ring]

    pieces: list[list[tuple[float, float]]] = []
    remaining = pts
    for _ in range(len(pts)):  # generous bound; each iteration removes >=1 crossing
        crossings = _crossings(remaining)
        if len(crossings) <= 1:
            break

        extracted = None
        for a in range(len(crossings)):
            b = (a + 1) % len(crossings)
            i1, d1 = crossings[a]
            i2, d2 = crossings[b]
            if d1 == d2:
                continue
            m = len(remaining)
            if i2 > i1:
                arc, rest = remaining[i1 + 1 : i2 + 1], remaining[: i1 + 1] + remaining[i2 + 1 :]
            else:
                arc = remaining[i1 + 1 :] + remaining[: i2 + 1]
                rest = remaining[i2 + 1 : i1 + 1]
            extracted = (arc, rest)
            break

        if extracted is None:
            return [ring]  # crossings never alternate direction; bail out safely
        arc, remaining = extracted
        pieces.append(_close_directly(arc))

    crossings = _crossings(remaining)
    if not crossings:
        pieces.append(_close_directly(remaining))
    else:
        i, _direction = crossings[0]
        chain = remaining[i + 1 :] + remaining[: i + 1]
        pieces.append(_close_through_pole(chain))

    return pieces


def _choose_seam_lon(e: BesselianElements) -> float:
    """A longitude, antipodal to the eclipse's greatest-eclipse point, to
    place the periodic grid's wraparound seam at. Real coverage is very
    unlikely to reach there, so the traced boundary rarely needs to cross
    it at all -- unlike the fixed antimeridian, which a shadow's own path
    can legitimately pass right through, producing extra crossings
    wherever the boundary happens to wiggle near it.
    """
    antipodal = e.lng_dd_ge + 180.0
    return (antipodal + 180.0) % 360.0 - 180.0


def _point_in_ring(pt: tuple[float, float], ring: list[tuple[float, float]]) -> bool:
    """Even-odd (ray casting) point-in-polygon test."""
    lat, lon = pt
    inside = False
    n = len(ring) - 1
    for i in range(n):
        lat1, lon1 = ring[i]
        lat2, lon2 = ring[i + 1]
        if (lon1 > lon) != (lon2 > lon):
            lat_intersect = lat1 + (lon - lon1) / (lon2 - lon1) * (lat2 - lat1)
            if lat < lat_intersect:
                inside = not inside
    return inside


def _ring_area(ring: list[tuple[float, float]]) -> float:
    total = 0.0
    for (lat1, lon1), (lat2, lon2) in zip(ring, ring[1:]):
        total += lon1 * lat2 - lon2 * lat1
    return abs(total) / 2.0


def _group_rings_with_holes(
    rings: list[list[tuple[float, float]]],
) -> list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]]:
    """Group rings into (outer_ring, [hole_ring, ...]) polygons.

    Marching squares on a coverage grid traces both the outer edge of a
    covered region and the edge of any uncovered gap within it (a "hole")
    as separate, equally-valid closed rings -- nothing about the trace
    itself marks which is which. A ring is classified as a hole of
    whichever other ring most tightly contains it (tested via a single
    point on it), rather than always being filled as if it were its own
    separate covered region.
    """
    n = len(rings)
    sample_points = [r[0] for r in rings]
    containers: list[list[int]] = [
        [j for j in range(n) if j != i and _point_in_ring(sample_points[i], rings[j])]
        for i in range(n)
    ]
    parent = [
        min(c, key=lambda j: _ring_area(rings[j])) if c else None for c in containers
    ]

    groups: dict[int, list[int]] = {i: [] for i in range(n) if parent[i] is None}
    for i in range(n):
        if parent[i] is not None:
            groups.setdefault(parent[i], []).append(i)

    return [(rings[outer], [rings[h] for h in holes]) for outer, holes in groups.items()]


def penumbral_coverage_rings(
    e: BesselianElements,
    radius_fn,
    resolution_deg: float = 0.15,
    step_minutes: float = 1.0,
    coarse_resolution_deg: float = 1.0,
    coarse_step_minutes: float = 5.0,
) -> list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]]:
    """Every closed boundary (outer_ring, [hole_ring, ...]) polygon of the
    region ever within the shadow circle (per `radius_fn`) on Earth's
    sunlit hemisphere, across the eclipse's tmin..tmax. Usually a single
    polygon with no holes; can be several, and/or have holes, for a
    shadow path that passes close to a pole. Empty if the shadow never
    touches Earth's visible surface at all.
    """
    bounds = _find_bounds(e, radius_fn, coarse_resolution_deg, coarse_step_minutes)
    if bounds is None:
        return []
    lat_min, lat_max, lon_min, lon_max, periodic = bounds
    seam_lon = _choose_seam_lon(e) if periodic else -180.0
    lats, lons, covered = _coverage_mask(
        e,
        radius_fn,
        lat_min,
        lat_max,
        lon_min,
        lon_max,
        resolution_deg,
        step_minutes,
        periodic,
        seam_lon,
    )

    envelope_min_lat, envelope_max_lat = _envelope_lat_bounds(e, radius_fn)
    covered &= (lats[:, None] <= envelope_max_lat + _ENVELOPE_CEILING_MARGIN_DEG) & (
        lats[:, None] >= envelope_min_lat - _ENVELOPE_CEILING_MARGIN_DEG
    )

    rings = _trace_contours(lats, lons, covered, periodic)
    split_rings = [piece for ring in rings for piece in _split_at_antimeridian(ring)]
    return _group_rings_with_holes(split_rings)
