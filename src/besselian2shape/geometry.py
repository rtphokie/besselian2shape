"""Solar eclipse path geometry from Besselian elements.

Implements the standard "Besselian elements" method (Meeus, *Astronomical
Algorithms*, ch. 54; Espenak & Meeus, *Five Millennium Canon of Solar
Eclipses*) for projecting the Moon's shadow onto the Earth's oblate surface:

- The Moon's shadow axis pierces a plane through the Earth's center,
  perpendicular to the axis (the "fundamental plane"). Besselian elements
  x(t), y(t) give the axis's position on that plane, in units of the
  Earth's equatorial radius, as polynomials in t = hours from t0 (TDT).
- l1(t)/l2(t) are the radii, in the same units, of the penumbral/umbral
  (or antumbral, for an annular eclipse) shadow circles on that plane.
- d(t)/mu(t) are the declination and Greenwich hour angle of the point
  where the shadow axis crosses the fundamental plane.

Given (x, y, d, mu) at an instant, the geographic point(s) on Earth's
ellipsoid where the fundamental-plane coordinates equal (x, y) are found
by inverting the standard geocentric-parallax reduction (the same
reduction used to compute topocentric parallax elsewhere in positional
astronomy), solved here in closed form via a quadratic rather than
iteratively.

This module computes the central line and umbral/antumbral (path of
totality/annularity) limits. The penumbral (partial-visibility) boundary
is computed separately, by direct coverage sampling -- see `raster.py`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .elements import BesselianElements

# b = sqrt(1 - e^2) for the reference ellipsoid (f = 1/298.257) used by
# NASA's Besselian element datasets -- the ratio of polar to equatorial
# radius.
_B = 0.99664719


@dataclass(frozen=True, slots=True)
class ShadowState:
    """Fundamental-plane shadow position/geometry at one instant."""

    t: float
    x: float
    y: float
    d: float  # radians
    mu: float  # radians
    l1: float
    l2: float
    xp: float  # dx/dt
    yp: float  # dy/dt
    l1p: float  # dl1/dt
    l2p: float  # dl2/dt


def _poly(c0: float, c1: float, c2: float, c3: float, t: float) -> float:
    return c0 + t * (c1 + t * (c2 + t * c3))


def _poly_deriv(c1: float, c2: float, c3: float, t: float) -> float:
    return c1 + t * (2 * c2 + t * 3 * c3)


def shadow_state(e: BesselianElements, t: float) -> ShadowState:
    """Evaluate the Besselian polynomials at t hours from t0 (TDT)."""
    x = _poly(e.x0, e.x1, e.x2, e.x3, t)
    y = _poly(e.y0, e.y1, e.y2, e.y3, t)
    d = math.radians(_poly(e.d0, e.d1, e.d2, 0.0, t))
    mu = math.radians(_poly(e.mu0, e.mu1, e.mu2, 0.0, t))
    l1 = _poly(e.l10, e.l11, e.l12, 0.0, t)
    l2 = _poly(e.l20, e.l21, e.l22, 0.0, t)
    xp = _poly_deriv(e.x1, e.x2, e.x3, t)
    yp = _poly_deriv(e.y1, e.y2, e.y3, t)
    l1p = _poly_deriv(e.l11, e.l12, 0.0, t)
    l2p = _poly_deriv(e.l21, e.l22, 0.0, t)
    return ShadowState(t, x, y, d, mu, l1, l2, xp, yp, l1p, l2p)


def invert(x: float, y: float, d: float, mu: float) -> tuple[float, float] | None:
    """Invert fundamental-plane coordinates (x, y) to geographic (lat, lon)
    in degrees, given the shadow axis's declination `d` and Greenwich hour
    angle `mu` (both radians), at sea level (h=0).

    Returns None if (x, y) does not correspond to any point on the visible
    (sunward) hemisphere of the Earth's ellipsoid at this instant.

    Derivation: writing A = rho*cos(phi') and B = rho*sin(phi') for the
    geocentric parallax factors (so A = cos(u), B = b*sin(u) with u the
    reduced latitude), the forward relations
        x = A sin(H),  y = B cos(d) - A cos(H) sin(d)
    combined with the ellipse constraint A^2 + (B/b)^2 = 1 (H = mu + lon)
    reduce to a quadratic in B, solved here directly rather than via
    iteration.
    """
    sin_d = math.sin(d)
    cos_d = math.cos(d)

    a_coef = cos_d * cos_d + (sin_d * sin_d) / (_B * _B)
    b_coef = -2.0 * y * cos_d
    c_coef = x * x * sin_d * sin_d + y * y - sin_d * sin_d

    disc = b_coef * b_coef - 4 * a_coef * c_coef
    if disc < 0:
        return None
    sqrt_disc = math.sqrt(disc)

    best: tuple[float, float, float] | None = None  # (zeta, lat, lon)
    for b_val in ((-b_coef + sqrt_disc) / (2 * a_coef), (-b_coef - sqrt_disc) / (2 * a_coef)):
        a_sq = 1.0 - (b_val / _B) ** 2
        if a_sq < 0:
            continue
        a_val = math.sqrt(a_sq)  # cos(u) >= 0 branch (u in (-90, 90))

        if a_val < 1e-12:
            sin_H, cos_H = 0.0, 1.0
        else:
            sin_H = max(-1.0, min(1.0, x / a_val))
            if abs(sin_d) < 1e-12:
                cos_H = math.sqrt(max(0.0, 1.0 - sin_H * sin_H))
            else:
                cos_H = (b_val * cos_d - y) / (a_val * sin_d)

        h_angle = math.atan2(sin_H, cos_H)
        zeta = b_val * sin_d + a_val * math.cos(h_angle) * cos_d
        if zeta < 0:
            continue

        u = math.atan2(b_val / _B, a_val)
        phi = math.atan2(math.tan(u), _B)
        lon = math.degrees(h_angle) - math.degrees(mu)
        lon = (lon + 180.0) % 360.0 - 180.0
        lat = math.degrees(phi)

        if best is None or zeta > best[0]:
            best = (zeta, lat, lon)

    if best is None:
        return None
    return best[1], best[2]


def central_line_point(e: BesselianElements, t: float) -> tuple[float, float] | None:
    """Geographic (lat, lon) of the shadow axis (central line) at time t,
    or None if the axis misses the Earth at this instant.
    """
    s = shadow_state(e, t)
    return invert(s.x, s.y, s.d, s.mu)


def envelope_points(
    s: ShadowState, radius: float, radius_deriv: float
) -> tuple[tuple[float, float] | None, tuple[float, float] | None]:
    """The two points of the envelope of the moving, radius-varying shadow
    circle at this instant -- i.e. the points where the circle at time t is
    tangent to its immediate neighbors in the swept family, which is the
    standard definition of the northern/southern limit curves (Meeus,
    *Astronomical Algorithms*, ch. 54).

    For a circle of center (x, y), radius r, moving with velocity (xp, yp),
    a point at angle theta (measured from the circle's center) lies on the
    envelope iff cos(theta - psi) = -r'/speed, where psi is the direction
    of travel and speed = |(xp, yp)|. This has two solutions (the north and
    south branches); for a constant radius this reduces to the simpler
    "perpendicular to the direction of travel" rule.

    Returns (point_a, point_b) in no particular north/south order, each
    None if that envelope point isn't on Earth's visible surface, or if
    the shadow circle isn't moving.
    """
    speed = math.hypot(s.xp, s.yp)
    if speed == 0:
        return None, None
    psi = math.atan2(s.yp, s.xp)
    cos_offset = max(-1.0, min(1.0, -radius_deriv / speed))
    delta = math.acos(cos_offset)

    points = []
    for theta in (psi + delta, psi - delta):
        points.append(
            invert(s.x + radius * math.cos(theta), s.y + radius * math.sin(theta), s.d, s.mu)
        )
    return points[0], points[1]


def _lat_at_angle(s: ShadowState, radius: float, theta: float) -> float | None:
    p = invert(s.x + radius * math.cos(theta), s.y + radius * math.sin(theta), s.d, s.mu)
    return p[0] if p is not None else None


def _bisect_visibility_edge(
    s: ShadowState, radius: float, visible_theta: float, invisible_theta: float, iters: int = 30
) -> float:
    """Refine the angle, between a known-visible and a known-invisible angle
    on the circle, where it crosses Earth's visible-disk edge."""
    a, b = visible_theta, invisible_theta
    for _ in range(iters):
        mid = (a + b) / 2
        if _lat_at_angle(s, radius, mid) is not None:
            a = mid
        else:
            b = mid
    return a


def _find_visible_arc(s: ShadowState, radius: float) -> tuple[float, float] | None:
    """The single contiguous range of angles theta (theta_hi > theta_lo,
    possibly > 2*pi if the arc wraps through 0) over which a point on the
    shadow circle is on Earth's visible surface. None if no point is
    visible; (0, 2*pi) if the whole circle is.

    Anchored at the circle's point farthest from the fundamental plane's
    origin (the weakest point for visibility) and its opposite, closest
    point (the strongest, so visible whenever any point on the circle is),
    then bisects outward from the weak point in each direction to find the
    gap's two edges precisely. Assumes a single contiguous invisible gap
    around the weak point.
    """
    weak_theta = math.atan2(s.y, s.x)
    if _lat_at_angle(s, radius, weak_theta) is not None:
        return 0.0, 2 * math.pi

    strong_theta = weak_theta + math.pi
    if _lat_at_angle(s, radius, strong_theta) is None:
        return None

    theta_lo = _bisect_visibility_edge(s, radius, weak_theta + math.pi, weak_theta)
    theta_hi = _bisect_visibility_edge(s, radius, weak_theta - math.pi, weak_theta)
    if theta_hi < theta_lo:
        theta_hi += 2 * math.pi
    return theta_lo, theta_hi


def sample_times(tmin: float, tmax: float, step_minutes: float = 0.5) -> list[float]:
    """Evenly spaced sample times (hours from t0) across [tmin, tmax]."""
    step = step_minutes / 60.0
    n = max(2, int(round((tmax - tmin) / step)) + 1)
    return [tmin + i * (tmax - tmin) / (n - 1) for i in range(n)]


def is_central(e: BesselianElements) -> bool:
    """True if this eclipse produces a central (total/annular/hybrid) path,
    as opposed to a partial-only eclipse.
    """
    return not e.eclipse_type.startswith("P")


def central_line(e: BesselianElements, step_minutes: float = 0.5) -> list[tuple[float, float]]:
    """Central line of a total/annular/hybrid eclipse, as (lat, lon) points
    ordered west-to-east (increasing t). Empty if the eclipse is partial-only
    or the axis never reaches the Earth.
    """
    if not is_central(e):
        return []
    points = []
    for t in sample_times(e.tmin, e.tmax, step_minutes):
        p = central_line_point(e, t)
        if p is not None:
            points.append(p)
    return points


@dataclass(frozen=True, slots=True)
class LimitTrack:
    """Parallel north/south limit curves for a shadow-circle radius swept
    over time.
    """

    north: list[tuple[float, float]]
    south: list[tuple[float, float]]


def _touches_visible_disk(s: ShadowState, radius: float) -> bool:
    """Whether any point of the shadow circle is on Earth's visible surface,
    checked only at the single point of the circle closest to the
    fundamental plane's origin -- provably the first to enter, and last to
    leave, Earth's visible disk as the circle approaches/recedes from it.
    """
    if radius <= 0:
        return False
    theta = math.atan2(-s.y, -s.x)
    return _lat_at_angle(s, radius, theta) is not None


def _tangent_point(s: ShadowState, radius: float) -> tuple[float, float] | None:
    """Geographic point where the shadow circle is tangent to Earth's
    visible-disk edge (the closest-to-origin point on the circle), used as
    the single start/end point of a limit track so it tapers to a point at
    first/last contact.
    """
    theta = math.atan2(-s.y, -s.x)
    return invert(s.x + radius * math.cos(theta), s.y + radius * math.sin(theta), s.d, s.mu)


def _find_contact_time(
    e: BesselianElements, radius_fn, t_not_touching: float, t_touching: float, iters: int = 50
) -> float:
    """Bisect between a time the shadow circle is known not to touch Earth
    and a time it's known to, to precisely locate the contact instant."""
    for _ in range(iters):
        t_mid = (t_not_touching + t_touching) / 2
        s = shadow_state(e, t_mid)
        radius, _ = radius_fn(s)
        if _touches_visible_disk(s, radius):
            t_touching = t_mid
        else:
            t_not_touching = t_mid
    return t_touching


def _visible_window(
    e: BesselianElements, radius_fn, scan_step_minutes: float = 1.0
) -> tuple[float, float] | None:
    """The precise [t_start, t_end] over which the shadow circle touches
    Earth's visible disk at all, found by a coarse scan across [tmin, tmax]
    to bracket the first/last contact, refined by bisection. None if the
    circle never touches Earth in that range.
    """
    times = sample_times(e.tmin, e.tmax, scan_step_minutes)
    touching = []
    for t in times:
        s = shadow_state(e, t)
        radius, _ = radius_fn(s)
        touching.append(_touches_visible_disk(s, radius))

    if not any(touching):
        return None

    first_idx = next(i for i, v in enumerate(touching) if v)
    last_idx = len(touching) - 1 - next(i for i, v in enumerate(reversed(touching)) if v)

    t_start = (
        _find_contact_time(e, radius_fn, times[first_idx - 1], times[first_idx])
        if first_idx > 0
        else times[first_idx]
    )
    t_end = (
        _find_contact_time(e, radius_fn, times[last_idx + 1], times[last_idx])
        if last_idx < len(times) - 1
        else times[last_idx]
    )
    return t_start, t_end


_EARTH_RADIUS_KM = 6371.0


def _geo_distance_km(p1: tuple[float, float], p2: tuple[float, float]) -> float:
    """Great-circle (haversine) distance in km."""
    lat1, lon1 = math.radians(p1[0]), math.radians(p1[1])
    lat2, lon2 = math.radians(p2[0]), math.radians(p2[1])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * _EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(a)))


# Near first/last contact, the shadow circle's chord of intersection with
# Earth's visible disk grows like sqrt(t - t_contact): an infinite initial
# slope, so no fixed time step avoids a long first segment there.
# Chord-length-based recursive bisection (down to this resolution) is used
# instead of a finer fixed step, so the curve is smooth everywhere without
# oversampling the well-behaved middle of the path.
_MAX_CHORD_KM = 40.0
_MAX_REFINE_DEPTH = 20

# The four mathematically distinct boundary-point formulas of a shadow
# circle at an instant: the two analytic-envelope solutions (valid
# wherever that specific point individually lands on Earth's visible
# surface, independent of whether some unrelated part of the same circle
# is clipped by Earth's edge elsewhere), and the two points where the
# circle crosses Earth's visible-disk edge (only present once the circle
# is clipped somewhere). Each is an unambiguous, continuous function of t.
_SOURCE_KEYS = ("A", "B", "lo", "hi")


def _all_source_points(
    e: BesselianElements, radius_fn, t: float
) -> dict[str, tuple[float, float]]:
    """Every one of `_SOURCE_KEYS` that's currently valid, at time t."""
    s = shadow_state(e, t)
    radius, radius_deriv = radius_fn(s)
    if radius <= 0:
        return {}

    points: dict[str, tuple[float, float]] = {}
    a, b = envelope_points(s, radius, radius_deriv)
    if a is not None:
        points["A"] = a
    if b is not None:
        points["B"] = b

    arc = _find_visible_arc(s, radius)
    if arc is not None:
        theta_lo, theta_hi = arc
        if theta_hi - theta_lo < 2 * math.pi - 1e-9:
            p_lo = invert(
                s.x + radius * math.cos(theta_lo), s.y + radius * math.sin(theta_lo), s.d, s.mu
            )
            p_hi = invert(
                s.x + radius * math.cos(theta_hi), s.y + radius * math.sin(theta_hi), s.d, s.mu
            )
            if p_lo is not None:
                points["lo"] = p_lo
            if p_hi is not None:
                points["hi"] = p_hi
    return points


def _resolve_source(
    points: dict[str, tuple[float, float]],
    current_source: str | None,
    other_source: str | None,
    ref: tuple[float, float],
) -> tuple[str, tuple[float, float]] | None:
    """This track's (source, point) at an instant where `points` are known
    valid: stays on `current_source` if it's still valid there (a track
    only switches source when it actually has to, not just because another
    formula also happens to be independently valid nearby). Otherwise picks
    whichever other available source is nearest to `ref` (excluding
    `other_source`, so the two tracks don't collide onto the same source).
    """
    if current_source is not None and current_source in points:
        return current_source, points[current_source]
    candidates = {k: v for k, v in points.items() if k != other_source} or points
    if not candidates:
        return None
    best = min(candidates, key=lambda k: _geo_distance_km(candidates[k], ref))
    return best, candidates[best]


_TrackState = tuple[tuple[str, tuple[float, float]], tuple[str, tuple[float, float]]]

# A source hand-off should land within roughly a normal sampling step's
# distance of where the track already was. A jump much larger than that
# means the old source vanished at a point unrelated to the other three
# formulas (e.g. where the disk-edge points converge and vanish together
# as a circle finishes becoming fully visible -- a pinch point
# structurally identical to the eclipse's overall first/last-contact
# tangent point, occurring mid-path). Such a jump is treated the same way
# as the very first divergence: re-seeded from whichever two
# currently-available sources are farthest from each other.
_MAX_PLAUSIBLE_HANDOFF_KM = 400.0


def _diverge(
    points: dict[str, tuple[float, float]], north_ref: tuple[float, float]
) -> _TrackState:
    """Seed north/south from whichever two available sources are farthest
    from each other (labeling by proximity to `north_ref` for output
    continuity) -- used both for the initial divergence away from the
    tangent point, and to restart after a pinch point mid-path.
    """
    keys = list(points)
    if len(keys) == 1:
        only = (keys[0], points[keys[0]])
        return only, only
    k1, k2 = max(
        ((a, b) for i, a in enumerate(keys) for b in keys[i + 1 :]),
        key=lambda pair: _geo_distance_km(points[pair[0]], points[pair[1]]),
    )
    p1, p2 = points[k1], points[k2]
    if _geo_distance_km(p1, north_ref) <= _geo_distance_km(p2, north_ref):
        return (k1, p1), (k2, p2)
    return (k2, p2), (k1, p1)


def _eval_tracks(
    e: BesselianElements,
    radius_fn,
    t: float,
    north_source: str | None,
    south_source: str | None,
    north_ref: tuple[float, float],
    south_ref: tuple[float, float],
) -> _TrackState | None:
    """((north_source, north_point), (south_source, south_point)) at time
    t, or None if the shadow circle doesn't touch Earth's visible surface
    at all here.
    """
    points = _all_source_points(e, radius_fn, t)
    if not points:
        return None

    if north_source is None and south_source is None and north_ref == south_ref:
        return _diverge(points, north_ref)

    north = _resolve_source(points, north_source, south_source, north_ref)
    south = (
        _resolve_source(points, south_source, north[0], south_ref) if north is not None else None
    )
    if (
        north is None
        or south is None
        or _geo_distance_km(north[1], north_ref) > _MAX_PLAUSIBLE_HANDOFF_KM
        or _geo_distance_km(south[1], south_ref) > _MAX_PLAUSIBLE_HANDOFF_KM
    ):
        return _diverge(points, north_ref)

    return north, south


def _refine_tracks(
    e: BesselianElements,
    radius_fn,
    t_lo: float,
    state_lo: _TrackState,
    t_hi: float,
    state_hi: _TrackState,
    depth: int,
    out: list[_TrackState],
) -> None:
    (_, n_lo), (_, s_lo) = state_lo
    (_, n_hi), (_, s_hi) = state_hi
    gap = max(_geo_distance_km(n_lo, n_hi), _geo_distance_km(s_lo, s_hi))
    if gap <= _MAX_CHORD_KM or depth >= _MAX_REFINE_DEPTH:
        out.append(state_hi)
        return

    (nsrc_lo, _), (ssrc_lo, _) = state_lo
    t_mid = (t_lo + t_hi) / 2
    state_mid = _eval_tracks(e, radius_fn, t_mid, nsrc_lo, ssrc_lo, n_lo, s_lo)
    if state_mid is None:
        out.append(state_hi)
        return

    _refine_tracks(e, radius_fn, t_lo, state_lo, t_mid, state_mid, depth + 1, out)
    _refine_tracks(e, radius_fn, t_mid, state_mid, t_hi, state_hi, depth + 1, out)


def _limit_track(e: BesselianElements, radius_fn, step_minutes: float) -> LimitTrack:
    """Build a LimitTrack by walking the shadow circle's boundary forward in
    time, tracking north and south each as a "sticky" choice among the four
    boundary-point formulas (`_all_source_points`): each stays on whichever
    formula it's currently using for as long as that formula remains valid,
    and only re-resolves, by nearest-to-last-position, at the point one
    actually stops being valid (`_resolve_source`). Chord-length-based
    recursive bisection (`_refine_tracks`) fills in the sampling wherever
    consecutive points end up further apart than `_MAX_CHORD_KM`.

    The track is sampled only within the shadow circle's true visible
    window (`_visible_window`, not simply [tmin, tmax]) and starts/ends
    exactly at the tangent point of first/last contact so it tapers to a
    point there.
    """
    window = _visible_window(e, radius_fn)
    if window is None:
        return LimitTrack([], [])
    t_start, t_end = window

    start_state = shadow_state(e, t_start)
    end_state = shadow_state(e, t_end)
    start_radius, _ = radius_fn(start_state)
    end_radius, _ = radius_fn(end_state)
    start_point = _tangent_point(start_state, start_radius)
    end_point = _tangent_point(end_state, end_radius)
    if start_point is None or end_point is None:
        return LimitTrack([], [])

    north = [start_point]
    south = [start_point]

    times = [t for t in sample_times(t_start, t_end, step_minutes) if t_start < t < t_end]
    prev_t = t_start
    prev_state: _TrackState = ((None, start_point), (None, start_point))  # type: ignore[arg-type]
    for t in times + [t_end]:
        (nsrc, _), (ssrc, _) = prev_state
        (_, n_ref), (_, s_ref) = prev_state
        state = _eval_tracks(e, radius_fn, t, nsrc, ssrc, n_ref, s_ref)
        if state is None:
            prev_t = t
            continue

        refined: list[_TrackState] = []
        _refine_tracks(e, radius_fn, prev_t, prev_state, t, state, 0, refined)
        for (_, n_pt), (_, s_pt) in refined:
            north.append(n_pt)
            south.append(s_pt)

        prev_state = refined[-1] if refined else state
        prev_t = t

    north.append(end_point)
    south.append(end_point)
    return LimitTrack(north, south)


def umbral_limits(e: BesselianElements, step_minutes: float = 0.5) -> LimitTrack:
    """Northern and southern limits of the path of totality/annularity.
    Empty if the eclipse is partial-only.
    """
    if not is_central(e):
        return LimitTrack([], [])

    def radius_fn(s: ShadowState) -> tuple[float, float]:
        sign = 1.0 if s.l2 >= 0 else -1.0
        return abs(s.l2), sign * s.l2p

    return _limit_track(e, radius_fn, step_minutes)


def limits_to_ring(track: LimitTrack) -> list[tuple[float, float]]:
    """Combine a LimitTrack into a single closed polygon ring: north
    west-to-east, then south east-to-west back to the start. Empty if the
    shadow circle never touched the Earth over the sampled time range.
    """
    if not track.north or not track.south:
        return []
    ring = list(track.north) + list(reversed(track.south))
    ring.append(ring[0])
    return ring
