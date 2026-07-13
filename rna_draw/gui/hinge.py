"""VARNA-faithful loop redistribution for the interactive helix hinge.

When an editor rotates a selected helix about its supporting loop's center,
the helix's rigid slice moves to a new angular position on that loop's
circle -- but the loop's own unpaired nucleotides must then re-space evenly
around the circle so the loop stays a clean, compact ring instead of
distorting (the bug this module fixes). This mirrors what VARNA does on an
interactive helix rotation.

VARNA model (confirmed from `fr.orsay.lri.varna.models.rna.RNA` bytecode in
VARNA.jar): the interactive rotation entry point `rotateEverything` calls
`rotateHelix(center, i, j, angle)` -- a rigid rotation of ONLY the dragged
helix's slice `[i, j]` about the loop center (subsequent sibling helices are
NOT carried) -- and then `fixUnpairedPositions` -> `distributeUnpaired`,
which re-places the unpaired nucleotides flanking the helix EVENLY BY ANGLE
along the loop circle (`center + radius * (cos angle, sin angle)`) in the
arcs between the helix's new endpoints and the neighbouring anchors. So only
the dragged helix moves; the loop's unpaired members redistribute; every
other anchor (the loop's closing pair and every sibling helix's base pair)
stays put.

`redistribute_loop` reproduces that redistribution. It reuses the same
circle geometry the engine's radial layout uses (`render_rna.
setup_coords_recursive` places every loop member at an equal arc-length
position on `center + radius * (cos, sin)`; the constructive engine's
`geometry_helpers.loop_member_point` is the same formula): unpaired members
are placed on the loop circle at angles linearly interpolated between their
two bounding anchors. The helix slice is expected to have been rotated about
`center` already, so its endpoints are still on the loop circle at their new
angles when this runs.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from rna_draw.layout.structure_tree import Loop


def _normalize_pi(angle: float) -> float:
    """Wrap an angle to the half-open interval ``(-pi, pi]``.

    Args:
        angle: An angle in radians.

    Returns:
        The equivalent angle in ``(-pi, pi]``.
    """
    wrapped = math.fmod(angle, 2.0 * math.pi)
    if wrapped <= -math.pi:
        wrapped += 2.0 * math.pi
    elif wrapped > math.pi:
        wrapped -= 2.0 * math.pi
    return wrapped


def _ring_points(loop: Loop) -> list[tuple[int, bool]]:
    """Loop members in geometric ring order, tagged anchor vs. unpaired.

    A loop's boundary is a ring of points on its circle. `Loop.members`
    lists them in structure order but records only the FIRST index of each
    child branch (its closing partner belongs to the branch's rigid slice);
    geometrically both endpoints of a child's base pair sit on the loop
    circle, adjacent in ring order. This expands each child branch's member
    into its two on-circle endpoints and tags every point:

    * anchors (``True``) -- fixed: the loop's own closing-pair endpoints and
      each child branch's two base-pair endpoints. Redistribution never
      moves these.
    * unpaired (``False``) -- the loop's own free nucleotides, the only
      points redistribution re-places.

    Args:
        loop: An interior loop (``closing_pair is not None``).

    Returns:
        ``(index, is_anchor)`` pairs in ring order from the closing pair's
        5' endpoint to its 3' endpoint.
    """
    closing = loop.closing_pair
    child_start_to_end = {b.start: b.end for b in loop.children}
    ring: list[tuple[int, bool]] = []
    for m in loop.members:
        if closing is not None and (m == closing[0] or m == closing[1]):
            ring.append((m, True))
        elif m in child_start_to_end:
            ring.append((m, True))
            ring.append((child_start_to_end[m], True))
        else:
            ring.append((m, False))
    return ring


def _winding_sign(
    ring: Sequence[tuple[int, bool]],
    x: Sequence[float],
    y: Sequence[float],
    center: tuple[float, float],
) -> float:
    """Direction (+1 CCW / -1 CW) the ring sweeps from first anchor to last.

    Summing the shortest signed angular step between consecutive anchors
    (in ring order) recovers the loop's winding: an interior loop is an open
    chain from its closing pair's 5' endpoint round to its 3' endpoint, so
    this sum is close to ``+-(2 pi - reserved sector)`` and its sign is a
    robust read of the traversal direction.

    Args:
        ring: The loop's ring points (see `_ring_points`).
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        center: ``(cx, cy)``, the loop circle center.

    Returns:
        ``+1.0`` for counter-clockwise, ``-1.0`` for clockwise (defaults to
        ``+1.0`` when the sum is exactly zero, e.g. fewer than two anchors).
    """
    cx, cy = center
    anchor_angles = [
        math.atan2(y[i] - cy, x[i] - cx) for i, is_anchor in ring if is_anchor
    ]
    total = 0.0
    for a, b in zip(anchor_angles, anchor_angles[1:]):
        total += _normalize_pi(b - a)
    return -1.0 if total < 0.0 else 1.0


def _directed_sweep(start_angle: float, end_angle: float, sign: float) -> float:
    """Positive angular distance from ``start_angle`` to ``end_angle`` along ``sign``.

    Args:
        start_angle: Angle of the anchor beginning the arc (radians).
        end_angle: Angle of the anchor ending the arc (radians).
        sign: Winding direction (`+1` CCW, `-1` CW) from `_winding_sign`.

    Returns:
        The non-negative sweep (in ``(0, 2 pi]``) to travel from
        ``start_angle`` to ``end_angle`` in the ring's winding direction.
    """
    raw = (end_angle - start_angle) * sign
    swept = math.fmod(raw, 2.0 * math.pi)
    if swept <= 0.0:
        swept += 2.0 * math.pi
    return swept


def _mean_anchor_radius(
    ring: Sequence[tuple[int, bool]],
    x: Sequence[float],
    y: Sequence[float],
    center: tuple[float, float],
) -> float:
    """Mean anchor distance to the loop center -- the loop circle's radius.

    Anchors (base-pair endpoints) define the loop circle, so their mean
    distance to the center is the radius unpaired members are placed at.

    Args:
        ring: The loop's ring points (see `_ring_points`).
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        center: ``(cx, cy)``, the loop circle center.

    Returns:
        The mean anchor-to-center distance, or ``0.0`` if there are no
        anchors.
    """
    cx, cy = center
    dists = [
        math.hypot(x[i] - cx, y[i] - cy) for i, is_anchor in ring if is_anchor
    ]
    if not dists:
        return 0.0
    return sum(dists) / len(dists)


def redistribute_loop(
    loop: Loop,
    x: Sequence[float],
    y: Sequence[float],
    center: tuple[float, float],
) -> tuple[list[float], list[float]]:
    """Re-space a loop's unpaired members evenly on its circle (VARNA hinge).

    Places every unpaired member of ``loop`` on the loop circle
    (``center + radius * (cos, sin)``) at an angle linearly interpolated
    between its two bounding anchors, so each run of unpaired nucleotides is
    spread evenly along the arc between the two helices/closing-pair
    endpoints that bracket it. Anchors -- the loop's closing pair and every
    child helix's base pair -- are left exactly where they are, so no
    sibling helix is moved or distorted and the caller's already-rotated
    helix keeps its new position. See the module docstring for the VARNA
    model this reproduces.

    Args:
        loop: The loop to redistribute. Must be an interior loop
            (``closing_pair is not None``); an exterior loop has no circle
            and is returned unchanged.
        x: Current nucleotide x-coordinates (helix slice already rotated).
        y: Current nucleotide y-coordinates.
        center: ``(cx, cy)``, the loop circle center -- the SAME pivot the
            helix was rotated about, so the rotated helix's endpoints still
            lie on this circle.

    Returns:
        ``(nx, ny)``: new coordinate lists, equal to the inputs everywhere
        except at ``loop``'s own unpaired members.
    """
    nx, ny = list(x), list(y)
    if loop.closing_pair is None:
        return nx, ny

    ring = _ring_points(loop)
    radius = _mean_anchor_radius(ring, x, y, center)
    if radius <= 0.0:
        return nx, ny
    cx, cy = center
    sign = _winding_sign(ring, x, y, center)

    # Walk the ring; for each maximal run of unpaired points bracketed by an
    # anchor before and an anchor after, spread the run evenly by angle.
    i = 0
    n = len(ring)
    while i < n:
        idx, is_anchor = ring[i]
        if is_anchor:
            i += 1
            continue
        run_start = i
        while i < n and not ring[i][1]:
            i += 1
        run_end = i  # exclusive; ring[i] is the trailing anchor (if any)
        before = run_start - 1
        if before < 0 or run_end >= n:
            # No bracketing anchor on one side (only happens for a
            # malformed/exterior ring); leave this run untouched.
            continue
        a_idx = ring[before][0]
        b_idx = ring[run_end][0]
        start_angle = math.atan2(y[a_idx] - cy, x[a_idx] - cx)
        end_angle = math.atan2(y[b_idx] - cy, x[b_idx] - cx)
        sweep = _directed_sweep(start_angle, end_angle, sign)
        count = run_end - run_start
        for k in range(count):
            frac = (k + 1) / (count + 1)
            angle = start_angle + sign * sweep * frac
            m = ring[run_start + k][0]
            nx[m] = cx + radius * math.cos(angle)
            ny[m] = cy + radius * math.sin(angle)
    return nx, ny


__all__ = ["redistribute_loop"]
