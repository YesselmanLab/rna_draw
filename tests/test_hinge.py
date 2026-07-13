"""Tests for the VARNA loop-redistribution hinge (`rna_draw.gui.hinge` +
`EditorModel.rotate_selection`).

The hinge rotates a selected helix about its supporting loop's center and
then re-spaces that loop's unpaired members evenly on the loop circle, so the
loop stays a clean, compact ring (as VARNA does) rather than distorting. These
tests assert: the loop stays round after rotation, sibling helices stay rigid,
only the dragged helix + the loop's own unpaired members move, and the
never-silent overlap flag still reflects the real native checker.
"""

from __future__ import annotations

import math

import pytest

from rna_draw.gui.hinge import redistribute_loop
from rna_draw.gui.model import EditorModel
from rna_draw.overlap_native import check_overlaps_native

# Three-way junction: click index 7 selects the whole hairpin helix [7, 18]
# whose parent loop is the central multiloop (unpaired members 4,5,6 / 19,20,21
# / 34,35,36).
_SS = "((((...((((....))))...((((....))))...))))"
_CLICK = 7
_HELIX = (7, 18)
_SIBLING = (22, 33)  # the other multiloop child helix


def _roundness(model, loop) -> float:
    """Std/mean of loop-member radii about the pivot (0 == perfect circle)."""
    cx, cy = model.pivot
    rs = [math.hypot(model._x[m] - cx, model._y[m] - cy) for m in loop.members]
    mean = sum(rs) / len(rs)
    var = sum((r - mean) ** 2 for r in rs) / len(rs)
    return math.sqrt(var) / mean


def _internal_distances(x, y, start, end):
    """All pairwise distances within a slice (a rigidity fingerprint)."""
    idx = range(start, end + 1)
    return [
        math.hypot(x[i] - x[j], y[i] - y[j]) for i in idx for j in idx if i < j
    ]


def test_redistribute_keeps_loop_round():
    # After a substantial rotation the redistributing hinge keeps the loop at
    # least as round as the rigid-only method (unpaired members land exactly
    # on the loop circle).
    angle = math.radians(35)

    m_old = EditorModel.from_ss(_SS)
    m_old.select(_CLICK)
    m_old.rotate_selection(angle, redistribute=False)
    round_old = _roundness(m_old, m_old._sel_loop)

    m_new = EditorModel.from_ss(_SS)
    m_new.select(_CLICK)
    m_new.rotate_selection(angle, redistribute=True)
    round_new = _roundness(m_new, m_new._sel_loop)

    assert round_new <= round_old + 1e-9
    # and the redistributed loop is genuinely tight (small radial spread)
    assert round_new < 0.05


def test_unpaired_members_land_on_the_circle():
    # Every unpaired member of the loop sits at (near) the mean anchor radius.
    m = EditorModel.from_ss(_SS)
    m.select(_CLICK)
    m.rotate_selection(math.radians(40), redistribute=True)
    loop = m._sel_loop
    cx, cy = m.pivot
    anchors = {loop.closing_pair[0], loop.closing_pair[1]}
    for b in loop.children:
        anchors.update((b.start, b.end))
    radius = sum(
        math.hypot(m._x[a] - cx, m._y[a] - cy) for a in anchors
    ) / len(anchors)
    for u in (4, 5, 6, 19, 20, 21, 34, 35, 36):
        r = math.hypot(m._x[u] - cx, m._y[u] - cy)
        assert r == pytest.approx(radius, rel=0.02)


def test_sibling_helix_stays_rigid():
    # The other child helix's internal shape is exactly preserved (not moved,
    # not distorted) by rotating + redistributing around it.
    m0 = EditorModel.from_ss(_SS)
    before = _internal_distances(m0._x, m0._y, *_SIBLING)

    m = EditorModel.from_ss(_SS)
    m.select(_CLICK)
    m.rotate_selection(math.radians(35), redistribute=True)
    after = _internal_distances(m._x, m._y, *_SIBLING)

    assert after == pytest.approx(before, abs=1e-9)
    # the sibling did not translate either -- its anchor is fixed
    assert (m._x[_SIBLING[0]], m._y[_SIBLING[0]]) == pytest.approx(
        (m0._x[_SIBLING[0]], m0._y[_SIBLING[0]]), abs=1e-9
    )


def test_selected_helix_stays_rigid():
    # The dragged helix rotates rigidly: its internal pairwise distances are
    # unchanged (redistribution never touches a helix's own nucleotides).
    m0 = EditorModel.from_ss(_SS)
    before = _internal_distances(m0._x, m0._y, *_HELIX)

    m = EditorModel.from_ss(_SS)
    m.select(_CLICK)
    m.rotate_selection(math.radians(35), redistribute=True)
    after = _internal_distances(m._x, m._y, *_HELIX)
    assert after == pytest.approx(before, abs=1e-9)


def test_only_helix_and_loop_unpaired_move():
    # Redistribution moves ONLY the dragged helix slice and the loop's own
    # unpaired members; every other nucleotide is untouched.
    m = EditorModel.from_ss(_SS)
    m.select(_CLICK)
    before = list(zip(m._x, m._y))
    m.rotate_selection(math.radians(30), redistribute=True)
    after = list(zip(m._x, m._y))

    allowed_to_move = set(range(_HELIX[0], _HELIX[1] + 1)) | {
        4, 5, 6, 19, 20, 21, 34, 35, 36
    }
    for i, (b, a) in enumerate(zip(before, after)):
        if i in allowed_to_move:
            continue
        assert a == pytest.approx(b, abs=1e-9), f"nt {i} moved but should not"


def test_never_silent_flag_matches_real_checker():
    # A large rotation drives helices into a collision; the model's flag must
    # equal the real native checker's verdict at the model's scaled params.
    m = EditorModel.from_ss(_SS)
    m.select(_CLICK)
    result = m.rotate_selection(math.radians(60), redistribute=True)
    params = m._scaled_params()
    count = check_overlaps_native(m._x, m._y, m._pair_map, params)
    assert result.flagged == (count > 0)
    assert m.flagged == (count > 0)


def test_redistribute_noop_without_unpaired_members():
    # Selecting a stacked inner pair gives a stack-continuation parent loop
    # with no unpaired members; redistribution must be a pure no-op there, so
    # the rotation stays exactly rigid.
    from rna_draw.layout.postpass import rotate_range

    m = EditorModel.from_ss(_SS)
    start, end, pivot = m.select(8)  # inner stacked pair -> degenerate parent
    nx, ny = rotate_range(m._x, m._y, start, end, pivot[0], pivot[1], 0.05)
    m.rotate_selection(0.05, redistribute=True)
    assert m._x == pytest.approx(nx)
    assert m._y == pytest.approx(ny)


def test_redistribute_loop_exterior_is_noop():
    # An exterior loop has no closing pair / circle; redistribute_loop returns
    # the coordinates unchanged.
    m = EditorModel.from_ss(_SS)
    exterior = m._tree.exterior
    nx, ny = redistribute_loop(exterior, m._x, m._y, (0.0, 0.0))
    assert nx == m._x
    assert ny == m._y
