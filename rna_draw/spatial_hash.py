"""Uniform-grid spatial hash for candidate overlap-pair generation.

Correctness argument: two footprints overlap implies their (already
inflated, i.e. tolerance/half-width-padded) axis-aligned bounding boxes
overlap, which implies some grid cell is covered by both boxes, which
implies both items were inserted into that cell (``insert`` adds an item
to every cell its AABB touches), which implies the pair is emitted by
``candidate_pairs``. Therefore ``candidate_pairs`` is a *superset* of every
true overlap for ANY positive ``cell_size`` -- the cell size can only ever
affect speed, never correctness. Expected performance is close to O(n) for
RNA layouts because nucleotide density is locally bounded (nucleotides are
~``PRIMARY_SPACE`` apart), so each cell holds few items.

``insert_segment`` extends this argument to long, thin footprints (e.g. a
base-pair chord drawn as a diagonal capsule): it splits the segment into
pieces short enough that each piece's own padded AABB is small, and
inserts each piece separately instead of one big AABB spanning the whole
diagonal. The superset argument above applies unchanged to each piece --
any point within ``pad`` of the full segment is within ``pad`` of the
particular piece nearest to it (each coordinate of that point differs from
the piece's segment-coordinate bounds by at most ``pad``, by the same
per-axis distance bound the whole-segment case relies on) -- so the union
of the pieces' cells is still a superset of the true footprint's cells.
This turns a diagonal segment from touching ``O((length / cell_size) ** 2)``
cells (one big square AABB over a thin diagonal) into ``O(length /
cell_size)`` cells (a thin tiled strip), which is what made large-structure
checks with long chords (e.g. a circle-layout fallback) slow.
"""

from math import ceil, floor, hypot


class SpatialHash:
    """Buckets items by the grid cells their AABBs touch.

    Args:
        cell_size: Side length of a square grid cell. A performance knob
            only (see module docstring); does not affect correctness.
    """

    def __init__(self, cell_size: float) -> None:
        self._cell_size = cell_size
        self._cells: dict[tuple[int, int], list[int]] = {}

    def insert(self, item_index: int, aabb: tuple[float, float, float, float]) -> None:
        """Insert an item into every grid cell its AABB touches.

        Args:
            item_index: Identifier for the item (e.g. its index in a
                primitive list).
            aabb: Bounding box as ``(minx, miny, maxx, maxy)``.
        """
        minx, miny, maxx, maxy = aabb
        fx0, fx1 = floor(minx / self._cell_size), floor(maxx / self._cell_size)
        fy0, fy1 = floor(miny / self._cell_size), floor(maxy / self._cell_size)
        for cx in range(fx0, fx1 + 1):
            for cy in range(fy0, fy1 + 1):
                self._cells.setdefault((cx, cy), []).append(item_index)

    def insert_segment(
        self, item_index: int, x0: float, y0: float, x1: float, y1: float, pad: float
    ) -> None:
        """Insert an item as a segment, tiled into ``pad``-padded pieces.

        Splits the ``(x0, y0)``-``(x1, y1)`` segment into
        ``ceil(length / cell_size)`` equal-length pieces and inserts each
        piece's own bounding box, padded by ``pad`` on every side (the same
        style of inflation ``insert``'s callers use for a whole AABB). See
        the module docstring for why this stays a superset of the true
        footprint's cells while avoiding one huge AABB over a long
        diagonal segment.

        Args:
            item_index: Identifier for the item.
            x0: Segment start x-coordinate.
            y0: Segment start y-coordinate.
            x1: Segment end x-coordinate.
            y1: Segment end y-coordinate.
            pad: Inflation applied to each piece's bounding box, e.g. a
                capsule's half-width.
        """
        length = hypot(x1 - x0, y1 - y0)
        if length == 0.0:
            self.insert(item_index, (x0 - pad, y0 - pad, x0 + pad, y0 + pad))
            return
        steps = max(1, ceil(length / self._cell_size))
        dx, dy = (x1 - x0) / steps, (y1 - y0) / steps
        for k in range(steps):
            px0, py0 = x0 + k * dx, y0 + k * dy
            px1, py1 = px0 + dx, py0 + dy
            aabb = (
                min(px0, px1) - pad,
                min(py0, py1) - pad,
                max(px0, px1) + pad,
                max(py0, py1) + pad,
            )
            self.insert(item_index, aabb)

    def candidate_pairs(self) -> set[tuple[int, int]]:
        """Return every ``(i, j)`` with ``i < j`` sharing at least one cell.

        Returns:
            Set of candidate index pairs; a superset of all truly
            AABB-overlapping pairs (see module docstring).
        """
        pairs: set[tuple[int, int]] = set()
        for items in self._cells.values():
            n = len(items)
            for a in range(n):
                for b in range(a + 1, n):
                    i, j = items[a], items[b]
                    pairs.add((i, j) if i < j else (j, i))
        return pairs
