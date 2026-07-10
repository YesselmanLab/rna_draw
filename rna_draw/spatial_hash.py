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
"""

from math import floor


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
