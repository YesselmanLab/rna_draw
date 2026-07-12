"""Multi-bracket dot-bracket parsing + stacked-run stem grouping (Phase 1).

`render_rna.get_pairmap_from_secstruct` (frozen) stays the arbiter of the
`()`-only view used everywhere else in `rna_draw`; this module is the ONLY
place that reads `[]{}<>` as base pairs, feeding the pseudoknot layout path
(`rna_draw.layout.pseudoknot.engine`).
"""

from __future__ import annotations

from dataclasses import dataclass

BRACKET_PAIRS: tuple[tuple[str, str], ...] = (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"))
_ALLOWED_CHARS = frozenset(".") | {char for pair in BRACKET_PAIRS for char in pair}


@dataclass(frozen=True)
class Stem:
    """A maximal run of stacked base pairs -- one physical helix.

    Args:
        i: 5' index of the stem's outermost pair.
        j: 3' index of the stem's outermost pair.
        length: Number of stacked rungs; rung `k` (`0 <= k < length`) pairs
            nucleotide `i + k` with `j - k`.
    """

    i: int
    j: int
    length: int

    def pairs(self) -> list[tuple[int, int]]:
        """Every base pair in this stem, outermost rung first."""
        return [(self.i + k, self.j - k) for k in range(self.length)]


def parse_all_pairs(secstruct: str) -> list[tuple[int, int]]:
    """Parse every base pair across all four dot-bracket types.

    Args:
        secstruct: Dot-bracket string using `.` plus any of `()[]{}<>`.

    Returns:
        Every `(i, j)`, `i < j`, sorted by `i`. An unmatched bracket (a
        malformed input) is simply left unpaired.

    Raises:
        ValueError: If `secstruct` contains a character outside
            `.()[]{}<>`.
    """
    bad = set(secstruct) - _ALLOWED_CHARS
    if bad:
        raise ValueError(f"secstruct contains unsupported character(s): {sorted(bad)!r}")
    pairs: list[tuple[int, int]] = []
    for opener, closer in BRACKET_PAIRS:
        pairs.extend(_parse_one_type(secstruct, opener, closer))
    return sorted(pairs)


def _parse_one_type(secstruct: str, opener: str, closer: str) -> list[tuple[int, int]]:
    """Stack-match one bracket type's openers against its closers."""
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for idx, char in enumerate(secstruct):
        if char == opener:
            stack.append(idx)
        elif char == closer and stack:
            pairs.append((stack.pop(), idx))
    return pairs


def group_stems(pairs: list[tuple[int, int]]) -> list[Stem]:
    """Group base pairs into maximal stacked runs (stems).

    Consecutive pairs `(i, j)` and `(i + 1, j - 1)` stack into one stem
    regardless of which bracket type produced them -- a stem is a purely
    index-geometric notion: consecutively-shrinking partners form one
    physical helix. (In practice a stem's rungs always share one bracket
    type in real input, since crossing a bracket-type boundary breaks the
    `partner[i] == j` chain; see the H-type example in the module tests.)

    Args:
        pairs: Every base pair (e.g. from `parse_all_pairs`), each
            `i < j`.

    Returns:
        One `Stem` per maximal stacked run, sorted by `i`.
    """
    partner = dict(pairs)
    used: set[int] = set()
    stems = []
    for i, j in sorted(pairs):
        if i in used:
            continue
        stems.append(_walk_stack(partner, used, i, j))
    return stems


def _walk_stack(partner: dict[int, int], used: set[int], i: int, j: int) -> Stem:
    """Follow a stacked run of pairs starting at `(i, j)`, marking each used."""
    start_i, start_j = i, j
    length = 0
    while partner.get(i) == j:
        used.add(i)
        length += 1
        i, j = i + 1, j - 1
    return Stem(i=start_i, j=start_j, length=length)


def stem_pairs(stems: list[Stem]) -> list[tuple[int, int]]:
    """Flatten every pair from a list of stems.

    Args:
        stems: Stems to flatten.

    Returns:
        Every `(i, j)` pair from every stem, in stem order then
        outermost-rung-first within a stem.
    """
    return [pair for stem in stems for pair in stem.pairs()]


def nested_secstruct(n: int, retained_pairs: list[tuple[int, int]]) -> str:
    """Emit a `()`-only dot-bracket string from a nested (non-crossing) pair set.

    Args:
        n: Total structure length.
        retained_pairs: Pairwise non-crossing base pairs (e.g. every pair
            from `max_nested_subset`'s `retained` stems, via `stem_pairs`).

    Returns:
        A dot-bracket string of length `n` using only `()`. Feeding it
        through `render_rna.get_pairmap_from_secstruct` reproduces
        `retained_pairs` exactly (Phase 1's round-trip guarantee, enforced
        by `pseudoknot.engine._assert_round_trip`).
    """
    pair_map = [-1] * n
    for i, j in retained_pairs:
        pair_map[i] = j
        pair_map[j] = i
    return "".join(_char_for(idx, partner) for idx, partner in enumerate(pair_map))


def _char_for(idx: int, partner: int) -> str:
    """The dot-bracket character for one `pair_map` entry."""
    if partner == -1:
        return "."
    return "(" if partner > idx else ")"


__all__ = [
    "BRACKET_PAIRS",
    "Stem",
    "parse_all_pairs",
    "group_stems",
    "stem_pairs",
    "nested_secstruct",
]
