"""Multi-format RNA secondary-structure file reading for the editor's Open.

`parse_structure_file(path)` auto-detects and parses the common RNA structure
file formats into a `(ss, seq, name)` triple the editor can lay out:

* **dot-bracket** (`.dbn`, `.dot`, `.txt`, `.dat`, `.ss`) -- bpRNA style:
  `#`-comment / header lines are skipped, then a sequence line (optional) and
  the dot-bracket line are read. Pseudoknot brackets ``[]{}<>`` are preserved.
* **FASTA-with-structure** (`.fasta`, `.fa`) -- a ``>name`` header, one or more
  sequence lines, then a structure line.
* **connectivity table** (`.ct`) -- a header (length [+ name]) then rows
  ``i base prev next pair k``; the partner column is turned into dot-bracket.
* **base-pair sequence** (`.bpseq`) -- rows ``i base j``; the partner column is
  turned into dot-bracket.

For the table formats, crossing (pseudoknotted) pairs are emitted across the
four bracket families ``()[]{}<>`` (greedy non-crossing assignment); pairs too
knotted to place in any of the four families are dropped (rare in practice).

The dot-bracket string this returns is fed straight to
`EditorModel.from_ss` / `MainWindow.load_ss`, so it uses exactly the alphabet
`render_rna.get_pairmap_from_secstruct` and the pseudoknot parser understand.
"""

from __future__ import annotations

import os

# Bracket families used to render pseudoknots when rebuilding dot-bracket from a
# pair list (CT / bpseq). Same four families the pseudoknot parser reads.
_FAMILIES: tuple[tuple[str, str], ...] = (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"))

# The full dot-bracket alphabet (all four families plus the unpaired dot).
_DB_CHARS = frozenset(".") | {c for pair in _FAMILIES for c in pair}

# Nucleotide letters a sequence line may contain (upper/lower, RNA/DNA, gaps).
_SEQ_CHARS = frozenset("ACGUTNRYSWKMBDHVacgutnryswkmbdhv-")


def _is_dot_bracket(line: str) -> bool:
    """True when every character of ``line`` is a dot-bracket symbol."""
    return bool(line) and set(line) <= _DB_CHARS


def _crosses(p: tuple[int, int], q: tuple[int, int]) -> bool:
    """True when base pairs ``p`` and ``q`` cross (are pseudoknotted)."""
    i, j = p
    a, b = q
    return (i < a < j < b) or (a < i < b < j)


def pairs_to_dot_bracket(n: int, pairs: list[tuple[int, int]]) -> str:
    """Render a set of base pairs as a length-``n`` dot-bracket string.

    Non-crossing pairs use ``()``; crossing (pseudoknotted) pairs spill over
    into ``[]``, ``{}`` then ``<>`` via a greedy first-fit assignment (each
    pair goes to the first family holding no pair it crosses). A pair that
    crosses a placed pair in all four families is dropped.

    Args:
        n: Total structure length.
        pairs: Base pairs as ``(i, j)`` (either order); duplicates/self-pairs
            are ignored.

    Returns:
        A dot-bracket string of length ``n``.
    """
    clean = sorted({(min(i, j), max(i, j)) for i, j in pairs if i != j and 0 <= min(i, j)})
    chars = ["."] * n
    families: list[list[tuple[int, int]]] = [[] for _ in _FAMILIES]
    for pair in clean:
        i, j = pair
        if j >= n:
            continue
        for placed, (opener, closer) in zip(families, _FAMILIES):
            if all(not _crosses(pair, other) for other in placed):
                placed.append(pair)
                chars[i], chars[j] = opener, closer
                break
    return "".join(chars)


def parse_structure_file(path: str) -> tuple[str, str | None, str | None]:
    """Parse an RNA structure file into ``(ss, seq, name)``.

    The format is chosen from the file extension, then confirmed / recovered
    from the content (so a mislabeled file still parses when its body is
    recognizable). ``seq`` and ``name`` are ``None`` when the file carries
    neither.

    Args:
        path: Path to the structure file.

    Returns:
        ``(ss, seq, name)``: the dot-bracket structure, the sequence (or
        ``None``), and a name/title (or ``None``).

    Raises:
        ValueError: When no structure can be extracted from the file.
    """
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    ext = os.path.splitext(path)[1].lower()

    if ext == ".ct":
        return _parse_ct(text)
    if ext == ".bpseq":
        return _parse_bpseq(text)

    # Everything else is dot-bracket / FASTA text; fall back to the table
    # parsers if the body clearly is a CT / bpseq that was mislabeled.
    try:
        return _parse_dot_bracket(text)
    except ValueError:
        for parser in (_parse_ct, _parse_bpseq):
            try:
                return parser(text)
            except ValueError:
                continue
        raise ValueError(f"could not parse a structure from {path!r}")


def _parse_dot_bracket(text: str) -> tuple[str, str | None, str | None]:
    """Parse bpRNA-style dot-bracket / FASTA-with-structure text.

    Skips ``#`` comments (capturing a ``Name:`` if present) and ``>`` FASTA
    headers, accumulates sequence line(s), and stops at the first dot-bracket
    line (the structure). Multi-line FASTA sequences are concatenated.
    """
    name: str | None = None
    seq_parts: list[str] = []
    struct: str | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is None:
                name = line[1:].strip() or None
            continue
        if line.startswith("#"):
            low = line.lower()
            if name is None and "name" in low and ":" in line:
                name = line.split(":", 1)[1].strip() or None
            continue
        if _is_dot_bracket(line):
            struct = line
            break
        # A plain letter line: part of the sequence.
        if set(line) <= _SEQ_CHARS:
            seq_parts.append(line)
    if struct is None:
        raise ValueError("no dot-bracket structure line found")
    seq = "".join(seq_parts) or None
    if seq is not None and len(seq) != len(struct):
        # Length mismatch means the letter line was not this structure's
        # sequence (stray annotation); drop it rather than mislabel.
        seq = None
    return struct, seq, name


def _iter_data_rows(text: str) -> list[list[str]]:
    """Whitespace-tokenized non-comment lines (``#``/``>`` stripped)."""
    rows: list[list[str]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith(">"):
            continue
        rows.append(line.split())
    return rows


def _parse_ct(text: str) -> tuple[str, str | None, str | None]:
    """Parse a connectivity-table (.ct) file into ``(ss, seq, name)``.

    Header: ``<length> [name/energy ...]``. Data rows:
    ``index base prev next pair k`` (1-based ``index``/``pair``; ``pair == 0``
    means unpaired). Only the first ``length`` data rows are read (some CT
    files concatenate several structures).
    """
    rows = _iter_data_rows(text)
    if not rows:
        raise ValueError("empty connectivity table")
    header = rows[0]
    try:
        length = int(header[0])
    except (ValueError, IndexError):
        raise ValueError("connectivity table header lacks a length")
    name = " ".join(header[1:]).strip() or None
    seq_chars: list[str] = []
    pairs: list[tuple[int, int]] = []
    count = 0
    for row in rows[1:]:
        if len(row) < 6:
            continue
        try:
            i = int(row[0])
            partner = int(row[4])
        except ValueError:
            continue
        seq_chars.append(row[1])
        if partner > 0 and partner > i:
            pairs.append((i - 1, partner - 1))
        count += 1
        if count >= length:
            break
    if count == 0:
        raise ValueError("connectivity table has no data rows")
    ss = pairs_to_dot_bracket(count, pairs)
    seq = "".join(seq_chars) or None
    return ss, seq, name


def _parse_bpseq(text: str) -> tuple[str, str | None, str | None]:
    """Parse a base-pair-sequence (.bpseq) file into ``(ss, seq, None)``.

    Data rows: ``index base pair`` (1-based; ``pair == 0`` means unpaired).
    """
    rows = _iter_data_rows(text)
    seq_chars: list[str] = []
    pairs: list[tuple[int, int]] = []
    n = 0
    for row in rows:
        if len(row) < 3:
            continue
        try:
            i = int(row[0])
            partner = int(row[2])
        except ValueError:
            continue
        if len(row[1]) != 1 or row[1] not in _SEQ_CHARS:
            continue
        seq_chars.append(row[1])
        n += 1
        if partner > 0 and partner > i:
            pairs.append((i - 1, partner - 1))
    if n == 0:
        raise ValueError("base-pair sequence has no data rows")
    ss = pairs_to_dot_bracket(n, pairs)
    seq = "".join(seq_chars) or None
    return ss, seq, None


# File-dialog filter fragment shared by the editor's Open action.
OPEN_FILTER = (
    "RNA documents (*.rnadoc.json *.rnadoc);;"
    "RNA structures (*.dbn *.dot *.ct *.bpseq *.fasta *.fa *.txt *.dat *.ss);;"
    "All files (*)"
)


__all__ = ["parse_structure_file", "pairs_to_dot_bracket", "OPEN_FILTER"]
