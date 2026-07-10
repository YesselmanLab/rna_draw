"""Shared test fixtures/helpers for the `rna_draw.layout` test suite."""

from __future__ import annotations

import random


def random_structure(seed: int, n: int) -> str:
    """Generate a random, pseudoknot-free dot-bracket structure of length `n`.

    Walks left to right maintaining the invariant "every open pair can
    still close within the remaining positions" (`open_count <= slots
    remaining`); whichever of `.`, `(`, `)` keep that invariant true are
    legal at each step, and a step that would violate it (only possible
    for `)`, once every remaining slot is needed to close existing opens)
    is forced. The result is always well-nested and exactly length `n`,
    using only `().`.

    Args:
        seed: RNG seed for reproducibility.
        n: Desired structure length.

    Returns:
        A random pseudoknot-free dot-bracket string of length `n`.
    """
    rng = random.Random(seed)
    chars: list[str] = []
    open_count = 0
    for i in range(n):
        remaining_after = n - i - 1
        if open_count > remaining_after:
            choice = ")"  # forced: this slot must close an existing open
        else:
            options = ["."]
            if open_count > 0:
                options.append(")")
            if open_count + 1 <= remaining_after:
                options.append("(")
            choice = rng.choice(options)
        if choice == "(":
            open_count += 1
        elif choice == ")":
            open_count -= 1
        chars.append(choice)
    return "".join(chars)
