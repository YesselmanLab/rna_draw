"""The `constructive` layout engine: a bottom-up, checker-verified layout.

Every layout is checked against the frozen checker (`rna_draw.overlap`)
before it is returned, never just hoped to be clean. See
`.claude/plans/constructive-engine-runbook.md` for the full staged plan;
this package currently implements M1 (a single multiloop of hairpins) --
see `rna_draw.layout.constructive.engine` for scope.
"""

from __future__ import annotations

from .engine import ConstructiveEngine

__all__ = ["ConstructiveEngine"]
