# Layout-rebuild status report generator

Generates the HTML for the published status Artifact
(https://claude.ai/code/artifact/eb56f105-7ace-427b-817b-ef9011f02fc7):
a build-over-time overlap trend, the current best algorithm, before/after
renders, and every attempt tried.

## Regenerate

```
micromamba run -n py3 python docs/status/generate.py
# writes docs/status/status.html
```

Then republish `status.html` to the Artifact (same file path keeps the URL).

## Updating with a new build

When a build changes the hard-set clean rate / overlap count, append a row to
`BUILDS` in `generate.py`, rerun, republish, and update the progression table
in `docs/layout_algorithm.md`. `_renders/` (rendered PNGs) is a scratch dir.
