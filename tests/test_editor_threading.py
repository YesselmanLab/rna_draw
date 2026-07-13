"""Tests for the background (off-UI-thread) layout in `MainWindow`.

`load_ss(threaded=True)` runs `EditorModel.from_ss` on a `QThread` so a large
structure never freezes the window. These tests cover: the worker's pure
`from_ss` call + emitted model; the main-thread apply-slot adopting the model;
a live threaded `load_ss` that eventually applies (event loop pumped); the
Draw/Open re-entry lock while a load is in flight; and that the close guard
neutralizes a late finish. All run headless under ``offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.main_window import _LayoutWorker  # noqa: E402
from rna_draw.gui.model import DEMO_SS, EditorModel  # noqa: E402

_SS = "((((....))))((((....))))"


@pytest.fixture
def win():
    _, w = build_app(ss=DEMO_SS, seq=None)
    yield w
    w.close()


def test_worker_run_emits_model():
    """The worker builds the model off-thread and emits it via `finished`."""
    worker = _LayoutWorker(_SS, seq=None)
    received: list[EditorModel] = []
    worker.finished.connect(received.append)
    worker.run()  # drive directly (pure, no Qt objects built)
    assert len(received) == 1
    model = received[0]
    assert isinstance(model, EditorModel)
    assert len(model.scene()["nucleotides"]) == len(_SS)


def test_worker_run_emits_failure_when_from_ss_raises(monkeypatch):
    """When `from_ss` raises, the worker emits `failed`, not `finished`."""

    def _boom(*_a, **_k):
        raise ValueError("bad structure")

    monkeypatch.setattr(EditorModel, "from_ss", classmethod(lambda cls, *a, **k: _boom()))
    worker = _LayoutWorker("....", seq=None)
    errors: list[str] = []
    models: list[EditorModel] = []
    worker.failed.connect(errors.append)
    worker.finished.connect(models.append)
    worker.run()
    assert models == []
    assert errors == ["bad structure"]


def test_apply_slot_adopts_model(win):
    """`_on_layout_done` adopts a worker-produced model into the view."""
    model = EditorModel.from_ss(_SS)
    win._on_layout_done(model, _SS, None)
    assert win._model is model
    assert len(win._view.scene().items()) > 0
    assert win._engine_label.text().startswith(f"{len(_SS)} nt")
    assert win._loading is False


def test_threaded_load_ss_eventually_applies(win):
    """A live threaded `load_ss` applies the model once the thread finishes."""
    app = QtCore.QCoreApplication.instance()
    win.load_ss(_SS, None, threaded=True)
    assert win._loading is True  # in flight
    assert win._draw_btn.isEnabled() is False  # re-entry locked

    deadline = QtCore.QDeadlineTimer(10_000)
    while win._loading and not deadline.hasExpired():
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 50)

    assert win._loading is False
    assert win._model is not None
    assert len(win._model.scene()["nucleotides"]) == len(_SS)
    assert win._draw_btn.isEnabled() is True  # re-enabled


def test_reentrant_threaded_load_is_ignored(win):
    """A second threaded load while one is in flight is ignored, not raced."""
    app = QtCore.QCoreApplication.instance()
    win.load_ss(_SS, None, threaded=True)
    first_thread = win._layout_thread
    win.load_ss("(((...)))", None, threaded=True)  # should be ignored
    assert win._layout_thread is first_thread

    deadline = QtCore.QDeadlineTimer(10_000)
    while win._loading and not deadline.hasExpired():
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 50)
    # the first (ignored-second) load applied its own structure
    assert len(win._model.scene()["nucleotides"]) == len(_SS)


def test_close_guard_neutralizes_late_finish(win):
    """After close, the finished-slot is a no-op (no touching a dead view)."""
    win._closing = True
    model = EditorModel.from_ss(_SS)
    win._on_layout_done(model, _SS, None)
    assert win._model is not model  # not adopted
