"""The picture an agent takes, and the two things that can go wrong cheaply."""
from __future__ import annotations

import struct

import pytest

from protein_inspector import capture


def _header(width: int, height: int, colour: int) -> bytes:
    """The first 26 bytes of a PNG: signature, chunk length, IHDR, fields."""
    return (b"\x89PNG\r\n\x1a\n" + struct.pack(">I", 13) + b"IHDR"
            + struct.pack(">IIBB", width, height, 8, colour))


def test_the_report_reads_the_size_and_alpha_out_of_the_png_itself():
    # Not reported by the browser: the whole point of the export is that it is
    # transparent, and a caller should be able to confirm that without an
    # image library.
    shape = capture._png_shape(_header(3809, 2369, 6))
    assert shape == {"png": True, "width": 3809, "height": 2369,
                     "has_alpha": True, "colour_type": 6}


def test_an_opaque_export_is_reported_as_opaque():
    # colour type 2 is truecolour with no alpha - what the provenance thumbnail
    # produces, and what a figure does NOT want.
    assert capture._png_shape(_header(900, 900, 2))["has_alpha"] is False


def test_something_that_is_not_a_png_is_not_described_as_one():
    assert capture._png_shape(b"<!doctype html><html>") == {"png": False}


def test_a_missing_bundle_is_refused_before_a_browser_is_started():
    pytest.importorskip("playwright", reason="capture is an optional extra")
    with pytest.raises(FileNotFoundError):
        capture.capture_bundle("no_such_bundle.html", "out.png")


def test_asking_for_a_picture_without_the_extra_says_how_to_get_it(monkeypatch):
    # The failure a user meets if they pip install the package and not the
    # extra. An ImportError that only says "playwright" sends them to the wrong
    # place; it has to name the extra and the browser step.
    import builtins

    real_import = builtins.__import__

    def refuse(name, *args, **kwargs):
        if name.startswith("playwright"):
            raise ImportError("No module named 'playwright'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", refuse)
    with pytest.raises(ImportError) as caught:
        capture.capture_bundle("anything.html", "out.png")
    message = str(caught.value)
    assert "protein-inspector[capture]" in message
    assert "playwright install chromium" in message
