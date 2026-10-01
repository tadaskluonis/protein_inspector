"""Take the picture, without a human clicking Save.

A bundle is for a reader, but an agent writing a report wants a still image it
can put in the report, and until now the only way to get one was to open the
page and press two buttons. This drives the bundle's OWN export - the same
`saveImage` the Save panel calls, which forces transparency and renders to an
offscreen canvas at the requested dpi - so the PNG an agent gets is the PNG a
reader would have got, not a second implementation that drifts from it.

Needs a browser, because the picture is drawn by WebGL and there is nothing to
screenshot without one:

    pip install "protein-inspector[capture]"
    playwright install chromium

Already have a Chrome or Chromium? Point `chrome=` at it and skip the download.
"""
from __future__ import annotations

import base64
import struct
from pathlib import Path
from typing import Any

DEFAULT_DPI = 300

# Headless Chromium has no GPU. Without these the canvas comes back blank and
# the export is a transparent rectangle, which looks like a bug in the bundle.
_GL_ARGS = ["--enable-unsafe-swiftshader", "--use-gl=angle",
            "--use-angle=swiftshader"]

# The viewer registers itself on the page as window.py2dmol_viewers, so this
# works on any py2Dmol bundle, not only ones this package wrote. _triggerDownload
# is replaced rather than the download intercepted at the browser level: no
# download permission, no temp directory, no filename to guess, and it works the
# same inside a sandboxed frame.
_EXPORT_JS = r"""
async (dpi) => {
  const reg = window.py2dmol_viewers;
  if (!reg) throw new Error(
    'no py2Dmol viewer on this page - is it an inspection bundle?');
  const wanted = (window.viewerConfig || {}).viewer_id;
  const keys = Object.keys(reg);
  const entry = (wanted && reg[wanted]) || (keys.length === 1 ? reg[keys[0]] : null);
  if (!entry) throw new Error(
    'this page has ' + keys.length + ' viewers and none is named by '
    + 'viewerConfig, so there is no single picture to export');
  const r = entry.renderer;
  if (!r || typeof r.saveImage !== 'function') throw new Error(
    'this viewer build has no saveImage - nothing here can export it');

  let blob = null, failed = null;
  const original = r._triggerDownload;
  r._triggerDownload = (b) => { blob = b; };
  try {
    r.saveImage({ format: 'png', dpi: dpi });
    // saveImage finishes inside canvas.toBlob, which is asynchronous.
    const deadline = Date.now() + 30000;
    while (!blob && !failed && Date.now() < deadline) {
      await new Promise((res) => setTimeout(res, 50));
    }
  } catch (e) {
    failed = String(e && e.message || e);
  } finally {
    r._triggerDownload = original;
  }
  if (failed) throw new Error('the viewer refused to export: ' + failed);
  if (!blob) throw new Error(
    'the viewer never produced a PNG - it was still drawing, or the canvas'
    + ' is not ready');

  const b64 = await new Promise((res, rej) => {
    const fr = new FileReader();
    fr.onload = () => res(String(fr.result).split(',')[1]);
    fr.onerror = () => rej(new Error('could not read the exported blob'));
    fr.readAsDataURL(blob);
  });
  return { b64: b64, canvas_css_px: r.canvas ? r.canvas.width : null };
}
"""


def _png_shape(data: bytes) -> dict[str, Any]:
    """Width, height and whether there is an alpha channel, from the IHDR.

    Read from the bytes rather than reported by the browser, because the point
    of the export is that it is transparent and a caller should be able to see
    that it is without installing an image library.
    """
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        return {"png": False}
    width, height, _depth, colour = struct.unpack(">IIBB", data[16:26])
    return {"png": True, "width": width, "height": height,
            # 6 = truecolour with alpha, 4 = greyscale with alpha
            "has_alpha": colour in (4, 6), "colour_type": colour}


def capture_bundle(bundle: str, out_png: str, dpi: int = DEFAULT_DPI,
                   chrome: str | None = None, settle_ms: int = 4000,
                   timeout_s: int = 120, width: int = 1600,
                   height: int = 1000) -> dict[str, Any]:
    """Write `bundle`'s picture to `out_png`. Returns what was written.

    `dpi` is the export resolution the Save panel offers; 300 gives roughly a
    3500 px wide image from a 1600 px window. The background is transparent,
    which is what the bundle's own Save does and what a figure wants.

    `settle_ms` is how long to let the viewer draw before exporting. A big
    structure or a morph needs longer than a small one; if the image comes back
    empty, raise it before suspecting anything else.
    """
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - depends on the extra
        raise ImportError(
            "taking a picture needs a browser driver, which is an optional"
            ' extra: pip install "protein-inspector[capture]" and then'
            " playwright install chromium"
        ) from exc

    source = Path(bundle).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    target = Path(out_png).expanduser()

    launch: dict[str, Any] = {"args": list(_GL_ARGS)}
    if chrome:
        launch["executable_path"] = str(Path(chrome).expanduser())

    with sync_playwright() as pw:
        browser = pw.chromium.launch(**launch)
        try:
            page = browser.new_page(viewport={"width": width, "height": height})
            problems: list[str] = []
            page.on("pageerror", lambda e: problems.append(str(e)))
            page.goto(source.as_uri(), wait_until="load",
                      timeout=timeout_s * 1000)
            page.wait_for_timeout(settle_ms)
            result = page.evaluate(_EXPORT_JS, dpi)
        finally:
            browser.close()

    data = base64.b64decode(result["b64"])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)

    report = {"path": str(target.resolve()), "bytes": len(data), "dpi": dpi,
              "source": str(source)}
    report.update(_png_shape(data))
    if problems:
        report["page_errors"] = problems
    return report
