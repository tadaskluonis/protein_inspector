"""Render immutable local mmCIF plus a Protein Inspector inspection manifest.

The renderer deliberately has no remote-fetch API.  It accepts bytes already on
disk, verifies the caller-provided digest, and emits a self-contained inspection
package.  py2Dmol supplies the interactive molecular canvas; the extra manifest
and projected snapshots make the scientific annotations independently auditable.
"""

from __future__ import annotations

import hashlib
import html
import json
import math
import re
import struct
import zlib
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import MMCIFParser

import py2Dmol


INSPECTOR_VERSION = "protein-inspector-1.14"
MANIFEST_SCHEMA = "protein-inspector-manifest-1"
UPSTREAM_REVISION = "78c2d489d0b5c5d19accd9eeeef878c2868f5271"
UPSTREAM_REPOSITORY = "https://github.com/sokrypton/py2Dmol"
# RETAINED IN EVERY EXPORT, because that is the whole of what the licence asks
# and minification removes it. Each bundle carries ~1 MB of py2Dmol, and the
# bundle is the thing that travels -- emailed, dropped in Slack, opened on a
# machine that has none of this checked out. A notice that lives only in the
# repository is a notice that does not reach the person holding the file.
UPSTREAM_LICENSE = (
    '"THE BEER-WARE LICENSE" (Revision 42): <so3@mit.edu> wrote this file. '
    "As long as you retain this notice you can do whatever you want with this "
    "stuff. If we meet some day, and you think this stuff is worth it, you can "
    "buy me a beer in return. Sergey Ovchinnikov"
)
INSPECTOR_REPOSITORY = "https://github.com/tadaskluonis/protein_inspector"
# The additions are on the SAME terms as the viewer, so a reader who is allowed
# to reuse py2Dmol is allowed to reuse the whole bundle rather than having to
# work out where one ends and the other begins.
INSPECTOR_LICENSE = (
    '"THE BEER-WARE LICENSE" (Revision 42): Tadas Kluonis '
    "<https://github.com/tadaskluonis> wrote the Protein Inspector "
    "additions. As long as you retain this notice you can do whatever you want "
    "with this stuff. If we meet some day, and you think this stuff is worth "
    "it, you can buy me a beer in return. Tadas Kluonis"
)

_ALLOWED_LAYER_KINDS = {
    "prediction_confidence",
    "domain",
    "topology",
    "missing_residue",
    "mutation",
    "noncanonical_residue",
    "modeled_region",
    "ptm",
    "membrane_annotation",
    "membrane_geometry",
    "partner_contact",
    "ligand",
    "ligand_contact",
    "accessibility",
    "candidate_hotspot",
    "prohibited_surface",
    "native_interface",
    "clash",
    "structural_uncertainty",
    "experimental_overlay",
    "custom",
}

_HIGHLIGHT_STYLES = {"color", "halo"}
_LAYER_COLORS = ("#dc2626", "#2563eb", "#059669", "#9333ea", "#d97706", "#0891b2")
_HEX_COLOR = re.compile(r"^#[0-9a-fA-F]{6}$")

# The neutral backdrop every bundle paints before its layers. Not a py2Dmol
# colour mode -- see _color_annotations for why a mode cannot do this.
BASE_COLOR = "#b9bfc7"

SIZE_TARGET_BYTES = 20 * 1024 * 1024
"""Keep a bundle emailable. Advisory: exceeding it warns, it never fails."""

_UNSAFE_BLOCK = re.compile(r"(?is)<\s*(script|style|iframe|object|embed)\b[^>]*>.*?<\s*/\s*\1\s*>")
_UNSAFE_TAG = re.compile(r"(?is)<\s*/?\s*(script|style|iframe|object|embed|link|meta|base|form)\b[^>]*>")
_UNSAFE_ATTR = re.compile(r"""(?is)\s(on\w+|srcdoc)\s*=\s*("[^"]*"|'[^']*'|[^\s>]+)""")
# data:image/... is allowed so a figure can be INLINED rather than linked --
# a bundle is one file, so a linked figure is a figure that does not arrive.
# Every other data: payload (notably data:text/html) is blocked along with the
# script pseudo-protocols.
_UNSAFE_URL = re.compile(
    r"""(?is)(href|src)\s*=\s*(["']?)\s*(?:javascript:|vbscript:|data:(?!image/(?:png|jpeg|gif|webp|svg\+xml);))""")


def _clean_about(body: str) -> str:
    """Accept freeform HTML for an About tab, minus anything that executes.

    The body is agent-authored and often quotes fetched text, so it is treated
    as untrusted: script/iframe/object/embed/link/meta/base/form elements, all
    on* handlers, and javascript:/data:/vbscript: URLs are removed. Everything
    else -- headings, tables, lists, links, code -- passes through, because the
    point of this field is not to be constrained.

    A body with no tags at all is taken as plain text and blank-line-separated
    blocks become paragraphs, so the simple case needs no markup.
    """
    text = str(body)
    if "<" not in text:
        blocks = [html.escape(b.strip()) for b in re.split(r"\n\s*\n", text) if b.strip()]
        return "".join(f"<p>{b}</p>" for b in blocks)
    text = _UNSAFE_BLOCK.sub("", text)   # drop the element AND its contents
    text = _UNSAFE_TAG.sub("", text)     # then any unpaired survivor
    text = _UNSAFE_ATTR.sub("", text)
    text = _UNSAFE_URL.sub(r"\1=\2#blocked:", text)
    return text


def _about_tabs(value: Any) -> list[dict[str, str]]:
    """Normalise the manifest's `about` into [{title, body_html}, ...].

    Accepts a string (one tab titled "About"), a {title: body} mapping, or a
    list of {"title", "body"} dicts. Several tabs are deliberately easy to
    produce: a bundle is ONE file, so extra context belongs in another tab
    rather than another file next to it.
    """
    if not value:
        return []
    items: list[tuple[str, str]]
    if isinstance(value, str):
        items = [("About", value)]
    elif isinstance(value, dict):
        items = [(str(k), str(v)) for k, v in value.items()]
    elif isinstance(value, list):
        items = []
        for entry in value:
            if not isinstance(entry, dict) or not entry.get("body"):
                raise ValueError("each about entry needs a 'body' (and optionally a 'title')")
            items.append((str(entry.get("title") or "About"), str(entry["body"])))
    else:
        raise ValueError("about must be a string, a {title: body} mapping, or a list of {title, body}")
    return [{"title": t, "body_html": _clean_about(b)} for t, b in items]



def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe_output(root: Path, name: str) -> Path:
    if not name or any(part in {"", ".", ".."} for part in Path(name).parts):
        raise ValueError("output_name must be a simple relative name")
    path = (root / name).resolve()
    if root.resolve() not in path.parents:
        raise ValueError("output_name escapes the output directory")
    return path


def _validate_manifest(value: dict[str, Any]) -> dict[str, Any]:
    if value.get("schema_version") != MANIFEST_SCHEMA:
        raise ValueError(f"inspection manifest must use {MANIFEST_SCHEMA}")
    annotations = value.get("annotations")
    if not isinstance(annotations, list):
        raise ValueError("inspection manifest annotations must be a list")
    ids: set[str] = set()
    for item in annotations:
        required = {"annotation_id", "kind", "label", "resolved", "evidence_ids", "method"}
        missing = required - set(item)
        if missing:
            raise ValueError(f"annotation is missing fields: {sorted(missing)}")
        if not all(isinstance(item[field], str) and item[field] for field in ("annotation_id", "kind", "label", "method")):
            raise ValueError("annotation_id, kind, label, and method must be non-empty strings")
        if not isinstance(item["resolved"], bool) or not isinstance(item["evidence_ids"], list) or not all(isinstance(value, str) for value in item["evidence_ids"]):
            raise ValueError("resolved must be a boolean and evidence_ids must be a list of strings")
        if item["kind"] not in _ALLOWED_LAYER_KINDS:
            raise ValueError(f"unsupported annotation kind: {item['kind']}")
        if item["annotation_id"] in ids:
            raise ValueError(f"duplicate annotation ID: {item['annotation_id']}")
        ids.add(item["annotation_id"])
        residue = item.get("residue")
        if item["resolved"]:
            required_residue_fields = {"component_id", "canonical_position", "chain_id", "author_residue_number"}
            if not isinstance(residue, dict) or any(residue.get(field) is None for field in required_residue_fields):
                raise ValueError(f"resolved annotation {item['annotation_id']} lacks a complete residue reference")
        for field in ("layer_id", "layer_label"):
            if field in item and (not isinstance(item[field], str) or not item[field].strip()):
                raise ValueError(f"annotation {item['annotation_id']} has an invalid {field}")
        if "layer_visible" in item and not isinstance(item["layer_visible"], bool):
            raise ValueError(f"annotation {item['annotation_id']} layer_visible must be a boolean")
        if "color" in item and (not isinstance(item["color"], str) or not _HEX_COLOR.fullmatch(item["color"])):
            raise ValueError(f"annotation {item['annotation_id']} color must be a #RRGGBB value")
        if item["kind"] == "partner_contact" and not item.get("partner_id"):
            raise ValueError("partner_contact annotations require partner_id")
    # ONE visual channel, never two. Painting the residues and haloing the same
    # residues encodes the same fact twice and makes both harder to read, so
    # `highlight` chooses between them for the whole bundle.
    if value.get("highlight", "color") not in _HIGHLIGHT_STYLES:
        raise ValueError(f"highlight must be one of {sorted(_HIGHLIGHT_STYLES)}")
    _about_tabs(value.get("about"))          # validate early, render later
    base_color = value.get("base_color", BASE_COLOR)
    if not isinstance(base_color, str) or not _HEX_COLOR.fullmatch(base_color):
        raise ValueError("base_color must be a #RRGGBB value")
    # base_mode is NOT validated against py2Dmol's list here: the page owns that
    # list (getAllValidColorModes registers custom modes such as 'ss' at load)
    # and the control is rebuilt from it at runtime.
    if not isinstance(value.get("base_mode", "chain"), str):
        raise ValueError("base_mode must be a string")
    return value


def _layer_key(annotation: dict[str, Any]) -> str:
    if annotation.get("layer_id"):
        return annotation["layer_id"]
    kind = annotation["kind"]
    if kind == "partner_contact":
        return f"partner_contact:{annotation['partner_id']}"
    if kind in {"ligand", "ligand_contact"} and annotation.get("component_id"):
        return f"{kind}:{annotation['component_id']}"
    return kind


def _layers(annotations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for item in annotations:
        key = _layer_key(item)
        layer = grouped.setdefault(key, {"annotation_ids": [], "label": item.get("layer_label"), "color": item.get("color"),
                                         "visible": bool(item.get("layer_visible", True))})
        if item.get("layer_label") and layer["label"] != item["layer_label"]:
            raise ValueError(f"layer {key!r} has conflicting labels")
        if item.get("color") and layer["color"] != item["color"]:
            raise ValueError(f"layer {key!r} has conflicting colors")
        if "layer_visible" in item and layer["visible"] != bool(item["layer_visible"]):
            raise ValueError(f"layer {key!r} disagrees about whether it starts visible")
        layer["annotation_ids"].append(item["annotation_id"])
    return [
        {
            "layer_id": key,
            "label": layer["label"] or key.replace("_", " ").replace(":", " — "),
            "color": layer["color"] or _LAYER_COLORS[index % len(_LAYER_COLORS)],
            "annotation_ids": layer["annotation_ids"],
            "visible": layer["visible"],
        }
        for index, (key, layer) in enumerate(grouped.items())
    ]


def _color_annotations(viewer: Any, annotations: list[dict[str, Any]], layers: list[dict[str, Any]],
                       base_color: str = BASE_COLOR, paint_base: bool = True,
                       mirror_chains: list[str] | None = None) -> None:
    """Paint a neutral base over the whole chain, then each layer on top.

    `mirror_chains` repeats every layer onto chains that are COPIES of the
    annotated one -- the other protomers of a homo-oligomer. An annotation
    names one chain because that is where the analysis ran, but a trimer's
    three subunits are one target presenting the same site three times, and
    painting one of them makes the reader hunt for the other two. A residue
    missing from a copy is skipped rather than raised on: a crystal often
    models a loop in one subunit and not the next, and that is not an error
    in the annotation.

    The base is painted EXPLICITLY rather than left to the viewer's colour
    mode. py2Dmol's mode list is {auto, chain, rainbow, plddt, deepmind,
    entropy, object, hydrophobicity, ss} — there is no flat-grey member, and an
    unrecognised name is not an error: `ui.js` silently falls back to `auto`,
    which on a single-chain object is rainbow. A rainbow backdrop competes with
    every layer colour, which is the one channel here that carries meaning.

    Painting the base also fixes the all-layers-off case. The per-position map
    is dropped entirely when nothing is painted (`object.color = null`), so
    unticking every layer used to reveal the rainbow rather than a neutral
    structure.
    """
    frame = viewer.objects[-1]["frames"][0]
    chains = frame.get("chains") or []
    positions: dict[tuple[str, int], list[int]] = {}
    for index, (chain, number) in enumerate(zip(chains, frame.get("residue_numbers") or [])):
        positions.setdefault((str(chain), int(number)), []).append(index)
    # ONLY IN CUSTOM MODE. py2Dmol's contract is that an explicit per-position
    # colour beats the colour mode and "the mode only decides the ones nobody
    # spoke for" (src/parts/embed.js). Seeding every position therefore speaks
    # for all of them and silently disables rainbow/plddt/chain/ss entirely --
    # which is what shipped in 1.7. A flat base is one CHOICE among the modes,
    # not a floor under them.
    if paint_base and chains:
        viewer.set_color(base_color, position=list(range(len(chains))))
    by_id = {item["annotation_id"]: item for item in annotations}
    for layer in layers:
        indices: list[int] = []
        for annotation_id in layer["annotation_ids"]:
            annotation = by_id[annotation_id]
            if not annotation["resolved"]:
                continue
            residue = annotation["residue"]
            matched = positions.get((str(residue["chain_id"]), int(residue["author_residue_number"])), [])
            if not matched:
                raise ValueError(f"resolved annotation {annotation_id} does not match a residue in the mmCIF")
            indices.extend(matched)
            for copy in mirror_chains or ():
                indices.extend(positions.get((str(copy), int(residue["author_residue_number"])), []))
        if indices:
            viewer.set_color(layer["color"], position=indices)


def _partner_copies(frame: dict[str, Any], partner_chains: set[str],
                    mode: bool | str = "auto") -> list[str]:
    """Which partner chains are copies of the annotated target, by sequence.

    A homo-oligomer is one target wearing several chain letters, so a layer
    computed on chain A is a fact about chains B and C too. A receptor, a Fab
    or a ligand-bearing chain is not, and painting an epitope onto it would be
    a claim nobody made. The two cases are told apart here rather than asked
    of the caller: a partner counts as a copy when it shares most of the
    target's author residue numbers AND agrees with it on the residue name at
    essentially all of them.

    `mode="auto"` applies that test, `True` takes every partner chain
    whatever it is (for numbering you have made equivalent yourself), and
    `False` returns nothing. Chains are reported under the names they carry
    in the rendered object, which is what `_color_annotations` addresses and
    what a collision may have renamed.
    """
    if mode is False or not partner_chains:
        return []
    chains = [str(c) for c in (frame.get("chains") or [])]
    numbers = [int(n) for n in (frame.get("residue_numbers") or [])]
    if mode is True:
        return sorted({c for c in chains if c in partner_chains})
    names = list(frame.get("position_names") or [])
    if len(names) != len(chains):
        # Nothing to compare: a copy cannot be told from a binding partner,
        # and guessing here would paint a receptor with the target's epitope.
        return []
    residues: dict[str, dict[int, str]] = {}
    for chain, number, name in zip(chains, numbers, names):
        residues.setdefault(chain, {})[number] = str(name)
    targets = [c for c in residues if c not in partner_chains]
    copies = []
    for partner in sorted(c for c in residues if c in partner_chains):
        other = residues[partner]
        for target in targets:
            ref = residues[target]
            shared = set(ref) & set(other)
            if len(shared) < 10 or len(shared) < 0.8 * len(ref):
                continue
            agree = sum(1 for n in shared if ref[n] == other[n])
            if agree >= 0.95 * len(shared):
                copies.append(partner)
                break
    return copies


def _trace(cif_path: Path) -> tuple[np.ndarray, list[str], list[int]]:
    structure = MMCIFParser(QUIET=True).get_structure(cif_path.stem, str(cif_path))
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError("inspection requires exactly one prepared coordinate model")
    coords: list[list[float]] = []
    chains: list[str] = []
    numbers: list[int] = []
    for chain in models[0]:
        for residue in chain:
            atom = residue.child_dict.get("CA") or residue.child_dict.get("C4'")
            if atom is None or residue.id[0] != " ":
                continue
            coords.append([float(value) for value in atom.coord])
            chains.append(chain.id)
            numbers.append(int(residue.id[1]))
    if len(coords) < 2:
        raise ValueError("prepared coordinate artifact has no renderable polymer trace")
    return np.asarray(coords, dtype=float), chains, numbers


def _project(coords: np.ndarray, size: int = 900, margin: int = 40) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    xy = centered @ vh[:2].T
    span = np.maximum(np.ptp(xy, axis=0), 1e-6)
    scale = min((size - margin * 2) / span[0], (size - margin * 2) / span[1])
    xy = xy * scale
    xy[:, 0] += size / 2
    xy[:, 1] = size / 2 - xy[:, 1]
    return xy


def _snapshot_svg(path: Path, xy: np.ndarray, chains: list[str], numbers: list[int]) -> None:
    palette = ["#2563eb", "#dc2626", "#059669", "#9333ea", "#d97706", "#0891b2"]
    chain_colors = {chain: palette[index % len(palette)] for index, chain in enumerate(dict.fromkeys(chains))}
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="900" viewBox="0 0 900 900">',
        '<rect width="900" height="900" fill="white"/>',
    ]
    for index in range(1, len(xy)):
        if chains[index] != chains[index - 1]:
            continue
        parts.append(f'<line x1="{xy[index-1,0]:.2f}" y1="{xy[index-1,1]:.2f}" x2="{xy[index,0]:.2f}" y2="{xy[index,1]:.2f}" stroke="{chain_colors[chains[index]]}" stroke-width="5" stroke-linecap="round"/>')
    for index, (point, chain, number) in enumerate(zip(xy, chains, numbers)):
        if index % 10 == 0:
            parts.append(f'<circle cx="{point[0]:.2f}" cy="{point[1]:.2f}" r="4" fill="{chain_colors[chain]}"><title>{chain}:{number}</title></circle>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def _png_chunk(kind: bytes, payload: bytes) -> bytes:
    return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)


def _snapshot_png(path: Path, xy: np.ndarray, chains: list[str], size: int = 900) -> None:
    rgb = bytearray([255] * (size * size * 3))
    palette = [(37, 99, 235), (220, 38, 38), (5, 150, 105), (147, 51, 234), (217, 119, 6), (8, 145, 178)]
    chain_colors = {chain: palette[index % len(palette)] for index, chain in enumerate(dict.fromkeys(chains))}

    def draw(x0: float, y0: float, x1: float, y1: float, color: tuple[int, int, int]) -> None:
        steps = max(1, int(max(abs(x1 - x0), abs(y1 - y0))))
        for step in range(steps + 1):
            x = int(round(x0 + (x1 - x0) * step / steps))
            y = int(round(y0 + (y1 - y0) * step / steps))
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    px, py = x + dx, y + dy
                    if 0 <= px < size and 0 <= py < size:
                        offset = (py * size + px) * 3
                        rgb[offset : offset + 3] = bytes(color)

    for index in range(1, len(xy)):
        if chains[index] == chains[index - 1]:
            draw(*xy[index - 1], *xy[index], chain_colors[chains[index]])
    scanlines = b"".join(b"\x00" + bytes(rgb[row * size * 3 : (row + 1) * size * 3]) for row in range(size))
    png = b"\x89PNG\r\n\x1a\n" + _png_chunk(b"IHDR", struct.pack(">IIBBBBB", size, size, 8, 2, 0, 0, 0)) + _png_chunk(b"IDAT", zlib.compress(scanlines, 9)) + _png_chunk(b"IEND", b"")
    path.write_bytes(png)


def _html(viewer_html: str, state: dict[str, Any]) -> str:
    payload = json.dumps(state, separators=(",", ":")).replace("</", "<\\/")
    # THE TWO ROWS THE HORIZONTAL SEAM SPLITS, seeded from what the caller
    # asked for rather than from a constant: --pinsp-view starts at the
    # `display` height, so overriding the viewer's own inline canvas height
    # (see #canvasContainer below) changes nothing until somebody drags.
    size = ((state.get("viewer") or {}).get("config") or {}).get("display") or {}
    view_height = int((size.get("size") or (720, 720))[1] or 720)
    seq_height = 152
    css = (
        "<style>"
        f".protein-inspector{{--pinsp-view:{view_height}px;--pinsp-seq:{seq_height}px}}"
        "html,body{margin:0;background:#fff}body{padding:14px;box-sizing:border-box}"
        ".protein-inspector{font:14px/1.45 system-ui,sans-serif;display:grid;"
        "grid-template-columns:minmax(240px,1fr) 11px var(--pinsp-panel,340px);"
        # ROW ONE IS A HEIGHT, NOT `auto`.
        #
        # It was `auto`, and the Layers panel lives in it: on a bundle with a
        # long residue list the panel is the tallest thing in the row, pinned
        # to its own max-height, so the row's height was the PANEL's. Shrinking
        # the canvas then moved nothing -- the seam stayed exactly where it
        # was while the picture quietly resized inside it, which is precisely
        # what the reader reported. Naming the row's height makes the seam the
        # authority over it, the canvas fills what is left of the stage, and
        # the panel takes the row and scrolls inside it.
        f"grid-template-rows:var(--pinsp-view,{view_height}px) 11px "
        f"var(--pinsp-seq,{seq_height}px);gap:0;"
        "color:#0f172a;align-items:stretch}"
        # PLACED, NOT AUTO-FLOWED. With a third zone the grid has rows that are
        # not all `auto`, and auto-placement hands each item the next cell --
        # so the narrow layout, where the vertical seam is display:none and
        # therefore placed nowhere, slid the Layers panel into the 11px seam
        # row. Every zone names its own cell instead.
        ".pinsp-stage{grid-column:1;grid-row:1}"
        ".pinsp-split{grid-column:2;grid-row:1}"
        ".pinsp-panel{grid-column:3;grid-row:1}"
        ".pinsp-hsplit{grid-column:1/-1;grid-row:2}"
        ".pinsp-seq{grid-column:1/-1;grid-row:3}"
        # THE TWO PANES TOUCH, AND THE SEAM IS THE CONTROL. A picture that can
        # be sized independently of the panel beside it can always be sized
        # into overlapping it, or into leaving a band of nothing between them.
        # One draggable boundary removes both states: the grid owns the split,
        # the stage takes whatever is left of it, and the canvas follows.
        ".pinsp-split{cursor:col-resize;position:relative;touch-action:none;"
        "align-self:stretch;min-height:120px}"
        ".pinsp-split::before{content:'';position:absolute;top:2px;bottom:2px;left:4px;"
        "width:3px;border-radius:2px;background:#e2e8f0}"
        ".pinsp-split:hover::before,.pinsp-split[data-drag='1']::before{background:#94a3b8}"
        ".pinsp-split:focus-visible{outline:2px solid #2563eb;outline-offset:-2px}"
        # THE TOOLS BELONG INSIDE THE PICTURE. py2Dmol lays the viewer out as
        # a flex row -- canvas, then a 340px column of Orient/Focus/Rotate/
        # Style/Clip/Capture -- and pins .py2dmol-viewer-instance to 948px to
        # fit both. Floating that column over the top-right of the canvas
        # frees the width, puts the controls on the thing they act on, and
        # leaves the grid's second column for the Layers panel.
        # A COLUMN, so the picture can take "the rest of the row" rather than a
        # number. The conformation ring and the capture bar keep their natural
        # heights and #mainContainer absorbs the remainder, which is what makes
        # one custom property on the row enough to resize the canvas.
        ".pinsp-stage{position:relative;min-width:0;min-height:0;overflow:hidden;"
        "display:flex;flex-direction:column}"
        # The viewer's own root div carries an id and NO class in an export,
        # so a `.py2dmol-viewer-instance` selector silently matches nothing --
        # it is reached from JS below instead, by walking up from
        # #mainContainer. Left here for the builds that do carry the class.
        ".pinsp-stage .py2dmol-viewer-instance,.pinsp-stage #viewerWrapper"
        "{display:contents!important}"
        # STILL A REAL POSITIONED BOX. The floating Orient/Focus/Style cluster
        # is its child and a SIBLING of #canvasContainer, so collapsing this
        # detaches the tools onto whatever sits above the picture. It is now
        # also the flex child that takes the leftover height.
        ".pinsp-stage #mainContainer{display:block!important;position:relative;"
        "width:auto!important;max-width:none!important;padding:0!important;"
        "flex:1 1 auto;min-height:0}"
        # WIDTH FLOWS ONE WAY ONLY. py2Dmol's ResizeObserver answers a
        # container resize by writing the observed width back onto
        # #viewerWrapper as an inline style. Once #canvasContainer is fluid,
        # that closes a loop -- wrapper width from container, container width
        # from wrapper -- and each pass loses the container's border, so the
        # picture starts at the right size and then walks itself narrow. The
        # !important here makes the observer's write inert: the wrapper takes
        # its width from the grid column and nothing else.
        # THE DISPLAY AREA STOPS AT ITS OWN COLUMN. #canvasContainer carries
        # `resize: both` and an explicit pixel width, so without a ceiling the
        # reader can drag the picture out from under the Layers panel -- or
        # simply be handed a bundle whose canvas is wider than the frame it
        # opens in. The cap is the column, so the two never overlap.
        # THE PICTURE IS AS WIDE AS ITS COLUMN. setupViewport writes an inline
        # pixel width on #canvasContainer, which makes the display area a fixed
        # box that is either narrower than the space it has or wide enough to
        # slide under the Layers panel. Overriding it with !important (inline
        # styles lose to that, and only to that) hands the width back to the
        # grid, so the right edge lands just short of the panel at any frame
        # size. A ResizeObserver on this element re-renders, so the canvas
        # follows. Height stays the caller's, and stays draggable.
        # HEIGHT IS THE SEAM'S, the way width is the grid column's -- but taken
        # from the STAGE rather than from the custom property directly. The
        # viewer writes an inline pixel height here from setupViewport, and
        # py2Dmol's own ResizeObserver watches this element for height as well
        # as width (parts/viewport.js: it re-sizes the canvas, calls
        # _updateCanvasDimensions and renders). So overriding the inline value
        # with `100%` of a flex-sized parent -- inline styles lose to
        # !important, and only to that -- is the whole of making the picture
        # follow the horizontal seam, and it leaves room for the conformation
        # ring and the capture bar instead of overflowing past them.
        ".pinsp-stage #canvasContainer{display:block!important;width:auto!important;"
        "max-width:none!important;margin:0!important;resize:none!important;"
        "height:100%!important;position:relative}"
        ".pinsp-stage #canvasContainer .resize-handle{display:none!important}"
        ".pinsp-stage #canvasContainer canvas{max-width:100%}"
        ".pinsp-stage #rightPanelContainer{position:absolute!important;top:10px;left:10px;"
        "z-index:5;width:336px;max-width:calc(100% - 20px);max-height:calc(100% - 20px);"
        "overflow-y:auto;background:rgba(255,255,255,.94);border:1px solid #e2e8f0;"
        "border-radius:10px;padding:8px;box-shadow:0 8px 24px rgba(15,23,42,.14)}"
        # THE LAYERS PANEL STAYS BESIDE THE STRUCTURE. Ticking a layer and
        # watching the structure change is the whole interaction, and it does
        # not survive the panel being pushed below the fold. It only stacks
        # where a side-by-side would leave the canvas unusably narrow.
        ".pinsp-stage[data-morph='1'] #controlsContainer{display:none!important}"
        "@media (max-width:820px){.protein-inspector{grid-template-columns:minmax(0,1fr);"
        "grid-template-rows:var(--pinsp-view,760px) auto 11px var(--pinsp-seq,152px)}"
        ".pinsp-split{display:none}.pinsp-panel{max-height:60vh;margin-top:12px}"
        ".pinsp-stage{grid-column:1;grid-row:1}.pinsp-panel{grid-column:1;grid-row:2}"
        ".pinsp-hsplit{grid-row:3}.pinsp-seq{grid-row:4}}"
        # THE RING. One segment per conformation in a single pill, so the set
        # of states is visible at a glance and the current one is obvious
        # without reading a label.
        # max-width was `calc(100% - 120px)`, reserving room for a button that
        # no longer exists, and the ring could not wrap -- so on a narrow
        # stage the later conformations ran off the edge, were clipped by the
        # stage's overflow, and their clicks landed on the Layers header
        # behind them. A control you can see the label of but cannot press is
        # worse than one that has wrapped onto a second line.
        ".pinsp-morph{display:flex;flex-wrap:wrap;align-items:center;gap:8px;"
        "padding:0 0 8px;max-width:100%}"
        ".bm-ring{display:inline-flex;flex-wrap:wrap;align-items:center;gap:2px;padding:3px;"
        "border:1px solid #cbd5e1;border-radius:999px;background:#f8fafc;max-width:100%}"
        ".bm-btn{font:11.5px system-ui;padding:5px 13px;border:0;background:none;"
        "border-radius:999px;cursor:pointer;color:#475569;white-space:nowrap;"
        "max-width:100%;overflow:hidden;text-overflow:ellipsis}"
        ".bm-btn:hover:not([aria-pressed='true']){background:#fff;color:#0f172a}"
        ".bm-btn[aria-pressed='true']{background:#0f172a;color:#fff}"
        ".pinsp-morph[data-busy='1'] .bm-btn{cursor:default;opacity:.55}"
        ".bm-solo{border:1px solid #cbd5e1;background:#f8fafc;padding:6px 14px}"
        ".bm-solo[aria-pressed='true']{background:#0f172a;color:#fff;border-color:#0f172a}"
        ".bm-solo[disabled]{opacity:.45;cursor:default;background:#f8fafc;color:#94a3b8;"
        "border-color:#e2e8f0}"
        # A CAPTURE THE READER CAN ACTUALLY KEEP. See the script: a download
        # anchor is dropped without error in a frame that lacks
        # allow-downloads, and the viewer reports success regardless.
        ".pinsp-capture{margin:8px 0 0;padding:9px 11px;border:1px solid #cbd5e1;"
        "border-radius:8px;background:#f8fafc;font:12px system-ui}"
        ".pinsp-capture[hidden]{display:none}"
        ".pinsp-capture a{color:#1d4ed8}"
        ".pinsp-capture .bp-btn{margin:0 6px 0 0}"
        ".pinsp-capture a.bp-btn{text-decoration:none;display:inline-block}"
        # THE PREVIEW SITS ON A CHECKERBOARD. saveImage forces
        # isTransparent for the export, so the PNG has no background at
        # all -- and on white it looked exactly like a white background
        # baked in, which is the one thing the reader needs to know it
        # is not.
        ".pinsp-capture img{display:block;margin:8px 0 0;max-width:100%;max-height:260px;"
        "border:1px solid #e2e8f0;border-radius:6px;"
        "background-image:linear-gradient(45deg,#eef2f7 25%,transparent 25%,transparent 75%,#eef2f7 75%),"
        "linear-gradient(45deg,#eef2f7 25%,transparent 25%,transparent 75%,#eef2f7 75%);"
        "background-size:18px 18px;background-position:0 0,9px 9px}"
        # THE PANEL TAKES THE ROW AND SCROLLS INSIDE IT. It was `max-height:88vh`
        # plus `position:sticky`, which is how it came to DICTATE the row's
        # height on any bundle with a long residue list -- a legend taller than
        # the picture beside it, and a seam that could not move because the
        # panel was holding the row open. `min-height:0` is what lets a flex/grid
        # child be shorter than its content; the internal overflow then does the
        # work, and the sticky `.bp-head` rows keep Layers and Residues visible
        # while it scrolls.
        ".pinsp-panel{border:1px solid #e2e8f0;border-radius:10px;min-height:0;overflow:auto;"
        "display:flex;flex-direction:column;background:#fff;"
        "box-shadow:0 8px 24px rgba(15,23,42,.06)}"
        ".bp-head{position:sticky;top:0;z-index:2;background:#fff;display:flex;justify-content:space-between;"
        "align-items:center;gap:8px;padding:10px 12px 6px;border-bottom:1px solid #f1f5f9}"
        ".bp-head h2{font-size:12px;letter-spacing:.04em;text-transform:uppercase;color:#475569;margin:0}"
        ".bp-btns{display:flex;gap:5px}"
        ".bp-btn{font:11px system-ui;padding:3px 9px;border:1px solid #cbd5e1;background:#fff;"
        "border-radius:6px;cursor:pointer;color:#0f172a}"
        ".bp-btn:hover{background:#f1f5f9;border-color:#94a3b8}"
        ".bp-layers{padding:6px 12px 10px}"
        ".bp-layer{display:flex;align-items:center;gap:7px;padding:2px 0}"
        ".bp-layer label{display:flex;align-items:center;gap:7px;flex:1;min-width:0;cursor:pointer}"
        ".bp-sw{width:13px;height:13px;border-radius:3px;border:1px solid rgba(0,0,0,.15);flex:none}"
        ".bp-lab{font-size:12.5px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}"
        ".bp-n{font:11px ui-monospace,monospace;color:#64748b;flex:none}"
        ".bp-only{font:10px system-ui;padding:1px 6px;border:1px solid #e2e8f0;background:#fff;"
        "border-radius:5px;cursor:pointer;color:#64748b;flex:none;visibility:hidden}"
        ".bp-layer:hover .bp-only{visibility:visible}.bp-only:hover{border-color:#94a3b8;color:#0f172a}"
        "#pinsp-filter{margin:0 12px 8px;padding:5px 8px;font:12px system-ui;border:1px solid #cbd5e1;"
        "border-radius:6px}"
        ".bp-card{margin:0 12px 10px;border:1px solid #cbd5e1;border-radius:8px;background:#f8fafc;"
        "padding:9px 11px}"
        ".bp-card-head{display:flex;justify-content:space-between;align-items:baseline;gap:8px}"
        ".bp-card-head b{font-size:15px}"
        ".bp-addr{font:10.5px ui-monospace,monospace;color:#64748b;margin:3px 0 7px}"
        ".bp-card ul{list-style:none;margin:0;padding:0}"
        ".bp-card li{display:flex;gap:7px;align-items:flex-start;padding:3px 0;font-size:12px;"
        "border-top:1px solid #e9eef4}"
        ".bp-card li:first-child{border-top:0}"
        ".bp-meta{font-size:10.5px;color:#64748b;display:block}"
        ".bp-x{font:14px system-ui;line-height:1;border:0;background:none;cursor:pointer;color:#64748b;padding:0 2px}"
        ".bp-x:hover{color:#0f172a}"
        ".bp-rows{padding:0 12px 12px}"
        ".bp-group{margin-bottom:6px}"
        ".bp-group>summary{cursor:pointer;font-size:11.5px;color:#475569;padding:3px 0;list-style:none;"
        "display:flex;align-items:center;gap:6px}"
        ".bp-group>summary::-webkit-details-marker{display:none}"
        ".bp-group>summary::before{content:'▸';font-size:9px;color:#94a3b8}"
        ".bp-group[open]>summary::before{content:'▾'}"
        ".pinsp-residue{cursor:pointer;padding:2px 6px;border-radius:5px;font-size:12px;"
        "border-left:3px solid transparent}"
        ".pinsp-residue:hover{background:#f1f5f9}"
        ".pinsp-residue[data-on='1']{background:#fef9c3;border-left-color:#eab308}"
        ".pinsp-warning{color:#9a3412}"
        ".bp-muted{font-size:11px;color:#64748b}"
        ".pinsp-tabs{display:flex;gap:2px;border-bottom:1px solid #e2e8f0;margin:0 0 14px;font:14px system-ui,sans-serif}"
        ".pinsp-tab{font:13px system-ui;padding:8px 16px;border:0;background:none;cursor:pointer;color:#64748b;border-bottom:2px solid transparent;margin-bottom:-1px}"
        ".pinsp-tab:hover{color:#0f172a}"
        '.pinsp-tab[aria-selected="true"]{color:#0f172a;font-weight:600;border-bottom-color:#2563eb}'
        ".pinsp-about{max-width:52em;font:15px/1.6 system-ui,sans-serif;color:#0f172a;padding:0 4px 28px}"
        ".pinsp-about h2{font-size:17px;margin:26px 0 6px;border-bottom:1px solid #e2e8f0;padding-bottom:4px}"
        ".pinsp-about h3{font-size:15px;margin:20px 0 4px}"
        ".pinsp-about table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px}"
        ".pinsp-about td,.pinsp-about th{border-bottom:1px solid #e2e8f0;padding:6px 8px;text-align:left;vertical-align:top}"
        ".pinsp-about code{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:13px}"
        ".pinsp-about .sw{display:inline-block;width:14px;height:14px;border-radius:3px;border:1px solid #cbd5e1;vertical-align:-2px}"
        "#pinsp-residue-details{display:none}"
        # ---------------------------------------------------------------
        # THE SEQUENCE IS THE THIRD ZONE, AND IT IS FULL WIDTH.
        #
        # A strip belongs under BOTH the picture and the panel rather than
        # beside either: it is one long line, so a column would wrap it into
        # an unreadable block, and it is the index to both of its neighbours
        # -- the letters carry the layer colours the structure carries, and
        # clicking one is the same act as clicking the residue. Its own seam
        # sets its height, so the reader decides how much of it is on screen
        # without the picture's width moving.
        ".pinsp-hsplit{cursor:row-resize;position:relative;touch-action:none;min-width:40px}"
        ".pinsp-hsplit::before{content:'';position:absolute;left:2px;right:2px;top:4px;"
        "height:3px;border-radius:2px;background:#e2e8f0}"
        ".pinsp-hsplit:hover::before,.pinsp-hsplit[data-drag='1']::before{background:#94a3b8}"
        ".pinsp-hsplit:focus-visible{outline:2px solid #2563eb;outline-offset:-2px}"
        ".pinsp-seq{min-height:0;display:flex;flex-direction:column;overflow:hidden;"
        "border:1px solid #e2e8f0;border-radius:10px;background:#fff;"
        "box-shadow:0 8px 24px rgba(15,23,42,.06)}"
        ".pinsp-seqbar{display:flex;align-items:center;gap:9px;flex-wrap:wrap;"
        "padding:7px 11px;border-bottom:1px solid #f1f5f9}"
        ".pinsp-seqbar h2{font-size:12px;letter-spacing:.04em;text-transform:uppercase;"
        "color:#475569;margin:0;flex:none}"
        ".pinsp-seqbar #pinsp-sel{min-width:0;overflow:hidden;text-overflow:ellipsis;"
        "white-space:nowrap}"
        ".pinsp-tools{margin-left:auto;display:flex;align-items:center;gap:6px;flex:none}"
        ".pinsp-tools label{font:11px system-ui;color:#64748b}"
        "#pinsp-around-a{width:50px;font:11px system-ui;padding:3px 4px;"
        "border:1px solid #cbd5e1;border-radius:5px}"
        ".bp-btn[aria-pressed='true']{background:#0f172a;color:#fff;border-color:#0f172a}"
        ".pinsp-pair{display:inline-flex}"
        ".pinsp-pair .pinsp-sw:first-child{border-radius:6px 0 0 6px}"
        ".pinsp-pair .pinsp-sw:last-child{border-radius:0 6px 6px 0;margin-left:-1px}"
        ".bp-btn[disabled]{opacity:.45;cursor:default;background:#f8fafc;color:#94a3b8;"
        "border-color:#e2e8f0}"
        # The folded half of the Style panel (see foldStylePanel). Deliberately
        # quiet: a summary that reads as a label rather than a button, because
        # what is behind it is authoring controls and not a second panel. The
        # rows inside keep every class and id the viewer's own stylesheet and
        # parts/ui.js already target, so the skin inside the fold is unchanged.
        ".pinsp-fine{margin-top:2px;border-top:1px solid #e2e8f0;padding-top:4px}"
        ".pinsp-fine>summary{font:11px system-ui;color:#64748b;cursor:pointer;"
        "list-style:none;padding:2px 0;user-select:none}"
        ".pinsp-fine>summary::-webkit-details-marker{display:none}"
        ".pinsp-fine>summary::before{content:'\\25b8 ';color:#94a3b8}"
        ".pinsp-fine[open]>summary::before{content:'\\25be '}"
        ".pinsp-fine>summary:hover{color:#0f172a}"
        ".pinsp-seqbody{flex:1;min-height:0;overflow:auto;padding:15px 11px 10px}"
        ".pinsp-seqchain{display:flex;align-items:flex-start;gap:9px;margin:0 0 14px}"
        ".pinsp-seqlab{font:11px ui-monospace,monospace;padding:2px 7px;flex:none;"
        "border:1px solid #cbd5e1;border-radius:5px;background:#f8fafc;color:#0f172a;"
        "cursor:pointer;position:sticky;left:0}"
        ".pinsp-seqlab:hover{border-color:#94a3b8}"
        # ROW-GAP, BECAUSE THE RULER SITS ABOVE THE LETTER. Every tenth cell
        # prints its author number in the space over it; without a row gap a
        # wrapped line's numbers land on the line above's letters.
        ".pinsp-seqcells{display:flex;flex-wrap:wrap;gap:13px 1px;min-width:0}"
        ".pinsp-cell{position:relative;width:13px;height:17px;flex:none;"
        "font:11px/17px ui-monospace,monospace;text-align:center;cursor:pointer;"
        "background:#f1f5f9;color:#334155;border-radius:2px;"
        "-webkit-user-select:none;user-select:none}"
        ".pinsp-cell[data-num]::before{content:attr(data-num);position:absolute;"
        "bottom:17px;left:0;font:9px system-ui;color:#94a3b8;white-space:nowrap}"
        ".pinsp-cell[data-lig='1']{border-radius:8px}"
        ".pinsp-cell:hover,.pinsp-cell[data-hover='1']{outline:1px solid #94a3b8}"
        ".pinsp-cell[data-sel='1']{box-shadow:inset 0 0 0 2px #eab308;font-weight:700}"
        ".pinsp-seqnote{margin:0;font-size:11.5px;color:#64748b}"
        # A RULE ACROSS THE FOOT OF THE SEQUENCE SECTION, not a block under it.
        # `flex:none` so it never takes space from the strip, and one line that
        # ellipsises rather than wrapping into a second: the seam's arithmetic
        # assumes this costs a constant height, and a credit that reflowed at
        # some window widths would quietly move the strip's floor with it.
        ".pinsp-credit{margin:0;flex:none;padding:5px 11px 6px;"
        "border-top:1px solid #e2e8f0;font:11px system-ui;color:#94a3b8;"
        "white-space:nowrap;overflow:hidden;text-overflow:ellipsis}"
        ".pinsp-credit a{color:#64748b}"
        "</style>"
    )
    script = r"""<script>(function(){
var root=document.getElementById('pinsp-inspection-state');if(!root)return;
var s=JSON.parse(root.textContent);
var byId={};s.annotations.forEach(function(a){byId[a.annotation_id]=a;});
var layerOf={};s.layers.forEach(function(l){layerOf[l.layer_id]=l;});
var PLACEHOLDER='Click a residue to see every annotation on it.';
var applying=false,selectedId=null,query='';
var baseMode=s.base_mode||'chain',baseColor=s.base_color||'#b9bfc7';
function esc(t){return String(t==null?'':t).replace(/[&<>"]/g,function(c){
  return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];});}
function renderer(){var reg=window.py2dmol_viewers;if(!reg)return null;
  var id=(s.viewer&&s.viewer.config&&s.viewer.config.viewer_id)||(window.viewerConfig&&window.viewerConfig.viewer_id);
  var e=id?reg[id]:null;if(!e){var k=Object.keys(reg);if(k.length===1)e=reg[k[0]];}
  return e&&e.renderer?e.renderer:null;}
function indicesFor(res){var r=renderer(),out=[];if(!r||!res)return out;
  for(var i=0;i<(r.residueNumbers||[]).length;i++){
    if(r.chains&&r.chains[i]===res.chain_id&&r.residueNumbers[i]===res.author_residue_number)out.push(i);}
  return out;}
function sameResidue(a,b){return a&&b&&a.component_id===b.component_id&&a.copy_index===b.copy_index&&
  a.canonical_position===b.canonical_position;}
function box(id){return document.querySelector('[data-layer="'+CSS.escape(id)+'"]');}
function layerVisible(id){var b=box(id);return !b||b.checked;}
function rowMatches(a){if(!query)return true;
  return (a.label+' '+a.annotation_id+' '+(layerOf[a.layer_id]?layerOf[a.layer_id].label:'')).toLowerCase().indexOf(query)>=0;}
function syncVisibleLayers(){
  var r=renderer(),selected={},shown=0;
  for(var li=0;li<s.layers.length;li++){var l=s.layers[li];var vis=layerVisible(l.layer_id);
    for(var ai=0;ai<l.annotation_ids.length;ai++){var id=l.annotation_ids[ai];var a=byId[id];
      var row=document.querySelector('[data-annotation="'+CSS.escape(id)+'"]');
      var on=vis&&rowMatches(a);
      if(row)row.hidden=!on;
      if(on)shown++;
      if(vis&&a&&a.resolved){var ix=indicesFor(a.residue);for(var k=0;k<ix.length;k++)selected[ix[k]]=1;}}
    var g=document.querySelector('[data-group="'+CSS.escape(l.layer_id)+'"]');if(g)g.hidden=!vis;}
  var c=document.getElementById('pinsp-count');
  if(c)c.textContent=shown+' of '+s.annotations.length;
  if(!r)return false;
  if(s.highlight!=='halo'){
    var name=r.currentObjectName,obj=r.objectsData&&r.objectsData[name];
    if(obj){var off=0;if(typeof r.localRangeOf==='function'){var w=r.localRangeOf(name);
        if(w&&typeof w.off==='number')off=w.off;}
      var paint={},painted=false;
      // A flat base is one CHOICE among the colour modes, not a floor under
      // them: an explicit per-position colour beats the mode, so seeding every
      // position would disable rainbow/plddt/chain/ss outright.
      if(baseMode==='custom'){var ntot=(r.residueNumbers||[]).length;
        for(var bi=0;bi<ntot;bi++){paint[bi-off]=baseColor;painted=true;}}
      for(var li2=0;li2<s.layers.length;li2++){var l2=s.layers[li2];if(!layerVisible(l2.layer_id))continue;
        for(var aj=0;aj<l2.annotation_ids.length;aj++){var a2=byId[l2.annotation_ids[aj]];
          if(!a2||!a2.resolved||!a2.color)continue;
          var ix2=indicesFor(a2.residue);for(var m=0;m<ix2.length;m++){paint[ix2[m]-off]=a2.color;painted=true;}}}
      obj.color=painted?{type:'advanced',value:{position:paint}}:null;
      r.colorsNeedUpdate=true;r.plddtColorsNeedUpdate=true;}}
  var set=new Set();for(var key in selected)set.add(Number(key));
  applying=true;try{if(s.highlight==='halo')r.setResidueSelection(set);
    r.render('Protein Inspector inspection layer visibility');}finally{applying=false;}
  // THE STRIP IS THE SAME LEGEND. Repainted from the same visible-layer set,
  // so a layer cannot be off in the picture and lit in the sequence.
  seqPaintLayers();seqPaintSelection();
  return true;}
// THE CARD IS ABOUT A RESIDUE, NOT ABOUT AN ANNOTATION.
//
// It was reachable only from a manifest row, so it could only ever open on a
// residue somebody had already annotated -- and a reader who clicks the
// structure is usually asking about one nobody spoke for. Same card, seeded
// either way: an address, and every annotation on it, which is sometimes none.
function showCard(res,related,headline){
  res=res||{};
  var rows=document.querySelectorAll('[data-annotation]');
  for(var i=0;i<rows.length;i++)rows[i].removeAttribute('data-on');
  for(var j=0;j<related.length;j++){
    var el=document.querySelector('[data-annotation="'+CSS.escape(related[j].annotation_id)+'"]');
    if(el&&!el.hidden)el.setAttribute('data-on','1');}
  var d=document.getElementById('pinsp-residue-details');
  if(d)d.textContent=headline||PLACEHOLDER;
  var card=document.getElementById('pinsp-card');if(!card)return;
  var title=res.component_id?(res.component_id+' '+res.author_residue_number)
    :((headline||'').split(' ')[0]||res.author_residue_number);
  var items='';
  for(var k=0;k<related.length;k++){var x=related[k];var lay=layerOf[x.layer_id]||{};
    var head=lay.label||x.layer_id||x.kind;
    var meta=(x.label&&x.label!==head)?'<span class="bp-meta">'+esc(x.label)+'</span>':'';
    items+='<li><span class="bp-sw" style="background:'+esc(x.color||'#94a3b8')+'"></span><span>'+
      esc(head)+meta+'</span></li>';}
  if(!items)items='<li><span class="bp-meta">No annotation on this residue.</span></li>';
  card.innerHTML='<div class="bp-card-head"><b>'+esc(title)+'</b>'+
    '<button class="bp-x" id="pinsp-card-close" title="clear (Esc)">&times;</button></div>'+
    '<div class="bp-addr" title="position '+esc(res.canonical_position)+
    ' of this chain in the file as written">chain '+esc(res.chain_id)+
    ' &middot; residue '+esc(res.author_residue_number)+'</div><ul>'+items+'</ul>';
  card.hidden=false;
  var xb=document.getElementById('pinsp-card-close');
  if(xb)xb.addEventListener('click',function(ev){ev.stopPropagation();clearDetails();});}
function showDetails(a){
  var related=s.annotations.filter(function(x){return sameResidue(x.residue,a.residue);});
  if(!related.length)related=[a];
  showCard(a.residue,related,a.label);}
// ...and the same card from a position index, which is what a click in the
// viewer and a letter in the strip both have.
function showIndex(i){
  var hits=annotationsAt(i);
  if(hits.length){selectedId=hits[0].annotation_id;showDetails(hits[0]);return;}
  selectedId=null;
  var res=residueAt(i);
  showCard(res,[],res?(res.component_id+' '+res.author_residue_number):PLACEHOLDER);}
// DOM ONLY. Clearing the card because the renderer's selection has just gone
// empty must not write the selection back -- see the selection-change handler.
function clearCard(){
  var rows=document.querySelectorAll('[data-annotation]');
  for(var i=0;i<rows.length;i++)rows[i].removeAttribute('data-on');
  var d=document.getElementById('pinsp-residue-details');if(d)d.textContent=PLACEHOLDER;
  var card=document.getElementById('pinsp-card');if(card){card.hidden=true;card.innerHTML='';}
  selectedId=null;}
function clearDetails(){
  clearCard();
  var r=renderer();
  if(r){applying=true;try{r.setResidueSelection(new Set());
    r.render('Protein Inspector clear residue selection');}finally{applying=false;}}
  seqPaintSelection();syncSelTools();}
function selectInStructure(res){var r=renderer();if(!r)return;
  applying=true;try{r.setResidueSelection(new Set(indicesFor(res)));
    r.render('Protein Inspector manifest residue selection');}finally{applying=false;}
  seqPaintSelection();syncSelTools();}
// ===========================================================================
// THE SEQUENCE, AND SELECTING WITH IT
//
// The Residues panel lists what was ANNOTATED. A reader asking "what is next
// to this one", "what is residue 118", or "show me that side chain" has no way
// in from there -- and until now no way in from the picture either, because
// the notebook bundle ships with canvas picking switched OFF (parts/ui.js
// turns it on for Focus mode alone).
//
// WRITTEN HERE RATHER THAN BORROWED. py2Dmol has a sequence panel,
// src/panels/seq.js, and it is not in the bundle this file inlines: it is a
// canvas strip wired to the web app's own selection module (src/app/
// selection.js, also absent), it keys its colours off the viewer's colour
// mode, and it knows nothing about inspection layers. Pulling it in means
// adding two web-only modules and rebuilding the bundle. What follows reads
// the renderer's own arrays -- the same arrays the canvas draws from -- and
// paints each letter in its LAYER colour, which is the one thing that makes
// the strip an index to this page rather than a second structure viewer.
// ===========================================================================
var AA1={ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLU:'E',GLN:'Q',GLY:'G',HIS:'H',
  ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',
  VAL:'V',SEC:'U',PYL:'O',
  // A modified or D-residue reads as its parent's letter, as gemmi's
  // one_letter_code does; without these a phosphopeptide is a row of X.
  MSE:'M',HSD:'H',HSE:'H',HID:'H',HIE:'H',HIP:'H',PTR:'Y',SEP:'S',TPO:'T',FME:'M',
  HYP:'P',PCA:'E',ALY:'K',DAL:'A',DAR:'R',DSG:'N',DAS:'D',DCY:'C',DGN:'Q',DGL:'E',
  DHI:'H',DIL:'I',DLE:'L',DLY:'K',MED:'M',DPN:'F',DPR:'P',DSN:'S',DTH:'T',DTR:'W',
  DTY:'Y',DVA:'V',
  DA:'A',DT:'T',DC:'C',DG:'G',DU:'U',A:'A',U:'U',C:'C',G:'G',T:'T',PSU:'U'};
var SEQ_MAX=12000;
var seqBody=null,seqCells=[],seqRuns=[],seqBuilt=false,seqAnchor=-1,
    seqDragging=false,seqHoverAt=-1,seqAddr={},
    aroundBtn=null,aroundInput=null,scOn=null,scOff=null,selClearBtn=null,selOut=null,
    seqMode='add',seqBase=null;
function setList(set){var out=[];if(set&&set.forEach)set.forEach(function(v){out.push(v);});
  return out;}
function currentSelection(){var r=renderer();
  return (r&&r.residueSelection&&r.residueSelection.forEach)?r.residueSelection:new Set();}
function seqNames(r){
  if(r&&r.positionNames&&r.positionNames.length)return r.positionNames;
  var o=(s.viewer&&s.viewer.objects&&s.viewer.objects[0])||null;
  var f=o&&o.frames&&o.frames[0];
  return (f&&f.position_names)||[];}
function residueAt(i){var r=renderer();if(!r)return null;
  var nm=seqNames(r),n=(r.residueNumbers||[])[i];
  return {component_id:(nm[i]==null?'':String(nm[i])),
          chain_id:String((r.chains||[])[i]),
          author_residue_number:(n==null?i:n),
          canonical_position:i+1};}
function annotationsAt(i){return s.annotations.filter(function(a){
  return a.resolved&&indicesFor(a.residue).indexOf(i)>=0;});}
// --- building the strip ----------------------------------------------------
// ONE ROW PER CHAIN OF ONE OBJECT. chainKeyAt, not the bare letter: with a
// partner merged in, chain B of the partner and chain B of the target are two
// molecules and one row each.
function seqChainRuns(r){
  var ch=r.chains||[],out=[],cur=null;
  for(var i=0;i<ch.length;i++){
    var key=(typeof r.chainKeyAt==='function')?r.chainKeyAt(i):ch[i];
    if(!cur||cur.key!==key){cur={key:key,label:String(ch[i]),from:i,to:i};out.push(cur);}
    else cur.to=i;}
  return out;}
function seqBuild(){
  if(seqBuilt)return true;
  var r=renderer();if(!r||!seqBody)return false;
  var ch=r.chains||[];if(!ch.length)return false;
  if(ch.length>SEQ_MAX){
    // A capsid is not a sequence strip. The Residues panel and the structure
    // itself still work; drawing 300,000 letters does not.
    seqBody.innerHTML='<p class="pinsp-seqnote">'+ch.length+' positions \u2014 too many '
      +'to draw as a sequence. Select in the structure instead; the tools above '
      +'still apply.</p>';
    seqBuilt=true;syncSelTools();return true;}
  var num=r.residueNumbers||[],nm=seqNames(r),ty=r.positionTypes||[];
  seqRuns=seqChainRuns(r);seqAddr={};
  var html='';
  for(var k=0;k<seqRuns.length;k++){var run=seqRuns[k],cells='';
    for(var i=run.from;i<=run.to;i++){
      var name=String(nm[i]==null?'':nm[i]).trim().toUpperCase(),
          n=(num[i]==null)?i:num[i],
          lig=(ty[i]==='L'),
          letter=lig?(name.charAt(0)||'*'):(AA1[name]||'X'),
          tick=(!lig&&typeof n==='number'&&n%10===0)?(' data-num="'+n+'"'):'';
      seqAddr[run.label+'|'+n]=i;
      cells+='<span class="pinsp-cell" data-i="'+i+'"'+(lig?' data-lig="1"':'')+tick
        +' title="'+esc(run.label+' \u00b7 '+(name||'?')+' '+n)+'">'+esc(letter)+'</span>';}
    html+='<div class="pinsp-seqchain"><button type="button" class="pinsp-seqlab" '
      +'data-run="'+k+'" title="Select every residue of this chain">'+esc(run.label)
      +'</button><div class="pinsp-seqcells">'+cells+'</div></div>';}
  seqBody.innerHTML=html;
  seqCells=[];
  var found=seqBody.querySelectorAll('.pinsp-cell');
  for(var c=0;c<found.length;c++)seqCells[Number(found[c].getAttribute('data-i'))]=found[c];
  seqBuilt=true;seqPaintLayers();seqPaintSelection();syncSelTools();
  return true;}
// --- what the letters look like -------------------------------------------
// THE SAME COLOUR THE RESIDUE HAS IN 3D, from the same visible-layer set, so
// the strip cannot disagree with the picture. Unannotated letters are left to
// the stylesheet rather than painted, for the same reason the structure's base
// is: a colour on every position is a colour that says nothing.
function seqDark(hex){var n=parseInt(String(hex).slice(1),16);
  if(!isFinite(n))return false;
  return (0.299*((n>>16)&255)+0.587*((n>>8)&255)+0.114*(n&255))<150;}
function seqPaintLayers(){
  if(!seqBuilt||!seqCells.length)return;
  var colors={};
  for(var li=0;li<s.layers.length;li++){var l=s.layers[li];
    if(!layerVisible(l.layer_id))continue;
    for(var ai=0;ai<l.annotation_ids.length;ai++){var a=byId[l.annotation_ids[ai]];
      if(!a||!a.resolved||!a.color)continue;
      var ix=indicesFor(a.residue);
      for(var m=0;m<ix.length;m++)colors[ix[m]]=a.color;}}
  for(var i=0;i<seqCells.length;i++){var cell=seqCells[i];if(!cell)continue;
    var hex=colors[i];
    cell.style.background=hex||'';
    cell.style.color=hex?(seqDark(hex)?'#ffffff':'#0f172a'):'';}}
function seqPaintSelection(){
  if(!seqBuilt||!seqCells.length)return;
  var sel=currentSelection();
  for(var i=0;i<seqCells.length;i++){var cell=seqCells[i];if(!cell)continue;
    if(sel.has(i))cell.setAttribute('data-sel','1');else cell.removeAttribute('data-sel');}}
function seqHover(i){
  if(i===seqHoverAt)return;
  var prev=seqCells[seqHoverAt];if(prev)prev.removeAttribute('data-hover');
  seqHoverAt=i;
  var cell=seqCells[i];if(cell)cell.setAttribute('data-hover','1');}
// THE VIEWER ALREADY LOOKS FOR A STRIP. core/mol.js calls
// window.SEQ.setHoveredResidue on mousemove when one is present and skips it
// when it is not, so hovering the structure lights the letter with no change
// on the viewer's side -- and buildView is its "the structure changed" signal.
// Installed only once the strip exists, so the pick-per-move cost is never
// paid by a page that has nothing to show for it.
function installSeqBridge(){
  if(window.SEQ&&window.SEQ.pinsp)return;
  window.SEQ={pinsp:true,
    setHoveredResidue:function(info){
      if(!info){seqHover(-1);return;}
      var i=seqAddr[String(info.chain)+'|'+info.resSeq];
      seqHover(i===undefined?-1:i);},
    setHoverAtoms:function(){},
    clear:function(){seqHover(-1);},
    buildView:function(){seqBuilt=false;seqBuild();},
    buildViewDeferred:function(){requestAnimationFrame(function(){
      seqBuilt=false;seqBuild();});}};}
// --- selecting -------------------------------------------------------------
function seqRunOf(i){for(var k=0;k<seqRuns.length;k++)
  if(i>=seqRuns[k].from&&i<=seqRuns[k].to)return seqRuns[k];
  return null;}
// A DRAG STAYS IN ONE CHAIN. Dragging from the end of chain A into chain B
// would otherwise select the gap between two molecules as though it were a
// stretch of sequence.
function seqRange(a,b){var run=seqRunOf(a),lo=Math.min(a,b),hi=Math.max(a,b),out=new Set();
  if(run){lo=Math.max(lo,run.from);hi=Math.min(hi,run.to);}
  for(var i=lo;i<=hi;i++)out.add(i);
  return out;}
function applySelection(set){
  var r=renderer();if(!r)return;
  applying=true;try{r.setResidueSelection(set);
    r.render('Protein Inspector residue selection');}finally{applying=false;}
  seqPaintSelection();syncSelTools();
  var list=setList(set);
  if(list.length===1)showIndex(list[0]);else clearCard();}
function describeSelection(list){
  if(!list.length){
    // THE PER-RESIDUE STATE OUTLIVES THE SELECTION, so say so. Which residues
    // draw their side chains is a fact about the residues, not about what is
    // picked -- dropping the selection does not undraw them, and a bar that
    // went blank made it look as though it had.
    var shown=scShownSet(),n=shown?shown.size:0;
    return n?('Nothing selected \u2014 '+n+' residue'+(n===1?'':'s')
              +' still drawing side chains.')
            :'Click the structure or a letter to select a residue.';}
  var r=renderer();if(!r)return list.length+' residues selected';
  var ch=r.chains||[],num=r.residueNumbers||[],
      sorted=list.slice().sort(function(a,b){return a-b;}),parts=[],run=null;
  function flush(){if(run)parts.push(run.c+' '+(run.first===run.last?run.first
    :(run.first+'\u2013'+run.last)));}
  for(var k=0;k<sorted.length;k++){var i=sorted[k],c=String(ch[i]),
      n=(num[i]==null)?i:Number(num[i]);
    if(run&&run.c===c&&i===run.lasti+1&&n===run.last+1){run.last=n;run.lasti=i;continue;}
    flush();run={c:c,first:n,last:n,lasti:i};}
  flush();
  // Joined with the count because the ranges alone do not say how many, and
  // the count alone does not say which.
  return parts.join(', ')+'  \u00b7  '+sorted.length+' residue'+(sorted.length===1?'':'s');}
// --- side chains of the selection -----------------------------------------
// WHICH RESIDUES SHOW THEIRS is the renderer's own set, in the OBJECT's
// numbering, so the button reads it back rather than keeping a second copy
// that the Style panel's own Focus mode could silently disagree with.
// WHICH RESIDUES DRAW THEIRS, in MERGED indices. `shownSidechainSet` is the
// renderer's own answer and already merged, which is what the selection is in;
// the fallback reads the object's own set, which is in the OBJECT's numbering
// and so needs the window offset.
function scShownSet(){
  var r=renderer();if(!r)return null;
  if(typeof r.shownSidechainSet==='function'){
    var merged=r.shownSidechainSet();
    if(merged&&merged.has)return merged;}
  var obj=r.objectsData&&r.objectsData[r.currentObjectName];
  if(!obj||!obj.sidechains||!obj.sidechains.has)return null;
  var off=0;
  if(typeof r.localRangeOf==='function'){
    var w=r.localRangeOf(r.currentObjectName);
    if(w&&typeof w.off==='number')off=w.off;}
  if(!off)return obj.sidechains;
  var out=new Set();
  obj.sidechains.forEach(function(i){out.add(i+off);});
  return out;}
// true: every one of them draws. false: none does. null: they disagree -- the
// state a single latch could not say, which is why this is a pair.
function scStateOf(list){
  if(!list.length)return false;
  var set=scShownSet(),hit=0;
  if(!set)return false;
  for(var i=0;i<list.length;i++)if(set.has(list[i]))hit++;
  if(!hit)return false;
  return hit===list.length?true:null;}
// ABSOLUTE, NOT A TOGGLE. The button says which state to put the selection in,
// so a mixed selection resolves to that state in one press instead of flipping
// to whichever one the majority happened not to be in.
function setSidechains(on){
  var r=renderer(),list=setList(currentSelection());
  if(!r||!list.length||typeof r.showSidechains!=='function')return;
  try{if(on)r.showSidechains({positions:list});
    else r.hideSidechains({positions:list});}
  catch(err){if(selOut)selOut.textContent=String((err&&err.message)||err);}
  syncSelTools();}
function selectAround(){
  var r=renderer();if(!r)return;
  var list=setList(currentSelection());
  if(!list.length||typeof r.residuesWithin!=='function')return;
  var a=parseFloat(aroundInput&&aroundInput.value);
  if(!isFinite(a)||a<=0)a=5;
  // ATOM TO ATOM, and the seed comes back with the answer -- PyMOL's
  // `byres (all within N of sele)`. A shell that drops what it was measured
  // from is a different question, and one more press of the same button
  // cannot un-ask this one.
  applySelection(r.residuesWithin(new Set(list),a));}
function syncSelTools(){
  var r=renderer(),list=setList(currentSelection());
  if(selOut)selOut.textContent=describeSelection(list);
  if(aroundBtn)aroundBtn.disabled=!(list.length&&r&&typeof r.residuesWithin==='function');
  if(selClearBtn)selClearBtn.disabled=!list.length;
  if(scOn&&scOff){
    var able=!!s.sidechains&&!!r&&typeof r.showSidechains==='function',
        // Nothing selected fills NEITHER, like a disagreeing selection: there
        // is no answer to "are these drawn" when there are no these.
        state=(able&&list.length)?scStateOf(list):null;
    scOn.disabled=scOff.disabled=!(able&&list.length);
    scOn.setAttribute('aria-pressed',String(state===true));
    scOff.setAttribute('aria-pressed',String(state===false));
    var why=able?(list.length
        ?(state===null?'These residues disagree \u2014 Show draws all of them'
                      :'Draw the selected residues\u2019 side chains')
        :'Select residues first; what is already drawn stays drawn')
      :'This bundle carries no side-chain atoms';
    scOn.title=why;scOff.title=why;}}
// DETAIL SURVIVES A STYLE SWITCH.
//
// `display={"detail": N}` seeds config.rendering.detail, and the constructor
// honours it -- but only for as long as the viewer stays in the style it
// opened in. core/mol.js's _applyLookDefaults assigns `cartoonDetail =
// d.detail` UNCONDITIONALLY, with no "did a person choose this" latch of the
// kind width and thickness have, and all three entries in LOOK_DEFAULTS carry
// detail: DETAIL_DEFAULT_N (4). A bundle opens as `tube`, so the first time a
// reader picks Richardson from the Style dropdown the sampling drops from 8 to
// 4 -- and Richardson is the style they pick it for. The bundle shipped with
// detail 8 that nobody could ever see, which is worse than shipping 4.
//
// That is not upstream's bug to fix: re-asserting every style-owned control is
// deliberate there, and the comment at the top of LOOK_DEFAULTS explains what
// it is for ("a high Outline set in cartoon silently survived into richardson,
// whose panel hides the slider, leaving no way to see or fix it"). Detail is
// simply not style-owned as far as a bundle is concerned -- it is a property
// of the document, chosen once by whoever rendered it.
//
// So the wrapper re-asserts it after the switch, and the reader still outranks
// us: dragging the Detail slider makes their value the one that is kept.
var wantDetail=(function(){
  var d=s.viewer&&s.viewer.config&&s.viewer.config.rendering
    &&s.viewer.config.rendering.detail;
  var n=parseInt(d,10);
  return (isFinite(n)&&n>=2&&n<=8)?n:null;})();
function keepDetail(){
  var r=renderer();if(!r||wantDetail===null)return;
  if(r.cartoonDetail===wantDetail)return;
  r.cartoonDetail=wantDetail;
  var sl=document.getElementById('detailSlider');
  if(sl)sl.value=String(wantDetail);
  r.render('Protein Inspector detail');}
function installDetailKeeper(){
  var r=renderer();if(!r||wantDetail===null)return false;
  if(r._pinspDetailKept)return true;
  r._pinspDetailKept=true;
  // WRAPPED RATHER THAN LISTENED FOR. A change listener on #styleSelect would
  // have to run after py2Dmol's own, which is an ordering assumption; and it
  // would miss setStyle called from anywhere else -- the loader picks tube
  // over cartoon by itself on a large structure.
  ['setStyle','setPreset'].forEach(function(name){
    var fn=r[name];if(typeof fn!=='function')return;
    r[name]=function(){var out=fn.apply(r,arguments);keepDetail();return out;};});
  var sl=document.getElementById('detailSlider');
  if(sl)['input','change'].forEach(function(ev){
    sl.addEventListener(ev,function(){
      var n=parseInt(sl.value,10);
      if(isFinite(n))wantDetail=n;});});
  keepDetail();
  return true;}

// THE STYLE PANEL, FOLDED -- NOT REWRITTEN.
//
// The panel comes from src/parts/panel.js, whose own header calls it "ONE
// PANEL, TWO SKINS": the website and this bundle build their controls from the
// SAME table of rows, and only the stylesheet differs. So "the site has fewer
// controls" is not true and there is nothing to copy; the site is a wide left
// sidebar with room for seventeen controls, and ours is a small panel floating
// over the canvas, where the same seventeen read as clutter.
//
// Four of them answer questions a reader of a bundle actually has -- which
// representation (Style), what are the colours (Color), how is my selection
// marked (Sele), and how finely is it sampled (Detail, the one that decides
// whether a magnified residue facets). Those stay out front. The rest are
// authoring controls; they go into a closed <details> and are one click away.
//
// Whole ROWS move, not individual controls: parts/ui.js hides per-style
// controls by setting `hidden` on the row that carries data-style, so a control
// lifted out of its row would be a slider that survives into a style that
// cannot honour it. Detail is the exception and is safe -- it is a `.half`
// cell, which carries its own data-style, and it is moving INTO the always-
// visible Style row, where it is visible under tube too (tube ignores it).
//
// Idempotent, and it returns false until the viewer has built the panel, so
// the init loop below can keep asking.
function foldStylePanel(){
  var panel=document.getElementById('stylePanel');
  if(!panel)return false;
  if(panel.querySelector('.pinsp-fine'))return true;
  var styleRow=null,colorRow=null,rows=[];
  var kids=panel.children;
  for(var i=0;i<kids.length;i++){
    var row=kids[i];
    if(row.className&&String(row.className).indexOf('toggle-item')<0)continue;
    if(row.querySelector('#styleSelect'))styleRow=row;
    else if(row.querySelector('#colorSelect'))colorRow=row;
    else rows.push(row);}
  if(!styleRow)return false;
  // Detail rides up beside Style, out of the row it shares with Pencil, Ink,
  // Shadow and Ortho -- all four of which are about to be folded away.
  var det=panel.querySelector('#detailSlider');
  var cell=det&&det.parentNode;
  if(cell&&cell.className&&String(cell.className).indexOf('half')>=0)styleRow.appendChild(cell);
  var fold=document.createElement('details');
  fold.className='pinsp-fine';
  var sum=document.createElement('summary');
  sum.textContent='Fine tuning';
  sum.title='Width, outline, thickness, flattening, lighting and grain.'
    +' Nothing here is hidden from you \\u2014 it is folded because a reader'
    +' rarely needs it.';
  fold.appendChild(sum);
  // THE PAINTER, AS A CONTROL -- because the two do not draw the same picture.
  //
  // parts/ui.js deliberately leaves this out of the notebook panel ("useGPU is
  // a rendering BACKEND"), and that is right for a notebook, where the build
  // carries one painter and the flag would be a request for a file that is not
  // there. A bundle is the build that carries BOTH, and the choice is real:
  // reported on a side chain drawn close up, the 2D painter's ink pass removes
  // the edges of a stick box that lies inside another, leaving one silhouette,
  // while the WebGL2 path showed the interior edges through the faces. Whoever
  // is looking at the picture is the only one who can say which is right for
  // it, and they cannot say it if the only way to switch is to re-render the
  // file.
  //
  // Only built when both painters are present and WebGL2 is actually available
  // -- py2Dmol's own rule from src/app/main.js, where a control that cannot do
  // anything is worse than none.
  var G=window.py2dmolCartoonGPU;
  if(G&&G.available&&G.available()&&typeof window.py2dmolCartoonPaint==='function'){
    var prow=document.createElement('div');prow.className='toggle-item';
    var lab=document.createElement('label');lab.className='btn-toggle';
    lab.title='WebGL2 painter. Off re-draws on the CPU every frame, which is'
      +' slower to turn but carves one silhouette per side chain instead of'
      +' showing the edges of each stick box.';
    var box=document.createElement('input');box.type='checkbox';
    box.id='pinsp-gpu';
    var rr=renderer();box.checked=!!(rr&&rr.useGPU);
    var cap=document.createElement('span');cap.textContent='GPU painter';
    lab.appendChild(box);lab.appendChild(cap);prow.appendChild(lab);
    box.addEventListener('change',function(){
      var r=renderer();if(!r)return;
      r.useGPU=box.checked;
      // The mesh was built from state that may have moved on while this was
      // off, so ask for a fresh one rather than trusting it -- main.js does
      // the same on its own checkbox.
      if(G.invalidate)G.invalidate();
      r.render('Protein Inspector painter switch');});
    fold.appendChild(prow);}
  for(var k=0;k<rows.length;k++)fold.appendChild(rows[k]);
  panel.appendChild(fold);
  if(colorRow)panel.insertBefore(colorRow,fold);
  return true;}
// SELECTING NEVER DRAWS ANYTHING. Focus mode is the exception, so it is made
// to stop being one.
//
// Nothing on this page draws a side chain except the Show/Hide pair -- but
// py2Dmol's Focus button, which is in the bundle's own control bar, is a MODE
// that redefines what a selection means. `enterFocusMode` clears every
// object's `sidechains` on the way in ("side chains they turned on by hand
// look like ones focus drew"), and `focusOn` then ASSIGNS the set from the
// neighbourhood of whatever was picked, writing every object on screen. It is
// wrapped around setResidueSelection rather than subscribed to it, so with the
// mode latched EVERY selection draws atoms: a click in the picture, a letter in
// the strip, `sel` on a layer, `Select around`. Which is exactly the
// complaint -- selecting a residue showed its side chains before anyone asked.
//
// Upstream already has the opt-out: focusOn's options carry `sidechains`, and
// `o.sidechains !== false` is the only gate on that whole block. So the camera
// half of the mode is kept -- it still flies to what you picked and clips to
// it, which is what the button is for -- and the drawn set is left to the one
// control that is about the drawn set. The snapshot/restore on the way in and
// out is neutered the same way: with nothing of its own to put back, there is
// nothing for it to overwrite.
//
// Wrapped, not re-implemented: the mode is upstream's and this says one thing
// about it.
function tameFocus(){
  var r=renderer();if(!r||r._pinspFocusTamed)return false;
  if(typeof r.focusOn!=='function')return true;   // a build without the part
  r._pinspFocusTamed=true;
  var focusOn=r.focusOn;
  r.focusOn=function(sel,opts){
    var o={};if(opts)for(var k in opts)if(Object.prototype.hasOwnProperty.call(opts,k))o[k]=opts[k];
    o.sidechains=false;
    return focusOn.call(r,sel,o);};
  // ...AND THE WIPE ON THE WAY IN. enterFocusMode drops every object's set
  // before the first click, so pressing the button undrew side chains the
  // reader had asked for. Put them back: the mode has no side-chain business
  // of its own any more, so there is nothing legitimate in that clearing.
  if(typeof r.enterFocusMode==='function'){
    var enter=r.enterFocusMode;
    r.enterFocusMode=function(){
      var keep={},data=r.objectsData||{},name;
      for(name in data)if(data[name]&&data[name].sidechains)
        keep[name]=new Set(data[name].sidechains);
      var out=enter.call(r),restored=false;
      for(name in keep)if(data[name]){data[name].sidechains=keep[name];restored=true;}
      if(restored){
        if(r._invalidateSegmentCache)r._invalidateSegmentCache();
        if(r.reloadDrawn)r.reloadDrawn(true);
        r.render('Protein Inspector focus keeps the side chains');}
      return out;};}
  return true;}
function enableClicks(){
  var r=renderer();if(!r)return false;
  // The notebook bundle ships with picking off, so an exported page's canvas
  // was inert: a reader could see a residue and had no way to ask what it is.
  // Single click selects one, shift-click extends, double-click takes the
  // chain, a click on the background clears -- all of it already written in
  // core/mol.js behind this one flag.
  r.selectionEnabled=true;return true;}
(function(){
  seqBody=document.getElementById('pinsp-seqbody');
  aroundBtn=document.getElementById('pinsp-around');
  aroundInput=document.getElementById('pinsp-around-a');
  scOn=document.getElementById('pinsp-sc-on');
  scOff=document.getElementById('pinsp-sc-off');
  selClearBtn=document.getElementById('pinsp-sel-clear');
  selOut=document.getElementById('pinsp-sel');
  if(aroundBtn)aroundBtn.addEventListener('click',selectAround);
  if(scOn)scOn.addEventListener('click',function(){setSidechains(true);});
  if(scOff)scOff.addEventListener('click',function(){setSidechains(false);});
  if(selClearBtn)selClearBtn.addEventListener('click',clearDetails);
  if(!seqBody)return;
  // ONE LISTENER FOR EVERY LETTER. A page with 20,000 cells cannot afford a
  // handler each, and the cells are leaves, so the event target IS the cell.
  function cellAt(t){var v=(t&&t.getAttribute)?t.getAttribute('data-i'):null;
    return v==null?-1:Number(v);}
  function runAt(t){var v=(t&&t.getAttribute)?t.getAttribute('data-run'):null;
    return v==null?-1:Number(v);}
  // THE STRIP ADDS AND SUBTRACTS; IT NEVER REPLACES.
  //
  // A click used to throw the selection away and start again, which makes the
  // strip useless for the thing it is for: building a set a few residues at a
  // time, out of two chains, around what the layers already gave you. What the
  // residue under the pointer ALREADY IS decides the gesture -- start on an
  // unselected letter and the drag adds, start on a selected one and it takes
  // away, which is one rule for both directions and needs no modifier.
  //
  // Each move recomputes from the selection as it was at mousedown, not from
  // the last frame's result, so dragging back over your own path shrinks the
  // range instead of ratcheting it.
  function beginSeqDrag(i,remove){
    seqBase=new Set(setList(currentSelection()));
    seqMode=remove?'remove':'add';
    seqAnchor=i;seqDragging=true;}
  function applySeqMode(range){
    var next=new Set(seqBase?setList(seqBase):[]);
    range.forEach(function(k){
      if(seqMode==='remove')next.delete(k);else next.add(k);});
    applySelection(next);}
  seqBody.addEventListener('mousedown',function(e){
    var i=cellAt(e.target);
    if(i>=0){
      if(e.preventDefault)e.preventDefault();
      // Shift keeps its usual meaning: extend from the last letter pressed,
      // in whichever direction that gesture was already going.
      if(e.shiftKey&&seqAnchor>=0){
        seqBase=new Set(setList(currentSelection()));
        applySeqMode(seqRange(seqAnchor,i));return;}
      beginSeqDrag(i,currentSelection().has(i));
      applySeqMode(seqRange(i,i));return;}
    var k=runAt(e.target);
    if(k>=0&&seqRuns[k]){
      // A chain button reads the same way: if the whole chain is in, it comes
      // out; otherwise it goes in. Nothing else in the selection is touched.
      var run=seqRuns[k],all=new Set(),cur=currentSelection(),every=true;
      for(var j=run.from;j<=run.to;j++){all.add(j);if(!cur.has(j))every=false;}
      seqBase=new Set(setList(cur));seqMode=every?'remove':'add';
      seqAnchor=run.from;applySeqMode(all);}});
  seqBody.addEventListener('mousemove',function(e){
    var i=cellAt(e.target);
    if(seqDragging){if(i>=0)applySeqMode(seqRange(seqAnchor,i));return;}
    seqHover(i);});
  seqBody.addEventListener('mouseleave',function(){seqHover(-1);});
  window.addEventListener('mouseup',function(){seqDragging=false;});
})();
// ONE SEAM, TWO ROWS, CONSTANT SUM -- which is what the vertical seam has
// always been, on the other axis.
//
// THE FIRST VERSION MEASURED AGAINST THE GRID'S OWN BOTTOM EDGE:
//
//     setSeq(grid.getBoundingClientRect().bottom - e.clientY)
//
// and that edge MOVES when the strip resizes. Every pointermove therefore read
// back its own previous result, so a 2px drag upward grew the strip, which
// pushed `bottom` down, which made the next event ask for more again -- the
// seam ran away from the pointer to its ceiling in a few events. The vertical
// seam is immune to the same expression because it reads `right`, and the grid
// is full width: changing the panel's width does not move the page's edge.
//
// So: ANCHOR THE GESTURE. The pointer's travel since pointerdown is the only
// quantity that matters, and nothing it writes can feed back into it.
//
// AND THE PICTURE GIVES UP WHAT THE STRIP TAKES. The first version only grew
// the strip, so the whole page got taller and the only thing a reader could
// see move was the credit line under it -- which is exactly what was reported.
// The two rows share a fixed sum instead: drag up and the canvas is shorter
// and the strip is taller, the page's height does not change, and nothing
// below it moves. The canvas follows because py2Dmol's own ResizeObserver
// watches #canvasContainer for height as well as width (parts/viewport.js),
// so one custom property is the whole mechanism.
// Set by the seam below, called by anything that changes how much room the
// grid has: a window resize, and coming back from an About tab, whose pane
// was display:none while the grid could not measure itself.
var pinspFit=null;
(function(){
  var grid=document.querySelector('.protein-inspector'),
      split=document.getElementById('pinsp-hsplit');
  if(!grid||!split)return;
  var MIN_SEQ=64,MIN_VIEW=200,dragging=false,
      anchorY=0,anchorSeq=0,anchorSum=0,homeSeq=0,homeSum=0;
  function cssPx(name,fallback){
    var v=parseFloat(getComputedStyle(grid).getPropertyValue(name));
    return isFinite(v)&&v>0?v:fallback;}
  function seqPx(){return cssPx('--pinsp-seq',152);}
  function viewPx(){return cssPx('--pinsp-view',720);}
  function apply(seq,sum){
    var most=Math.max(MIN_SEQ,sum-MIN_VIEW),
        value=Math.round(Math.min(Math.max(seq,MIN_SEQ),most));
    grid.style.setProperty('--pinsp-seq',value+'px');
    grid.style.setProperty('--pinsp-view',Math.round(sum-value)+'px');
    split.setAttribute('aria-valuenow',String(value));}
  homeSeq=seqPx();homeSum=homeSeq+viewPx();
  // AND THE STACK ENDS WHERE THE WINDOW DOES.
  //
  // Both rows were fixed pixel heights -- the picture's from the caller's
  // `display` height, the strip's from its own default -- so the page was
  // exactly as tall as the sum of them plus the body's padding, and on any
  // window taller than that the sequence section stopped short of the bottom
  // with a band of blank page under it. Moving the credit inside the section
  // did not help, because the section itself was the thing ending early.
  //
  // So the SUM is the viewport's, measured from where the grid actually
  // starts (the tab strip sits above it) less the padding the page writes
  // itself. The strip keeps the height it was given and the PICTURE takes the
  // slack, which is the right way round: the extra space goes to the
  // structure, and the seam still decides the division.
  //
  // The seam's own contract is untouched -- it reads the sum fresh at
  // pointerdown, so "two rows, constant sum" is now constant at the window's
  // height rather than at the caller's. The one thing `display={"height"}`
  // loses is the final say on the canvas, which it did not have anyway: the
  // display area has always been fluid in width, and this makes it fluid in
  // the other axis for the same reason.
  var PAGE_PAD=14;        // body{padding:14px}, written by this page
  function availPx(){
    if(typeof grid.getBoundingClientRect!=='function')return 0;
    var box=grid.getBoundingClientRect(),h=window.innerHeight;
    if(!box||!isFinite(box.top)||!isFinite(h))return 0;
    return Math.round(h-box.top-PAGE_PAD);}
  function fitStack(){
    var avail=availPx();if(avail<=0)return;
    var sum=avail-(split.offsetHeight||11),floor=MIN_SEQ+MIN_VIEW;
    if(sum<floor)sum=floor;
    homeSum=sum;
    apply(seqPx(),sum);}
  pinspFit=fitStack;
  fitStack();
  window.addEventListener('resize',fitStack);
  function fromPointer(e){if(dragging)apply(anchorSeq+(anchorY-e.clientY),anchorSum);}
  split.addEventListener('pointerdown',function(e){
    dragging=true;anchorY=e.clientY;anchorSeq=seqPx();anchorSum=anchorSeq+viewPx();
    split.setAttribute('data-drag','1');
    if(split.setPointerCapture)try{split.setPointerCapture(e.pointerId);}catch(err){}
    e.preventDefault();});
  split.addEventListener('pointermove',fromPointer);
  window.addEventListener('pointermove',fromPointer);
  function stop(){if(!dragging)return;dragging=false;split.removeAttribute('data-drag');}
  split.addEventListener('pointerup',stop);
  window.addEventListener('pointerup',stop);
  window.addEventListener('pointercancel',stop);
  split.addEventListener('keydown',function(e){
    var step=e.shiftKey?48:16,sum=seqPx()+viewPx();
    if(e.key==='ArrowUp'){apply(seqPx()+step,sum);e.preventDefault();}
    else if(e.key==='ArrowDown'){apply(seqPx()-step,sum);e.preventDefault();}
    else if(e.key==='Home'){apply(homeSeq,homeSum);e.preventDefault();}});
})();
function setAllLayers(on){for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
    if(b)b.checked=on;}syncVisibleLayers();}
// A LAYER IS A SET OF RESIDUES, SO IT CAN BE A SELECTION.
//
// Everything the selection tools do -- grow it by distance, draw its side
// chains -- was reachable only for residues picked one at a time, which is the
// wrong unit: the interesting question is almost always "this whole epitope"
// or "every hotspot", and those are already named, counted and coloured in the
// Layers panel. Selecting a layer is therefore not a new concept, it is the
// one that was missing between the two halves of the page.
function layerPositions(id){
  var layer=layerOf[id],out=[];
  if(!layer)return out;
  for(var i=0;i<layer.annotation_ids.length;i++){
    var a=byId[layer.annotation_ids[i]];
    if(!a||!a.resolved)continue;
    var ix=indicesFor(a.residue);
    for(var k=0;k<ix.length;k++)out.push(ix[k]);}
  return out;}
function selectLayer(id,add){
  var next=add?new Set(setList(currentSelection())):new Set(),
      list=layerPositions(id);
  for(var i=0;i<list.length;i++)next.add(list[i]);
  if(!next.size)return;
  applySelection(next);}
// ...AND EVERY VISIBLE LAYER AT ONCE, which is what the reader means by "show
// the side chains for these": the layers they have left ticked. Reads the
// same visibility the structure is painted from, so what gets selected is
// what is on screen and nothing else.
function selectVisibleLayers(){
  var next=new Set();
  for(var i=0;i<s.layers.length;i++){
    var id=s.layers[i].layer_id;
    if(!layerVisible(id))continue;
    var list=layerPositions(id);
    for(var k=0;k<list.length;k++)next.add(list[k]);}
  if(!next.size)return;
  applySelection(next);}
function onlyLayer(id){for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
    if(b)b.checked=(s.layers[i].layer_id===id);}syncVisibleLayers();}
for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
  if(b)b.addEventListener('change',syncVisibleLayers);}
var onlys=document.querySelectorAll('[data-only]');
for(var o=0;o<onlys.length;o++)(function(btn){btn.addEventListener('click',function(ev){
  ev.preventDefault();ev.stopPropagation();onlyLayer(btn.getAttribute('data-only'));});})(onlys[o]);
var picks=document.querySelectorAll('[data-select]');
for(var q=0;q<picks.length;q++)(function(btn){btn.addEventListener('click',function(ev){
  ev.preventDefault();ev.stopPropagation();
  // Shift adds, so several layers can be selected together -- the same
  // modifier the strip and the canvas use for the same meaning.
  selectLayer(btn.getAttribute('data-select'),!!(ev&&ev.shiftKey));});})(picks[q]);
var rws=document.querySelectorAll('[data-annotation]');
for(var w=0;w<rws.length;w++)(function(row){row.addEventListener('click',function(){
  var id=row.dataset?row.dataset.annotation:row.getAttribute('data-annotation');
  if(selectedId===id){clearDetails();return;}
  selectedId=id;var a=byId[id];if(!a)return;showDetails(a);selectInStructure(a.residue);});})(rws[w]);
var bAll=document.getElementById('pinsp-layers-all');
if(bAll)bAll.addEventListener('click',function(){setAllLayers(true);});
var bNone=document.getElementById('pinsp-layers-none');
if(bNone)bNone.addEventListener('click',function(){setAllLayers(false);});
var bPick=document.getElementById('pinsp-layers-select');
if(bPick)bPick.addEventListener('click',selectVisibleLayers);
var bClear=document.getElementById('pinsp-clear-residue');
if(bClear)bClear.addEventListener('click',clearDetails);
var filt=document.getElementById('pinsp-filter');
if(filt)filt.addEventListener('input',function(){query=(filt.value||'').toLowerCase().trim();
  syncVisibleLayers();});
document.addEventListener('keydown',function(e){if(e.key==='Escape')clearDetails();});
// EVERY ROUTE IN LANDS HERE. The canvas (click, shift-click, double-click,
// click-on-background), the strip, the tools and the manifest rows all go
// through setResidueSelection, which dispatches this -- so the strip, the
// readout and the three tool buttons are refreshed in one place rather than
// at four call sites that can fall out of step.
document.addEventListener('py2dmol-residue-selection-change',function(){
  seqPaintSelection();syncSelTools();
  // `applying` means this page asked for the selection and has already said
  // what the card should show. Only a gesture in the CANVAS gets past here.
  if(applying)return;
  var r=renderer();if(!r)return;
  // AND THE MARK HAS TO BE DRAWN. core/mol.js's mouseup handler sets the
  // selection and then renders only `isLargeMolecule`, on the reasoning that
  // a small structure's next frame is along shortly -- which is true on the
  // website, where a hover readout and a sequence canvas are repainting
  // anyway, and false in a still bundle. So the yellow outline appeared on
  // the first frame AFTER the click: nudge the structure and there it was.
  // One render here covers every canvas route in, including the
  // clearResidueSelection a click on the background makes.
  r.render('Protein Inspector selection mark');
  var list=setList(r.residueSelection);
  if(!list.length){clearCard();return;}
  // The card is about ONE residue. A chain or a neighbourhood is reported by
  // the readout instead, which can say what a card cannot: which ranges.
  if(list.length===1){showIndex(list[0]);return;}
  clearCard();});
window.proteinInspector={syncVisibleLayers:syncVisibleLayers,clearDetails:clearDetails,
  setAllLayers:setAllLayers,onlyLayer:onlyLayer,
  selectAnnotation:function(id){var a=byId[id];if(a){selectedId=id;showDetails(a);
    selectInStructure(a.residue);}},
  // The three new verbs, on the same object, so a reader who wants a figure
  // scripted rather than clicked has the same selection the panel has.
  selectResidues:function(positions){applySelection(new Set(positions||[]));},
  selectAround:function(angstroms){
    if(angstroms!==undefined&&aroundInput)aroundInput.value=String(angstroms);
    selectAround();},
  setSidechains:setSidechains,
  sidechainState:function(){return scStateOf(setList(currentSelection()));},
  selectLayer:selectLayer,
  selectVisibleLayers:selectVisibleLayers,
  selection:function(){return setList(currentSelection());}};
// Colour of the un-annotated structure is the VIEWER's business, set from its
// own Style panel; there is no duplicate control here. The only thing this page
// does is honour a manifest that asked for a flat custom base, and stand down
// the moment the Style panel says otherwise -- otherwise the base would
// silently override whatever the user just picked.
document.addEventListener('py2dmol-color-change',function(){
  var r=renderer();if(!r||!r.colorMode||baseMode===r.colorMode)return;
  baseMode=r.colorMode;syncVisibleLayers();});
// THE PICTURE FILLS ITS PANE, BY LAYOUT. The grid column owns the width;
// #canvasContainer is a flow child of the stage and takes all of it. The
// only thing JS still has to do is the one element CSS cannot reach: the
// viewer's root div carries an id and no class in an export, and it is an
// inline-block, so it would shrink-wrap its content and latch narrow. One
// write, once, at startup -- nothing here measures anything, so there is no
// path by which a width can feed back into itself.
(function(){
  var main=document.querySelector('#pinsp-stage #mainContainer');
  if(!main)return;
  var instance=main.closest('#pinsp-stage')===main.parentElement?null:main.parentElement;
  while(instance&&instance.id!=='pinsp-stage'){
    instance.style.setProperty('display','contents','important');
    instance=instance.parentElement;}
})();
// ONE BOUNDARY, DRAGGED. The panel width is a custom property on the grid,
// so moving the seam re-lays out both panes at once and they cannot come
// apart or overlap. The canvas takes its width from the column it is in and
// py2Dmol's own ResizeObserver redraws it.
(function(){
  var grid=document.querySelector('.protein-inspector'),
      split=document.getElementById('pinsp-split');
  if(!grid||!split)return;
  var MIN_PANEL=240,MIN_STAGE=320,dragging=false;
  function panelWidth(){
    var v=parseFloat(getComputedStyle(grid).getPropertyValue('--pinsp-panel'));
    return isFinite(v)&&v>0?v:340;}
  function setPanel(px){
    var box=grid.getBoundingClientRect(),
        most=Math.max(MIN_PANEL,box.width-MIN_STAGE);
    grid.style.setProperty('--pinsp-panel',
      Math.round(Math.min(Math.max(px,MIN_PANEL),most))+'px');
    split.setAttribute('aria-valuenow',String(Math.round(panelWidth())));}
  function fromPointer(e){setPanel(grid.getBoundingClientRect().right-e.clientX);}
  split.addEventListener('pointerdown',function(e){
    dragging=true;split.setAttribute('data-drag','1');
    if(split.setPointerCapture)try{split.setPointerCapture(e.pointerId);}catch(err){}
    e.preventDefault();});
  split.addEventListener('pointermove',function(e){if(dragging)fromPointer(e);});
  window.addEventListener('pointermove',function(e){if(dragging)fromPointer(e);});
  function stop(){if(!dragging)return;dragging=false;split.removeAttribute('data-drag');}
  split.addEventListener('pointerup',stop);
  window.addEventListener('pointerup',stop);
  window.addEventListener('pointercancel',stop);
  // A seam that only a mouse can move is a seam some readers cannot move.
  split.addEventListener('keydown',function(e){
    var step=e.shiftKey?48:16;
    if(e.key==='ArrowLeft'){setPanel(panelWidth()+step);e.preventDefault();}
    else if(e.key==='ArrowRight'){setPanel(panelWidth()-step);e.preventDefault();}
    else if(e.key==='Home'){setPanel(340);e.preventDefault();}});
})();
// A CAPTURE THE READER CAN KEEP, WHEREVER THE PAGE IS EMBEDDED.
//
// Save Image builds a Blob, clicks a <a download>, revokes the URL and then
// writes "Saved PNG: ...x..., N dpi" to its status line -- it never checks
// whether anything was saved. In an iframe without allow-downloads (a chat
// artifact tile, a docs embed, a notebook output cell) the click is dropped in
// silence, so the reader is told the file is on disk when no file exists
// anywhere. That is the worst possible failure for a figure someone is about
// to put in a talk.
//
// Measured, on this bundle, clicking Save in the capture panel:
//
//   top-level file:// tab                      downloads (344,768 B)
//   plain <iframe>                             downloads
//   sandboxed <iframe>, no allow-downloads     BLOCKED, nothing written
//   sandboxed <iframe> + allow-downloads       downloads
//
// So the block is real but narrow, and the first version of this bar warned
// about it unconditionally -- telling three readers out of four that a save
// which had just succeeded might not have. Worse, the escape hatch it offered
// was "right-click the link", which an Electron/WebView2 shell does not have a
// context menu for, on a blob: URL from a file:// page (opaque origin) that is
// the least reliable thing to Save-link-as even where it does.
//
// So: COPY comes first. The clipboard needs neither a download permission nor
// a context menu, which is precisely what an embedded app shell withholds, and
// PNG carries the alpha channel through. The link stays, upgraded to a data:
// URL that survives right-click and drag-out. The warning appears only when
// the page is actually embedded. Revoking is deferred rather than skipped, so
// the URL stays valid long enough to use and is still collected.
(function(){
  if(!window.URL||!URL.createObjectURL)return;
  var makeUrl=URL.createObjectURL.bind(URL),dropUrl=URL.revokeObjectURL.bind(URL),lastBlob=null;
  URL.createObjectURL=function(blob){lastBlob=blob;return makeUrl(blob);};
  URL.revokeObjectURL=function(url){setTimeout(function(){try{dropUrl(url);}catch(e){}},600000);};
  var made=document.createElement.bind(document);
  document.createElement=function(tag){
    var el=made(tag);
    if(String(tag).toLowerCase()!=='a')return el;
    var realClick=el.click.bind(el);
    el.click=function(){
      var name=el.getAttribute&&el.getAttribute('download');
      if(name&&lastBlob){
        // Rename BEFORE the click, or the file the reader actually receives
        // keeps the upstream project's prefix and the bar's copy is a second,
        // differently-named download of the same picture.
        var nice=exportName(name);
        try{el.setAttribute('download',nice);}catch(e){}
        showCapture(nice,lastBlob);}
      try{realClick();}catch(e){}};
    return el;};
  function exportName(name){
    // py2Dmol names its own exports; after the rename that prefix is somebody
    // else's project on the file that lands in the reader's Downloads folder.
    try{
      var el=document.getElementById('pinsp-inspection-state');
      var st=JSON.parse(el.textContent.split('<\/').join('</'));
      var src=(st.source&&st.source.path)||'';
      var base=src.split('/').pop().split('\\').pop().replace(/\.(cif|mmcif)$/i,'');
      if(base)return base+name.replace(/^.*?(_\d{4}-\d\d-\d\dT[\d-]+)?(\.[a-z0-9]+)$/i,'$1$2');
    }catch(e){}
    return String(name).replace(/^py2dmol_/,'');}
  function asDataUrl(blob,cb){
    try{var r=new FileReader();
      r.onload=function(){cb(String(r.result));};
      r.onerror=function(){cb(null);};
      r.readAsDataURL(blob);}catch(e){cb(null);}}
  function showCapture(name,blob){
    var box=document.getElementById('pinsp-capture');if(!box)return;
    var url=makeUrl(blob),mb=(blob.size/1048576).toFixed(2),nice=name;
    box.textContent='';
    var line=document.createElement('span');
    line.textContent='Capture ready \u2014 '+mb+' MB, transparent background. ';
    box.appendChild(line);
    var note=document.createElement('div');
    note.className='bp-muted';note.style.marginTop='6px';
    function say(t){note.textContent=t;}
    // COPY FIRST: the only route that needs neither a download permission nor
    // a context menu. PNG on the clipboard keeps the alpha channel.
    var copy=document.createElement('button');
    copy.type='button';copy.className='bp-btn';copy.textContent='Copy image';
    // LEGACY PATH FIRST WHEN THE MODERN ONE IS ABSENT OR REFUSED. Selecting an
    // <img> in a contenteditable and running execCommand('copy') happens
    // synchronously inside the click, so it needs no permission grant and no
    // promise -- which is exactly the ground the async API loses on inside an
    // embedded frame.
    function copyViaSelection(){
      try{
        var holder=document.createElement('div');
        holder.contentEditable='true';
        holder.style.cssText='position:fixed;left:-9999px;top:0;opacity:0';
        var im=document.createElement('img');im.src=url;
        holder.appendChild(im);document.body.appendChild(holder);
        var rng=document.createRange();rng.selectNode(im);
        var sel=window.getSelection();sel.removeAllRanges();sel.addRange(rng);
        var ok=document.execCommand('copy');
        sel.removeAllRanges();document.body.removeChild(holder);
        return ok;
      }catch(e){return false;}}
    copy.addEventListener('click',function(){
      var settle=function(text,msg){
        copy.textContent=text;copy.disabled=false;say(msg||'');
        if(text==='Copied')setTimeout(function(){copy.textContent='Copy image';},4000);};
      copy.disabled=true;copy.textContent='Copying\u2026';
      var viaApi=window.ClipboardItem&&navigator.clipboard&&navigator.clipboard.write;
      if(!viaApi){
        settle(copyViaSelection()?'Copied':'Copy image',
          copyViaSelection?'':'This browser has no clipboard image support \u2014 '
            +'use Download, or right-click the picture below.');
        return;}
      try{
        navigator.clipboard.write([new ClipboardItem({'image/png':blob})]).then(
          function(){settle('Copied');},
          function(err){
            if(copyViaSelection()){settle('Copied');return;}
            settle('Copy image','Clipboard refused ('+((err&&err.name)||'unknown')
              +'). This usually means the frame is not allowed to write to the '
              +'clipboard \u2014 open the file in a browser tab, or drag the picture '
              +'below onto your desktop.');});
      }catch(e){
        if(copyViaSelection()){settle('Copied');return;}
        settle('Copy image','Clipboard refused ('+((e&&e.name)||'unknown')+').');}});
    box.appendChild(copy);
    var link=document.createElement('a');
    link.className='bp-btn';link.href=url;link.textContent='Download';
    link.setAttribute('download',nice);
    box.appendChild(link);
    // A data: URL right-clicks and drags out where a blob: from an opaque
    // origin does not. Swapped in once it is ready; the blob URL works until.
    asDataUrl(blob,function(d){if(d)link.href=d;});
    var nm=document.createElement('span');
    nm.className='bp-muted';nm.textContent=nice;
    box.appendChild(nm);
    box.appendChild(note);
    // ONLY WARN WHERE IT CAN BE TRUE. A top-level page can always download.
    if(window.top!==window.self)
      say('This page is embedded in a frame. If Download does nothing the frame '
        +'blocks downloads \u2014 use Copy image, or open the file in a browser tab.');
    if(/^image\//.test(blob.type||'')){
      var img=document.createElement('img');img.src=url;img.alt=nice;box.appendChild(img);}
    box.hidden=false;}
})();
// PARTNERS ARE DRAWN, NOT ANNOTATED, AND THEY BELONG TO A STATE. Every
// conformation's partners sit in the same object as their own block of
// positions; only the current conformation's are shown. So changing state
// swaps the partner too -- the tethered receptor does not keep wearing the
// Fab that was bound to the extended one, which would be a composite passed
// off as an observation. `owner` maps each partner chain to the conformation
// it came from; everything not in it is the target and is always visible.
var partners=(s.partners&&s.partners.owner)?s.partners:null,partnersOn=true;
function applyPartners(showFor){
  var state=(showFor===undefined)?morphAt:showFor;
  var r=renderer();if(!r||!partners)return false;
  var list=r.chains||[];if(!list.length)return false;
  var pos=new Set(),chs=new Set(),shown=0;
  for(var i=0;i<list.length;i++){
    var of=partners.owner[list[i]];
    if(of!==undefined){
      if(!partnersOn||of!==state)continue;
      shown++;}
    pos.add(i);
    if(typeof r.chainKeyAt==='function')chs.add(r.chainKeyAt(i));}
  // THE BUTTON IS INDEPENDENT OF THE CONFORMATION. A visibility patch sends
  // the renderer back to the first frame, so showing or hiding a partner was
  // an unasked-for jump to the reference state. `morphAt` -- not whatever the
  // renderer currently holds -- is the authority on which conformation is on
  // screen, and the frame index of conformation k is k. Re-assert it after
  // the patch AND on the next animation frame, because the drop does not
  // always land inside the setVisibility call.
  var keep=(morph&&morphReady)?morphAt
    :((typeof r.currentFrame==='number'&&r.currentFrame>=0)?r.currentFrame:0);
  if(keep<0)keep=0;
  r.setVisibility(chs.size?{positions:pos,chains:chs}:{positions:pos});
  function holdFrame(){
    var live=renderer();if(!live)return;
    if(typeof live.setFrame==='function'&&live.currentFrame!==keep)live.setFrame(keep);
    else live.render('Protein Inspector partners');}
  holdFrame();
  requestAnimationFrame(holdFrame);
  return shown;}
function hasPartnersHere(){
  if(!partners)return false;
  for(var c in partners.owner)if(partners.owner[c]===morphAt)return true;
  return false;}
// THE BUTTON NAMES WHAT IS ON SCREEN. Partners belong to a state, so the
// word for them does too: the apo trimer's are other protomers, the bound
// state's is a receptor. `labels` carries one per conformation; the global
// label is the fallback for a state that did not name its own.
function partnerLabelAt(st){
  var m=partners&&partners.labels;
  if(m){var v=m[String(st)];if(v)return v;}
  return (partners&&partners.label)||'partners';}
function syncPartnerButton(){
  var btn=document.getElementById('pinsp-partners');if(!btn||!partners)return;
  var here=hasPartnersHere(),name=partnerLabelAt(morphAt);
  btn.disabled=!here;
  if(here)btn.textContent=name;
  btn.setAttribute('aria-pressed',String(!!(here&&partnersOn)));
  btn.title=here?(partnersOn?'Hide ':'Show ')+name
    :'No partners in this conformation';}
var partnerBtn=document.getElementById('pinsp-partners');
if(partnerBtn)partnerBtn.addEventListener('click',function(){
  if(partnerBtn.disabled)return;
  partnersOn=!partnersOn;applyPartners();syncPartnerButton();});
// MORPH, INTERPOLATED ON DEMAND. The file carries one frame per conformation
// and nothing between them: the in-between coordinates are a straight line, so
// the page works them out as it draws. That is what makes ANY pair reachable
// directly -- a chain of stored intermediates can only be walked in order, so
// going from the first state to the last had to travel through every state in
// between, which is a claim about a pathway the data does not make. Cartesian
// interpolation is a depiction of two endpoints; the frames between them are
// not physical and bond geometry is not preserved there.
//
// One scratch frame, appended past the conformations, is the animation buffer.
// _loadFrameData re-reads object.frames[i] on every setFrame with no caching,
// so rewriting that one frame's coords and asking for it again is the whole
// mechanism.
var morph=(s.morph&&s.morph.mode==='browser')?s.morph:null,
    morphSteps=(morph&&morph.animation_steps)||18,morphAt=0,morphBusy=false,morphReady=false,
    morphScratch=-1;
function morphNow(){return (window.performance&&window.performance.now)
  ?window.performance.now():Date.now();}
function morphObject(){var r=renderer();if(!r)return null;
  var name=r.currentObjectName||'prepared-target';
  return (r.objectsData&&r.objectsData[name])||null;}
function expandMorph(){
  if(!morph||morphReady)return;
  var r=renderer(),obj=morphObject(),n=morph.conformers.length;
  if(!r||!obj||!obj.frames||obj.frames.length<n)return;
  if(obj.frames.length===n){
    var seed={},base=obj.frames[0];
    for(var key in base)seed[key]=base[key];
    seed.coords=base.coords.map(function(p){return [p[0],p[1],p[2]];});
    seed.pae=undefined;
    obj.frames.push(seed);
    if(typeof r.updateUIControls==='function')r.updateUIControls();}
  morphScratch=n;morphReady=true;setMorphButtons(morphAt);
  applyPartners();syncPartnerButton();}
function setMorphButtons(active){
  // [data-conf] MATTERS. The partners control wears .bm-btn too, so that it
  // looks like it belongs beside the ring -- and a bare .bm-btn query swept
  // it into the ring's own handlers. Number(null) is 0, so pressing
  // "partners" also asked to morph to the first conformation.
  var bs=document.querySelectorAll('.bm-btn[data-conf]');
  for(var i=0;i<bs.length;i++)bs[i].setAttribute('aria-pressed',
    String(Number(bs[i].getAttribute('data-conf'))===active));
  var bar=document.getElementById('pinsp-morph');
  if(bar)bar.setAttribute('data-busy',morphBusy?'1':'0');
  var note=document.getElementById('pinsp-morph-note');
  if(note&&morph){var c=morph.conformers[active];
    note.textContent=(c&&c.rmsd_to_reference_A)
      ?c.rmsd_to_reference_A+' \u00c5 C\u03b1 RMSD from '+morph.conformers[0].label
      :'reference';}}
function goMorph(target){
  var r=renderer();if(!morph||!r||morphBusy||target===morphAt)return;
  expandMorph();if(!morphReady)return;
  var obj=morphObject(),from=obj.frames[morphAt],to=obj.frames[target],
      buf=obj.frames[morphScratch];
  if(!from||!to||!buf)return;
  // THE SIDE CHAINS TRAVEL WITH THE STATE BEING LEFT. Each real frame has its
  // own side-chain table -- coefficients in the residue's backbone frame, so
  // they follow interpolated coordinates without a transform -- and the buffer
  // borrows whichever one the animation starts from rather than always the
  // reference's. The rotamer is frozen for the length of the animation, which
  // is the same thing the interpolation itself is: a depiction of two
  // endpoints, not a pathway. Landing on the target's own frame picks up the
  // target's own rotamers.
  buf.sidechains=from.sidechains;
  morphBusy=true;setMorphButtons(morphAt);
  // The outgoing partner goes the moment the target starts moving: leaving it
  // on through the animation shows it bound to coordinates it was never
  // solved against.
  applyPartners(-1);syncPartnerButton();
  var a=from.coords,b=to.coords,out=buf.coords,
      t0=morphNow(),dur=Math.max(320,Math.min(900,morphSteps*36));
  (function step(){
    var p=Math.min(1,(morphNow()-t0)/dur),e=p<0.5?2*p*p:-1+(4-2*p)*p;
    for(var k=0;k<out.length&&k<a.length;k++){
      out[k][0]=a[k][0]+(b[k][0]-a[k][0])*e;
      out[k][1]=a[k][1]+(b[k][1]-a[k][1])*e;
      out[k][2]=a[k][2]+(b[k][2]-a[k][2])*e;}
    r.setFrame(morphScratch);
    if(p<1){requestAnimationFrame(step);return;}
    // Land on the conformation's own frame, so the viewer is showing a real
    // structure and not the buffer once the animation stops.
    r.setFrame(target);
    morphAt=target;morphBusy=false;setMorphButtons(target);
    applyPartners();syncPartnerButton();})();}
// A WATCHDOG, BECAUSE THE FRAME MOVES BEHIND OUR BACK. Applying a visibility
// patch makes py2Dmol jump the frame -- sometimes inside the call, sometimes
// on a later tick, and sometimes to the LAST frame (the interpolation buffer)
// rather than the first. Chasing each path with a targeted restore kept
// missing one, so this watches the renderer's own frame-change event and puts
// the frame back whenever it disagrees with the conformation the reader
// chose. setFrame re-fires the event, but the guard makes the correction
// idempotent, so it settles in one pass.
document.addEventListener('py2dmol-frame-change',function(){
  if(!morph||!morphReady||morphBusy)return;
  var r=renderer();if(!r||typeof r.setFrame!=='function')return;
  if(r.currentFrame!==morphAt)r.setFrame(morphAt);});
var morphBtns=document.querySelectorAll('.bm-btn[data-conf]');
for(var mb=0;mb<morphBtns.length;mb++)(function(btn){
  btn.addEventListener('click',function(){goMorph(Number(btn.getAttribute('data-conf')));});
})(morphBtns[mb]);
var tabs=document.querySelectorAll('[data-tab]');
for(var t=0;t<tabs.length;t++)(function(btn){btn.addEventListener('click',function(){
  var want=btn.getAttribute('data-tab');
  for(var i=0;i<tabs.length;i++)tabs[i].setAttribute('aria-selected',
    tabs[i].getAttribute('data-tab')===want?'true':'false');
  var panes=document.querySelectorAll('[data-pane]');
  for(var j=0;j<panes.length;j++)panes[j].hidden=(panes[j].getAttribute('data-pane')!==want);
  // The canvas is sized on layout; coming back from a hidden pane needs a nudge.
  if(want==='structure'){if(pinspFit)pinspFit();
    var r=renderer();if(r)r.render('Protein Inspector tab shown');}});})(tabs[t]);
if(partners&&!morph){var pt=0;(function pwait(){
  if(applyPartners()!==false||++pt>600){syncPartnerButton();return;}
  requestAnimationFrame(pwait);})();}
var tries=0;(function wait(){var ok=syncVisibleLayers(),seqOk=false,foldOk=foldStylePanel();
  installDetailKeeper();tameFocus();
  if(ok){expandMorph();enableClicks();seqOk=seqBuild();if(seqOk)installSeqBridge();}
  if((ok&&seqOk&&foldOk&&(!morph||morphReady))||++tries>600)return;
  requestAnimationFrame(wait);})();
})();</script>"""
    panel = (css
             + '<script id="pinsp-inspection-state" type="application/json">__STATE__</script>'
             + script).replace("__STATE__", payload)

    counts = {item["layer_id"]: len(item["annotation_ids"]) for item in state["layers"]}
    layers = "".join(
        '<div class="bp-layer">'
        f'<label><input type="checkbox"{" checked" if item.get("visible", True) else ""} data-layer="{html.escape(item["layer_id"], quote=True)}">'
        f'<span class="bp-sw" style="background:{item["color"]}"></span>'
        f'<span class="bp-lab" title="{html.escape(item["label"], quote=True)}">{html.escape(item["label"])}</span>'
        '</label>'
        f'<span class="bp-n">{counts[item["layer_id"]]}</span>'
        f'<button type="button" class="bp-only" data-only="{html.escape(item["layer_id"], quote=True)}">only</button>'
        # SELECTING THE LAYER IS WHAT CONNECTS THE TWO HALVES OF THE PAGE: the
        # side-chain and neighbourhood tools act on a selection, and a layer
        # is already a named set of residues. Hover-revealed beside `only`,
        # the same affordance in the same place.
        '<button type="button" class="bp-only" data-select='
        f'"{html.escape(item["layer_id"], quote=True)}" '
        'title="Select this layer\u2019s residues \u2014 shift-click to add to '
        'the selection">sel</button>'
        '</div>'
        for item in state["layers"]
    )

    by_id = {a["annotation_id"]: a for a in state["annotations"]}
    groups = []
    for item in state["layers"]:
        inner = "".join(
            f'<div class="pinsp-residue{"" if by_id[aid]["resolved"] else " pinsp-warning"}" '
            f'data-annotation="{html.escape(aid, quote=True)}">{html.escape(by_id[aid]["label"])}</div>'
            for aid in item["annotation_ids"] if aid in by_id
        )
        groups.append(
            f'<details class="bp-group" data-group="{html.escape(item["layer_id"], quote=True)}">'
            f'<summary><span class="bp-sw" style="background:{item["color"]}"></span>'
            f'{html.escape(item["label"])} <span class="bp-n">{counts[item["layer_id"]]}</span></summary>'
            f'{inner}</details>'
        )
    rows = "".join(groups)

    # ONE BUTTON PER CONFORMATION, not a scrub bar. A slider asks the reader to
    # find the endpoints; a button says where it is going and gets there.
    # THE TOGGLE SITS BY THE PICTURE, not in the Layers panel: partners are
    # part of what is drawn, not an annotation over it.
    partners = state.get("partners") or {}
    partnerbar = ""
    if partners.get("owner"):
        # The opening state is the reference, so the button opens wearing the
        # reference's word for its partners; syncPartnerButton swaps it with
        # the conformation from there.
        opening = (partners.get("labels") or {}).get("0") or partners.get("label") or "Partners"
        partnerbar = ('<button type="button" class="bm-btn bm-solo" id="pinsp-partners" '
                      f'aria-pressed="true">{html.escape(opening)}'
                      '</button>')

    morph = state.get("morph") or {}
    morphbar = ""
    if morph.get("mode") == "browser":
        buttons = "".join(
            '<button type="button" class="bm-btn" '
            f'data-conf="{index}" aria-pressed="{"true" if index == 0 else "false"}">'
            f'{html.escape(item["label"])}</button>'
            for index, item in enumerate(morph.get("conformers", []))
        )
        morphbar = ('<div class="pinsp-morph" id="pinsp-morph" data-busy="0">'
                    f'<div class="bm-ring" role="group" aria-label="Conformation">{buttons}</div>'
                    f'{partnerbar}'
                    '<span class="bp-muted" id="pinsp-morph-note"></span></div>')
    elif partnerbar:
        morphbar = f'<div class="pinsp-morph" id="pinsp-morph">{partnerbar}</div>'
    stage_has_ring = morph.get("mode") == "browser"

    about = state.get("about") or []
    tabstrip = ""
    panes_open, panes_close = "", ""
    if about:
        buttons = ['<button type="button" class="pinsp-tab" data-tab="structure" '
                   'aria-selected="true">Structure</button>']
        buttons += [f'<button type="button" class="pinsp-tab" data-tab="about-{i}" '
                    f'aria-selected="false">{html.escape(t["title"])}</button>'
                    for i, t in enumerate(about)]
        tabstrip = f'<nav class="pinsp-tabs">{"".join(buttons)}</nav>'
        panes_open = '<div data-pane="structure">'
        panes_close = "</div>" + "".join(
            f'<div data-pane="about-{i}" hidden><article class="pinsp-about">{t["body_html"]}</article></div>'
            for i, t in enumerate(about))

    credit = (
        '<!--\n  This file is free to reuse. Both notices below are BEER-WARE\n'
        '  (Revision 42) and both ask the same one thing: keep them here.\n\n'
        f'  Viewer: py2Dmol by Sergey Ovchinnikov -- {UPSTREAM_REPOSITORY}\n'
        f'  at revision {UPSTREAM_REVISION}, inlined in this file.\n'
        f'  {UPSTREAM_LICENSE}\n\n'
        f'  Annotation layers, morph and export: Protein Inspector\n'
        f'  {INSPECTOR_VERSION} -- {INSPECTOR_REPOSITORY}\n'
        f'  {INSPECTOR_LICENSE}\n-->'
    )
    # ONE LINE, ON THE BOTTOM EDGE OF THE SEQUENCE SECTION.
    #
    # This was four sentences in a band of its own below the strip, which left
    # a stripe of empty page under the only zone that already fills its width.
    # The long form said three things twice: that the viewer is inlined (it is
    # a self-contained file -- there is nowhere else it could be), which parts
    # are whose (the repository links say it), and where the full notice lives
    # (it is a comment at the top of the file anyone reading this far can
    # find). What is left is the two names, the two versions, and the licence.
    #
    # The rest goes in `title`, so it is one hover away rather than gone --
    # and the notice itself is not weakened: the BEER-WARE comment at the top
    # of the file is the copy that travels.
    credit_line = (
        '<footer class="pinsp-credit" title="Both are BEER-WARE'
        ' (Revision 42) \u2014 free to reuse, keep the notice. The full text of'
        ' both is in the comment at the top of this file.">'
        f'<a href="{UPSTREAM_REPOSITORY}">py2Dmol</a> by Sergey Ovchinnikov'
        f' rev&nbsp;{UPSTREAM_REVISION[:7]} \u00b7 '
        f'<a href="{INSPECTOR_REPOSITORY}">Protein Inspector</a> '
        # ...the NUMBER, not the slug. INSPECTOR_VERSION is
        # "protein-inspector-1.14", which after the link text read "Protein
        # Inspector protein-inspector-1.14 by Tadas Kluonis". The full slug
        # still travels in the comment at the top of the file and in the
        # manifest, where it is an identifier rather than a label.
        f'{INSPECTOR_VERSION.removeprefix("protein-inspector-")}'
        ' by Tadas Kluonis \u00b7 BEER-WARE (Rev&nbsp;42)'
        '</footer>'
    )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        f'<title>Protein Inspector</title></head><body>{credit}'
        f'{tabstrip}{panes_open}'
        '<main class="protein-inspector">'
        '<section class="pinsp-stage" id="pinsp-stage" data-controls="1" '
        f'data-morph="{"1" if stage_has_ring else "0"}">'
        f'{morphbar}{viewer_html}'
        '<div class="pinsp-capture" id="pinsp-capture" hidden></div>'
        '</section>'
        '<div class="pinsp-split" id="pinsp-split" role="separator" '
        'aria-orientation="vertical" tabindex="0" '
        'aria-label="Resize the panel"></div>'
        '<aside class="pinsp-panel">'
        '<div class="bp-head"><h2>Layers</h2><div class="bp-btns">'
        '<button type="button" class="bp-btn" id="pinsp-layers-all">All</button>'
        '<button type="button" class="bp-btn" id="pinsp-layers-none">None</button>'
        # ...and the whole legend as one selection, which is what "show the
        # side chains for these" means when several layers are ticked. The
        # header's buttons act on every layer, so this belongs beside them
        # rather than as a fourth control on the sequence bar.
        '<button type="button" class="bp-btn" id="pinsp-layers-select" '
        'title="Select every residue in the ticked layers">Select</button>'
        '</div></div>'
        f'<div class="bp-layers">{layers}</div>'
        '<div class="bp-head"><h2>Residues</h2>'
        '<span class="bp-muted" id="pinsp-count"></span></div>'
        '<input id="pinsp-filter" type="search" placeholder="Filter residues \u2014 try 433, His, glycan">'
        '<div class="bp-card" id="pinsp-card" hidden></div>'
        '<pre id="pinsp-residue-details">Click a residue to see every annotation on it.</pre>'
        '<button type="button" class="bp-btn" id="pinsp-clear-residue" style="display:none">Clear</button>'
        f'<div class="bp-rows">{rows}</div>'
        '</aside>'
        # THE SEQUENCE, UNDER BOTH OF THEM, with its own seam above it. The
        # three tools sit on this bar rather than in the Layers panel: what
        # they act on is the SELECTION, which is made here and in the picture,
        # and neither of those is a layer.
        '<div class="pinsp-hsplit" id="pinsp-hsplit" role="separator" '
        'aria-orientation="horizontal" tabindex="0" aria-valuenow="152" '
        'aria-label="Resize the sequence strip"></div>'
        '<section class="pinsp-seq" id="pinsp-seq">'
        '<div class="pinsp-seqbar"><h2>Sequence</h2>'
        '<span class="bp-muted" id="pinsp-sel">'
        'Click the structure or a letter to select a residue.</span>'
        '<span class="pinsp-tools">'
        '<label for="pinsp-around-a">within</label>'
        '<input id="pinsp-around-a" type="number" min="0.5" max="30" step="0.5" '
        'value="5" aria-label="Distance in angstroms">'
        '<label for="pinsp-around">&#197;</label>'
        '<button type="button" class="bp-btn" id="pinsp-around" disabled '
        'title="Add every residue with an atom within this distance of the '
        'selection">Select around</button>'
        # A SHOW/HIDE PAIR, NOT A LATCH. One question, two buttons, and the
        # state on their faces -- which is how src/app/selection.js puts it
        # upstream, and the reason is the THIRD state: a selection where some
        # residues draw their side chains and some do not. A single latch can
        # only render that as "off", so it said what it would do and left the
        # reader to work out what it had done. `Show` fills when every selected
        # residue draws, `Hide` when none does, neither when they disagree.
        '<span class="bp-muted">side chains</span>'
        '<span class="pinsp-pair" role="group" aria-label="Side chains of the selection">'
        '<button type="button" class="bp-btn pinsp-sw" id="pinsp-sc-on" disabled '
        'aria-pressed="false">Show</button>'
        '<button type="button" class="bp-btn pinsp-sw" id="pinsp-sc-off" disabled '
        'aria-pressed="false">Hide</button>'
        '</span>'
        '<button type="button" class="bp-btn" id="pinsp-sel-clear" disabled>'
        'Clear</button>'
        '</span></div>'
        '<div class="pinsp-seqbody" id="pinsp-seqbody"></div>'
        # ...AND THE CREDIT ON ITS BOTTOM EDGE, so the sequence section reaches
        # the foot of the page instead of floating above a band of nothing. It
        # is `flex:none` inside the section's column, which means the strip
        # keeps every pixel the seam gives it minus one line -- and the seam
        # arithmetic is untouched, because the row it measures is the section.
        f'{credit_line}'
        '</section>'
        '</main>'
        # THE CREDIT IS INSIDE THE SEQUENCE SECTION NOW, above this line -- and
        # therefore still inside the STRUCTURE PANE, which is the part that
        # mattered. Emitted here, after `panes_close`, it lived outside every
        # pane and reappeared under each custom About tab: a footer about a
        # viewer, stamped under prose that has nothing to do with the viewer.
        f'{panes_close}'
        f'{panel}</body></html>'
    )

def render_inspection_bundle(
    *,
    mmcif_path: str,
    mmcif_sha256: str,
    inspection_manifest: dict[str, Any],
    output_dir: str,
    output_name: str = "inspection",
    display_options: dict[str, Any] | None = None,
    extras: bool = False,
    conformers: list[dict[str, Any]] | None = None,
    partner_chains: list[str] | None = None,
    partner_label: str = "Partners",
    layers_on_partners: bool | str = "auto",
    morph_mapping: str = "exact",
    morph_reference_label: str = "reference",
    morph_to: list[dict[str, Any]] | None = None,
    morph_steps: int = 12,
) -> dict[str, Any]:
    """Render ONE self-contained HTML bundle (plus, on request, audit files).

    A bundle is a single file on purpose. It gets emailed, dropped in Slack and
    opened on a machine that has none of this checked out, so anything that
    lives beside it is a thing that arrives detached or not at all. Put extra
    context in an `about` tab (see `_about_tabs`) rather than in a second file:
    `about` takes a string, a {title: body} mapping, or a list of
    {"title", "body"}, and each entry becomes a tab next to "Structure".

    Pass `extras=True` to also write `.viewer.json`, `.svg`, `.png` and
    `.manifest.json`. Those are for auditing your own annotations -- the viewer
    state is where a wrong `canonical_position` shows up -- not for the reader.

    Bundles above SIZE_TARGET_BYTES (20 MB) still render, and the returned
    manifest carries a `size_warning`.
    """
    source = Path(mmcif_path).resolve()
    if not source.is_file() or source.suffix.lower() not in {".cif", ".mmcif"}:
        raise ValueError("mmcif_path must identify an existing local .cif/.mmcif file")
    actual_hash = _sha256(source)
    if actual_hash != mmcif_sha256:
        raise ValueError("mmCIF hash does not match the immutable artifact reference")
    manifest = _validate_manifest(dict(inspection_manifest))
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    base = _safe_output(root, output_name)
    options = dict(display_options or {})
    allowed_options = {"width", "height", "style", "color", "background", "chain_palette",
                       "sidechains", "gpu", "detail"}
    unknown_options = set(options) - allowed_options
    if unknown_options:
        raise ValueError(f"unsupported display options: {sorted(unknown_options)}")
    capture_sidechains = bool(options.get("sidechains", True))

    viewer = py2Dmol.view(
        size=(int(options.get("width", 720)), int(options.get("height", 720))),
        style=options.get("style", "tube"),
        # Chain identity in neutral greys (see chainPaletteFor in core/mol.js), so
        # hue stays free for the annotation layers. Recoloured from the viewer's
        # own Style panel rather than from a bespoke control in the side panel.
        color=options.get("color", "chain"),
        chain_palette=options.get("chain_palette", "greys"),
        bg=options.get("background", "white"),
        # THE GPU PAINTER, BECAUSE THE BUNDLE ALREADY CARRIES IT.
        #
        # This was False, and it is what made a `richardson` or `cartoon`
        # bundle feel slow next to py2dmol.solab.org on the same machine: the
        # website's default is gpu=True (viewer.py's DEFAULT_CONFIG, line 63)
        # and this was handing the same geometry to the CPU painter. What that
        # costs is structural rather than a benchmark -- paintgl.js's own
        # header states it: the drawing is "resident on the GPU so that turning
        # the model costs one draw call instead of one full repaint", where the
        # 2D painter runs cartoon/geom.js every frame (core/mol.js, line
        # 10276). No numbers are quoted here because none were measured on this
        # machine, and upstream does not publish a head-to-head.
        #
        # NOTHING IS GIVEN UP, because the notebook bundle is the one build
        # with BOTH painters in it (see BUNDLES in tools/bundle.py) and
        # core/mol.js only honours this flag when both are present:
        #   - No WebGL2, or a lost context: `_gpuWillTake` declines the frame
        #     and the 2D painter draws it. Per frame, not per session.
        #   - SVG export: paintgl refuses any context carrying
        #     getSerializedSvg, so vector output still comes off the 2D
        #     painter, which is the only one that can produce it.
        #   - `protein-inspector capture`: headless Chromium has no GPU, and
        #     capture.py already launches it with SwiftShader for exactly this
        #     reason.
        # `display={"gpu": False}` forces the CPU painter back.
        gpu=bool(options.get("gpu", True)),
        # CARTOON SAMPLING, RAISED TO THE TOP OF ITS RANGE.
        #
        # py2Dmol's default is 4, and 4 is a reasonable default for a viewer
        # someone is orbiting at whole-protein zoom. It is the wrong default
        # for a bundle, and the difference is the point of a bundle: this file
        # exists to be zoomed into and read residue by residue, and `detail`
        # is the one control whose cost was deliberately made independent of
        # zoom. Upstream retired an adaptive term that used to target a fixed
        # on-screen chord length -- it was the largest term in the frame and it
        # made a picture depend on the zoom it was first built at -- so
        # sampling is now exactly what this number says at every zoom, and
        # cartoon/geom.js states the consequence plainly: "Magnified curves
        # facet rather than resample."
        #
        # At 4 a residue gets 4 subdivisions in a helix (HELIX_SUB 8) and 3 in
        # a strand or loop (SHEET_SUB / SUB, both 6), scaled by detail/8. Close
        # in, that is three flat plates per residue of loop, each carrying its
        # own outline stroke -- read by the first person to zoom a bundle in as
        # malformed side-chain bonds, which is how this was found. 8 is the top
        # of the range and reproduces the sampling the retired adaptive term
        # used to reach (~8.7 per residue on a 600px canvas).
        #
        # Affordable here precisely because the GPU painter is on above: the
        # doubled station count is paid once, at capture, rather than on every
        # frame of every rotation. The Detail slider in the Style
        # panel still spans 2-8 live, and `display={"detail": 4}` restores
        # py2Dmol's own default.
        detail=min(8, max(2, int(options.get("detail", 8)))),
        controls=True,
        # CARRY THE SIDE-CHAIN ATOMS, so a reader who clicks a residue can ask
        # to see it. Without this flag py2Dmol ships a Cα trace, the table that
        # `showSidechains` needs is never built, and the verb raises -- which
        # is the one thing a reader cannot work around in an exported file.
        #
        # The atoms are stored once per residue, not per frame, so the cost is
        # a function of the structure and not of the number of conformations.
        # `display={"sidechains": False}` declines it; the Side chains button
        # then says why it is disabled rather than pretending the capability is
        # there. See CHANGELOG for the measured cost on a real structure.
        sidechains=capture_sidechains,
    )
    # INLINE THE LIBRARY IN EVERY EXPORT, ALWAYS.
    #
    # py2Dmol's default is notebook semantics: the first view() of a process
    # writes the ~1 MB library and offers it over a BroadcastChannel, and every
    # later view() writes a short borrow stub instead. That is right for cells
    # in one document and wrong for a file on disk. A bundle is opened on its
    # own, often on another machine, with no lender anywhere on the page — the
    # stub then fails with "the viewer library never loaded. Re-run the first
    # cell in this notebook that created a viewer", which is advice that cannot
    # be followed for an exported file.
    #
    # The failure is silent at render time and depends on how many bundles the
    # process rendered before this one, so only the FIRST export of a session
    # worked. tests/test_renderer.py::test_every_export_is_self_contained
    # renders twice in one process and fails if the second borrows.
    viewer._share_library = False
    morph_report = None
    partner_state = None
    if conformers and morph_to:
        raise ValueError("pass either conformers= or morph_to=, not both")
    if conformers:
        # ENDPOINTS ONLY. The intermediates are a straight line between two
        # coordinate sets, so shipping them is shipping a number the reader's
        # own machine can work out: N conformers cost N frames here and the
        # page expands them to N + (N-1)*steps on load. At 600 residues that
        # is the difference between 15 frames and 3 in the file.
        # A PARTNER BELONGS TO A STATE. The reference's partner is where the
        # reference put it; a conformer that brings its own has it somewhere
        # else entirely. So each conformation's partners are their own block of
        # positions, present in every frame and revealed only while that
        # conformation is the one on screen -- which is the only honest way to
        # show a receptor whose bound partner differs between states.
        partner_spec = [list(partner_chains or ())]
        for entry in conformers:
            partner_spec.append(list(entry.get("partner_chains") or ()))
        excluded = sorted({c for group in partner_spec for c in group})
        keys, names, stacks, labels, morph_report, fits = _conformer_stacks(
            source, conformers, morph_mapping, morph_reference_label,
            exclude_chains=excluded)
        chain_ids = [k[0] for k in keys]
        residue_numbers = [k[1] for k in keys]
        paths = [source] + [Path(e["path"]).resolve() for e in conformers]
        owner, taken, blocks = {}, set(chain_ids), []
        partner_sidechains: list[list[Any]] = []
        for index, (path, group) in enumerate(zip(paths, partner_spec)):
            if not group:
                continue
            trace = _ca_trace(path)
            wanted = set(group)
            extra = sorted(k for k in trace if k[0] in wanted)
            if not extra:
                continue
            # Two conformations may both carry a chain "B", and one object
            # cannot hold the same (chain, number) address twice -- the
            # annotation lookup and the selection both key on it. Rename on
            # collision rather than silently merging two different molecules.
            renamed = {}
            for chain in sorted(wanted):
                label = chain
                suffix = 1
                while label in taken:
                    label = f"{chain}{suffix}"
                    suffix += 1
                taken.add(label)
                renamed[chain] = label
                owner[label] = index
            # THE PARTNER MOVES WITH ITS RECEPTOR. Conformer i's target was
            # rotated onto the reference; its partner has to take the same
            # trip or it is left floating where the original crystal put it.
            block = fits[index](np.array([trace[k][1] for k in extra], dtype=float))
            # ...and so do its side chains, under the same transform and at the
            # offset this block lands on. Read before chain_ids grows, because
            # that length IS the offset.
            if capture_sidechains:
                partner_sidechains += _sidechain_rows(path, extra, len(chain_ids),
                                                      fits[index])
            stacks = [np.vstack([stack, block]) for stack in stacks]
            chain_ids += [renamed[k[0]] for k in extra]
            residue_numbers += [k[1] for k in extra]
            names += [trace[k][0] for k in extra]
            blocks.append({"conformer": index, "label": labels[index],
                           "chains": [renamed[c] for c in sorted(wanted)],
                           "renamed_from": {renamed[c]: c for c in sorted(wanted)
                                            if renamed[c] != c} or None,
                           "n_positions": len(extra)})
        if blocks:
            morph_report["partner_blocks"] = blocks
            # ONE LABEL PER STATE. The button names what is on screen, and
            # what is on screen changes with the conformation: the apo trimer
            # shows its other protomers, the bound state shows a receptor.
            # One global string made the button say "other protomers" over a
            # receptor, which is a caption contradicting its own picture.
            per_state = {0: partner_label}
            for index, entry in enumerate(conformers, start=1):
                per_state[index] = entry.get("partner_label") or partner_label
            partner_state = {"label": partner_label, "owner": owner,
                             "labels": {str(k): v for k, v in per_state.items()
                                        if k in set(owner.values())},
                             "per_conformer": True}
        # ONE SIDE-CHAIN ROW SET PER CONFORMATION. `add_pdb` gets these for
        # free; this path builds its frames from Cα coordinates, so without
        # them a morph bundle could show the sequence, select, and grow a
        # neighbourhood, and the Side chains button was dead -- on exactly the
        # bundles that are most worth interrogating. See `_sidechain_rows` for
        # why the table survives both the fit and the interpolation: it is
        # coefficients in each residue's own backbone frame, rebuilt per frame.
        #
        # The partner rows are shared, because a partner block holds the same
        # coordinates in every frame; only the target's rotamers differ.
        sidechain_rows: list[list[Any] | None] = [None] * len(stacks)
        if capture_sidechains:
            sidechain_rows = [
                (_sidechain_rows(path, keys, 0, fits[index]) + partner_sidechains) or None
                for index, path in enumerate(paths)
            ]
        # ...AND THE TYPES, WHICH IS WHAT MAKES THE SIDE CHAINS DRAW AS ATOMS.
        #
        # For a protein trace this says exactly what core/mol.js's own default
        # would have said ('P' everywhere) -- and sending it anyway is
        # load-bearing, because the default is reached by a DIFFERENT route. _materialiseSidechains appends one position per
        # side-chain atom and types each 'L', building the new array as
        # `(data.position_types || []).slice()` plus one push per atom. With no
        # types in the frame that starts EMPTY, so the finished array is as
        # long as the appended atoms rather than as long as the coordinates;
        # _setDataField takes a per-position array whose length does not match
        # the coordinate count and silently replaces it with
        # `Array(n).fill('P')`. Every side-chain atom then claimed to be a
        # protein backbone position, and cartoon/geom.js partitions on exactly
        # that -- backbone drawn as cartoon, everything else as sticks -- so
        # the atoms were swept into the ribbon and drawn as slabs running
        # between them. The same bundle rendered through add_pdb was right,
        # because add_pdb sends types.
        #
        # Reported as "the boxes connecting the atoms look wrong in between",
        # and invisible to every check that had been run: the side-chain table
        # is chemically exact (559/559 residues, every ring closure) and its
        # geometry survives the fit intact (CA-CB median 1.530 A in both
        # frames). Nothing was wrong with the atoms. They were mislabelled.
        position_types = _position_types(names)
        for index, coords in enumerate(stacks):
            viewer.add(np.asarray(coords, dtype=float), chains=chain_ids,
                       residue_numbers=residue_numbers, position_names=names,
                       position_types=position_types,
                       name="prepared-target", align=(index == 0),
                       allow_reflection=False,
                       sidechain_atoms=sidechain_rows[index])
        morph_report["mode"] = "browser"
        morph_report["animation_steps"] = int(morph_steps)
        morph_report["frames_in_file"] = len(stacks)
        morph_report["stored_intermediates"] = 0
    elif morph_to:
        # A MORPH IS ONE OBJECT WITH MANY FRAMES, which is what keeps the
        # colouring still. Layer colour is written per POSITION onto the object
        # (see _color_annotations), so it is shared by every frame -- the
        # structure moves and the annotation stays on the residue it names.
        keys, names, frames, morph_report = _morph_frames(source, morph_to, morph_steps)
        chain_ids = [k[0] for k in keys]
        residue_numbers = [k[1] for k in keys]
        # Types here too. This path carries no side-chain atoms, so nothing is
        # appended and the array cannot fall out of step -- but a frame that
        # describes its own positions is the contract the other paths keep, and
        # the asymmetry is what hid the bug above for three rounds.
        position_types = _position_types(names)
        for index, coords in enumerate(frames):
            viewer.add(np.asarray(coords, dtype=float), chains=chain_ids, residue_numbers=residue_numbers,
                       position_names=names, position_types=position_types,
                       name="prepared-target",
                       align=(index == 0), allow_reflection=False)
    else:
        viewer.add_pdb(str(source), use_biounit=False, filter_additives=False, load_ligands=True, name="prepared-target")
        if partner_chains:
            # No conformations to follow, so every partner belongs to the one
            # state there is. Same button, same visibility patch.
            partner_state = {"label": partner_label,
                             "owner": {chain: 0 for chain in partner_chains},
                             "labels": {"0": partner_label},
                             "per_conformer": False}
    if not any(item.get("frames") for item in viewer.objects):
        raise ValueError("py2Dmol could not load a renderable structure")
    layers = _layers(manifest["annotations"])
    highlight = manifest.get("highlight", "color")
    # THE OTHER PROTOMERS WEAR THE ANNOTATION TOO. A copy of the target is
    # the target: hiding it is still one button, but while it is on screen it
    # carries the same colours rather than sitting next to them in grey.
    mirror_chains: list[str] = []
    if partner_state:
        mirror_chains = _partner_copies(viewer.objects[-1]["frames"][0],
                                        set(partner_state["owner"]),
                                        layers_on_partners)
        partner_state["painted"] = mirror_chains or None
    if highlight == "color":
        _color_annotations(viewer, manifest["annotations"], layers,
                           base_color=manifest.get("base_color", BASE_COLOR),
                           paint_base=manifest.get("base_mode", "chain") == "custom",
                           mirror_chains=mirror_chains)
    # WHETHER THE ATOMS ARE ACTUALLY IN THE FILE, which is not the same as
    # whether they were asked for: the conformer and morph paths build their
    # frames from Cα coordinates, so there is nothing to capture there however
    # the flag is set. The page reads this to decide whether the Side chains
    # button is a control or an explanation.
    has_sidechains = any(
        bool(frame.get("sidechain_atoms"))
        for item in viewer.objects for frame in (item.get("frames") or ())
    )
    state = {
        "sidechains": has_sidechains,
        "highlight": highlight,
        "base_color": manifest.get("base_color", BASE_COLOR),
        "base_mode": manifest.get("base_mode", "chain"),
        "schema_version": "protein-inspector-state-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream_revision": UPSTREAM_REVISION,
        "source": {"path": str(source), "sha256": actual_hash},
        "layers": layers,
        "annotations": manifest["annotations"],
        "overlays": manifest.get("overlays", []),
        "about": _about_tabs(manifest.get("about")),
        "morph": morph_report,
        "partners": partner_state,
        "viewer": {"config": viewer.config, "objects": viewer.objects},
    }
    html_path = base.with_suffix(".html")
    state_path = base.with_suffix(".viewer.json")
    svg_path = base.with_suffix(".svg")
    png_path = base.with_suffix(".png")
    html_path.write_text(_html(viewer._display_viewer(static_data=viewer.objects), state))
    wanted = [("html", html_path, "text/html")]
    if extras:
        state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
        coords, chains, numbers = _trace(source)
        xy = _project(coords)
        _snapshot_svg(svg_path, xy, chains, numbers)
        _snapshot_png(png_path, xy, chains)
        wanted += [
            ("viewer_state", state_path, "application/json"),
            ("svg", svg_path, "image/svg+xml"),
            ("png", png_path, "image/png"),
        ]

    artifacts = []
    for kind, path, media_type in wanted:
        artifacts.append({"kind": kind, "path": str(path), "sha256": _sha256(path), "size_bytes": path.stat().st_size, "media_type": media_type})
    compact = {
        "schema_version": "protein-inspector-artifact-manifest-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream": {"repository": "https://github.com/sokrypton/py2Dmol", "revision": UPSTREAM_REVISION, "license": "BEER-WARE Revision 42"},
        "source_sha256": actual_hash,
        "layer_count": len(state["layers"]),
        "annotation_count": len(state["annotations"]),
        "about_tabs": [item["title"] for item in state.get("about") or []],
        "morph": morph_report,
        "partners": ({"labels": partner_state.get("labels"),
                      "chains": sorted(partner_state["owner"]),
                      "painted": partner_state.get("painted")}
                     if partner_state else None),
        "artifacts": artifacts,
    }
    html_bytes = html_path.stat().st_size
    if html_bytes > SIZE_TARGET_BYTES:
        compact["size_warning"] = (
            f"bundle is {html_bytes / 1024 / 1024:.1f} MB, above the "
            f"{SIZE_TARGET_BYTES // 1024 // 1024} MB target for something that gets emailed; "
            "consider fewer annotations or a smaller structure"
        )
    if extras:
        compact_path = base.with_suffix(".manifest.json")
        compact_path.write_text(json.dumps(compact, indent=2, sort_keys=True) + "\n")
        compact["manifest"] = {"path": str(compact_path), "sha256": _sha256(compact_path), "size_bytes": compact_path.stat().st_size, "media_type": "application/json"}
    return compact


def read_inspection_bundle(html_path: str) -> dict[str, Any]:
    """Recover the full inspection state from a rendered bundle. The inverse of render.

    A bundle is one self-contained HTML file, which is right for a reader and
    unhelpful for a program unless the machine-readable part is a supported
    surface rather than something to scrape. It is: the whole state -- layers,
    annotations, residue coordinates, About bodies, the source digest -- is one
    JSON object in `<script id="pinsp-inspection-state">`, and this reads it
    back without a browser, an HTML parser, or any third-party package.

    Returns a dict with:
      inspector_version, schema_version, source {path, sha256}, highlight,
      base_mode, base_color
      layers       [{layer_id, label, color, visible, n_annotations}]
      about        [{title, body_html}]
      annotations  the manifest's annotations, verbatim
      residues     one flat row per modelled residue -- chain_id,
                   author_residue_number, canonical_position, residue_name,
                   x, y, z, plddt, color, layers[], labels[] -- which is the
                   table most callers actually want and is CSV-ready as-is.

    >>> state = read_inspection_bundle("inspection.html")
    >>> [r for r in state["residues"] if "ph_anchor" in r["layers"]]
    """
    text = Path(html_path).read_text(encoding="utf-8")
    match = re.search(
        r'<script id="pinsp-inspection-state" type="application/json">(.*?)</script>',
        text, re.S)
    if not match:
        raise ValueError(f"{html_path} carries no protein_inspector inspection state")
    # `_html` escapes "</" as "<\/" so the payload cannot close its own tag.
    state = json.loads(match.group(1).replace("<\\/", "</"))

    by_residue: dict[tuple[str, int], dict[str, Any]] = {}
    for item in state.get("annotations", []):
        if not item.get("resolved"):
            continue
        res = item["residue"]
        key = (str(res["chain_id"]), int(res["author_residue_number"]))
        slot = by_residue.setdefault(key, {"layers": [], "labels": [], "color": None})
        slot["layers"].append(item.get("layer_id") or item["kind"])
        slot["labels"].append(item["label"])
        if item.get("color"):
            slot["color"] = item["color"]      # last layer wins, as on screen

    # A COPY WEARS THE ANNOTATION IT WAS PAINTED WITH. The annotations name
    # one chain, but `layers_on_partners` may have painted the same colours
    # onto the other protomers of a homo-oligomer, and this table is meant to
    # report what the reader sees rather than what the manifest said. Mirror
    # the rows the same way the colour pass did; a residue absent from a copy
    # simply never appears in the loop below.
    partners = state.get("partners") or {}
    owners = set(partners.get("owner") or {})
    for copy in partners.get("painted") or ():
        for (chain, number), slot in list(by_residue.items()):
            if chain in owners:
                continue
            by_residue.setdefault((str(copy), number), slot)

    frames = (state.get("viewer", {}).get("objects") or [{}])[0].get("frames") or [{}]
    frame = frames[0]
    chains = frame.get("chains") or []
    numbers = frame.get("residue_numbers") or []
    coords = frame.get("coords") or []
    names = frame.get("position_names") or []
    plddts = frame.get("plddts") or []

    residues = []
    for index in range(len(chains)):
        key = (str(chains[index]), int(numbers[index]))
        extra = by_residue.get(key, {})
        xyz = coords[index] if index < len(coords) else [None, None, None]
        residues.append({
            "chain_id": key[0],
            "author_residue_number": key[1],
            "canonical_position": index + 1,
            "residue_name": names[index] if index < len(names) else None,
            "x": xyz[0], "y": xyz[1], "z": xyz[2],
            "plddt": plddts[index] if index < len(plddts) else None,
            "color": extra.get("color"),
            "layers": extra.get("layers", []),
            "labels": extra.get("labels", []),
        })

    counts: dict[str, int] = {}
    for item in state.get("annotations", []):
        key = item.get("layer_id") or item["kind"]
        counts[key] = counts.get(key, 0) + 1
    layers = [dict(layer, n_annotations=counts.get(layer["layer_id"], 0))
              for layer in state.get("layers", [])]
    for layer in layers:
        layer.pop("annotation_ids", None)

    return {
        "inspector_version": state.get("inspector_version"),
        "schema_version": state.get("schema_version"),
        "source": state.get("source"),
        "highlight": state.get("highlight"),
        "base_mode": state.get("base_mode"),
        "base_color": state.get("base_color"),
        "layers": layers,
        "about": state.get("about", []),
        "annotations": state.get("annotations", []),
        "residues": residues,
    }


#: Residue names whose trace atom is C4' rather than CA -- `_ca_trace` accepts
#: either, so a conformer stack can hold nucleic positions, and py2Dmol types
#: them 'D' (DNA) or 'R' (RNA) rather than 'P'. Getting this wrong is not
#: cosmetic: core/mol.js picks localFrame's step range from the type, and a base
#: rebuilt through the peptide's 3.8 A range fails the frame and has its atoms
#: dropped in silence -- 347 of them, on the one occasion upstream hit it.
_DNA_RESIDUES = frozenset({"DA", "DC", "DG", "DT", "DU", "DI"})
_RNA_RESIDUES = frozenset({"A", "C", "G", "U", "I"})


def _position_types(names: list[str]) -> list[str]:
    """One py2Dmol position type per trace position, from the residue name.

    'P' for protein, 'D' for DNA, 'R' for RNA. Sent explicitly rather than left
    to the renderer's default -- see the long note at the conformer `add` call
    for why a type array that is merely ABSENT is not the same as one that says
    'P', once side-chain atoms are appended to it.
    """
    out = []
    for name in names:
        key = (name or "").strip().upper()
        out.append("D" if key in _DNA_RESIDUES
                   else "R" if key in _RNA_RESIDUES else "P")
    return out


def _ca_trace(cif_path: Path) -> dict[tuple[str, int], tuple[str, list[float]]]:
    """{(chain, author_resnum): (residue_name, xyz)} for CA/C4' only.

    Raises on a duplicate key: an insertion code or an altloc that collapses two
    residues onto one address makes the residue mapping ambiguous, and a morph
    built on an ambiguous mapping silently interpolates the wrong pairs.
    """
    structure = MMCIFParser(QUIET=True).get_structure(cif_path.stem, str(cif_path))
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError(f"{cif_path.name}: a morph conformer needs exactly one model")
    out: dict[tuple[str, int], tuple[str, list[float]]] = {}
    for chain in models[0]:
        for residue in chain:
            atom = residue.child_dict.get("CA") or residue.child_dict.get("C4'")
            if atom is None or residue.id[0] != " ":
                continue
            key = (chain.id, int(residue.id[1]))
            if key in out:
                raise ValueError(
                    f"{cif_path.name}: residue {key} appears twice, so the mapping between "
                    "conformers is ambiguous; renumber or split the file before morphing")
            out[key] = (residue.get_resname(), [float(v) for v in atom.coord])
    return out


def _sidechain_rows(cif_path: Path, keys: list[tuple[str, int]], offset: int = 0,
                    transform=None) -> list[list[Any]]:
    """`[[position_index, residue_name, [[atom, x, y, z, element], ...]], ...]`.

    WHY A MORPH CAN CARRY SIDE CHAINS AT ALL. The browser does not store these
    atoms in world space: `buildSidechainTable` (src/io/sidechains.js) keeps
    each one as three coefficients in its residue's OWN backbone frame, built
    from the Cα trace, and rebuilds that frame from the final positions at draw
    time. So the table is invariant under translation and rotation, and
    `parts/ui.js` builds a SEPARATE one per frame from that frame's backbone --
    "building once and inheriting would freeze the rotamer of whichever frame
    carried the atoms". One row set per conformer is therefore both necessary
    and sufficient: each conformation shows its own rotamers.

    THE FIT HAS TO BE APPLIED. The coefficients are measured against the
    coordinates in the frame, and a conformer's frame holds its Kabsch-fitted
    trace while its file holds the original crystal coordinates. Mixing the two
    would measure every side chain against a backbone that has rotated out from
    under it, so `transform` is the same callable the trace went through.

    Keyed by (chain, author number) so it lines up with `_ca_trace`, and read
    with gemmi so the atom selection is `py2Dmol`'s own rather than a second
    copy of the rule. A key gemmi does not have is skipped rather than guessed:
    the residue then simply has no side chain to draw, as glycine has none.
    """
    import gemmi

    from py2Dmol.viewer import _sidechain_atoms

    structure = gemmi.read_structure(str(cif_path))
    structure.remove_hydrogens()
    structure.remove_alternative_conformations()
    by_key: dict[tuple[str, int], Any] = {}
    for chain in structure[0]:
        for residue in chain:
            by_key.setdefault((str(chain.name), int(residue.seqid.num)), residue)
    rows: list[list[Any]] = []
    for index, key in enumerate(keys):
        residue = by_key.get((str(key[0]), int(key[1])))
        if residue is None:
            continue
        atoms = _sidechain_atoms(residue)
        if not atoms:
            continue
        if transform is not None:
            moved = transform(np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float))
            atoms = [[a[0], round(float(p[0]), 2), round(float(p[1]), 2),
                      round(float(p[2]), 2), a[4]] for a, p in zip(atoms, moved)]
        rows.append([offset + index, residue.name, atoms])
    return rows


def _kabsch(mobile: np.ndarray, target: np.ndarray):
    """Return the transform that moves `mobile` onto `target`. No reflection.

    Returned as a callable rather than as coordinates, because everything that
    travels with the mobile body -- a bound partner, a ligand -- has to move by
    the SAME transform. Applying the fit to the receptor and leaving its
    partner in the original frame leaves the partner floating in space beside
    a receptor that has rotated out from under it.
    """
    mobile_center, target_center = mobile.mean(0), target.mean(0)
    mc, tc = mobile - mobile_center, target - target_center
    v, _, wt = np.linalg.svd(mc.T @ tc)
    d = np.sign(np.linalg.det(v @ wt))
    rotation = v @ np.diag([1.0, 1.0, d]) @ wt

    def apply(points: np.ndarray) -> np.ndarray:
        return (np.asarray(points, dtype=float) - mobile_center) @ rotation + target_center
    return apply


def _superpose(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Kabsch, reflection forbidden. Returns `mobile` moved onto `target`."""
    return _kabsch(mobile, target)(mobile)


def _conformer_stacks(reference: Path, conformers: list[dict[str, Any]],
                      mapping: str = "exact", reference_label: str = "reference",
                      exclude_chains: list[str] | None = None):
    """Superpose any number of conformers onto a shared residue mapping.

    Returns (keys, residue_names, stacks, labels, report) where `stacks` holds
    one (n_residues, 3) array per conformer, reference first.

    `mapping="exact"` is the default and refuses anything but an identical
    (chain, author number) set across every file: a morph over a quietly
    intersected mapping looks just as smooth while interpolating the wrong
    pairs. `mapping="intersection"` opts in to the common residues -- which is
    what comparing a crystal chain against a full-length prediction needs --
    and the report then names how many residues each file lost, so a dropped
    epitope cannot pass unnoticed.
    """
    if mapping not in {"exact", "intersection"}:
        raise ValueError('morph_mapping must be "exact" or "intersection"')
    skip = set(exclude_chains or ())
    ref = {k: v for k, v in _ca_trace(reference).items() if k[0] not in skip}
    traces, labels, digests = [ref], [reference_label], [_sha256(reference)]
    for entry in conformers:
        path = Path(entry["path"]).resolve()
        digest = _sha256(path)
        if entry.get("sha256") and entry["sha256"] != digest:
            raise ValueError(f"{path.name}: morph conformer hash does not match")
        traces.append({k: v for k, v in _ca_trace(path).items() if k[0] not in skip})
        labels.append(entry.get("label") or path.stem)
        digests.append(digest)

    common = set(ref)
    for trace in traces[1:]:
        common &= set(trace)
    dropped = []
    for label, trace in zip(labels, traces):
        missing = len(set(trace) - common)
        dropped.append(missing)
        if mapping == "exact" and (missing or len(trace) != len(common)):
            only_here = sorted(set(trace) - common)[:6]
            raise ValueError(
                f"{label}: residue mapping is not exact -- {missing} residues are not "
                f"shared by every conformer (e.g. {only_here}). Prepare the files over "
                'the same residue range, or pass morph_mapping="intersection".')
    if len(common) < 4:
        raise ValueError("conformers share fewer than 4 residues; nothing to superpose")

    keys = sorted(common)
    names = [ref[k][0] for k in keys]
    base = np.array([ref[k][1] for k in keys], dtype=float)
    stacks, transforms = [base], [lambda points: np.asarray(points, dtype=float)]
    for trace in traces[1:]:
        moved = np.array([trace[k][1] for k in keys], dtype=float)
        fit = _kabsch(moved, base)
        stacks.append(fit(moved))
        transforms.append(fit)
    report = {
        "mapping": mapping,
        "n_residues": len(keys),
        "residue_range": [int(keys[0][1]), int(keys[-1][1])],
        "conformers": [
            {"label": label, "sha256": digest, "residues_dropped": drop,
             "rmsd_to_reference_A": round(float(np.sqrt(((stack - base) ** 2).sum(1).mean())), 2)}
            for label, digest, drop, stack in zip(labels, digests, dropped, stacks)
        ],
    }
    return keys, names, stacks, labels, report, transforms


def _morph_frames(reference: Path, conformers: list[dict[str, Any]], steps: int):
    """Pre-expanded Cartesian morph: every intermediate written into the file.

    Kept for callers that want the frames on disk. `conformers=` on
    `render_inspection_bundle` is the cheaper route -- it ships the endpoints
    only and the page interpolates on demand, which is the same picture for a
    fraction of the bytes. Cartesian interpolation is a depiction of the
    endpoints, NOT a pathway: intermediates are not physical and bond geometry
    is not preserved.
    """
    keys, names, stacks, labels, report, _ = _conformer_stacks(reference, conformers, "exact")
    frames, rmsds = [], []
    for index in range(len(stacks) - 1):
        start, end = stacks[index], stacks[index + 1]
        rmsds.append(float(np.sqrt(((start - end) ** 2).sum(1).mean())))
        for step in range(steps):
            t = step / float(steps)
            frames.append(start * (1.0 - t) + end * t)
    frames.append(stacks[-1])
    report = {"conformers": labels, "steps_between": steps, "n_frames": len(frames),
              "n_residues": len(keys), "endpoint_rmsd_A": [round(v, 2) for v in rmsds]}
    return keys, names, frames, report
