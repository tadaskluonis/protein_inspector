"""Render immutable local mmCIF plus a BindOS inspection manifest.

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


INSPECTOR_VERSION = "bindos-inspector-1.2"
MANIFEST_SCHEMA = "bindos-inspection-manifest-1"
UPSTREAM_REVISION = "78c2d489d0b5c5d19accd9eeeef878c2868f5271"

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
        if "color" in item and (not isinstance(item["color"], str) or not _HEX_COLOR.fullmatch(item["color"])):
            raise ValueError(f"annotation {item['annotation_id']} color must be a #RRGGBB value")
        if item["kind"] == "partner_contact" and not item.get("partner_id"):
            raise ValueError("partner_contact annotations require partner_id")
    # ONE visual channel, never two. Painting the residues and haloing the same
    # residues encodes the same fact twice and makes both harder to read, so
    # `highlight` chooses between them for the whole bundle.
    if value.get("highlight", "color") not in _HIGHLIGHT_STYLES:
        raise ValueError(f"highlight must be one of {sorted(_HIGHLIGHT_STYLES)}")
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
        layer = grouped.setdefault(key, {"annotation_ids": [], "label": item.get("layer_label"), "color": item.get("color")})
        if item.get("layer_label") and layer["label"] != item["layer_label"]:
            raise ValueError(f"layer {key!r} has conflicting labels")
        if item.get("color") and layer["color"] != item["color"]:
            raise ValueError(f"layer {key!r} has conflicting colors")
        layer["annotation_ids"].append(item["annotation_id"])
    return [
        {
            "layer_id": key,
            "label": layer["label"] or key.replace("_", " ").replace(":", " — "),
            "color": layer["color"] or _LAYER_COLORS[index % len(_LAYER_COLORS)],
            "annotation_ids": layer["annotation_ids"],
            "visible": True,
        }
        for index, (key, layer) in enumerate(grouped.items())
    ]


def _color_annotations(viewer: Any, annotations: list[dict[str, Any]], layers: list[dict[str, Any]]) -> None:
    """Apply each declared layer color to its resolved residue positions."""
    frame = viewer.objects[-1]["frames"][0]
    positions: dict[tuple[str, int], list[int]] = {}
    for index, (chain, number) in enumerate(zip(frame.get("chains") or [], frame.get("residue_numbers") or [])):
        positions.setdefault((str(chain), int(number)), []).append(index)
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
        if indices:
            viewer.set_color(layer["color"], position=indices)


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
    panel = """
<style>.bindos-inspector{font:14px system-ui;display:grid;grid-template-columns:minmax(0,1fr) 320px;gap:12px}.bindos-panel{padding:12px;border:1px solid #ddd;border-radius:8px;max-height:860px;overflow:auto}.bindos-layer{display:block;margin:5px 0}.bindos-residue{cursor:pointer;padding:2px 4px}.bindos-residue:hover{background:#fff7cc}.bindos-warning{color:#9a3412}#bindos-residue-details{white-space:pre-wrap;background:#f8fafc;padding:8px;border-radius:6px}</style>
<script id="bindos-inspection-state" type="application/json">__STATE__</script>
<script>(function(){const s=JSON.parse(document.getElementById('bindos-inspection-state').textContent);const byId=Object.fromEntries(s.annotations.map(a=>[a.annotation_id,a]));let applyingLayerSelection=false;function renderer(){var reg=window.py2dmol_viewers;if(!reg)return null;var id=(s.viewer&&s.viewer.config&&s.viewer.config.viewer_id)||(window.viewerConfig&&window.viewerConfig.viewer_id);var entry=id?reg[id]:null;if(!entry){var keys=Object.keys(reg);if(keys.length===1)entry=reg[keys[0]];}return entry&&entry.renderer?entry.renderer:null;}function sameResidue(a,b){return a&&b&&a.component_id===b.component_id&&a.copy_index===b.copy_index&&a.canonical_position===b.canonical_position;}function indicesFor(residue){const r=renderer(),out=[];if(!r||!residue)return out;for(let i=0;i<(r.residueNumbers||[]).length;i++){const chain=r.chains&&r.chains[i];if(chain===residue.chain_id&&r.residueNumbers[i]===residue.author_residue_number)out.push(i);}return out;}function showDetails(annotation){for(const x of document.querySelectorAll('[data-annotation]'))x.style.outline='';const related=s.annotations.filter(x=>sameResidue(x.residue,annotation.residue));for(const x of related){const match=document.querySelector('[data-annotation="'+CSS.escape(x.annotation_id)+'"]');if(match&&!match.hidden)match.style.outline='2px solid #eab308';}document.getElementById('bindos-residue-details').textContent=JSON.stringify({residue:annotation.residue,annotations:related.length?related:[annotation]},null,2);}function syncVisibleLayers(){const r=renderer(),selected=new Set();for(const l of s.layers){const box=document.querySelector('[data-layer="'+CSS.escape(l.layer_id)+'"]');const visible=!box||box.checked;for(const id of l.annotation_ids){const row=document.querySelector('[data-annotation="'+CSS.escape(id)+'"]');if(row)row.hidden=!visible;const a=byId[id];if(visible&&a&&a.resolved)for(const index of indicesFor(a.residue))selected.add(index);}}if(!r)return false;var halo=s.highlight==='halo';if(!halo){var name=r.currentObjectName,object=r.objectsData&&r.objectsData[name];if(object){var off=0;if(typeof r.localRangeOf==='function'){var win=r.localRangeOf(name);if(win&&typeof win.off==='number')off=win.off;}var paint={},painted=false;for(const l of s.layers){const box=document.querySelector('[data-layer="'+CSS.escape(l.layer_id)+'"]');if(box&&!box.checked)continue;for(const id of l.annotation_ids){const a=byId[id];if(!a||!a.resolved||!a.color)continue;for(const index of indicesFor(a.residue)){paint[index-off]=a.color;painted=true;}}}object.color=painted?{type:'advanced',value:{position:paint}}:null;r.colorsNeedUpdate=true;r.plddtColorsNeedUpdate=true;}}applyingLayerSelection=true;try{if(halo)r.setResidueSelection(selected);r.render('BindOS inspection layer visibility');}finally{applyingLayerSelection=false;}return true;}function selectInStructure(residue){const r=renderer();if(!r)return;applyingLayerSelection=true;try{r.setResidueSelection(new Set(indicesFor(residue)));r.render('BindOS manifest residue selection');}finally{applyingLayerSelection=false;}}for(const l of s.layers){const box=document.querySelector('[data-layer="'+CSS.escape(l.layer_id)+'"]');if(box)box.addEventListener('change',syncVisibleLayers);}for(const row of document.querySelectorAll('[data-annotation]'))row.addEventListener('click',()=>{const a=byId[row.dataset.annotation];showDetails(a);selectInStructure(a.residue);});document.addEventListener('py2dmol-residue-selection-change',()=>{if(applyingLayerSelection)return;const r=renderer();if(!r||!r.residueSelection||!r.residueSelection.size)return;const selected=s.annotations.find(a=>indicesFor(a.residue).some(i=>r.residueSelection.has(i)));if(selected)showDetails(selected);});window.bindosInspection={syncVisibleLayers:syncVisibleLayers,selectAnnotation:function(id){const a=byId[id];if(a){showDetails(a);selectInStructure(a.residue);}}};var tries=0;(function waitForViewer(){if(syncVisibleLayers()||++tries>600)return;requestAnimationFrame(waitForViewer);})();})();</script>
""".replace("__STATE__", payload)
    layers = "".join(
        f'<label class="bindos-layer"><input type="checkbox" checked data-layer="{html.escape(item["layer_id"], quote=True)}"> '
        f'<span style="color:{item["color"]}">&#9679;</span> {html.escape(item["label"])}'
        f'</label>'
        for item in state["layers"]
    )
    rows = "".join(
        f'<div class="bindos-residue{(" bindos-warning" if not item["resolved"] else "")}" data-annotation="{html.escape(item["annotation_id"], quote=True)}"><b>{html.escape(item["label"])}</b><br><small>{html.escape(str(item.get("method", "unknown")))} · {html.escape(", ".join(item.get("evidence_ids", [])) or "no evidence")}</small></div>'
        for item in state["annotations"]
    )
    return f'<!doctype html><html><head><meta charset="utf-8"><title>BindOS structure inspection</title></head><body><main class="bindos-inspector"><section>{viewer_html}</section><aside class="bindos-panel"><h2>Inspection layers</h2>{layers}<h2>Selected residue</h2><pre id="bindos-residue-details">Click an annotated residue to inspect numbering, evidence, and warnings.</pre><h2>Residues and warnings</h2>{rows}</aside></main>{panel}</body></html>'


def render_inspection_bundle(
    *,
    mmcif_path: str,
    mmcif_sha256: str,
    inspection_manifest: dict[str, Any],
    output_dir: str,
    output_name: str = "inspection",
    display_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Create HTML, viewer state, SVG/PNG snapshots, and a compact manifest."""
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
    allowed_options = {"width", "height", "style", "color", "background"}
    unknown_options = set(options) - allowed_options
    if unknown_options:
        raise ValueError(f"unsupported display options: {sorted(unknown_options)}")

    viewer = py2Dmol.view(
        size=(int(options.get("width", 720)), int(options.get("height", 720))),
        style=options.get("style", "cartoon"),
        color=options.get("color", "chain"),
        bg=options.get("background", "white"),
        gpu=False,
        controls=True,
    )
    viewer.add_pdb(str(source), use_biounit=False, filter_additives=False, load_ligands=True, name="prepared-target")
    if not any(item.get("frames") for item in viewer.objects):
        raise ValueError("py2Dmol could not load a renderable structure")
    layers = _layers(manifest["annotations"])
    highlight = manifest.get("highlight", "color")
    if highlight == "color":
        _color_annotations(viewer, manifest["annotations"], layers)
    state = {
        "highlight": highlight,
        "schema_version": "bindos-viewer-state-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream_revision": UPSTREAM_REVISION,
        "source": {"path": str(source), "sha256": actual_hash},
        "layers": layers,
        "annotations": manifest["annotations"],
        "overlays": manifest.get("overlays", []),
        "viewer": {"config": viewer.config, "objects": viewer.objects},
    }
    html_path = base.with_suffix(".html")
    state_path = base.with_suffix(".viewer.json")
    svg_path = base.with_suffix(".svg")
    png_path = base.with_suffix(".png")
    state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    html_path.write_text(_html(viewer._display_viewer(static_data=viewer.objects), state))
    coords, chains, numbers = _trace(source)
    xy = _project(coords)
    _snapshot_svg(svg_path, xy, chains, numbers)
    _snapshot_png(png_path, xy, chains)

    artifacts = []
    for kind, path, media_type in [
        ("html", html_path, "text/html"),
        ("viewer_state", state_path, "application/json"),
        ("svg", svg_path, "image/svg+xml"),
        ("png", png_path, "image/png"),
    ]:
        artifacts.append({"kind": kind, "path": str(path), "sha256": _sha256(path), "size_bytes": path.stat().st_size, "media_type": media_type})
    compact = {
        "schema_version": "bindos-inspection-artifact-manifest-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream": {"repository": "https://github.com/sokrypton/py2Dmol", "revision": UPSTREAM_REVISION, "license": "BEER-WARE Revision 42"},
        "source_sha256": actual_hash,
        "layer_count": len(state["layers"]),
        "annotation_count": len(state["annotations"]),
        "artifacts": artifacts,
    }
    compact_path = base.with_suffix(".manifest.json")
    compact_path.write_text(json.dumps(compact, indent=2, sort_keys=True) + "\n")
    compact["manifest"] = {"path": str(compact_path), "sha256": _sha256(compact_path), "size_bytes": compact_path.stat().st_size, "media_type": "application/json"}
    return compact
