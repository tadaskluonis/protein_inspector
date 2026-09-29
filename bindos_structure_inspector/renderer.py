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


INSPECTOR_VERSION = "bindos-inspector-1.12"
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


def _color_annotations(viewer: Any, annotations: list[dict[str, Any]], layers: list[dict[str, Any]],
                       base_color: str = BASE_COLOR, paint_base: bool = True) -> None:
    """Paint a neutral base over the whole chain, then each layer on top.

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
    css = (
        "<style>"
        ".bindos-inspector{font:14px/1.45 system-ui,sans-serif;display:grid;"
        "grid-template-columns:minmax(0,1fr) 360px;gap:14px;color:#0f172a}"
        ".bindos-panel{border:1px solid #e2e8f0;border-radius:10px;max-height:880px;overflow:auto;"
        "display:flex;flex-direction:column;background:#fff}"
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
        "#bindos-filter{margin:0 12px 8px;padding:5px 8px;font:12px system-ui;border:1px solid #cbd5e1;"
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
        ".bindos-residue{cursor:pointer;padding:2px 6px;border-radius:5px;font-size:12px;"
        "border-left:3px solid transparent}"
        ".bindos-residue:hover{background:#f1f5f9}"
        ".bindos-residue[data-on='1']{background:#fef9c3;border-left-color:#eab308}"
        ".bindos-warning{color:#9a3412}"
        ".bp-muted{font-size:11px;color:#64748b}"
        ".bindos-tabs{display:flex;gap:2px;border-bottom:1px solid #e2e8f0;margin:0 0 14px;font:14px system-ui,sans-serif}"
        ".bindos-tab{font:13px system-ui;padding:8px 16px;border:0;background:none;cursor:pointer;color:#64748b;border-bottom:2px solid transparent;margin-bottom:-1px}"
        ".bindos-tab:hover{color:#0f172a}"
        '.bindos-tab[aria-selected="true"]{color:#0f172a;font-weight:600;border-bottom-color:#2563eb}'
        ".bindos-about{max-width:52em;font:15px/1.6 system-ui,sans-serif;color:#0f172a;padding:0 4px 28px}"
        ".bindos-about h2{font-size:17px;margin:26px 0 6px;border-bottom:1px solid #e2e8f0;padding-bottom:4px}"
        ".bindos-about h3{font-size:15px;margin:20px 0 4px}"
        ".bindos-about table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px}"
        ".bindos-about td,.bindos-about th{border-bottom:1px solid #e2e8f0;padding:6px 8px;text-align:left;vertical-align:top}"
        ".bindos-about code{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:13px}"
        ".bindos-about .sw{display:inline-block;width:14px;height:14px;border-radius:3px;border:1px solid #cbd5e1;vertical-align:-2px}"
        "#bindos-residue-details{display:none}"
        "</style>"
    )
    script = r"""<script>(function(){
var root=document.getElementById('bindos-inspection-state');if(!root)return;
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
  var c=document.getElementById('bindos-count');
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
    r.render('BindOS inspection layer visibility');}finally{applying=false;}
  return true;}
function showDetails(a){
  var rows=document.querySelectorAll('[data-annotation]');
  for(var i=0;i<rows.length;i++)rows[i].removeAttribute('data-on');
  var related=s.annotations.filter(function(x){return sameResidue(x.residue,a.residue);});
  if(!related.length)related=[a];
  for(var j=0;j<related.length;j++){
    var el=document.querySelector('[data-annotation="'+CSS.escape(related[j].annotation_id)+'"]');
    if(el&&!el.hidden)el.setAttribute('data-on','1');}
  var d=document.getElementById('bindos-residue-details');
  if(d)d.textContent=a.label;
  var card=document.getElementById('bindos-card');if(!card)return;
  var res=a.residue||{};
  var title=(a.label||'').split(' ')[0]||res.author_residue_number;
  var items='';
  for(var k=0;k<related.length;k++){var x=related[k];var lay=layerOf[x.layer_id]||{};
    items+='<li><span class="bp-sw" style="background:'+esc(x.color||'#94a3b8')+'"></span><span>'+
      esc(lay.label||x.layer_id||x.kind)+'<span class="bp-meta">'+esc(x.label)+'</span></span></li>';}
  card.innerHTML='<div class="bp-card-head"><b>'+esc(title)+'</b>'+
    '<button class="bp-x" id="bindos-card-close" title="clear (Esc)">&times;</button></div>'+
    '<div class="bp-addr">chain '+esc(res.chain_id)+' &middot; author '+esc(res.author_residue_number)+
    ' &middot; model index '+esc(res.canonical_position)+'</div><ul>'+items+'</ul>';
  card.hidden=false;
  var xb=document.getElementById('bindos-card-close');
  if(xb)xb.addEventListener('click',function(ev){ev.stopPropagation();clearDetails();});}
function clearDetails(){
  var rows=document.querySelectorAll('[data-annotation]');
  for(var i=0;i<rows.length;i++)rows[i].removeAttribute('data-on');
  var d=document.getElementById('bindos-residue-details');if(d)d.textContent=PLACEHOLDER;
  var card=document.getElementById('bindos-card');if(card){card.hidden=true;card.innerHTML='';}
  selectedId=null;
  var r=renderer();
  if(r){applying=true;try{r.setResidueSelection(new Set());
    r.render('BindOS clear residue selection');}finally{applying=false;}}}
function selectInStructure(res){var r=renderer();if(!r)return;
  applying=true;try{r.setResidueSelection(new Set(indicesFor(res)));
    r.render('BindOS manifest residue selection');}finally{applying=false;}}
function setAllLayers(on){for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
    if(b)b.checked=on;}syncVisibleLayers();}
function onlyLayer(id){for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
    if(b)b.checked=(s.layers[i].layer_id===id);}syncVisibleLayers();}
for(var i=0;i<s.layers.length;i++){var b=box(s.layers[i].layer_id);
  if(b)b.addEventListener('change',syncVisibleLayers);}
var onlys=document.querySelectorAll('[data-only]');
for(var o=0;o<onlys.length;o++)(function(btn){btn.addEventListener('click',function(ev){
  ev.preventDefault();ev.stopPropagation();onlyLayer(btn.getAttribute('data-only'));});})(onlys[o]);
var rws=document.querySelectorAll('[data-annotation]');
for(var w=0;w<rws.length;w++)(function(row){row.addEventListener('click',function(){
  var id=row.dataset?row.dataset.annotation:row.getAttribute('data-annotation');
  if(selectedId===id){clearDetails();return;}
  selectedId=id;var a=byId[id];if(!a)return;showDetails(a);selectInStructure(a.residue);});})(rws[w]);
var bAll=document.getElementById('bindos-layers-all');
if(bAll)bAll.addEventListener('click',function(){setAllLayers(true);});
var bNone=document.getElementById('bindos-layers-none');
if(bNone)bNone.addEventListener('click',function(){setAllLayers(false);});
var bClear=document.getElementById('bindos-clear-residue');
if(bClear)bClear.addEventListener('click',clearDetails);
var filt=document.getElementById('bindos-filter');
if(filt)filt.addEventListener('input',function(){query=(filt.value||'').toLowerCase().trim();
  syncVisibleLayers();});
document.addEventListener('keydown',function(e){if(e.key==='Escape')clearDetails();});
document.addEventListener('py2dmol-residue-selection-change',function(){
  if(applying)return;var r=renderer();if(!r||!r.residueSelection||!r.residueSelection.size)return;
  var hit=s.annotations.find(function(a){return indicesFor(a.residue).some(function(i){
    return r.residueSelection.has(i);});});
  if(hit){selectedId=hit.annotation_id;showDetails(hit);}});
window.bindosInspection={syncVisibleLayers:syncVisibleLayers,clearDetails:clearDetails,
  setAllLayers:setAllLayers,onlyLayer:onlyLayer,
  selectAnnotation:function(id){var a=byId[id];if(a){selectedId=id;showDetails(a);
    selectInStructure(a.residue);}}};
// Colour of the un-annotated structure is the VIEWER's business, set from its
// own Style panel; there is no duplicate control here. The only thing this page
// does is honour a manifest that asked for a flat custom base, and stand down
// the moment the Style panel says otherwise -- otherwise the base would
// silently override whatever the user just picked.
document.addEventListener('py2dmol-color-change',function(){
  var r=renderer();if(!r||!r.colorMode||baseMode===r.colorMode)return;
  baseMode=r.colorMode;syncVisibleLayers();});
var tabs=document.querySelectorAll('[data-tab]');
for(var t=0;t<tabs.length;t++)(function(btn){btn.addEventListener('click',function(){
  var want=btn.getAttribute('data-tab');
  for(var i=0;i<tabs.length;i++)tabs[i].setAttribute('aria-selected',
    tabs[i].getAttribute('data-tab')===want?'true':'false');
  var panes=document.querySelectorAll('[data-pane]');
  for(var j=0;j<panes.length;j++)panes[j].hidden=(panes[j].getAttribute('data-pane')!==want);
  // The canvas is sized on layout; coming back from a hidden pane needs a nudge.
  if(want==='structure'){var r=renderer();if(r)r.render('BindOS tab shown');}});})(tabs[t]);
var tries=0;(function wait(){if(syncVisibleLayers()||++tries>600)return;
  requestAnimationFrame(wait);})();
})();</script>"""
    panel = (css
             + '<script id="bindos-inspection-state" type="application/json">__STATE__</script>'
             + script).replace("__STATE__", payload)

    counts = {item["layer_id"]: len(item["annotation_ids"]) for item in state["layers"]}
    layers = "".join(
        '<div class="bp-layer">'
        f'<label><input type="checkbox" checked data-layer="{html.escape(item["layer_id"], quote=True)}">'
        f'<span class="bp-sw" style="background:{item["color"]}"></span>'
        f'<span class="bp-lab" title="{html.escape(item["label"], quote=True)}">{html.escape(item["label"])}</span>'
        '</label>'
        f'<span class="bp-n">{counts[item["layer_id"]]}</span>'
        f'<button type="button" class="bp-only" data-only="{html.escape(item["layer_id"], quote=True)}">only</button>'
        '</div>'
        for item in state["layers"]
    )

    by_id = {a["annotation_id"]: a for a in state["annotations"]}
    groups = []
    for item in state["layers"]:
        inner = "".join(
            f'<div class="bindos-residue{"" if by_id[aid]["resolved"] else " bindos-warning"}" '
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

    about = state.get("about") or []
    tabstrip = ""
    panes_open, panes_close = "", ""
    if about:
        buttons = ['<button type="button" class="bindos-tab" data-tab="structure" '
                   'aria-selected="true">Structure</button>']
        buttons += [f'<button type="button" class="bindos-tab" data-tab="about-{i}" '
                    f'aria-selected="false">{html.escape(t["title"])}</button>'
                    for i, t in enumerate(about)]
        tabstrip = f'<nav class="bindos-tabs">{"".join(buttons)}</nav>'
        panes_open = '<div data-pane="structure">'
        panes_close = "</div>" + "".join(
            f'<div data-pane="about-{i}" hidden><article class="bindos-about">{t["body_html"]}</article></div>'
            for i, t in enumerate(about))

    return (
        '<!doctype html><html><head><meta charset="utf-8">'
        '<title>BindOS structure inspection</title></head><body>'
        f'{tabstrip}{panes_open}'
        '<main class="bindos-inspector">'
        f'<section>{viewer_html}</section>'
        '<aside class="bindos-panel">'
        '<div class="bp-head"><h2>Layers</h2><div class="bp-btns">'
        '<button type="button" class="bp-btn" id="bindos-layers-all">All</button>'
        '<button type="button" class="bp-btn" id="bindos-layers-none">None</button>'
        '</div></div>'
        f'<div class="bp-layers">{layers}</div>'
        '<div class="bp-head"><h2>Residues</h2>'
        '<span class="bp-muted" id="bindos-count"></span></div>'
        '<input id="bindos-filter" type="search" placeholder="Filter residues \u2014 try 433, His, glycan">'
        '<div class="bp-card" id="bindos-card" hidden></div>'
        '<pre id="bindos-residue-details">Click a residue to see every annotation on it.</pre>'
        '<button type="button" class="bp-btn" id="bindos-clear-residue" style="display:none">Clear</button>'
        f'<div class="bp-rows">{rows}</div>'
        '</aside></main>'
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
    allowed_options = {"width", "height", "style", "color", "background", "chain_palette"}
    unknown_options = set(options) - allowed_options
    if unknown_options:
        raise ValueError(f"unsupported display options: {sorted(unknown_options)}")

    viewer = py2Dmol.view(
        size=(int(options.get("width", 720)), int(options.get("height", 720))),
        style=options.get("style", "tube"),
        # Chain identity in neutral greys (see chainPaletteFor in core/mol.js), so
        # hue stays free for the annotation layers. Recoloured from the viewer's
        # own Style panel rather than from a bespoke control in the side panel.
        color=options.get("color", "chain"),
        chain_palette=options.get("chain_palette", "greys"),
        bg=options.get("background", "white"),
        gpu=False,
        controls=True,
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
    # worked. tests/test_bindos_renderer.py::test_every_export_is_self_contained
    # renders twice in one process and fails if the second borrows.
    viewer._share_library = False
    morph_report = None
    if morph_to:
        # A MORPH IS ONE OBJECT WITH MANY FRAMES, which is what keeps the
        # colouring still. Layer colour is written per POSITION onto the object
        # (see _color_annotations), so it is shared by every frame -- the
        # structure moves and the annotation stays on the residue it names.
        keys, names, frames, morph_report = _morph_frames(source, morph_to, morph_steps)
        chain_ids = [k[0] for k in keys]
        residue_numbers = [k[1] for k in keys]
        for index, coords in enumerate(frames):
            viewer.add(np.asarray(coords, dtype=float), chains=chain_ids, residue_numbers=residue_numbers,
                       position_names=names, name="prepared-target",
                       align=(index == 0), allow_reflection=False)
    else:
        viewer.add_pdb(str(source), use_biounit=False, filter_additives=False, load_ligands=True, name="prepared-target")
    if not any(item.get("frames") for item in viewer.objects):
        raise ValueError("py2Dmol could not load a renderable structure")
    layers = _layers(manifest["annotations"])
    highlight = manifest.get("highlight", "color")
    if highlight == "color":
        _color_annotations(viewer, manifest["annotations"], layers,
                           base_color=manifest.get("base_color", BASE_COLOR),
                           paint_base=manifest.get("base_mode", "chain") == "custom")
    state = {
        "highlight": highlight,
        "base_color": manifest.get("base_color", BASE_COLOR),
        "base_mode": manifest.get("base_mode", "chain"),
        "schema_version": "bindos-viewer-state-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream_revision": UPSTREAM_REVISION,
        "source": {"path": str(source), "sha256": actual_hash},
        "layers": layers,
        "annotations": manifest["annotations"],
        "overlays": manifest.get("overlays", []),
        "about": _about_tabs(manifest.get("about")),
        "morph": morph_report,
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
        "schema_version": "bindos-inspection-artifact-manifest-1",
        "inspector_version": INSPECTOR_VERSION,
        "upstream": {"repository": "https://github.com/sokrypton/py2Dmol", "revision": UPSTREAM_REVISION, "license": "BEER-WARE Revision 42"},
        "source_sha256": actual_hash,
        "layer_count": len(state["layers"]),
        "annotation_count": len(state["annotations"]),
        "about_tabs": [item["title"] for item in state.get("about") or []],
        "morph": morph_report,
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
    JSON object in `<script id="bindos-inspection-state">`, and this reads it
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
        r'<script id="bindos-inspection-state" type="application/json">(.*?)</script>',
        text, re.S)
    if not match:
        raise ValueError(f"{html_path} carries no bindos inspection state")
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


def _superpose(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Kabsch, reflection forbidden. Returns `mobile` moved onto `target`."""
    mc, tc = mobile - mobile.mean(0), target - target.mean(0)
    v, _, wt = np.linalg.svd(mc.T @ tc)
    d = np.sign(np.linalg.det(v @ wt))
    rotation = v @ np.diag([1.0, 1.0, d]) @ wt
    return mc @ rotation + target.mean(0)


def _morph_frames(reference: Path, conformers: list[dict[str, Any]], steps: int):
    """Interpolate between superposed conformers that share an exact residue mapping.

    Returns (keys, residue_names, frames, report). Cartesian interpolation is a
    depiction of the endpoints, NOT a pathway: intermediates are not physical
    and bond geometry is not preserved. It is here because seeing domain II
    swing out says more about why an epitope is or is not reachable than two
    static pictures side by side.

    The residue mapping must be EXACT -- same (chain, author number) set in
    every conformer. Anything else is refused rather than silently intersected,
    because a morph over a quiet intersection looks just as smooth while
    interpolating the wrong pairs.
    """
    ref = _ca_trace(reference)
    keys = sorted(ref)
    names = [ref[k][0] for k in keys]
    base = np.array([ref[k][1] for k in keys], dtype=float)

    stacks, labels = [base], ["reference"]
    for entry in conformers:
        path = Path(entry["path"]).resolve()
        digest = _sha256(path)
        if entry.get("sha256") and entry["sha256"] != digest:
            raise ValueError(f"{path.name}: morph conformer hash does not match")
        other = _ca_trace(path)
        if set(other) != set(ref):
            only_ref = sorted(set(ref) - set(other))[:6]
            only_other = sorted(set(other) - set(ref))[:6]
            raise ValueError(
                f"{path.name}: residue mapping is not exact -- "
                f"{len(set(ref) - set(other))} residues only in the reference "
                f"(e.g. {only_ref}), {len(set(other) - set(ref))} only here (e.g. {only_other}). "
                "Prepare both files over the same residue range before morphing.")
        stacks.append(_superpose(np.array([other[k][1] for k in keys], dtype=float), base))
        labels.append(entry.get("label") or path.stem)

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
