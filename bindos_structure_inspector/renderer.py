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


INSPECTOR_VERSION = "bindos-inspector-1.6"
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

    return (
        '<!doctype html><html><head><meta charset="utf-8">'
        '<title>BindOS structure inspection</title></head><body>'
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
        style=options.get("style", "tube"),
        color=options.get("color", "gray"),
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
