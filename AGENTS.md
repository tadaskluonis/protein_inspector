# BindOS Agent Guide

Use `bindos_structure_inspector.render_inspection_bundle` to create a local,
self-contained inspection bundle from a prepared single-model mmCIF file. Do
not fetch structures, edit JavaScript, or modify generated HTML for ordinary
figures.

## Workflow

1. Verify the mmCIF SHA-256 before rendering.
2. Build an inspection manifest with one annotation per residue-level claim.
3. Group arbitrary categories with `layer_id`; categories are user-defined.
4. Give each layer one `#RRGGBB` color and render the bundle.
5. Check the returned artifact hashes and inspect the generated HTML.

Every resolved annotation must identify the residue with `component_id`,
`canonical_position`, `chain_id`, and `author_residue_number`. Do not guess an
address: an unresolved claim should use `"resolved": false` instead.

## Layer Example

```python
from bindos_structure_inspector import render_inspection_bundle

manifest = {
    "schema_version": "bindos-inspection-manifest-1",
    "annotations": [
        {
            "annotation_id": "pocket-42",
            "kind": "candidate_hotspot",
            "label": "Predicted binding hotspot",
            "layer_id": "binding-hotspots",
            "layer_label": "Binding hotspots",
            "color": "#dc2626",
            "resolved": True,
            "residue": {
                "component_id": "target",
                "canonical_position": 42,
                "chain_id": "A",
                "author_residue_number": 42,
            },
            "evidence_ids": ["model-run-17"],
            "method": "classifier",
        },
    ],
}

result = render_inspection_bundle(
    mmcif_path="prepared.cif",
    mmcif_sha256="...",
    inspection_manifest=manifest,
    output_dir="inspection-output",
)
```

`kind` must be a supported scientific annotation kind; use `custom` when none
applies. `layer_id` controls the visualization category and may be any non-empty
string. Annotations in the same layer must use the same optional `layer_label`
and `color`.

## One file, with tabs (inspector 1.8)

`render_inspection_bundle` writes **one self-contained HTML file** and nothing
else. A bundle gets emailed, dropped in a channel and opened on a machine that
has none of this checked out, so anything sitting beside it arrives detached or
not at all.

**Extra context goes in a tab, never in a second file.** The manifest's
optional top-level `about` takes a string, a `{title: body}` mapping, or a list
of `{"title", "body"}`; each entry becomes a tab beside "Structure". Bodies are
freeform HTML — headings, tables, lists, links, code all pass through — so an
overview, a legend and a provenance note are three tabs, not three documents.
Bodies are sanitised (script/style/iframe/object/embed and their contents, all
`on*` handlers, and `javascript:`/`data:` URLs are stripped), because the text
is usually agent-authored and often quotes fetched sources. A body with no tags
at all is treated as plain text and blank-line blocks become paragraphs.

Keep a bundle under `SIZE_TARGET_BYTES` (20 MB) so it stays sendable; above it
the render still succeeds and the returned manifest carries a `size_warning`.

Pass `extras=True` only when you need the audit set — `.viewer.json`, `.svg`,
`.png`, `.manifest.json`. The viewer state is where a wrong
`canonical_position` shows up, so use it while checking your own annotations,
then re-render without it for the copy you hand over.

## The inspection panel (inspector 1.4)

The panel is emitted by `render_inspection_bundle`. Do not re-implement or
hand-edit any of it in generated HTML; it is covered by
`tests/bindos_inspection_dom.js`, which clicks the controls rather than merely
asserting their strings appear.

- **Layers** section, at the top: `All` / `None` buttons, then one row per layer
  with a colour swatch, the label, its annotation count, and a hover-revealed
  `only` button that isolates that layer.
- **Residues** section: a live count, a filter box (matches residue label, layer
  label and annotation id), and the residue list grouped into one collapsible
  `<details>` per layer rather than one flat list of every annotation.
- **Residue card**: clicking a row opens a formatted card — residue identity,
  chain / author number / model index, and one line per annotation on that
  residue with its layer swatch. It is dismissed by the card's `×`, by clicking
  the same row again, or by Escape; all three also drop the 3D selection.
  `#bindos-residue-details` still exists as a hidden plain-text mirror for
  programmatic checks.
- `window.bindosInspection` exposes `syncVisibleLayers`, `clearDetails`,
  `setAllLayers`, `onlyLayer` and `selectAnnotation`.

## Colour (inspector 1.7)

The structure is painted a neutral grey (`BASE_COLOR`) and layer colour is the
only hue that carries meaning. **Do not try to get a flat base out of
`display_options={"color": ...}`** — py2Dmol's modes are `auto`, `chain`,
`rainbow`, `plddt`, `deepmind`, `entropy`, `object`, `hydrophobicity`, `ss`;
there is no grey, and an unrecognised name is NOT rejected: `ui.js` falls back
to `auto`, which is rainbow on a single chain. The base therefore comes from
the per-position colour map, in `_color_annotations` and again in the page's
`syncVisibleLayers`, and it is present when every layer is off — dropping the
map to `null` there is what used to expose the rainbow.

Halo mode is the deliberate exception and paints no base, because it exists to
leave a pLDDT- or chain-coloured structure intact and annotate on top of it.

The default render style is **`tube`**. Pass `display_options={"style": ...}`
for `cartoon` / `richardson` / `ribbon` / `3d` when a specific figure needs it.

## One visual channel, not two

The manifest's optional top-level `"highlight"` picks how layers mark residues:

- `"color"` (default) paints the residues in their layer color. **Prefer this.**
- `"halo"` leaves the structure its default color and rings the layer's residues
  with the selection highlight instead.

Do not try to get both. Painting residues and haloing the same residues encodes
one fact twice, and the halo is the weaker signal of the two -- it is thin, it
reads poorly against a cartoon, and it competes with the color underneath.
`highlight` is bundle-wide for exactly this reason: there is no per-layer form,
because mixing the two channels across layers has the same problem. Reach for
`"halo"` only when the structure's own coloring is the subject -- a pLDDT or
chain-colored view you need to keep intact -- and annotate on top of it.

Clicking a residue row always halos that one residue, in either mode. That is a
transient "you are here" marker for a single residue, not a layer channel, and
it does not count as the second encoding this section warns about.

Layer colors are applied to the 3D residues. Unchecking a generated layer
checkbox removes that layer's color from the structure, hides its panel rows,
and drops it from the selection highlight; re-checking restores it. With every
layer unchecked the structure returns to its default coloring. If layers
overlap, the later annotation layer in the manifest is the visible residue
color -- and that is re-resolved on every toggle, so hiding the layer on top
reveals the one underneath rather than leaving a gap.

`residue.canonical_position` is the 1-based index of the residue within the
modeled chain, which is NOT the author numbering: for a chain whose first
residue is `21`, `author_residue_number: 98` is `canonical_position: 78`. Both
fields are required and a wrong `canonical_position` colors the wrong residue
without raising. The exported `inspection.viewer.json` is the place to check
your work -- `viewer.objects[0].color.value.position` maps 0-based object index
(`canonical_position - 1`) to color.

The exported `inspection.svg` and `inspection.png` are geometry-only projections
and carry no layer color. Use `inspection.html` for anything a reader will look
at.

`render_inspection_bundle` needs `gemmi` and `IPython` installed alongside
`biopython` and `numpy`.
