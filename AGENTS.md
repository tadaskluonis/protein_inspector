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

Layer colors are applied to the 3D residues. The generated layer checkboxes
control the current selection highlight and panel rows; they do not erase the
saved layer colors. If layers overlap, the later annotation layer in the
manifest is the visible residue color.
