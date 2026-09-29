---
name: bindos-inspector
description: "Build a single self-contained interactive HTML structure viewer with annotated, toggleable residue layers — epitopes, hotspots, glycans, domains, mutations, interface contacts — and optional button-driven morphing between any number of conformations. Use when visualizing a binder-design target, showing which residues an analysis picked out on a protein, sharing a structure figure someone must be able to rotate and click rather than read off a PNG, or checking a design's contacts against a reference structure. Renders inline in chat as an artifact; no external viewer needed."
---

# BindOS structure inspector

One structure, any number of named residue layers, one HTML file. The file is
self-contained — it opens with no server, no network and no install, and in
Claude Science it renders **inline** as an artifact, so the reader never leaves
the conversation.

Use it whenever an analysis produces a *set of residues* that someone has to
look at on the structure: a predicted epitope, hotspot residues, a glycan
shadow, conserved-vs-divergent positions between orthologues, contacts in a
docked pose, liabilities in a design.

## Setup

The engine is a separate package. Once per environment:

```bash
pip install git+https://github.com/profdocpizza/bindos-structure-inspector@v1.13
```

If the user has a local checkout instead, set `BINDOS_INSPECTOR_HOME` to it —
`inspector_engine()` also probes `~/code/bindos-structure-inspector`.
`gemmi` is needed only when you hand it a `.pdb`.

## Workflow

1. **Prepare one mmCIF per structure.** One model, one chain per file is
   simplest. Author residue numbering is what you will annotate against, so
   renumber to UniProt positions *before* rendering if that is the numbering
   your analysis speaks.
2. **Write the layer spec** — an ordered list, background first, since later
   layers paint over earlier ones.
3. **Render**, then `save_artifacts` the HTML and embed it inline.
4. **Read it back** with `inspection_table()` to verify the residues actually
   landed where you meant.

```python
report = inspect_structure(
    "target.cif",
    layers=[
        {"id": "domain-III", "label": "Domain III", "color": "#8fa8bd",
         "residues": range(310, 481)},
        {"id": "epitope", "label": "Cetuximab epitope", "color": "#dc2626",
         "residues": {356: "Q384 rim", 384: "core contact", 468: "hot spot"}},
    ],
    about={"Why this target": "...prose, HTML allowed..."},
    out="target_inspection.html",
)
report["path"], report["unmapped_residues"]
```

`residues` takes a list/range of author residue numbers, or a
`{number: note}` mapping when each residue deserves its own caption.
Always check `report["unmapped_residues"]` — residues absent from the model are
skipped silently otherwise.

## Morphing between conformations

Pass `conformers=` to compare states — apo/holo, tethered/extended, crystal vs
prediction, wild-type vs design. Any number of structures; **one button per
conformation**, no slider hunting.

```python
inspect_structure(
    "6ARU_A.cif",                                  # the reference
    layers=layers,
    conformers=[{"path": "1NQL_A.cif",  "label": "Tethered (1NQL)"},
                {"path": "AF-P00533.cif", "label": "AlphaFold"}],
    morph_mapping="intersection",                  # or "exact"
)
```

Only the endpoint traces are written to the file; the page interpolates the
in-between coordinates on load, so N conformations cost N frames rather than
`N + (N-1)×steps`. The animation is **Cartesian interpolation — a depiction of
the endpoints, not a pathway**; intermediates are not physical and bond geometry
is not preserved. Say so in an `about` tab whenever the morph is part of an
argument.

`morph_mapping="exact"` (default in the engine) refuses anything but an
identical residue set across all files. `"intersection"` opts in to the shared
residues — needed to compare a crystal chain against a full-length prediction —
and `report["morph"]["conformers"][i]["residues_dropped"]` says what each file
lost. Check it: a dropped epitope is a silent hole in the comparison.

## Delivering it

```python
save_artifacts(files=[report["path"]], language="python")
```

Then embed the returned version id inline. Keep bundles under ~20 MB
(`report["size_warning"]` appears above that) — a big structure with a morph is
the usual cause; trim the chain or drop a conformer.

## Verifying

`inspection_table(path)` returns one row per modelled residue —
`chain_id, author_residue_number, canonical_position, residue_name, x/y/z,
plddt, color, layers[], labels[]` — recovered from the rendered file itself, so
it proves what the reader will actually see. There is also a no-import CLI:
`python -m bindos_structure_inspector bundle.html residues`.

## Choosing layers

- **Six or fewer layers** on screen at once. Beyond that the legend is the
  figure and the structure is decoration.
- **One fact per layer.** "Epitope" and "glycan shadow" are two layers, not one
  two-colour layer.
- **Tiling layers hide the base.** Four domain layers covering the whole chain
  means every residue is painted; add them last and expect the reader to untick
  them to see anything underneath.
- Let colour carry annotation meaning only — leave the unannotated chain in the
  engine's grey and don't also set `style`/`color` display options.
