---
name: bindos-inspector
description: "Turn a structure and a set of residues into one self-contained interactive HTML page: named colour layers the reader can toggle, optional partner chains, and optional button-driven morphing between any number of conformations. Reach for it whenever an analysis has picked out residues that someone needs to see ON the structure rather than read as a list — epitopes, hotspots, pockets, conserved or divergent positions, mutations, contacts, confidence, anything per-residue. Renders inline in chat as an artifact and opens anywhere with no install."
---

# Structure inspection bundles

One structure in, one HTML file out. It opens with no server, no network and
no install, and in Claude Science it renders **inline** as an artifact, so the
reader never leaves the conversation.

The unit it works in is a **layer**: a named, coloured set of residues with an
optional note on each. What a layer *means* is entirely yours — the tool has
no opinion about biology. A layer has been a crystallographic epitope, a
docking hotspot set, columns of a conservation alignment, residues above a
pLDDT cutoff, positions a mutagenesis scan called dead, the output of a
clustering run, or just "the twelve residues I want you to look at".

Use it whenever the result of some work is *which residues*, and a list of
numbers would make the reader open a viewer themselves. That covers far more
than design work: showing a collaborator what your script found, checking your
own annotation landed where you meant, putting a rotatable figure in a report,
or handing someone a structure they can interrogate without your environment.

## The one thing to get right

**Show less than you have.** A bundle with every layer you could compute is a
legend with a structure behind it. Pick the few that carry the argument; the
reader can always be sent a second bundle. Six layers visible at once is
plenty, and two is often the whole point.

## Setup

```bash
pip install git+https://github.com/profdocpizza/bindos-structure-inspector@v1.13
```

Or set `BINDOS_INSPECTOR_HOME` to a local checkout (`~/code/bindos-structure-inspector`
is probed by default). `gemmi` is needed only for `.pdb` input.

## Making one

```python
report = inspect_structure(
    "target.cif",
    layers=[
        {"id": "core", "label": "Hydrophobic core", "color": "#8fa8bd",
         "residues": range(40, 78)},
        {"id": "picked", "label": "Residues from the scan", "color": "#dc2626",
         "residues": {56: "largest effect", 60: "second", 91: "third"}},
    ],
    about={"What this is": "...prose, HTML allowed, becomes a tab..."},
    out="inspection.html",
)
save_artifacts(files=[report["path"]], language="python")
```

`residues` takes a range, a list, or a `{number: note}` mapping when each
residue deserves its own caption. Numbers are **author residue numbers** in
the file, so renumber before rendering if your analysis speaks a different
numbering. Later layers paint over earlier ones — order background first.

Always read `report["unmapped_residues"]`: residues absent from the model are
skipped, and silently if you don't look.

## What else it can do

**Partner chains.** Anything in the same file that isn't the target —
`partner_chains=["B", "C"]` — gets one toggle button. Any chain: a bound
antibody, a ligand-bearing chain, the other half of a dimer, a docked design.

**Conformations.** `conformers=[{"path": ..., "label": ...}]`, any number, one
button each, and any state morphs straight to any other. Only the endpoints
are stored; the page interpolates. Each conformer may bring its own
`partner_chains`, and a partner is shown only while its own state is on
screen — a partner solved against one conformation, left draped over
another's coordinates, is a composite passed off as an observation.

Two warnings worth repeating to your reader in an `about` tab:

- The morph is **Cartesian interpolation between endpoints, not a pathway**.
  Intermediates are not physical and bond geometry is not preserved.
- `morph_mapping="intersection"` (as against `"exact"`) drops residues not
  shared by every file. Check
  `report["morph"]["conformers"][i]["residues_dropped"]` — a silently dropped
  region is a hole in the comparison.

**About tabs.** `about=` takes a string, a `{title: body}` mapping, or a list
of `{"title", "body"}`; each becomes a tab beside the structure. HTML is
allowed and images can be inlined as data URIs. Context belongs here rather
than in a second file — a bundle is one file precisely so that nothing arrives
detached from it.

## Checking it

`inspection_table(path)` reads the rendered file back and returns one row per
modelled residue — chain, author number, canonical position, residue name,
coordinates, pLDDT, colour, layers, labels. It reports what the reader will
actually see, so it is the honest check that an annotation landed on the
residue you named. There is also a CLI:
`python -m bindos_structure_inspector bundle.html residues`.

## Notes

Keep bundles under ~20 MB so they stay emailable; `report["size_warning"]`
appears above that. Save Image writes a capture bar under the viewer with a
right-click-saveable PNG, because the browser download is silently dropped in
some embedded frames.

Built on **py2Dmol** by Sergey Ovchinnikov
(<https://github.com/sokrypton/py2Dmol>), inlined into every bundle. Layers,
morph, partners and export are by profdocpizza. Both BEER-WARE (Revision 42) —
free to reuse, keep the notice. Every bundle carries both; don't strip the
credit footer from one you pass on.
