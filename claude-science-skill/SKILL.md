---
name: protein-inspector
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

A layer can also be built and listed but start switched **off**, with
`"visible": False` on the layer dict. Use it for the supporting set that would
crowd the opening view: the reader still has it, one click away, and the first
thing they see is the thing you brought them for.

## Setup

```bash
pip install protein-inspector
```

Or set `PROTEIN_INSPECTOR_HOME` to a local checkout (`~/code/protein-inspector`
is probed by default). `gemmi` is installed with it.

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
residue deserves its own caption. `"visible": False` keeps a layer in the
legend but unchecked at open. Numbers are **author residue numbers** in
the file, so renumber before rendering if your analysis speaks a different
numbering. Later layers paint over earlier ones — order background first.

Always read `report["unmapped_residues"]`: residues absent from the model are
skipped, and silently if you don't look.

## What else it can do

**Partner chains.** Anything in the same file that isn't the target —
`partner_chains=["B", "C"]` — gets one toggle button. Any chain: a bound
antibody, a ligand-bearing chain, the other half of a dimer, a docked design.

**Conformations — show them when you have them.** This is the feature that does
the most for a reader and the one most often left unused. One conformation
shows a shape; two show a mechanism — which domain swings, what closes over
what, where the hinge is. That is the big picture, and it is exactly what a
list of residue numbers cannot carry. Apo and bound, open and closed, wild type
and mutant, prediction and experiment, design before and after relaxation: when
the comparison is the result, put both states in.

`conformers=[{"path": ..., "label": ...}]`, any number, one named button each,
and any state morphs straight to any other — not a fixed tour. Each button
carries the **Cα RMSD to the reference**, so the size of the change is on
screen rather than guessed from the animation. Only the endpoints
are stored; the page interpolates, so states are cheap — three conformations of
a 148-residue protein is about 0.64 MB. Each conformer may bring its own
`partner_chains`, and a partner is shown only while its own state is on
screen — a partner solved against one conformation, left draped over
another's coordinates, is a composite passed off as an observation.

Two things for you to know, not to write a tab about:

- The morph is **Cartesian interpolation between endpoints, not a pathway**.
  Intermediates are not physical and bond geometry is not preserved. If that
  matters to the argument, say it in one clause of your context paragraph. It
  does not need a tab of its own, and a reader who has to open a tab to learn
  it will have formed the wrong impression before they get there.
- `morph_mapping="intersection"` (as against `"exact"`) drops residues not
  shared by every file. Check
  `report["morph"]["conformers"][i]["residues_dropped"]`. A silently dropped
  region is a hole in the comparison — and if something was dropped, that
  belongs in the same paragraph, not in a caveats tab.

**About tabs, and how few to write.** `about=` takes a string, a
`{title: body}` mapping, or a list of `{"title", "body"}`; each becomes a tab
beside the structure. HTML is allowed and images can be inlined as data URIs.
Context belongs here rather than in a second file — a bundle is one file
precisely so that nothing arrives detached from it.

**Then compress it.** Tabs are where a bundle goes wrong. The reader came to
look at a structure; every tab is somewhere else they have to go and something
they have to carry back. Most bundles want **one** tab: the paragraph a
colleague needs in order to read the picture — what this is, where the numbers
came from, what to look at. Write that, and fold everything else into it. A
second tab has to earn itself: a long table that would drown the paragraph, a
derivation, a methods block someone will actually check.

Four things that look like tabs and are not. Caveats about the whole figure
belong in the one paragraph. The legend is the Layers panel. The list of
residues is the Residues panel. How the controls work is not your reader's
problem — they can click.

Your tabs are yours and start empty; nothing is added to them.

## Putting a still image in a report

The bundle is for a reader, but a report wants a figure. Take one without
opening the page:

```python
shot = capture_bundle("inspection.html", "figure.png", dpi=300)
save_artifacts(files=[shot["path"]], language="python")
```

The background is transparent and the view is rendered at the requested dpi
rather than scaled up from the screen, so 300 dpi gives roughly 3800 px across.
It drives the bundle's own export, so the image is the one a reader would have
saved — and it captures the view as the page has it, so a layer marked
`visible: False` is off in the figure too. Orient the structure in the page
first if the default view is not the one that makes the point.

Needs a browser, which is an optional extra:
`pip install "protein-inspector[capture]"` then `playwright install chromium`,
or pass `chrome=` to use one already on the machine. Send the reader the HTML
as well as the figure — the still cannot be rotated, and the rotating is the
reason to use this tool at all.

## Checking it

`inspection_table(path)` reads the rendered file back and returns one row per
modelled residue — chain, author number, canonical position, residue name,
coordinates, pLDDT, colour, layers, labels. It reports what the reader will
actually see, so it is the honest check that an annotation landed on the
residue you named. There is also a CLI:
`python -m protein_inspector bundle.html residues`.

## Notes

The display area is fluid: it fills the column beside the Layers panel at
whatever width the page opens in, so `display={"height": ...}` is the size
worth setting and width is ignored in practice.

Keep bundles under ~20 MB so they stay emailable; `report["size_warning"]`
appears above that. Save exports a **transparent** PNG at the
chosen dpi and offers Copy image beside it, which is the route that still works
inside a frame that blocks downloads.

`extras=True` is **not** a way to get that figure. It writes a `.viewer.json`,
an `.svg` and a `.png`, but that PNG is a 900x900 single-colour C-alpha wire
trace with no alpha channel — a provenance thumbnail, not the rendered view. Use
it to prove what was loaded, never as the picture you hand someone.

Built on **py2Dmol** by Sergey Ovchinnikov
(<https://github.com/sokrypton/py2Dmol>), inlined into every bundle. Layers,
morph, partners and export are by Tadas Kluonis. Both BEER-WARE (Revision 42) —
free to reuse, keep the notice. Every bundle carries both; don't strip the
credit footer from one you pass on.
