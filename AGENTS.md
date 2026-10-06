# Protein Inspector — agent guide

Turn a structure and a set of residues into **one self-contained interactive
HTML page**: named colour layers the reader can toggle, optional partner
chains, and optional morphing between any number of conformations. It opens
with no server, no network and no install — from an email attachment, on a
machine that has none of your environment.

Reach for it whenever the result of some work is *which residues*, and a list
of numbers would make the reader open a viewer themselves: epitopes, hotspots,
pockets, conserved or divergent positions, mutations, contacts, confidence,
anything per-residue. It is not only for design work — showing a collaborator
what your script found, checking your own annotation landed where you meant,
or handing someone a structure they can interrogate all qualify.

```bash
pip install protein-inspector
```

## The call

```python
from protein_inspector import inspect_structure, inspection_table

report = inspect_structure(
    "target.cif",
    layers=[
        {"id": "iface", "label": "Interface, <4.5 A", "color": "#9aa7b4",
         "residues": [31, 34, 35, 38, 42, 45]},
        {"id": "hot", "label": "Hotspots from the scan", "color": "#dc2626",
         "residues": {56: "largest effect", 60: "second", 91: "third"}},
    ],
    partner_chains=["B"],
    about={"What this is": "Contacts at 4.5 A from the deposited coordinates."},
    out="inspection.html",
)
print(report["path"], report["unmapped_residues"])
```

`.cif` or `.mmcif` in; `.pdb` is converted for you. `residues` takes a range,
a list, or a `{number: note}` mapping when each residue deserves its own
caption. Numbers are **author residue numbers** in the file — renumber before
rendering if your analysis speaks a different numbering. Later layers paint
over earlier ones, so order background first.

## Show the morph

**If you have two or more conformations, show them.** This is the feature that
does the most for a reader and the one most often left on the table.

A reader looking at a single conformation sees a shape. A reader who can watch
it move sees the mechanism — which domain swings, what closes over what, where
the hinge is. That is the big picture, and it is precisely what a list of
residue numbers can never carry. If your work involves an apo and a bound
form, an open and a closed state, a wild type and a mutant, a predicted model
and its experimental counterpart, or a design before and after relaxation, the
comparison *is* the result. Put both in.

```python
report = inspect_structure(
    "closed.cif",
    layers=[...],
    conformers=[
        {"path": "open.cif",  "label": "open (4AKE)"},
        {"path": "bound.cif", "label": "bound to AP5A (1AKE)", "partner_chains": ["B"]},
    ],
    reference_label="closed (1ANK)",
    morph_mapping="intersection",
    out="inspection.html",
)
```

What the reader gets: one named button per state, any state morphing straight
to any other — not a fixed tour — and the **Cα RMSD to the reference printed
beside each button**, so the size of the change is on screen rather than
inferred from an animation. Only the endpoints are stored and the page
interpolates, so states are cheap: three conformations of a 148-residue
protein is about 0.64 MB.

Two things to know yourself rather than write a tab about:

- The morph is **Cartesian interpolation between endpoints, not a pathway**.
  Intermediates are not physical and bond geometry is not preserved. If that
  matters to the argument, say it in one clause of your context paragraph.
- `morph_mapping="intersection"` drops residues not shared by every file;
  `"exact"` requires them to match. Check
  `report["morph"]["conformers"][i]["residues_dropped"]` — a silently dropped
  region is a hole in the comparison, and if something went it belongs in that
  same paragraph.
- **The mapping is by address, not by sequence.** Residues are matched on
  (chain, author number) and their NAMES are never compared, so a point
  substitution between two entries passes `"exact"` in silence. Morphing a
  wild type against a mutant is legitimate, but you should know it is what you
  are doing: compare the residue names yourself before rendering. It now shows
  up in the picture, because each state draws its own side chains while the
  strip and the labels read the reference's identity.

## Partner chains

Anything in the same file that is not the target — `partner_chains=["B", "C"]`
— gets one toggle button: a bound antibody, a ligand-bearing chain, the other
half of a dimer, a docked design.

Each conformer may bring its own `partner_chains`, and a partner is shown only
while its own state is on screen. This is enforced, not merely documented: a
partner solved against one conformation and left draped over another's
coordinates is a composite passed off as an observation.

It may bring its own `partner_label` too, and should whenever the states hold
different things — the button is a caption, and one string over an apo trimer
and a receptor complex is wrong for one of them.

A partner that is a COPY of the target is painted with the target's layers:
`layers_on_partners="auto"` (the default) tells copies from binding partners
by sequence, so a homo-oligomer comes out annotated on every subunit while a
receptor or Fab in the same file keeps its own colour. The copies are still
partners — the toggle hides them with the rest. `True` paints every partner
chain, `False` none; `report["partners"]["painted"]` says which it did.

## Show less than you have

A bundle with every layer you could compute is a legend with a structure
behind it. Pick the few that carry the argument; you can always send a second
bundle. Six visible layers is plenty, and two is often the whole point.

**Write one About tab, not five.** `about=` takes a string, a `{title: body}`
mapping, or a list of `{"title", "body"}`; each becomes a tab beside the
structure, HTML is allowed, and images can be inlined as data URIs. But every
tab is somewhere the reader has to go and something they have to carry back.
Most bundles want one: the paragraph a colleague needs in order to read the
picture — what this is, where the numbers came from, what to look at. Fold
everything else into it. A second tab has to earn itself: a long table that
would drown the paragraph, a derivation, a methods block someone will check.

Four things that look like tabs and are not. Whole-figure caveats belong in
that paragraph. The legend is the Layers panel. The residue list is the
Residues panel. How the controls work is not your reader's problem.

Your tabs are yours and open empty; nothing is added to them.

## Check it before you hand it over

**Always read `report["unmapped_residues"]`.** Residues absent from the model
are skipped, and silently if you do not look.

`inspection_table(path)` reads the rendered file back and returns one row per
modelled residue — chain, author number, canonical position, residue name,
coordinates, pLDDT, colour, layers, labels. It reports what the reader will
actually *see*, which makes it the honest check that an annotation landed on
the residue you named, and it is CSV-ready as it stands.

Keep bundles under ~20 MB so they stay emailable; `report["size_warning"]`
appears above that.

## What the reader gets

Emitted by the renderer — do not re-implement or hand-edit it in generated
HTML; `tests/inspection_dom.js` clicks these controls rather than asserting
their strings appear.

- **Layers**: `All` / `None` / `Select`, then one row per layer with a colour
  swatch, the label, its annotation count, and two hover-revealed buttons:
  `only` isolates the layer, `sel` makes its residues the selection
  (shift-click adds, so several layers can be selected together). `Select` in
  the header takes every **ticked** layer at once. This is how the selection
  tools reach a whole epitope rather than one residue: `sel` then
  `Side chains`, or `sel` then `Select around`.
- **Residues**: a live count, a filter box matching residue label, layer label
  and annotation id, and the list grouped into one collapsible section per
  layer.
- **Residue card**: clicking a row — or a residue in the structure, or a
  letter in the sequence — opens a card with the residue identity, chain /
  author number, and one line per annotation on it, which is sometimes none.
  Dismissed by `×`, by clicking the row again, or by Escape; all three drop
  the 3D selection.
- **The structure is clickable**: one click selects a residue, shift-click
  extends, double-click takes its chain, a click on the background clears.
  py2Dmol ships with picking switched off in a notebook build; bundles turn it
  on.
- **Sequence**: a full-width strip under both panes, one letter per modelled
  residue in file order, wearing the layer colours the structure wears and
  printing the author number over every tenth residue. **The strip adds and
  subtracts; it never replaces.** A click or a drag that starts on an
  unselected letter adds that letter or range to what is already selected; one
  that starts on a selected letter removes it, and dragging back over your own
  path shrinks the change rather than ratcheting it. Shift extends from the
  last letter pressed, and the chain button takes the whole chain unless it is
  already all in, in which case it drops it. Hovering the structure lights the
  corresponding letter. The seam above it trades height between the picture
  and the strip — one seam, two rows, constant sum — so the page's height does
  not change and nothing below it moves. The vertical seam still owns width.
- **Selection tools**, on the sequence bar, acting on whatever is selected
  however it was selected:
  - `within N Å` + `Select around` adds every residue with an atom inside that
    distance — the renderer's own atom-to-atom search, so side chains count,
    and the selection it grew from is kept (PyMOL's
    `byres (all within N of sele)`).
  - `side chains` is a **Show / Hide pair, not a latch**. Selecting changes
    nothing about what is drawn; you press afterwards, and the press is
    absolute rather than a toggle. The pair reads back the state of the
    selection as well as setting it: Show fills when every selected residue is
    drawn, Hide when none is, and **neither** fills when they disagree or when
    nothing is selected — so one press always resolves a mixed selection in a
    known direction instead of inverting each residue separately. Deselecting
    does not undo anything: what is drawn stays drawn, and the readout says how
    many residues are still drawing side chains. Works on a morph, where each
    conformation carries its own rotamers and the animation shows those of the
    state it is leaving. Needs the atoms: both halves disabled, with the reason
    in their tooltip, on a bundle rendered `display={"sidechains": False}` or
    from coordinates that carry none.
  - `Clear` drops the selection.
- **Save**: Capture opens the panel; Save writes a **transparent** PNG at the
  chosen dpi, with Copy image beside it — the route that still works inside a
  frame that blocks downloads.
- `window.proteinInspector` exposes `syncVisibleLayers`, `clearDetails`,
  `setAllLayers`, `onlyLayer`, `selectAnnotation`, `selectResidues(positions)`,
  `selectAround(angstroms)`, `setSidechains(on)`, `sidechainState()`,
  `selectLayer(id, add)`,
  `selectVisibleLayers()` and `selection()` — the selection tools, for a
  figure that is scripted rather than clicked. `positions` are 0-based indices
  into what is drawn, not author residue numbers.

**A frame must declare its own `position_types`.** Not optional, and not
cosmetic: `_materialiseSidechains` appends one position per side-chain atom and
types each `'L'`, building the array from `(data.position_types || []).slice()`.
With no types in the frame that starts empty, the finished array is shorter than
the coordinate count, and `_setDataField` silently replaces a length-mismatched
per-position array with `Array(n).fill('P')` — so every side-chain atom claims
to be protein backbone, and `cartoon/geom.js` draws it as cartoon rather than as
a stick. That is one missing keyword argument between a correct picture and
ribbon slabs running between a residue's atoms. Derive the value with
`_position_types(names)`; never fill `'P'`, because `_ca_trace` accepts `C4'` as
well as `CA` and a nucleic position typed `'P'` loses its atoms without a word.

## Colour

The panel's **Base colour** control chooses how the un-annotated structure is
coloured: `Custom colour` (the default — a flat grey with a picker beside it)
or any mode the page reports from `window.py2dmol_colorModes()`: `auto`,
`chain`, `rainbow`, `plddt`, `ss`, `hydrophobicity`, `entropy`, `deepmind`,
`object`. Set the starting point with `base_mode` / `base_color`.

**A flat base is one choice among the modes, not a floor under them.** An
explicit per-position colour beats the mode, and the mode only decides the
positions nobody spoke for — so painting every position speaks for all of them
and silently disables rainbow/plddt/chain/ss. The base is therefore seeded only
when the mode is `custom`. Two consequences: under `custom`, unticking every
layer leaves a clean flat structure rather than the rainbow `auto` resolves to
on a single chain; and changing the viewer's own colour dropdown drops the
panel out of `custom`, because otherwise the base would override what the user
just picked.

The default style is **`tube`**. Pass `display={"style": ...}` for `cartoon` /
`richardson` / `ribbon` / `3d` when a figure needs it.

`display` also takes `width` / `height` (the starting split for the horizontal
seam), `color`, `chain_palette`, `background`, and two switches you will rarely
want to change:

- `gpu` (default `True`) — the WebGL2 painter, which `paintgl.js` describes as
  keeping the drawing "resident on the GPU so that turning the model costs one
  draw call instead of one full repaint"; the CPU painter re-runs
  `cartoon/geom.js` every frame instead. Both are in every bundle, so `False`
  is for reproducing a rendering difference, not for compatibility: a machine
  without WebGL2 falls back on its own, per frame. No head-to-head timing is
  quoted here — none has been measured, and upstream publishes none.
- `sidechains` (default `True`) — carry the side-chain atoms, which is what the
  `side chains` pair draws. ~8% of the bundle; `False` disables both halves.
- `detail` (default **`8`**, py2Dmol's own default is 4) — cartoon subdivisions
  per residue, 2–8. This is the only thing that sets sampling: upstream retired
  an adaptive term that targeted a fixed on-screen chord length, so
  `cartoon/geom.js` now says plainly that "magnified curves facet rather than
  resample". At 4 a residue gets 4 subdivisions in a helix and 3 in a strand or
  loop, which close in is three flat plates per residue of loop, each with its
  own outline stroke — read by the first person to zoom a bundle in as
  malformed side-chain bonds. A bundle exists to be zoomed into, so it asks for
  the top of the range; the `Detail` slider still spans 2–8 live.

**The Style panel is folded, not rewritten.** `src/parts/panel.js` builds it
from one table of rows shared with py2dmol.solab.org — "ONE PANEL, TWO SKINS" —
so the website exposes the same seventeen controls and there is no simpler set
to borrow from it. What differs is the room: the site is a wide sidebar, ours
floats over the canvas. So `Style`, `Detail`, `Color` and `Sele` stay out
front and the rest go into a closed `Fine tuning` group — along with one
control the shared panel does **not** ship: **`GPU painter`**. `parts/ui.js`
leaves `useGPU` out deliberately, and that is right for a notebook, where the
build carries one painter and the flag would ask for a file that is absent. A
bundle carries both, and the two do not draw a close-up side chain the same
way: the 2D painter's ink pass drops the edges of a stick box that lies inside
another, leaving one silhouette per side chain, where the WebGL2 path has been
reported showing those interior edges through the faces. Only built when
WebGL2 actually answers (`py2dmolCartoonGPU.available()`), per py2Dmol's own
rule that a control which cannot do anything is worse than none.

Move whole **rows**
if you change this: `parts/ui.js` hides per-style controls by setting `hidden`
on the row carrying `data-style`, so a control lifted out of its row survives
into a style that cannot honour it. A `.half` cell carries its own
`data-style` and may move alone, which is how `Detail` rides up beside
`Style`.

## One visual channel, not two

`highlight` picks how layers mark residues, bundle-wide:

- `"color"` (default) paints the residues in their layer colour. **Prefer this.**
- `"halo"` leaves the structure its own colour and rings the residues with the
  selection highlight instead.

Do not try to get both. Painting and haloing the same residues encodes one
fact twice, and the halo is the weaker signal — thin, poor against a cartoon,
and competing with the colour underneath. There is no per-layer form for the
same reason. Reach for `"halo"` only when the structure's own colouring is the
subject (a pLDDT- or chain-coloured view you need to keep intact) and annotate
on top of it.

One more cost under `"halo"`: the halo *is* the selection channel, so a
layer and the reader's own selection compete for it — ticking a layer
replaces whatever they had selected. Selecting, `Select around` and
`Side chains` all still work, and the sequence strip still marks the
selection, but under `"color"` the two channels never collide.

## Other hosts

Nothing in the package knows which agent is calling it.

| host | entry |
| --- | --- |
| Claude Science / Claude Code | the `protein-inspector` skill; `skills/protein-inspector/SKILL.md` is its source |
| any shell | `protein-inspector render spec.json` — same keys as `inspect_structure`, prints the bundle path |
| | `protein-inspector read bundle.html residues` — the table above, as JSON |
| MCP clients | `uvx --from "protein-inspector[mcp]" protein-inspector-mcp` → `render_bundle`, `read_bundle` |

## The strict surface underneath

`render_inspection_bundle` takes a validated manifest in which every annotation
carries a `canonical_position` that is **not** the author residue number:

```python
{"schema_version": "protein-inspector-manifest-1",
 "highlight": "color", "base_mode": "chain",
 "annotations": [
   {"annotation_id": "pocket-42", "kind": "candidate_hotspot",
    "label": "Y42 — rim of the cleft", "layer_id": "pocket",
    "layer_label": "Binding pocket", "color": "#dc2626", "resolved": True,
    "method": "structure-derived annotation", "evidence_ids": [],
    "residue": {"component_id": "target", "canonical_position": 18,
                "chain_id": "A", "author_residue_number": 42}}]}
```

Getting `canonical_position` wrong colours the wrong residue and raises
nothing, which is why `inspect_structure` derives it from the file instead of
asking. Use the manifest form only when you need a field the simple form does
not expose, and verify the mmCIF SHA-256 before rendering. An unresolved claim
should carry `"resolved": false` rather than a guessed address.

Internals of the vendored viewer — file map, paint order, GPU lifecycle — are
in `docs/viewer_internals.md`, not here.
