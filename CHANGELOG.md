# Changelog

## Unreleased

### Added

**The other protomers wear the annotation, and the button names what is on
screen.** Two things a homo-oligomeric target got wrong, both of them the
same mistake: treating a copy of the target as if it were a stranger.

- **`layers_on_partners`**, default `"auto"`. A trimer is one target
  presenting the same site three times, but the layers were painted only on
  the chain the analysis named, leaving the reader to work out that the grey
  subunits carry the same epitope. `"auto"` paints the copies and only the
  copies: a partner qualifies when it shares most of the target's author
  residue numbers and agrees on the residue name at essentially all of them
  (and there are at least ten — two copies of a tetrapeptide are not
  evidence). A receptor, a Fab or a ligand chain in the same file fails that
  test and keeps its own colour, which is the point: painting an epitope onto
  a receptor would be a claim nobody made. `True` paints every partner chain,
  for numbering you have made equivalent yourself; `False` is the old
  behaviour. The chains actually painted come back as
  `report["partners"]["painted"]`, and `read_inspection_bundle` mirrors the
  rows the same way, so the read-back still reports what the reader sees
  rather than what the manifest said. A residue missing from one copy is
  skipped rather than raised on — a crystal models a loop in one subunit and
  not the next.
- **A partner label per conformation.** `partner_label` was one string for
  the whole bundle, so a morph between an apo trimer and a receptor complex
  had a button reading "other protomers" over a receptor — a caption
  contradicting its own picture. A conformer may now carry its own
  `"partner_label"`; the page swaps the button's text with the state and
  falls back to the global label for any state that did not name one.

**A sequence strip, a clickable structure, and selections made from
selections** (`protein-inspector-1.14`). Three things a reader of a bundle
could not do, all of which are the same question — *this residue, and the ones
around it*:

- **The structure is clickable.** py2Dmol gates canvas picking on
  `selectionEnabled`, and the notebook bundle leaves it `false` (`parts/ui.js`
  turns it on for Focus mode alone), so an exported page's canvas was inert:
  single click, shift-click, double-click-for-chain and clear-on-background
  were all already written in `core/mol.js` behind one flag that nothing set.
- **A sequence strip**, full width under both panes, with its own draggable
  seam for its height so it cannot take width from the picture. One letter per
  modelled position in file order, in the layer colours the structure is
  wearing — the strip is the Layers legend with an index — the author number
  over every tenth residue, click/drag/shift/ctrl to select, and a chain
  button for the whole chain. Hovering the structure lights the letter, via
  the `window.SEQ.setHoveredResidue` hook `core/mol.js` already calls when a
  strip is present.

  Written here rather than borrowed. `src/panels/seq.js` is not in the
  notebook bundle and could not be dropped in: it is a canvas strip wired to
  `src/app/selection.js` (also web-only), and it colours by the viewer's
  colour mode, which is the one thing that would make it disagree with the
  page it is on.
- **`Select around`**, with a distance box, calling the renderer's own
  `residuesWithin` — atom to atom on a grid it keeps between calls, so side
  chains count — and keeping the selection it grew from, as PyMOL's
  `byres (all within N of sele)` does.
- **`side chains`**, a Show / Hide pair drawing the selected residues' side
  chains. `parts/sidechains.js` had the verb and nothing in an
  exported page reached it. The atoms have to be in the file, so
  `view(sidechains=True)` is now set; `display={"sidechains": False}` declines
  it, and the button is then disabled with the reason in its tooltip rather
  than raising. Morph and conformer bundles build their frames from Cα
  coordinates and carry none either way, which the page reads off the state
  rather than guessing.

  What the atoms cost, measured on 3PTB (233 drawn positions, 223 of them
  with a side chain): **607,877 → 656,162 bytes, +7.9%**. They are stored once
  per residue rather than once per frame, so a morph does not multiply them.

**A layer is a selection.** The tools act on a selection and a layer is
already a named, counted, coloured set of residues, so the two halves of the
page were not connected: side chains could only be drawn for residues picked
one at a time, and the unit anyone actually wants is "this whole epitope".
`sel` on a layer row selects it (shift-click adds, so several layers go
together) and `Select` in the Layers header takes every **ticked** layer at
once — the same visibility the structure is painted from.

**Side chains work on a morph.** The conformer path builds its frames from Cα
coordinates, so the new button was dead on exactly the bundles most worth
interrogating. It works because the browser's table is not world-space: each
atom is three coefficients in its residue's own backbone frame, rebuilt from
the final positions at draw time, and `parts/ui.js` builds a separate table
per frame. So one row set per conformer gives every state its own rotamers,
and the interpolation buffer borrows the table of the state it is leaving —
frozen for the length of the animation, which is what the interpolation itself
already is. The atoms have to travel through the same Kabsch fit as the trace
they belong to, or every side chain is measured against a backbone that has
rotated out from under it.

Measured on the EGFR ECD bundle (609 residues, 559 with a side chain, two
conformations): 2,309,668 → 2,627,741 bytes, +13.8%.

The residue card opens on any residue now, not only an annotated one, and says
so when there is nothing on it.

`window.proteinInspector` gains `selectResidues`, `selectAround`,
`setSidechains`, `sidechainState` and `selection`.

**Selecting no longer draws anything, and selecting never erases.** Two
borrowings from `src/app/selection.js`, both of which the website had and an
exported page did not:

- The side-chain control is a **Show / Hide pair with three states**, and the
  press is absolute rather than a toggle. It reads the selection's state back
  per residue from `shownSidechainSet()`: Show fills when all of them are
  drawn, Hide when none is, and neither fills when they disagree or when
  nothing is selected. A single latch could only show that third state as a
  grey smear — it said what it would do and left you to work out what it had
  done — and on a mixed selection it inverted each residue separately, which
  is never the thing anyone wanted. Deselecting now undoes nothing: what is
  drawn stays drawn, and the readout counts it.
- **The strip adds and subtracts.** A click or drag beginning on an unselected
  letter adds; one beginning on a selected letter removes; the mode is decided
  at pointer-down and the range is recomputed from a snapshot of the selection
  taken there, so dragging back over your own path shrinks the change instead
  of ratcheting. A click never replaces what was selected.

**The cartoon was sampled for a reader who does not zoom in.** The first person
to zoom a bundle to where one residue fills the canvas reported the side-chain
bonds as malformed. They were not: the side-chain table checks out exactly on
all 559 residues of the EGFR bundle — every heavy-atom and bond count,
including every aromatic ring closure — and the stick section is the website's
own. What they were looking at was the cartoon.

`detail` is the only thing that sets cartoon sampling, and py2Dmol defaults it
to 4. Upstream retired an adaptive term that targeted a fixed on-screen chord
length — it was the largest single term in the frame, and it made a picture
depend on the zoom it was first built at — so `cartoon/geom.js` now states the
consequence plainly: *"Magnified curves facet rather than resample."* At 4, a
residue gets 4 subdivisions in a helix and 3 in a strand or loop, which close
in is three flat plates per residue of loop, each carrying its own outline
stroke. That is what read as extra bonds and rough box joins.

A bundle is the one viewer made to be zoomed into, so `display` takes `detail`
and defaults it to **8**, the top of the range and the sampling the retired
adaptive term used to reach (~8.7 per residue on a 600 px canvas). Affordable
because the GPU painter is on: the doubled station count is paid once, at
capture, instead of on every frame of every rotation.
`display={"detail": 4}` restores
py2Dmol's default, and the `Detail` slider still spans 2–8 live.

**The Style panel is folded.** `Style`, `Detail`, `Color` and `Sele` stay out
front; `Width`, `Outline`, `Thick`, `Flat`, `Hilite`, `Shade`, `Pencil`, `Ink`,
`Shadow`, `Ortho` and the two toggles go into a closed `Fine tuning` group.

**Which painter draws is a control now**, in that group, and it is the one
thing in the panel the shared table does not provide. `parts/ui.js` omits
`useGPU` on the grounds that it is a backend rather than a look — right for a
notebook, which carries one painter, and wrong for a bundle, which carries
both. The two do not draw a close-up side chain the same way: `cartoon/geom.js`
describes the 2D painter's rule for a stick box lying inside another — "one
that lies inside the other is behind its surface, so the ink pass removes it" —
and the WebGL2 path was reported showing those interior edges through the
faces, on a structure where the same residue on py2dmol.solab.org carried a
single silhouette. Whoever is looking at the picture is the only one who can
say which is right for it, and they cannot say it if switching means
re-rendering the file. Built only when `py2dmolCartoonGPU.available()` answers,
and it drops the GPU mesh on the way across rather than re-drawing geometry
built while the painter was off.

Nothing is removed, and nothing was borrowed from py2dmol.solab.org either —
there was nothing to borrow. `src/parts/panel.js` builds this panel from one
shared table of rows, and its own header calls it "ONE PANEL, TWO SKINS": the
website and the bundle get the same seventeen controls and differ only in
stylesheet. What differs is the room they have. The site is a wide left
sidebar; ours floats over the canvas, where seventeen controls read as clutter.

Whole rows move, because `parts/ui.js` hides per-style controls by setting
`hidden` on the row that carries `data-style` — a control lifted out of its row
would be a slider surviving into a style that cannot honour it. `Detail` is the
exception and is safe: a `.half` cell carries its own `data-style`, and it
moves into the row that is always visible.

**The cartoon is painted on the GPU.** `render_inspection_bundle` passed
`gpu=False`, while py2Dmol's own default — and the website's — is `True`, so a
`richardson` or `cartoon` bundle was handing the same geometry to the CPU
painter and felt slow next to py2dmol.solab.org on the same machine. The
difference is structural rather than a benchmark, and `paintgl.js`'s own header
states it: the drawing is "resident on the GPU so that turning the model costs
one draw call instead of one full repaint", where the 2D painter runs
`cartoon/geom.js` every frame. No head-to-head figures are quoted because
none were measured here and upstream publishes none.

Nothing is given up, because the notebook bundle is the one build that carries
BOTH painters and `core/mol.js` only honours the flag when both are present:
no WebGL2 or a lost context and `_gpuWillTake` declines the frame; paintgl
refuses any context carrying `getSerializedSvg`, so vector export still comes
off the 2D painter; and `protein-inspector capture` already launches headless
Chromium with SwiftShader. `display={"gpu": False}` forces the CPU painter.

### Fixed

**A conformer frame did not say what its positions were, and its side chains
were drawn as cartoon.** Reported as "the boxes connecting the atoms look wrong
in between" on a morph bundle, while the same structure rendered as a single
file through `add_pdb` was clean — with a byte-identical side-chain table.

`_materialiseSidechains` appends one position per side-chain atom and types
each `'L'`, building the array as `(data.position_types || []).slice()` plus one
push per atom. The conformer path called `viewer.add()` without
`position_types`, so that started **empty**: the finished array was as long as
the appended atoms rather than as long as the coordinates, and `_setDataField`
takes a per-position array whose length does not match the coordinate count and
silently replaces it with `Array(n).fill('P')`. Every side-chain atom then
claimed to be a protein backbone position — and `cartoon/geom.js` partitions on
exactly that, backbone drawn as cartoon and everything else as sticks. The
atoms were swept into the ribbon and drawn as slabs running between them.

Nothing was wrong with the atoms, which is why this survived every check made
while looking for it. The table is chemically exact — all 559 residues of the
EGFR bundle match expected heavy-atom and bond counts, every aromatic ring
closure included — and its geometry survives the Kabsch fit intact: CA–CB
median 1.530 Å and intra-side-chain bonds 1.488 Å in **both** frames, no
outliers. They were mislabelled, not misplaced.

Both conformer `add()` call sites now send types, derived per residue by
`_position_types` rather than filled with `'P'`: `_ca_trace` accepts `C4'` as
well as `CA`, so a stack can hold nucleic positions, and a base typed `'P'` is
rebuilt through the peptide's step range — `localFrame` fails and its atoms are
dropped in silence, which is how upstream once lost 347 of them.

Three explanations were proposed and withdrawn before this one: the GPU
painter, the cartoon sampling (a real defect, fixed above, but not this one),
and the load path. What localised it was rendering the same 6ARU through both
paths and diffing the frames.

**Detail did not survive the style switch it was rendered for.**
`display={"detail": N}` seeds the config and the constructor honours it, but
`_applyLookDefaults` re-asserts `cartoonDetail = d.detail` unconditionally —
no "did a person choose this" latch of the kind width and thickness carry — and
all three `LOOK_DEFAULTS` entries hold 4. A bundle opens as `tube`, so the first
pick of Richardson dropped the sampling from 8 to 4, and Richardson is the style
it was picked for: the file shipped a detail nobody could ever see. The page now
wraps `setStyle` / `setPreset` and re-asserts it afterwards — wrapped rather
than listening on the dropdown, because `setStyle` is also called
programmatically. Dragging the Detail slider makes the reader's value the kept
one.

**The sequence seam ran away from the pointer, and only the page got taller.**
Two faults with one cause. The drag was computed as
`grid.getBoundingClientRect().bottom - e.clientY`, and that edge moves when the
strip resizes — so every `pointermove` read back its own previous result and
the seam reached its ceiling in a few events. The vertical seam survives the
same expression because it reads `right`, and changing the panel's width does
not move the page's edge. The gesture is now anchored to `pointerdown` and the
two rows share a fixed sum, so the picture gives up exactly what the strip
takes, the page's height does not change, and the credit line underneath it
stops being the only thing that visibly moves.

**The selection mark appeared one frame late.** `core/mol.js`'s mouseup handler
sets the selection and then renders only `isLargeMolecule`, on the reasoning
that a small structure's next frame is along shortly — true on the website,
where a hover readout and a sequence canvas are repainting anyway, and false in
a still bundle. The yellow outline therefore appeared on the first frame AFTER
the click: nudge the structure and there it was. The page renders on every
canvas-driven selection change, which covers the background click too.

### Added (earlier in this cycle)

**`chain_palette=` on `view()`, and a colour per chain in the Style panel.**
Chain mode drew from one hard-coded palette — `chainColors`, a module-level
`const` in `core/mol.js` — so the only way to change what colour a chain came
out was to edit the source. Two ways in now:

- `view(color='chain', chain_palette='greys')` picks a named palette.
  `'pymol'` (the previous colours) is still the default, so nothing that does
  not ask changes. `window.py2dmol_chainPalettes()` lists them.
- Style → Chain shows one swatch per chain ID. Picking a colour pins that
  chain; double-clicking the swatch puts it back on the palette. Overrides are
  keyed by the renderer's chain key rather than the bare letter, because in a
  merged view two objects both have a chain `A`.

Colourblind mode still outranks both: it is an accessibility setting, not a
preference.

Everything that draws a chain — the viewport and the sequence strip — now goes
through one resolver, `chainColorHexFor`, so an override cannot be honoured in
one place and missed in the other. That was four separate palette lookups
before.

### Fixed

**`normalizeConfig` dropped unrecognised `color` keys silently.** It rebuilds
`config.color` from a named whitelist, so a key plumbed end to end through
`view()` and read back in `parts/ui.js` never arrived, with no error on the
way — the renderer simply used the default. Worth knowing when adding the
next colour option: the emitted JSON carrying a key proves nothing about
whether the code that consumes it kept it.

### Tests

`tests/test_chain_colors.py` executes the real `src/core/mol.js` under QuickJS
(new, test-only, `importorskip`-guarded) rather than asserting that strings
appear in the bundle. The system Node is 12.x and cannot parse the optional
chaining the source uses, which is how the whitelist bug survived review.

## 2.0.0

280 commits since `v1.6.5`. The major number is not for the size of it — it is
for the five things below that change what existing code does.

### Breaking

**`position_atoms` is gone.** `add()`, `replace()` and the payload no longer
take or carry a ligand atom's NAME. It was produced by both parsers, copied
through every field-by-field rebuilder and stored on the renderer, and read by
nothing — 2.7 KB a frame on 4HHB, 574 of whose 748 entries were blank, paid
again per frame per viewer in a notebook. `position_elements` stays and now
does two jobs: colour by element, and the per-pair bond thresholds. Passing
`position_atoms=` is now a `TypeError` rather than a silent no-op.

**Biological assemblies are built by default.** `add_pdb`, `from_pdb` and
`from_afdb` all took `use_biounit=False`; they take `use_biounit=True` now, and
the website and the embed do the same. A multimeric entry loads as the multimer
it is. Two consequences worth knowing before you upgrade:

- **The chains are renamed.** `gemmi.make_assembly` renames every chain it
  copies, so chain `A` comes back as `A1 A2 A3`. Any `set_color(..., chain='A')`
  or `add_contacts([['A', 10, 'A', 20, 1.0]])` against such an entry stops
  resolving — silently, because a selector that matches nothing is not an error.
  Pass `use_biounit=False` to keep the deposited chains.
- It is skipped when it would not expand the structure, compared by atom count,
  so a monomer is unaffected.

The reason it was off is that it never worked from an mmCIF:
`extractCIFBiounitOperations` asked the parser for loops without telling it what
to seek, and read the assembly only as a `loop_` when a file with one assembly
writes it as key-value items. The website's Load Biounit box has been ticked all
along and drawing the asymmetric unit.

**`set_color(chain=..., position=...)` together now raises.** It used to colour
the whole chain *and* those positions — a union — while the selector everywhere
else in py2Dmol reads the same pair as an intersection. Rather than let one word
mean two things, the combination is refused by name and the error says which
spelling you want. Either key alone is unchanged.

**Metal coordination is no longer drawn as a bond.** `_struct_conn` carries
`covale`, `disulf`, `hydrog` and `metalc`, and a coordination record is not a
bond: 7P1E declares Ca 506 chelated by both carboxylate oxygens of the ligand
K99, and those two sticks plus the carboxylate's own close a four-ring that
reads as a solid triangle. `metalc` is excluded.

**Lone-atom radii follow PyMOL.** `loneAtomRadiusA` is PyMOL's `ElementTable` —
Bondi where Bondi reaches, 1.80 for everything else. Eight elements change size,
calcium most visibly (2.31 → 1.80).

**`detect_cyclic` is now `cyclic`.** It was the only argument on `view()` named
for an action; everything else there names a state (`box`, `rotate`, `overlay`,
`multi`, `gpu`, `shadow`, `arrows`), and the toggle it drives is labelled
Cyclic. The verb-named family — `use_biounit`, `load_ligands`,
`filter_additives`, `ignore_ligands`, `allow_reflection` — all sit on the
loaders and describe a one-off instruction; this is a live setting. The config
key is `rendering.cyclic` and the checkbox is `#cyclicCheckbox`, with no alias:
a state file written before 2.0.0 loses that one setting on load and the toggle
comes back at its default.

**Some `view()` defaults moved**: `ortho` 1.0 → 0.5, `outline` `"full"` → `None`,
`width` 3.0 → `None` (the style decides). And `best_view`, `kabsch` and
`align_a_to_b` are gone from `py2Dmol.viewer` — the browser chooses the angle
now, so Python sends the *request* (`align`, `allow_reflection`) rather than the
result.

### New

- **The colour picker is organised by colour.** One row per family — reds,
  greens, blues, yellows, magentas, cyans, oranges, tints, grays — each running
  its own shades, which is how PyMOL's colour menu is laid out and where these
  values come from. It used to be the chain cycle: 40 colours in the order
  PyMOL hands them to chains, deliberately unlike itself from one to the next
  so neighbours contrast, which is right for chains and unreadable as a grid.
- **`view.show_sidechains()` / `view.hide_sidechains()`** — name the residues
  whose side chains are drawn, with the same selector as `clip` and `focus`.
  Relative, both of them: `show_sidechains(chain="A")` then
  `hide_sidechains(position=45)` is chain A without residue 45, and with
  nothing named either means every residue. Needs `view(sidechains=True)`,
  which is what carries the atoms; without it there is nothing to draw and the
  call says so. The verb itself moved out of `parts/embed.js` into
  `parts/sidechains.js`, so all three shells have it — the embed's JS API could
  draw a side chain and a notebook could not ask at all.
- **`view(style=...)`** — one flat list: `tube`, `richardson`, `ribbon`, `3d`,
  with `preset`, `smooth`, `thickness`, `sheet_flat`, `pencil`, `arrows`,
  `base_plates`, `detail`, `fade`, `highlight`, `outline_tint`, `shade`, `bg`
  and `ss_palette` beside it.
- **`view(gpu=...)`** — WebGL2 by default; `gpu=False` inlines the 2D painter
  instead, which is 46 KB smaller and the only build that can save an SVG. A
  notebook cannot fall back at runtime, so this chooses which bundle is written
  into the cell.
- **One notebook bundle, and `gpu` is a runtime setting again.** There were
  three — WebGL2, 2D, and a cartoon-less tube — because the library is inlined
  into the `.ipynb` once per `show()` cell. Sharing pays it once per document,
  so the notebook now ships both painters for 26 KB more: `gpu=False` reaches
  the 2D painter without a different file, a machine with no WebGL2 has a
  fallback, and **the cartoon can export an SVG from a notebook**. The embeds
  still ship one painter each — they are gzipped over HTTP, a different trade.
- **A PAE travels as base64, and the payload is written without spaces.** A
  PAE is N² numbers and it is inlined into the `.ipynb`: one 837×837 matrix was
  **72% of the demo notebook**. It is stored as a `Uint8Array` at 1/8 Å either
  way — that part is unchanged — but writing it as a JSON list of the scaled
  integers costs three characters and a comma each, where base64 of the same
  bytes costs 1.33. With compact JSON separators beside it (a megabyte of
  numbers was paying one space per element), the demo notebook goes from
  **4.32 MB to 2.03 MB**. Lossless, and `setData` still takes the three forms
  it always did, so an older payload still draws.

    | AF-Q5VSL9, 837×837 | |
    |---|---|
    | list of ints, `", "` separators | 3,048 KB |
    | …compact separators | 2,364 KB |
    | base64 of the bytes | 912 KB |
    | …resampled to the panel's 300px | **120 KB** |

- **A PAE carries no more resolution than the panel can draw.** The plot is an
  n×n image scaled into a canvas of `pae_size` pixels — 300 by default — so
  above that the browser was already throwing the detail away on every frame,
  and an 837-row matrix gave each residue 0.36 of a pixel. Doing the resample
  once, in Python, is the same picture. With the two changes above the demo
  notebook goes **4.32 MB → 1.22 MB**.

    The matrix side and the residue count are two numbers now, and `pae_n`
    carries the second: a box dragged on the plot is a range of *residues*
    handed to `setVisibility`, so it is scaled back out on the way. Selection
    edges on a resampled matrix land on a block of residues rather than one —
    on a plot where a residue was already a third of a pixel.

- **The notebook library is shared between cells where it can be.** Each
  `show()` used to write ~450 KB into its own output, because Colab gives every
  cell output its own iframe — ten viewers was 4.5 MB of `.ipynb`. The first
  viewer of a session now writes it and offers it on a `BroadcastChannel`;
  later ones ask, and keep a copy so they can lend to the next. Two viewers go
  from ~950 KB to ~505; ten from 4.5 MB to ~700 KB.

  There is no flag, and it is always on. In Jupyter every output shares one
  document, so a borrower finds the library already there and the channel is
  never used; in Colab the outputs are separate frames and the channel carries
  it. Eight viewers: **3,872 KB down to 645**.

  Python also pings the page at each `show()`, and a positive answer lets a
  fresh kernel borrow from a page that still has a lender. A negative one is
  ignored — it cannot be told apart from a question that never arrived.

  Re-running the cell that carries the library makes it carry the library
  again. Re-running a cell replaces its output, so the copy went with it while
  the kernel still believed it had lent one — the next viewer asked a lender
  that no longer existed and came up as an error box. A cell that creates
  several viewers is one execution, so a grid still writes the library once.

  The library is offered under a key that includes a **hash of its content**, so
  a cell can only borrow a library that matches the payload it is writing. A
  notebook is re-run cell by cell, and the cell holding the library can be from
  an older build than the cell now asking; today's payload handed to
  yesterday's renderer draws, silently missing whatever the two disagree about.
  A borrower that finds no matching lender inlines its own copy.
- **`view.focus(name=, chain=, position=, cutoff=)`**, a **Focus** button
  beside Style and `v.focus(sel)` on the embed — click a residue or a
  ligand and see what it is doing. Four things at once, each replaced by the
  next call: the residue is selected, the side chains of everything within 5 Å
  are drawn and the last focus's are taken away, the camera moves in, and a
  slab is cut around it. **It does not turn the structure** — only the centre
  and the zoom move, so focusing from one residue to the next walks through a
  structure rather than spinning it, which is the difference between this and
  `orient()`. Every step is an existing verb (`residuesWithin`, the object's
  side-chain set, `autoClip`); `parts/focus.js` is the composition, so all
  three shells get it from one place.

  A click in the **sequence strip** focuses too — it wrote the selection field
  directly and so went past everything hanging off a selection; it uses the
  renderer's setter now, as does Select all.

  The camera **moves over about a third of a second** rather than jumping, and
  clicking the background comes back out — an empty selection is the same
  signal the mode already reads, so the way back needs no second control.

  Side chains in a notebook need `view(sidechains=True)` — see below.
- **`view(sidechains=True)`** — a notebook can draw side chains. The payload
  carried one position per residue and nothing else, so `showSidechains`, the
  side-chain half of Focus, and any side-chain-to-side-chain measurement were
  web-only. Python now sends the raw atoms and the browser builds the table
  with the same code the website uses — `buildSidechainTable`, cut out of
  `src/io/parse.js` into `src/io/sidechains.js` so the notebook can have the
  chemistry (8.7 KB) without the parser (36 KB).

  **Off by default, and the reason is the payload.** Side-chain atoms are
  coordinates, so every frame carries its own: a 251-residue design goes from
  9.0 KB a frame to 37.2 — six frames is 54 KB against 223, and a hundred would
  be 2.8 MB inlined into the `.ipynb`. It is a viewer-wide setting rather than
  a per-load one because the table is per frame with no inheritance; mixing
  would make side chains appear and vanish as you step through frames.
- **`view.clip(name=, chain=, position=)`**, and a Clip button in the notebook
  and embed shells. `parts/clip.js` was in every bundle already — the slab, the
  tracking and the per-frame refit — and only the website could reach any of it.
  The depth is the selection's own depth along the view, so to cut deeper, clip
  to less; `clip()` with nothing turns it off.
- **`view(multi=True)` and `view.show_objects(...)`** — several objects in one
  picture, which is a different question from `overlay=True` (every FRAME of
  one object). The renderer has drawn several at once since the website grew
  its Multi button and the embed has exposed it as `v.showObjects()`; Python
  had no way to ask. `show_objects()` with no argument means every object
  loaded, resolved in Python at the moment of the call.
- **`view.orient(name=, chain=, position=, animate=)`** — turn the camera onto
  a selection. A viewer already does this once, unprompted, when the first
  frame lands; this is for afterwards.
- **`set_sse(sse, name=, chain=, position=)`** — force a region to helix, strand
  or coil, or pass `None` to return it to the automatic assignment.
- **Structural alignment** — TM-align, vendored from foldjs, running in a worker.
- **Cross-object contacts** — a contact whose two ends are in different objects.
- **An embeddable build**: `py2Dmol.embed.min.js` (453 KB, WebGL2) and
  `py2Dmol.embed.cpu.min.js` (414 KB, 2D and SVG-capable), documented by
  `embed.html`, with one selector grammar shared with the Python API.

### Fixed

- **Focus gives back a camera that shows what is drawn.** The mode snapshots
  the camera when you press Focus and restores it on every background click —
  a camera measured on the picture that was there at the time. Turn Multi on,
  or switch which objects are drawn, and that snapshot framed one object inside
  a scene of several: clicking the background put you back there every time,
  with Orient the only way out. It now widens to fit when the drawn set has
  changed, and the slab goes rather than coming back, since a near and a far
  along the view cut somewhere arbitrary once the picture is different.
- **The sequence strip forgets its hover when it is rebuilt.** Copy a
  selection, move the pointer up to a chain label, and part of the old
  selection stayed marked. The hover is a set of position indices and the only
  thing that took it back was a mousemove on the strip's canvas — which a
  rebuild destroys, so no `mouseleave` could ever arrive.
- **Flat ribbons get flat side chains on the GPU.** In `ribbon` — plain
  cartoon, whose whole look is flatness — side chains were drawn as fat 3-D
  sticks while the Thickness knob read 0, and dragging the knob fixed it until
  something rebuilt. The GPU floors the thickness at 0.05 Å before capturing so
  each ribbon piece is a closed solid and keeps its outline (at 0 the silhouette
  rule has one face to work with and the drawing loses almost every line:
  1,049 dark pixels against the 2D painter's 14,789). The stick rule reads that
  same field to ask whether the preset *wants* flat — `thickness === 0` — so the
  floored value said no, and every side chain took the ligand's own 0.5. The
  reader's value now survives the floor. The 2D painter never had the fault.
- **The PAE plot's selection box is drawn where the selection is.** A stored
  box is in residues and the canvas is laid out in cells, and only the mask had
  been converted — the black outline multiplied residue indices by a cell
  width, so on a matrix resampled for the panel (which happens in a notebook,
  above `pae.size`) it framed a region 1.2× out from the region it was
  highlighting. The chain boundary lines had the same fault from the other
  side, walking cells while indexing residues. `render()` converts once and
  hands cells to both.
- **Focus is a mode with one focus at a time.** Four things it got wrong:
  loading a structure, or anything else that emptied the selection, handed the
  Sele dropdown back to Highlight while the Focus button stayed lit; a click in
  the sequence strip ADDED to the focus rather than moving it, so each click
  focused the union and walked the camera to the centroid of everything ever
  clicked; switching objects and coming back parked you at the pocket you had
  focused with nothing marked, since the camera is per object and the selection
  was not; and merging objects left one structure's slab cutting through
  another. Clear All drops the mode too, rather than carrying a snapshot of
  objects that no longer exist into the next structure.
- **Anything that can be set can now be unset.** `set_color(None)`,
  `add_contacts([])` and `add_bonds([])` clear, and for `set_color` the clear
  reaches exactly as far as the selector that set it — the object, one chain,
  some positions, or one frame. Previously `set_color(None)` returned silently
  and `add_contacts([])` warned and refused, so a colour or a contact could be
  put on and never taken off.
- **The live path works in Colab.** Colab renders every cell output in its own
  iframe, so a later `add()` cannot reach the viewer directly and
  `BroadcastChannel` is the only bridge. It does not retain, and there was no
  handshake — `viewerReady` was posted by the viewer and listened for by nobody
  — so on a notebook **reopen** the update cells routinely posted before the
  viewer's channel existed and their frames were lost. They answer the
  announcement now, and because the replay arrives in whatever order the iframes
  ran, the viewer holds it briefly and applies it sorted.
- **`align=True` actually superposes a trajectory loaded with `add()` then
  `show()`.** The fitting moved to the browser, so the payload carries the
  request rather than the result — and `_display_viewer` builds its frames
  field by field and never named it, as did the static loader on the other
  side. Both ends were dropping it, each silently. `show()` then `add()` was
  unaffected, because the live path sends the frame whole. A helix and the same
  helix turned 90 degrees came out 7.07 A apart instead of 0.
- **`set_sse()` and `set_color(frame=N)` reached a live viewer at all.** The
  first was dropped on arrival by a second, drifted copy of the metadata
  applier; the second had no route, because a frame is delivered once.
- **The wheel ships what `viewer.py` opens.** `package_data` omitted the GPU
  renderer, which `viewer.py` reads unconditionally — a `FileNotFoundError` on
  the first `show()`, in the wheel only.
- 🔴 **Orient uses the whole canvas.** The best-view search already reads the
  window's aspect and lays the long axis along the long side; the framing then
  fitted the result into a square of side `min(width, height)`, so a 600×300
  viewer drew an elongated structure across 47% of its width. It now measures
  the shape under the rotation it chose. The same rod fills **94%**; a globular
  structure in the same window goes from ~50% to 99% of its height.

  The framing belongs to that rotation — turning a long structure end-on
  afterwards can push it past the edge, as in PyMOL. Press Orient again to
  reframe.
- 🔴 **An SVG of the tube kept its shading with the GPU on.** The CPU occlusion
  pass is skipped when the GPU is going to draw — it computes its own — but the
  question was asked of the renderer's state rather than of the context, and an
  SVG context is one the GPU refuses. So the export took the 2D path with a
  pass that had been skipped on its behalf, and `gpu` + `tube` + SVG came out
  flat.
- 🔴 **`bg` did nothing whenever `box` was off.** The box is the frame and `bg`
  is the paper: turning the frame off set the canvas background to transparent
  *after* the requested colour had been put on it, and told the renderer to
  clear transparent too. `py2Dmol.grid` defaults `box` to False, so
  `g.view(bg="black")` came out on white. White still means "not asked for" and
  still floats on the page, which is what `box=False` is for; any other colour
  is honoured. **`grid(bg=...)`** is a grid-wide default now, beside `size`,
  `controls` and `box`.
- **No white space under a viewer in Colab.** The stylesheet gave the canvas
  box 600×600 and a script corrected it, so the markup alone was 648px for a
  300px viewer. Colab inserts output HTML with `innerHTML` — which never runs a
  script — and sizes the output iframe from what it measures in that window: a
  2×2 grid of 300px viewers is ~1,220px of unsized markup, and the frame kept
  ~1,000px around a 644px page. The size is in the markup now.
- **A grid emits one output, not twenty-eight.** `Grid.view()` said "do not
  show yourself" by setting `_is_live` — and that flag also means "you are on
  the page, so send updates". The grid has not been emitted yet, so every
  `add()` during collection wrote an update addressed to a viewer that did not
  exist: four viewers came to twenty-seven of them, twenty from a single NMR
  ensemble, each an empty output element and together a band of white space
  under the cell. `_managed` says the first thing alone, and `Grid.show()`
  marks the viewers published afterwards — so an `add()` *after* the grid still
  reaches the viewer beside it, which is the one thing the old flag got right.
- **`rotate=True` turns the structure.** It reached the config, the
  constructor and the Rotate checkbox, and was then switched off by the
  viewer's own opening orient — which stopped the spin unconditionally, on the
  reasoning that a reader pressing Orient wants the view held. Nobody presses
  the automatic one. It passes `keepSpin` now; a deliberate Orient still stops
  the turn. Affected the notebook and the embed alike.
- **`py2dmol_scatter_loaded`** was dispatched on `document` and listened for on
  `window`; a bare `Event` does not bubble.
- An NMR ensemble kept one model of six when an assembly was built.
- A contact drifted when the view rotated: the GPU repack dropped `zBias` and
  `wA` from a line primitive, so the near-surface bias was applied in model
  space along whichever way the view happened to be pointing.
- 🔴 **A ribose-bearing cofactor is a ligand, not a nucleotide.** `add_pdb`
  promoted any residue carrying `C4'` + `O4'` + `C1'` to a nucleotide — a rule
  written for modified bases inside an RNA chain, which also caught SAM, SAH,
  ATP, NAD and FAD. Such a ligand collapsed onto a single position at its `C4'`
  and drew as one sphere, in every frame of a trajectory. The rule now applies
  only inside a chain that holds nucleic acid. The website was never affected:
  its parser asks a stricter question.
- Element colouring is on by default, and the parser infers an element from the
  atom name when columns 77–78 are blank.

### The interface

- **A Clip control in the notebook and the embed**, and `view.clip()` in Python
  — `parts/clip.js` was in every bundle and only the website could reach it.
  It shows whether it is on, which it previously did not.
- **The Cyclic toggle works outside the website.** The shared Style panel put
  it on screen in the notebook and the embed, where nothing wired it: it came
  up unticked while `cyclic` was true and the ring was closed, and
  clicking it did nothing.
- **An Orient control in the notebook.** The website and the embed both had one.
- **Draw is not offered where no 2D painter can honour it**, and SVG export of
  the cartoon is refused rather than writing an empty file. The tube exports a
  vector on every build, because it is stroked by the core rather than by a
  painter.
- **Capture lights up while its panel is open**, like Style beside it and Clip
  above it. The open cue in the notebook and embed shells was written as
  `#styleToggle[aria-expanded="true"]` — one button by name — so Capture put
  its panel up with its own button unlit and nothing said which of the two
  panels you had. Both shells key it on the state now, as `index.html` always
  did, and a latch (`aria-pressed`) and an open panel (`aria-expanded`) wear
  one skin.
- **A more compact viewer menu** in the notebook and the embed: one control
  height for the whole viewer, stated once, and one spacing rhythm. The column
  is 88px closed against 106, and an embed keeps its own spacing in a host page
  that styles bare elements.
- **The website is compacted and aligned**: the header block is 122px against
  161, the page is 79px shorter, and every block shares its edges — header,
  viewer row, sequence strip, control panel, PAE and scatter.
- **The object picker is not shown when there is one object** — where the
  picker is all its row holds, which is the notebook's shell and the embed's.
  The website keeps its row: there the picker sits beside Multi and the
  prev/next buttons, which stay useful with one object.
- **Chain mode in the sequence strip shows the selection.** It showed nothing
  at all, while the same click lit the structure.

**The selection mark follows the ribbon.** It joined consecutive residues with
a straight line, and a cartoon helix is a ribbon spiralling *through* those
residues — so the mark chorded the thing it was marking, by more than twice its
own width. `cartoon/geom.js` hands back the centre line it actually drew and
the mark traces that: on 4HHB's longest helix the path is 328 px against 295.8
of chords, which is the arc-over-chord ratio of a helical step. A tube is
unchanged, because a tube *is* the straight lines between its residues.

**Focus is a mode with a door at each end.** Entering remembers the selection,
every object's side chains, the slab and the camera's centre and zoom, then
clears the decorations so the session starts from the structure rather than
from whatever was left on screen — and focuses a selection that is already
there, since pressing Focus with something picked is asking to look at it
closer. Leaving puts it all back. Two things it deliberately does not take
back: an angle you turned to while inside (focus never rotates, so if it moved,
you moved it) and a selection mark you chose in there. The mark goes to
**Outline** while the mode is on, which is what that option was built for, and
the selection panel stays out of the way.

**How a selection is marked is a setting.** A `Sele` dropdown in the Style
panel, beside `Color` — **Highlight**, **Outline**, **None** — in the website,
the notebook and the embed alike. Highlight is the translucent band this has always drawn and
stays the default. Outline draws the same band and punches its middle out, so a
thin line traces the residue and the geometry inside is untouched: over a
yellow chain, where the band is a blot with the answer somewhere inside it,
that is the difference between marking a residue and covering it. None draws
nothing.

`docs/SELECTION_MARK.md` has the six treatments this was chosen from, with the
pictures and the constants, and the measurements that say the choice costs
nothing — 0.02 ms between the cheapest and the dearest, on an operation that
runs in a tenth of a millisecond.

### Faster

**A side-chain click costs half what it did.** (And half of that again: see
the three-part mesh below.) Showing side chains appends
positions, which changed every term of the GPU mesh signature and rebuilt the
whole cartoon — 8,514 ribbon faces recomputed to draw 182 new stick ones. The
mesh is built in two halves now, ribbon and sticks, concatenated into one
buffer so the draw path is untouched, and the ribbon half is reused whenever
its own faces are unchanged. On 4HHB with 40 side chains, toggled eight times:
**62 ms a click to 32.5**. It applies to every way side chains change — the
Style panel toggle, `showSidechains` from the embed, `view(sidechains=True)`,
a contact, and Focus, which is the one that does it most often.

The mesh is three parts now — the ribbon, the ligands and contacts, and the
side chains — and a click rebuilds only the last: on 4HHB that is 198 faces
instead of 10,336, taking the stick work from 12 ms to 2.5 and the whole click
to a **31 ms median**. A heme is 1,822 faces and was being rebuilt on every
side-chain click purely because it is made of sticks.

What it costs is 260 pixels of 357,604, all of them where a stick surface meets
the ribbon it grows out of: the halves arrive ribbon-then-stick rather than
interleaved by depth, and a tie under `depthFunc(LESS)` goes to whoever drew
first. The winner there was already arbitrary. Indistinguishable at 6x.

**A large structure allocates about half the garbage it used to.** The prim
list is dropped as it is read rather than one pass later, the per-edge objects
and per-face corner arrays are flat typed arrays, and the edge table's Map of
Maps is one map and a chain. On a ribosome the live peak in a build is 288 MB
against 259, and the heap as V8 actually lets it grow went from ~450 MB to
214-436. Nothing about the mesh changed: verified by an order-independent
digest of the fill and the outline, identical across five structures, because
two runs of the same code differ in a third of their pixels on this path and
images cannot settle the question.

### Internal

One translation per surface, not three. A selector becomes a slab in
`renderer.clipTo` and a camera move in `py2dmolOrient.orientTo`, and the
website, the notebook and the embed all call those - each had grown its own
copy, and the copies had already diverged over whether the Clip button gets
re-synced.

The JavaScript is 26 source files under `src/`, merged per target by
`tools/bundle.py` — one manifest, from which every other file list is derived.
`dev.html` is generated from it and `bundle.py check` fails if it has drifted.
The suite is three lanes and about 40 seconds: node checks, headless-Chrome
probes, and GPU probes that run alone because they measure time.

## 1.6.5 and earlier

See the commit history.
