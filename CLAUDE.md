# Protein Inspector — where to look

- **Using the tool** (layers, morph, partners, the About tab, reading a bundle
  back): `AGENTS.md`. That is the one to read first, whichever agent you are.
- **The Claude Science skill**: `claude-science-skill/SKILL.md` is the source of
  the published skill; `claude-science-skill/README.md` explains why the skill
  and the engine ship separately.
- **The renderer's own public surface**: docstrings on
  `protein_inspector.inspect_structure` and
  `protein_inspector.render_inspection_bundle`.
- **Internals of the vendored py2Dmol viewer** — the `src/` file map, bundle
  targets, paint order, GPU lifecycle, selection marks:
  `docs/viewer_internals.md`. You need it to change how the picture is drawn,
  and not before.

This file used to be the viewer file map itself, which meant the one document
an agent is pointed at by convention was 143 kB about how the renderer is
built, and nothing about what the tool does.
- Comments in `src/` cite `tests/*` files this repo does not ship; they are upstream py2Dmol's. See `docs/viewer_internals.md`.
