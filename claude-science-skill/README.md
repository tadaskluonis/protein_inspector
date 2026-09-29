# bindos-inspector — Claude Science skill

The Claude Science packaging of `bindos_structure_inspector`. Two files do the
work:

| file | role |
| --- | --- |
| `SKILL.md` | what the agent reads — workflow, layer conventions, morph caveats |
| `kernel.py` | sidecar helpers auto-loaded into the agent's kernel |

## Why it is split from the engine

The engine vendors py2Dmol: four minified bundles, 2.2 MB of JavaScript. A
published Claude Science skill is shown to the user as a card listing every
file in the bundle, and publish is refused for binary or oversized files. So
the skill stays small and text-only, and installs the engine:

```bash
pip install git+https://github.com/profdocpizza/bindos-structure-inspector@v1.13
```

For a local checkout, `BINDOS_INSPECTOR_HOME=/path/to/repo` or leave it at
`~/code/bindos-structure-inspector`, which `inspector_engine()` probes.

Pin a **tag**, not `main`: a published skill outlives the working copy, and an
unpinned install lets a later push change what an old conversation runs.

## What the sidecar adds over the raw engine

- `canonical_position` is derived from the file instead of being asked for.
  The engine treats it as authoritative and a wrong value colours the wrong
  residue without raising; this is the single most common way to get a bundle
  that is confidently incorrect.
- Layers are `{"id", "label", "color", "residues"}` with `residues` as a range,
  list or `{number: note}` mapping, rather than one fully-specified annotation
  dict per residue.
- `.pdb` in, mmCIF out, via gemmi.
- `report["unmapped_residues"]` names the residues that were not in the model,
  which the engine otherwise skips in silence.

## Publishing

```python
host.skills.edit("bindos-inspector", "SKILL.md", open("SKILL.md").read())
host.skills.edit("bindos-inspector", "kernel.py", open("kernel.py").read())
host.skills.publish("bindos-inspector")
```

## Tests

```bash
python -m pytest claude-science-skill/tests -q
```

They run against the sibling engine checkout, so they also catch a sidecar that
has drifted from the renderer's signature.
