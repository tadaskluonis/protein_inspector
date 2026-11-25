"""The API most callers want: layers in, one self-contained bundle out.

`render_inspection_bundle` in .renderer is the strict surface -- it takes a
validated manifest in which every annotation carries a `canonical_position`
that is NOT the author residue number. Getting that wrong colours the wrong
residue and raises nothing. So nothing here asks the caller for it:
`inspect_structure` reads the ordering out of the file.
"""
from __future__ import annotations

import os
from pathlib import Path

from .renderer import _ca_trace, _trace, read_inspection_bundle, render_inspection_bundle

DEFAULT_PALETTE = ("#dc2626", "#2563eb", "#059669", "#9333ea", "#d97706", "#0891b2",
                   "#db2777", "#65a30d", "#0ea5e9", "#f59e0b")


def to_mmcif(structure, out_dir=None):
    """Return a path to an .cif for `structure`, converting a .pdb if needed."""
    import gemmi
    path = os.path.abspath(os.path.expanduser(str(structure)))
    if path.lower().endswith((".cif", ".mmcif")):
        return path
    target_dir = out_dir or os.path.dirname(path) or "."
    out = os.path.join(target_dir, os.path.splitext(os.path.basename(path))[0] + ".cif")
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.make_mmcif_document().write_file(out)
    return out


def _chain_list(chain):
    """`None` / "A" / ["A", "B"] -> a list of chain ids, possibly empty."""
    if chain is None:
        return []
    if isinstance(chain, str):
        return [chain]
    return [str(c) for c in chain]


def residue_order(cif_path, chain=None):
    """[(chain_id, author_resnum)] in file order — the canonical_position basis.

    `chain` takes one chain id or several; `None` is the whole file.
    """
    coords, chains, numbers = _trace(Path(cif_path))
    del coords
    pairs = list(zip(chains, [int(n) for n in numbers]))
    wanted = _chain_list(chain)
    if wanted:
        pairs = [p for p in pairs if p[0] in set(wanted)]
    return pairs


def residue_components(cif_path, chain=None):
    """{(chain_id, author_resnum): component_id} -- MET, ALA, DA, HOH.

    component_id is a REQUIRED field of every resolved annotation, and it was
    being filled with the string "target" for every residue this API has ever
    annotated. The manifest validator only checks that the field is present,
    so nothing complained, and the residue card had the layer's name to show
    and never the residue's own identity.
    """
    out = {}
    for (chain_id, number), (name, _xyz) in _ca_trace(Path(cif_path)).items():
        if chain is not None and chain_id != chain:
            continue
        out[(chain_id, int(number))] = str(name)
    return out


def build_manifest(cif_path, layers, chain=None, highlight="color", about=None,
                   base_mode="chain", method="structure-derived annotation"):
    """Turn simple layer dicts into a validated inspection manifest.

    Each layer: {"id", "label", "residues", "color"?, "kind"?, "note"?,
    "description"?, "visible"?}. `description` is one sentence saying what
    the layer IS -- it becomes the row's tooltip and an info glyph that opens
    it inline, so a reader meeting "AS1" does not have to go to a tab to find
    out. `visible: False` means the layer is built and listed but
    starts unchecked, so the opening view can be the one thing that matters.
    `residues` is a list of author residue numbers, or {resnum: note}.
    Later layers paint over earlier ones, so order background-first.

    `chain` is the default a layer with no `"chain"` key of its own paints on:
    one chain id, or SEVERAL. A list means every one of them is the target and
    wears the layer -- which is what a homo-oligomer is, and what the file's
    other protomers need when they are the biological unit rather than
    something the reader should be able to switch off. Keep `partner_chains`
    for what it is for: chains that are genuinely not the target.
    """
    # EVERY CHAIN IN THE FILE, not just the target. This was built from the
    # target chain alone, so a layer naming any other chain -- the documented
    # {"chain": "B"} key -- found no address, landed in `unmapped_residues`
    # whole, and rendered a bundle byte-identical to one that had never been
    # given the key. A silent no-op is the expensive kind.
    order = residue_order(cif_path)
    components = residue_components(cif_path)
    canonical, seen = {}, {}
    for index, (chain_id, number) in enumerate(order):
        seen[chain_id] = seen.get(chain_id, 0) + 1
        canonical[(chain_id, number)] = seen[chain_id]
    present = {chain_id for chain_id, _ in order}
    default_chains = _chain_list(chain) or [order[0][0]]
    unknown = [c for c in default_chains if c not in present]
    if unknown:
        raise ValueError(
            "chain=%r names %s, which this file does not have (chains: %s)"
            % (chain, ", ".join(unknown), ", ".join(sorted(present))))
    annotations, missing = [], []
    for position, layer in enumerate(layers):
        residues = layer["residues"]
        items = residues.items() if hasattr(residues, "items") else [(r, None) for r in residues]
        layer_chains = _chain_list(layer.get("chain")) or default_chains
        # A CHAIN THE ENGINE WOULD IGNORE IS AN ERROR, not an empty result. A
        # typo'd or stale chain key used to come back as every residue of the
        # layer in `unmapped_residues`, which reads exactly like "the model
        # does not have those residues".
        missing_chains = [c for c in layer_chains if c not in present]
        if missing_chains:
            raise ValueError(
                "layer %r names chain %s, which this file does not have (chains: %s)"
                % (layer["id"], ", ".join(missing_chains), ", ".join(sorted(present))))
        # ONE ANNOTATION PER (CHAIN, RESIDUE). A layer that names several
        # chains is the same fact about each of them -- the protomers of a
        # homo-oligomer present the same site -- so each gets its own address
        # rather than the first one standing in for all.
        for target_chain in layer_chains:
            for number, note in sorted(items):
                key = (target_chain, int(number))
                if key not in canonical:
                    missing.append((layer["id"], target_chain, int(number)))
                    continue
                item = {
                    "annotation_id": "%s-%s%d" % (layer["id"], target_chain, int(number)),
                    # A layer can ask to start switched off: the reader still has
                    # it, but the first view is not all sixteen things at once.
                    # Stamped on the annotations because the manifest has no
                    # separate layer list.
                    **({"layer_visible": False} if layer.get("visible", True) is False else {}),
                    "kind": layer.get("kind", "custom"),
                    "label": note or layer.get("note") or layer["label"],
                    "layer_id": layer["id"],
                    "layer_label": layer["label"],
                    # WHAT THE LAYER IS, not only what it is called. "AS1" is a
                    # name; "the first amphipathic helix, which the transmembrane
                    # segment feeds into" is the thing the reader needs, and it
                    # belongs on the row rather than in a tab they have to find.
                    **({"layer_description": layer["description"]} if layer.get("description") else {}),
                    "color": layer.get("color") or DEFAULT_PALETTE[position % len(DEFAULT_PALETTE)],
                    "resolved": True,
                    "evidence_ids": list(layer.get("evidence_ids", [])),
                    "method": layer.get("method", method),
                    "residue": {"component_id": (layer.get("component_id")
                                                or components.get(key) or "target"),
                                "canonical_position": canonical[key],
                                "chain_id": target_chain,
                                "author_residue_number": int(number)},
                }
                if item["kind"] == "partner_contact":
                    item["partner_id"] = layer.get("partner_id", layer["id"])
                annotations.append(item)
    manifest = {"schema_version": "protein-inspector-manifest-1", "highlight": highlight,
                "base_mode": base_mode, "annotations": annotations}
    if about:
        manifest["about"] = about
    return manifest, missing


def inspect_structure(structure, layers=None, out="inspection.html", chain=None,
                      about=None, conformers=None, reference_label="reference",
                      morph_mapping="intersection", morph_steps=12,
                      partner_chains=None, partner_label="partners",
                      layers_on_partners="auto",
                      highlight="color", display=None, extras=False,
                      interactions=None):
    """Render one self-contained annotated HTML bundle. Returns the render report.

    `chain` names the target: one chain id, or a list of them. Every chain
    named is drawn always and carries the layers that do not name a chain of
    their own, which is what the protomers of a homo-oligomer are — use it
    rather than `partner_chains` when the assembly IS the thing, since a
    partner is something the reader can switch off. A layer's own `"chain"`
    key (one id or a list) overrides it.

    `conformers` is a list of {"path", "label"?, "partner_chains"?,
    "partner_label"?} — any number. Only the endpoints are written to the
    file; the page interpolates between them.

    `partner_chains` names chains in the SAME file that are binding partners
    rather than the target; they get a Show-partners toggle and are excluded
    from the morph residue mapping. The toggle is named by `partner_label`,
    or by the conformer's own `partner_label` while that state is on screen —
    use it when the states carry different partners ("other protomers" in the
    apo trimer, "TNFR1" in the bound one).

    `layers_on_partners` repeats the layer colours onto partner chains that
    are COPIES of the target: "auto" (default) detects them by sequence, so a
    homo-oligomer is painted on every subunit and a receptor or Fab is not,
    True takes every partner chain, False keeps the colours on the annotated
    chain alone. The painted chains come back as report["partners"]["painted"].
    """
    import hashlib
    out_dir = os.path.dirname(os.path.abspath(out)) or "."
    cif = to_mmcif(structure, out_dir)
    manifest, missing = build_manifest(cif, layers or [], chain=chain,
                                       highlight=highlight, about=about)
    prepared = None
    if conformers:
        prepared = [{"path": to_mmcif(c["path"], out_dir), "label": c.get("label"),
                     "partner_chains": list(c.get("partner_chains") or ()) or None,
                     "partner_label": c.get("partner_label")}
                    for c in conformers]
    report = render_inspection_bundle(
        mmcif_path=cif,
        mmcif_sha256=hashlib.sha256(open(cif, "rb").read()).hexdigest(),
        inspection_manifest=manifest,
        output_dir=out_dir,
        output_name=os.path.splitext(os.path.basename(out))[0],
        display_options=display or {"width": 900, "height": 760},
        conformers=prepared,
        morph_mapping=morph_mapping,
        morph_reference_label=reference_label,
        partner_chains=list(partner_chains) if partner_chains else None,
        partner_label=partner_label,
        layers_on_partners=layers_on_partners,
        morph_steps=morph_steps,
        extras=extras,
        interactions=interactions,
    )
    report["path"] = report["artifacts"][0]["path"]
    report["unmapped_residues"] = missing
    # ...AND PER STATE, because a morph's states need not model the same
    # residues. A layer is painted on the residues the conformation on screen
    # actually has; the ones it does not have are listed here under that
    # state's label, as (layer_id, chain, author_number) — so "the loop loses
    # 12 residues in state 2" is something a caller reads rather than
    # something a reader notices by toggling the button.
    if report.get("morph") and report["morph"].get("conformers"):
        report["unmapped_residues_by_state"] = {
            conformer["label"]: [tuple(entry) for entry
                                 in conformer.get("layer_residues_unmapped") or []]
            for conformer in report["morph"]["conformers"]
        }
    return report


def inspection_table(html_path):
    """Flat per-residue rows (chain, number, colour, layers, labels) from a bundle."""
    return read_inspection_bundle(html_path)["residues"]
