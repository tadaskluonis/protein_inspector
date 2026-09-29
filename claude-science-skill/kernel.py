"""Kernel helpers for the bindos-inspector skill.

The engine (`bindos_structure_inspector`) speaks a strict manifest schema with a
`canonical_position` that is NOT the author residue number.  Getting that wrong
colours the wrong residue without raising, so nothing here asks the caller for
it: `inspect_structure` reads the ordering out of the file.
"""
import os
import sys

DEFAULT_PALETTE = ("#dc2626", "#2563eb", "#059669", "#9333ea", "#d97706", "#0891b2",
                   "#db2777", "#65a30d", "#0ea5e9", "#f59e0b")
SEARCH_PATHS = ("~/code/bindos-structure-inspector", "~/bindos-structure-inspector")


def inspector_engine(path=None):
    """Import the renderer, adding a local checkout to sys.path if needed."""
    candidates = []
    if path:
        candidates.append(path)
    if os.environ.get("BINDOS_INSPECTOR_HOME"):
        candidates.append(os.environ["BINDOS_INSPECTOR_HOME"])
    candidates.extend(os.path.expanduser(p) for p in SEARCH_PATHS)
    for candidate in candidates:
        if candidate and os.path.isdir(candidate) and candidate not in sys.path:
            sys.path.insert(0, candidate)
    try:
        import bindos_structure_inspector as engine
    except ImportError as exc:
        raise ImportError(
            "bindos_structure_inspector not importable. Either "
            "`pip install git+https://github.com/profdocpizza/bindos-structure-inspector` "
            "or point BINDOS_INSPECTOR_HOME at a local checkout."
        ) from exc
    return engine


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


def residue_order(cif_path, chain=None):
    """[(chain_id, author_resnum)] in file order — the canonical_position basis."""
    engine = inspector_engine()
    from pathlib import Path
    coords, chains, numbers = engine.renderer._trace(Path(cif_path))
    del coords
    pairs = list(zip(chains, [int(n) for n in numbers]))
    if chain is not None:
        pairs = [p for p in pairs if p[0] == chain]
    return pairs


def build_manifest(cif_path, layers, chain=None, highlight="color", about=None,
                   base_mode="chain", method="structure-derived annotation"):
    """Turn simple layer dicts into a validated inspection manifest.

    Each layer: {"id", "label", "residues", "color"?, "kind"?, "note"?}.
    `residues` is a list of author residue numbers, or {resnum: note}.
    Later layers paint over earlier ones, so order background-first.
    """
    order = residue_order(cif_path, chain)
    canonical, seen = {}, {}
    for index, (chain_id, number) in enumerate(order):
        seen[chain_id] = seen.get(chain_id, 0) + 1
        canonical[(chain_id, number)] = seen[chain_id]
    annotations, missing = [], []
    for position, layer in enumerate(layers):
        residues = layer["residues"]
        items = residues.items() if hasattr(residues, "items") else [(r, None) for r in residues]
        target_chain = layer.get("chain", chain) or order[0][0]
        for number, note in sorted(items):
            key = (target_chain, int(number))
            if key not in canonical:
                missing.append((layer["id"], int(number)))
                continue
            item = {
                "annotation_id": "%s-%s%d" % (layer["id"], target_chain, int(number)),
                "kind": layer.get("kind", "custom"),
                "label": note or layer.get("note") or layer["label"],
                "layer_id": layer["id"],
                "layer_label": layer["label"],
                "color": layer.get("color") or DEFAULT_PALETTE[position % len(DEFAULT_PALETTE)],
                "resolved": True,
                "evidence_ids": list(layer.get("evidence_ids", [])),
                "method": layer.get("method", method),
                "residue": {"component_id": layer.get("component_id", "target"),
                            "canonical_position": canonical[key],
                            "chain_id": target_chain,
                            "author_residue_number": int(number)},
            }
            if item["kind"] == "partner_contact":
                item["partner_id"] = layer.get("partner_id", layer["id"])
            annotations.append(item)
    manifest = {"schema_version": "bindos-inspection-manifest-1", "highlight": highlight,
                "base_mode": base_mode, "annotations": annotations}
    if about:
        manifest["about"] = about
    return manifest, missing


def inspect_structure(structure, layers=None, out="inspection.html", chain=None,
                      about=None, conformers=None, reference_label="reference",
                      morph_mapping="intersection", morph_steps=12,
                      highlight="color", display=None, extras=False):
    """Render one self-contained annotated HTML bundle. Returns the render report.

    `conformers` is a list of {"path", "label"?} — any number. Only the
    endpoints are written to the file; the page interpolates between them.
    """
    import hashlib
    engine = inspector_engine()
    out_dir = os.path.dirname(os.path.abspath(out)) or "."
    cif = to_mmcif(structure, out_dir)
    manifest, missing = build_manifest(cif, layers or [], chain=chain,
                                       highlight=highlight, about=about)
    prepared = None
    if conformers:
        prepared = [{"path": to_mmcif(c["path"], out_dir), "label": c.get("label")}
                    for c in conformers]
    report = engine.render_inspection_bundle(
        mmcif_path=cif,
        mmcif_sha256=hashlib.sha256(open(cif, "rb").read()).hexdigest(),
        inspection_manifest=manifest,
        output_dir=out_dir,
        output_name=os.path.splitext(os.path.basename(out))[0],
        display_options=display or {"width": 900, "height": 760},
        conformers=prepared,
        morph_mapping=morph_mapping,
        morph_reference_label=reference_label,
        morph_steps=morph_steps,
        extras=extras,
    )
    report["path"] = report["artifacts"][0]["path"]
    report["unmapped_residues"] = missing
    return report


def inspection_table(html_path):
    """Flat per-residue rows (chain, number, colour, layers, labels) from a bundle."""
    engine = inspector_engine()
    return engine.read_inspection_bundle(html_path)["residues"]
