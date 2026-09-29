import hashlib
import json
from pathlib import Path

from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

from bindos_structure_inspector import render_inspection_bundle


def _fixture(path: Path) -> None:
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for index in range(1, 6):
        residue = Residue.Residue((" ", index, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", (float(index * 3), float(index % 2), 0.0), 0.0, 1.0, " ", "CA", index, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))


def test_local_bundle_has_partner_specific_layers_and_hashed_exports(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {"annotation_id": "x", "kind": "partner_contact", "label": "Partner X contact", "partner_id": "Partner X", "residue": {"component_id": "target", "canonical_position": 2, "chain_id": "A", "author_residue_number": 2}, "resolved": True, "evidence_ids": ["e1"], "method": "distance"},
            {"annotation_id": "y", "kind": "partner_contact", "label": "Partner Y contact", "partner_id": "Partner Y", "residue": {"component_id": "target", "canonical_position": 4, "chain_id": "A", "author_residue_number": 4}, "resolved": True, "evidence_ids": ["e2"], "method": "distance"},
        ],
    }
    result = render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=digest, inspection_manifest=manifest, output_dir=str(tmp_path / "out"))
    state = json.loads(Path(next(item["path"] for item in result["artifacts"] if item["kind"] == "viewer_state")).read_text())
    assert {item["layer_id"] for item in state["layers"]} == {"partner_contact:Partner X", "partner_contact:Partner Y"}
    assert {item["kind"] for item in result["artifacts"]} == {"html", "viewer_state", "svg", "png"}
    assert all(Path(item["path"]).is_file() and len(item["sha256"]) == 64 for item in result["artifacts"])
    html = Path(next(item["path"] for item in result["artifacts"] if item["kind"] == "html")).read_text()
    assert "setResidueSelection" in html
    assert "bindos-residue-details" in html
    assert "syncVisibleLayers" in html
    assert "py2dmol-residue-selection-change" in html
    assert "window.bindosInspection" in html
    # Partner layers remain independently addressable in both the annotation
    # list and the 3D selection overlay rather than collapsing by kind.
    assert 'data-layer="partner_contact:Partner X"' in html
    assert 'data-layer="partner_contact:Partner Y"' in html


def test_custom_layer_colors_resolved_residues(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {"annotation_id": "hotspot", "kind": "custom", "label": "Pocket hotspot", "layer_id": "binding-hotspots", "layer_label": "Binding hotspots", "color": "#ef4444", "residue": {"component_id": "target", "canonical_position": 3, "chain_id": "A", "author_residue_number": 3}, "resolved": True, "evidence_ids": [], "method": "model"},
        ],
    }
    result = render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(), inspection_manifest=manifest, output_dir=str(tmp_path / "out"))
    state = json.loads(Path(next(item["path"] for item in result["artifacts"] if item["kind"] == "viewer_state")).read_text())
    assert state["layers"] == [{"layer_id": "binding-hotspots", "label": "Binding hotspots", "color": "#ef4444", "annotation_ids": ["hotspot"], "visible": True}]
    # Advanced colour is {"type": "advanced", "value": {"position": {...}}} -- the
    # shape py2Dmol's renderer reads (src/core/mol.js) and renderer.setColor writes.
    # Position keys are 0-based object indices, so canonical_position 3 is key "2".
    assert state["viewer"]["objects"][0]["color"]["value"]["position"]["2"] == "#ef4444"


def test_resolved_annotation_requires_a_matching_residue_address(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {"annotation_id": "bad", "kind": "mutation", "label": "Bad address", "residue": {"component_id": "target", "canonical_position": 9, "chain_id": "A", "author_residue_number": 9}, "resolved": True, "evidence_ids": [], "method": "model"},
        ],
    }
    try:
        render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(), inspection_manifest=manifest, output_dir=str(tmp_path / "out"))
    except ValueError as error:
        assert "does not match a residue" in str(error)
    else:
        raise AssertionError("expected an invalid resolved residue to fail")


def _run_dom_harness(html_path, spec):
    """Execute the emitted behaviour script in Node against a stub DOM."""
    import subprocess
    out = subprocess.run(
        ["node", str(Path(__file__).parent / "bindos_inspection_dom.js"), str(html_path), json.dumps(spec)],
        capture_output=True, text=True)
    if out.returncode != 0:
        raise AssertionError(f"harness failed ({out.returncode}):\n{out.stderr}")
    return json.loads(out.stdout)


def _two_layer_bundle(tmp_path, highlight=None):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {"annotation_id": "a1", "kind": "custom", "label": "Layer A res 1", "layer_id": "layer-a",
             "layer_label": "Layer A", "color": "#111111", "resolved": True, "evidence_ids": [], "method": "rule",
             "residue": {"component_id": "target", "canonical_position": 1, "chain_id": "A", "author_residue_number": 1}},
            {"annotation_id": "a2", "kind": "custom", "label": "Layer A res 2", "layer_id": "layer-a",
             "layer_label": "Layer A", "color": "#111111", "resolved": True, "evidence_ids": [], "method": "rule",
             "residue": {"component_id": "target", "canonical_position": 2, "chain_id": "A", "author_residue_number": 2}},
            {"annotation_id": "b2", "kind": "custom", "label": "Layer B res 2", "layer_id": "layer-b",
             "layer_label": "Layer B", "color": "#222222", "resolved": True, "evidence_ids": [], "method": "rule",
             "residue": {"component_id": "target", "canonical_position": 2, "chain_id": "A", "author_residue_number": 2}},
        ],
    }
    if highlight is not None:
        manifest["highlight"] = highlight
    result = render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
                                      inspection_manifest=manifest, output_dir=str(tmp_path / "out"))
    return Path(next(i["path"] for i in result["artifacts"] if i["kind"] == "html"))


def test_layer_checkboxes_repaint_the_structure(tmp_path):
    """Unchecking a layer must remove its colour, not merely hide panel rows.

    Regression test for two faults that every substring assertion missed: the
    script resolved the renderer through a bare `viewerApi` (module-scoped in
    src/app/main.js, so always undefined from a separate script tag) and it
    never rewrote the colour map at all.
    """
    html_path = _two_layer_bundle(tmp_path)
    spec = {
        "chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5],
        "steps": [
            {"name": "hide_a", "set": {"layer-a": False}},
            {"name": "hide_both", "set": {"layer-a": False, "layer-b": False}},
            {"name": "show_both", "set": {"layer-a": True, "layer-b": True}},
        ],
    }
    report = _run_dom_harness(html_path, spec)
    base = "#b9bfc7"
    lit = {"0": "#111111", "1": "#222222", "2": base, "3": base, "4": base}
    # Manifest order decides overlap: layer-b is emitted second, so it owns residue 2.
    # Everything else carries the neutral base rather than the viewer's colour mode.
    assert report["initial"] == lit
    assert report["renders"] >= 1
    # Layer A off: its exclusive residue falls back to the base, layer B's survives.
    assert report["steps"]["hide_a"]["paint"] == {**lit, "0": base}
    # Everything off: a fully neutral structure. NOT None — dropping the map here
    # is what used to expose the rainbow that `auto` resolves to on one chain.
    assert report["steps"]["hide_both"]["paint"] == {str(i): base for i in range(5)}
    # ...and back, byte for byte.
    assert report["steps"]["show_both"]["paint"] == lit


def test_layers_repaint_when_the_viewer_starts_late(tmp_path):
    """The bundle builds its renderer on load, after this script runs."""
    report = _run_dom_harness(_two_layer_bundle(tmp_path), {
        "chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5], "rafDelay": 5})
    base = "#b9bfc7"
    assert report["initial"] == {"0": "#111111", "1": "#222222", "2": base, "3": base, "4": base}
    assert report["frames"] >= 5


def test_inspection_script_uses_the_public_viewer_registry(tmp_path):
    """`viewerApi` is module-scoped in the bundle and unreachable from here."""
    html = _two_layer_bundle(tmp_path).read_text()
    script = html[html.index("<script>(function(){"):]
    assert "py2dmol_viewers" in script
    assert "viewerApi" not in script


def test_color_mode_does_not_also_halo_the_layer_residues(tmp_path):
    """Colour and halo are alternatives; a layer must not use both channels."""
    report = _run_dom_harness(_two_layer_bundle(tmp_path), {
        "chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5],
        "steps": [{"name": "hide_a", "set": {"layer-a": False}}]})
    base = "#b9bfc7"
    assert report["initial"] == {"0": "#111111", "1": "#222222", "2": base, "3": base, "4": base}
    # Layer visibility must leave the selection overlay alone in colour mode.
    assert report["steps"]["hide_a"]["selection"] == 0


def test_halo_mode_selects_instead_of_painting(tmp_path):
    # Halo mode deliberately paints NO base: it exists to keep the structure's
    # own colouring (pLDDT, chain) intact and annotate on top of it.
    html_path = _two_layer_bundle(tmp_path, highlight="halo")
    report = _run_dom_harness(html_path, {
        "chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5],
        "steps": [{"name": "hide_a", "set": {"layer-a": False}},
                  {"name": "hide_both", "set": {"layer-a": False, "layer-b": False}}]})
    assert report["initial"] is None                       # nothing painted
    assert report["steps"]["hide_a"]["selection"] == 1     # only layer-b's residue
    assert report["steps"]["hide_both"]["selection"] == 0
    state = json.loads(Path(str(html_path).replace(".html", ".viewer.json")).read_text())
    assert state["highlight"] == "halo"
    assert state["viewer"]["objects"][0].get("color") in (None, {})


def test_unknown_highlight_style_is_rejected(tmp_path):
    import pytest
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    with pytest.raises(ValueError, match="highlight must be one of"):
        render_inspection_bundle(
            mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
            inspection_manifest={"schema_version": "bindos-inspection-manifest-1",
                                 "annotations": [], "highlight": "both"},
            output_dir=str(tmp_path / "out2"))


def test_every_export_is_self_contained(tmp_path):
    """Render TWICE in one process; both bundles must carry the library.

    py2Dmol shares its library across view() calls in a process — right for
    notebook cells, fatal for a file on disk, because the borrower looks for a
    lender that does not exist and dies with "the viewer library never loaded".
    Rendering once cannot catch it: the first export is always fine.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {
                "annotation_id": "a-1",
                "kind": "custom",
                "label": "probe",
                "resolved": True,
                "residue": {
                    "component_id": "t",
                    "canonical_position": 2,
                    "chain_id": "A",
                    "author_residue_number": 2,
                },
                "evidence_ids": [],
                "method": "test",
            }
        ],
    }
    sizes = []
    for index in (1, 2):
        out = tmp_path / f"out{index}"
        render_inspection_bundle(
            mmcif_path=str(cif),
            mmcif_sha256=digest,
            inspection_manifest=manifest,
            output_dir=str(out),
            output_name="inspection",
        )
        page = (out / "inspection.html").read_text(encoding="utf-8")
        assert "initializePy2DmolViewer" in page, f"export {index} has no viewer library"
        assert "BroadcastChannel('py2dmol_lib')" not in page, (
            f"export {index} borrows the library instead of inlining it"
        )
        assert "the viewer library never" in page, (
            f"export {index} lost the bootstrap diagnostic"
        )
        sizes.append(len(page))
    assert min(sizes) > 400_000, f"a bundle is too small to contain the library: {sizes}"


def test_base_color_is_painted_and_survives_unticking_every_layer(tmp_path):
    """A neutral base must be painted explicitly, not left to a colour mode.

    py2Dmol has no flat-grey mode and does not reject an unknown one — ui.js
    falls back to 'auto', i.e. rainbow on a single chain. So the backdrop has
    to come from the per-position map, and it has to be there when every layer
    is off, which is exactly when the map used to be dropped for null.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
        "annotations": [
            {
                "annotation_id": "a-1",
                "kind": "custom",
                "label": "probe",
                "color": "#dc2626",
                "resolved": True,
                "residue": {
                    "component_id": "t",
                    "canonical_position": 2,
                    "chain_id": "A",
                    "author_residue_number": 2,
                },
                "evidence_ids": [],
                "method": "test",
            }
        ],
    }
    result = render_inspection_bundle(
        mmcif_path=str(cif),
        mmcif_sha256=digest,
        inspection_manifest=manifest,
        output_dir=str(tmp_path / "out"),
    )
    state_path = next(i["path"] for i in result["artifacts"] if i["kind"] == "viewer_state")
    state = json.loads(Path(state_path).read_text())

    assert state["base_color"].startswith("#"), "no base_color published to the page"
    base = state["base_color"]

    paint = state["viewer"]["objects"][0]["color"]["value"]["position"]
    # every residue of the 5-residue fixture is painted, not only the annotated one
    assert len(paint) == 5, f"base not painted over the whole chain: {paint}"
    assert paint["1"] == "#dc2626", "layer colour did not land on its residue"
    assert {v for k, v in paint.items() if k != "1"} == {base}, "unannotated residues are not the base"

    # the viewer's own mode must be one py2Dmol actually honours
    from py2Dmol.viewer import VALID_COLOR_MODES
    assert state["viewer"]["config"]["color"]["mode"] in VALID_COLOR_MODES

    page = Path(next(i["path"] for i in result["artifacts"] if i["kind"] == "html")).read_text()
    assert "var BASE=s.base_color" in page, "JS repaint does not seed the base colour"
