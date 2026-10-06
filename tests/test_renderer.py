import hashlib
import json
import numpy as np
import pytest
from pathlib import Path

from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

from protein_inspector import read_inspection_bundle, render_inspection_bundle
from protein_inspector.renderer import _position_types


def _fixture(path: Path, n: int = 5, bend: float = 0.0) -> None:
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for index in range(1, n + 1):
        residue = Residue.Residue((" ", index, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", (float(index * 3), float(index % 2) + bend * index, 0.0), 0.0, 1.0, " ", "CA", index, element="C"))
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
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "x", "kind": "partner_contact", "label": "Partner X contact", "partner_id": "Partner X", "residue": {"component_id": "target", "canonical_position": 2, "chain_id": "A", "author_residue_number": 2}, "resolved": True, "evidence_ids": ["e1"], "method": "distance"},
            {"annotation_id": "y", "kind": "partner_contact", "label": "Partner Y contact", "partner_id": "Partner Y", "residue": {"component_id": "target", "canonical_position": 4, "chain_id": "A", "author_residue_number": 4}, "resolved": True, "evidence_ids": ["e2"], "method": "distance"},
        ],
    }
    result = render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=digest, inspection_manifest=manifest, output_dir=str(tmp_path / "out"), extras=True)
    state = json.loads(Path(next(item["path"] for item in result["artifacts"] if item["kind"] == "viewer_state")).read_text())
    assert {item["layer_id"] for item in state["layers"]} == {"partner_contact:Partner X", "partner_contact:Partner Y"}
    assert {item["kind"] for item in result["artifacts"]} == {"html", "viewer_state", "svg", "png"}
    assert all(Path(item["path"]).is_file() and len(item["sha256"]) == 64 for item in result["artifacts"])
    html = Path(next(item["path"] for item in result["artifacts"] if item["kind"] == "html")).read_text()
    assert "setResidueSelection" in html
    assert "pinsp-residue-details" in html
    assert "syncVisibleLayers" in html
    assert "py2dmol-residue-selection-change" in html
    assert "window.proteinInspector" in html
    # Partner layers remain independently addressable in both the annotation
    # list and the 3D selection overlay rather than collapsing by kind.
    assert 'data-layer="partner_contact:Partner X"' in html
    assert 'data-layer="partner_contact:Partner Y"' in html


def test_custom_layer_colors_resolved_residues(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "hotspot", "kind": "custom", "label": "Pocket hotspot", "layer_id": "binding-hotspots", "layer_label": "Binding hotspots", "color": "#ef4444", "residue": {"component_id": "target", "canonical_position": 3, "chain_id": "A", "author_residue_number": 3}, "resolved": True, "evidence_ids": [], "method": "model"},
        ],
    }
    result = render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(), inspection_manifest=manifest, output_dir=str(tmp_path / "out"), extras=True)
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
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "bad", "kind": "mutation", "label": "Bad address", "residue": {"component_id": "target", "canonical_position": 9, "chain_id": "A", "author_residue_number": 9}, "resolved": True, "evidence_ids": [], "method": "model"},
        ],
    }
    try:
        render_inspection_bundle(mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(), inspection_manifest=manifest, output_dir=str(tmp_path / "out"), extras=True)
    except ValueError as error:
        assert "does not match a residue" in str(error)
    else:
        raise AssertionError("expected an invalid resolved residue to fail")


def _run_dom_harness(html_path, spec):
    """Execute the emitted behaviour script in Node against a stub DOM."""
    import subprocess
    out = subprocess.run(
        ["node", str(Path(__file__).parent / "inspection_dom.js"), str(html_path), json.dumps(spec)],
        capture_output=True, text=True)
    if out.returncode != 0:
        raise AssertionError(f"harness failed ({out.returncode}):\n{out.stderr}")
    return json.loads(out.stdout)


def _two_layer_bundle(tmp_path, highlight=None):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
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
                                      inspection_manifest=manifest, output_dir=str(tmp_path / "out"), extras=True)
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
    lit = {"0": "#111111", "1": "#222222"}
    # Manifest order decides overlap: layer-b is emitted second, so it owns residue 2.
    # Unannotated residues are left to the VIEWER's colour mode (chain, in neutral
    # greys) rather than painted over: the mode only decides the ones nobody spoke
    # for, so speaking for all of them would disable it.
    assert report["initial"] == lit
    assert report["renders"] >= 1
    # Layer A off: its exclusive residue falls back to the mode, layer B's survives.
    assert report["steps"]["hide_a"]["paint"] == {"1": "#222222"}
    # Everything off: the map is dropped and the chain palette shows through.
    assert report["steps"]["hide_both"]["paint"] is None
    # ...and back, byte for byte.
    assert report["steps"]["show_both"]["paint"] == lit


def test_layers_repaint_when_the_viewer_starts_late(tmp_path):
    """The bundle builds its renderer on load, after this script runs."""
    report = _run_dom_harness(_two_layer_bundle(tmp_path), {
        "chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5], "rafDelay": 5})
    assert report["initial"] == {"0": "#111111", "1": "#222222"}
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
    assert report["initial"] == {"0": "#111111", "1": "#222222"}
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
            inspection_manifest={"schema_version": "protein-inspector-manifest-1",
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
        "schema_version": "protein-inspector-manifest-1",
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
    """base_mode='custom' is OPT-IN; when asked for, it must be complete.

    py2Dmol has no flat-grey mode and does not reject an unknown one — ui.js
    falls back to 'auto', i.e. rainbow on a single chain. So the backdrop has
    to come from the per-position map, and it has to be there when every layer
    is off, which is exactly when the map used to be dropped for null.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "base_mode": "custom",
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
        extras=True,
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
    assert "if(baseMode==='custom')" in page, "JS repaint does not gate the base on custom mode"
    # There is deliberately NO base-colour control in this panel: colouring the
    # un-annotated structure is the viewer's Style panel's job, and a second
    # control for it is a second source of truth.
    assert 'id="pinsp-base-mode"' not in page and 'id="pinsp-base-color"' not in page


def test_default_render_is_one_file_with_about_tabs(tmp_path):
    """A bundle is one file. Extra context becomes a TAB, never a sibling file.

    The About body is agent-authored and often quotes fetched text, so it is
    sanitised: the element and its contents go, not just the tags.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "about": [
            {"title": "Read me", "body": "<h2>Heading</h2><p>Body <code>x</code></p>"},
            {"title": "Methods", "body": "plain one\n\nplain two"},
        ],
        "annotations": [
            {
                "annotation_id": "a-1", "kind": "custom", "label": "probe", "resolved": True,
                "residue": {"component_id": "t", "canonical_position": 2, "chain_id": "A",
                            "author_residue_number": 2},
                "evidence_ids": [], "method": "test",
            }
        ],
    }
    out = tmp_path / "out"
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=digest,
        inspection_manifest=manifest, output_dir=str(out),
    )
    assert sorted(f.name for f in out.iterdir()) == ["inspection.html"], "default render is not one file"
    assert [a["kind"] for a in result["artifacts"]] == ["html"]
    assert result["about_tabs"] == ["Read me", "Methods"]
    assert "size_warning" not in result

    page = (out / "inspection.html").read_text(encoding="utf-8")
    assert 'data-tab="structure"' in page and 'data-tab="about-0"' in page
    assert 'data-pane="about-1"' in page
    assert "<h2>Heading</h2>" in page                    # markup preserved
    assert "<p>plain one</p><p>plain two</p>" in page    # bare text becomes paragraphs

    # ...and extras=True still writes the audit set alongside it
    out2 = tmp_path / "out2"
    render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=digest,
        inspection_manifest=manifest, output_dir=str(out2), extras=True,
    )
    assert sorted(f.suffix for f in out2.iterdir()) == [".html", ".json", ".json", ".png", ".svg"]


def test_about_body_cannot_smuggle_script(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "about": "<p>ok</p><script>steal()</script><a href=\'javascript:go()\' onclick=\'go()\'>z</a>",
        "annotations": [],
    }
    out = tmp_path / "out"
    render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest=manifest, output_dir=str(out),
    )
    about = (out / "inspection.html").read_text(encoding="utf-8")
    about = about[about.index('data-pane="about-0"'):]
    assert "<p>ok</p>" in about
    assert "steal()" not in about, "script contents survived into the page"
    assert "onclick" not in about
    assert "javascript:" not in about


def test_base_mode_other_than_custom_leaves_the_colour_mode_alone(tmp_path):
    """A flat base is one CHOICE among the modes, not a floor under them.

    py2Dmol: an explicit per-position colour beats the mode, and "the mode only
    decides the ones nobody spoke for". Painting every position — what 1.7 did —
    therefore disabled rainbow/plddt/chain/ss outright. Under a real mode only
    the annotated residues may be painted.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()

    def render(base_mode):
        manifest = {
            "schema_version": "protein-inspector-manifest-1",
            "base_mode": base_mode,
            "annotations": [
                {
                    "annotation_id": "a-1", "kind": "custom", "label": "probe",
                    "color": "#dc2626", "resolved": True,
                    "residue": {"component_id": "t", "canonical_position": 2,
                                "chain_id": "A", "author_residue_number": 2},
                    "evidence_ids": [], "method": "test",
                }
            ],
        }
        out = tmp_path / f"out-{base_mode}"
        result = render_inspection_bundle(
            mmcif_path=str(cif), mmcif_sha256=digest, inspection_manifest=manifest,
            output_dir=str(out), extras=True,
        )
        state_path = next(i["path"] for i in result["artifacts"] if i["kind"] == "viewer_state")
        return json.loads(Path(state_path).read_text())

    rainbow = render("rainbow")
    assert rainbow["base_mode"] == "rainbow"
    paint = rainbow["viewer"]["objects"][0]["color"]["value"]["position"]
    assert paint == {"1": "#dc2626"}, (
        f"under a colour mode only annotated residues may be painted, got {paint}"
    )

    custom = render("custom")
    paint = custom["viewer"]["objects"][0]["color"]["value"]["position"]
    assert len(paint) == 5 and paint["1"] == "#dc2626", "custom mode must still paint a full base"


def test_base_color_must_be_a_hex_value(tmp_path):
    import pytest
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    with pytest.raises(ValueError, match="base_color must be"):
        render_inspection_bundle(
            mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
            inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                                 "annotations": [], "base_color": "grey"},
            output_dir=str(tmp_path / "out"),
        )


def test_bundle_reads_back_cleanly_without_a_browser(tmp_path):
    """The bundle must be machine-readable as a supported call, not by scraping.

    One file is right for a reader and useless to a program unless the
    machine-readable part is an interface. Everything a caller needs comes back
    from the HTML alone: layers, annotations, About bodies, and a flat
    per-residue table with coordinates.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "about": "context that must survive the round trip",
        "annotations": [
            {"annotation_id": "a-1", "kind": "custom", "label": "probe two",
             "layer_id": "L", "layer_label": "Layer L", "color": "#dc2626", "resolved": True,
             "residue": {"component_id": "t", "canonical_position": 2, "chain_id": "A",
                         "author_residue_number": 2},
             "evidence_ids": ["e"], "method": "test"},
        ],
    }
    out = tmp_path / "out"
    render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest=manifest, output_dir=str(out),
    )
    state = read_inspection_bundle(str(out / "inspection.html"))

    assert state["inspector_version"].startswith("protein-inspector-")
    assert state["source"]["sha256"] == hashlib.sha256(cif.read_bytes()).hexdigest()
    assert state["base_mode"] == "chain"   # defaults defer to the viewer Style panel
    assert [l["layer_id"] for l in state["layers"]] == ["L"]
    assert state["layers"][0]["n_annotations"] == 1
    assert state["about"][0]["title"] == "About"
    assert "must survive" in state["about"][0]["body_html"]

    rows = state["residues"]
    assert len(rows) == 5, "every modelled residue should appear once"
    assert [r["canonical_position"] for r in rows] == [1, 2, 3, 4, 5]
    hit = [r for r in rows if r["author_residue_number"] == 2][0]
    assert hit["layers"] == ["L"] and hit["labels"] == ["probe two"]
    assert hit["color"] == "#dc2626"
    assert hit["residue_name"] == "ALA"
    assert all(isinstance(r["x"], (int, float)) for r in rows), "coordinates did not survive"
    assert [r for r in rows if r["author_residue_number"] == 1][0]["layers"] == []


def test_reader_rejects_a_page_that_is_not_a_bundle(tmp_path):
    import pytest
    stray = tmp_path / "other.html"
    stray.write_text("<html><body>not a bundle</body></html>")
    with pytest.raises(ValueError, match="no protein_inspector inspection state"):
        read_inspection_bundle(str(stray))


# --- the viewer config the bundle emits ----------------------------------
# Lives here rather than beside the other chain-colour tests because it needs
# the renderer: it checks the whole path in the direction it runs, from the
# config this package emits through the normaliser in core/mol.js that reads
# it. That normaliser rebuilds `color` from a named whitelist and silently
# dropped `chain_palette` for a whole revision.

def test_bundle_chain_palette_survives_the_viewers_own_config_normalizer(tmp_path):
    quickjs = pytest.importorskip("quickjs", reason="QuickJS runs the real mol.js")
    from tests.test_chain_colors import BROWSER_STUB, MOL_JS

    cif = tmp_path / "f.cif"
    _fixture(cif)
    render_inspection_bundle(
        mmcif_path=str(cif),
        mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "out"),
        extras=True,
    )
    color = json.loads((tmp_path / "out" / "inspection.viewer.json").read_text())
    color = color["viewer"]["config"]["color"]
    assert color["chain_palette"] == "greys", "renderer stopped asking for greys"

    ctx = quickjs.Context()
    ctx.eval(BROWSER_STUB)
    ctx.eval(MOL_JS.read_text())
    survived = json.loads(ctx.eval(
        "JSON.stringify(normalizeConfig(" + json.dumps({"color": color}) + ").color)"))
    assert survived["chain_palette"] == "greys", "config drops the key on the way in"
    assert survived["mode"] == "chain"


def _morph_bundle(tmp_path, n_conformers=3, mapping="exact", steps=6, ragged=False):
    """Reference plus n-1 conformers that differ by a real internal change.

    A rigid translation would be removed by the superposition, so each
    conformer bends instead. `ragged` shortens each conformer by one residue,
    which is what "exact" must refuse and "intersection" must report.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif, n=8)
    others = []
    for index in range(1, n_conformers):
        other = tmp_path / f"conf{index}.cif"
        _fixture(other, n=8 - (index if ragged else 0), bend=index * 2.0)
        others.append({"path": str(other), "label": f"State {index}"})
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"), conformers=others,
        morph_mapping=mapping, morph_steps=steps, morph_reference_label="Start")
    html = next(item["path"] for item in result["artifacts"] if item["kind"] == "html")
    return result, Path(html)


def test_conformers_ship_endpoints_only_and_the_page_interpolates(tmp_path):
    """N frames in the file and NO stored intermediates; the page draws them."""
    result, html = _morph_bundle(tmp_path, n_conformers=3, steps=6)
    assert result["morph"]["mode"] == "browser"
    assert result["morph"]["frames_in_file"] == 3
    assert result["morph"]["stored_intermediates"] == 0
    state = json.loads(html.read_text().split('id="pinsp-inspection-state" type="application/json">')[1]
                       .split("</script>")[0].replace("<\\/", "</"))
    assert len(state["viewer"]["objects"][0]["frames"]) == 3

    report = _run_dom_harness(html, {"chains": ["A"] * 8, "residueNumbers": [1, 2, 3, 4, 5, 6, 7, 8]})
    assert report["morph"]["framesInFile"] == 3
    # One scratch frame is appended as the animation buffer -- and only one.
    assert report["morph"]["framesAfterExpansion"] == 4
    # FIRST TO LAST GOES STRAIGHT THERE. Every drawn frame is the buffer, and
    # the animation lands on the target conformation's own frame; it must never
    # pass through the middle conformation's frame on the way.
    assert report["morph"]["visitedCount"] > 3
    assert report["morph"]["buffersOnly"] is True
    assert report["morph"]["landedOn"] == report["morph"]["expectedLanding"] == 2
    assert report["morph"]["pressedAfter"] == ["false", "false", "true"]
    # The buffer is a genuine blend of the two endpoints it was asked for: the
    # harness seeds conformer k at y = 10k, so halfway from 0 to 2 reads y = 10
    # -- which is also what the middle conformer sits at, so the test checks
    # the buffer got there by interpolation, not by visiting frame 1.
    assert 0 < report["morph"]["bufferMaxY"] <= 20
    assert report["morph"]["bufferEndY"] == 20


def test_reference_label_and_button_per_conformation(tmp_path):
    result, html = _morph_bundle(tmp_path, n_conformers=3)
    assert [c["label"] for c in result["morph"]["conformers"]] == ["Start", "State 1", "State 2"]
    text = html.read_text()
    for index in range(3):
        assert f'data-conf="{index}"' in text
    assert 'class="bm-ring"' in text
    # A ring of segments is the thing this replaces the slider with, and the
    # frame transport (Play included) is gone with it.
    assert 'class="bm-btn"' in text
    assert "data-morph=\"1\"" in text
    assert ".pinsp-stage[data-morph='1'] #controlsContainer{display:none!important}" in text


def test_exact_mapping_refuses_a_ragged_conformer(tmp_path):
    with pytest.raises(ValueError, match="not exact"):
        _morph_bundle(tmp_path, n_conformers=2, mapping="exact", ragged=True)


def test_intersection_mapping_reports_what_each_file_lost(tmp_path):
    result, _ = _morph_bundle(tmp_path, n_conformers=3, mapping="intersection", ragged=True)
    morph = result["morph"]
    assert morph["mapping"] == "intersection"
    assert morph["n_residues"] == 6
    dropped = {c["label"]: c["residues_dropped"] for c in morph["conformers"]}
    assert dropped == {"Start": 2, "State 1": 1, "State 2": 0}


def test_conformers_and_morph_to_are_mutually_exclusive(tmp_path):
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    with pytest.raises(ValueError, match="not both"):
        render_inspection_bundle(
            mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
            inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
            output_dir=str(tmp_path / "out"),
            conformers=[{"path": str(cif)}], morph_to=[{"path": str(cif)}])


def test_the_viewer_cannot_paint_over_the_layer_panel(tmp_path):
    """py2Dmol fixes .py2dmol-viewer-instance at 948px; layout, not JS, must beat it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    # A flex COLUMN now, so the picture can take "the rest of the row" and the
    # conformation ring and capture bar keep their natural heights.
    assert (".pinsp-stage{position:relative;min-width:0;min-height:0;overflow:hidden;"
            "display:flex;flex-direction:column}") in text
    # Orient/Focus/Rotate/Style/Clip/Capture float over the top-LEFT of the
    # canvas instead of sitting in a 340px column beside it.
    assert ".pinsp-stage #rightPanelContainer{position:absolute!important;top:10px;left:10px;" in text
    # WIDTH FLOWS ONE WAY: grid column -> #mainContainer -> #canvasContainer.
    # py2Dmol's ResizeObserver answers a container resize by writing the
    # observed width back onto #viewerWrapper as an INLINE style, which beats
    # a stylesheet rule however important it is. Overriding that write was the
    # old approach and it lost: container and wrapper sized from each other,
    # each pass shedding the container's border, and the picture walked itself
    # narrower a couple of pixels per animation frame until it hit the floor.
    # The write is now inert by construction rather than by cascade -- the
    # element it lands on does not generate a box at all.
    assert (".pinsp-stage .py2dmol-viewer-instance,.pinsp-stage #viewerWrapper"
            "{display:contents!important}") in text
    # #mainContainer stays a real box: the floating Orient/Focus/Style cluster
    # is its child and a SIBLING of #canvasContainer, so it is the positioning
    # context those tools resolve against. Collapse it too and they detach and
    # land on whatever sits above the picture.
    assert (".pinsp-stage #mainContainer{display:block!important;position:relative;"
            "width:auto!important;max-width:none!important;padding:0!important;"
            "flex:1 1 auto;min-height:0}") in text
    # setupViewport writes an INLINE pixel width on #canvasContainer, so the
    # override has to be !important or the picture stays a fixed box.
    assert ".pinsp-stage #canvasContainer{display:block!important;width:auto!important;" in text
    assert "resize:none!important" in text
    # NOTHING MEASURES ANYTHING. The JS width controller is gone; if it ever
    # comes back, so does the feedback loop.
    assert "fitViewer" not in text
    assert "box.style.setProperty('width'" not in text
    assert "new ResizeObserver(fitViewer)" not in text
    # The one write CSS cannot do: the viewer root carries an id and no class
    # in an export, so it is reached by walking up from #mainContainer -- once,
    # at startup, with no measurement.
    assert "instance.style.setProperty('display','contents','important')" in text
    assert ".pinsp-stage #canvasContainer canvas{max-width:100%}" in text
    # The Layers panel keeps its own column beside the structure, and there is
    # no collapse button to hide the viewer's controls with.
    assert "grid-template-columns:minmax(240px,1fr) 11px var(--pinsp-panel,340px)" in text
    assert "@media (max-width:820px)" in text
    assert 'id="pinsp-controls"' not in text
    assert "Hide controls" not in text


def test_a_blocked_download_still_leaves_the_capture_recoverable(tmp_path):
    """Save Image reports success without checking; an embedded frame may drop it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert 'id="pinsp-capture"' in text
    # The blob is kept and the revoke deferred, or there is nothing to offer.
    assert "URL.createObjectURL=function(blob){lastBlob=blob" in text
    assert "URL.revokeObjectURL=function(url){setTimeout(" in text
    assert "blocks downloads" in text


def _partner_bundle(tmp_path, morph=False):
    """One file: chain A is the target, chain P is a bound partner."""
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    for chain_id, shift in (("A", 0.0), ("P", 12.0)):
        chain = Chain.Chain(chain_id)
        for index in range(1, 5):
            residue = Residue.Residue((" ", index, " "), "ALA", " ")
            residue.add(Atom.Atom("CA", (float(index * 3), float(index % 2) + shift, 0.0),
                                  0.0, 1.0, " ", "CA", index, element="C"))
            chain.add(residue)
        model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    cif = tmp_path / "complex.cif"
    io.save(str(cif))

    conformers = None
    if morph:
        other = tmp_path / "bent.cif"
        _fixture(other, n=4, bend=1.5)
        conformers = [{"path": str(other), "label": "Bent"}]
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"), partner_chains=["P"], partner_label="Fab",
        conformers=conformers, morph_mapping="exact", morph_steps=6)
    html = next(item["path"] for item in result["artifacts"] if item["kind"] == "html")
    return result, Path(html)


def test_partner_chains_get_one_button_that_hides_exactly_them(tmp_path):
    _, html = _partner_bundle(tmp_path)
    text = html.read_text()
    # ONE BUTTON, not a checkbox.
    assert 'id="pinsp-partners"' in text
    assert 'class="bm-btn bm-solo"' in text
    # ...and it is not the old checkbox. (The py2Dmol library itself contains
    # the word, so the assertion has to name this control.)
    assert 'type="checkbox" id="pinsp-partners"' not in text
    report = _run_dom_harness(html, {"chains": ["A", "A", "A", "A", "P", "P", "P", "P"],
                                     "residueNumbers": [1, 2, 3, 4, 1, 2, 3, 4]})
    partners = report["partners"]
    assert partners["disabledAtStart"] is False
    assert partners["pressedAtStart"] == "true"
    assert partners["onAtStart"] == partners["expectedOn"] == [0, 1, 2, 3, 4, 5, 6, 7]
    assert partners["afterHide"] == partners["expectedOff"] == [0, 1, 2, 3]
    assert partners["pressedAfterHide"] == "false"
    assert partners["afterShow"] == [0, 1, 2, 3, 4, 5, 6, 7]


def test_a_partner_is_excluded_from_the_morph_and_held_still(tmp_path):
    """The partner exists in the reference only; interpolating it is meaningless."""
    result, _ = _partner_bundle(tmp_path, morph=True)
    morph = result["morph"]
    # The residue mapping is the target chain alone -- the conformer has no P.
    assert morph["n_residues"] == 4
    assert morph["mapping"] == "exact"
    # ...but the partner is still drawn, as its own block of positions.
    assert morph["partner_blocks"] == [
        {"conformer": 0, "label": "reference", "chains": ["P"],
         "renamed_from": None, "n_positions": 4}]


def test_every_export_retains_the_upstream_licence_and_credit(tmp_path):
    """BEER-WARE asks one thing: retain the notice. Minification removes it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert "THE BEER-WARE LICENSE" in text
    assert "Sergey Ovchinnikov" in text
    assert "https://github.com/sokrypton/py2Dmol" in text
    assert "78c2d489d0b5c5d19accd9eeeef878c2868f5271" in text
    # Not only in a comment: a reader of the page can see who wrote the viewer.
    assert 'class="pinsp-credit"' in text
    assert "py2Dmol</a> by Sergey Ovchinnikov" in text
    # BOTH notices, because the additions are on the same terms -- a reader
    # allowed to reuse py2Dmol should not have to work out where it ends.
    assert "Tadas Kluonis" in text
    assert "https://github.com/tadaskluonis/protein_inspector" in text
    assert "free to reuse" in text
    # ONE footer, on the bottom edge of the sequence section. It used to be
    # emitted twice over -- once inside the structure pane and once after
    # `panes_close` -- and it is inside the section now so that zone reaches
    # the foot of the page instead of floating above a band of nothing.
    assert text.count('class="pinsp-credit"') == 1
    seq = text.split('<section class="pinsp-seq"')[1].split("</section>")[0]
    assert 'class="pinsp-credit"' in seq
    assert seq.index('id="pinsp-seqbody"') < seq.index('class="pinsp-credit"')


def _two_partner_bundle(tmp_path):
    """Reference carries partner P; the second conformation carries its own, Q."""
    def complex_cif(path, partner_chain, shift, bend=0.0):
        structure = Structure.Structure("fixture")
        model = Model.Model(0)
        for chain_id, offset, bendy in (("A", 0.0, bend), (partner_chain, shift, 0.0)):
            chain = Chain.Chain(chain_id)
            for index in range(1, 5):
                residue = Residue.Residue((" ", index, " "), "ALA", " ")
                residue.add(Atom.Atom("CA", (float(index * 3),
                                             float(index % 2) + offset + bendy * index, 0.0),
                                      0.0, 1.0, " ", "CA", index, element="C"))
                chain.add(residue)
            model.add(chain)
        structure.add(model)
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(path))
        return path

    reference = complex_cif(tmp_path / "ref.cif", "P", 12.0)
    other = complex_cif(tmp_path / "alt.cif", "Q", -12.0, bend=1.5)
    result = render_inspection_bundle(
        mmcif_path=str(reference),
        mmcif_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"),
        partner_chains=["P"], partner_label="Bound partner",
        conformers=[{"path": str(other), "label": "Alt", "partner_chains": ["Q"]}],
        morph_mapping="exact", morph_steps=6)
    html = next(item["path"] for item in result["artifacts"] if item["kind"] == "html")
    return result, Path(html)


def test_each_conformation_shows_its_own_partner(tmp_path):
    """The tethered receptor must not keep wearing the extended form's Fab."""
    result, html = _two_partner_bundle(tmp_path)
    blocks = result["morph"]["partner_blocks"]
    assert [b["conformer"] for b in blocks] == [0, 1]
    assert [b["chains"] for b in blocks] == [["P"], ["Q"]]
    assert [b["n_positions"] for b in blocks] == [4, 4]
    # Target mapping is chain A alone in both files.
    assert result["morph"]["n_residues"] == 4

    report = _run_dom_harness(html, {
        "chains": ["A"] * 4 + ["P"] * 4 + ["Q"] * 4,
        "residueNumbers": [1, 2, 3, 4] * 3})
    partners = report["partners"]
    # State 0: target + P, no Q.
    assert partners["onAtStart"] == partners["expectedOn"] == [0, 1, 2, 3, 4, 5, 6, 7]
    # After morphing to state 1: target + Q, and P is gone.
    assert partners["afterMorph"] == partners["expectedAfterMorph"] \
        == [0, 1, 2, 3, 8, 9, 10, 11]
    assert partners["disabledAfterMorph"] is False


def test_a_colliding_partner_chain_is_renamed_not_merged(tmp_path):
    """Two conformations both carrying chain B are two molecules, not one."""
    def complex_cif(path, shift, bend=0.0):
        structure = Structure.Structure("fixture")
        model = Model.Model(0)
        for chain_id, offset, bendy in (("A", 0.0, bend), ("B", shift, 0.0)):
            chain = Chain.Chain(chain_id)
            for index in range(1, 5):
                residue = Residue.Residue((" ", index, " "), "ALA", " ")
                residue.add(Atom.Atom("CA", (float(index * 3),
                                             float(index % 2) + offset + bendy * index, 0.0),
                                      0.0, 1.0, " ", "CA", index, element="C"))
                chain.add(residue)
            model.add(chain)
        structure.add(model)
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(path))
        return path

    reference = complex_cif(tmp_path / "ref.cif", 12.0)
    other = complex_cif(tmp_path / "alt.cif", -12.0, bend=1.5)
    result = render_inspection_bundle(
        mmcif_path=str(reference),
        mmcif_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"),
        partner_chains=["B"],
        conformers=[{"path": str(other), "label": "Alt", "partner_chains": ["B"]}],
        morph_mapping="exact")
    blocks = result["morph"]["partner_blocks"]
    assert [b["chains"] for b in blocks] == [["B"], ["B1"]]
    assert blocks[1]["renamed_from"] == {"B1": "B"}


def test_a_conformers_partner_travels_with_its_receptor(tmp_path):
    """The fit applied to the target must be applied to its partner too.

    Otherwise the receptor rotates onto the reference and its partner stays
    where its own crystal put it -- floating in space beside nothing.
    """
    def complex_cif(path, rotate_by, partner_chain):
        angle = np.deg2rad(rotate_by)
        rot = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                        [np.sin(angle), np.cos(angle), 0.0],
                        [0.0, 0.0, 1.0]])
        structure = Structure.Structure("fixture")
        model = Model.Model(0)
        for chain_id, offset in (("A", 0.0), (partner_chain, 9.0)):
            chain = Chain.Chain(chain_id)
            for index in range(1, 7):
                xyz = np.array([float(index * 3), float(index % 2) + offset, 0.0]) @ rot.T
                residue = Residue.Residue((" ", index, " "), "ALA", " ")
                residue.add(Atom.Atom("CA", tuple(xyz), 0.0, 1.0, " ", "CA", index, element="C"))
                chain.add(residue)
            model.add(chain)
        structure.add(model)
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(path))
        return path

    # Same complex, one rotated 70 degrees in the plane. Superposing the
    # targets must bring the partners on top of each other too.
    reference = complex_cif(tmp_path / "ref.cif", 0.0, "P")
    turned = complex_cif(tmp_path / "turned.cif", 70.0, "Q")
    result = render_inspection_bundle(
        mmcif_path=str(reference),
        mmcif_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"), partner_chains=["P"],
        conformers=[{"path": str(turned), "label": "Turned", "partner_chains": ["Q"]}],
        morph_mapping="exact", extras=True)

    state = json.loads(Path(next(item["path"] for item in result["artifacts"]
                                 if item["kind"] == "viewer_state")).read_text())
    frame = state["viewer"]["objects"][0]["frames"][0]
    chains = frame["chains"]
    coords = np.array(frame["coords"], dtype=float)
    p_block = coords[[i for i, c in enumerate(chains) if c == "P"]]
    q_block = coords[[i for i, c in enumerate(chains) if c == "Q"]]
    assert len(p_block) == len(q_block) == 6
    # The turned copy's partner lands on the reference's partner, not 70
    # degrees away from it.
    assert np.abs(p_block - q_block).max() < 1e-6


def test_the_partner_button_does_not_move_the_conformation(tmp_path):
    """Showing or hiding a partner is not a request to change state.

    py2Dmol's own setVisibility drops the viewer back to the first frame, so
    the toggle became an unasked-for jump to the reference. The page has to
    put the frame back; the harness reproduces the drop so it is exercised.
    """
    _, html = _two_partner_bundle(tmp_path)
    report = _run_dom_harness(html, {
        "chains": ["A"] * 4 + ["P"] * 4 + ["Q"] * 4,
        "residueNumbers": [1, 2, 3, 4] * 3})
    partners = report["partners"]
    expected = partners["expectedFrameOnState"]
    assert partners["frameAfterHideOnState"] == expected
    assert partners["frameAfterShowOnState"] == expected
    # ...and on the first state too, where the bug is invisible.
    assert partners["frameAfterToggle"] == 0


def test_the_two_panes_share_one_draggable_boundary(tmp_path):
    """Sized independently, the picture can always overlap the panel or leave
    a gap. One seam, owned by the grid, removes both states."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert 'id="pinsp-split"' in text
    assert 'role="separator"' in text
    # gap:0 -- the panes touch; the seam IS the gutter. The columns and the
    # rows are separate declarations now that the sequence strip is a third
    # zone, so this asserts the two facts rather than one byte string that
    # happened to carry both.
    assert "grid-template-columns:minmax(240px,1fr) 11px var(--pinsp-panel,340px);" in text
    assert "gap:0;" in text
    # The picture is no longer independently resizable.
    assert "resize:none!important" in text
    assert ".pinsp-stage #canvasContainer .resize-handle{display:none!important}" in text

    report = _run_dom_harness(html, {"chains": ["A"] * 5,
                                     "residueNumbers": [1, 2, 3, 4, 5]})
    split = report["split"]
    assert split["dragAttr"] == "1"
    # Grid spans 0..1000, so a pointer at 600 asks for a 400px panel.
    assert split["at600"] == 400
    # Clamped: never wider than width - 320 stage, never below 240.
    assert split["clampedWide"] == 680
    assert split["clampedNarrow"] == 240
    assert split["dragAttrAfter"] is None
    # Released, the seam stops following the pointer.
    assert split["afterRelease"] == 240
    # Keyboard moves it too -- left widens the panel.
    assert split["afterArrowLeft"] == 256
    assert split["afterHome"] == 340


def _atom_fixture(path: Path, n: int = 12) -> None:
    """A short chain with real atoms, not a Cα trace.

    `_fixture` writes one CA per residue, which is the right fixture for the
    layer and morph tests and the wrong one here: a trace carries no
    side-chain atoms, so `view(sidechains=True)` has nothing to capture and
    every assertion about the Side chains button would pass against a bundle
    that cannot draw one. CB is what makes this a side chain; ALA and LEU
    alternate so the strip's letters are distinguishable.
    """
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for index in range(1, n + 1):
        residue = Residue.Residue((" ", index, " "), "ALA" if index % 3 else "LEU", " ")
        x = float(index * 3)
        y = float(index % 2)
        for name, offset, element in (("N", (-1.2, 0.0, 0.0), "N"),
                                      ("CA", (0.0, 0.0, 0.0), "C"),
                                      ("C", (1.2, 0.0, 0.0), "C"),
                                      ("O", (1.2, 1.1, 0.0), "O"),
                                      ("CB", (0.0, 1.5, 0.6), "C")):
            residue.add(Atom.Atom(name, (x + offset[0], y + offset[1], offset[2]),
                                  0.0, 1.0, " ", name, index * 5, element=element))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))


_SEQUENCE_SPEC = {
    "chains": ["A"] * 12,
    "residueNumbers": list(range(1, 13)),
    "positionNames": ["ALA", "ALA", "LEU"] * 4,
}


def _sequence_bundle(tmp_path, sidechains=True):
    """A 12-residue bundle with two layers, for the strip and the tools."""
    cif = tmp_path / "atoms.cif"
    _atom_fixture(cif)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "a2", "kind": "custom", "label": "Layer A res 2",
             "layer_id": "layer-a", "layer_label": "Layer A", "color": "#111111",
             "resolved": True, "evidence_ids": [], "method": "rule",
             "residue": {"component_id": "ALA", "canonical_position": 2,
                         "chain_id": "A", "author_residue_number": 2}},
            {"annotation_id": "b5", "kind": "custom", "label": "Layer B res 5",
             "layer_id": "layer-b", "layer_label": "Layer B", "color": "#dc2626",
             "resolved": True, "evidence_ids": [], "method": "rule",
             "residue": {"component_id": "ALA", "canonical_position": 5,
                         "chain_id": "A", "author_residue_number": 5}},
        ],
    }
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest=manifest, output_dir=str(tmp_path / "out"), extras=True,
        display_options={"width": 600, "height": 400, "sidechains": sidechains})
    return Path(next(i["path"] for i in result["artifacts"] if i["kind"] == "html"))


def test_the_sequence_strip_is_one_letter_per_residue_in_the_layer_colours(tmp_path):
    """The strip is the Layers legend with an index: same colours, same
    visible-layer set, one cell per drawn position in file order.

    Built from the renderer's own arrays rather than from the manifest, so a
    position the manifest never mentions still has a letter -- which is the
    whole point, since most of a structure is unannotated.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    seq = report["seq"]
    assert seq["cells"] == 12
    assert seq["letters"] == "AALAALAALAAL"
    # The author number is printed over every tenth residue and nowhere else.
    assert seq["ticks"] == ["10"]
    assert seq["chainLabels"] == ["A"]
    assert seq["firstTitle"] == "A \u00b7 ALA 1"
    # Only the two annotated positions are painted; the rest are left to the
    # stylesheet, for the same reason the structure's base is.
    assert seq["colours"][1] == "#111111"
    assert seq["colours"][4] == "#dc2626"
    assert [c for c in seq["colours"] if c] == ["#111111", "#dc2626"]


def test_clicking_the_viewer_selects_a_residue_and_opens_its_card(tmp_path):
    """py2Dmol ships with canvas picking OFF -- parts/ui.js turns it on for
    Focus mode alone -- so in an exported page the structure was inert.

    Turning it on is one flag; the card is the other half. It could only open
    on a residue that carried an annotation, and a reader clicking the
    structure is usually asking about one that does not.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    assert report["seq"]["selectionEnabled"] is True
    click = report["viewerClick"]
    # Position 7 carries no annotation. The card still names it.
    assert "ALA 8" in click["card"]
    assert "chain A &middot; residue 8" in click["card"]
    assert "No annotation on this residue." in click["card"]
    assert click["details"] == "ALA 8"
    assert click["marked"] == [7]
    assert click["readout"].startswith("A 8")
    # A double-click takes the whole chain, which is not one residue: the
    # readout says which, and the card stands down rather than picking one.
    assert click["chainCard"] is True
    assert click["chainReadout"].startswith("A 1\u20133")
    # Clicking the background clears every surface at once.
    assert click["afterBackground"]["cardHidden"] is True
    assert click["afterBackground"]["marked"] == []
    assert click["afterBackground"]["details"].startswith("Click a residue")


def test_the_strip_adds_and_subtracts_and_never_replaces(tmp_path):
    """A click used to throw the selection away and start again, which makes
    the strip useless for what it is for: building a set a few residues at a
    time, across two chains, on top of what the layers already gave you.

    What the residue under the pointer ALREADY IS decides the gesture -- start
    on an unselected letter and the drag adds, start on a selected one and it
    takes away. One rule, both directions, no modifier to remember.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    seq = report["seq"]
    assert seq["afterClick"]["selection"] == [2]
    assert seq["afterClick"]["marked"] == [2]
    assert seq["afterClick"]["readout"].startswith("A 3")
    assert "LEU 3" in seq["afterClick"]["card"]
    # Dragging from the anchor to another letter takes everything between.
    assert seq["afterDrag"] == [2, 3, 4, 5]
    # ...and a release ends it, so the next move is hover and not selection.
    assert seq["afterRelease"] == [2, 3, 4, 5]
    # An unselected letter JOINS what is already there.
    assert seq["afterAddClick"] == [2, 3, 4, 5, 8]
    # ...and pressing it again takes only it away.
    assert seq["afterRemoveClick"] == [2, 3, 4, 5]
    # A press on a SELECTED letter starts a subtracting drag.
    assert seq["beforeSubtract"] == [2, 4, 5]
    assert seq["afterSubtractDrag"] == [2]
    # ...and dragging back over your own path shrinks the range rather than
    # ratcheting it: every move recomputes from the mousedown snapshot, so 5
    # comes back the moment the range no longer covers it.
    assert seq["afterDragBack"] == [2, 5]
    # Shift extends from the last letter pressed, additively.
    assert seq["afterShift"] == [1, 2, 3, 4, 5]
    # The chain button reads the same way: all in -> out, otherwise in.
    assert seq["afterChainLabel"] == 12
    assert seq["afterChainLabelAgain"] == 0


def test_select_around_asks_the_renderers_own_neighbourhood_search(tmp_path):
    """A selection made FROM a selection, which is the other half of picking
    one: "what lines this pocket" is a question about what is near what you
    already have.

    The page must not measure anything itself -- `residuesWithin` is the
    renderer's own atom-to-atom search on a grid it keeps between calls -- so
    what is asserted here is that the distance box reaches it and whatever it
    answers becomes the selection, seed included, as PyMOL's byres does.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    around = report["around"]
    assert around["disabledWithSelection"] is False
    assert around["calls"] == [{"seed": [4], "cutoff": 2}]
    assert around["selection"] == [2, 3, 4, 5, 6]
    assert around["readout"].startswith("A 3\u20137")
    # A tool that acts on the selection is not pressable without one.
    assert around["clearedSelection"] == []
    assert around["disabledWhenEmpty"] is True


def test_the_side_chain_button_draws_exactly_the_selection(tmp_path):
    """Clicking a residue and asking to SEE it is the third thing the reader
    could not do: the verb exists in the bundle (parts/sidechains.js) and
    nothing in an exported page reached it.

    It is a latch, and its state is read back off the renderer's own set
    rather than kept here, so Focus mode cannot turn side chains on behind the
    button's back and leave it saying off.
    """
    html = _sequence_bundle(tmp_path)
    state = json.loads(Path(str(html).replace(".html", ".viewer.json")).read_text())
    # The atoms are in the file -- this is what the button's enabled state and
    # `showSidechains` both depend on.
    assert state["sidechains"] is True
    frame = state["viewer"]["objects"][0]["frames"][0]
    assert len(frame["sidechain_atoms"]) == 12

    report = _run_dom_harness(html, _SEQUENCE_SPEC)
    sc = report["sidechains"]
    assert sc["stateSaysPresent"] is True
    assert sc["disabled"] is False
    # A PAIR, AND THREE STATES. `Show` fills when every selected residue
    # draws, `Hide` when none does, NEITHER when they disagree -- the state a
    # single latch could only render as "off", which is what made it say what
    # it would do rather than what it had done (src/app/selection.js).
    assert sc["pressed"] == ["false", "true"]                 # none drawn yet
    assert sc["afterShow"]["calls"] == [{"on": True, "positions": [6]}]
    assert sc["afterShow"]["pressed"] == ["true", "false"]
    assert sc["afterHide"]["calls"][-1] == {"on": False, "positions": [6]}
    assert sc["afterHide"]["pressed"] == ["false", "true"]
    # One residue drawing, one not: neither button is filled.
    assert sc["mixedSelection"] == [6, 7]
    assert sc["mixedPressed"] == ["false", "false"]
    # ...and the buttons are ABSOLUTE, so one press resolves a disagreement
    # instead of flipping it to whichever state the majority was not in.
    assert sc["afterResolve"]["calls"] == [{"on": True, "positions": [6, 7]}]
    assert sc["afterResolve"]["pressed"] == ["true", "false"]
    # THE PER-RESIDUE STATE OUTLIVES THE SELECTION. Dropping the selection
    # undraws nothing, and the bar says so rather than going blank.
    assert sc["afterDeselect"]["shown"] == 2
    assert "2 residues still drawing side chains" in sc["afterDeselect"]["readout"]
    assert sc["afterDeselect"]["disabled"] is True
    assert sc["afterDeselect"]["pressed"] == ["false", "false"]


def test_without_side_chain_atoms_the_button_explains_itself(tmp_path):
    """`showSidechains` RAISES on a structure that carries none, and a Cα
    trace and every morph frame carry none. A control that throws is worse
    than one that is disabled and says why."""
    html = _sequence_bundle(tmp_path, sidechains=False)
    state = json.loads(Path(str(html).replace(".html", ".viewer.json")).read_text())
    assert state["sidechains"] is False
    assert "sidechain_atoms" not in state["viewer"]["objects"][0]["frames"][0]
    report = _run_dom_harness(html, _SEQUENCE_SPEC)
    assert report["sidechains"]["disabled"] is True
    # Pressed and nothing happened, rather than an exception in the console.
    assert report["sidechains"]["afterShow"]["calls"] == []
    assert "no side-chain atoms" in Path(html).read_text()


def test_hovering_the_structure_lights_the_letter(tmp_path):
    """core/mol.js already calls window.SEQ.setHoveredResidue on mousemove
    when a strip is present, and skips the pick entirely when it is not. The
    bridge is installed only once the strip exists, so the per-move cost is
    never paid by a page with nothing to show for it."""
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    assert report["hover"]["lit"] == [5]
    assert report["hover"]["afterLeave"] == 0


def test_the_sequence_seam_trades_height_between_the_picture_and_the_strip(tmp_path):
    """One seam, two rows, constant sum -- which is what the vertical seam has
    always been, on the other axis.

    Two faults, reported together and with the same cause. The gesture
    measured against the grid's own bottom edge, which MOVES when the strip
    resizes, so every pointermove read back its own previous result and the
    seam ran away from the pointer to its ceiling in a few events. And only
    the strip's row was written, so the page simply got taller: the one thing
    a reader could see move was the credit line underneath it.
    """
    html = _sequence_bundle(tmp_path)          # display height 400
    text = Path(html).read_text()
    assert 'id="pinsp-hsplit"' in text
    assert 'aria-orientation="horizontal"' in text
    # ROW ONE IS A HEIGHT, NOT `auto` -- see the test below for why.
    assert "grid-template-rows:var(--pinsp-view,400px) 11px var(--pinsp-seq,152px)" in text
    # Every zone names its own cell: with rows that are not all `auto`,
    # auto-placement put the Layers panel in the 11px seam row on narrow
    # screens, where the vertical seam is display:none and placed nowhere.
    assert ".pinsp-seq{grid-column:1/-1;grid-row:3}" in text
    # The picture's height is the property the seam writes, seeded from the
    # caller's `display` height so nothing moves until somebody drags.
    assert ".protein-inspector{--pinsp-view:400px;--pinsp-seq:152px}" in text
    # The picture's height comes from the STAGE, which the row sizes -- not
    # from the property directly, or the conformation ring and the capture bar
    # would be pushed out of a row that is now a fixed height.
    assert "height:100%!important" in text
    # The gesture is pointer-travel-since-pointerdown and nothing else. Said
    # positively: the comment above it quotes the expression it replaced, so a
    # "the old one is gone" assertion matches the explanation and passes
    # whatever the code does.
    assert "anchorSeq+(anchorY-e.clientY)" in text

    report = _run_dom_harness(html, _SEQUENCE_SPEC)
    hsplit = report["hsplit"]
    # THE SUM IS THE WINDOW'S, NOT THE CALLER'S. Both rows used to be fixed
    # pixel heights, so the page was exactly as tall as `display` asked for and
    # any taller window left a band of blank page under the sequence section --
    # the section stopped short of the bottom. The strip keeps the height it was
    # given and the picture takes the slack: the harness window is 900 px tall
    # with the grid at y=0, less the 14 px the page pads the body by and the
    # 11 px seam, so the pair sums to 875 and the view is 875 - 152.
    assert hsplit["start"] == {"seq": 152, "view": 723}
    assert hsplit["dragAttr"] == "1"
    total = hsplit["start"]["seq"] + hsplit["start"]["view"]
    assert total == 875
    # A pointermove at the anchor moves nothing. Under the old expression this
    # was the step that drifted, because the measurement had already changed.
    assert hsplit["atAnchor"] == {"seq": 152, "view": 723}
    # 40 px up: the strip takes 40, the picture gives up 40, the sum holds.
    assert hsplit["up40"] == {"seq": 192, "view": 683}
    # ...and the gesture is anchored, so coming back to the start comes back
    # to the start rather than to somewhere further along.
    assert hsplit["backToAnchor"] == {"seq": 152, "view": 723}
    assert hsplit["clampedShort"] == {"seq": 64, "view": total - 64}
    assert hsplit["clampedTall"] == {"seq": total - 200, "view": 200}
    assert hsplit["dragAttrAfter"] is None
    assert hsplit["afterRelease"] == hsplit["clampedTall"]   # stops following
    assert hsplit["afterHome"] == {"seq": 152, "view": 723}
    assert hsplit["afterArrowUp"] == {"seq": 168, "view": 707}
    assert hsplit["afterArrowDown"] == {"seq": 152, "view": 723}
    # Every step keeps the pair summing to the same total, which is what keeps
    # the page's height -- and the credit line -- still.
    for step in ("atAnchor", "up40", "backToAnchor", "clampedShort",
                 "clampedTall", "afterHome", "afterArrowUp", "afterArrowDown"):
        assert hsplit[step]["seq"] + hsplit[step]["view"] == total, step


def test_the_layers_panel_cannot_hold_the_seam_open(tmp_path):
    """The seam moved nothing while the picture resized inside it.

    Row one was `auto` and the Layers panel lives in it. On a bundle with a
    long residue list the panel is the tallest thing in that row -- it was
    capped at `max-height:88vh` and filled it -- so the row's height was the
    PANEL's, not the picture's. Shrinking the canvas then left the seam
    exactly where it was: the reader grabbed it, saw the structure resize, and
    saw the thing in their hand stay still.

    ASSERTED AS CSS, not as behaviour, and deliberately. A seam is a box-model
    fact: the DOM harness has no layout, so it reported this one working --
    the page wrote both custom properties and the harness read them back. What
    actually decides whether the seam moves is which element's height the grid
    row takes, and these five declarations are that decision. Run the page in
    a browser to watch it; run this to keep it.
    """
    text = Path(_sequence_bundle(tmp_path)).read_text()
    # The row is a height the seam owns, seeded from the caller's `display`.
    assert "grid-template-rows:var(--pinsp-view,400px) 11px var(--pinsp-seq,152px)" in text
    # The panel takes the row and scrolls inside it, rather than setting it.
    assert ".pinsp-panel{border:1px solid #e2e8f0;border-radius:10px;min-height:0;overflow:auto;" in text
    assert "max-height:88vh" not in text
    assert "position:sticky;top:14px" not in text
    # The picture takes what is left of the stage, so the conformation ring and
    # the capture bar are not pushed out of a row that is now a fixed height.
    assert ".pinsp-stage{position:relative;min-width:0;min-height:0;overflow:hidden;" in text
    assert "display:flex;flex-direction:column}" in text
    assert "height:100%!important;position:relative}" in text
    # ...and #mainContainer is still a positioned box, because the floating
    # Orient/Focus/Style cluster is its child and a sibling of the canvas.
    assert "flex:1 1 auto;min-height:0}" in text
    assert ".pinsp-stage #mainContainer{display:block!important;position:relative;" in text


def test_the_selection_mark_is_drawn_on_the_click_that_made_it(tmp_path):
    """core/mol.js's mouseup sets the selection and renders only when the
    molecule is large, on the reasoning that a small one's next frame is along
    shortly. True on the website, where a hover readout and a sequence canvas
    are repainting anyway; false in a still bundle, where the yellow outline
    appeared on the first frame AFTER the click -- nudge the structure and
    there it was."""
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    assert report["viewerClick"]["rendered"] >= 1
    assert report["viewerClick"]["marked"] == [7]


def test_the_cartoon_is_painted_on_the_gpu_by_default(tmp_path):
    """A `richardson` or `cartoon` bundle felt slow next to py2dmol.solab.org
    on the same machine because this passed `gpu=False` while the website's
    default is True -- the same geometry, handed to the CPU painter.

    The notebook bundle is the one build carrying BOTH painters, and
    core/mol.js only honours the flag when both are present, so nothing is
    given up: no WebGL2 or a lost context and `_gpuWillTake` declines the
    frame, and paintgl refuses any export context, so SVG still comes off the
    2D painter.
    """
    html = _sequence_bundle(tmp_path)
    state = json.loads(Path(str(html).replace(".html", ".viewer.json")).read_text())
    assert state["viewer"]["config"]["rendering"]["gpu"] is True

    cif = tmp_path / "atoms.cif"
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "cpu"), extras=True,
        display_options={"gpu": False})
    cpu = json.loads(Path(next(i["path"] for i in result["artifacts"]
                               if i["kind"] == "viewer_state")).read_text())
    assert cpu["viewer"]["config"]["rendering"]["gpu"] is False


def test_the_cartoon_is_sampled_for_a_reader_who_zooms_in(tmp_path):
    """py2Dmol samples the cartoon at `detail` subdivisions per residue and
    nothing else: upstream retired an adaptive term that targeted a fixed
    on-screen chord length, so a magnified curve now facets rather than
    resampling -- cartoon/geom.js says so in those words.

    At the library default of 4 that is 4 subdivisions per helix residue and 3
    per strand or loop (HELIX_SUB 8, SHEET_SUB and SUB both 6, scaled by
    detail/8). Zoomed to where a residue fills the canvas, three flat plates
    per residue of loop -- each carrying its own outline stroke -- was read by
    the first reader to zoom a bundle in as malformed side-chain bonds.

    A bundle is made to be zoomed into, so it asks for the top of the range.
    The Detail slider still spans 2-8 live and an explicit value still wins.
    """
    html = _sequence_bundle(tmp_path)
    state = json.loads(Path(str(html).replace(".html", ".viewer.json")).read_text())
    assert state["viewer"]["config"]["rendering"]["detail"] == 8

    cif = tmp_path / "atoms.cif"
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "coarse"), extras=True,
        display_options={"detail": 2})
    coarse = json.loads(Path(next(i["path"] for i in result["artifacts"]
                                  if i["kind"] == "viewer_state")).read_text())
    assert coarse["viewer"]["config"]["rendering"]["detail"] == 2


def test_a_conformer_frame_says_what_its_positions_are(tmp_path):
    """A frame that does not carry `position_types` cannot carry side chains.

    `_materialiseSidechains` appends one position per side-chain atom and types
    each `'L'`, building the array as `(data.position_types || []).slice()` plus
    one push per atom. With no types in the frame that starts EMPTY, so the
    finished array is as long as the appended atoms instead of as long as the
    coordinates -- and `_setDataField` takes a per-position array whose length
    does not match the coordinate count and silently replaces it with
    `Array(n).fill('P')`. Every side-chain atom then claimed to be a protein
    backbone position, and `cartoon/geom.js` partitions on exactly that
    (backbone drawn as cartoon, everything else as sticks), so the atoms were
    swept into the ribbon and drawn as slabs running between them. The same
    structure through `add_pdb` was right, because `add_pdb` sends types.

    Reported as "the boxes connecting the atoms look wrong in between", and
    invisible to everything that had been checked: the side-chain table is
    chemically exact and its geometry survives the Kabsch fit intact (CA-CB
    median 1.530 A in both frames of the EGFR bundle). The atoms were never
    wrong. They were mislabelled.

    The type is derived per residue rather than filled with 'P', because
    `_ca_trace` accepts C4' as well as CA: a nucleic position typed 'P' is
    rebuilt through the peptide's step range, `localFrame` fails, and the base's
    atoms are dropped without a word.
    """
    _result, html = _morph_bundle(tmp_path)
    state = json.loads(html.read_text()
                       .split('id="pinsp-inspection-state" type="application/json">')[1]
                       .split("</script>")[0].replace("<\\/", "</"))
    frames = state["viewer"]["objects"][0]["frames"]
    assert len(frames) >= 2
    for frame in frames:
        types = frame.get("position_types")
        assert types is not None, "a conformer frame must declare its types"
        assert len(types) == len(frame["coords"])
        assert set(types) <= {"P", "D", "R", "L"}

    assert _position_types(["ALA", "DA", "A", "HOH", "", "dg"]) == [
        "P", "D", "R", "P", "P", "D"]


def test_selecting_never_draws_a_side_chain(tmp_path):
    """Nothing on the page draws a side chain except the Show/Hide pair -- and
    py2Dmol's Focus button, which is in the bundle's own control bar, was the
    exception that made that untrue.

    Focus is a MODE that redefines what a selection means. `enterFocusMode`
    clears every object's `sidechains` on the way in, and `focusOn` then
    ASSIGNS the set from the neighbourhood of whatever was picked -- writing
    every object on screen, because "hide the last click's" is the point of it.
    It is wrapped around `setResidueSelection` rather than subscribed to it, so
    with the mode latched EVERY selection drew atoms: a click in the picture, a
    letter in the strip, `sel` on a layer, `Select around`. Reported as
    selecting a residue showing its side chains before anyone asked.

    Upstream already has the opt-out -- `focusOn`'s options carry `sidechains`
    and `o.sidechains !== false` is the only gate on that block -- so the
    camera half of the mode is kept and the drawn set is left to the one
    control that is about the drawn set.
    """
    text = Path(_sequence_bundle(tmp_path)).read_text()
    assert "o.sidechains=false" in text
    assert "r._pinspFocusTamed" in text
    # The wrap is on the renderer's own verbs, not a re-implementation of the
    # mode: a bundle should say one thing about Focus, not own a copy of it.
    assert "focusOn.call(r,sel,o)" in text
    assert "enter.call(r)" in text


def test_detail_survives_the_style_switch_it_was_rendered_for(tmp_path):
    """`display={"detail": N}` seeds the viewer's config, and the constructor
    honours it -- for as long as the viewer stays in the style it opened in.

    core/mol.js's `_applyLookDefaults` re-asserts every style-owned control
    from `LOOK_DEFAULTS`, and all three entries there carry `detail: 4`, with
    none of the "did a person choose this" latching that width and thickness
    get. A bundle opens as `tube`, so the first time a reader picks Richardson
    the sampling dropped from 8 back to 4 -- and Richardson is the style they
    picked it for. The bundle shipped a detail nobody could ever see, which is
    worse than shipping the default.

    Upstream is not wrong to re-assert: the comment above LOOK_DEFAULTS gives
    the reason (an Outline set in cartoon surviving into richardson, whose
    panel hides the slider, "leaving no way to see or fix it"). Detail is just
    not a style-owned quantity for a bundle -- it is a property of the
    document, chosen once by whoever rendered it. So the page re-asserts it
    after the switch, and a reader who drags the slider outranks the file.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    keep = report["detailKeeper"]
    assert keep["configured"] == 8
    assert keep["atStart"] == 8
    assert keep["afterRichardson"] == 8          # not the preset's 4
    assert keep["afterDragThenSwitch"] == 3      # the reader's, not ours


def test_the_style_panel_keeps_four_controls_and_folds_the_rest(tmp_path):
    """src/parts/panel.js builds this panel from ONE table of rows shared with
    py2dmol.solab.org -- its own header calls it "ONE PANEL, TWO SKINS" -- so
    the website does not expose fewer controls and there is no simpler set to
    borrow. What differs is the room: the site is a wide sidebar, ours floats
    over the canvas, and seventeen controls there read as clutter.

    So the four a reader of a bundle has a question for stay out front, and
    the rest fold into a closed <details>. NOTHING IS REMOVED: the assertion
    is that front plus folded is still every control, because a fold that
    dropped a row would look exactly like a fold that worked.

    Whole rows move. parts/ui.js hides per-style controls by setting `hidden`
    on the row carrying data-style, so a control lifted out of its row would
    survive into a style that cannot honour it. Detail is the exception and is
    safe: it is a `.half` with its own data-style, moving into the row that is
    always visible.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    fold = report["styleFold"]
    assert fold["front"] == ["styleSelect", "detailSlider",
                             "colorSelect", "selectionMarkSelect"]
    assert fold["summary"] == "Fine tuning"
    assert fold["closed"] is True        # no `open` attribute was set
    assert fold["foldCount"] == 1
    assert fold["detailIsFront"] is True
    everything = {
        "styleSelect", "lineWidthSlider", "outlineWidthSlider", "thicknessSlider",
        "sheetFlatSlider", "highlightSlider", "shadeSlider", "pencilSlider",
        "shadowSlider", "outlineTintSlider", "detailSlider", "orthoSlider",
        "colorSelect", "selectionMarkSelect", "cartoonArrowsToggle", "basePlatesToggle",
    }
    # ...plus the one control the panel does not ship: which painter draws.
    # parts/ui.js leaves it out because a notebook build carries one painter
    # and the flag would ask for a file that is absent; a bundle carries both,
    # and the two do not draw a close-up side chain the same way.
    assert set(fold["front"]) | set(fold["folded"]) == everything | {"pinsp-gpu"}
    assert not set(fold["front"]) & set(fold["folded"])
    assert "pinsp-gpu" in fold["folded"]
    assert fold["gpuBefore"] is True
    assert fold["gpuAfterOff"] is False
    assert fold["gpuAfterOn"] is True
    assert fold["invalidated"] == 1      # the mesh is dropped, not reused


def test_a_layer_is_a_selection_so_the_tools_can_act_on_it(tmp_path):
    """The tools act on a selection and a layer is already a named set of
    residues, so the two halves of the page were not connected: side chains
    could only be drawn for residues picked one at a time, and the interesting
    unit is "this whole epitope".

    `sel` on a row takes one layer, shift adds, and the header's `Select`
    takes every TICKED layer -- the same visibility the structure is painted
    from, so what gets selected is what is on screen.
    """
    report = _run_dom_harness(_sequence_bundle(tmp_path), _SEQUENCE_SPEC)
    pick = report["layerSelect"]
    # layer-a is residue 2 (position 1), layer-b residue 5 (position 4).
    assert pick["one"] == [1]
    assert pick["marked"] == [1]
    assert pick["readout"].startswith("A 2")
    assert pick["replaced"] == [4]              # a plain click replaces
    assert pick["added"] == [1, 4]              # shift adds
    assert pick["allVisible"] == [1, 4]
    assert pick["visibleAfterUntick"] == [4]    # an unticked layer is not in it
    # ...and that is the set the one Side chains button then draws.
    assert pick["sidechainsForLayers"] == [1, 4]

    text = Path(_sequence_bundle(tmp_path)).read_text()
    assert 'data-select="layer-a"' in text
    assert 'id="pinsp-layers-select"' in text


def test_a_morph_carries_its_side_chains_one_rotamer_set_per_state(tmp_path):
    """The conformer path builds its frames from Cα coordinates, so the Side
    chains button was dead on exactly the bundles most worth interrogating.

    It works because the browser's table is not world-space: each atom is
    three coefficients in its residue's own backbone frame, rebuilt from the
    final positions at draw time, and `parts/ui.js` builds a SEPARATE table
    per frame -- so one row set per conformer gives each state its own
    rotamers, and the interpolation buffer borrows the table of the state it
    is leaving.

    The fit has to be applied to the atoms: a conformer's frame holds its
    Kabsch-fitted trace while its file holds the original coordinates, and
    measuring a side chain against a backbone that has rotated out from under
    it puts it in the wrong place.
    """
    cif = tmp_path / "ref.cif"
    _atom_fixture(cif, n=8)
    other = tmp_path / "conf1.cif"
    _atom_fixture(other, n=8)
    # Rotate the conformer bodily: the superposition removes it, so the fitted
    # trace and the file's own coordinates are genuinely different frames --
    # which is what makes the transform load-bearing rather than decorative.
    import gemmi
    st = gemmi.read_structure(str(other))
    rot = gemmi.Mat33([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    v = rot.multiply(atom.pos)
                    atom.pos = gemmi.Position(v.x + 17.0, v.y - 9.0, v.z + 4.0)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(other))

    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "out"), extras=True,
        conformers=[{"path": str(other), "label": "rotated"}],
        morph_mapping="exact", morph_reference_label="start")
    state = json.loads(Path(next(i["path"] for i in result["artifacts"]
                                 if i["kind"] == "viewer_state")).read_text())
    assert state["sidechains"] is True
    frames = state["viewer"]["objects"][0]["frames"]
    assert len(frames) == 2
    # One row per residue with a side chain, in BOTH frames -- not inherited.
    for frame in frames:
        rows = frame["sidechain_atoms"]
        assert len(rows) == 8
        assert [row[0] for row in rows] == list(range(8))
        assert all(row[2] and row[2][0][0] == "CB" for row in rows)
    # THE FIT IS APPLIED. The conformer was rotated bodily and superposed
    # back, so its stored side-chain atoms must sit near its stored trace --
    # measured against the raw file they would be a whole rigid transform
    # away. Each CB is a known 1.5 A off its own CA in the fixture.
    import numpy as np
    for frame in frames:
        trace = np.asarray(frame["coords"], dtype=float)
        for row in frame["sidechain_atoms"]:
            cb = np.asarray(row[2][0][1:4], dtype=float)
            assert float(np.linalg.norm(cb - trace[row[0]])) < 2.0

    # ...and the page hands the animation buffer the table of the state it is
    # leaving, rather than always the reference's.
    html = Path(next(i["path"] for i in result["artifacts"] if i["kind"] == "html"))
    assert "buf.sidechains=from.sidechains" in html.read_text()


def test_the_conformation_ring_wraps_instead_of_running_off_the_stage(tmp_path):
    """A label you can read but cannot press is worse than a second line.

    The stage clips its overflow, so a ring that cannot wrap put the later
    conformations off the edge with their hit areas over the Layers panel.
    """
    _, html = _morph_bundle(tmp_path, n_conformers=3)
    text = html.read_text()
    assert ".bm-ring{display:inline-flex;flex-wrap:wrap;" in text
    assert "max-width:100%}" in text
    # The old reservation was for a control that has since been removed.
    assert "calc(100% - 120px)" not in text


def test_the_partners_button_is_not_swept_into_the_conformation_ring(tmp_path):
    """It wears .bm-btn so it looks of a piece with the ring; it must not be
    wired like one. `Number(null)` is 0, so a bare `.bm-btn` query made
    pressing "partners" also morph to the first conformation."""
    _, html = _two_partner_bundle(tmp_path)
    text = html.read_text()
    assert "document.querySelectorAll('.bm-btn')" not in text
    assert text.count("document.querySelectorAll('.bm-btn[data-conf]')") == 2
    # The partners button carries no data-conf, which is what excludes it.
    partners_tag = text.split('id="pinsp-partners"')[0].rsplit("<button", 1)[1]
    assert "data-conf" not in partners_tag


_SUBUNIT = ["ALA", "GLY", "SER", "VAL", "LEU", "THR", "ILE", "PRO", "GLU", "LYS",
            "ASP", "ARG"]
_RECEPTOR = ["TRP", "TYR", "PHE", "HIS", "CYS", "MET", "ASN", "GLN", "TRP", "TYR",
             "PHE", "HIS"]


def _oligomer_cif(path, sequences):
    """One file, one chain per {chain_id: [resname, ...]} entry, laid side by side."""
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    for offset, (chain_id, sequence) in enumerate(sequences.items()):
        chain = Chain.Chain(chain_id)
        for index, resname in enumerate(sequence, start=1):
            residue = Residue.Residue((" ", index, " "), resname, " ")
            residue.add(Atom.Atom("CA", (float(index * 3), float(index % 2) + 15.0 * offset,
                                         0.0), 0.0, 1.0, " ", "CA", index, element="C"))
            chain.add(residue)
        model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def _oligomer_bundle(tmp_path, sequences, partner_chains, **kwargs):
    """One annotation on chain A of a multi-chain file, rendered with partners."""
    cif = _oligomer_cif(tmp_path / "oligomer.cif", sequences)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "site", "kind": "custom", "label": "Epitope",
             "layer_id": "epitope", "layer_label": "Epitope", "color": "#16a34a",
             "residue": {"component_id": "target", "canonical_position": 3,
                         "chain_id": "A", "author_residue_number": 3},
             "resolved": True, "evidence_ids": [], "method": "model"},
        ],
    }
    result = render_inspection_bundle(
        mmcif_path=str(cif), mmcif_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
        inspection_manifest=manifest, output_dir=str(tmp_path / "out"),
        partner_chains=partner_chains, partner_label="other protomers",
        **kwargs)
    html = next(item["path"] for item in result["artifacts"] if item["kind"] == "html")
    return result, Path(html)


def test_a_homo_oligomer_wears_the_annotation_on_every_copy(tmp_path):
    """Three protomers are one target presenting the same site three times."""
    result, html = _oligomer_bundle(
        tmp_path, {"A": _SUBUNIT, "B": _SUBUNIT, "C": _SUBUNIT}, ["B", "C"])
    assert result["partners"]["painted"] == ["B", "C"]
    painted = {(row["chain_id"], row["author_residue_number"]): row
               for row in read_inspection_bundle(str(html))["residues"]
               if row["layers"]}
    # The analysis named chain A alone; the reader sees the site on all three.
    assert sorted(painted) == [("A", 3), ("B", 3), ("C", 3)]
    assert {row["color"] for row in painted.values()} == {"#16a34a"}
    assert {tuple(row["layers"]) for row in painted.values()} == {("epitope",)}


def test_a_binding_partner_is_not_painted_with_the_targets_epitope(tmp_path):
    """A receptor is not a copy of the target, so it keeps its own colour."""
    result, html = _oligomer_bundle(
        tmp_path, {"A": _SUBUNIT, "R": _RECEPTOR}, ["R"])
    assert result["partners"]["painted"] is None
    coloured = {row["chain_id"] for row in read_inspection_bundle(str(html))["residues"]
                if row["layers"]}
    assert coloured == {"A"}


def test_painting_the_copies_can_be_forced_and_refused(tmp_path):
    """`auto` decides by sequence; True and False take the decision back."""
    one, two = tmp_path / "forced", tmp_path / "refused"
    one.mkdir()
    two.mkdir()
    forced, _ = _oligomer_bundle(one, {"A": _SUBUNIT, "R": _RECEPTOR},
                                 ["R"], layers_on_partners=True)
    assert forced["partners"]["painted"] == ["R"]
    refused, _ = _oligomer_bundle(two, {"A": _SUBUNIT, "B": _SUBUNIT}, ["B"],
                                  layers_on_partners=False)
    assert refused["partners"]["painted"] is None


def test_the_partner_button_is_named_by_the_conformation_on_screen(tmp_path):
    """A bound receptor captioned 'other protomers' contradicts its own picture."""
    apo = _oligomer_cif(tmp_path / "apo.cif", {"A": _SUBUNIT, "B": _SUBUNIT})
    bound = _oligomer_cif(tmp_path / "bound.cif", {"A": _SUBUNIT, "R": _RECEPTOR})
    result = render_inspection_bundle(
        mmcif_path=str(apo), mmcif_sha256=hashlib.sha256(apo.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "out"),
        partner_chains=["B"], partner_label="other protomers",
        conformers=[{"path": str(bound), "label": "Bound", "partner_chains": ["R"],
                     "partner_label": "receptor"}],
        morph_mapping="exact", morph_steps=6)
    assert result["partners"]["labels"] == {"0": "other protomers", "1": "receptor"}
    html = Path(next(item["path"] for item in result["artifacts"]
                     if item["kind"] == "html")).read_text()
    # The button opens on the reference's word for its partners and is renamed
    # from the same map as the conformation changes.
    assert '>other protomers</button>' in html
    assert "function partnerLabelAt(" in html
    assert "btn.textContent=name" in html


def test_a_conformer_without_its_own_label_keeps_the_global_one(tmp_path):
    """Naming one state must not leave the others anonymous."""
    apo = _oligomer_cif(tmp_path / "apo.cif", {"A": _SUBUNIT, "B": _SUBUNIT})
    bound = _oligomer_cif(tmp_path / "bound.cif", {"A": _SUBUNIT, "R": _RECEPTOR})
    result = render_inspection_bundle(
        mmcif_path=str(apo), mmcif_sha256=hashlib.sha256(apo.read_bytes()).hexdigest(),
        inspection_manifest={"schema_version": "protein-inspector-manifest-1",
                             "annotations": []},
        output_dir=str(tmp_path / "out"),
        partner_chains=["B"], partner_label="partners",
        conformers=[{"path": str(bound), "label": "Bound", "partner_chains": ["R"]}],
        morph_mapping="exact", morph_steps=6)
    assert result["partners"]["labels"] == {"0": "partners", "1": "partners"}
