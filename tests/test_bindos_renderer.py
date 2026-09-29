import hashlib
import json
import pytest
from pathlib import Path

from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

from bindos_structure_inspector import read_inspection_bundle, render_inspection_bundle


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
        "schema_version": "bindos-inspection-manifest-1",
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
        "schema_version": "bindos-inspection-manifest-1",
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
        "schema_version": "bindos-inspection-manifest-1",
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
    assert 'id="bindos-base-mode"' not in page and 'id="bindos-base-color"' not in page


def test_default_render_is_one_file_with_about_tabs(tmp_path):
    """A bundle is one file. Extra context becomes a TAB, never a sibling file.

    The About body is agent-authored and often quotes fetched text, so it is
    sanitised: the element and its contents go, not just the tags.
    """
    cif = tmp_path / "fixture.cif"
    _fixture(cif)
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "bindos-inspection-manifest-1",
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
        "schema_version": "bindos-inspection-manifest-1",
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
            "schema_version": "bindos-inspection-manifest-1",
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
            inspection_manifest={"schema_version": "bindos-inspection-manifest-1",
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
        "schema_version": "bindos-inspection-manifest-1",
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

    assert state["inspector_version"].startswith("bindos-inspector-")
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
    with pytest.raises(ValueError, match="no bindos inspection state"):
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
        inspection_manifest={"schema_version": "bindos-inspection-manifest-1",
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
        inspection_manifest={"schema_version": "bindos-inspection-manifest-1", "annotations": []},
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
    state = json.loads(html.read_text().split('id="bindos-inspection-state" type="application/json">')[1]
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
    assert ".bindos-stage[data-morph='1'] #controlsContainer{display:none!important}" in text


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
            inspection_manifest={"schema_version": "bindos-inspection-manifest-1", "annotations": []},
            output_dir=str(tmp_path / "out"),
            conformers=[{"path": str(cif)}], morph_to=[{"path": str(cif)}])


def test_the_viewer_cannot_paint_over_the_layer_panel(tmp_path):
    """py2Dmol fixes .py2dmol-viewer-instance at 948px; the stage must clamp it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert ".bindos-stage{position:relative;min-width:0;overflow-x:auto" in text
    assert "@media (max-width:1340px)" in text
    assert 'id="bindos-controls"' in text
    # Orient/Focus/Rotate/Style/Clip/Capture float over the top-right of the
    # canvas instead of sitting in a 340px column beside it.
    assert ".bindos-stage #rightPanelContainer{position:absolute!important;top:10px;right:10px;" in text
    assert ".bindos-stage .py2dmol-viewer-instance{width:auto!important}" in text
    report = _run_dom_harness(html, {"chains": ["A"] * 5, "residueNumbers": [1, 2, 3, 4, 5]})
    assert report["controlsAfterClick"] == "0"


def test_a_blocked_download_still_leaves_the_capture_recoverable(tmp_path):
    """Save Image reports success without checking; an embedded frame may drop it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert 'id="bindos-capture"' in text
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
        inspection_manifest={"schema_version": "bindos-inspection-manifest-1", "annotations": []},
        output_dir=str(tmp_path / "out"), partner_chains=["P"], partner_label="Fab",
        conformers=conformers, morph_mapping="exact", morph_steps=6)
    html = next(item["path"] for item in result["artifacts"] if item["kind"] == "html")
    return result, Path(html)


def test_partner_chains_get_a_toggle_that_hides_exactly_them(tmp_path):
    _, html = _partner_bundle(tmp_path)
    text = html.read_text()
    assert 'id="bindos-partners"' in text
    assert "Show Fab" in text
    report = _run_dom_harness(html, {"chains": ["A", "A", "A", "A", "P", "P", "P", "P"],
                                     "residueNumbers": [1, 2, 3, 4, 1, 2, 3, 4]})
    assert report["partners"]["mode"] == "patch"
    assert report["partners"]["shown"] == [0, 1, 2, 3] == report["partners"]["expected"]
    assert report["partners"]["chainKeys"] == ["obj:A"]
    # Re-ticking restores everything rather than leaving a narrowed selection.
    assert report["partners"]["afterRetick"] == "all"


def test_a_partner_is_excluded_from_the_morph_and_held_still(tmp_path):
    """The partner exists in the reference only; interpolating it is meaningless."""
    result, _ = _partner_bundle(tmp_path, morph=True)
    morph = result["morph"]
    # The residue mapping is the target chain alone -- the conformer has no P.
    assert morph["n_residues"] == 4
    assert morph["mapping"] == "exact"
    # ...but the partner is still drawn, identically in every frame.
    assert morph["partner_positions_static"] == 4


def test_every_export_retains_the_upstream_licence_and_credit(tmp_path):
    """BEER-WARE asks one thing: retain the notice. Minification removes it."""
    html = _two_layer_bundle(tmp_path)
    text = Path(html).read_text()
    assert "THE BEER-WARE LICENSE" in text
    assert "Sergey Ovchinnikov" in text
    assert "https://github.com/sokrypton/py2Dmol" in text
    assert "78c2d489d0b5c5d19accd9eeeef878c2868f5271" in text
    # Not only in a comment: a reader of the page can see who wrote the viewer.
    assert 'class="bindos-credit"' in text
    assert "py2Dmol</a> by Sergey Ovchinnikov" in text
