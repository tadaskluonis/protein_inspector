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
            {"annotation_id": "x", "kind": "partner_contact", "label": "Partner X contact", "partner_id": "Partner X", "residue": {"component_id": "target", "canonical_position": 2}, "resolved": True, "evidence_ids": ["e1"], "method": "distance"},
            {"annotation_id": "y", "kind": "partner_contact", "label": "Partner Y contact", "partner_id": "Partner Y", "residue": {"component_id": "target", "canonical_position": 4}, "resolved": True, "evidence_ids": ["e2"], "method": "distance"},
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
