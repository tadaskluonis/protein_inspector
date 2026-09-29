"""The sidecar's contract: simple layers in, correct canonical positions out."""
import ast
import importlib.util
import sys
from pathlib import Path

import pytest
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

SKILL_DIR = Path(__file__).resolve().parents[1]
REPO = SKILL_DIR.parent


def _load():
    sys.path.insert(0, str(REPO))
    spec = importlib.util.spec_from_file_location("bindos_skill_kernel", SKILL_DIR / "kernel.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _cif(path, first=41, n=6, bend=0.0):
    """Author numbering deliberately does NOT start at 1."""
    structure = Structure.Structure("f")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for offset in range(n):
        number = first + offset
        residue = Residue.Residue((" ", number, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", (float(offset * 3), float(offset % 2) + bend * offset, 0.0),
                              0.0, 1.0, " ", "CA", offset + 1, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def test_kernel_py_passes_the_sidecar_structural_gate():
    """Loading a skill may only DEFINE names: no top-level calls, no decorators."""
    tree = ast.parse((SKILL_DIR / "kernel.py").read_text())
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            assert not any(a.name == "*" for a in getattr(node, "names", []))
            continue
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            assert not node.decorator_list, f"{node.name} is decorated"
            for default in node.args.defaults + [d for d in node.args.kw_defaults if d]:
                assert isinstance(default, ast.Constant), f"{node.name} has a non-literal default"
            assert not node.name.startswith("_"), node.name
            continue
        if isinstance(node, ast.Assign):
            assert all(isinstance(t, ast.Name) and not t.id.startswith("_") for t in node.targets)
            assert isinstance(node.value, (ast.Constant, ast.Tuple)), "non-literal top-level assign"
            continue
        assert isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant), \
            f"disallowed top-level node: {type(node).__name__}"


def test_canonical_position_is_derived_not_the_author_number(tmp_path):
    """Author 43 is canonical 3. Asking the caller for this is how bundles lie."""
    kernel = _load()
    cif = _cif(tmp_path / "t.cif")
    manifest, missing = kernel.build_manifest(
        str(cif), [{"id": "hot", "label": "Hot", "color": "#dc2626", "residues": [41, 43, 46]}])
    assert missing == []
    positions = {a["residue"]["author_residue_number"]: a["residue"]["canonical_position"]
                 for a in manifest["annotations"]}
    assert positions == {41: 1, 43: 3, 46: 6}


def test_residues_absent_from_the_model_are_reported_not_dropped(tmp_path):
    kernel = _load()
    cif = _cif(tmp_path / "t.cif")
    manifest, missing = kernel.build_manifest(
        str(cif), [{"id": "hot", "label": "Hot", "residues": [41, 999]}])
    assert len(manifest["annotations"]) == 1
    assert missing == [("hot", 999)]


def test_notes_and_ranges_both_work_as_residue_input(tmp_path):
    kernel = _load()
    cif = _cif(tmp_path / "t.cif")
    manifest, _ = kernel.build_manifest(str(cif), [
        {"id": "dom", "label": "Domain", "residues": range(41, 44)},
        {"id": "ep", "label": "Epitope", "residues": {45: "rim contact"}},
    ])
    labels = {a["annotation_id"]: a["label"] for a in manifest["annotations"]}
    assert labels["dom-A41"] == "Domain"
    assert labels["ep-A45"] == "rim contact"


def test_end_to_end_bundle_reads_back_with_the_layers_it_was_given(tmp_path):
    kernel = _load()
    cif = _cif(tmp_path / "t.cif")
    report = kernel.inspect_structure(
        str(cif),
        layers=[{"id": "ep", "label": "Epitope", "color": "#dc2626", "residues": [43, 44]}],
        out=str(tmp_path / "out.html"))
    assert report["unmapped_residues"] == []
    rows = kernel.inspection_table(report["path"])
    painted = {r["author_residue_number"]: r["color"] for r in rows if "ep" in r["layers"]}
    assert painted == {43: "#dc2626", 44: "#dc2626"}


def test_conformers_ship_endpoints_only(tmp_path):
    kernel = _load()
    reference = _cif(tmp_path / "a.cif")
    report = kernel.inspect_structure(
        str(reference), layers=[],
        conformers=[{"path": str(_cif(tmp_path / "b.cif", bend=1.5)), "label": "Bent"}],
        morph_steps=10, out=str(tmp_path / "m.html"))
    assert report["morph"]["frames_in_file"] == 2
    assert report["morph"]["stored_intermediates"] == 0
    assert report["morph"]["animation_steps"] == 10
    assert report["morph"]["conformers"][1]["rmsd_to_reference_A"] > 0


def test_pdb_input_is_converted(tmp_path):
    pytest.importorskip("gemmi")
    kernel = _load()
    import gemmi
    cif = _cif(tmp_path / "t.cif")
    pdb = tmp_path / "t.pdb"
    gemmi.read_structure(str(cif)).write_pdb(str(pdb))
    assert kernel.to_mmcif(str(pdb)).endswith(".cif")
