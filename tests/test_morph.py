"""Conformational morphs: residue mapping, superposition, interpolation.

The morph is a DEPICTION OF ITS ENDPOINTS, not a pathway -- intermediates are
Cartesian interpolations, so they are not physical and bond geometry is not
preserved. Everything here tests that the endpoints are exactly the structures
they claim to be and that the residues paired between them are the right ones.

The residue-mapping guard carries most of the weight. A morph built over a
quiet intersection of two residue sets looks exactly as smooth as a correct
one while interpolating the wrong pairs, so a mismatch has to be refused
loudly rather than worked around.
"""

import hashlib
import json
import pathlib

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue, Structure
from Bio.PDB.mmcifio import MMCIFIO

from protein_inspector import render_inspection_bundle
from protein_inspector.renderer import _ca_trace, _morph_frames, _superpose


def _write(path, coords, chain_id="A", numbers=None, resname="ALA", extra_atom=False):
    """A single-chain CA-only structure at the given coordinates."""
    numbers = numbers if numbers is not None else range(1, len(coords) + 1)
    structure = Structure.Structure("s")
    model = Model.Model(0)
    chain = Chain.Chain(chain_id)
    for number, xyz in zip(numbers, coords):
        residue = Residue.Residue((" ", int(number), " "), resname, " ")
        residue.add(Atom.Atom("CA", tuple(float(v) for v in xyz), 0.0, 1.0, " ",
                              "CA", int(number), element="C"))
        if extra_atom:
            residue.add(Atom.Atom("CB", tuple(float(v) + 1.0 for v in xyz), 0.0, 1.0,
                                  " ", "CB", int(number) + 1000, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def _helix(n=12, twist=0.6, rise=1.5, radius=4.0, phase=0.0):
    t = np.arange(n) * twist + phase
    return np.stack([radius * np.cos(t), radius * np.sin(t), np.arange(n) * rise], 1)


# --- residue mapping -----------------------------------------------------

def test_ca_trace_keys_on_chain_and_author_number(tmp_path):
    trace = _ca_trace(_write(tmp_path / "a.cif", _helix(4), chain_id="B",
                             numbers=[10, 11, 12, 13]))
    assert sorted(trace) == [("B", 10), ("B", 11), ("B", 12), ("B", 13)]
    assert all(name == "ALA" and len(xyz) == 3 for name, xyz in trace.values())


def test_ca_trace_takes_one_atom_per_residue(tmp_path):
    """Sidechain atoms must not enter the trace, or the frames carry more
    points than the colour array has residues."""
    assert len(_ca_trace(_write(tmp_path / "a.cif", _helix(5), extra_atom=True))) == 5


def test_ca_trace_refuses_a_duplicate_residue_address(tmp_path):
    """An insertion code or altloc collapsing two residues onto one address
    makes the pairing ambiguous, and ambiguity here is silent."""
    path = tmp_path / "dup.cif"
    coords = _helix(3)
    structure = Structure.Structure("s")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    # 3 and 3A are different residues to Biopython and the same address to us.
    for number, icode, xyz in [(1, " ", coords[0]), (3, " ", coords[1]),
                               (3, "A", coords[2])]:
        residue = Residue.Residue((" ", number, icode), "ALA", " ")
        residue.add(Atom.Atom("CA", tuple(float(v) for v in xyz), 0.0, 1.0, " ",
                              "CA", 1, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    with pytest.raises(ValueError, match="twice|ambiguous"):
        _ca_trace(path)


def test_morph_refuses_a_mismatched_residue_set_and_names_the_offenders(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(8), numbers=range(1, 9))
    short = _write(tmp_path / "short.cif", _helix(6, phase=0.3), numbers=range(1, 7))
    with pytest.raises(ValueError) as excinfo:
        _morph_frames(ref, [{"path": str(short)}], steps=4)
    message = str(excinfo.value)
    assert "not exact" in message
    # The specific residues, not just a count: "they differ" is not actionable.
    assert "('A', 7)" in message and "('A', 8)" in message


def test_morph_refuses_a_conformer_whose_hash_does_not_match(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(6))
    other = _write(tmp_path / "other.cif", _helix(6, phase=0.4))
    with pytest.raises(ValueError, match="hash does not match"):
        _morph_frames(ref, [{"path": str(other), "sha256": "0" * 64}], steps=2)


def test_morph_accepts_a_matching_hash(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(6))
    other = _write(tmp_path / "other.cif", _helix(6, phase=0.4))
    digest = hashlib.sha256(other.read_bytes()).hexdigest()
    _, _, frames, _ = _morph_frames(ref, [{"path": str(other), "sha256": digest}], steps=2)
    assert len(frames) == 3


# --- superposition -------------------------------------------------------

def test_superpose_recovers_a_rigid_transform(tmp_path):
    target = _helix(20)
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1.0]])
    moved = target @ rotation.T + np.array([13.0, -4.0, 7.5])
    assert np.allclose(_superpose(moved, target), target, atol=1e-6)


def test_superpose_forbids_reflection(tmp_path):
    """A mirrored chain is not the same chain. Kabsch without the determinant
    correction will happily flip one onto the other and report RMSD 0."""
    target = _helix(20)
    mirrored = target * np.array([1.0, 1.0, -1.0])
    rmsd = np.sqrt(((_superpose(mirrored, target) - target) ** 2).sum(1).mean())
    assert rmsd > 1.0, "reflection was allowed, so a mirror image superposes perfectly"


# --- interpolation -------------------------------------------------------

def test_endpoints_are_the_structures_not_interpolations(tmp_path):
    """First and last frame must be the inputs exactly (the last after
    superposition), or the morph never actually shows either conformer."""
    a, b = _helix(10), _helix(10, phase=0.9)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    keys, names, frames, report = _morph_frames(ref, [{"path": str(other)}], steps=5)
    # mmCIF carries three decimals, so the round trip is exact to ~1e-3 A,
    # not to machine precision. The point is that the endpoint is the
    # structure rather than an interpolation of it.
    assert np.allclose(frames[0], a, atol=5e-3)
    assert np.allclose(frames[-1], _superpose(b, a), atol=5e-3)
    assert report["n_frames"] == len(frames) == 5 + 1
    assert report["n_residues"] == len(keys) == len(names) == 10


def test_interpolation_is_linear_in_time(tmp_path):
    a, b = _helix(8), _helix(8, phase=1.1)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    _, _, frames, _ = _morph_frames(ref, [{"path": str(other)}], steps=2)
    assert np.allclose(frames[1], (frames[0] + frames[2]) / 2.0, atol=1e-6)


def test_multiple_conformers_chain_end_to_end(tmp_path):
    a, b, c = _helix(9), _helix(9, phase=0.5), _helix(9, phase=1.0)
    ref = _write(tmp_path / "a.cif", a)
    paths = [{"path": str(_write(tmp_path / "b.cif", b)), "label": "mid"},
             {"path": str(_write(tmp_path / "c.cif", c))}]
    _, _, frames, report = _morph_frames(ref, paths, steps=4)
    assert report["conformers"] == ["reference", "mid", "c"]
    assert report["n_frames"] == len(frames) == 4 * 2 + 1
    assert len(report["endpoint_rmsd_A"]) == 2


def test_reported_rmsd_is_the_superposed_endpoint_distance(tmp_path):
    a, b = _helix(15), _helix(15, phase=0.8)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    _, _, _, report = _morph_frames(ref, [{"path": str(other)}], steps=3)
    expected = np.sqrt(((_superpose(b, a) - a) ** 2).sum(1).mean())
    assert report["endpoint_rmsd_A"] == [round(float(expected), 2)]


# --- the property the viewer depends on ----------------------------------

def test_morph_does_not_disturb_the_colouring(tmp_path):
    """THE POINT OF THE FEATURE. Colour is per-residue and frames carry only
    coordinates, so the annotation layers must paint identically with a morph
    attached and without one. If this fails, the morph is repainting."""
    a, b = _helix(6), _helix(6, phase=0.7)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "p", "kind": "custom", "label": "probe",
             "color": "#dc2626",
             "residue": {"component_id": "target", "canonical_position": 3,
                         "chain_id": "A", "author_residue_number": 3},
             "resolved": True, "evidence_ids": ["e"], "method": "manual"},
        ],
    }
    digest = hashlib.sha256(ref.read_bytes()).hexdigest()

    def build(name, morph):
        render_inspection_bundle(
            mmcif_path=str(ref), mmcif_sha256=digest, inspection_manifest=manifest,
            output_dir=str(tmp_path / name), output_name="b", extras=True,
            **({"morph_to": morph, "morph_steps": 3} if morph else {}))
        state = json.loads((tmp_path / name / "b.viewer.json").read_text())
        return state["viewer"]["objects"][0]

    plain = build("plain", None)
    morphed = build("morphed", [{"path": str(other), "label": "open"}])

    assert morphed["color"] == plain["color"], "morph changed the colour assignment"
    assert len(plain["frames"]) == 1 and len(morphed["frames"]) == 4
    # ...and the extra frames are coordinates only: same residues, same order.
    assert morphed["frames"][0]["chains"] == plain["frames"][0]["chains"]
    assert (morphed["frames"][0]["residue_numbers"]
            == plain["frames"][0]["residue_numbers"])
